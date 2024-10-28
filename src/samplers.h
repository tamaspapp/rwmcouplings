#include <cmath>
#include <Eigen/Core>
#include <Eigen/Cholesky>

#include "targets.h"
#include "rng.h"

using Eigen::Lower;
using Eigen::Upper;

using Eigen::MatrixBase;
using Eigen::MatrixXd;
using Eigen::VectorXd;

/************* Helper functions -- linear algebra *************/

// Cholesky factorization
inline Eigen::MatrixXd chol_l(const Eigen::Ref<const Eigen::MatrixXd> &Sigma)
{
    if (Sigma.cols() == 1)
        return Sigma.cwiseSqrt();
    else
        return Sigma.llt().matrixL();
}

// Multiply by dense matrix
inline Eigen::VectorXd precondition(const Eigen::Ref<const Eigen::MatrixXd> &Sigma, const Eigen::Ref<const Eigen::VectorXd> &x)
{
    if (Sigma.cols() == 1)
        return Sigma.cwiseProduct(x);
    else
        return Sigma * x;
}

// Multiply by triangular matrices
inline Eigen::VectorXd precondition_L(const Eigen::Ref<const Eigen::MatrixXd> &L, const Eigen::Ref<const Eigen::VectorXd> &x)
{
    if (L.cols() == 1)
        return L.cwiseProduct(x);
    else
        return L.triangularView<Lower>() * x;
}
inline Eigen::VectorXd precondition_Ltop(const Eigen::Ref<const Eigen::MatrixXd> &L, const Eigen::Ref<const Eigen::VectorXd> &x)
{
    if (L.cols() == 1)
        return L.cwiseProduct(x);
    else
        return L.triangularView<Lower>().transpose() * x;
}

// Multiply by inverses of triangular matrices
inline Eigen::VectorXd precondition_Linverse(const Eigen::Ref<const Eigen::MatrixXd> &L, const Eigen::Ref<const Eigen::VectorXd> &x)
{
    if (L.cols() == 1)
        return x.cwiseProduct(L.cwiseInverse());
    else
        return L.triangularView<Lower>().solve(x); // Uses Gaussian elimination, as fast as multiplying by the actual inverse
}

// Reflect z in e (where e is assumed to be normalized)
inline Eigen::VectorXd reflect(const Eigen::Ref<const Eigen::VectorXd> &z, const Eigen::Ref<const Eigen::VectorXd> &e)
{
    return z - (2. * z.dot(e)) * e;
}


/************* Helper functions -- MCMC *************/

// Gaussian proposal, using lower triangular Cholesky factor
inline Eigen::VectorXd gaussian_proposal(const Eigen::Ref<const Eigen::VectorXd> &x,
                                         const Eigen::Ref<const Eigen::MatrixXd> &L,
                                         const Eigen::Ref<const Eigen::VectorXd> &z)
{
    return x + precondition_L(L, z);
}

// Proposal log density difference log q(x | xp) - log q(xp | x) for MALA.
//      - See Proposition 1 in Titsias (2023, NeurIPS).
//      - This refinement saves two matmuls, at the expense of storing two O(d) temporaries x_mean and xp_mean.
inline double mala_proposal_logdensity_difference(const Eigen::Ref<const Eigen::VectorXd> &x, const Eigen::Ref<const Eigen::VectorXd> &x_mean, const Eigen::Ref<const Eigen::VectorXd> &x_grad,
                                                  const Eigen::Ref<const Eigen::VectorXd> &xp, const Eigen::Ref<const Eigen::VectorXd> &xp_mean, const Eigen::Ref<const Eigen::VectorXd> &xp_grad)
{
    return 0.5 * (x - 0.5 * (xp + xp_mean)).dot(xp_grad) - 0.5 * (xp - 0.5 * (x + x_mean)).dot(x_grad);
}


/************* Single-chain MCMC *************/

template <typename m1>
void rwm(const Target &target,
         const Eigen::Ref<const Eigen::VectorXd> &x0,
         const Eigen::Ref<const Eigen::MatrixXd> &Sigma, // Cholesky factor of proposal covariance. Can also input a vector for diagonal covariances.
         const int &iter,
         const int &thin,
         Eigen::MatrixBase<m1> &xs,
         double &acceptance_rate)
{
    int d = x0.size();
    MatrixXd Sigma_chol_l = chol_l(Sigma);

    Eigen::VectorXd x = x0;
    double logpi_x = target.LogDensity(x);
    xs.row(0) = x;

    int accepted = 0;
    for (int it = 1; it <= iter; it++)
    {
        Eigen::VectorXd xp = gaussian_proposal(x, Sigma_chol_l, RNG::GaussianVector(d));
        double logpi_xp = target.LogDensity(xp);

        // Metropolis correction
        double logHR_x = logpi_xp - logpi_x;
        double log_u = log(RNG::runif(RNG::rng));
        if (logHR_x > 0 || log_u < logHR_x)
        {
            ++accepted;
            x = xp;
            logpi_x = logpi_xp;
        }

        if (it % thin == 0)
            xs.row(it / thin) = x;
    }
    acceptance_rate = (double)accepted / (double)iter;
}

template <typename m1>
void mala(const Target &target,
          const Eigen::Ref<const Eigen::VectorXd> &x0,
          const Eigen::Ref<const Eigen::MatrixXd> &Sigma,
          const int &iter,
          const int &thin,
          Eigen::MatrixBase<m1> &xs,
          double &acceptance_rate)
{
    int d = x0.size();
    MatrixXd Sigma_chol_l = chol_l(Sigma);

    Eigen::VectorXd x = x0;
    Eigen::VectorXd x_grad = target.GradLogDensity(x);
    Eigen::VectorXd x_mean = x + 0.5 * precondition(Sigma, x_grad);
    double logpi_x = target.LogDensity(x);
    xs.row(0) = x;

    int accepted = 0;
    for (int it = 1; it <= iter; it++)
    {
        Eigen::VectorXd xp = gaussian_proposal(x_mean, Sigma_chol_l, RNG::GaussianVector(d));
        Eigen::VectorXd xp_grad = target.GradLogDensity(xp);
        Eigen::VectorXd xp_mean = xp + 0.5 * precondition(Sigma, xp_grad);
        double logpi_xp = target.LogDensity(xp);

        // Metropolis correction
        double logHR_x = logpi_xp - logpi_x + mala_proposal_logdensity_difference(x, x_mean, x_grad, xp, xp_mean, xp_grad);
        double log_u = log(RNG::runif(RNG::rng));
        if (logHR_x > 0 || log_u < logHR_x)
        {
            ++accepted;
            x = xp;
            x_grad = xp_grad;
            x_mean = xp_mean;
            logpi_x = logpi_xp;
        }

        if (it % thin == 0)
            xs.row(it / thin) = x;
    }
    acceptance_rate = (double)accepted / (double)iter;
}

/************* Coupled MCMC, same marginal kernel **********/

// Two-scale GCRN coupling, with reflection-maximal coupling when: thresh > || L.inverse() * (x - y) ||^2.
void rwm_twoscaleGCRN(Target &target,
                      const Eigen::Ref<const Eigen::VectorXd> &x0,
                      const Eigen::Ref<const Eigen::VectorXd> &y0,
                      const Eigen::Ref<const Eigen::MatrixXd> &Sigma,
                      const int &iter,
                      const double &thresh,
                      Eigen::Ref<Eigen::VectorXd> sumdiff, // Fishy function for asymptotic variance calculation
                      int &tau)
{
    int d = x0.size();
    MatrixXd L = chol_l(Sigma);

    // X-chain
    Eigen::VectorXd x = x0, xp(d);
    double logpi_x = target.LogDensity(x), logpi_xp;
    bool update_gx = true;
    // xs.row(0) = x;
    int accepted_x = 0;

    // Y-chain
    Eigen::VectorXd y = y0, yp(d);
    double logpi_y = target.LogDensity(y), logpi_yp;
    bool update_gy = true;
    //ys.row(0) = y;
    int accepted_y = 0;

    // Will need these for the coupling.
    Eigen::VectorXd e_gx(d), e_gy(d), e(d);

    // Fishy function
    sumdiff = Eigen::VectorXd::Zero(d);

    //int iter_done = iter;
    for (int it = 1; it <= iter; it++)
    {
        sumdiff += (x-y);

        bool same_proposal = false, accepted_x_proposal = false, accepted_y_proposal = false;

        e = precondition_Linverse(L, x - y);
        if (e.squaredNorm() < thresh) // Reflection-maximal
        {
            Eigen::VectorXd zx = RNG::GaussianVector(d);
            xp = gaussian_proposal(x, L, zx);
            logpi_xp = target.LogDensity(xp);

            double log_probcpl = -0.5 * (zx + e).squaredNorm() + 0.5 * zx.squaredNorm();
            double log_ucpl = log(RNG::runif(RNG::rng));
            if (log_probcpl > 0 || log_ucpl < log_probcpl) // Coalesce proposals
            {
                yp = xp;
                logpi_yp = logpi_xp;
                same_proposal = true;
            }
            else // Reflect
            {
                e.normalize();
                Eigen::VectorXd zy = reflect(zx, e);

                yp = gaussian_proposal(y, L, zy);
                logpi_yp = target.LogDensity(yp);
            }
        }
        else // GCRN
        {
            if (update_gx)
            {
                e_gx = precondition_Ltop(L, target.GradLogDensity(x));
                e_gx.normalize();
                update_gx = false;
            }
            if (update_gy)
            {
                e_gy = precondition_Ltop(L, target.GradLogDensity(y));
                e_gy.normalize();
                update_gy = false;
            }

            Eigen::VectorXd z = RNG::GaussianVector(d);
            double z1 = RNG::rnorm(RNG::rng);

            Eigen::VectorXd zx = z; zx += (z1 - zx.dot(e_gx)) * e_gx;
            Eigen::VectorXd zy = z; zy += (z1 - zy.dot(e_gy)) * e_gy;

            xp = gaussian_proposal(x, L, zx);
            yp = gaussian_proposal(y, L, zy);

            logpi_xp = target.LogDensity(xp);
            logpi_yp = target.LogDensity(yp);
        }

        double logHR_x = logpi_xp - logpi_x;
        double logHR_y = logpi_yp - logpi_y;

        double log_u = log(RNG::runif(RNG::rng));
        if (logHR_x > 0 || log_u < logHR_x)
        {
            ++accepted_x;
            x = xp;
            logpi_x = logpi_xp;
            accepted_x_proposal = update_gx = true;
        }
        if (logHR_y > 0 || log_u < logHR_y)
        {
            ++accepted_y;
            y = yp;
            logpi_y = logpi_yp;
            accepted_y_proposal = update_gy = true;
        }

        // Break out of the loop if we've coalesced
        if (same_proposal && accepted_x_proposal && accepted_y_proposal)
        {
            tau = it;
            break;
        }
    }
}

// Two-scale GCRefl coupling, with reflection-maximal coupling when: thresh > || L.inverse() * (x - y) ||^2.
void rwm_twoscaleGCRefl(Target &target,
                        const Eigen::Ref<const Eigen::VectorXd> &x0,
                        const Eigen::Ref<const Eigen::VectorXd> &y0,
                        const Eigen::Ref<const Eigen::MatrixXd> &Sigma,
                        const int &iter,
                        const double &thresh,
                        Eigen::Ref<Eigen::VectorXd> sumdiff,
                        int &tau)
{

    int d = x0.size();
    MatrixXd L = chol_l(Sigma);

    // X-chain
    Eigen::VectorXd x = x0, xp(d);
    double logpi_x = target.LogDensity(x), logpi_xp;
    bool update_gx = true;
    // xs.row(0) = x;
    int accepted_x = 0;

    // Y-chain
    Eigen::VectorXd y = y0, yp(d);
    double logpi_y = target.LogDensity(y), logpi_yp;
    bool update_gy = true;
    //ys.row(0) = y;
    int accepted_y = 0;

    // Will need these for the coupling.
    Eigen::VectorXd e(d), e_gx(d), e_gy(d), c_x(d), c_y(d);

    // Fishy function
    sumdiff = Eigen::VectorXd::Zero(d);

    //int iter_done = iter;
    for (int it = 1; it <= iter; it++)
    {
        sumdiff += (x-y);

        bool same_proposal = false, accepted_x_proposal = false, accepted_y_proposal = false;

        e = precondition_Linverse(L, x - y);
        if (e.squaredNorm() < thresh) // Reflection-maximal
        {
            Eigen::VectorXd zx = RNG::GaussianVector(d);
            xp = gaussian_proposal(x, L, zx);
            logpi_xp = target.LogDensity(xp);

            double log_probcpl = -0.5 * (zx + e).squaredNorm() + 0.5 * zx.squaredNorm();
            double log_ucpl = log(RNG::runif(RNG::rng));
            if (log_probcpl > 0 || log_ucpl < log_probcpl) // Coalesce proposals
            {
                same_proposal = true;

                yp = xp;
                logpi_yp = logpi_xp;
            }
            else // Reflect
            {
                e.normalize();
                Eigen::VectorXd zy = reflect(zx, e);

                yp = gaussian_proposal(y, L, zy);
                logpi_yp = target.LogDensity(yp);
            }
        }
        else // GCRefl
        {
            if (update_gx)
            {
                e_gx = precondition_Ltop(L, target.GradLogDensity(x));
                e_gx.normalize();
                update_gx = false;
            }
            if (update_gy)
            {
                e_gy = precondition_Ltop(L, target.GradLogDensity(y));
                e_gy.normalize();
                update_gy = false;
            }

            e.normalize();
            c_x = (e_gx - e_gx.dot(e) * e).normalized(); 
            c_y = (e_gy - e_gy.dot(e) * e).normalized();

            Eigen::VectorXd z = RNG::GaussianVector(d);
            double z1 = RNG::rnorm(RNG::rng);

            Eigen::VectorXd zx = z;             zx += (z1 - zx.dot(c_x)) * c_x;
            Eigen::VectorXd zy = reflect(z, e); zy += (z1 - zy.dot(c_y)) * c_y;

            xp = gaussian_proposal(x, L, zx);
            yp = gaussian_proposal(y, L, zy);

            logpi_xp = target.LogDensity(xp);
            logpi_yp = target.LogDensity(yp);
        }

        double logHR_x = logpi_xp - logpi_x;
        double logHR_y = logpi_yp - logpi_y;

        double log_u = log(RNG::runif(RNG::rng));
        if (logHR_x > 0 || log_u < logHR_x)
        {
            ++accepted_x;
            x = xp;
            logpi_x = logpi_xp;
            accepted_x_proposal = update_gx = true;
        }
        if (logHR_y > 0 || log_u < logHR_y)
        {
            ++accepted_y;
            y = yp;
            logpi_y = logpi_yp;
            accepted_y_proposal = update_gy = true;
        }

        // Break out of the loop if we've coalesced
        if (same_proposal && accepted_x_proposal && accepted_y_proposal)
        {
            tau = it;
            break;
        }
    }
}

// Two-scale CRN coupling, with reflection-maximal coupling when: thresh > || L.inverse() * (x_mean - y_mean) ||^2, where L is the Cholesky factor of the proposal covariance.
void mala_twoscaleCRN(const Target &target,
                      const Eigen::Ref<const Eigen::VectorXd> &x0,
                      const Eigen::Ref<const Eigen::VectorXd> &y0,
                      const Eigen::Ref<const Eigen::MatrixXd> &Sigma,
                      const int &iter,
                      const double &thresh,
                      Eigen::Ref<Eigen::VectorXd> sumdiff,
                      int &tau)
{
    int d = x0.size();
    MatrixXd Sigma_chol_l = chol_l(Sigma);

    // X-chain
    Eigen::VectorXd x = x0, xp(d);
    Eigen::VectorXd x_grad = target.GradLogDensity(x), xp_grad(d);
    Eigen::VectorXd x_mean = x + 0.5 * precondition(Sigma, x_grad), xp_mean(d);
    double logpi_x = target.LogDensity(x), logpi_xp;
    int accepted_x = 0;

    // Y-chain
    Eigen::VectorXd y = y0, yp(d);
    Eigen::VectorXd y_grad = target.GradLogDensity(y), yp_grad(d);
    Eigen::VectorXd y_mean = y + 0.5 * precondition(Sigma, y_grad), yp_mean(d);
    double logpi_y = target.LogDensity(y), logpi_yp;
    int accepted_y = 0;

    // Fishy function
    sumdiff = Eigen::VectorXd::Zero(d);

    for (int it = 1; it <= iter; it++)
    {
        sumdiff += (x-y);

        Eigen::VectorXd e = precondition_Linverse(Sigma_chol_l, x_mean - y_mean);
        bool same_proposal = false, accepted_x_proposal = false, accepted_y_proposal = false;
        if (e.squaredNorm() < thresh) // Reflection-maximal
        {
            Eigen::VectorXd zx = RNG::GaussianVector(d);

            xp = gaussian_proposal(x_mean, Sigma_chol_l, zx);
            logpi_xp = target.LogDensity(xp);
            xp_grad = target.GradLogDensity(xp);
            xp_mean = xp + 0.5 * precondition(Sigma, xp_grad);

            double log_probcpl = -0.5 * (zx + e).squaredNorm() + 0.5 * zx.squaredNorm();
            double log_ucpl = log(RNG::runif(RNG::rng));
            if (log_probcpl > 0 || log_ucpl < log_probcpl) // Coalesce proposals
            {
                same_proposal = true;
                
                yp = xp;
                yp_grad = xp_grad;
                yp_mean = xp_mean;
                logpi_yp = logpi_xp;
            }
            else // Reflect
            {
                e.normalize();
                Eigen::VectorXd zy = reflect(zx, e);

                yp = gaussian_proposal(y_mean, Sigma_chol_l, zy);
                logpi_yp = target.LogDensity(yp);
                yp_grad = target.GradLogDensity(yp);
                yp_mean = yp + 0.5 * precondition(Sigma, yp_grad);
            }
        }
        else // CRN
        {
            Eigen::VectorXd Lz = precondition_L(Sigma_chol_l, RNG::GaussianVector(d));

            xp = x_mean + Lz;
            logpi_xp = target.LogDensity(xp);
            xp_grad = target.GradLogDensity(xp);
            xp_mean = xp + 0.5 * precondition(Sigma, xp_grad);

            yp = y_mean + Lz;
            logpi_yp = target.LogDensity(yp);
            yp_grad = target.GradLogDensity(yp);
            yp_mean = yp + 0.5 * precondition(Sigma, yp_grad);
        }

        double logHR_x = logpi_xp - logpi_x + mala_proposal_logdensity_difference(x, x_mean, x_grad, xp, xp_mean, xp_grad);
        double logHR_y = logpi_yp - logpi_y + mala_proposal_logdensity_difference(y, y_mean, y_grad, yp, yp_mean, yp_grad);

        double log_u = log(RNG::runif(RNG::rng));
        if (logHR_x > 0 || log_u < logHR_x)
        {
            accepted_x_proposal = true;
            ++accepted_x;

            x = xp;
            logpi_x = logpi_xp;
            x_grad = xp_grad;
            x_mean = xp_mean;
        }
        if (logHR_y > 0 || log_u < logHR_y)
        {
            ++accepted_y;
            accepted_y_proposal = true;

            y = yp;
            logpi_y = logpi_yp;
            y_grad = yp_grad;
            y_mean = yp_mean;
        }

        // Break out of the loop if we've coalesced
        if (same_proposal && accepted_x_proposal && accepted_y_proposal)
        {
            tau = it;
            break;
        }
    }
}

