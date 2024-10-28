// Samplers for the stochastic volatility model

#include <RcppEigen.h>
#include "rng.h"
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

//' Sample the latent variables from the model
//'
//' @export
//[[Rcpp::export(rng = false)]]
Rcpp::NumericVector SampleLatentVariables(int T, double sig, double phi) {
  
    Rcpp::NumericVector x(T);
  
    x(0) = sqrt(sig * sig / (1. - phi * phi)) * RNG::rnorm(RNG::rng);

    for(int i = 1; i < T; i++)
    {
        x(i) = phi * x(i - 1) + sig * RNG::rnorm(RNG::rng);
    }
  
  return x;
}

// Target log-density
inline double logpi_nograd_(const Eigen::Ref<const Eigen::ArrayXd> &x,
                            const Eigen::Ref<const Eigen::ArrayXd> &y_sq,
                            const double &inv_beta_sq,
                            const double &phi, const double &phi_sq, const double &inv_phi_sq,
                            const double &inv_sigma_sq,
                            const int &T)
{
    // "temp" is (a constant shift off) from twice the negative log-posterior
    double temp = x.sum() + inv_beta_sq  * (y_sq * exp(-x)).sum()
                          + inv_sigma_sq * (phi * x.segment(0, T - 2) -  x.segment(1, T - 1)).square().sum()
                          + (1. - phi_sq) * inv_sigma_sq * x(0) * x(0);

    return (-0.5) * temp;
}

// Target log-density
inline double logpi_(const Eigen::Ref<const Eigen::ArrayXd> &x,
                    const Eigen::Ref<const Eigen::ArrayXd> &exp_minus_x,
                    const Eigen::Ref<const Eigen::ArrayXd> &x_diff, // Entries are: phi * x_{t} - x_{t+1}
                    const Eigen::Ref<const Eigen::ArrayXd> &y_sq,
                    const double &inv_beta_sq,
                    const double &phi, const double &phi_sq, const double &inv_phi_sq,
                    const double &inv_sigma_sq,
                    const int &T)
{
    // "temp" is (a constant off) from twice the negative log-posterior
    double temp = x.sum() + inv_beta_sq  * (y_sq * exp_minus_x).sum();

    temp += inv_sigma_sq * x_diff.square().sum();

    temp += (1 - phi_sq) * inv_sigma_sq * x(0) * x(0);

    return (-0.5) * temp;
}

// Target log-gradient
inline Eigen::ArrayXd gradlogpi_(const Eigen::Ref<const Eigen::ArrayXd> &x,
                         const Eigen::Ref<const Eigen::ArrayXd> &exp_minus_x,
                         const Eigen::Ref<const Eigen::ArrayXd> &x_diff, // Entries are: phi * x_{t} - x_{t+1}
                         const Eigen::Ref<const Eigen::ArrayXd> &y_sq,
                         const double &inv_beta_sq,
                         const double &phi, const double &phi_sq, const double &inv_phi_sq,
                         const double &inv_sigma_sq,
                         const int &T)
{
    // "temp" is the log-posterior-gradient
    Eigen::ArrayXd temp = 0.5 * inv_beta_sq * (y_sq * exp_minus_x);

    temp(0) -= (1 - phi_sq) * inv_sigma_sq * x(0);

    temp.segment(0, T - 2) -= phi * inv_sigma_sq * x_diff;

    temp.segment(1, T - 1) += inv_sigma_sq * x_diff;

    return (temp-0.5);
}

//' RWM algorithm
//'
//' @export
//[[Rcpp::export(rng = false)]]
Rcpp::List svmRWM(const Eigen::Map<Eigen::ArrayXd> x0,
                  const Eigen::Map<Eigen::ArrayXd> y,
                  const double &beta,
                  const double &sigma,
                  const double &phi,
                  const double &h,
                  const int &iter)
{
    // Declare constants
    double beta_sq     = beta * beta;
    double inv_beta_sq = 1 / beta_sq;

    double phi_sq     = phi * phi;
    double inv_phi_sq = 1 / phi_sq;

    double sigma_sq     = sigma * sigma;
    double inv_sigma_sq = 1 / sigma_sq;

    int T = x0.rows();
    Eigen::ArrayXd y_sq = y.square();

    // Declare target as lamda
    auto logpi = [&](Eigen::ArrayXd t) { return logpi_nograd_(t, y_sq, inv_beta_sq, phi, phi_sq, inv_phi_sq, inv_sigma_sq, T); };

    // Declare temporary variables
    double u;

    double logpi_x, logpi_xp, logHR_x;
    Eigen::ArrayXd x(T), xp(T);
    Eigen::ArrayXd x_dot(T);

    // Intialize storage
    int accepted = 0;

    x = x0;
    logpi_x = logpi(x);

    // "iter" iterations of MCMC
    for (int i = 1; i <= iter; i++)
    {
        // Proposal noise
        for (int j = 0; j < T; j++)
        {
            x_dot(j) = RNG::rnorm(RNG::rng);
        }
        xp = x + h * x_dot;

        // Compute log of Hastings ratio
        logpi_xp = logpi(xp);
        logHR_x = logpi_xp - logpi_x;

        // Perform acceptance step
        u = RNG::runif(RNG::rng);
        if (logHR_x > 0 || log(u) < logHR_x)
        {
            ++accepted;
            x = xp;
            logpi_x = logpi_xp;
        }
    }

    return Rcpp::List::create(Rcpp::Named("x") = x,
                              Rcpp::Named("acc_rate") = (double) accepted / (double) iter);
}


//' RWM algorithms, with CRN coupling
//'
//' @export
//[[Rcpp::export(rng = false)]]
Rcpp::List svmCRNRWM(const Eigen::Map<Eigen::ArrayXd> x0,
                     const Eigen::Map<Eigen::ArrayXd> y0,
                     const Eigen::Map<Eigen::ArrayXd> y_data,
                     const double &beta,
                     const double &sigma,
                     const double &phi,
                     const double &h,
                     const int &iter,
                     const int &thin)
{
    // Declare constants
    double beta_sq     = beta * beta;
    double inv_beta_sq = 1 / beta_sq;

    double phi_sq     = phi * phi;
    double inv_phi_sq = 1 / phi_sq;

    double sigma_sq     = sigma * sigma;
    double inv_sigma_sq = 1 / sigma_sq;

    int T = x0.rows();
    Eigen::ArrayXd y_data_sq = y_data.square();

    // Declare target as lamda
    auto logpi = [&](Eigen::ArrayXd t) { return logpi_nograd_(t, y_data_sq, inv_beta_sq, phi, phi_sq, inv_phi_sq, inv_sigma_sq, T); };

    // Declare temporary variables
    double log_u;

    double logpi_x, logpi_xp, logHR_x;
    Eigen::ArrayXd x(T), xp(T), x_dot(T);

    double logpi_y, logpi_yp, logHR_y;
    Eigen::ArrayXd y(T), yp(T);
    
    // Intialize storage
    std::vector<double> xy_square_dist(iter/thin + 1, 0.);
    
    // Start off MCMC
    x = x0;
    logpi_x = logpi(x);
    y = y0;
    logpi_y = logpi(y);

    xy_square_dist[0] = (x - y).square().sum();

    // Do "iter" iterations of MCMC
    for (int i = 1; i <= iter; i++)
    {
        // Sample proposal noise
        for (int j = 0; j < T; j++)
        {
            x_dot(j) = RNG::rnorm(RNG::rng);
        }
        xp = x + h * x_dot;
        // CRN coupling for Y
        yp = y + h * x_dot;
        
        // Compute log of Hastings ratio
        logpi_xp = logpi(xp);
        logpi_yp = logpi(yp);
        logHR_x = logpi_xp - logpi_x;
        logHR_y = logpi_yp - logpi_y;

        // Accept-reject
        log_u = log(RNG::runif(RNG::rng));
        if (logHR_x > 0 || log_u < logHR_x)
        {
            x = xp;
            logpi_x = logpi_xp;
        }
        if (logHR_y > 0 || log_u < logHR_y)
        {
            y = yp;
            logpi_y = logpi_yp;
        }

        if (i % thin == 0)
        {
            xy_square_dist[i/thin] = (x - y).square().sum();
        }
    }

    return Rcpp::List::create(Rcpp::Named("squaredist") = Rcpp::wrap(xy_square_dist));
}

//' RWM algorithms, with reflection-maximal coupling
//'
//' @export
//[[Rcpp::export(rng = false)]]
Rcpp::List svmReflMaxRWM(const Eigen::Map<Eigen::ArrayXd> x0,
                         const Eigen::Map<Eigen::ArrayXd> y0,
                         const Eigen::Map<Eigen::ArrayXd> y_data,
                         const double &beta,
                         const double &sigma,
                         const double &phi,
                         const double &h,
                         const int &iter,
                         const int &thin)
{
    // Declare constants
    double beta_sq     = beta * beta;
    double inv_beta_sq = 1 / beta_sq;

    double phi_sq     = phi * phi;
    double inv_phi_sq = 1 / phi_sq;

    double sigma_sq     = sigma * sigma;
    double inv_sigma_sq = 1 / sigma_sq;

    int T = x0.rows();
    Eigen::ArrayXd y_data_sq = y_data.square();

    // Declare target as lambda
    auto logpi = [&](Eigen::ArrayXd t) { return logpi_nograd_(t, y_data_sq, inv_beta_sq, phi, phi_sq, inv_phi_sq, inv_sigma_sq, T); };

    // Declare temporary variables
    double log_u, logcpl, ucpl, z_sqnorm;
    bool coupled = false;
    Eigen::ArrayXd z(T);

    double logpi_x, logpi_xp;
    double logHR_x;
    bool acc_x;

    Eigen::ArrayXd x(T), xp(T);
    Eigen::ArrayXd x_dot(T);

    double logpi_y, logpi_yp;
    double logHR_y;
    bool acc_y;

    Eigen::ArrayXd y(T), yp(T);
    Eigen::ArrayXd y_dot(T);
    
    // Intialize storage
    int tau = -1;
    std::vector<double> xy_square_dist(iter/thin + 1, 0.);

    x = x0;
    logpi_x = logpi(x);
    y = y0;
    logpi_y = logpi(y);

    xy_square_dist[0] = (x-y).square().sum();

    // Iterations until coupling or max iteration count
    for (int i = 1; i <= iter; i++)
    {
        // Proposal noise
        for (int j = 0; j < T; j++)
        {
            x_dot(j) = RNG::rnorm(RNG::rng);
        }
        xp = x + h * x_dot;
        logpi_xp = logpi(xp);

        // Reflection-maximal coupling for Y
        z = (x - y) / h;
        z_sqnorm = z.square().sum();

        logcpl = -0.5 * (x_dot + z).square().sum() + 0.5 * x_dot.square().sum();
        ucpl = RNG::runif(RNG::rng);

        coupled = false;
        if (logcpl > 0 || log(ucpl) < logcpl) // Maximal coupling of proposals
        {
            coupled = true;
            yp = xp;
            logpi_yp = logpi_xp;
        }
        else // Reflect
        {
            yp = y + h * (x_dot - 2 / z_sqnorm * (z * x_dot).sum() * z);
            logpi_yp = logpi(yp);
        }

        // Compute log of Hastings ratio
        logHR_x = logpi_xp - logpi_x;
        logHR_y = logpi_yp - logpi_y;
        
        // Accept-reject
        log_u = log(RNG::runif(RNG::rng));
        if (logHR_x > 0 || log_u < logHR_x)
        {
            x = xp;
            logpi_x = logpi_xp;
            acc_x = true;
        } else { acc_x = false; }        
        if (logHR_y > 0 || log_u < logHR_y)
        {
            y = yp;
            logpi_y = logpi_yp;
            acc_y = true;
        } else { acc_y = false; }

        if (i % thin == 0)
        {
            xy_square_dist[i/thin] = z_sqnorm*h*h;
        }

        // Stop early if coupled
        if (acc_x * acc_y * coupled == true)
        {
            tau = i;
            break;
        }
    }

    return Rcpp::List::create(Rcpp::Named("squaredist") = Rcpp::wrap(xy_square_dist),
                              Rcpp::Named("tau") = tau);
}

//' RWM algorithms, with two-scale coupling using GCRN when (|x-y|^2 >= thresh) and using reflection-maximal coupling when (|x-y|^2 < thresh)
//'
//' @export
//[[Rcpp::export(rng = false)]]
Rcpp::List svmTwoScaleGCRNRWM(const Eigen::Map<Eigen::ArrayXd> x0,
                              const Eigen::Map<Eigen::ArrayXd> y0,
                              const Eigen::Map<Eigen::ArrayXd> y_data,
                              const double &beta,
                              const double &sigma,
                              const double &phi,
                              const double &h,
                              const int &iter,
                              const int &thin,
                              const double &thresh)
{
    // Declare constants
    double h_sq     = h * h;
    double inv_h_sq = 1 / h_sq;
    
    double beta_sq     = beta * beta;
    double inv_beta_sq = 1 / beta_sq;

    double phi_sq     = phi * phi;
    double inv_phi_sq = 1 / phi_sq;

    double sigma_sq     = sigma * sigma;
    double inv_sigma_sq = 1 / sigma_sq;

    int T = x0.rows();
    Eigen::ArrayXd y_data_sq = y_data.square();

    // Declare target and gradient as lambda
    auto logpi     = [&](Eigen::ArrayXd t1, Eigen::ArrayXd t2, Eigen::ArrayXd t3) { return logpi_(t1, t2, t3, y_data_sq, inv_beta_sq, phi, phi_sq, inv_phi_sq, inv_sigma_sq, T); };
    auto gradlogpi = [&](Eigen::ArrayXd t1, Eigen::ArrayXd t2, Eigen::ArrayXd t3) { return gradlogpi_(t1, t2, t3, y_data_sq, inv_beta_sq, phi, phi_sq, inv_phi_sq, inv_sigma_sq, T); };
    auto logpi_nograd = [&](Eigen::ArrayXd t) { return logpi_nograd_(t, y_data_sq, inv_beta_sq, phi, phi_sq, inv_phi_sq, inv_sigma_sq, T); };

    // Declare temporary variables
    double log_u, ucpl, logcpl, squaredist, z1, z_sqnorm;
    Eigen::ArrayXd z(T);

        // X-chain
    double logpi_x, logpi_xp;
    double logHR_x;
    bool acc_x = true;

    Eigen::ArrayXd x(T), xp(T);
    Eigen::ArrayXd x_dot(T);
    Eigen::ArrayXd n_x(T);
    Eigen::ArrayXd exp_minus_x(T), exp_minus_xp(T);
    Eigen::ArrayXd x_diff(T - 1), xp_diff(T - 1);

        // Y-chain
    double logpi_y, logpi_yp;
    double logHR_y;
    bool acc_y = true;

    Eigen::ArrayXd y(T), yp(T);
    Eigen::ArrayXd y_dot(T);
    Eigen::ArrayXd n_y(T);
    Eigen::ArrayXd exp_minus_y(T), exp_minus_yp(T);
    Eigen::ArrayXd y_diff(T - 1), yp_diff(T - 1);
    
    // Intialize storage
    int tau = -1;
    std::vector<double> xy_square_dist(iter/thin + 1, 0.);

    // Start MCMC
    x = x0;
    exp_minus_x = (-x).exp();
    x_diff = phi * x.segment(0, T - 2) -  x.segment(1, T - 1);
    logpi_x = logpi(x, exp_minus_x, x_diff);
    n_x = gradlogpi(x, exp_minus_x, x_diff);
    n_x /= sqrt(n_x.square().sum());  // Normalize

    y = y0;
    exp_minus_y = (-y).exp();
    y_diff = phi * y.segment(0, T - 2) -  y.segment(1, T - 1);
    logpi_y = logpi(y, exp_minus_y, y_diff);
    n_y = gradlogpi(y, exp_minus_y, y_diff);
    n_y /= sqrt(n_y.square().sum());

    squaredist = (x-y).square().sum();
    xy_square_dist[0] = squaredist;


    bool coupled = false, refl = false;
    // Iterations until coupling or max iteration count
    for (int i = 1; i <= iter; i++)
    {
        if(squaredist < thresh) // Do reflection-maximal coupling
        {
            refl = true;

            // Proposal noises
            x_dot = Eigen::ArrayXd::NullaryExpr(T,[&](){ return RNG::rnorm(RNG::rng); }); 
            xp = x + h * x_dot;
            logpi_xp = logpi_nograd(xp);

            // Reflection-maximal coupling for Y
            z = (x - y) / h;
            z_sqnorm = z.square().sum();

            logcpl = -0.5 * (x_dot + z).square().sum() + 0.5 * x_dot.square().sum();
            ucpl = RNG::runif(RNG::rng);
            if (logcpl > 0 || log(ucpl) < logcpl) // Maximal coupling of proposals
            {
                coupled = true;
                yp = xp;
                logpi_yp = logpi_xp;
            }
            else // Reflect
            {
                coupled = false;
                yp = y + h * (x_dot - 2. / z_sqnorm * (z * x_dot).sum() * z);
                logpi_yp = logpi_nograd(yp);
            }

                // Compute log of Hastings ratio
            logHR_x = logpi_xp - logpi_x;
            logHR_y = logpi_yp - logpi_y;

            // Accept-reject
            log_u = log(RNG::runif(RNG::rng));
            if (logHR_x > 0 || log_u < logHR_x)
            {
                x = xp;
                logpi_x = logpi_xp;
                acc_x = true;
            } else { acc_x = false; }
            if (logHR_y > 0 || log_u < logHR_y)
            {
                y = yp;
                logpi_y = logpi_yp;
                acc_y = true;
            } else { acc_y = false; }

            squaredist = z_sqnorm*h*h;
            if (i % thin == 0)
            {
                xy_square_dist[i/thin] = squaredist;
            }
            
            // Stop early if coupled
            if (acc_x * acc_y * coupled == true)
            {
                tau = i;
                break;
            } 
        }
        else 
        {
            if(refl == true) // If did a reflection coupling previously, refresh the gradients
            {
                refl = false;

                exp_minus_x = (-x).exp();
                x_diff = phi * x.segment(0, T - 2) -  x.segment(1, T - 1);
                exp_minus_y = (-y).exp();
                y_diff = phi * y.segment(0, T - 2) -  y.segment(1, T - 1);

                n_x = gradlogpi(x, exp_minus_x, x_diff); n_x /= sqrt(n_x.square().sum());
                n_y = gradlogpi(y, exp_minus_y, y_diff); n_y /= sqrt(n_y.square().sum());
                
            } 
            else // Refresh the gradients only if needed
            {
                if(acc_x == true)
                {
                    n_x = gradlogpi(x, exp_minus_x, x_diff);
                    n_x /= sqrt(n_x.square().sum());
                }
                if(acc_y == true)
                {
                    n_y = gradlogpi(y, exp_minus_y, y_diff);
                    n_y /= sqrt(n_y.square().sum());
                }
            }
            
            // Proposal noises
            double z1 = RNG::rnorm(RNG::rng);
            z = Eigen::ArrayXd::NullaryExpr(T,[&](){ return RNG::rnorm(RNG::rng); }); 
            x_dot = z + (z1 - (z * n_x).sum()) * n_x;
            y_dot = z + (z1 - (z * n_y).sum()) * n_y;
            
            xp = x + h * x_dot;
            exp_minus_xp = (-xp).exp();
            xp_diff = phi * xp.segment(0, T - 2) -  xp.segment(1, T - 1);
            logpi_xp = logpi(xp, exp_minus_xp, xp_diff);

            yp = y + h * y_dot;
            exp_minus_yp = (-yp).exp();
            yp_diff = phi * yp.segment(0, T - 2) -  yp.segment(1, T - 1);
            logpi_yp = logpi(yp, exp_minus_yp, yp_diff);

            // Compute log of Hastings ratio
            logHR_x = logpi_xp - logpi_x;
            logHR_y = logpi_yp - logpi_y;

            // Accept-reject
            log_u = log(RNG::runif(RNG::rng));
            if (logHR_x > 0 || log_u < logHR_x)
            {
                x = xp;
                logpi_x = logpi_xp;
                exp_minus_x = exp_minus_xp;
                x_diff = xp_diff;
                acc_x = true;
            } else { acc_x = false; }
            if (logHR_y > 0 || log_u < logHR_y)
            {
                y = yp;
                logpi_y = logpi_yp;
                exp_minus_y = exp_minus_yp;
                y_diff = yp_diff;
                acc_y = true;
            } else { acc_y = false; }

            squaredist = (x-y).square().sum();
            if (i % thin == 0)
            {
                xy_square_dist[i/thin] = squaredist;
            }
        }
    }
    
    return Rcpp::List::create(Rcpp::Named("squaredist") = Rcpp::wrap(xy_square_dist),
                              Rcpp::Named("tau") = tau);
}






/* Two-scale coupling, with corrected reflection and reflection. */

// log-density
inline double logpi_svm(const Eigen::Ref<const Eigen::VectorXd> &x,
                        const Eigen::Ref<const Eigen::VectorXd> &y_sq,
                        const double &beta,
                        const double &sigma,
                        const double &phi)
{
    double temp = x.sum() + y_sq.dot((-x).array().exp().matrix()) / (beta*beta)
                          + (phi * x.head(x.size()-1) - x.tail(x.size()-1)).squaredNorm() / (sigma*sigma)
                          + (1. - phi*phi) * x(0) * x(0) / (sigma*sigma);
    return (-0.5) * temp;
}
// score
inline Eigen::VectorXd gradlogpi_svm(const Eigen::Ref<const Eigen::VectorXd> &x,
                                     const Eigen::Ref<const Eigen::VectorXd> &y_sq,
                                     const double &beta,
                                     const double &sigma,
                                     const double &phi)
{
    Eigen::VectorXd temp = Eigen::VectorXd::Constant(x.size(),-0.5); 
    temp += 0.5 * y_sq.cwiseProduct((-x).array().exp().matrix()) / (beta*beta);
    temp(0) -= (1 - phi*phi) * x(0) / (sigma*sigma);
    temp.head(temp.size()-1) -= phi * (phi * x.head(x.size()-1) - x.tail(x.size()-1)) / (sigma*sigma);
    temp.tail(temp.size()-1) += (phi * x.head(x.size()-1) - x.tail(x.size()-1)) / (sigma*sigma);
    return temp;
}

//' Two-scale coupling, with corrected reflection and reflection.
//' @export
//[[Rcpp::export(rng = false)]]
Rcpp::List svmTwoScaleGCReflRWM(const Eigen::Map<Eigen::VectorXd> x0,
                                const Eigen::Map<Eigen::VectorXd> y0,
                                const Eigen::Map<Eigen::VectorXd> y_data,
                                const double &beta,
                                const double &sigma,
                                const double &phi,
                                const double &h,
                                const int &iter,
                                const int &thin,
                                const double &thresh) // Note that this is used in the check: |x-y|^2 < thresh
{   
    // Declare logdensity and score as lambdas
    Eigen::VectorXd y_data_sq = y_data.cwiseProduct(y_data);
    auto logpi = [&](Eigen::VectorXd x) { return logpi_svm(x, y_data_sq, beta, sigma, phi); };
    auto gradlogpi = [&](Eigen::VectorXd x) { return  gradlogpi_svm(x, y_data_sq, beta, sigma, phi); };

    // Declare constants
    int d = x0.rows();
    auto inv_h = 1./h;

    // Declare temporary variables
    Eigen::VectorXd x(d), xp(d), e_gx(d);
    Eigen::VectorXd y(d), yp(d), e_gy(d);
    Eigen::VectorXd x_dot(d);
    Eigen::VectorXd y_dot(d);
    Eigen::VectorXd z(d), e(d), c_x(d), c_y(d);
    double z1;
    double u;
    double e_sqnorm;
    double logpi_x, logpi_xp, logHR_x;
    double logpi_y, logpi_yp, logHR_y;
    int accepted_x = 0, accepted_y = 0;
    bool acc_x_flag = true, acc_y_flag = true, cpl_flag = false;
    double log_refl, u_refl;
    int tau = -1, iter_done = iter;

    // Do MCMC
    x = x0;
    logpi_x = logpi(x);
    e_gx = gradlogpi(x).normalized();

    y = y0;
    logpi_y = logpi(y);
    e_gy = gradlogpi(y).normalized();

    Eigen::VectorXd xy_square_dist = Eigen::VectorXd::Zero(iter/thin + 1);
    double squaredist = (x - y).squaredNorm();
    xy_square_dist(0) = squaredist;

    // "iter" iterations of MCMC
    for (int i = 1; i <= iter; i++)
    {
        cpl_flag = false;
        e = (x - y) * inv_h;
        if (squaredist < thresh) // Reflection-maximal
        {
            for (int j = 0; j < d; j++)
            {
              x_dot(j) = RNG::rnorm(RNG::rng);
            }
            xp = x + x_dot * h;
            logpi_xp = logpi(xp);

            log_refl = -0.5 * (x_dot + e).squaredNorm() + 0.5 * x_dot.squaredNorm();
            u_refl = RNG::runif(RNG::rng);
            if (log_refl > 0 || log(u_refl) < log_refl) // Couple the proposals
            {
                cpl_flag = true;
                yp = xp;
                logpi_yp = logpi_xp;
            }
            else // Reflect
            {
                e.normalize();
                y_dot = x_dot - (2 * x_dot.dot(e)) * e;
                yp = y + y_dot * h;
                logpi_yp = logpi(yp);
            }
        }
        else // Corrected reflection
        {
            e.normalize();
            c_x = (e_gx - e_gx.dot(e) * e).normalized();
            c_y = (e_gy - e_gy.dot(e) * e).normalized();

            // Proposals
            for (int j = 0; j < d; j++)
            {
                z(j) = RNG::rnorm(RNG::rng);
            }
            z1 = RNG::rnorm(RNG::rng);

            x_dot = z - z.dot(c_x) * c_x + z1 * c_x;
            y_dot = z - 2 * z.dot(e) * e - z.dot(c_y) * c_y + z1 * c_y;

            xp = x + x_dot * h;
            logpi_xp = logpi(xp);

            yp = y + y_dot * h;
            logpi_yp = logpi(yp);
        }

        // Acceptances
        logHR_x = logpi_xp - logpi_x;
        logHR_y = logpi_yp - logpi_y;

        u = RNG::runif(RNG::rng);
        if (logHR_x > 0 || log(u) < logHR_x)
        {
            acc_x_flag = true;
            ++accepted_x;
            x = xp;
            logpi_x = logpi_xp;
            e_gx = gradlogpi(x).normalized(); // Update gradient if accepted; could update gradients less often with more complicated control flow.
        }
        else{acc_x_flag = false;}

        if (logHR_y > 0 || log(u) < logHR_y)
        {
            acc_y_flag = true;
            ++accepted_y;
            y = yp;
            logpi_y = logpi_yp;
            e_gy = gradlogpi(y).normalized();
        }
        else{acc_y_flag = false;}

        if (cpl_flag && acc_x_flag && acc_y_flag) // Coalescence!
        {
            tau = i;
            iter_done = tau;
            break;
        }

        squaredist = (x - y).squaredNorm();
        if (i % thin == 0)
        {
            xy_square_dist(i/thin) = squaredist;
        }
    }
    return Rcpp::List::create(Rcpp::Named("x") = x,
                              Rcpp::Named("y") = y,
                              Rcpp::Named("squaredist") = xy_square_dist,
                              Rcpp::Named("tau") = tau,
                              Rcpp::Named("acc_x") = (double)accepted_x / (double) iter_done,
                              Rcpp::Named("acc_y") = (double)accepted_y / (double) iter_done);
}
