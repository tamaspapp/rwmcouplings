// Stochastic volatility model, bias of Laplace approximation

#include <RcppEigen.h>
#include "rng.h"
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]


// Targets and gradients

// SVM log-density
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

// SVM log-density gradient
inline Eigen::VectorXd gradlogpi_svm(const Eigen::Ref<const Eigen::VectorXd> &x,
                                     const Eigen::Ref<const Eigen::VectorXd> &y_sq,
                                     const double &beta,
                                     const double &sigma,
                                     const double &phi)
{
    Eigen::VectorXd temp = Eigen::VectorXd::Constant(x.size(),-0.5); 
    temp += 0.5 * y_sq.cwiseProduct((-x).array().exp().matrix())  / (beta*beta);
    temp(0) -= (1 - phi*phi) * x(0) / (sigma*sigma);
    temp.head(temp.size()-1) -= phi * (phi * x.head(x.size()-1) - x.tail(x.size()-1)) / (sigma*sigma);
    temp.tail(temp.size()-1) += (phi * x.head(x.size()-1) - x.tail(x.size()-1)) / (sigma*sigma);
    return temp;
}

// Gaussian log-density
inline double logpi_gauss(const Eigen::Ref<const Eigen::VectorXd> &x,
                          const Eigen::Ref<const Eigen::VectorXd> &mu,
                          const Eigen::Ref<const Eigen::SparseMatrix<double>> &Omega_chol) // Upper triangular part of Cholesky decomposition of precision matrix Omega
{

    return (-0.5) * (Omega_chol * (x - mu)).squaredNorm();
}

inline Eigen::VectorXd gradlogpi_gauss(const Eigen::Ref<const Eigen::VectorXd> &x,
                                       const Eigen::Ref<const Eigen::VectorXd> &mu,
                                       const Eigen::Ref<const Eigen::SparseMatrix<double>> &Omega)
{
    return Omega * (mu - x);
}


//' GCRN-coupled RWM algorithms targeting SVM and a Laplace approximation
//'
//' @export
//[[Rcpp::export(rng = false)]]
Rcpp::List svmLaplaceGCRNRWM(const Eigen::Map<Eigen::VectorXd> x0, // SVM
                             const Eigen::Map<Eigen::VectorXd> y0, // Laplace approximation
                             const Eigen::Map<Eigen::ArrayXd> y_data,
                             const double &beta,
                             const double &sigma,
                             const double &phi,
                             const Eigen::Map<Eigen::VectorXd> mu,                     // Mean for the Laplace approximation
                             const Eigen::Map<Eigen::SparseMatrix<double>> Omega,      // Precision matrix of the Laplace approximation
                             const Eigen::Map<Eigen::SparseMatrix<double>> Omega_chol, // Cholesky factor (upper triangular) for the precision matrix of the Laplace approximation
                             const double &h,
                             const int &iter,
                             const int &thin)
{
    // Declare constants
    Eigen::VectorXd y_data_sq = y_data.square();
    int T = y_data.size();

    // Declare targets andgradients as lambdas
    auto logpi     = [&](Eigen::VectorXd t) { return logpi_svm(t, y_data_sq, beta, sigma, phi); };
    auto gradlogpi = [&](Eigen::VectorXd t) { return gradlogpi_svm(t, y_data_sq, beta, sigma, phi); };
    auto logpi_laplace     = [&](Eigen::VectorXd t) { return logpi_gauss(t,mu,Omega_chol); };
    auto gradlogpi_laplace = [&](Eigen::VectorXd t) { return gradlogpi_gauss(t,mu,Omega); };

    // Declare temporary variables
    double log_u, z1;
    Eigen::VectorXd z(T);

        // X-chain
    double logpi_x, logpi_xp;
    double logHR_x;

    Eigen::VectorXd x(T), xp(T);
    Eigen::VectorXd x_dot(T);
    Eigen::VectorXd n_x(T);

        // Y-chain
    double logpi_y, logpi_yp;
    double logHR_y;

    Eigen::VectorXd y(T), yp(T);
    Eigen::VectorXd y_dot(T);
    Eigen::VectorXd n_y(T);
    
    // Intialize storage
    std::vector<double> xy_square_dist(iter/thin + 1);

    // Start MCMC
    x = x0;
    logpi_x = logpi(x);
    n_x = gradlogpi(x).normalized();

    y = y0;
    logpi_y = logpi_laplace(y);
    n_y = gradlogpi_laplace(y).normalized();

    xy_square_dist[0] = (x-y).squaredNorm();

    // "iter" iterations
    for (int i = 1; i <= iter; i++)
    {
        // Proposals
        double z1 = RNG::rnorm(RNG::rng);
        for (int j = 0; j < T; j++)
        {
            z(j) = RNG::rnorm(RNG::rng);
        }
        x_dot = z + (z1 - z.dot(n_x)) * n_x;
        y_dot = z + (z1 - z.dot(n_y)) * n_y;
            
        xp = x + h * x_dot;
        logpi_xp = logpi(xp);
        yp = y + h * y_dot;
        logpi_yp = logpi_laplace(yp);

        // Hastings ratio
        logHR_x = logpi_xp - logpi_x;
        logHR_y = logpi_yp - logpi_y;

        // Accept-reject
        log_u = log(RNG::runif(RNG::rng));
        if (logHR_x > 0 || log_u < logHR_x)
        {
            x = xp;
            logpi_x = logpi_xp;
            n_x = gradlogpi(x).normalized();
        }
        if (logHR_y > 0 || log_u < logHR_y)
        {
            y = yp;
            logpi_y = logpi_yp;
            n_y = gradlogpi_laplace(y).normalized();
        }

        // Storage
        if (i % thin == 0)
        {
            xy_square_dist[i/thin] = (x-y).squaredNorm();
        }
    }
    
    return Rcpp::List::create(Rcpp::Named("x") = Rcpp::wrap(x),
                              Rcpp::Named("y") = Rcpp::wrap(y),
                              Rcpp::Named("squaredist") = Rcpp::wrap(xy_square_dist));
}


//' CRN-coupled RWM algorithms targeting SVM and a Laplace approximation
//'
//' @export
//[[Rcpp::export(rng = false)]]
Rcpp::List svmLaplaceCRNRWM(const Eigen::Map<Eigen::VectorXd> x0, // SVM
                            const Eigen::Map<Eigen::VectorXd> y0, // Laplace approximation
                            const Eigen::Map<Eigen::ArrayXd> y_data,
                            const double &beta,
                            const double &sigma,
                            const double &phi,
                            const Eigen::Map<Eigen::VectorXd> mu,                     // Mean for the Laplace approximation
                            const Eigen::Map<Eigen::SparseMatrix<double>> Omega,      // Precision matrix of the Laplace approximation
                            const Eigen::Map<Eigen::SparseMatrix<double>> Omega_chol, // Cholesky factor (upper triangular) for the precision matrix of the Laplace approximation
                            const double &h,
                            const int &iter,
                            const int &thin)
{
    // Declare constants
    Eigen::VectorXd y_data_sq = y_data.square();
    int T = y_data.size();

    // Declare targets andgradients as lambdas
    auto logpi     = [&](Eigen::VectorXd t) { return logpi_svm(t, y_data_sq, beta, sigma, phi); };
    auto logpi_laplace     = [&](Eigen::VectorXd t) { return logpi_gauss(t,mu,Omega_chol); };

    // Declare temporary variables
    double log_u;
    Eigen::VectorXd z(T);

        // X-chain
    double logpi_x, logpi_xp;
    double logHR_x;

    Eigen::VectorXd x(T), xp(T);

        // Y-chain
    double logpi_y, logpi_yp;
    double logHR_y;

    Eigen::VectorXd y(T), yp(T);
    
    // Intialize storage
    std::vector<double> xy_square_dist(iter/thin + 1);

    // Start MCMC
    x = x0;
    logpi_x = logpi(x);

    y = y0;
    logpi_y = logpi_laplace(y);

    xy_square_dist[0] = (x-y).squaredNorm();

    // "iter" iterations
    for (int i = 1; i <= iter; i++)
    {
        // Proposals
        for (int j = 0; j < T; j++)
        {
            z(j) = RNG::rnorm(RNG::rng);
        }
            
        xp = x + h * z;
        logpi_xp = logpi(xp);
        yp = y + h * z;
        logpi_yp = logpi_laplace(yp);

        // Hastings ratio
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

        // Storage
        if (i % thin == 0)
        {
            xy_square_dist[i/thin] = (x-y).squaredNorm();
        }
    }
    
    return Rcpp::List::create(Rcpp::Named("x") = Rcpp::wrap(x),
                              Rcpp::Named("y") = Rcpp::wrap(y),
                              Rcpp::Named("squaredist") = Rcpp::wrap(xy_square_dist));
}


//' Reflection-coupled RWM algorithms targeting SVM and a Laplace approximation
//'
//' @export
//[[Rcpp::export(rng = false)]]
Rcpp::List svmLaplaceReflRWM(const Eigen::Map<Eigen::VectorXd> x0, // SVM
                            const Eigen::Map<Eigen::VectorXd> y0, // Laplace approximation
                            const Eigen::Map<Eigen::ArrayXd> y_data,
                            const double &beta,
                            const double &sigma,
                            const double &phi,
                            const Eigen::Map<Eigen::VectorXd> mu,                     // Mean for the Laplace approximation
                            const Eigen::Map<Eigen::SparseMatrix<double>> Omega,      // Precision matrix of the Laplace approximation
                            const Eigen::Map<Eigen::SparseMatrix<double>> Omega_chol, // Cholesky factor (upper triangular) for the precision matrix of the Laplace approximation
                            const double &h,
                            const int &iter,
                            const int &thin)
{
    // Declare constants
    Eigen::VectorXd y_data_sq = y_data.square();
    int T = y_data.size();

    // Declare targets andgradients as lambdas
    auto logpi     = [&](Eigen::VectorXd t) { return logpi_svm(t, y_data_sq, beta, sigma, phi); };
    auto logpi_laplace     = [&](Eigen::VectorXd t) { return logpi_gauss(t,mu,Omega_chol); };

    // Declare temporary variables
    double log_u, squaredist;
    Eigen::VectorXd e(T);

        // X-chain
    double logpi_x, logpi_xp;
    double logHR_x;

    Eigen::VectorXd x(T), xp(T), xdot(T);

        // Y-chain
    double logpi_y, logpi_yp;
    double logHR_y;

    Eigen::VectorXd y(T), yp(T), ydot(T);
    
    // Intialize storage
    std::vector<double> xy_square_dist(iter/thin + 1);

    // Start MCMC
    x = x0;
    logpi_x = logpi(x);

    y = y0;
    logpi_y = logpi_laplace(y);

    e = x - y;
    squaredist = e.squaredNorm();
    xy_square_dist[0] = squaredist;

    // "iter" iterations
    for (int i = 1; i <= iter; i++)
    {
        // Proposals
        for (int j = 0; j < T; j++)
        {
            xdot(j) = RNG::rnorm(RNG::rng);
        }
            
        xp = x + h * xdot;
        logpi_xp = logpi(xp);
        yp = y + h * (xdot - (2. / squaredist) * e.dot(xdot) * e);
        logpi_yp = logpi_laplace(yp);

        // Hastings ratio
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

        e = x - y;
        squaredist = e.squaredNorm();
        // Storage
        if (i % thin == 0)
        {
            xy_square_dist[i/thin] = squaredist;
        }
    }
    
    return Rcpp::List::create(Rcpp::Named("x") = Rcpp::wrap(x),
                              Rcpp::Named("y") = Rcpp::wrap(y),
                              Rcpp::Named("squaredist") = Rcpp::wrap(xy_square_dist));
}
