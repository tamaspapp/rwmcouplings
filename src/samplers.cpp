#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]

#include "samplers.h"

using Eigen::ArrayXd;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;
using Eigen::Ref;
using Rcpp::as;
using std::string;

/******************** Single-chain MCMC ********************/

//' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List rwm_cpp(const Rcpp::List &target_params,
                   const Rcpp::List &sampler_params,
                   const Eigen::Map<Eigen::VectorXd> &theta0,
                   const int &iter,
                   const int &thin)
{
    // Sampler parameters
    Map<MatrixXd> Sigma = Rcpp::as<Map<MatrixXd>>(sampler_params["Sigma"]);  // Proposal covariance. Can be fed in as a vector for diagonal covariances.

    // Output
    MatrixXd xs(iter / thin + 1, theta0.size());
    double acc_rate_x;

    string target_type = target_params["target_type"];
    if (target_type == "logistic_regression")
    {
        Map<MatrixXd> yX = as<Map<MatrixXd>>(target_params["yX"]);
        double lambda = target_params["lambda"];
        LogisticRegression pi(yX, lambda);

        rwm(pi, theta0, Sigma, iter, thin, xs, acc_rate_x);
    }
    else
    {
        std::cout << "Target type not implemented."
                  << "\n";
        assert(false);
    }

    return Rcpp::List::create(Rcpp::Named("xs") = xs,
                              Rcpp::Named("acc_rate_x") = acc_rate_x);
}

//' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List mala_cpp(const Rcpp::List &target_params,
                    const Rcpp::List &sampler_params,
                    const Eigen::Map<Eigen::VectorXd> &theta0,
                    const int &iter,
                    const int &thin)
{
    // Sampler parameters
    Map<MatrixXd> Sigma = Rcpp::as<Map<MatrixXd>>(sampler_params["Sigma"]); // Proposal covariance. Can be fed in as a vector for diagonal covariances.
    // Output
    MatrixXd xs(iter / thin + 1, theta0.size());
    double acc_rate_x;

    string target_type = target_params["target_type"];

    if (target_type == "logistic_regression")
    {
        Map<MatrixXd> yX = as<Map<MatrixXd>>(target_params["yX"]);
        double lambda = target_params["lambda"];
        LogisticRegression pi(yX, lambda);

        mala(pi, theta0, Sigma, iter, thin, xs, acc_rate_x);
    }
    else
    {
        std::cout << "Target type not implemented."
                  << "\n";
        assert(false);
    }

    return Rcpp::List::create(Rcpp::Named("xs") = xs,
                              Rcpp::Named("acc_rate_x") = acc_rate_x);
}

// /************** Coupled MCMC, same kernel ***************/

//' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List rwm_twoscalegcrn_cpp(const Rcpp::List &target_params,
                                const Rcpp::List &sampler_params, // Proposal covariance. Can be fed in as a vector for diagonal covariances.
                                const Eigen::Map<Eigen::VectorXd> &x0,
                                const Eigen::Map<Eigen::VectorXd> &y0,
                                const int &iter)
{
    // Sampler parameters
    Map<MatrixXd> Sigma = Rcpp::as<Map<MatrixXd>>(sampler_params["Sigma"]);
    double thresh = sampler_params["thresh"];

    // Output
    VectorXd fishy(x0.size());
    int tau = -1;
    int grad_evals = -1;

    string target_type = target_params["target_type"];
    if (target_type == "logistic_regression")
    {
        Map<MatrixXd> yX = as<Map<MatrixXd>>(target_params["yX"]);
        double lambda = target_params["lambda"];
        LogisticRegression pi(yX, lambda);

        rwm_twoscaleGCRN(pi, x0, y0, Sigma, iter, thresh, fishy, tau);
        grad_evals = pi.NumGradientEvaluations();
    }
    else
    {
        std::cout << "Target type not implemented."
                  << "\n";
        assert(false);
    }

    return Rcpp::List::create(Rcpp::Named("tau") = tau,
                              Rcpp::Named("fishy") = fishy,
                              Rcpp::Named("grad_evals") = grad_evals);
}


//' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List rwm_twoscalegcrefl_cpp(const Rcpp::List &target_params,
                                  const Rcpp::List &sampler_params, // Proposal covariance. Can be fed in as a vector for diagonal covariances.
                                  const Eigen::Map<Eigen::VectorXd> &x0,
                                  const Eigen::Map<Eigen::VectorXd> &y0,
                                  const int &iter)
{
    // Sampler parameters
    Map<MatrixXd> Sigma = Rcpp::as<Map<MatrixXd>>(sampler_params["Sigma"]);
    double thresh = sampler_params["thresh"];

    // Output
    VectorXd fishy(x0.size());
    int tau = -1;
    int grad_evals = -1;

    string target_type = target_params["target_type"];
    if (target_type == "logistic_regression")
    {
        Map<MatrixXd> yX = as<Map<MatrixXd>>(target_params["yX"]);
        double lambda = target_params["lambda"];
        LogisticRegression pi(yX, lambda);

        rwm_twoscaleGCRefl(pi, x0, y0, Sigma, iter, thresh, fishy, tau);
        grad_evals = pi.NumGradientEvaluations();
    }
    else
    {
        std::cout << "Target type not implemented."
                  << "\n";
        assert(false);
    }

    return Rcpp::List::create(Rcpp::Named("tau") = tau,
                              Rcpp::Named("fishy") = fishy,
                              Rcpp::Named("grad_evals") = grad_evals);
}

//' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List mala_twoscalecrn_cpp(const Rcpp::List &target_params,
                                const Rcpp::List &sampler_params,
                                const Eigen::Map<Eigen::VectorXd> &x0,
                                const Eigen::Map<Eigen::VectorXd> &y0,
                                const int &iter)
{
    // Sampler parameters
    Map<MatrixXd> Sigma = Rcpp::as<Map<MatrixXd>>(sampler_params["Sigma"]);
    double thresh = sampler_params["thresh"];

    // Output
    VectorXd fishy(x0.size());
    int tau = -1;

    string target_type = target_params["target_type"];
    if (target_type == "logistic_regression")
    {
        Map<MatrixXd> yX = as<Map<MatrixXd>>(target_params["yX"]);
        double lambda = target_params["lambda"];
        LogisticRegression pi(yX, lambda);

        mala_twoscaleCRN(pi, x0, y0, Sigma, iter, thresh, fishy, tau);
    }
    else
    {
        std::cout << "Target type not implemented."
                  << "\n";
        assert(false);
    }

    return Rcpp::List::create(Rcpp::Named("tau") = tau,
                              Rcpp::Named("fishy") = fishy);
}
