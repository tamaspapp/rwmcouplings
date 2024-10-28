#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include "targets.h"

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Rcpp::as;
using std::cout;
using std::string;

/************* Target evaluation **********/

//' @export
// [[Rcpp::export(rng = false)]]
double logpi_cpp(const Rcpp::List &target_params, 
                 const Eigen::Map<Eigen::VectorXd> &theta)
{
    string target_type = target_params["target_type"];

    if (target_type == "logistic_regression")
    {
        Map<MatrixXd> yX = as<Map<MatrixXd>>(target_params["yX"]);
        double lambda = target_params["lambda"];
        LogisticRegression pi(yX, lambda);

        return pi.LogDensity(theta);
    }
    else
    {
        cout << "Target type not implemented."
             << "\n";
        return 0;
    }
}

//' @export
// [[Rcpp::export(rng = false)]]
Eigen::VectorXd gradlogpi_cpp(const Rcpp::List &target_params, 
                              const Eigen::Map<Eigen::VectorXd> &theta)
{
    string target_type = target_params["target_type"];

    if (target_type == "logistic_regression")
    {
        Map<MatrixXd> yX = as<Map<MatrixXd>>(target_params["yX"]);
        double lambda = target_params["lambda"];
        LogisticRegression pi(yX, lambda);

        return pi.GradLogDensity(theta);
    }
    else
    {
        cout << "Target type not implemented."
             << "\n";
        return Eigen::VectorXd::Zero(theta.size());
    }
}

//' @export
// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd hesslogpi_cpp(const Rcpp::List &target_params, 
                              const Eigen::Map<Eigen::VectorXd> &theta)
{
    string target_type = target_params["target_type"];

    if (target_type == "logistic_regression")
    {
        Map<MatrixXd> yX = as<Map<MatrixXd>>(target_params["yX"]);
        double lambda = target_params["lambda"];
        LogisticRegression pi(yX, lambda);

        return pi.HessLogDensity(theta);
    }
    else
    {
        cout << "Target type not implemented."
             << "\n";
        return Eigen::MatrixXd::Zero(theta.size(), theta.size());
    }
}
