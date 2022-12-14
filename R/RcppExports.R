# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Hug and Hop algorithm, CRN for Hug and GCRN for Hop
#' 
#' Standard Gaussian momentum
#' 
#' @export
GCRNHugHop <- function(x0, y0, hug_params, hop_params, iter, thin) {
    .Call('_rwmcouplings_GCRNHugHop', PACKAGE = 'rwmcouplings', x0, y0, hug_params, hop_params, iter, thin)
}

#' @export
SetSeed_cpp <- function(seed, stream = 0L) {
    invisible(.Call('_rwmcouplings_SetSeed_cpp', PACKAGE = 'rwmcouplings', seed, stream))
}

#' Sample the latent variables from the model
#'
#' @export
SampleLatentVariables <- function(T, sig, phi) {
    .Call('_rwmcouplings_SampleLatentVariables', PACKAGE = 'rwmcouplings', T, sig, phi)
}

#' RWM algorithm
#'
#' @export
svmRWM <- function(x0, y, beta, sigma, phi, h, iter) {
    .Call('_rwmcouplings_svmRWM', PACKAGE = 'rwmcouplings', x0, y, beta, sigma, phi, h, iter)
}

#' RWM algorithms, with CRN coupling
#'
#' @export
svmCRNRWM <- function(x0, y0, y_data, beta, sigma, phi, h, iter, thin) {
    .Call('_rwmcouplings_svmCRNRWM', PACKAGE = 'rwmcouplings', x0, y0, y_data, beta, sigma, phi, h, iter, thin)
}

#' RWM algorithms, with reflection-maximal coupling
#'
#' @export
svmReflMaxRWM <- function(x0, y0, y_data, beta, sigma, phi, h, iter, thin) {
    .Call('_rwmcouplings_svmReflMaxRWM', PACKAGE = 'rwmcouplings', x0, y0, y_data, beta, sigma, phi, h, iter, thin)
}

#' RWM algorithms, with two-scale coupling using GCRN when (|x-y|^2 >= thresh) and using reflection-maximal coupling when (|x-y|^2 < thresh)
#'
#' @export
svmTwoScaleGCRNRWM <- function(x0, y0, y_data, beta, sigma, phi, h, iter, thin, thresh) {
    .Call('_rwmcouplings_svmTwoScaleGCRNRWM', PACKAGE = 'rwmcouplings', x0, y0, y_data, beta, sigma, phi, h, iter, thin, thresh)
}

#' GCRN-coupled RWM algorithms targeting SVM and a Laplace approximation
#'
#' @export
svmLaplaceGCRNRWM <- function(x0, y0, y_data, beta, sigma, phi, mu, Omega, Omega_chol, h, iter, thin) {
    .Call('_rwmcouplings_svmLaplaceGCRNRWM', PACKAGE = 'rwmcouplings', x0, y0, y_data, beta, sigma, phi, mu, Omega, Omega_chol, h, iter, thin)
}

#' CRN-coupled RWM algorithms targeting SVM and a Laplace approximation
#'
#' @export
svmLaplaceCRNRWM <- function(x0, y0, y_data, beta, sigma, phi, mu, Omega, Omega_chol, h, iter, thin) {
    .Call('_rwmcouplings_svmLaplaceCRNRWM', PACKAGE = 'rwmcouplings', x0, y0, y_data, beta, sigma, phi, mu, Omega, Omega_chol, h, iter, thin)
}

#' Reflection-coupled RWM algorithms targeting SVM and a Laplace approximation
#'
#' @export
svmLaplaceReflRWM <- function(x0, y0, y_data, beta, sigma, phi, mu, Omega, Omega_chol, h, iter, thin) {
    .Call('_rwmcouplings_svmLaplaceReflRWM', PACKAGE = 'rwmcouplings', x0, y0, y_data, beta, sigma, phi, mu, Omega, Omega_chol, h, iter, thin)
}

