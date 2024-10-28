########################################################################### 
# Calculate Gelbrich (1990) lower bound on W2^2 for Laplace approximation #
###########################################################################
# - We do many MCMC runs of an efficient algorithm.
# - We calculate a sample mean and covariance from each of them.
# - We then average the sample means and covariances.
# - This gives us the posterior mean and covariance, plus error bars since
# we're i.i.d. over the replicates.
#########################################################################

library(rwmcouplings)
library(Rcpp)
library(Matrix)
library(doParallel)
library(doRNG)
library(reshape2)
library(ggplot2)
library(latex2exp)
library(ggpubr)

# Sample using Hug and Hop
hughop_ <- function(x0,Time,Bounces,lam,kap,iter,logpi,gradlogpi){
  
  normalize <- function(x) {return((1/sqrt(sum(x^2)))*x)}
  mu = sqrt(lam * kap);
  
  d = length(x0);
  
  delta = Time / Bounces;
  
  x = x0;
  
  logpi_x = logpi(x);
  
  acc_x_hug = 0; acc_x_hop = 0;
  xs = matrix(NA, nrow = iter+1, ncol = d);
  xs[1,] = x;
  for(i in 1:iter) {
    # Hug
    # Propose: bounce B times
    z = rnorm(d); # X-velocity
    # Initial position half-step
    xp = x + 0.5 * delta * z;
    g_xp = normalize(gradlogpi(xp));
    
    for(b in 1:Bounces)
    {
      # Reflect velocity in gradient
      z = z -  2 * sum(z * g_xp) * g_xp;
      # Full step
      xp = xp + delta * z;
      g_xp = normalize(gradlogpi(xp));
    }
    # Went too far
    xp = xp - 0.5 * delta * z;
    
    # Accept-reject
    logpi_xp = logpi(xp);
    logHR_x = logpi_xp - logpi_x;
    
    log_u = log(runif(1));
    if (logHR_x > 0 || log_u < logHR_x)
    {
      x = xp;
      logpi_x = logpi_xp;
      acc_x_hug = acc_x_hug + 1;
    }
    
    # Hop
    g_x = gradlogpi(x); gnx = sum(g_x^2);
    
    z = rnorm(d);
    z1 = sum(z * g_x)/gnx *g_x;
    xp <- x + lam/sqrt(gnx) * z1 + mu/sqrt(gnx) *(z - z1)
    
    ###
    # Accept-reject
    ###
    g_xp = gradlogpi(xp); gnxp = sum(g_xp^2);
    logpi_xp = logpi(xp);
    logHR_x = logpi_xp - logpi_x + 0.5*d*log(gnxp/gnx) - 0.5/mu^2*sum((xp-x)^2)*(gnxp - gnx) - 0.5*(1/lam^2 - 1/mu^2)*(sum((xp - x)*g_xp)^2 - sum((x - xp)*g_x)^2);
    
    log_u = log(runif(1))
    if (logHR_x > 0 || log_u < logHR_x) {
      x = xp;
      logpi_x = logpi_xp;
      acc_x_hop = acc_x_hop + 1;
    }
    xs[i + 1,] = x;
  }
  return(list("x" = x,
              "xs" = xs,
              "acc_hop" = acc_x_hop/iter, 
              "acc_hug" = acc_x_hug/iter))
}

#####
# Squared Wasserstein distance lower bound: Gelbrich (1990)
#####
# - Partly adapted from from: https://github.com/niloyb/BoundWasserstein/blob/main/estimators.R
# - Matrix square roots are calculated using Singular Value Decomposition (SVD).
#

# Do some matrix operations in C++
Rcpp::cppFunction('
Eigen::MatrixXd cpp_crossprod(Eigen::Map<Eigen::MatrixXd> X){
  return (X.transpose())*X;
}', depends = 'RcppEigen')
Rcpp::cppFunction('
Eigen::MatrixXd cpp_prod(Eigen::Map<Eigen::MatrixXd> X, 
                         Eigen::Map<Eigen::MatrixXd> Y){
  return X*Y;
}', depends = 'RcppEigen')

w2sq_gelbrich <- function(mu1, Sigma1, mu2, Sigma2){ 
  svd1 <- svd(Sigma1)
  svd2 <- svd(Sigma2)
  
  sqrtSigma1 <- cpp_crossprod(cpp_prod(diag(svd1$d^0.25),t(svd1$u)))
  sqrtSigma2 <- cpp_crossprod(cpp_prod(diag(svd2$d^0.25),t(svd2$u)))
  
  trace_Sigma1 <- sum(diag(Sigma1))
  trace_Sigma2 <- sum(diag(Sigma2))
  
  # Cross term is trace((M^T * M)^(1/2)), where M = Sigma1^(1/2) * Sigma2^(1/2)
  # Use trick: if SVD of M is U * D * V^T, then M^T * M = V * D^2 * V^T
  # So that (M^T * M)^(1/2) = V^(1/2) * D * {V^(1/2)}^T
  # So that the trace of trace((M^T * M)^(1/2)) = trace(D)
  M <- cpp_prod(sqrtSigma1, sqrtSigma2)
  CrossSigma_svd <- svd(M)
  trace_sqrtCrossSigma <- sum(CrossSigma_svd$d)
  
  out <- sum((mu1-mu2)^2) + trace_Sigma1 + trace_Sigma2 - 2 * trace_sqrtCrossSigma
  
  return(abs(out))
} # Sanity check: have compared with implementation using explicit square-roots. 



####################

# Load in model data 
if(!file.exists("svmData.RData")) {source("svm-prelim.R") } # Ensure we run this first
load(file = "svmData.RData")

# Replicates; number of iterations for each replicate
R_gelb <- 50
iter_gelb <- 5e4


# Parallel computing setup
ncores <- 6
cl <- parallel::makeCluster(ncores) 
doParallel::registerDoParallel(cl)

# RNG setup
seed <- 12345; set.seed(seed); SetSeed_cpp(seed) 

# MCMC parameters
Time <- 0.5; Bounces <- 10 # Hug params
lam <- 20; kap <- 1 # Hop params

x0_gelb <- vector(mode = "list", R_gelb)
for(i in 1:R_gelb){x0_gelb[[i]] <- mu + as.vector(Sigma_cholLT%*%rnorm(t))}

gelb_out <- 
  foreach(x0 = x0_gelb, 
          .combine='comb',.multicombine=TRUE,.init=list(list(), list())) %dorng% {
  out <- hughop_(x0,Time,Bounces,lam,kap,iter_gelb,logpi,gradlogpi)
  list(cov(out$xs), colMeans(out$xs))
}; names(gelb_out) <- c("cov", "mean")
parallel::stopCluster(cl)



# Sample means
pi_cov  <- Reduce("+", gelb_out$cov)  / length(gelb_out$cov)
pi_mean <- Reduce("+", gelb_out$mean) / length(gelb_out$mean)

w2sq_gelbrich_lb <- w2sq_gelbrich(pi_mean,pi_cov,laplace_mean,laplace_cov)
w2sq_gelbrich_lb


# Jackknife bias and standard error
# NB: could make this faster by only computing the SVD of cov2 once.
w2sq_jack <- function(means, covs, mean2, cov2) {
  w2sqs <- rep(NA,length(covs))
  for(i in 1:length(covs)){
    pi_mean <- Reduce("+", means[-i]) / length(means[-i])
    pi_cov  <- Reduce("+", covs[-i])  / length(covs[-i])
    w2sqs[i] <- w2sq_gelbrich(pi_mean,pi_cov,mean2,cov2)
  }
  w2sqs
}
w2sqs <- w2sq_jack(gelb_out$mean, gelb_out$cov,laplace_mean,laplace_cov)


w2sq_gelbrich_jackbias <- (R_gelb-1) * (mean(w2sqs) - w2sq_gelbrich_lb)
w2sq_gelbrich_jackse <- sqrt((R_gelb-1)^2/R_gelb) * sd(w2sqs)
# NB: can obtain explicit standard errors: Rippl. et al (2016), equation (8).

w2sq_gelbrich_jackbias
w2sq_gelbrich_jackse

save(w2sq_gelbrich_lb, 
     w2sq_gelbrich_jackbias,
     w2sq_gelbrich_jackse, file = "laplaceBiasGelbrich.RData")

