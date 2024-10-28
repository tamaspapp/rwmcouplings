library(doRNG)
library(doParallel)
library(rwmcouplings)

if(!file.exists("binreg-prelim.RData")) {source("binreg-prelim.R")}
load("binreg-prelim.RData")

seed <- 144169
SetSeed_cpp(seed)
set.seed(seed)

######
# Meeting times, on a grid of step sizes and thresholds
######
ncores <- 7
cl <- parallel::makeCluster(ncores) # Start parallel processing
doParallel::registerDoParallel(cl)

R_abl <- 50

# Sample initial values from Gaussian approximation to target
L <- t(chol(Sigma))
x0s <- foreach(i = 1:R_abl) %dorng% {mu + as.vector(L %*% rnorm(d))} 
y0s <- foreach(i = 1:R_abl) %dorng% {mu + as.vector(L %*% rnorm(d))}

run_ablation_globalscale <- function(setting, algorithm, step_sizes, threshholds,
                                     seed_ = seed, target_ = target, x0s_ = x0s, y0s_ = y0s, iter_, R_ = R_abl) {
  # Switching thresholds on the scale of delta ~ sum(|| L^(-1) (x_mean - y_mean)||^2), 
  # where proposal covariance is Sigma = h^2 L^\top L.
  
  foreach(l = step_sizes, .combine = "rbind") %:%
  foreach(thresh = threshholds, .combine = "rbind") %:%
  foreach(i = 1:R_, x0_ = x0s_, y0_ = y0s_, .packages = "rwmcouplings", .combine = "rbind") %dopar% {
    SetSeed_cpp(seed_, i)
    tau <- algorithm(target_, list("Sigma" = l^2 * setting$cov, "thresh" = thresh / l^2), x0_, y0_, iter_)$tau
    return(data.frame("l" = l, "thresh" = thresh, "tau" = tau))
    }
}

run_ablation_proposalscale <- function(setting, algorithm, step_sizes, threshholds,
                                       seed_ = seed, target_ = target, x0s_ = x0s, y0s_ = y0s, iter_, R_ = R_abl) {
  # Switching thresholds on the scale of delta ~ sum(|| (hL)^(-1) (x_mean - y_mean)||^2),
  # where proposal covariance is Sigma = (hL)^\top (hL).
  foreach(l = step_sizes, .combine = "rbind") %:%
    foreach(thresh = threshholds, .combine = "rbind") %:%
    foreach(i = 1:R_, x0_ = x0s_, y0_ = y0s_, .packages = "rwmcouplings", .combine = "rbind") %dopar% {
      SetSeed_cpp(seed_, i)
      tau <- algorithm(target_, list("Sigma" = l^2 * setting$cov, "thresh" = thresh), x0_, y0_, iter_)$tau
      return(data.frame("l" = l, "thresh" = thresh, "tau" = tau))
    }
}

# Diagonal proposal variance
iter_max <- 1e7 # Reduce to e.g. 2e5 for reasonable execution times
mala_abl       <- run_ablation_globalscale(DiagVar, mala_twoscalecrn_cpp, DiagVar$l_mala, c(10, 100, Inf), iter_ = iter_max)
rwm_abl_gcrefl <- run_ablation_proposalscale(DiagVar, rwm_twoscalegcrefl_cpp, DiagVar$l_rwm, c(0.1, 1, 10, 100), iter_ = iter_max) # Our intuition is that GCRefl is a strict improvement over ReflMax unless the chains are close enough to meet
rwm_abl_gcrn   <- run_ablation_globalscale(DiagVar, rwm_twoscalegcrn_cpp, DiagVar$l_rwm, c(10, 100, Inf), iter_ = iter_max)

save(mala_abl, rwm_abl_gcrefl, rwm_abl_gcrn, iter_max,
     file = "binreg-ablation-diagvar.RData")

# Full proposal variance
# iter_max <- 1e7
mala_abl       <- run_ablation_globalscale(FullVar, mala_twoscalecrn_cpp, FullVar$l_mala, c(1, 10, Inf), iter_ = iter_max)
rwm_abl_gcrefl <- run_ablation_globalscale(FullVar, rwm_twoscalegcrefl_cpp, FullVar$l_rwm, c(0.01, 0.1, 1, 10, Inf), iter_ = iter_max)
rwm_abl_gcrn   <- run_ablation_globalscale(FullVar, rwm_twoscalegcrn_cpp, FullVar$l_rwm, c(1, 10, Inf), iter_ = iter_max)

save(mala_abl, rwm_abl_gcrefl, rwm_abl_gcrn, iter_max,
     file = "binreg-ablation-fullvar.RData")

parallel::stopCluster(cl) # Stop parallel processing

