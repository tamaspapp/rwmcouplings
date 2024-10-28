# This code was run on a compute cluster. #

library(doRNG)
library(doParallel)
library(rwmcouplings)

if(!file.exists("binreg-prelim.RData")) {source("binreg-prelim.R")}
load("binreg-prelim.RData")

# "Combine" functions for "foreach": get separate lists for separate outputs
comb_targetsample <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}
comb_coupling <- function(x, ...) {  
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}

seed <- 144169
ncores <- 56
R <- 1.12e5

######
# 1. Obtain samples from "the target" and calculate posterior mean and variance.
# We use several long enough, warm-started, MCMC runs.
######
SetSeed_cpp(seed)
set.seed(seed)

iter   <- 5e3
burnin <- 1e3

sample_target <- function(R, iter, burnin, mala_proposal_var, target, mu, Sigma, seed, ncores_ = ncores){
  
  cl <- parallel::makeCluster(ncores) 
  doParallel::registerDoParallel(cl)
  
  d <- length(mu)
  
  # Warm-start from Gaussian approximation to posterior
  L <- t(chol(Sigma))
  x0s <- foreach(i = 1:R) %dorng% {mu + as.vector(L %*% rnorm(d))}
  
  par_out <- 
    foreach(i = 1:R, x0_ = x0s,
            .init = list(list(),list(),list()), .combine = comb_targetsample, 
            .packages = "rwmcouplings") %dopar% {
              SetSeed_cpp(seed, i)
              samples <- mala_cpp(target, list("Sigma" =  mala_proposal_var), x0_, iter, thin = 1)$xs
              list(samples[iter + 1, ], colMeans(samples[-(0:burnin + 1), ]), apply(samples[-(0:burnin + 1), ],2,var))
            }
  
  
  # Store samples from the target
  xinfinitys <- par_out[[1]]
  
  # Get posterior summaries
  means <- do.call(rbind, par_out[[2]])
  vars <- do.call(rbind, par_out[[3]])
  posterior <- data.frame("mean" = colMeans(means),
                          "mean.stderr" = apply(means, 2, sd) / sqrt(R),
                          "var" = colMeans(vars), 
                          "var.stderr" = apply(vars, 2, sd) / sqrt(R),
                          "coord" = 1:d)
  
  parallel::stopCluster(cl)
  
  return(list("xinfinitys" = xinfinitys, 
              "posterior" = posterior))
}

t_in <- Sys.time()
targetsample_out <- sample_target(R, iter, burnin, 0.65^2 * FullVar$cov, target, mu, Sigma, seed)
t_diff <- Sys.time() - t_in

xinfinitys <- targetsample_out$xinfinitys
posterior <- targetsample_out$posterior

save(xinfinitys, posterior,
     R, t_diff,
     file = "binreg-targetsample.RData")

######
# 2. Get meeting times (in terms of iterations and wall-time) and near-unbiased
# asymptotic variance estimates.
######
SetSeed_cpp(seed)
set.seed(seed)

getAsymptoticVariance <- function(x0, fishy0, posterior) {
  # EPAVE with t=0, equation 3.2 of (Douc et al, 2022, "Solving the Poisson equation with coupled Markov chains")
  2*(x0 - posterior$mean) * fishy0 - posterior$var 
}

y0 <- posterior$mean # Start Y-chain from the posterior mean
iter_max <- 2e6 # The replicate count R is specified above

run_coupling_globalscale <- function(setting, step_sizes, algorithm, thresh,
                                     target_ = target, x0s_ = xinfinitys, y0_ = y0, iter_ = iter_max, R_ = R, 
                                     comb_ = comb_coupling, seed_ = seed, asympvar = getAsymptoticVariance,
                                     posterior_ = posterior, ncores_ = ncores){
  
  cl <- parallel::makeCluster(ncores_)
  doParallel::registerDoParallel(cl)
  
  foreach_out <-
    foreach(l = step_sizes, .combine = comb_) %:%
    foreach(i = 1:R_, x0_ = x0s_,
            .packages = "rwmcouplings", .combine = comb_, .multicombine=TRUE) %dopar% {
              
              SetSeed_cpp(seed_, stream = i + R_)
              t0  <- Sys.time()
              out <- algorithm(target_, list("Sigma" = l^2 * setting$cov, "thresh" = thresh / l^2), x0_, y0_, iter_)
              t1  <- Sys.time()
              
              list(data.frame("l" = l,"tau" = out$tau,"time" = t1 - t0),
                   data.frame("l" = l,"coord" = 1:length(x0_),"asympvar" = asympvar(x0_, out$fishy, posterior_)))
            }
  parallel::stopCluster(cl)
  
  return(foreach_out)  
}

run_coupling_proposalscale <- function(setting, step_sizes, algorithm, thresh,
                                       target_ = target, x0s_ = xinfinitys, y0_ = y0, iter_ = iter_max, R_ = R, 
                                       comb_ = comb_coupling, seed_ = seed, asympvar = getAsymptoticVariance,
                                       posterior_ = posterior, ncores_ = ncores){
  cl <- parallel::makeCluster(ncores_)
  doParallel::registerDoParallel(cl)
  
  foreach_out <-
    foreach(l = step_sizes, .combine = comb_) %:%
    foreach(i = 1:R_, x0_ = x0s_,
            .packages = "rwmcouplings", .combine = comb_, .multicombine=TRUE) %dopar% {

              SetSeed_cpp(seed_, stream = i + R_)
              t0  <- Sys.time()
              out <- algorithm(target_, list("Sigma" = l^2 * setting$cov, "thresh" = thresh), x0_, y0_, iter_)
              t1  <- Sys.time()
              
              list(data.frame("l" = l,"tau" = out$tau,"time" = t1 - t0),
                   data.frame("l" = l,"coord" = 1:length(x0_),"asympvar" = asympvar(x0_, out$fishy, posterior_)))
            }
  
  parallel::stopCluster(cl)
  
  return(foreach_out)   
}

# Diagonal variance matrix #####
t_in <- Sys.time()

rwm_cpl <- run_coupling_proposalscale(DiagVar, DiagVar$l_rwm, rwm_twoscalegcrefl_cpp, 10)
rwm_cpl_meet <- rwm_cpl[[1]]
rwm_cpl_meet$algorithm <- "RWM"
rwm_cpl_asympvar <- rwm_cpl[[2]]
rwm_cpl_asympvar$algorithm <- "RWM"

mala_cpl <- run_coupling_globalscale(DiagVar, DiagVar$l_mala, mala_twoscalecrn_cpp, Inf)
mala_cpl_meet <- mala_cpl[[1]]
mala_cpl_meet$algorithm <- "MALA"
mala_cpl_asympvar <- mala_cpl[[2]]
mala_cpl_asympvar$algorithm <- "MALA"

cpl_meet <- rbind(rwm_cpl_meet, mala_cpl_meet)
cpl_asympvar <- rbind(rwm_cpl_asympvar, mala_cpl_asympvar)

t_diff <- Sys.time() - t_in

save(cpl_meet, cpl_asympvar, t_diff, file = "binreg-coupling-diagvar.RData")
#####


# Full variance matrix #####
t_in <- Sys.time()

rwm_cpl <- run_coupling_globalscale(FullVar, FullVar$l_rwm, rwm_twoscalegcrefl_cpp, thresh = 1, R_ = min(1.12e4, R))
rwm_cpl_meet <- rwm_cpl[[1]]
rwm_cpl_meet$algorithm <- "RWM"
rwm_cpl_asympvar <- rwm_cpl[[2]]
rwm_cpl_asympvar$algorithm <- "RWM"

mala_cpl <- run_coupling_globalscale(FullVar, FullVar$l_mala, mala_twoscalecrn_cpp, thresh = Inf, R_ = min(1.12e4, R))
mala_cpl_meet <- mala_cpl[[1]]
mala_cpl_meet$algorithm <- "MALA"
mala_cpl_asympvar <- mala_cpl[[2]]
mala_cpl_asympvar$algorithm <- "MALA"

cpl_meet <- rbind(rwm_cpl_meet, mala_cpl_meet)
cpl_asympvar <- rbind(rwm_cpl_asympvar, mala_cpl_asympvar)

t_diff <- Sys.time() - t_in

save(cpl_meet, cpl_asympvar, t_diff, file = "binreg-coupling-fullvar.RData")
#####
