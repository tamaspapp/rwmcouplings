library(doParallel)
library(rwmcouplings)

######
# 1. Preprocess data, set prior
######

# Preprocess the data (follow Gelman et al (2008) - "A weakly informative default prior...")
#   Center covariates, then scale to 0.5 standard deviation
#   Add intercept to the covariates
preprocess <- function(data_frame) {
  # Center covariates, scale them to 0.5 standard deviation
  X <- as.matrix(data_frame[, -ncol(data_frame)]) # Remove response "Y"
  X <- 0.5 * scale(X, center = T, scale = T) # 
  # Add intercept
  X <- as.matrix(cbind(rep(1, nrow(X)), X))  
  
  # Response
  Y <- as.matrix(data_frame[, ncol(data_frame)]) 
  return(list("X" = X, "Y" = Y))
}
######

data.p <- preprocess(read.csv("sonar.csv", header = F))

X <- data.p$X
Y <- data.p$Y
yX <- sweep(X, 1, (2*Y-1), "*");  colnames(yX) <- NULL

# "lambda" is the standard deviation of the spherical Gaussian prior
lambda <- 5
d <- ncol(yX)
n <- nrow(yX)

target <- 
  list(
    "target_type" = "logistic_regression",
    "yX" = yX,
    "lambda" = lambda
  )


# RNG setup
seed <- 144169
SetSeed_cpp(seed)
set.seed(seed)

# Two settings: full and diagonal proposal covariance matrices
DiagVar <- list()
FullVar <- list()

######
# 2. Obtain a diagonal preconditioner.
######
# Use a preliminary run, MALA with spherical proposals
iter <- 5e5
burnin <- 1e4 # Discard initialization and some burn-in

mala_out <- mala_cpp(target, list("Sigma" = rep(0.27^2, d)), rnorm(d), iter, 1); mala_out$acc_rate

# # ESS calculation for preliminary run
# library(coda)
# summary(effectiveSize(mala_out$xs))

# Sample from the target
x_infty <- mala_out$xs[iter + 1,]

# Target mean and covariance
remove_burnin <- -(0:burnin)-1
mu    <- colMeans(mala_out$xs[remove_burnin, ])
Sigma <- cov(mala_out$xs[remove_burnin, ])

# Proposal covariances
DiagVar$cov <- diag(Sigma)
FullVar$cov <- Sigma

######
# 3. Check the variation in the Hessian term of the acceptance ratio
######
# Eigenvalues of the precision matrix, accounting for preconditioning
D_inv <- diag(1/sqrt(DiagVar$cov))
v <- 1 / eigen(D_inv %*% Sigma %*% D_inv, symmetric = T, only.values = T)$values
print("Ratio (max eigenvalue)/(trace), for precision matrix:")
print(max(v) / sum(v))


######
# 4. Consider a grid of step sizes, calculate acceptance rates on the grid
######
ncores <- 4 # Set up parallel processing
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

iter_acc <- 1e6

find_acceptance_rate <- function(setting, algorithm, step_sizes, 
                                 seed_ = seed, target_ = target, x0_ = x_infty, iter_ = iter_acc){
  acc_rates <- 
    foreach(l_ = step_sizes, .packages = "rwmcouplings", .combine = "c") %dopar% {
      SetSeed_cpp(seed_); algorithm(target_, list("Sigma" = l_^2 * setting$cov), x0_, iter_, 1)$acc_rate
    }
  print("Acceptance rates:")
  print(acc_rates)
  
  return(acc_rates)
}

DiagVar$l_rwm <- c(0.045, 0.05, 0.06, 0.07, 0.08, 0.09)
DiagVar$acc_rate_rwm <- find_acceptance_rate(DiagVar, rwm_cpp, DiagVar$l_rwm)
DiagVar$l_mala <- c(0.07, 0.08, 0.09, 0.1, 0.11, 0.12)
DiagVar$acc_rate_mala <- find_acceptance_rate(DiagVar, mala_cpp, DiagVar$l_mala)

FullVar$l_rwm <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35)
FullVar$acc_rate_rwm <- find_acceptance_rate(FullVar, rwm_cpp, FullVar$l_rwm)
FullVar$l_mala <- c(0.45, 0.55, 0.65, 0.75, 0.8)
FullVar$acc_rate_mala <- find_acceptance_rate(FullVar, mala_cpp, FullVar$l_mala)

parallel::stopCluster(cl) # Stop parallel processing

compile_acceptance_rates <- function(setting){
  acc_rates <- rbind(data.frame("l" = setting$l_mala, "algorithm" = "MALA", "acc" = setting$acc_rate_mala),
                     data.frame("l" = setting$l_rwm, "algorithm" = "RWM", "acc" = setting$acc_rate_rwm))
  return(acc_rates)
}

FullVar$acc_rates <- compile_acceptance_rates(FullVar)
DiagVar$acc_rates <- compile_acceptance_rates(DiagVar)

######
# 5. Timing
######
iter_time <- 2e6

time_per_iteration <- function(setting, algorithm, step_size, 
                               target_ = target, x0_ = x_infty, iter_ = iter_time){
  t0 <- Sys.time()
  out <- algorithm(target_, list("Sigma" = step_size^2 * setting$cov), x0_, iter_, 1)
  time <- as.double(Sys.time() - t0) / iter_
  
  print("Time per iteration (in seconds):")
  print(time)
  
  return(time)
}

DiagVar$time_rwm <- time_per_iteration(DiagVar, rwm_cpp, 0.09)
DiagVar$time_mala <- time_per_iteration(DiagVar, mala_cpp, 0.11)

FullVar$time_rwm <- time_per_iteration(FullVar, rwm_cpp, 0.6)
FullVar$time_mala <- time_per_iteration(FullVar, mala_cpp, 0.8)

# Save data
save(mu, Sigma, # Target summaries
     DiagVar, FullVar, # Two settings for MCMC
     target, yX, lambda, d, n,
     file = "binreg-prelim.RData")
