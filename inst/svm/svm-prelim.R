########
# Preliminary steps: params, RNG seeding, data generation, Laplace approximation
########
library(rwmcouplings)
library(Matrix)

####
# Auxiliary functions
####
# "Combine" function for "foreach", to get distances and meeting times separately
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

# Calculate the TVD bound at given iterations (assumes all Y-chains meet before their time L)
getTVDBound <- function(iters, taus) {
  
  out <- rep(NA, length(iters))
  for(i in 1:length(iters)) {
    out[i] <- mean(taus>iters[i])
  }
  out
}
####



########
# Set up model and score
########

####
# Evaluate log-density and score for the stochastic volatility model
svmLoglike <- function(x, y_data, beta = 0.65, phi = 0.98, sig = 0.15) {
  -0.5 * (sum(x) + 
            (1/beta^2) * sum(y_data^2 * exp(-x)) + 
            (1/sig^2) * sum((phi * x[-length(x)] - x[-1])^2) + 
            (1-phi^2) * x[1]^2 / sig^2)
}
svmGradLoglike <- function(x, y_data, beta = 0.65, phi = 0.98, sig = 0.15) {
  diff_x <- phi * x[-length(x)] - x[-1]
  temp <- 0.5 * (1/beta^2) * y_data^2 * exp(-x)
  temp[1] <- temp[1] - ((1-phi^2) / sig^2) * x[1]
  temp[-1] <- temp[-1] + (1/sig^2) * diff_x
  temp[-length(temp)] <- temp[-length(temp)] - (phi / sig^2) * diff_x
  return(temp-0.5)
}
####

####
# Model parameters
t <- 360L # Dimension
beta <- 0.65
sig  <- 0.15
phi  <- 0.98

# Generate model data 
seed <- 12345; set.seed(seed); SetSeed_cpp(seed)
y_data <- beta * rnorm(t) * exp(0.5 * SampleLatentVariables(t, sig, phi)) # Model data

# Save log-density and score
logpi <- function(x){svmLoglike(x, y_data, beta, phi, sig)}
gradlogpi <- function(x){svmGradLoglike(x, y_data, beta, phi, sig)}
####


########
# Fit Laplace approximation
########
x0_optim <- SampleLatentVariables(t, sig, phi)
optim_out <- optim(par=x0_optim, fn=logpi, gr=gradlogpi,
                   method = "L-BFGS-B",
                   hessian = TRUE,
                   control = list(fnscale = -1,  maxit = 2000))

optim_out$hessian[abs(optim_out$hessian)<1e-6] <- 0 # The Hessian is sparse; make the approximation sparse too.
Omega_optim <- as(as(as(-optim_out$hessian, "dMatrix"), "symmetricMatrix"), "CsparseMatrix")
Sigma_optim <- solve(Omega_optim)

# Laplace approximation parameters
mu <- optim_out$par # Mean
Omega <- as(Omega_optim, "dgCMatrix") # Precision
Omega_chol <- as(chol(Omega_optim), "dgCMatrix") # Upper triangular Cholesky factor of precision
Sigma_cholLT <- t(chol(Sigma_optim)) # Lower triangular Cholesky factor of covariance

# Average traces of precision and covariance matrices
avg_prec <- mean(diag(Omega_optim))
avg_cov  <- mean(diag(Sigma_optim))

print(paste("Eccentricity of Laplace approximation is:", signif(avg_prec*avg_cov)))

laplace_cov  <- as.matrix(Sigma_optim)
laplace_mean <- as.vector(mu)
########
# Save output
########
save(t, beta, sig, phi, y_data, svmLoglike, svmGradLoglike, logpi, gradlogpi, 
     mu, Omega, Omega_chol, Sigma_cholLT,
     laplace_cov, laplace_mean,
     avg_prec, avg_cov,
     comb, getTVDBound,
     file = "svmData.RData")

