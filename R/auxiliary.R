

#####
# Expectation E[1 ^ exp(lx^{1/2} Z_1 - l^2/2) ^ exp(ly^{1/2} Z_2 - l^2/2)],
# under (Z_1, Z_2) ~ BvN(rho).
#####
g <- function(x,y,rho,l) {
  
  # Deal with rho very close to 1 or very close to -1. Otherwise, the calculations seem numerically stable.
  if(rho > 1-1e-12) {
    return(2*pnorm(-0.5*l))
  } else if(rho < -1+1e-12) {
    return(g(x,y,-1+1e-12,l))
  }
  
  # Owen (1980)'s (10,011.1)
  owen_10010.1 <- function(Y,a,b){
    Sigma <- diag(2); Sigma[1,2] <- Sigma[2,1] <- (-b)/sqrt(1+b^2)
    mvtnorm::pmvnorm(upper = c(a/sqrt(1+b^2), Y), corr = Sigma)
  }
  
  h <- function(x,y,rho,l) {
    b <- -(sqrt(x/y) - rho)/sqrt(1 - rho^2)
    a <- b * l * sqrt(x)
    U <- 0.5 * l / sqrt(x) - l * sqrt(x)
    exp(0.5*l^2*(x-1)) * owen_10010.1(U, a, b)
  }
  
  Sigma <- diag(2); Sigma[1,2] <- Sigma[2,1] <- rho
  lb <- 0.5*l/sqrt(c(x,y))
  mvtnorm::pmvnorm(lower = lb, corr = Sigma) + h(x,y,rho,l) + h(y,x,rho,l) 
}


#####
# Return an estimate of the squared distance between chains under reflection and
# CRN couplings. Accurate if the curvature of the target is approximately constant.
#####
## Parameters:
#
# l = step length of RWM algorithm with proposal variance l^2/d
# avg_cov = estimate of average eigenvalue of the covariance
# avg_precision = estimate of average eigenvalue of precision matrix
# rho_lim = grid of correlations of gradient projections (in [0,1)) to search over with "uniroot"
#
#' @export
squareDistAsymptote <- function(l, d, avg_cov, avg_prec, 
                                rho_lim = c(0, 1-1e-5), plot = TRUE, rho_plot = seq(0, 1, 0.001)) {
  l1sq <- l^2 * avg_prec # Natural parameter
  ellipt <- 1 / (avg_cov * avg_prec) # Measure of ellipticity, 1 if spherically symmetric
  
  b_crn  <- 2 * pnorm(-0.5*sqrt(l1sq))
  b_refl <- 2 / (1 - ellipt) * pnorm(-0.5*sqrt(l1sq))
  
  fun_crn  <- function(rho) {g(1,1,rho,sqrt(l1sq))[1] - b_crn * rho}
  fun_refl <- function(rho) {g(1,1,rho,sqrt(l1sq))[1] - b_refl * (rho - ellipt) }
  
  # Get solution for CRN
  rho_crn <- uniroot(fun_crn, rho_lim, extendInt = c("downX"),tol = 1e-9, maxiter = 1e6)$root
  v_crn <- rho_crn

  if(ellipt == 1) {
    print("Spherically symmetric target, reflection-coupled chains asymptote at 0.")
    v_refl <- 1
  } else {  
    rho_refl <- uniroot(fun_refl, rho_lim, extendInt = c("downX"),tol = 1e-9, maxiter = 1e6)$root
    v_refl <-  (rho_refl - ellipt) / (1 - ellipt)
  }
  
  if(plot){
    gs <- rep(NA, length(rho_plot))
    for(i in seq_along(rho_plot)){
      gs[i] <- g(1,1, rho_plot[i],sqrt(l1sq))
    }
    
    if(ellipt != 1) {
      plot(rho_plot, gs, xlab = "Correlation of gradients", main = "Diagnostic plot: stable fixed point of ODE limit", ylab = "Functions to intersect", type = "l")
      lines(rho_plot, b_crn * rho_plot, col = "red")
      abline(v = rho_crn, lty = 2, col = "red")
      legend(x = "topleft", legend = c("CRN", "Reflection", "uniroot solution"), col = c("red", "blue", "black"), lty = c(1,1,2))
      lines(rho_plot, b_refl * (rho_plot - ellipt), col = "blue")
      abline(v = rho_refl, lty = 2, col = "blue")
      
    } else {
      plot(rho_plot, gs, xlab = "Correlation of gradients", main = "Diagnostic plot: stable fixed point of ODE limit", ylab = "Functions to intersect", type = "l")
      lines(rho_plot, b_crn * rho_plot, col = "red")
      abline(v = rho_crn, lty = 2, col = "red")
      legend(x = "topleft", legend = c("CRN", "uniroot solution"), col = c("red", "black"), lty = c(1, 2))
    }
  }

  return(list("crn"  = 2*avg_cov*d*(1 - v_crn), 
              "refl" = 2*avg_cov*d*(1 - v_refl))) 
}



## Parameters:
#
# l = step length of RWM algorithm with proposal variance l^2/d
# avg_cov = estimate of average eigenvalue of the covariance
# avg_precision = estimate of average eigenvalue of precision matrix
#
# rho = grid of correlations of gradient projections (in [0,1]) to search over.
# n = Monte Carlo sample size for numerical evaluation of expectation
# plot = logical, whether to output diagnostic plots
# correct = logical, whether to correct for error at rho = 1 or not
#
# squareDistAsymptote <- function(l, d, avg_cov, avg_prec, 
#                                 rho = c(seq(0, 1, 0.01)), n = 1e7, 
#                                 plot = TRUE, correct = TRUE) {
#   l1sq <- l^2 * avg_prec
#   
#   MCavg <- function(rho,n) {
#     z1 <- rnorm(n)
#     z <- rnorm(n)
#     
#     out <- rep(NA, length(rho))
#     for (r in seq_along(rho)){
#       z2 <- rho[r] * z1 + sqrt(1 - rho[r]^2) * z
#       out[r] <- mean(exp(pmin(0, -sqrt(l1sq)*z1-0.5 * l1sq, -sqrt(l1sq)*z2-0.5 * l1sq)))
#     }
#     # Correct by forcing LHS to take the analytically correct value when the correlation is 1
#     offset <- 2 * pnorm(-0.5*sqrt(l1sq)) - mean(exp(pmin(0, -sqrt(l1sq)*z1-0.5 * l1sq)))
#     out <- out + offset
#     
#     return(out)
#   }
#   lhs <- MCavg(rho, n)
#   
#   # lhs <- MCavg(rho, l1sq, n) # Called from C++ code
#   
#   # Get solution for CRN
#   rhs_crn  <- 2 * rho * pnorm(-0.5*sqrt(l1sq))
#   if(sum(lhs >= rhs_crn) == length(rhs_crn)) {
#     rho_crn <- 1
#   } else {
#     rho_crn <- rho[Position(isFALSE, lhs > rhs_crn, right=FALSE) - 1] # Stop at last "true"
#   }
#   
#   # Get solution for reflection coupling
#   ellipt <- 1/ (avg_cov * avg_prec) # Measure of ellipticity
#   
#   if(ellipt == 1) {
#     print("Spherically symmetric target, reflection-coupled chains asymptote at 0.")
#     v_refl <- 1
#   } else {
#     
#     rhs_refl <- 2 * (rho - ellipt) / (1 - ellipt) * pnorm(-0.5*sqrt(l1sq))
#     
#     # Get solution for reflection
#     if(sum(lhs >= rhs_refl) == length(rhs_refl)) {
#       rho_refl <- 1
#     } else {
#       rho_refl <- rho[Position(isFALSE, lhs > rhs_refl, right=FALSE) - 1]
#     }
#     v_refl <-  (rho_refl - ellipt) / (1 - ellipt)
#   }
#   
#   if(plot){
#     if(ellipt != 1) {
#       par(mfrow = c(1,2))
#       plot(rho, lhs, xlab = "Correlation of gradients", main = "CRN coupling", ylab = "Functions to intersect")
#       abline(v = rho_crn, lty = 2, col = "red")
#       lines(rho, rhs_crn, col = "red")
#       
#       plot(rho, lhs, xlab = "Correlation of gradients", main = "Reflection coupling", ylab = "Functions to intersect")
#       abline(v = rho_refl, lty = 2, col = "red")
#       lines(rho, rhs_refl, col = "red")
#       par(mfrow = c(1,1))
#     } else {
#       plot(rho, lhs, xlab = "Correlation of gradients", main = "CRN coupling", ylab = "Functions to intersect") 
#       abline(v = rho_crn, lty = 2, col = "red")
#       lines(rho, rhs_crn, col = "red")
#     }
#     
#   }
#   
#   return(list("crn"  = 2*avg_cov*d*(1 - rho_crn), 
#               "refl" = 2*avg_cov*d*(1 - v_refl))) 
# }


