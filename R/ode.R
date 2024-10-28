#####
# Compute ODE limits numerically
#####


# Auxiliary functions ####
a <- function(x,l) {
  pnorm(-0.5 * l / sqrt(x)) + (1-2*x) * exp(0.5 * l^2 *(x-1)) * pnorm(0.5 * l/sqrt(x) - l * sqrt(x))
}

# ODE limits ####
# Parameters:
#
# ode_times  = sequence from 0 to "tmax" of times at which to compute the ODE limit
# Other parameters are self-explanatory.

# GCRN coupling ###
b_gcrn <- function(x,y,v,l) {
  m <- pmin(x,y)
  M <- pmax(x,y)
  
  ( pnorm(-0.5 * l / sqrt(m)) +  exp(0.5* l^2 *(m-1)) * ( pnorm(0.5 * l/sqrt(m) - l * sqrt(m)) - pnorm(- l * sqrt(m)) ) + exp(0.5* l^2 *(M-1)) * pnorm(- l * sqrt(M)) ) - 
  #g(x,y,1,l) - 
    v * ( exp(0.5* l^2 *(x-1)) * pnorm(0.5 * l/sqrt(x) - l * sqrt(x)) + 
            exp(0.5* l^2 *(y-1)) * pnorm(0.5 * l/sqrt(y) - l * sqrt(y)) )
} 

#' @export
ode_gcrn <- function(x0,y0,v0,l,ode_times){
  
  limiting_ode <- function(t, y, parms) {
    with(as.list(c(y, parms)), {
      dy1 <- l^2 * a(y[1],l)
      dy2 <- l^2 * a(y[2],l)
      dy3 <- l^2 * b_gcrn(y[1], y[2], y[3],l)
      return(list(c(dy1, dy2, dy3)))
    })
  }
  
  parms <- c(l = l)
  yini <- c(y1 = x0, 
            y2 = y0, 
            y3 = v0)
  
  out <- deSolve::ode(times = ode_times, 
                      func = limiting_ode,
                      y = yini, 
                      parms = parms)
  
  return(list("xsq" = out[, 2],
              "ysq" = out[, 3],
              "innerprod" = out[, 4]))
}



# CRN coupling ###
b_crn <- function(x,y,v,l) {
  if(abs(x + y - 2*v) < 1e-12) {
    rho <- 1 
  } else {
    rho <- v / sqrt(x*y) #pmax(pmin(, 1), -1)
  }
  g(x,y,rho,l) - 
    v * ( exp(0.5* l^2 *(x-1)) * pnorm(0.5 * l/sqrt(x) - l * sqrt(x)) + 
            exp(0.5* l^2 *(y-1)) * pnorm(0.5 * l/sqrt(y) - l * sqrt(y)) )
}

#' @export
ode_crn <- function(x0,y0,v0,l,ode_times){

  limiting_ode <- function(t, y, parms) {
    with(as.list(c(y, parms)), {
      dy1 <- l^2 * a(y[1],l)
      dy2 <- l^2 * a(y[2],l)
      dy3 <- l^2 * b_crn(y[1], y[2], y[3],l)
      return(list(c(dy1, dy2, dy3)))
    })
  }

  parms <- c(l = l)
  yini <- c(y1 = x0,
            y2 = y0,
            y3 = v0)

  out <- deSolve::ode(times = ode_times,
                      func = limiting_ode,
                      y = yini,
                      parms = parms)

  return(list("xsq" = out[, 2],
              "ysq" = out[, 3],
              "innerprod" = out[, 4]))
}



# Reflection coupling ###
b_refl <- function(x,y,v,l) { # ,z1,z,offset) {
  if(abs(x + y - 2*v) < 1e-12) {
    rho <- 1 
  } else {
    rho <- (2*x*y - (x+y)*v) / sqrt(x*y) / (x+y-2*v) #pmax(pmin(, 1), -1)
  }
g(x,y,rho,l) - 
  v * ( exp(0.5* l^2 *(x-1)) * pnorm(0.5 * l/sqrt(x) - l * sqrt(x)) + 
          exp(0.5* l^2 *(y-1)) * pnorm(0.5 * l/sqrt(y) - l * sqrt(y)) )
}

#' @export
ode_refl <- function(x0,y0,v0,l,ode_times) {
  
  limiting_ode <- function(t, y, parms) {
    with(as.list(c(y, parms)), {
      dy1 <- l^2 * a(y[1],l)
      dy2 <- l^2 * a(y[2],l)
      dy3 <- l^2 * b_refl(y[1], y[2], y[3],l)
      return(list(c(dy1, dy2, dy3)))
    })
  }
  
  parms <- c(l = l)
  yini <- c(y1 = x0, 
            y2 = y0, 
            y3 = v0)
  
  out <- deSolve::ode(times = ode_times, 
                      func = limiting_ode,
                      y = yini, 
                      parms = parms)
  
  return(list("xsq" = out[, 2],
              "ysq" = out[, 3],
              "innerprod" = out[, 4]))
}



# Asymptotically optimal Markovian coupling ###
b_optimal <- function(x,y,v,l) {
  
  pmin( pnorm(-0.5 * l / sqrt(x)) + exp(0.5* l^2 *(x-1)) * pnorm(0.5 * l/sqrt(x) - l * sqrt(x)),
        pnorm(-0.5 * l / sqrt(y)) + exp(0.5* l^2 *(y-1)) * pnorm(0.5 * l/sqrt(y) - l * sqrt(y)) ) - 
    v * ( exp(0.5* l^2 *(x-1)) * pnorm(0.5 * l/sqrt(x) - l * sqrt(x)) + 
            exp(0.5* l^2 *(y-1)) * pnorm(0.5 * l/sqrt(y) - l * sqrt(y)) )
}

#' @export
ode_optimal <- function(x0,y0,v0,l,ode_times) {
  
  limiting_ode <- function(t, y, parms) {
    with(as.list(c(y, parms)), {
      dy1 <- l^2 * a(y[1],l)
      dy2 <- l^2 * a(y[2],l)
      dy3 <- l^2 * b_optimal(y[1], y[2], y[3],l)
      return(list(c(dy1, dy2, dy3)))
    })
  }
  
  parms <- c(l = l)
  yini <- c(y1 = x0, 
            y2 = y0, 
            y3 = v0)
  
  out <- deSolve::ode(times = ode_times, 
                      func = limiting_ode,
                      y = yini, 
                      parms = parms)
  
  return(list("xsq" = out[, 2],
              "ysq" = out[, 3],
              "innerprod" = out[, 4]))
}



# ode_crn <- function(x0,y0,v0,l,ode_times,n_MC = 1e5){
# 
#   z1 <- rnorm(n_MC, sd = l)
#   z <- rnorm(n_MC, sd = l)
#   offset <- 0
# 
#   b_crn <- function(x,y,v,l,z1,z,offset) {
#     rho <- pmax(pmin(v / sqrt(x*y), 1), -1)
#     z2 <- rho * z1 + sqrt(1 - rho^2) * z
#     # MCavg2(rho, l, z1, z, offset)- # Called from C++
#     mean(exp(pmin(0, sqrt(x) * z1 - 0.5 * l^2, sqrt(y) * z2 - 0.5 * l^2))) +
#     offset -
#       v * ( exp(0.5* l^2 *(x-1)) * pnorm(0.5 * l/sqrt(x) - l * sqrt(x)) +
#               exp(0.5* l^2 *(y-1)) * pnorm(0.5 * l/sqrt(y) - l * sqrt(y)) )
#   }
# 
#   limiting_ode <- function(t, y, parms) {
#     with(as.list(c(y, parms)), {
#       dy1 <- l^2 * a(y[1],l)
#       dy2 <- l^2 * a(y[2],l)
#       dy3 <- l^2 * b_crn(y[1], y[2], y[3],l,z1,z,offset)
#       return(list(c(dy1, dy2, dy3)))
#     })
#   }
# 
#   parms <- c(l = l, z = z, z1 = z1, offset = offset)
#   yini <- c(y1 = x0,
#             y2 = y0,
#             y3 = v0)
# 
#   out <- deSolve::ode(times = ode_times,
#                       func = limiting_ode,
#                       y = yini,
#                       parms = parms, method = "euler")
# 
#   return(list("xsq" = out[, 2],
#               "ysq" = out[, 3],
#               "innerprod" = out[, 4]))
# 
# }


# ode_refl <- function(x0,y0,v0,l,ode_times,n_MC = 1e5){
#   
#   z1 <- rnorm(n_MC, sd = l)
#   z <- rnorm(n_MC, sd = l)
#   offset <- 2 * pnorm(-0.5 * l) - mean(exp(pmin(0, z1 - 0.5 * l^2, z1 - 0.5 * l^2)))
#   
#   b_refl <- function(x,y,v,l,z1,z,offset) {
#     if(abs(x + y - 2*v) < sqrt(.Machine$double.eps)) {
#       rho <- 1 
#     } else {
#       rho <- pmax(pmin((2*x*y - (x+y)*v) / sqrt(x*y) / (x+y-2*v), 1), -1)
#     }
#     # z2 <- rho * z1 + sqrt(1 - rho^2) * z
#     # mean(exp(pmin(0, z1 - 0.5 * l^2, z2 - 0.5 * l^2))) + offset
#     
#     MCavg2(rho, l1, z1, z, offset) - # Called from C++
#      v * ( exp(0.5* l^2 *(x-1)) * pnorm(0.5 * l/sqrt(x) - l * sqrt(x)) +
#             exp(0.5* l^2 *(y-1)) * pnorm(0.5 * l/sqrt(y) - l * sqrt(y)) )
#   }
#   
#   limiting_ode <- function(t, y, parms) {
#     with(as.list(c(y, parms)), {
#       dy1 <- l^2 * a(y[1],l)
#       dy2 <- l^2 * a(y[2],l)
#       dy3 <- l^2 * b_refl(y[1], y[2], y[3],l,z1,z,offset)
#       return(list(c(dy1, dy2, dy3)))
#     })
#   }
#   
#   parms <- c(l = l, z = z, z1 = z1, offset = offset)
#   yini <- c(y1 = x0, 
#             y2 = y0, 
#             y3 = v0)
#   
#   out <- deSolve::ode(times = ode_times, 
#                       func = limiting_ode,
#                       y = yini, 
#                       parms = parms)
#   
#   return(list("xsq" = out[, 2],
#               "ysq" = out[, 3],
#               "innerprod" = out[, 4]))
# }




# x <- 1
# y <- 1
# l <- 1
# rhos <- seq(-0.999,0.9999, 0.0001)
# gs <- rep(NA, length(rhos))
# 
# for(i in seq_along(rhos)) {
#   gs[i] <- g(x,y,rhos[i],l)[1]
# }
# plot(gs)
# abline(h = 2 * pnorm(-l/2), col = "red")
# 
# 
# l1sq <- l^2
# MCavg <- function(x,y,rho,n) {
#   z1 <- rnorm(n)
#   z <- rnorm(n)
# 
#   out <- rep(NA, length(rho))
#   for (r in seq_along(rho)){
#     z2 <- rho[r] * z1 + sqrt(1 - rho[r]^2) * z
#     out[r] <- mean(exp(pmin(0, sqrt(l1sq*x)*z1-0.5 * l1sq, sqrt(l1sq*y)*z2-0.5 * l1sq)))
#   }
#   # Correct by forcing LHS to take the analytically correct value when the correlation is 1
#   offset <- 0
#   out <- out + offset
# 
#   return(out)
# }
# 
# lines(MCavg(x,y,rhos, 1e4), col = "red")


#####################
#####
# Check that the Owen table is correct
#####
#

# owen_10010.1 <- function(x,a,b){
#   Sigma <- diag(2); Sigma[1,2] <- Sigma[2,1] <- (-b)/sqrt(1+b^2)
#   
#   xs <- rep(NA, length(x))
#   for(i in seq_along(xs)) {
#     xs[i] <- mvtnorm::pmvnorm(upper = c(a/sqrt(1+b^2), x[i]), corr = Sigma)
#   }
#   return(xs)
# }
# 
# owen_10010.1_MC  <- function(x,a,b,n_MC){
#   z <- rnorm(n_MC)
#   
#   xs <- rep(NA, length(x))
#   for(i in seq_along(xs)) {
#     xs[i] <- mean( (z < x[i]) * pnorm(a+b*z))
#   }
#  return(xs)
# }

# a <- 1
# b <- 1
# n_MC <- 1e6
# 
# x <- seq(-2,2,0.1)
# 
# analytical <- owen_10010.1(x,a,b)
# numerical  <- owen_10010.1_MC(x,a,b,n_MC)
# 
# analytical
# numerical
# plot(numerical - analytical)

#####################

# mean(exp(pmin(0, l * rnorm(1e6) - 0.5 * l^2, l * rnorm(1e6) - 0.5 * l^2)))
# 
# l <- 2.38
# x <- 1
# y <- 1
# rho <- 0
# 
# h(x,y,rho,l)[1]
# 
# n <- 1e5
# z1 <- rnorm(n)
# z2 <- rho * z1 + sqrt(1- rho^2) * rnorm(n)
# mean((z1 < l/(2* sqrt(x))) * (sqrt(x) *z1 < sqrt(y) *z2) * exp(l*sqrt(x) *z1 - l^2/2))
# 
# z_ <- rnorm(1e4)
# pnorm(-(sqrt(x/y) - rho)/sqrt(1 - rho^2) * z1[1])
# mean( (sqrt(x) *z1[1] < sqrt(y) * (rho * z1[1] + sqrt(1-rho^2)* z_)))

