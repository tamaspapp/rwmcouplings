
library(ggplot2)
# library(mvtnorm)
# # Auxiliary functions ####
# a <- function(x,l) {
#   pnorm(-0.5 * l / sqrt(x)) + (1-2*x) * exp(0.5 * l^2 *(x-1)) * pnorm(0.5 * l/sqrt(x) - l * sqrt(x))
# }
# 
# g <- function(x,y,rho,l) {
#   # Deal with rho very close to 1 or very close to -1. Otherwise, the calculations seem numerically stable.
#   if(rho > 1-1e-12) {
#     return(2*pnorm(-0.5*l))
#   } else if(rho < -1+1e-12) {
#     return(g(x,y,-1+1e-12,l))
#   }
#   
#   # Owen (1980)'s (10,011.1)
#   owen_10010.1 <- function(Y,a,b){
#     Sigma <- diag(2); Sigma[1,2] <- Sigma[2,1] <- (-b)/sqrt(1+b^2)
#     mvtnorm::pmvnorm(upper = c(a/sqrt(1+b^2), Y), corr = Sigma)
#   }
#   
#   h <- function(x,y,rho,l) {
#     b <- -(sqrt(x/y) - rho)/sqrt(1 - rho^2)
#     a <- b * l * sqrt(x)
#     U <- 0.5 * l / sqrt(x) - l * sqrt(x)
#     exp(0.5*l^2*(x-1)) * owen_10010.1(U, a, b)
#   }
#   
#   Sigma <- diag(2); Sigma[1,2] <- Sigma[2,1] <- rho
#   lb <- 0.5*l/sqrt(c(x,y))
#   mvtnorm::pmvnorm(lower = lb, corr = Sigma) + h(x,y,rho,l) + h(y,x,rho,l) 
# }
# 
# # Reflection coupling ###
# b_refl <- function(x,y,v,l) { 
#   if(abs(x + y - 2*v) < 1e-12) {
#     rho <- 1 
#   } else {
#     rho <- (2*x*y - (x+y)*v) / sqrt(x*y) / (x+y-2*v) #pmax(pmin(, 1), -1)
#   }
#   g(x,y,rho,l) - 
#     v * ( exp(0.5* l^2 *(x-1)) * pnorm(0.5 * l/sqrt(x) - l * sqrt(x)) + 
#             exp(0.5* l^2 *(y-1)) * pnorm(0.5 * l/sqrt(y) - l * sqrt(y)) )
# }
# 
# ####
# # Heatmap of optimal step size, varying position of one chain relative to a stationary chain
# #### 
# drift_refl <- function(x,y,rho,l) {
#   v <- sqrt(x*y)*rho
#   l^2 * (a(x,l) + a(y,l) - 2 * b_refl(x,y,v,l))
# }
# drift_refl_vec <- Vectorize(drift_refl,"l")

res <- 5e-2 # Resolution of grid

x <- seq(res,2,res)
y <- 1
rho <- seq(-1,1-res,res)

l_heatmap <- expand.grid(x,rho)

ls <- seq(1e-2,5,1e-2)

library(doParallel)

# Set up parallel processing
ncores <- 7 # Change as needed!
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

xs_grid <- l_heatmap[,1]
rhos_grid <- l_heatmap[,2]

l_refl <- foreach(x_ = xs_grid, rho_ = rhos_grid, 
        .combine = "c", .packages = "mvtnorm") %dopar% {
          
          
          # Auxiliary functions ####
          a <- function(x,l) {
            pnorm(-0.5 * l / sqrt(x)) + (1-2*x) * exp(0.5 * l^2 *(x-1)) * pnorm(0.5 * l/sqrt(x) - l * sqrt(x))
          }
          
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
          
          # Reflection coupling ###
          b_refl <- function(x,y,v,l) { 
            if(abs(x + y - 2*v) < 1e-12) {
              rho <- 1 
            } else {
              rho <- (2*x*y - (x+y)*v) / sqrt(x*y) / (x+y-2*v) #pmax(pmin(, 1), -1)
            }
            g(x,y,rho,l) - 
              v * ( exp(0.5* l^2 *(x-1)) * pnorm(0.5 * l/sqrt(x) - l * sqrt(x)) + 
                      exp(0.5* l^2 *(y-1)) * pnorm(0.5 * l/sqrt(y) - l * sqrt(y)) )
          }
          
          ####
          # Heatmap of optimal step size, varying position of one chain relative to a stationary chain
          #### 
          drift_refl <- function(x,y,rho,l) {
            v <- sqrt(x*y)*rho
            l^2 * (a(x,l) + a(y,l) - 2 * b_refl(x,y,v,l))
          }
          drift_refl_vec <- Vectorize(drift_refl,"l")
          
  fun <- function(l){drift_refl_vec(x_,1,rho_,l)}
  vals <- fun(ls)
  ls[which.min(vals)]
}

parallel::stopCluster(cl)

l_heatmap <- as.data.frame(cbind(l_heatmap, l_refl))
names(l_heatmap) <- c("y","rho","l")


ggplot(l_heatmap, aes(y, rho, fill = l)) +
  geom_raster(interpolate = T) +
  scale_fill_gradient(low = "blue",high="yellow") +
  scale_x_continuous(expand=c(0,0), limits = c(0,2)) +
  scale_y_continuous(expand=c(0,0), limits = c(-1,1))


  