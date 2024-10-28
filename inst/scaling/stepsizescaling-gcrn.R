library(latex2exp)
library(ggplot2)
library(ggpubr)
library(Cairo)
library(doParallel)

# Auxiliary functions related to the drift ####
a <- function(x,l) {
  l^2 * pnorm(-0.5 * l / sqrt(x)) + 
    l^2 *(1-2*x) * exp(0.5 * l^2 *(x-1)) * pnorm(0.5 * l/sqrt(x) - l * sqrt(x))
}
b_gcrn <- function(x,y,v,l) {
  g <- function(x,y,v,l) {
    m <- pmin(x,y)
    M <- pmax(x,y)
    
    l^2 * pnorm(-0.5 * l / sqrt(m)) +  
      l^2 * exp(0.5* l^2 *(M-1)) * pnorm(- l * sqrt(M)) +
      l^2 * exp(0.5* l^2 *(m-1)) * ( pnorm(0.5 * l/sqrt(m) - l * sqrt(m)) - pnorm(- l * sqrt(m)) )
  }
  p <- function(x,l){l^2 * exp(0.5* l^2 *(x-1)) * pnorm(0.5 * l/sqrt(x) - l * sqrt(x))}
  
  g(x,y,v,l) - v * (p(x,l) + p(y,l)) 
}

# Drift of GCRN coupling, parametrized in terms of cosine similarity rho (= correlation)
drift_gcrn <- function(x,y,rho,l) {
  v <- sqrt(x*y)*rho
  a(x,l) + a(y,l) - 2 * b_gcrn(x,y,v,l)
}

####
# Heatmap of absolute value of relative drift 
# |a-relative(x,l)| against x and l 
#### 

res <- 1e-2 # resolution of grid
l <- seq(res, 4, res)
x <- c(1e-9 + seq(0,2,res))

# 1. Get optimal step size "l" and maximum of absolute value of drift "a" ###
l_opt <- a_abs_max <- rep(NA,length(x))
for(i in seq_along(x)) {
  opt <- optim(2,function(l){-abs(a(x[i],l))},method = "Brent", lower = 0.5, upper = 4)
  l_opt[i]  <- opt$par
  a_abs_max[i] <- -opt$value
}
max_df <- as.data.frame(cbind(x, l_opt)); names(max_df) <- c("x", "l")
max_df <- max_df[max_df$x != 1,]

# 1.1. How far off is the stationary-optimal step size? ### 
min(abs(a(x, 2.38)) / a_abs_max)
# 1.2. How far off is the mode-optimal step size? ### 
min(abs(a(x, sqrt(2))) / a_abs_max)

# 2. Calculate relative drift ###
a_abs <- abs(outer(x,l, FUN = a))
a_rel <- diag(1/a_abs_max) %*% a_abs
a_heatmap <- cbind(expand.grid(x,l), as.vector(a_rel))
names(a_heatmap) <- c("x", "l", "a")

# 3. Plot the heatmap and optimal step sizes ###
convergence_scaling <- 
ggplot(a_heatmap, aes(x, l, fill = a)) +
  geom_raster(interpolate = T) + 
  scale_fill_gradient(low = "blue",high="yellow", limits = c(0,1)) +
  geom_line(data=max_df, mapping=aes(x, l, fill = NULL), color = "black", lwd = 0.6, lty=2) +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  labs(x = TeX("$x$"), y = "\u2113", fill = TeX(paste0("$a_{rel}(x$, ","\u2113","$)$"))) + # \u2113 == "\ell"
  theme_bw() +
  theme(axis.title.y = element_text(size = 14),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))
convergence_scaling

acc_rate <- function(x,l) {pnorm(-0.5 * l / sqrt(x)) + exp(0.5 * l^2 *(x-1)) * pnorm(0.5 * l/sqrt(x) - l * sqrt(x))}

# 4. What is the acceptance rate pointwise at the pointwise optimal step size?
quantile(acc_rate(x,l_opt),c(0,0.1,0.5,0.9,1))

# 5. What is the acceptance rate pointwise at a step size of 2.38?
quantile(acc_rate(x,2.38),c(0,0.1,0.5,0.9,1))

####
# Heatmap of optimal step size, varying position of one chain relative to a stationary chain
#### 

# 1. Create the grid of values ###
y <- 1

# Reparametrize in terms of cosine similarity
res1 <- 1e-2 # resolution of grid
x <- seq(res1,2,res1)
rho <- seq(0,1-res1,res1) # rho <- seq(-1,1-res1,res1)
xrho <- expand.grid(x,rho)

xs <- xrho[, 1]
rhos <- xrho[, 2]

s_fun <- function(x,y,rho){x + y - 2 * sqrt(x*y) * rho} 
s <- sapply(1:nrow(xrho),function(i){s_fun(xrho[i,1],y,xrho[i,2])})

# 2. Search for the optimum over this grid ###
ls <- seq(res1,4,res1)

ncores <- 4
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

l_heatmap <- 
  foreach(x_ = xs, rho_ = rhos, .combine = "rbind", .multicombine = T) %dopar% {
      drifts <- drift_gcrn(x_,y,rho_,ls)
  
      idx_optim <- which.min(drifts)
      drift_optim <- drifts[idx_optim]
  
      # If the optimal drift is non-negative, it shows up as undefined in the heatmap
      l_gcrn <- NA
      if(drift_optim < 0) {
        l_gcrn <- ls[idx_optim]
      }
      
      # 2.1 Also check what interval we're 80% optimal in
      l_left <- NA
      l_right <- NA
      if(drift_optim < 0) {
        drift_optim  <- drifts[idx_optim]
        drifts_left  <- drifts[1:idx_optim]
        drifts_right <- drifts[idx_optim:length(drifts)]

        idx_left  <- base::Position(function(x){x < 0.8 * drift_optim}, drifts_left)
        idx_right <- base::Position(function(x){x >= 0.8 * drift_optim}, drifts_right)
        l_left  <- ls[idx_left]
        l_right <- ls[idx_right + idx_optim - 1]
      }
    return(c(l_gcrn, l_left, l_right))
  }

parallel::stopCluster(cl)

l_heatmap <- as.data.frame(cbind(xrho, l_heatmap))
names(l_heatmap) <- c("y","rho","l", "l_left, l_right")

gcrn_plot_params <- list(
  geom_raster(interpolate = F),
  scale_fill_gradient(low = "blue",high="yellow", limits = c(NA,2.5)),
  scale_x_continuous(expand=c(0,0), limits = c(0,2)),
  scale_y_continuous(expand=c(0,0), limits = c(0,1)), # Restrict to rho \in [0,1]
  theme_bw(),
  theme(axis.title.y = element_text(size = 14),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))
  )

gcrn_scaling <- 
  ggplot(l_heatmap, aes(y, rho, fill = l)) +
  gcrn_plot_params +
  labs(x = TeX("$y$"), y = TeX("$\\rho$"), 
       fill = TeX(paste0("$","\u2113","_{gcrn}(y, \\rho)$")))
gcrn_scaling

# # How insensitive is the coupling to the choice of step size (within 80% of optimum)?
# ggplot(l_heatmap, aes(y, rho, fill = l_left)) +
#   gcrn_plot_params +
#   labs(x = TeX("$y$"), y = TeX("$\\rho$"), fill = TeX("$l_{left}$"))
# 
# ggplot(l_heatmap, aes(y, rho, fill = l_right)) +
#   gcrn_plot_params +
#   labs(x = TeX("$y$"), y = TeX("$\\rho$"), fill = TeX("$l_{right}$"))


comb_plot <- ggarrange(convergence_scaling, 
                       gcrn_scaling, 
                       ncol=2, nrow=1, common.legend = FALSE)
comb_plot

ggsave(filename = "scaling.pdf",
       plot = comb_plot, device=cairo_pdf, width = 30, height = 12, units = "cm", bg = "transparent")








# l_gcrn <- l_left <- l_right <- rep(NA, nrow(xrho))
# for(j in 1:nrow(xrho)){
#   x <- xrho[j,1]
#   rho <- xrho[j,2]
#   drifts <- drift_gcrn(x,y,rho,ls)
#   
#   idx_optim <- which.min(drifts)
#   drift_optim <- drifts[idx_optim]
#   
#   # print(drift_optim)
#   # If the optimal drift is non-negative, it shows up as undefined in the heatmap
#   if(drift_optim < 0) { 
#     l_gcrn[j] <- ls[idx_optim]
#   }
#   
#   # # 2.1 Also check what interval we're 80% optimal in
#   # if(drift_optim < 0) {
#   #   drift_optim  <- drifts[idx_optim]
#   #   drifts_left  <- drifts[1:idx_optim]
#   #   drifts_right <- drifts[idx_optim:length(drifts)]
#   #   
#   #   idx_left  <- base::Position(function(x){x < 0.8 * drift_optim}, drifts_left)
#   #   idx_right <- base::Position(function(x){x >= 0.8 * drift_optim}, drifts_right)
#   #   l_left[j]  <- ls[idx_left]
#   #   l_right[j] <- ls[idx_right + idx_optim - 1]
#   # }
#   # # if(j == 1){plot(drifts, type = "l", ylim = c(-4, 1))} else{lines(drifts, col = j)} # Looks unimodal
# }

  