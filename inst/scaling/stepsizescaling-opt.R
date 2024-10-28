library(latex2exp)
library(ggplot2)
library(ggpubr)
library(Cairo)

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
  q <- function(x,l){l^2 * exp(0.5* l^2 *(x-1)) * pnorm(0.5 * l/sqrt(x) - l * sqrt(x))}
  
  g(x,y,v,l) - v * (q(x,l) + q(y,l))
}
b_opt <- function(x,y,v,l) {
  p <- function(x,l) {l^2 * pnorm(-0.5 * l / sqrt(x)) + l^2 * exp(0.5 * l^2 *(x-1)) * pnorm(0.5 * l/sqrt(x) - l * sqrt(x))}
  q <- function(x,l) {l^2 * exp(0.5 * l^2 *(x-1)) * pnorm(0.5 * l/sqrt(x) - l * sqrt(x))}
  
  pmin(p(x,l),p(y,l)) - v * (q(x,l) + q(y,l))
} 

drift_gcrn <- function(x,y,rho,l) {
  v <- sqrt(x*y)*rho
  a(x,l) + a(y,l) - 2 * b_gcrn(x,y,v,l)
}
drift_opt <- function(x,y,rho,l) {
  v <- sqrt(x*y)*rho
  a(x,l) + a(y,l) - 2 * b_opt(x,y,v,l)
}



res1 <- 1e-2

x <- seq(res1,2,res1)
y <- 1
rho <- seq(-1,1-res1,res1)

l_heatmap <- expand.grid(x,rho)

ls <- seq(1e-2,5,1e-2)

l_opt <- rep(NA, nrow(l_heatmap))
drifts <- rep(NA, nrow(l_heatmap))
  
for(j in 1:nrow(l_heatmap)){
  x <- l_heatmap[j,1]
  rho <- l_heatmap[j,2]
  fun <- function(l){drift_opt(x,y,rho,l)}
  
  l_opt[j] <- ls[which.min(fun(ls))]
  drifts[j] <- fun(l_opt[j])
}

plot(-drifts, log ="y")


l_heatmap <- as.data.frame(cbind(l_heatmap, l_opt))
names(l_heatmap) <- c("y","rho","l")


plot_params <- list(
  geom_raster(interpolate = T),
  scale_fill_gradient(low = "blue",high="yellow", limits = c(NA,NA)),
  scale_x_continuous(expand=c(0,0), limits = c(0,2)),
  scale_y_continuous(expand=c(0,0), limits = c(0,1)), # Limit to rho in [0,1] here
  theme_bw(),
  theme(axis.title.y = element_text(size = 14),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))
)

opt <-
  ggplot(l_heatmap, aes(y, rho, fill = l)) +
    plot_params +
    scale_fill_gradient(low = "blue",high="yellow", limits = c(NA,2.6)) +
    labs(x = TeX("$y$"), y = TeX("$\\rho$"), 
       fill = TeX(paste0("$","\u2113","_{opt}(y, \\rho)$")))
opt


####
# Relative efficiency of GCRN vs optimal Markovian
#### 

# 1. Create the grid of values ###
y <- 1

# Reparametrize in terms of cosine similarity
res1 <- 1e-2 # resolution of grid
x <- seq(res1,2,res1)
rho <- seq(-1,1-res1,res1)
xrho <- expand.grid(x,rho)

s_fun <- function(x,y,rho){x + y - 2 * sqrt(x*y) * rho} 
s <- sapply(1:nrow(xrho),function(i){s_fun(xrho[i,1],y,xrho[i,2])})

# 2. Search for the optimum over this grid ###
ls <- seq(res1,5,res1)
drifts_gcrn <- drifts_opt <- rep(NA, nrow(xrho))

for(j in 1:nrow(xrho)){
  x <- xrho[j,1]
  rho <- xrho[j,2]
  
  drifts <- drift_gcrn(x,y,rho,ls)
  idx_l_gcrn <- which.min(drifts)
  drifts_gcrn[j] <- drifts[idx_l_gcrn]
  
  fun <- function(l){drift_opt(x,y,rho,l)}
  idx_l_opt <- which.min(fun(ls))
  drifts_opt[j] <- fun(ls[idx_l_opt])
}

eff_heatmap <- as.data.frame(cbind(xrho, drifts_opt, drifts_gcrn))
names(eff_heatmap) <- c("y","rho","opt","gcrn")

eff <- 
  ggplot(eff_heatmap, aes(y, rho, fill =  gcrn/opt)) +
  plot_params +
  labs(x = TeX("$y$"), y = TeX("$\\rho$"), 
       fill = "Relative efficiency")
eff

comb_plot <- ggarrange(opt, 
                       eff, 
                       ncol=2, nrow=1, common.legend = FALSE)
comb_plot

ggsave(filename = "opt-vs-gcrn.pdf",
       plot = comb_plot, device=cairo_pdf, width = 30, height = 12, units = "cm", bg = "transparent")




#

q <- function(x,l = 1){l^2 * exp(0.5* l^2 *(x-1)) * pnorm(0.5 * l/sqrt(x) - l * sqrt(x))}

# plot(x, q(x), type = "l", main = "Function q(x)")
# abline(v = 1, lty = 2)

x <- seq(0, 20, 0.1)
l <- c(0.001, 1, 2.38, 10)

for(l_ in l) plot(x, sqrt(x) * q(x, l_), type = "l", main = paste0("Function sqrt(x) * q(x,",l_, ")"))
# Numerical error for large (l,x) since we don't work on the log scale

