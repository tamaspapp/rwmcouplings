
library(rwmcouplings)
library(Matrix)
library(doParallel)
library(doRNG)
library(reshape2)
library(ggplot2)
library(latex2exp)
library(ggpubr)
library(reshape2)

ncores <- 4

################################################################################
# Hug and Hop
################################################################################
# Load in model data 
if(!file.exists("svmData.RData")) {source("svm-prelim.R")} # Ensure we run this first
load(file = "svmData.RData")

seed <- 12345; set.seed(seed); SetSeed_cpp(seed) # Seed RNGs

B    <- 10 # Hug and HMC params
T_hh <- 0.35
T_hmc <- 0.225

lam <- 20 # Hop params
kap <- 1
thresh_hop <- 1e-5

rwm_proportion <- 0.05 # RWM params
l <- 2.38/sqrt(avg_prec)
h <- l / sqrt(t)
thresh_rwm <- 0.1

# # Acceptance rate of Hug and Hop
# temp1 <- hughop(SampleLatentVariables(t, sig, phi),T_hh,B,lam,kap,1e4,logpi,gradlogpi)
# temp1$acc_hug
# temp1$acc_hop
# 
# # Acceptance rates of HMC
# temp2 <- rwmhmc(SampleLatentVariables(t, sig, phi),T_hmc,B,rwm_proportion,h,1e4,logpi,gradlogpi)
# temp2$acc_x_rwm
# temp2$acc_x_hmc
# ###

# Parallel MCMC setup
R_hh <- 100 # Replicates
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

x0 <- y0 <- vector("list", R_hh)
for(i in 1:R_hh) { x0[[i]] <- SampleLatentVariables(t, sig, phi) }
for(i in 1:R_hh) { y0[[i]] <- SampleLatentVariables(t, sig, phi) }


# First L iterations
L <- 2.5e3 # Lag parameter

# HH
xL_hh <- foreach(x0_ = x0, .packages = "rwmcouplings") %dorng% {
  firstL <- hughop(x0_,T_hh,B,lam,kap,L,logpi,gradlogpi)
  firstL$x
}
# HMC
xL_hmc <- foreach(x0_ = x0, .packages = "rwmcouplings") %dorng% {
  firstL <- rwmhmc(x0_,T_hmc,B,rwm_proportion,h,L,logpi,gradlogpi)
  firstL$x
}

# Coupled MCMC iterations
iter <- 1e5 # Cap at 100,000 iterations

hh <- 
  foreach(xL_ = xL_hh, y0_ = y0,.packages = "rwmcouplings",
          .combine='comb',.multicombine=TRUE,.init=list(list(), list())) %dorng% {
  out <- cplhughop(xL_,y0_,T_hh,B,lam,kap,thresh_hop,iter,logpi,gradlogpi)
  list(out$squaredist, out$tau)
  }; names(hh) <- c("squaredist", "tau")


hmc <- 
  foreach(xL_ = xL_hmc, y0_ = y0,.packages = "rwmcouplings",
          .combine='comb',.multicombine=TRUE,.init=list(list(), list())) %dorng% {
    out <- twoscalehmc(xL_,y0_,T_hmc,B,rwm_proportion,h,thresh_rwm,iter,logpi,gradlogpi)
    list(out$squaredist, out$tau)
  }; names(hmc) <- c("squaredist", "tau")

parallel::stopCluster(cl) # Stop cluster
save(hh, hmc, file = "hhout.RData")


#####
# Plotting ###
load(file = "hhout.RData")

iter_plot <- 2.5e3
iters <- seq(0, iter_plot)

# Squared Wasserstein distance
# NB: assumes all chains have met before time L, which is the case with our setup.
w2sq_df <- data.frame("hh" = Reduce("+", hh$squaredist)[0:iter_plot + 1] / R_hh,
                      "hmc" = Reduce("+", hmc$squaredist)[0:iter_plot + 1] / R_hh,
                      "iter" = iters)
w2sq_df <- melt(w2sq_df, id.vars = "iter")


# Define custom transformation for y-axis: log(0.1 + y)
log.const.p <- function(const = 0.1, name = "log.const.p") {
  scales::trans_new(name, 
                    function(x) {log10(x + const)},
                    function(y) {10^(y) - const},
                    domain = c(-const,Inf))
}
log0.1p <- log.const.p(0.1)
  
hh_w2sq <-
  ggplot(w2sq_df, aes(x = iter, y = value, color = variable)) +
  geom_hline(yintercept = 0,linetype = 2) +
  geom_line() +
  scale_y_continuous(trans = log0.1p, breaks = c(0,0.1,1,10,100, 500), minor_breaks = NULL) +
  #coord_trans(y = log0.1p) + 
  #scale_y_continuous(breaks = c(0,1,10,100, 500), minor_breaks = NULL) +
  labs(x = "Iteration", y = "Squared Wasserstein distance", color = "Algorithm") +
  scale_color_discrete(labels = list("hmc" = "HMC", "hh" = "Hug & Hop")) +
  theme_bw()
hh_w2sq

# Total variation distance
tvd_df <- data.frame("hh" = getTVDBound(iters, unlist(hh$tau)),
                     "hmc" = getTVDBound(iters, unlist(hmc$tau)), 
                     "iter" = iters)
tvd_df <- melt(tvd_df, id.vars = "iter")  

hh_tvd <-
  ggplot(tvd_df, aes(x = iter, y = value, color = variable)) +
  geom_hline(yintercept = 0,linetype = 2) +
  geom_line() +
  labs(x = "Iteration", y = "Total variation distance", color = "Algorithm") +
  scale_color_discrete(labels = list("hmc" = "HMC", "hh" = "Hug & Hop")) +
  theme_bw()
hh_tvd

# Join plots into one, with common legend ###
comb_plot <- ggarrange(hh_w2sq, hh_tvd, 
                       ncol=2, nrow=1, common.legend = T, legend = "right")
comb_plot

ggsave(filename = "hh_convergence.pdf",
       plot = comb_plot,
       device = "pdf",  width = 24, height = 8, units = "cm", bg = "transparent")




