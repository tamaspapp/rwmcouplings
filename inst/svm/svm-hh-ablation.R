
library(rwmcouplings)
library(Matrix)
library(doParallel)
library(doRNG)
library(reshape2)
library(ggplot2)
library(latex2exp)
library(egg)

crnhmc <- function(x0,y0,T_hmc,L_hmc,iter,gradlogpi) {
  # The zeros below mean that we don't mix in the RWM at all
  twoscalehmc(x0,y0,T_hmc,L_hmc,0,0,0,iter,logpi,gradlogpi)
}

gcrnhh <- function(x0,y0,T_hug,B_hug,lam,kap,iter,gradlogpi) {
  # The -1 below means that only GCRN coupling is used for Hop
  cplhughop(x0,y0,T_hug,B_hug,lam,kap,-1,iter,logpi,gradlogpi)
}


# Load in model data 
if(!file.exists("svmData.RData")) {source("svm-prelim.R") } # Ensure we run this first
load(file = "svmData.RData")

# Hop parameters, fixed a priori
lam <- 20
kap <- 1

# Parallel setup
ncores <- 6

########
# 1. Contractivity, varying the integration time and fixing a small step size 
########
# - There is a sharp phase transition for the contractivity of HMC:
#   we need a short enough integration time.
# - Hug does not have this issue.

# Replicates and iterations / parameter setting
R <- 5
iter <- 1e3

# HMC and Hug parameters
eps <- 1e-3
Ts  <- seq(0.1,0.5,0.05)

# Set RNG
seed <- 12345; set.seed(seed); SetSeed_cpp(seed)
rng <- RNGseq(R, seed) # For nested foreach "doRNG" loops

# Start coupled chains from independent draws from the Laplace approximation
x0 <- y0 <- vector(mode = "list", R)
for(i in 1:R) { x0[[i]] <- mu + as.vector(Sigma_cholLT%*%rnorm(t)) }
for(i in 1:R) { y0[[i]] <- mu + as.vector(Sigma_cholLT%*%rnorm(t)) }

cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

hmc_out <- 
  foreach(T_  = Ts, .combine = "rbind") %:%     
  foreach(i=1:R, x0_ = x0, y0_ = y0, r = rng, # Need i=1:R to count replicates
          .combine='rbind',.multicombine=T,.packages = "rwmcouplings") %dopar% {
    rngtools::setRNG(r)
    B <- as.integer(T_/eps)
    cbind(crnhmc(x0_,y0_,T_,B,iter,gradlogpi)$squaredist,
          0:iter,T_,i)
  }
hmc_out <- as.data.frame(hmc_out)
names(hmc_out) <- c("squaredist", "iteration", "T", "replicate")

hh_out <- 
  foreach(T_  = Ts, .combine = "rbind") %:%     
  foreach(i=1:R, x0_ = x0, y0_ = y0, r = rng,
          .combine='rbind',.multicombine=T,.packages = "rwmcouplings") %dopar% {
    rngtools::setRNG(r)
    B <- as.integer(T_/eps)
    cbind(gcrnhh(x0_,y0_,T_,B,lam,kap,iter,gradlogpi)$squaredist,
          0:iter,T_,i)
  }
hh_out <- as.data.frame(hh_out)
names(hh_out) <- c("squaredist", "iteration", "T", "replicate")

parallel::stopCluster(cl)
save(hmc_out,hh_out,file = "hhIntegrationTime.RData")

# Plotting
load(file = "hhIntegrationTime.RData")

hmc_inttime <- 
ggplot(hmc_out, 
       aes(y = squaredist, x = iteration, 
           color = T, group = interaction(replicate,T)))+ 
  geom_line(alpha = 0.7) +
  scale_y_log10(breaks = scales::breaks_log(7)) +
  coord_cartesian(ylim = c(1e-12,1e2)) +
  scale_color_gradient2(limits = c(min(Ts),max(Ts)),low = "red", mid = "yellow", high = "blue") +
  labs(x = "Iteration", 
       y = "Squared distance", 
       title = "HMC") +
  theme_bw()
# hmc_inttime

hh_inttime <-
ggplot(hh_out, 
       aes(y = squaredist, x = iteration, 
           color = T, group = interaction(T,replicate)))+ 
  geom_line(alpha = 0.7) +
  scale_y_log10(breaks = scales::breaks_log(7)) +
  coord_cartesian(ylim = c(1e-12,1e2)) +
  scale_color_gradient2(limits = c(min(Ts),max(Ts)),low = "red", mid = "yellow", high = "blue") +
  labs(x = "Iteration", 
       y = "Squared distance",
       title = "Hug and Hop") +
  theme_bw()
# hh_inttime

# # Avg. squared distance after 1000 iterations
# aggregate(squaredist ~ T, hh_out[hh_out$iteration == 1000, ], mean)
# aggregate(squaredist ~ T, hmc_out[hmc_out$iteration == 1000, ], mean)

inttime <- egg::ggarrange(hmc_inttime + theme(legend.position = "none"), 
                          hh_inttime + theme(axis.text.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.title.y = element_blank() ), 
                     ncol=2, nrow=1)
inttime

ggsave(filename = "hhIntegrationTime.pdf",
       plot = inttime,
       device = "pdf",  width = 24, height = 8, units = "cm", bg = "transparent")




########
# 2. Search for parameter settings
########

# Replicates and iterations / parameter setting
R <- 20
iter <- 2e3

# HMC and Hug parameters
B <- c(10,20,30)
T_hmc  <- seq(0.15,0.25,0.025)
T_hug  <- seq(0.2,0.5,0.05)

# Set RNG
seed <- 12345; set.seed(seed); SetSeed_cpp(seed)
rng <- RNGseq(R, seed) # For nested foreach "doRNG" loops

# Start coupled chains from independent draws from the Laplace approximation
x0 <- y0 <- vector(mode = "list", R)
for(i in 1:R) { x0[[i]] <- mu + as.vector(Sigma_cholLT%*%rnorm(t)) }
for(i in 1:R) { y0[[i]] <- mu + as.vector(Sigma_cholLT%*%rnorm(t)) }

cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

hmc_out <- 
  foreach(T_  = T_hmc, .combine = "rbind") %:%
  foreach(B_  = B, .combine = "rbind") %:%     
  foreach(i=1:R, x0_ = x0, y0_ = y0, r = rng,
          .combine='rbind',.multicombine=T,.packages = "rwmcouplings") %dopar% {
            rngtools::setRNG(r)
            cbind(crnhmc(x0_,y0_,T_,B_,iter,gradlogpi)$squaredist,
                  0:iter,T_,B_,i)
          }
hmc_out <- as.data.frame(hmc_out)
names(hmc_out) <- c("squaredist", "iteration", "T", "B", "replicate")

hh_out <- 
  foreach(T_  = T_hug, .combine = "rbind") %:%
  foreach(B_  = B, .combine = "rbind") %:%     
  foreach(i=1:R, x0_ = x0, y0_ = y0, r = rng,
          .combine='rbind',.multicombine=T,.packages = "rwmcouplings") %dopar% {
            rngtools::setRNG(r)
            cbind(gcrnhh(x0_,y0_,T_,B_,lam,kap,iter,gradlogpi)$squaredist,
                  0:iter,T_,B_,i)
          }
hh_out <- as.data.frame(hh_out)
names(hh_out) <- c("squaredist", "iteration", "T", "B", "replicate")

parallel::stopCluster(cl)
save(hmc_out,hh_out,file = "hhParamTuning.RData")

# Plotting
load(file = "hhParamTuning.RData")

# Plot using facets
hmc_out$algorithm <- "HMC"
hh_out$algorithm <- "Hug and Hop"
paramtuning_out <- rbind(hmc_out, hh_out)
paramtuning_out <- paramtuning_out[paramtuning_out$iteration == c(1000,2000),]


prepender <- function(string, prefix = "Iteration: ") paste0(prefix, string)
iter_labels <- as_labeller(prepender)

paramtuning <- 
  ggplot(paramtuning_out,
         aes(x = as.factor(T), color = as.factor(B), y = squaredist))+
    geom_boxplot() +
    facet_grid(vars(iteration), vars(algorithm), scales = "free", labeller = labeller(iteration = iter_labels)) +
    scale_y_log10(breaks = scales::breaks_log(9)) +
    labs(x = "T (Integration time)",
         y = "Squared distance",
         color = "B") +
    theme_bw()
paramtuning

ggsave(filename = "hhParamTuning.pdf",
       plot = paramtuning,
       device = "pdf",  width = 24, height = 16, units = "cm", bg = "transparent")


########
# 2.1. Acceptance rates at selected parameter tunings
########

B_final     <- 10
T_hh_final  <- 0.35
T_hmc_final <- 0.225

hh_accrate <- hughop(x0[[1]], T_hh_final, B_final, lam, kap, 5e4, logpi, gradlogpi)
hh_accrate$acc_hug
hh_accrate$acc_hop

hmc_accrate <- rwmhmc(x0[[1]], T_hh_final, B_final, 0, 0, 5e4, logpi, gradlogpi)
hmc_accrate$acc_x_hmc


########
# 3. Choice of threshold delta in two-scale Hop coupling
########
# Data wrangling for threshold plots
threshDataWrangle <- function(taus, thresh){
  taus <- data.frame(taus[, order(thresh)])
  names(taus) <- signif(sort(thresh),3)
  taus <- reshape2::melt(taus, value.name = "tau", variable.name = "thresh")
}

cl <- parallel::makeCluster(ncores) # Start cluster
doParallel::registerDoParallel(cl)

iter_max <- 5e4
R <- 100
rng <- RNGseq(R, seed) # For nested foreach "doRNG" loops
x0 <- y0 <- vector(mode = "list", R)
for(i in 1:R) { x0[[i]] <- mu + as.vector(Sigma_cholLT%*%rnorm(t)) }
for(i in 1:R) { y0[[i]] <- mu + as.vector(Sigma_cholLT%*%rnorm(t)) }

thresh_hh <- c(1e-9,1e-7,1e-5,1e-4,1e-3,5e-3)
taus_hh <- 
  foreach(thresh = thresh_hh, .combine='cbind') %:% 
  foreach(x0_ = x0, y0_ = y0, r = rng,
          .combine='c', .packages = "rwmcouplings") %dopar% {
    rngtools::setRNG(r)
    cplhughop(x0_,y0_,T_hh_final,B_final,lam,kap,thresh,iter_max,logpi,gradlogpi)$tau
  }
taus_hh <- threshDataWrangle(taus_hh, thresh_hh)

parallel::stopCluster(cl) # Stop cluster
save(taus_hh, file = "hhThresh.RData")

###
# Visualize the data
###
load(file = "hhThresh.RData")

thresh_hh_plot <- 
  ggplot(taus_hh, aes(x = thresh, y = tau)) +
  geom_boxplot() +
  xlab(TeX("Threshold $\\delta$ (squared distance)")) +
  ylab("Meeting time") + 
  coord_trans(y = "log10") +
  # #scale_x_discrete(labels = scales::scientific) +
  scale_y_continuous(labels = scales::scientific,
                     breaks = scales::breaks_log(7),#c(1,3,10,30,100,300,1000)*1e2,
                     minor_breaks = NULL) + #breaks = c(5e5, 1e6, 1.5e6, 2e6, 2.5e6, 3e6),
  theme_bw()
thresh_hh_plot

ggsave(filename = "hh_thresh.pdf",
       plot = thresh_hh_plot,
       device = "pdf",  width = 16, height = 8, units = "cm", bg = "transparent")
