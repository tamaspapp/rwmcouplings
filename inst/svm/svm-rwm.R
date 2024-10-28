
library(rwmcouplings)
library(Matrix)
library(doParallel)
library(doRNG)
library(reshape2)
library(ggplot2)
library(latex2exp)
library(ggpubr)

ncores <- 6
cl <- parallel::makeCluster(ncores) # Parallel computing setup
doParallel::registerDoParallel(cl)

########
# 1. Convergence of RWM
########
# Load in model data 
if(!file.exists("svmData.RData")) {source("svm-prelim.R")} # Ensure we run this first
load(file = "svmData.RData")

# Set seed for RNGs
seed <- 12345; set.seed(seed); SetSeed_cpp(seed)

# Step size
l <- 2.38/sqrt(avg_prec)
h <- l / sqrt(t)
# Acceptance rate
print(paste("RWM acceptance rate is:", signif(svmRWM(SampleLatentVariables(t, sig, phi),y_data,beta,sig,phi,h,1e6)$acc_rate, 3)))

# Predict asymptotes
asymptote_out <- squareDistAsymptote(l,t,avg_cov,avg_prec)

# Parameters for parallel MCMC
R <- 100    # Number of replicates
L <- 1.5e6    # Lag
iter <- 2e6 # Coupled iterations
thin <- 1e3

# First L iterations are in common for all couplings
xL <- foreach(i = 1:R,.packages = "rwmcouplings") %dopar% {
  SetSeed_cpp(seed, i)
  x0 <- SampleLatentVariables(t, sig, phi)
  firstL <- svmRWM(x0,y_data,beta,sig,phi,h,L)
  firstL$x
}

####
# Coupled iterations
# CRN coupling
crn <- foreach(i=1:2,xL_ = xL,.packages = "rwmcouplings") %dopar% {
  SetSeed_cpp(seed, i + R)
  y0 <- SampleLatentVariables(t, sig, phi)
  crn_out <- svmCRNRWM(xL_,y0,y_data,beta,sig,phi,h,iter,thin)
  crn_out$squaredist
}

# Reflection coupling
refl <- foreach(i=1:2,xL_ = xL,.packages = "rwmcouplings", 
                .combine='comb',.multicombine=TRUE,.init=list(list(), list())) %dopar% {
                  SetSeed_cpp(seed, i + R)
                  y0 <- SampleLatentVariables(t, sig, phi)
                  refl_out <- svmReflMaxRWM(xL_,y0,y_data,beta,sig,phi,h,iter,thin)
                  refl_out
                }; names(refl) <- c("squaredist", "tau")

# GCRN coupling
thresh_gcrn <- 0.1
gcrn <- foreach(i=1:R, xL_ = xL,.packages = "rwmcouplings",
                .combine='comb',.multicombine=TRUE,.init=list(list(), list())) %dopar% {
                  SetSeed_cpp(seed, i + R)
                  y0 <- SampleLatentVariables(t, sig, phi)
                  gcrn_out <- svmTwoScaleGCRNRWM(xL_,y0,y_data,beta,sig,phi,h,iter,thin,thresh_gcrn)
                  gcrn_out
                }; names(gcrn) <- c("squaredist", "tau")

# GCRefl coupling
thresh_gcrefl <- 0.001
gcrefl <- foreach(i=1:R,xL_ = xL,.packages = "rwmcouplings",
                  .combine='comb',.multicombine=TRUE,.init=list(list(), list())) %dopar% {
                    SetSeed_cpp(seed, i + R)
                    y0 <- SampleLatentVariables(t, sig, phi)
                    gcrefl_out <- svmTwoScaleGCReflRWM(xL_,y0,y_data,beta,sig,phi,h,iter,thin,thresh_gcrefl)
                    list(gcrefl_out$squaredist, gcrefl_out$tau)
                  }; names(gcrefl) <- c("squaredist", "tau")

####
# Save data
save(xL, crn, refl, gcrn, gcrefl, file = "svmRWMconvergence.RData")
####

parallel::stopCluster(cl)



####
# Plot bounds
load(file = "svmRWMconvergence.RData")
####
# W2^2 bound; assumes all chains have met before time L, which is the case with this setup

# Some data wrangling to do
iters <- seq(0, iter, thin)
w2sq_df <- data.frame("GCRN"   = Reduce("+", gcrn$squaredist) / R, 
                      "GCRefl" = Reduce("+", gcrefl$squaredist) / R, 
                      "CRN"    = crn[[1]], 
                      "Reflection" = refl$squaredist[[1]], 
                      "iter" = iters)
w2sq_df <- melt(w2sq_df, id.vars = "iter", variable.name = "Coupling", value.name = "squaredist")
w2sq_df$Linetype <- "MCMC output"
levels(w2sq_df$Coupling) <- c("GCRN", "GCRefl", "CRN", "Reflection")

asympt_df <- data.frame("Coupling" = c("GCRN", "GCRefl", "CRN", "Reflection"), 
                        "squaredist" = c(NA,NA, asymptote_out$crn, asymptote_out$refl))
asympt_df$Linetype <- "Predicted asymptote"
asympt_df$Coupling <- c("GCRN", "GCRefl", "CRN", "Reflection")
# End data wrangling

svm_w2sq <-
  ggplot(w2sq_df, aes(x = iter, y = squaredist, color = Coupling, linetype = Linetype)) +
  geom_line() +
  geom_hline(data = asympt_df, aes(yintercept = squaredist, 
                                   color = Coupling, 
                                   linetype = Linetype)) +
  scale_linetype_manual(values = c(1,2)) +
  coord_trans(y = "log1p", xlim = c(0, 1.5e6)) + 
  scale_y_continuous(breaks = c(0,1,10,100,500), minor_breaks = NULL) +
  ylab("Squared Wasserstein distance") + xlab("Iteration") +
  theme_bw() +
  theme(legend.position = "bottom")
svm_w2sq
####

# TVD bound
tvd_df <- data.frame("GCRN"       = getTVDBound(iters, unlist(gcrn$tau)),
                     "GCRefl"     = getTVDBound(iters, unlist(gcrefl$tau)),
                     "CRN"        = rep(NA, length(iters)),
                     "Reflection" = rep(NA, length(iters)),
                     "iter" = iters)
tvd_df <- melt(tvd_df, id.vars = "iter", variable.name = "Coupling", value.name = "TVD")

svm_tvd <-
  ggplot(tvd_df, aes(x = iter, y = TVD, color = Coupling)) +
  geom_line() +
  # geom_hline(aes(yintercept = 0), linetype = 2, color = "black") +
  ylab("Total variation distance") + xlab("Iteration") +
  coord_cartesian(xlim = c(0, 1.5e6)) +
  theme_bw()
svm_tvd


# Join plots into one, with common legend ###
comb_plot <- ggarrange(svm_w2sq, svm_tvd, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
comb_plot

ggsave(filename = "svm_convergence.pdf",
       plot = comb_plot,
       device = "pdf",  width = 24, height = 8, units = "cm", bg = "transparent")






################################################################################
# 2. Bias of Laplace approximation
################################################################################

# Load in model data 
if(!file.exists("laplaceBiasGelbrich.RData")) {source("svm-gelbrich.R")} # Ensure we run this first
load(file = "laplaceBiasGelbrich.RData")

load(file = "svmData.RData")
load(file = "svmRWMconvergence.RData")
seed <- 12345; set.seed(seed); SetSeed_cpp(seed)

# Replicates, iteration count, thinning
R_bias <- R # 100
iter_bias <- 1e6
thin_bias <- 1e2

# Parallel setup
ncores <- 6
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

y0 <- vector(mode = "list", R_bias) # Starting values, chain targeting Laplace approximation
for(i in 1:R_bias) { y0[[i]] <- mu + as.vector(Sigma_cholLT%*%rnorm(t)) }

gcrn_bias <- foreach(i=1:R_bias, xL_ = xL, y0_ = y0, .packages = c("rwmcouplings")) %dopar% {
  SetSeed_cpp(seed, i+R_bias)
  svmLaplaceGCRNRWM(xL_,y0_,y_data,beta,sig,phi,mu,Omega,Omega_chol,h,iter_bias,thin_bias)$squaredist
}

crn_bias <- foreach(i=1:R_bias, xL_ = xL, y0_ = y0, .packages = c("rwmcouplings")) %dopar% {
  SetSeed_cpp(seed, i+R_bias)
  svmLaplaceCRNRWM(xL_,y0_,y_data,beta,sig,phi,mu,Omega,Omega_chol,h,iter_bias,thin_bias)$squaredist
}

refl_bias <- foreach(i=1:R_bias, xL_ = xL, y0_ = y0, .packages = c("rwmcouplings")) %dopar% {
  SetSeed_cpp(seed, i+R_bias)
  svmLaplaceReflRWM(xL_,y0_,y_data,beta,sig,phi,mu,Omega,Omega_chol,h,iter_bias,thin_bias)$squaredist
}

gcrefl_bias <- foreach(i=1:R_bias, xL_ = xL, y0_ = y0, .packages = c("rwmcouplings")) %dopar% {
  SetSeed_cpp(seed, i+R_bias)
  svmLaplaceGCReflRWM(xL_,y0_,y_data,beta,sig,phi,mu,Omega,Omega_chol,h,iter_bias,thin_bias)$squaredist
}

save(gcrn_bias, crn_bias, refl_bias, gcrefl_bias, file = "svmMCMCbiasOutput.RData")
parallel::stopCluster(cl)

#################
# Produce plots #
#################  
load(file = "svmMCMCbiasOutput.RData")
iters <- seq(0,iter_bias,thin_bias)

##############
# Pointwise averages w/ error bands
##############

# Iteration-wise means and standard deviations
meansAndSds <- function(lst, iters,  name) {
  mn   <- Reduce("+", lst) / length(lst)
  mnsq <- Reduce("+", lapply(lst, "^", 2)) / length(lst)
  stderr <- sqrt(mnsq - mn^2) / sqrt(length(lst))
  df_ <- data.frame("mean" = mn, "se" = stderr, "iter" = iters)
  df_$coupling <- name
  df_
}

gcrn_df   <- meansAndSds(gcrn_bias, iters, "GCRN")
crn_df    <- meansAndSds(crn_bias, iters, "CRN")
refl_df   <- meansAndSds(refl_bias, iters, "Reflection")
gcrefl_df <- meansAndSds(gcrefl_bias, iters, "GCRefl")

plot_df <- rbind(gcrn_df, crn_df, refl_df, gcrefl_df)

svm_laplace_avgtrace <- 
  ggplot(plot_df, aes(x = iter, y = mean, color = coupling, fill = coupling)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean - 2*se, ymax = mean + 2*se, color = NULL), alpha = 0.3) +
  #coord_trans(y = "log1p") +
  xlab("Iteration") + ylab("Squared Wasserstein distance") +
  labs(color = "Coupling", fill = "Coupling") +
  geom_hline(yintercept = w2sq_gelbrich_lb, linetype = "dashed")+
  annotate('ribbon', x = c(-Inf, Inf), # Add Gelbrich lower bound
           ymin = w2sq_gelbrich_lb - 2 * w2sq_gelbrich_jackse, 
           ymax = w2sq_gelbrich_lb + 2 * w2sq_gelbrich_jackse, 
           alpha = 0.3, fill = "black") +
  theme_bw()
svm_laplace_avgtrace

##############
# Estimates and error bars
##############
meltDf <- function(list_, iters, name) {
  df_ <- data.frame(list_)
  names(df_) <- seq(1,ncol(df_))
  df_ <- cbind(df_, iters)
  df_ <- melt(df_, id.vars = c("iters"), variable = c("R"), value.name = c("squaredist"))
  df_$coupling <- name
  return(df_)
}
gcrn_df   <- meltDf(gcrn_bias, iters, "GCRN")
crn_df    <- meltDf(crn_bias, iters, "CRN")
refl_df   <- meltDf(refl_bias, iters, "Reflection")
gcrefl_df <- meltDf(gcrefl_bias, iters, "GCRefl")

burnin <- 3e5

bias_df <- rbind(gcrn_df, crn_df, refl_df, gcrefl_df)
bias_df_noburn <- aggregate(squaredist ~  R + coupling, bias_df[bias_df$iters >= burnin,], mean)

getMeanSd <- function(df_, R) {
  mn_ <- aggregate(squaredist ~ coupling, df_, mean)
  sd_ <- aggregate(squaredist ~ coupling, df_, sd)
  data.frame("coupling" = mn_$coupling, "mean" = mn_$squaredist, "se" = sd_$squaredist / sqrt(R))
}
plot_df_noburn <- getMeanSd(bias_df_noburn, R_bias)

y_breaks <- c(0, 1, 6.8, 36.5) #signif(plot_df_noburn$mean, 2)

svm_laplace_boxerror <- 
  ggplot(plot_df_noburn, aes(x = coupling, y = mean, color = coupling, fill = coupling)) +
  geom_point() +
  geom_errorbar(aes(ymin=mean-2*se, ymax=mean+2*se)) +
  coord_trans(y = "log1p", ylim = c(0, 40)) +
  scale_y_continuous(breaks = y_breaks, minor_breaks = NULL) +
  xlab("Coupling") + ylab("Squared Wasserstein distance") +
  # Add Gelbrich lower bound
  geom_hline(yintercept = w2sq_gelbrich_lb, linetype = "dashed")+
  annotate('ribbon', x = c(-Inf, Inf), 
           ymin = w2sq_gelbrich_lb - 2 * w2sq_gelbrich_jackse, 
           ymax = w2sq_gelbrich_lb + 2 * w2sq_gelbrich_jackse, 
           alpha = 0.3, fill = 'black') +
  theme_bw()
svm_laplace_boxerror

comb_plot_bias <- ggarrange(svm_laplace_avgtrace, 
                            svm_laplace_boxerror, 
                            ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
comb_plot_bias
ggsave(filename = "svm_laplace.pdf",
       plot = comb_plot_bias,
       device = "pdf",  width = 24, height = 8, units = "cm", bg = "transparent")



####
# Bonus: can we get a non-trivial upper bound on the total variation distance?
####

cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

gcrefl_bias_tvd <- foreach(i=1:R_bias, xL_ = xL, y0_ = y0, .packages = c("rwmcouplings"), .combine = "c") %dopar% {
  SetSeed_cpp(seed, i+R_bias)
  out <- svmLaplacetwoscaleGCReflRWM(xL_,y0_,y_data,beta,sig,phi,mu,Omega,Omega_chol,h,iter_bias,thin_bias,thresh_gcrefl)
  mean(!out$coalesced[burnin:(iter_bias + 1)])
}

parallel::stopCluster(cl)

mean(gcrefl_bias_tvd)
2 * sd(gcrefl_bias_tvd) / sqrt(R_bias)
