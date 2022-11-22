
library(rwmcouplings)
library(Matrix)
library(doParallel)
library(doRNG)
library(reshape2)
library(ggplot2)
library(latex2exp)
library(ggpubr)
########
# Auxiliary functions
########
# Evaluate log-density, and its gradient, for the stochastic volatility model
svmLoglike <- function(x, y_data, beta = 0.65, phi = 0.98, sig = 0.15) {
  -0.5 * (sum(x) + 
            (1/beta^2) * sum(y_data^2 * exp(-x)) + 
            (1/sig^2) * sum((phi * x[-length(x)] - x[-1])^2) + 
            (1-phi^2) * x[1]^2 / sig^2)
}
svmGradLoglike <- function(x, y_data, beta = 0.65, phi = 0.98, sig = 0.15) {
  diff_x <- phi * x[-length(x)] - x[-1]
  temp <- 0.5 * (1/beta^2) * y_data^2 * exp(-x)
  temp[1] <- temp[1] - ((1-phi^2) / sig^2) * x[1]^2
  temp[-1] <- temp[-1] + (1/sig^2) * diff_x
  temp[-length(temp)] <- temp[-length(temp)] - (phi / sig^2) * diff_x
  return(temp-0.5)
}

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

########

########
# Preliminary steps: params, RNG seeding, data generation, Laplace approximation
########

####
# Model parameters
t <- 360L # dimension of the problem

beta <- 0.65
sig  <- 0.15
phi  <- 0.98
####

###
# Set seed for RNGs
seed <- 12345
set.seed(seed)
SetSeed_cpp(seed)
###

###
# Make and save model data
y_data <- beta * rnorm(t) * exp(0.5 * SampleLatentVariables(t, sig, phi))
save(y_data, file = "svmData.RData")
###

###
# Run "optim" to find Laplace approximation. Used to predict where CRN and Reflection couplings will asymptote.
x0_optim <- SampleLatentVariables(t, sig, phi)

optim_out <- optim(par=x0_optim, fn=svmLoglike, gr=svmGradLoglike,
                   y_data = y_data,
                   beta = 0.65, phi = 0.98, sig = 0.15,
                   method = "L-BFGS-B",
                   hessian = TRUE,
                   control = list(fnscale = -1,  maxit = 2000))
optim_out$hessian[abs(optim_out$hessian)<1e-6] <- 0 # The Hessian is sparse, make the approximation sparse too.

Omega_optim <- as(-optim_out$hessian, "dsCMatrix")
Sigma_optim <- solve(Omega_optim)

avg_prec <- mean(diag(Omega_optim)) # Trace of precision matrix at mode
avg_cov  <- mean(diag(Sigma_optim)) # Trace of covariance matrix at mode
########
#Frobenius norm 
# sum(expm::sqrtm(Sigma_optim)^2)

########
# 1. Quantifying the convergence of the RWM: comparison of couplings
########

# Step size
l <- 2.38/sqrt(avg_prec) # Optimal step size for mixing. Is about 0.25, confirming what was found in Papp&Sherlock(2022) by optimizing the acceptance rate.
h <- l / sqrt(t)


# Number of iterations
L <- 2e6 #  Preliminary L iterations for X-chain
iter <- 2e6 # Coupled iterations
thin <- 1e3

# Predict asymptotes
asymptote_out <- squareDistAsymptote(l,t,avg_cov,avg_prec)


######
# Do MCMC in parallel
######
R <- 100 # Number of replicates

# Parallel setup
ncores <- 7
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)


# First L iterations are in common for all couplings
xL <- foreach(i = 1:R, .packages = "rwmcouplings") %dopar% {
  SetSeed_cpp(seed, i)
  x0 <- SampleLatentVariables(t, sig, phi)
  firstL <- svmRWM(x0,y_data,beta,sig,phi,h,L)
  firstL$x
}

###
# Coupled iterations
###

# CRN coupling ###
crn <- foreach(i=1:2,.packages = "rwmcouplings", xL_ = xL) %dopar% {
  SetSeed_cpp(seed, i + R)
  y0 <- SampleLatentVariables(t, sig, phi)
  crn_out <- svmCRNRWM(xL_,y0,y_data,beta,sig,phi,h,iter,thin)
  crn_out$squaredist
}
####


# Reflection coupling ###
refl <- foreach(i=1:2,.packages = "rwmcouplings", xL_ = xL,
                .combine='comb',
                .multicombine=TRUE,
                .init=list(list(), list())) %dopar% {
  SetSeed_cpp(seed, i + R)
  y0 <- SampleLatentVariables(t, sig, phi)
  refl_out <- svmReflMaxRWM(xL_,y0,y_data,beta,sig,phi,h,iter,thin)
  refl_out
}; names(refl) <- c("squaredist", "tau")
####


# GCRN coupling ###

########
# 1.1 Check meeting times vs. the threshold in the two-scale coupling

# Do MCMC
R_thresh_test <- 100

iter_thresh_test <- 1e7

thresh_test <- c(0.001, 0.01, 0.1, 1)
prob_thresh_test <- 2 * pnorm(sqrt(thresh_test), sd = h, lower.tail = FALSE)
prob_thresh_test

taus_thresh_test <- matrix(NA, nrow = R_thresh_test, ncol = length(thresh_test))
for(j in 1:length(thresh_test)) {
  out <- foreach(i=1:R_thresh_test,.packages = "rwmcouplings", xL_ = xL,
                 .combine='comb',
                 .multicombine=TRUE,
                 .init=list(list(), list())) %dopar% {
                   SetSeed_cpp(seed, i + R)
                   y0 <- SampleLatentVariables(t, sig, phi)
                   gcrn_out <- svmTwoScaleGCRNRWM(xL_,y0,y_data,beta,sig,phi,h,iter_thresh_test,1e6,thresh_test[j])
                   gcrn_out
                 }; names(out) <- c("squaredist", "tau")
  taus_thresh_test[,j] <- unlist(out$tau)
}
save(thresh_test, taus_thresh_test, file = "svmThresh.RData")


load(file = "svmThresh.RData")

thresh_plot_df <- data.frame(taus_thresh_test)
names(thresh_plot_df) <- thresh_test
thresh_plot_df <- melt(thresh_plot_df, value.name = "meet", variable.name = "thresh")

# Plot
thresh_plot <- ggplot(thresh_plot_df, aes(x = thresh, y = meet)) +
  geom_boxplot() +
  xlab("Threshold (squared distance)") +
  ylab("Meeting time") + 
  coord_trans(y = "log10") +
  scale_y_continuous(breaks = c(3e5, 5e5, 1e6, 3e6, 5e6), minor_breaks = NULL,
                     labels = scales::scientific) +
  theme_bw()
thresh_plot

ggsave(filename = "svm_thresh.pdf", 
       plot = thresh_plot,
       device = "pdf",  width = 16, height = 8, units = "cm", bg = "transparent")

########
# Actually do GCRN coupling

thresh <- 0.1

gcrn <- foreach(i=1:R,.packages = "rwmcouplings", xL_ = xL,
                               .combine='comb',
                               .multicombine=TRUE,
                               .init=list(list(), list())) %dopar% {
  SetSeed_cpp(seed, i + R)
  y0 <- SampleLatentVariables(t, sig, phi)
  gcrn_out <- svmTwoScaleGCRNRWM(xL_,y0,y_data,beta,sig,phi,h,iter,thin,thresh)
  gcrn_out
}; names(gcrn) <- c("squaredist", "tau")
####

# Stop parallel clusters
parallel::stopCluster(cl)

# Save data
save(xL, crn, refl, gcrn, file = "svmMCMCOutput.RData")

######
# Produce plots
######
load(file = "svmMCMCOutput.RData")
# Data wrangling ###
iters <- seq(0, iter, thin)

w2sq_gcrn <- Reduce("+", gcrn$squaredist) / R
sqdist_crn  <- crn[[1]]
sqdist_refl <- refl$squaredist[[1]]

w2sq_df <- data.frame("GCRN" = w2sq_gcrn, "CRN" = sqdist_crn, "Reflection" = sqdist_refl, "iter" = iters)
w2sq_df <- melt(w2sq_df, id.vars = "iter", variable.name = "Coupling", value.name = "squaredist")
w2sq_df$Linetype <- "MCMC output"
levels(w2sq_df$Coupling) <- c("GCRN (2-scale)", "CRN", "Reflection")

asympt_df <- data.frame("Coupling" = c("GCRN", "CRN", "Reflection"), "squaredist" = c(0, asymptote_out$crn, asymptote_out$refl))
asympt_df$Linetype <- "Predicted asymptote"
asympt_df$Coupling <- c("GCRN (2-scale)", "CRN", "Reflection")

svm_w2sq <-
  ggplot(w2sq_df, aes(x = iter, y = squaredist, color = Coupling, linetype = Linetype)) +
  geom_line() +
  geom_hline(data = asympt_df, aes(yintercept = squaredist, 
                                   color = Coupling, 
                                   linetype = Linetype)) +
  scale_linetype_manual(values = c(1,2)) +
  coord_trans(y = "log1p") + 
  scale_y_continuous(breaks = c(0, 1,10,100, 500), minor_breaks = NULL) +
  ylab("Squared Wasserstein distance") + xlab("Iteration") +
  theme_bw() +
  theme(legend.position = "bottom")
svm_w2sq

taus_gcrn <- unlist(gcrn$tau)
tvd_gcrn  <- getTVDBound(iters, taus_gcrn)
tvd_df <- data.frame("TVD" = tvd_gcrn, "iter" = iters)
svm_tvd <-
  ggplot(tvd_df, aes(x = iter, y = TVD)) +
  geom_line(color = "#F8766D") +
  geom_hline(aes(yintercept = 0), linetype = 2, color = "#F8766D") +
  ylab("Total variation distance") + xlab("Iteration") +
  theme_bw()
svm_tvd

# Join plots into one, with common legend ###
comb_plot <- ggarrange(svm_w2sq, svm_tvd, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
comb_plot

ggsave(filename = "svm_convergence.pdf",
       plot = comb_plot,
       device = "pdf",  width = 24, height = 8, units = "cm", bg = "transparent")


################################################################################
# 2. Quantifying the bias of Laplace approximation to the SVM 
################################################################################

# Laplace approximation parameters:
mu <- optim_out$par # Mean,
Omega <- Omega_optim #precision matrix
Omega_chol <- chol(Omega_optim) # upper triangular Cholesky factor of precision matrix
Sigma_cholLT <- t(chol(Sigma_optim)) # lower triangular Cholesky factor of covariance matrix

# Correct sparse types for RcppEigen
Omega <- as(Omega, "dgCMatrix")
Omega_chol <- as(Omega_chol, "dgCMatrix")

# sum(mu^2); avg_cov # Relevant scales of Laplace approximation

#Step size
l <- 2.38 / sqrt(avg_prec)
h <- l / sqrt(t)
  
# RNG seed
set.seed(seed)
SetSeed_cpp(seed)

# Replicates, iteration count, thinning
R_bias <- 100
iter_bias <- 1e6
thin_bias <- 1e2

# Parallel setup
ncores <- 7

cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

# Starting values for Laplace approximation
y0 <- vector(mode = "list", R_bias)
for(i in 1:R_bias) { y0[[i]] <- mu + as.vector(Sigma_cholLT%*%rnorm(t)) }

gcrn_bias <- foreach(i=1:R_bias,.packages = c("rwmcouplings"), xL_ = xL, y0_ = y0) %dopar% {
  SetSeed_cpp(seed, i+R_bias)
  return(svmLaplaceGCRNRWM(xL_,y0_,y_data,beta,sig,phi,mu,Omega,Omega_chol,h,iter_bias,thin_bias)$squaredist)
}

crn_bias <- foreach(i=1:R_bias,.packages = c("rwmcouplings"), xL_ = xL, y0_ = y0) %dopar% {
  SetSeed_cpp(seed, i+R_bias)
  return(svmLaplaceCRNRWM(xL_,y0_,y_data,beta,sig,phi,mu,Omega,Omega_chol,h,iter_bias,thin_bias)$squaredist)
}

refl_bias <- foreach(i=1:R_bias,.packages = c("rwmcouplings"), xL_ = xL, y0_ = y0) %dopar% {
  SetSeed_cpp(seed, i+R_bias)
  return(svmLaplaceReflRWM(xL_,y0_,y_data,beta,sig,phi,mu,Omega,Omega_chol,h,iter_bias,thin_bias)$squaredist)
}

# Stop parallel clusters
parallel::stopCluster(cl)

save(gcrn_bias, crn_bias, refl_bias, file = "svmMCMCbiasOutput.RData")


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
  mn   <-  Reduce("+", lst) / length(lst)
  mnsq <- Reduce("+", lapply(lst, "^", 2)) / length(lst)
  stddev <- sqrt(mnsq - mn^2) / sqrt(length(lst))
  df_ <- data.frame("mean" = mn, "sd" = stddev, "iter" = iters)
  df_$coupling <- name
  df_
}

gcrn_df <- meansAndSds(gcrn_bias, iters, "GCRN")
crn_df  <- meansAndSds(crn_bias, iters, "CRN")
refl_df <- meansAndSds(refl_bias, iters, "Reflection")
plot_df <- rbind(gcrn_df, crn_df, refl_df)

svm_laplace_avgtrace <- 
  ggplot(plot_df, aes(x = iter, y = mean, color = coupling, fill = coupling)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean - 2*sd, ymax = mean + 2*sd, color = NULL), alpha = 0.3) +
  #coord_trans(y = "log10") +
  xlab("Iteration") + ylab("Squared Wasserstein distance") +
  labs(color = "Coupling", fill = "Coupling") +
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
gcrn_df <- meltDf(gcrn_bias, iters, "GCRN")
crn_df <- meltDf(crn_bias, iters, "CRN")
refl_df <- meltDf(refl_bias, iters, "Reflection")

bias_df <- rbind(gcrn_df, crn_df, refl_df)
bias_df_noburn <- aggregate(squaredist ~  R + coupling, bias_df[bias_df$iters >= 3e5,], mean)

getMeanSd <- function(df_, R) {
  mn_ <- aggregate(squaredist ~ coupling, df_, mean)
  sd_ <- aggregate(squaredist ~ coupling, df_, sd)
  data.frame("coupling" = mn_$coupling, "mean" = mn_$squaredist, "sd" = sd_$squaredist / sqrt(R))
}
plot_df_noburn <- getMeanSd(bias_df_noburn, R_bias)

y_breaks <- c(0, 1.28, 3, 10, 36.5)

svm_laplace_boxerror <- 
  ggplot(plot_df_noburn, aes(x = coupling, y = mean, color = coupling, fill = coupling)) +
  geom_point() +
  geom_errorbar(aes(ymin=mean-2*sd, ymax=mean+2*sd)) +
  coord_trans(y = "log1p", ylim = c(0, 40)) +
  scale_y_continuous(breaks = y_breaks, minor_breaks = NULL) +
  xlab("Coupling") + ylab("Squared Wasserstein distance") +
  theme_bw()
svm_laplace_boxerror

comb_plot_bias <- ggarrange(svm_laplace_avgtrace, 
                       svm_laplace_boxerror, 
                       ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
comb_plot_bias
ggsave(filename = "svm_laplace.pdf",
       plot = comb_plot_bias,
       device = "pdf",  width = 24, height = 8, units = "cm", bg = "transparent")





################################################################################
# 3. Hug and Hop
################################################################################

logpi <-function(x){svmLoglike(x, y_data)}
gradlogpi <-function(x){svmGradLoglike(x, y_data)}

L <- 6e3

Time <- 0.5
Bounces <- 10

lam <- 20
kap <- 1

seed <- 12345
set.seed(seed)
SetSeed_cpp(seed)

# Parallel setup
R <- 100
ncores <- 7
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

x0 <- vector("list", R)
for(i in 1:R) {
  SetSeed_cpp(seed, i)
  x0[[i]] <- SampleLatentVariables(t, sig, phi)
}

# First L iterations
xL_hh <- foreach(i = 1:R, .packages = "rwmcouplings", x0_ = x0) %dorng% {
  firstL <- hughop(x0_,Time,Bounces,lam,kap,L,logpi,gradlogpi)
  firstL$x
}

# save(xL_hh, file = "hhPrelim.RData")
# load(file = "hhPrelim.RData")

###
# Coupled iterations
###
# Check meeting times vs. the threshold in the two-scale coupling ########
R_thresh <- R
iter_thresh <- 1e5
thresh <- c(1e-8,1e-6,1e-2,1e-4)

taus_thresh <- matrix(NA, nrow = R_thresh, ncol = length(thresh))
for(j in 1:length(thresh)) {
  out <- foreach(i=1:R_thresh,.packages = "rwmcouplings", xL_ = xL_hh,
                 .combine='comb',
                 .multicombine=TRUE,
                 .init=list(list(), list())) %dorng% {
                   SetSeed_cpp(seed, i + R)
                   y0 <- SampleLatentVariables(t, sig, phi)
                   gcrn_out <- cplhughop(xL_,y0,Time,Bounces,lam,kap,thresh[j],iter_thresh,logpi,gradlogpi)
                   gcrn_out 
                 }; names(out) <- c("squaredist", "tau")
                 taus_thresh[,j] <- unlist(out$tau)
  if(thresh[j] == 1e-4) {hh <- out}
}
# Stop parallel clusters
parallel::stopCluster(cl)

# Save data
save(xL_hh, hh, file = "hh.RData")
save(thresh, taus_thresh, file = "hhThresh.RData")

thresh_plot_df <- data.frame(taus_thresh)
thresh_plot_df <- thresh_plot_df[, order(thresh)]
names(thresh_plot_df) <- sort(thresh)
thresh_plot_df <- melt(thresh_plot_df, value.name = "meet", variable.name = "thresh")

# Plot
thresh_plot <- 
  ggplot(thresh_plot_df, aes(x = thresh, y = meet)) +
  geom_boxplot() +
  xlab("Threshold (squared distance)") +
  ylab("Meeting time") + 
  coord_trans(y = "log10", ylim = c(1e2,1e5)) +
  #scale_x_discrete(labels = scales::scientific) +
  scale_y_continuous(labels = scales::scientific,
                     breaks = c(1,3,10,30,100,300,1000)*1e2,
                     minor_breaks = NULL) + #breaks = c(5e5, 1e6, 1.5e6, 2e6, 2.5e6, 3e6), 
  theme_bw()
thresh_plot

ggsave(filename = "hh_thresh.pdf", 
       plot = thresh_plot,
       device = "pdf",  width = 16, height = 8, units = "cm", bg = "transparent")

#####
# Plotting ###
iter_plot <- 5e3
iters <- seq(0, iter_plot)
w2sq_df <- data.frame("squaredist" = Reduce("+", hh$squaredist)[0:iter_plot+1] / R_thresh, "iter" = iters)

hh_w2sq <-
  ggplot(w2sq_df, aes(x = iter, y = squaredist)) +
  geom_line() +
  geom_hline(yintercept = 0,linetype = 2) +
  coord_trans(y = "log1p") + 
  scale_y_continuous(breaks = c(0, 1,10,100, 500), minor_breaks = NULL) +
  ylab("Squared Wasserstein distance") + xlab("Iteration") +
  theme_bw() +
  theme(legend.position = "bottom")
hh_w2sq

taus_gcrn <- unlist(hh$tau)
tvd_gcrn  <- getTVDBound(iters, taus_gcrn)
tvd_df <- data.frame("TVD" = tvd_gcrn, "iter" = iters)
hh_tvd <-
  ggplot(tvd_df, aes(x = iter, y = TVD)) +
  geom_line() +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  ylab("Total variation distance") + xlab("Iteration") +
  theme_bw()
hh_tvd

# Join plots into one, with common legend ###
comb_plot <- ggarrange(hh_w2sq, hh_tvd, ncol=2, nrow=1)
comb_plot

ggsave(filename = "hh_convergence.pdf",
       plot = comb_plot,
       device = "pdf",  width = 24, height = 8, units = "cm", bg = "transparent")



###
# EXTRA: check that two-scale CRN coupling performs poorly for H&H
###

load(file = "hh.RData")

set.seed(seed)
ncores <- 7
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

out <- foreach(i=1:10,.packages = "rwmcouplings", xL_ = xL_hh,
               .combine='comb',
               .multicombine=TRUE,
               .init=list(list(), list())) %dorng% {
                 SetSeed_cpp(seed, i + R)
                 y0 <- SampleLatentVariables(t, sig, phi)
                 gcrn_out <- crnhughop(xL_,y0,Time,Bounces,lam,kap,1e-4,1e5,logpi,gradlogpi)
                 gcrn_out 
               }; names(out) <- c("squaredist", "tau")
unlist(out$tau)

parallel::stopCluster(cl)