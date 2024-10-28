
library(rwmcouplings)
library(Matrix)
library(doParallel)
library(doRNG)
library(reshape2)
library(ggplot2)
library(latex2exp)
library(ggpubr)

# Parallel setup
ncores <- 20
seed <- 12345; set.seed(seed); SetSeed_cpp(seed)

########
# 0. Setup
########

# Probability of RWM proposals coalescing in 1 iteration, under maximal coupling.
probProposalsCoalesced <- function(delta, h) {
  2 * pnorm(sqrt(delta), sd = h, lower.tail = FALSE)
}

# Data wrangling for threshold plots
threshDataWrangle <- function(taus, thresh){
  taus <- data.frame(taus[, order(thresh)])
  names(taus) <- signif(sort(thresh),3)
  taus <- reshape2::melt(taus, value.name = "tau", variable.name = "thresh")
}

# Load in model data 
if(!file.exists("svmData.RData")) {source("svm-prelim.R") } # Ensure we run this first
load(file = "svmData.RData")
########


########
# 1. RWM two-scale couplings: find optimal threshold by grid search
########
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

# Step size
l <- 2.38/sqrt(avg_prec)
h <- l / sqrt(t)


# Parallel MCMC setup
R_thresh <- 100
iter_thresh <- 2e7
thin_thresh <- 1e7 # Pick large; we don't use MCMC output aside from the meeting time

# Start the coupled chains from independent draws from the Laplace approximation
y0_thresh <- x0_thresh <- vector(mode = "list", R_thresh)
for(i in 1:R_thresh) { y0_thresh[[i]] <- mu + as.vector(Sigma_cholLT%*%rnorm(t)) }
for(i in 1:R_thresh) { x0_thresh[[i]] <- mu + as.vector(Sigma_cholLT%*%rnorm(t)) }


# GCRN
thresh_gcrn <- c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 5e-1, 1)
taus_gcrn <- 
  foreach(thresh=thresh_gcrn, .combine='cbind') %:% 
  foreach(i=1:R_thresh, x0_ = x0_thresh, y0_ = y0_thresh, .combine='c', .packages = "rwmcouplings") %dopar% {
    SetSeed_cpp(seed, i + R_thresh)
    svmTwoScaleGCRNRWM(x0_,y0_,y_data,beta,sig,phi,h,iter_thresh,thin_thresh,thresh)$tau
  }
taus_gcrn <- threshDataWrangle(taus_gcrn, thresh_gcrn)
taus_gcrn$Coupling <- "GCRN"

# GCRefl
thresh_gcrefl <- c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 5e-1) # Don't use the last threshold because the meeting times are an order of magnitude larger
taus_gcrefl <- 
  foreach(thresh=thresh_gcrefl, .combine='cbind') %:% 
  foreach(i=1:R_thresh, x0_ = x0_thresh, y0_ = y0_thresh, .combine='c', .packages = "rwmcouplings") %dopar% {
    SetSeed_cpp(seed, i + R_thresh)
    svmTwoScaleGCReflRWM(x0_,y0_,y_data,beta,sig,phi,h,iter_thresh,thin_thresh,thresh)$tau
  }
taus_gcrefl <- threshDataWrangle(taus_gcrefl, thresh_gcrefl)
taus_gcrefl$Coupling <- "GCRefl"

parallel::stopCluster(cl)


taus_rwm <- rbind(taus_gcrn, taus_gcrefl)
save(taus_rwm, file = "rwmThresh.RData")


###
# Visualize the data
###
load(file = "rwmThresh.RData")

probs <- probProposalsCoalesced(thresh_gcrn, h)
print("Probability of coalcescing proposals, at threshold delta:")
print(paste0("delta=", signif(thresh_gcrn,3), ": prob=", signif(probs,3)))

thresh_rwm_plot <- 
  ggplot(taus_rwm, aes(x = thresh, y = tau, fill = Coupling)) + 
  geom_boxplot() +
  labs(x = TeX("Threshold $\\delta$ (squared distance)"), y = "Meeting time") +
  theme_bw() +
  stat_summary(fun.y = mean, geom = "errorbar", 
               aes(ymax = after_stat(y), ymin = after_stat(y), group = factor(Coupling)),
               width = 0.65, linetype = "dashed", position = position_dodge(0.75), linewidth = 0.75)
thresh_rwm_plot

ggsave(filename = "rwm_thresh.pdf",
       plot = thresh_rwm_plot,
       device = "pdf",  width = 24, height = 12, units = "cm", bg = "transparent")

# labs(x = "Threshold (squared distance)", y = "Meeting time") +
#coord_trans(y = "log10") +
#scale_y_continuous(breaks = c(5e4, 1e5, 3e5, 5e5, 1e6, 3e6, 5e6), minor_breaks = NULL,
#                 labels = scales::scientific) +










########
# 2. HMC: contractivity with long and short integration times
########
# There is a sharp phase transition for the contractivity of HMC:
# we need a short enough integration time.
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

R_hmc <- 2

T_hmc    <- c(0.9,0.7,0.5,0.3,0.25,0.2,0.1,0.05)
L_hmc    <- 100
iter_hmc <- 1e3

rng_hmc <- RNGseq(R_hmc, seed + 1) # Set RNG manually to use nested foreach loops

hmc_intTime <- 
  foreach(j=1:length(T_hmc), .combine = "rbind") %:%
  foreach(i=1:R_hmc, x0_ = x0_thresh, y0_ = y0_thresh, r=rng_hmc,.combine='rbind',.packages = "rwmcouplings") %dopar% {
    rngtools::setRNG(r)
    # The zeros below are so that we don't mix in the RWM at all
    cbind(twoscalehmc(x0_,y0_,T_hmc[j],L_hmc,0,0,0,iter_hmc,logpi,gradlogpi)$squaredist,T_hmc[j],i,0:iter_hmc) 
  }

parallel::stopCluster(cl) # Stop parallel clusters

hmc_intTime <- as.data.frame(hmc_intTime)
names(hmc_intTime) <- c("squaredist", "time", "replicate", "iteration")

ggplot(hmc_intTime, aes(y = squaredist, x = iteration, color = factor(time), linetype = factor(replicate))) + 
  geom_line(alpha = 0.5) +
  labs(x = "Iteration", y = "Squared distance", color = "Integration time", title = "HMC", subtitle = "Contraction with CRN coupling") +
  scale_y_log10() +
  theme_bw()



########
# 3. Hug + Hop: select Hug parameters for optimal contractivity.
########
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)


# Fix Hop parameters a priori.
lam <- 20 # Hop params
kap <- 1


R_hh <- 4




T_hh    <- c(0.6,0.5,0.4,0.3,0.2)
B_hh    <- c(10,20,30,40)
iter_hh <- 1e3



rng_hug <- RNGseq(R_hh, seed + 1) # Set RNG manually to use nested foreach loops

hh <- 
  foreach(B_ = B_hh, .combine = "rbind") %:%
  foreach(T_ = T_hh, .combine = "rbind") %:%
  foreach(i=1:R_hh, x0_ = x0_thresh, y0_ = y0_thresh, r=rng_hug,.combine='rbind',.packages = "rwmcouplings") %dopar% {
    rngtools::setRNG(r)
    # Zero below means we do not use maximal coupling
    cbind(cplhughop(x0_,y0_,T_,B_,lam,kap,0,iter_hh,logpi,gradlogpi)$squaredist,
          T_,
          B_,
          i,
          0:iter_hh) 
  }
parallel::stopCluster(cl) # Stop parallel clusters

hh <- as.data.frame(hh)
names(hh) <- c("squaredist", "time", "bounces", "replicate", "iteration")


hh_ <- hh[hh$bounces == 40,]

ggplot(hh_, aes(y = squaredist, x = iteration, 
                color = factor(time), linetype = factor(replicate))) + 
  geom_line(alpha = 0.5) +
  labs(x = "Iteration", y = "Squared distance", color = "Integration time",
       title = "Hug+Hop", subtitle = "Contraction with CRN+GCRN coupling") +
  scale_y_log10() +
  theme_bw()











########
# CHANGE #. Hug and Hop: select threshold for two-scale Hop coupling
########
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

R_hh <- 100
iter_hh <- 1e5

Time <- 0.5; Bounces <- 10 # Hug params
lam <- 20; kap <- 1 # Hop params

thresh_hh <- c(1e-8,1e-6,1e-4,1e-2)
taus_hh <- 
  foreach(thresh=thresh_hh, .combine='cbind') %:% 
  foreach(i=1:R_hh, x0_ = x0_thresh, y0_ = y0_thresh, .combine='c', .packages = "rwmcouplings") %dopar% {
    cplhughop(x0_,y0_,Time,Bounces,lam,kap,thresh,iter_hh,logpi,gradlogpi)$tau
  }
taus_hh <- threshDataWrangle(taus_hh, thresh_hh)
save(thresh_hh, taus_hh, file = "hhThresh.RData")


parallel::stopCluster(cl) # Stop parallel clusters


###
# Visualize the data
###
load(file = "hhThresh.RData")

thresh_hh_plot <- 
  ggplot(taus_hh, aes(x = thresh, y = tau)) +
  geom_boxplot() +
  xlab("Threshold (squared distance)") +
  ylab("Meeting time") + 
  coord_trans(y = "log10", ylim = c(1e2,1e5)) +
  #scale_x_discrete(labels = scales::scientific) +
  scale_y_continuous(labels = scales::scientific,
                     breaks = c(1,3,10,30,100,300,1000)*1e2,
                     minor_breaks = NULL) + #breaks = c(5e5, 1e6, 1.5e6, 2e6, 2.5e6, 3e6), 
  theme_bw()
thresh_hh_plot

ggsave(filename = "hh_thresh.pdf",
       plot = thresh_hh_plot,
       device = "pdf",  width = 16, height = 8, units = "cm", bg = "transparent")




########
# 3. Hug and Hop: two-scale CRN coupling performs poorly
########
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

R_crnhh <- 10
iter_crnhh   <- 1e5 
thresh_crnhh <- 1e-4

crnhh <- foreach(i=1:R_crnhh, x0_ = x0_thresh, y0_ = y0_thresh, .combine='c',.packages = "rwmcouplings") %dorng% {
                 crnhughop(x0_,y0_,Time,Bounces,lam,kap,thresh_crnhh,iter_crnhh,logpi,gradlogpi)$tau
               }
print("Hug and Hop -- meeting times with two-scale CRN coupling for Hop:")
print(as.vector(crnhh))

parallel::stopCluster(cl) # Stop parallel clusters






########
# 5. Hug: vary the integration time
########
crnhug <- function(x0,y0,Time,Bounces,iter,logpi,gradlogpi) {
  
  normalize <- function(x) {return((1/sqrt(sum(x^2)))*x)}
  d = length(x0);
  
  delta = Time / Bounces;
  
  x = x0;
  y = y0;
  
  logpi_x = logpi(x);
  logpi_y = logpi(y);

  squaredist = rep(0, iter+1);
  squaredist[1] = sum((x-y)^2);
  
  acc_x_hug = 0;
  tau = -1;
  for(i in 1:iter) {
    # Hug
    # Propose: bounce B times
    z = rnorm(d); # X-velocity
    z_y = z;  # Y-velocity
    
    # Initial position half-step
    xp = x + 0.5 * delta * z;
    yp = y + 0.5 * delta * z_y;
    g_xp = normalize(gradlogpi(xp));
    g_yp = normalize(gradlogpi(yp));
    
    for(b in 1:Bounces)
    {
      # Reflect velocity in gradient
      z = z -  2 * sum(z * g_xp) * g_xp;
      z_y = z_y -  2 * sum(z_y * g_yp) * g_yp;
      
      # Full step
      xp = xp + delta * z;
      yp = yp + delta * z_y;
      g_xp = normalize(gradlogpi(xp));
      g_yp = normalize(gradlogpi(yp));
    }
    # Went too far
    xp = xp - 0.5 * delta * z;
    yp = yp - 0.5 * delta * z_y;
    
    # Accept-reject
    logpi_xp = logpi(xp);
    logpi_yp = logpi(yp);
    
    logHR_x = logpi_xp - logpi_x;
    logHR_y = logpi_yp - logpi_y;
    
    log_u = log(runif(1));
    if (logHR_x > 0 || log_u < logHR_x)
    {
      x = xp;
      logpi_x = logpi_xp;
      acc_x_hug = acc_x_hug + 1;
    }
    if (logHR_y > 0 || log_u < logHR_y)
    {
      y = yp;
      logpi_y = logpi_yp;
    }
    squaredist[i+1] = sum((x-y)^2);
  }
  
  return(list("squaredist" = squaredist,
              "acc_x_hug" = acc_x_hug / i,
              "tau" = tau))
}

cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

R_hug <- 3

T_hug    <- c(0.5,0.25,0.2)
B_hug    <- 100
iter_hug <- 1e3

rng_hug <- RNGseq(R_hug, seed + 1) # Set RNG manually to use nested foreach loops

hug_intTime <- 
  foreach(T_=T_hug, .combine = "rbind") %:%
  foreach(i=1:R_hug, x0_ = x0_thresh, y0_ = y0_thresh, r=rng_hug,.combine='rbind') %dopar% {
    rngtools::setRNG(r)
    cbind(crnhug(x0_,y0_,T_,B_hug,iter_hug,logpi,gradlogpi)$squaredist,
          T_,
          i,
          0:iter_hug) 
  }

parallel::stopCluster(cl) # Stop parallel clusters

hug_intTime <- as.data.frame(hug_intTime)
names(hug_intTime) <- c("squaredist", "time", "replicate", "iteration")

ggplot(hug_intTime, aes(y = squaredist, x = iteration, color = factor(time), linetype = factor(replicate))) + 
  geom_line(alpha = 0.5) +
  labs(x = "Iteration", y = "Squared distance", color = "Integration time", title = "Hug", subtitle = "Contraction with CRN coupling") +
  scale_y_log10() +
  theme_bw()








