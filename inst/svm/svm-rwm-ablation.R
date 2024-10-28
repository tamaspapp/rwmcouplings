
library(rwmcouplings)
library(Matrix)
library(doParallel)
library(doRNG)
library(reshape2)
library(ggplot2)
library(latex2exp)
library(ggpubr)

# Parallel setup
ncores <- 6
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
  ggplot(taus_rwm, aes(x = thresh, y = tau, color = Coupling)) + 
  geom_boxplot() +
  labs(x = TeX("Threshold $\\delta$ (squared distance)"), y = "Meeting time") +
  theme_bw() +
  coord_trans(y = "log10")+
  scale_y_continuous(breaks = scales::breaks_log(7)) +
  stat_summary(fun.y = mean, geom = "errorbar", 
               aes(ymax = after_stat(y), ymin = after_stat(y), group = factor(Coupling)),
               width = 0.65, linetype = "dashed", position = position_dodge(0.75), linewidth = 0.75)
thresh_rwm_plot

ggsave(filename = "rwm_thresh.pdf",
       plot = thresh_rwm_plot,
       device = "pdf",  width = 24, height = 12, units = "cm", bg = "transparent")
