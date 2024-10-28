
#####
# Plots for squared distance between chains over the average covariance of the target, elliptical case
#####

library(Matrix)
library(ggplot2)
library(ggh4x)
library(rwmcouplings)
library(reshape2)
library(latex2exp)

# For optimal Markovian coupling
library(emg)
library(TruncatedNormal)

seed <- 12345

####
# Preliminary functions
####

# Kac-Murdock-Szego matrix and inverse constructors.
KMSMatrix <- function(rho, d) {
  toeplitz(rho^(0:(d-1)))
}
KMSMatrixInverse <- function(rho, d) {
  if(d == 1) {
    M <- 1
  } else {
    M <- Matrix::bandSparse(d, 
                            k = c(0,1), 
                            diag = list( rep((1 + rho^2)/ (1-rho^2),d), rep(-rho/ (1-rho^2), d-1)), 
                            symm=TRUE)
    M[1,1] <- M[d,d] <- 1 / (1-rho^2)
  }
  return(M)
}


# Do MCMC, compute ODE limits, save as long format dataframe for plotting
getPlotDf <- function(x,y,x0,y0,v0,l,avg_cov,d,tmax,length_plot,logpi,gradlogpi) {
  
  h <- l / sqrt(d)
  
  ode_times <- seq(0, tmax, length.out = length_plot)
  
  iter <- tmax * d
  thinned_mcmc <- 1 + floor(seq(0, iter, length.out = length_plot))
  
  # GCRN ####
  set.seed(seed+1) # Same seed across couplings
  out_mcmc <- RWM_gcrn(x,y,h,iter,logpi,gradlogpi)
  squaredist_mcmc <- out_mcmc$squaredist[thinned_mcmc]
  ode_asympt<- rep(NA, length_plot)
  df_gcrn <- data.frame("scaled.squaredist" = squaredist_mcmc/d/avg_cov, "ode.asymptote" = ode_asympt, "iter/d" = ode_times, "Coupling" = "GCRN")
  
  # CRN ####
  set.seed(seed+1)
  out_mcmc <- RWM_crn(x,y,h,iter,logpi,gradlogpi)
  squaredist_mcmc <- out_mcmc$squaredist[thinned_mcmc]
  ode_asympt <- rep(NA, length_plot)
  df_crn <- data.frame("scaled.squaredist" = squaredist_mcmc/d/avg_cov,"ode.asymptote" = ode_asympt, "iter/d" = ode_times, "Coupling" = "CRN")
  
  # Reflection ####
  set.seed(seed+1)
  out_mcmc <- RWM_refl(x,y,h,iter,logpi,gradlogpi)
  squaredist_mcmc <- out_mcmc$squaredist[thinned_mcmc]
  ode_asympt<- rep(NA, length_plot)
  df_refl <- data.frame("scaled.squaredist" = squaredist_mcmc/d/avg_cov, "ode.asymptote" = ode_asympt, "iter/d" = ode_times, "Coupling" = "Reflection")
  
  # Optimal Markovian ####
  set.seed(seed+1)
  out_mcmc <- RWM2targets_optimal(x,y,h,iter,logpi,gradlogpi,logpi,gradlogpi)
  squaredist_mcmc <- out_mcmc$squaredist[thinned_mcmc]
  ode_asympt<- rep(NA, length_plot)
  df_opt <- data.frame("scaled.squaredist" = squaredist_mcmc/d/avg_cov, "ode.asymptote" = ode_asympt, "iter/d" = ode_times, "Coupling" = "Optimal Markovian")
  
  # Join dataframes into a single one, for plotting ####
  df_ <- rbind(df_gcrn,df_crn, df_refl, df_opt)
  df_ <- melt(df_, id.var = c("iter.d", "Coupling"), variable.name = "Linetype", value.name = "scaled.squaredist")
  levels(df_$Linetype) <- c("MCMC", "Predicted asymptote")
  df_ <- na.omit(df_)
  return(df_)
}

d <- 2000
length_plot <- 2001

####
# Target setups
####
ar1_setup <- function(l1, seed, rho = 0.5) {
  set.seed(seed)
  
  # rho <<- 0.5 # AR(1) correlation
  tmax <<- 20
  
  Omega <<- KMSMatrixInverse(rho, d)
  Omega_chol <<- chol(Omega)
  Sigma <<- KMSMatrix(rho, d)
  Sigma_cholLT <<- t(chol(Sigma))
  
  logpi <<- function(x) { - 0.5 * sum(as.vector((Omega_chol %*% x))^2)}
  gradlogpi <<- function(x) { - as.vector(Omega %*% x)}
  
  avg_cov <<- mean(diag(Sigma))  # 1
  avg_prec <<- mean(diag(Omega)) # (1 + rho^2)/ (1-rho^2)
  
  l <<- l1/sqrt(avg_prec)
  
  # MCMC starting values
  x <<- as.vector(Sigma_cholLT%*%rnorm(d))
  y <<- as.vector(Sigma_cholLT%*%rnorm(d))
  
  asymptote_out <<- squareDistAsymptote(l, d, avg_cov, avg_prec) # Asymptote for CRN and reflection
  asymptote_crn <<- asymptote_out$crn/avg_cov/d
  asymptote_refl <<- asymptote_out$refl/avg_cov/d
}

chisq_setup <- function(l1, seed, nu = 3) {
  set.seed(seed)
  
  tmax <<- 40
  # nu <- 3 # Chi-square degrees of freedom, assumed nu>2
  Sigma_diagonal_chisq <- rchisq(d, df = nu)
  
  Sigma <<- Matrix::Diagonal(x = Sigma_diagonal_chisq)
  Sigma_chol <<- Sigma^0.5
  Omega <<- Matrix::Diagonal(x = 1/Sigma_diagonal_chisq)
  Omega_chol <<- Omega^0.5
  
  avg_cov <<- mean(diag(Sigma)) #nu
  avg_prec <<- mean(diag(Omega)) #1/(nu - 2)
  
  l <<- l1/sqrt(avg_prec)
  
  asymptote_out <<- squareDistAsymptote(l, d, avg_cov, avg_prec) # Asymptote for CRN and reflection
  asymptote_crn <<- asymptote_out$crn/avg_cov/d
  asymptote_refl <<- asymptote_out$refl/avg_cov/d
  
  logpi <<- function(x) { - 0.5 * sum((diag(Omega_chol) * x)^2)}
  gradlogpi <<- function(x) { - diag(Omega) * x}
  
  # MCMC starting values
  x <<- as.vector(Sigma_chol%*%rnorm(d))
  y <<- as.vector(Sigma_chol%*%rnorm(d))
}

# Diagonal with alternating eigenvalues Diag(1, sig_sq, 1, sig_sq, ...)
alternating_setup <- function(l1, seed, sig_sq = 24) {
  set.seed(seed)
  
  tmax <<- 80
  # sig_sq <- 24 # Eigenvalue that's not 1
  Sigma_diagonal <- rep(c(1,sig_sq), ceiling(d/2))[1:d]
  
  Sigma <<- Matrix::Diagonal(x = Sigma_diagonal)
  Sigma_chol <<- Sigma^0.5
  Omega <<- Matrix::Diagonal(x = 1/Sigma_diagonal)
  Omega_chol <<- Omega^0.5
  
  avg_cov <<- mean(diag(Sigma)) #(1 + sig_sq)/2
  avg_prec <<- mean(diag(Omega)) #(1 + 1/sig_sq)/2
  
  l <<- l1/sqrt(avg_prec)
  
  asymptote_out <<- squareDistAsymptote(l, d, avg_cov, avg_prec) # Asymptote for CRN and reflection
  asymptote_crn <<- asymptote_out$crn/avg_cov/d
  asymptote_refl <<- asymptote_out$refl/avg_cov/d
  
  logpi <<- function(x) { - 0.5 * sum((diag(Omega_chol) * x)^2)}
  gradlogpi <<- function(x) { - diag(Omega) * x}
  
  # MCMC starting values
  x <<- as.vector(Sigma_chol%*%rnorm(d))
  y <<- as.vector(Sigma_chol%*%rnorm(d))
}


####
# Setup: stationary-optimal step size
####
l1 <- 2.38

# AR(1) target
ar1_setup(l1, seed)

plot_df1 <- getPlotDf(x,y,x0,y0,v0,l,avg_cov,d,tmax,length_plot,logpi,gradlogpi)
plot_df1$l <- "Step size (a)"
plot_df1$target <- "Target (i)"

plot_df_asympt1 <- data.frame("Coupling" = c("CRN", "GCRN", "Reflection", "Optimal Markovian"), "scaled.squaredist" = c(asymptote_crn, 0, asymptote_refl, 0))
plot_df_asympt1$Linetype <- "Predicted asymptote"
plot_df_asympt1$l <- "Step size (a)"
plot_df_asympt1$target <- "Target (i)"


# Eigenvalues sampled from Chi-square
chisq_setup(l1, seed)

plot_df2 <- getPlotDf(x,y,x0,y0,v0,l,avg_cov,d,tmax,length_plot,logpi,gradlogpi)
plot_df2$l <- "Step size (a)"
plot_df2$target <- "Target (ii)"

plot_df_asympt2 <- data.frame("Coupling" = c("CRN", "GCRN", "Reflection", "Optimal Markovian"), "scaled.squaredist" = c(asymptote_crn, 0, asymptote_refl, 0))
plot_df_asympt2$Linetype <- "Predicted asymptote"
plot_df_asympt2$l <- "Step size (a)"
plot_df_asympt2$target <- "Target (ii)"


# Diag(1, 24, 1, 24,...)
alternating_setup(l1, seed)

plot_df3 <- getPlotDf(x,y,x0,y0,v0,l,avg_cov,d,tmax,length_plot,logpi,gradlogpi)
plot_df3$l <- "Step size (a)"
plot_df3$target <- "Target (iii)"

plot_df_asympt3 <- data.frame("Coupling" = c("CRN", "GCRN", "Reflection", "Optimal Markovian"), "scaled.squaredist" = c(asymptote_crn, 0, asymptote_refl, 0))
plot_df_asympt3$Linetype <- "Predicted asymptote"
plot_df_asympt3$l <- "Step size (a)"
plot_df_asympt3$target <- "Target (iii)"


# Gather data
plot_df_bigl <- rbind(plot_df1, plot_df2, plot_df3)
plot_df_bigl_asympt <- rbind(plot_df_asympt1, plot_df_asympt2, plot_df_asympt3)


####
# Setup: smaller step size, close to optimal for transient phase
####
l1 <- sqrt(2)

# AR(1) target
ar1_setup(l1, seed)

plot_df1 <- getPlotDf(x,y,x0,y0,v0,l,avg_cov,d,tmax,length_plot,logpi,gradlogpi)
plot_df1$l <- "Step size (b)"
plot_df1$target <- "Target (i)"

plot_df_asympt1 <- data.frame("Coupling" = c("CRN", "GCRN", "Reflection", "Optimal Markovian"), "scaled.squaredist" = c(asymptote_crn, 0, asymptote_refl, 0))
plot_df_asympt1$Linetype <- "Predicted asymptote"
plot_df_asympt1$l <- "Step size (b)"
plot_df_asympt1$target <- "Target (i)"


# Eigenvalues sampled from Chi-square
chisq_setup(l1, seed)

plot_df2 <- getPlotDf(x,y,x0,y0,v0,l,avg_cov,d,tmax,length_plot,logpi,gradlogpi)
plot_df2$l <- "Step size (b)"
plot_df2$target <- "Target (ii)"

plot_df_asympt2 <- data.frame("Coupling" = c("CRN", "GCRN", "Reflection", "Optimal Markovian"), "scaled.squaredist" = c(asymptote_crn, 0, asymptote_refl, 0))
plot_df_asympt2$Linetype <- "Predicted asymptote"
plot_df_asympt2$l <- "Step size (b)"
plot_df_asympt2$target <- "Target (ii)"


# Diagonal with alternating eigenvalues, e.g. (1,2,1,2,...)
alternating_setup(l1, seed)

plot_df3 <- getPlotDf(x,y,x0,y0,v0,l,avg_cov,d,tmax,length_plot,logpi,gradlogpi)
plot_df3$l <- "Step size (b)"
plot_df3$target <- "Target (iii)"

plot_df_asympt3 <- data.frame("Coupling" = c("CRN", "GCRN", "Reflection", "Optimal Markovian"), "scaled.squaredist" = c(asymptote_crn, 0, asymptote_refl, 0))
plot_df_asympt3$Linetype <- "Predicted asymptote"
plot_df_asympt3$l <- "Step size (b)"
plot_df_asympt3$target <- "Target (iii)"


# Gather data
plot_df_smalll <- rbind(plot_df1, plot_df2, plot_df3)
plot_df_smalll_asympt <- rbind(plot_df_asympt1, plot_df_asympt2, plot_df_asympt3)




##########
# Plot ###

plot_df <- rbind(plot_df_bigl, plot_df_smalll)
plot_df_asympt <- rbind(plot_df_bigl_asympt, plot_df_smalll_asympt)

save(plot_df, plot_df_asympt, file = "elliptical.RData")

load("elliptical.RData")

out_plot <- 
  ggplot(plot_df, aes(x = iter.d, y = scaled.squaredist)) +
  geom_line(aes(color = Coupling, linetype = Linetype), size = 0.75) +
  facet_grid(cols = vars(target), rows = vars(l), scales = "free_x") +
  geom_hline(data = plot_df_asympt, aes(yintercept = scaled.squaredist, 
                                        color = Coupling, 
                                        linetype = Linetype), size = 0.75) +
  theme_bw() +
  theme(legend.position = "bottom", 
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text=element_text(size = 11)) +
  xlab(paste("Iteration", TeX("$t/d$"))) +
  ylab("Scaled squared distance") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 1), 
         linetype = guide_legend(override.aes = list(linewidth = 0.75), order = 2))
out_plot 

ggsave(filename = "elliptical.pdf", 
       plot = out_plot,
       device = "pdf",  width = 30, height = 16, units = "cm", bg = "transparent")



############
# Contractivity and long-time behaviour in finite dimension

library(rwmcouplings)
library(foreach)
library(doParallel)
library(ggplot2)

algs <- c(RWM_gcrn)
cpls <- c("GCRN")
ds <- c(2,4,6,8,10,100,1000,10000)
lambda <- 2.38
iters <- ds * 1000
reps <- 3

# Omegas <- lapply(ds, function(d) rchisq(d,df=3)) # Precision matrices
Omegas <- lapply(ds, function(d) rep(c(1, 1/24),ceiling(d/2))[1:d])  # Precision matrices

cl <- makeCluster(3)
registerDoParallel(cl)
contractivity_df <-
  foreach(r = 1:reps, .combine = "rbind") %:%
  foreach(alg = algs, cpl = cpls, .combine = "rbind") %:%
  foreach(d = ds, Omega = Omegas, iter = iters, .combine = "rbind", .packages = "rwmcouplings") %dopar% {

    thin <- ceiling(d/2); thinned_iter <- seq(0,iter,thin)
    h <- lambda / sqrt(sum(Omega))
    
    x <- rnorm(d, sd = sqrt(1/Omega))
    y <- rnorm(d, sd = sqrt(1/Omega))
    
    elliptical <- function(x){-0.5 * sum(Omega * x^2)}
    elliptical_grad <- function(x){- Omega * x}
    
    squaredist <- alg(x,y,h,iter,elliptical,elliptical_grad)$squaredist[thinned_iter + 1]
    
    return(data.frame("iter" = thinned_iter, "squaredist" = squaredist, "h" = h^2, "d" = d, "coupling" = cpl, "replicate" = r))
  }
stopCluster(cl)

contractivity_df$replicate <- as.factor(contractivity_df$replicate)

ggplot(contractivity_df, aes(x = iter/d, y = squaredist, color = replicate)) +
  facet_grid(coupling~d, scales = "free", labeller = label_both) +
  geom_line() +
  geom_hline(aes(yintercept = h^2), color = "black", linetype = "dashed") +
  geom_text(aes(x=0,y=h^2,label="y = h^2"), stat = "unique", vjust="inward",hjust="inward",color = "black") +
  scale_y_log10() +
  theme_bw() +
  labs(x  = "Iteration (times d)", 
       y = "Squared distance between chains", 
       title = "Contractivity in finite dimension",
       subtitle = "Elliptical Gaussian target")
