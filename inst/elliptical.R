
#####
# Plots for squared distance between chains over the average covariance of the target, elliptical case
#####

library(Matrix)
library(ggplot2)
library(ggh4x)
library(rwmcouplings)
library(reshape2)
library(latex2exp)
set.seed(12345)

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

####
# Do MCMC, compute ODE limits, save as long format dataframe for plotting
getPlotDf <- function(x,y,x0,y0,v0,l,avg_cov,d,tmax,length_plot,logpi,gradlogpi) {
  
  h <- l / sqrt(d)
  
  ode_times <- seq(0, tmax, length.out = length_plot)
  
  iter <- tmax * d
  thinned_mcmc <- 1 + floor(seq(0, iter, length.out = length_plot))
  
  # GCRN ####
  out_mcmc <- RWM_gcrn(x,y,h,iter,logpi,gradlogpi)
  squaredist_mcmc <- out_mcmc$squaredist[thinned_mcmc]
  ode_asympt<- rep(NA, length_plot)
  df_gcrn <- data.frame("scaled.squaredist" = squaredist_mcmc/d/avg_cov, "ode.asymptote" = ode_asympt, "iter/d" = ode_times, "Coupling" = "GCRN")
  
  # CRN ####
  out_mcmc <- RWM_crn(x,y,h,iter,logpi,gradlogpi)
  squaredist_mcmc <- out_mcmc$squaredist[thinned_mcmc]
  ode_asympt <- rep(NA, length_plot)
  df_crn <- data.frame("scaled.squaredist" = squaredist_mcmc/d/avg_cov,"ode.asymptote" = ode_asympt, "iter/d" = ode_times, "Coupling" = "CRN")
  
  # Reflection ####
  out_mcmc <- RWM_refl(x,y,h,iter,logpi,gradlogpi)
  squaredist_mcmc <- out_mcmc$squaredist[thinned_mcmc]
  ode_asympt<- rep(NA, length_plot)
  df_refl <- data.frame("scaled.squaredist" = squaredist_mcmc/d/avg_cov, "ode.asymptote" = ode_asympt, "iter/d" = ode_times, "Coupling" = "Reflection")
  
  # Join dataframes into a single one, for plotting ####
  df_ <- rbind(df_gcrn,df_crn, df_refl)
  df_ <- melt(df_, id.var = c("iter.d", "Coupling"), variable.name = "Linetype", value.name = "scaled.squaredist")
  levels(df_$Linetype) <- c("MCMC", "Predicted asymptote")
  df_ <- na.omit(df_)
  return(df_)
}

d <- 2000
length_plot <- 2001

####
# Setup: stationary-optimal step size
####
l1 <- 2.38

####
# AR(1) target
####
tmax <- 40
rho <- 0.5 # AR(1) correlation

avg_cov <- 1
avg_prec <- (1 + rho^2)/ (1-rho^2)
  
l <- l1/sqrt(avg_prec)

asymptote_out  <- squareDistAsymptote(l, d, avg_cov, avg_prec) # Asymptote for CRN and reflection
asymptote_crn  <- asymptote_out$crn/avg_cov/d
asymptote_refl <- asymptote_out$refl/avg_cov/d

Omega <- KMSMatrixInverse(rho, d)
Omega_chol <- chol(Omega)

logpi <- function(x) { - 0.5 * sum(as.vector((Omega_chol %*% x))^2)}
gradlogpi <- function(x) { - as.vector(Omega %*% x)}

# MCMC starting values
Sigma <- KMSMatrix(rho, d)
Sigma_cholLT <- t(chol(Sigma)) 

x <- as.vector(Sigma_cholLT%*%rnorm(d))
y <- as.vector(Sigma_cholLT%*%rnorm(d))

# Compute  ###
plot_df1 <- getPlotDf(x,y,x0,y0,v0,l,avg_cov,d,tmax,length_plot,logpi,gradlogpi)
plot_df1$l <- "Step size (a)"
plot_df1$target <- "Target (i)"

plot_df_asympt1 <- data.frame("Coupling" = c("CRN", "GCRN", "Reflection"), "scaled.squaredist" = c(asymptote_crn, 0, asymptote_refl))
plot_df_asympt1$Linetype <- "Predicted asymptote"
plot_df_asympt1$l <- "Step size (a)"
plot_df_asympt1$target <- "Target (i)"


####
# Eigenvalues sampled from Chi-square
####
tmax <- 80
nu <- 3 # Chi-square degrees of freedom, assumed nu>2
Sigma_diagonal_chisq <- rchisq(d, df = nu)

avg_cov <- nu
avg_prec <- 1/(nu - 2)

l <- l1/sqrt(avg_prec)
  
asymptote_out  <- squareDistAsymptote(l, d, avg_cov, avg_prec) # Asymptote for CRN and reflection
asymptote_crn  <- asymptote_out$crn/avg_cov/d
asymptote_refl <- asymptote_out$refl/avg_cov/d

Sigma <- Matrix::Diagonal(x = Sigma_diagonal_chisq)
Sigma_chol <- Sigma^0.5
Omega <- Matrix::Diagonal(x = 1/Sigma_diagonal_chisq)
Omega_chol <- Omega^0.5

logpi <- function(x) { - 0.5 * sum((diag(Omega_chol) * x)^2)}
gradlogpi <- function(x) { - diag(Omega) * x}
  
# MCMC starting values
x <- as.vector(Sigma_chol%*%rnorm(d))
y <- as.vector(Sigma_chol%*%rnorm(d))

# Compute  ###
plot_df2 <- getPlotDf(x,y,x0,y0,v0,l,avg_cov,d,tmax,length_plot,logpi,gradlogpi)
plot_df2$l <- "Step size (a)"
plot_df2$target <- "Target (ii)"

plot_df_asympt2 <- data.frame("Coupling" = c("CRN", "GCRN", "Reflection"), "scaled.squaredist" = c(asymptote_crn, 0, asymptote_refl))
plot_df_asympt2$Linetype <- "Predicted asymptote"
plot_df_asympt2$l <- "Step size (a)"
plot_df_asympt2$target <- "Target (ii)"



####
# Diagonal with repeating eigenvalues, e.g. (1,2,1,2,...)
####
tmax <- 160
sig_sq <- 24 # Eigenvalue that's not 1
Sigma_diagonal <- rep(c(1,sig_sq), ceiling(d/2))[1:d]

avg_cov  <- (1 + sig_sq)/2
avg_prec <- (1 + 1/sig_sq)/2

l <- l1/sqrt(avg_prec)

asymptote_out  <- squareDistAsymptote(l, d, avg_cov, avg_prec) # Asymptote for CRN and reflection
asymptote_crn  <- asymptote_out$crn/avg_cov/d
asymptote_refl <- asymptote_out$refl/avg_cov/d

Sigma <- Matrix::Diagonal(x = Sigma_diagonal)
Sigma_chol <- Sigma^0.5
Omega <- Matrix::Diagonal(x = 1/Sigma_diagonal)
Omega_chol <- Omega^0.5

logpi <- function(x) { - 0.5 * sum((diag(Omega_chol) * x)^2)}
gradlogpi <- function(x) { - diag(Omega) * x}

# MCMC starting values
x <- as.vector(Sigma_chol%*%rnorm(d))
y <- as.vector(Sigma_chol%*%rnorm(d))

# Compute  ###
plot_df3 <- getPlotDf(x,y,x0,y0,v0,l,avg_cov,d,tmax,length_plot,logpi,gradlogpi)
plot_df3$l <- "Step size (a)"
plot_df3$target <- "Target (iii)"

plot_df_asympt3 <- data.frame("Coupling" = c("CRN", "GCRN", "Reflection"), "scaled.squaredist" = c(asymptote_crn, 0, asymptote_refl))
plot_df_asympt3$Linetype <- "Predicted asymptote"
plot_df_asympt3$l <- "Step size (a)"
plot_df_asympt3$target <- "Target (iii)"


# Gather data ####
plot_df_bigl <- rbind(plot_df1, plot_df2, plot_df3)
plot_df_bigl_asympt <- rbind(plot_df_asympt1, plot_df_asympt2, plot_df_asympt3)




# 
# ####
# # Setup: smaller step size, close to optimal for transient phase
# ####
l1 <- sqrt(2)

####
# AR(1) target
####
tmax <- 40
rho <- 0.5 # AR(1) correlation

avg_cov <- 1
avg_prec <- (1 + rho^2)/ (1-rho^2)

l <- l1/sqrt(avg_prec)

asymptote_out  <- squareDistAsymptote(l, d, avg_cov, avg_prec) # Asymptote for CRN and reflection
asymptote_crn  <- asymptote_out$crn/avg_cov/d
asymptote_refl <- asymptote_out$refl/avg_cov/d

Omega <- KMSMatrixInverse(rho, d)
Omega_chol <- chol(Omega)
logpi <- function(x) { - 0.5 * sum(as.vector((Omega_chol %*% x))^2)}
gradlogpi <- function(x) { - as.vector(Omega %*% x)}

# MCMC starting values
x <- as.vector(solve(Omega_chol, rnorm(d)))
y <- as.vector(solve(Omega_chol, rnorm(d)))

# Compute  ###
plot_df1 <- getPlotDf(x,y,x0,y0,v0,l,avg_cov,d,tmax,length_plot,logpi,gradlogpi)
plot_df1$l <- "Step size (b)"
plot_df1$target <- "Target (i)"

plot_df_asympt1 <- data.frame("Coupling" = c("CRN", "GCRN", "Reflection"), "scaled.squaredist" = c(asymptote_crn, 0, asymptote_refl))
plot_df_asympt1$Linetype <- "Predicted asymptote"
plot_df_asympt1$l <- "Step size (b)"
plot_df_asympt1$target <- "Target (i)"


####
# Eigenvalues sampled from Chi-square
####
tmax <- 80
# nu <- 3 # Chi-square degrees of freedom, assumed nu>2
# Sigma_diagonal <- rchisq(d, df = nu)

avg_cov <- nu
avg_prec <- 1/(nu - 2)

l <- l1/sqrt(avg_prec)

asymptote_out  <- squareDistAsymptote(l, d, avg_cov, avg_prec) # Asymptote for CRN and reflection
asymptote_crn  <- asymptote_out$crn/avg_cov/d
asymptote_refl <- asymptote_out$refl/avg_cov/d

Sigma <- Matrix::Diagonal(x = Sigma_diagonal_chisq)
Sigma_chol <- Sigma^0.5
Omega <- Matrix::Diagonal(x = 1/Sigma_diagonal_chisq)
Omega_chol <- Omega^0.5

logpi <- function(x) { - 0.5 * sum((diag(Omega_chol) * x)^2)}
gradlogpi <- function(x) { - diag(Omega) * x}

# MCMC starting values
x <- as.vector(Sigma_chol%*%rnorm(d))
y <- as.vector(Sigma_chol%*%rnorm(d))

# Compute  ###
plot_df2 <- getPlotDf(x,y,x0,y0,v0,l,avg_cov,d,tmax,length_plot,logpi,gradlogpi)
plot_df2$l <- "Step size (b)"
plot_df2$target <- "Target (ii)"

plot_df_asympt2 <- data.frame("Coupling" = c("CRN", "GCRN", "Reflection"), "scaled.squaredist" = c(asymptote_crn, 0, asymptote_refl))
plot_df_asympt2$Linetype <- "Predicted asymptote"
plot_df_asympt2$l <- "Step size (b)"
plot_df_asympt2$target <- "Target (ii)"



####
# Diagonal with repeating eigenvalues, e.g. (1,2,1,2,...)
####
tmax <- 160
sig_sq <- 24 # Eigenvalue that's not 1
Sigma_diagonal <- rep(c(1,sig_sq), ceiling(d/2))[1:d]

avg_cov  <- (1 + sig_sq)/2
avg_prec <- (1 + 1/sig_sq)/2

l <- l1/sqrt(avg_prec)

asymptote_out  <- squareDistAsymptote(l, d, avg_cov, avg_prec) # Asymptote for CRN and reflection
asymptote_crn  <- asymptote_out$crn/avg_cov/d
asymptote_refl <- asymptote_out$refl/avg_cov/d

Sigma <- Matrix::Diagonal(x = Sigma_diagonal)
Sigma_chol <- Sigma^0.5
Omega <- Matrix::Diagonal(x = 1/Sigma_diagonal)
Omega_chol <- Omega^0.5

logpi <- function(x) { - 0.5 * sum((diag(Omega_chol) * x)^2)}
gradlogpi <- function(x) { - diag(Omega) * x}

# MCMC starting values
x <- as.vector(Sigma_chol%*%rnorm(d))
y <- as.vector(Sigma_chol%*%rnorm(d))

# Compute  ###
plot_df3 <- getPlotDf(x,y,x0,y0,v0,l,avg_cov,d,tmax,length_plot,logpi,gradlogpi)
plot_df3$l <- "Step size (b)"
plot_df3$target <- "Target (iii)"

plot_df_asympt3 <- data.frame("Coupling" = c("CRN", "GCRN", "Reflection"), "scaled.squaredist" = c(asymptote_crn, 0, asymptote_refl))
plot_df_asympt3$Linetype <- "Predicted asymptote"
plot_df_asympt3$l <- "Step size (b)"
plot_df_asympt3$target <- "Target (iii)"



# Gather data ####
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
  ylab("Scaled squared distance") + guides(linetype = guide_legend(override.aes = list(size = 0.75)),
                                           color = guide_legend(override.aes = list(size = 2)))
out_plot 

ggsave(filename = "elliptical.pdf", 
       plot = out_plot,
       device = "pdf",  width = 30, height = 16, units = "cm", bg = "transparent")




