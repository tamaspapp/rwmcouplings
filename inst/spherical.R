
library(ggplot2)
library(ggh4x)
library(rwmcouplings)
library(reshape2)
library(latex2exp)
set.seed(12345)

####
# Preliminary functions
####

####
# Do MCMC, compute ODE limits, save as long format dataframe for plotting
getPlotDf <- function(x,y,x0,y0,v0,l,d,tmax,length_plot,asymptote_crn) {
  
  h <- l / sqrt(d)
  
  ode_times <- seq(0, tmax, length.out = length_plot)
  
  iter <- tmax * d
  thinned_mcmc <- 1 + floor(seq(0, iter, length.out = length_plot))
  
  
  # GCRN ####
  out_mcmc <- RWM_gcrn(x,y,h,iter,spherical,spherical_grad)
  out_ode  <- ode_gcrn(x0,y0,v0,l,ode_times)
  
  squaredist_mcmc <- out_mcmc$squaredist[thinned_mcmc]
  squaredist_ode <- out_ode$xsq + out_ode$ysq - 2 * out_ode$innerprod
  ode_asympt<- rep(NA, length_plot)
  
  
  df_gcrn <- data.frame("squaredist/d" = squaredist_mcmc/d, "ode" = squaredist_ode, "ode.asymptote" = ode_asympt, "iter/d" = ode_times, "type" = "GCRN")
  
  # CRN ####
  out_mcmc <- RWM_crn(x,y,h,iter,spherical,spherical_grad)
  out_ode  <- ode_crn(x0,y0,v0,l,ode_times)
  
  # out_ode$innerprod / sqrt(out_ode$xsq * out_ode$ysq)
  
  squaredist_mcmc <- out_mcmc$squaredist[thinned_mcmc]
  squaredist_ode <- out_ode$xsq + out_ode$ysq - 2 * out_ode$innerprod
  ode_asympt <- rep(NA, length_plot)#rep(asymptote_crn, length_plot)
  
  df_crn <- data.frame("squaredist/d" = squaredist_mcmc/d, "ode" = squaredist_ode, "ode.asymptote" = ode_asympt, "iter/d" = ode_times, "type" = "CRN")
  
  # Reflection ####
  out_mcmc <- RWM_refl(x,y,h,iter,spherical,spherical_grad)
  out_ode  <- ode_refl(x0,y0,v0,l,ode_times)
  
  squaredist_mcmc <- out_mcmc$squaredist[thinned_mcmc]
  squaredist_ode <- out_ode$xsq + out_ode$ysq - 2 * out_ode$innerprod
  ode_asympt<- rep(NA, length_plot)
  
  df_refl <- data.frame("squaredist/d" = squaredist_mcmc/d, "ode" = squaredist_ode, "ode.asymptote" = ode_asympt, "iter/d" = ode_times, "type" = "Reflection")
  
  # Optimal ####
  out_ode <- ode_optimal(x0,y0,v0,l,ode_times)
  
  squaredist_mcmc <- rep(NA, length(thinned_mcmc))
  squaredist_ode <- out_ode$xsq + out_ode$ysq - 2 * out_ode$innerprod
  ode_asympt<- rep(NA, length_plot)
  
  df_optimal <- data.frame("squaredist/d" = squaredist_mcmc/d, "ode" = squaredist_ode, "ode.asymptote" = ode_asympt, "iter/d" = ode_times, "type" = "Optimal Markovian")
  
  # Join dataframes into a single one, for plotting ####
  df_ <- rbind(df_gcrn,df_crn, df_refl, df_optimal)
  df_ <- melt(df_, id.var = c("iter.d", "type"), variable.name = "ode")
  levels(df_$ode) <- c("MCMC", "ODE limit", "ODE asymptote")
  names(df_) <- c("iter.d", "Coupling", "Linetype", "squaredist.d")
  df_ <- na.omit(df_)
  return(df_)
}

d <- 1000
length_plot <- 1001
tmax <- 10

####
# Setup: stationary-optimal step size
####
l <- 2.38
asymptote_out <- squareDistAsymptote(l, d, 1, 1) # Asymptote for CRN and reflection
asymptote_crn <- asymptote_out$crn/d

####
# Independent start from stationary
####

# Parameters and starting values ###

rho <- 0

# ODE starting values
x0 <- 1
y0 <- 1
v0 <- sqrt(x0*y0)*rho

z <- rnorm(d)
# MCMC starting values
x <- sqrt(x0) * z
y <- sqrt(y0) * (rho * z + sqrt(1-rho^2) * rnorm(d))

# Compute  ###
plot_df1 <- getPlotDf(x,y,x0,y0,v0,l,d,tmax,length_plot, asymptote_crn)


####
# Positively correlated start, stationary + stationary
####

# Parameters and starting values ###
rho <- 0.9

# ODE starting values
x0 <- 1
y0 <- 1
v0 <- sqrt(x0*y0)*rho

# MCMC starting values
z <- rnorm(d)

x <- sqrt(x0) * z
y <- sqrt(y0) * (rho * z + sqrt(1-rho^2) * rnorm(d))

# Compute ###
plot_df2 <- getPlotDf(x,y,x0,y0,v0,l,d,tmax,length_plot,asymptote_crn)


####
# Independent start, underdispersed + overdispersed
####

rho <- 0

# ODE starting values
x0 <- 2.25
y0 <- 0.4
v0 <- sqrt(x0*y0)*rho

# MCMC starting values
z <- rnorm(d)

x <- sqrt(x0) * z
y <- sqrt(y0) * (rho * z + sqrt(1-rho^2) * rnorm(d))

# Compute ###
plot_df3 <- getPlotDf(x,y,x0,y0,v0,l,d,tmax,length_plot,asymptote_crn)


####
# Negatively correlated start, underdispersed + underdispersed
####

rho <- -0.5

# ODE starting values
x0 <- 0.01
y0 <- 0.25
v0 <- sqrt(x0*y0)*rho

# MCMC starting values
z <- rnorm(d)

x <- sqrt(x0) * z
y <- sqrt(y0) * (rho * z + sqrt(1-rho^2) * rnorm(d))

# Compute ###
plot_df4 <- getPlotDf(x,y,x0,y0,v0,l,d,tmax,length_plot,asymptote_crn)

# Gather data ####
plot_df1$start <- "Start (i)" 
plot_df2$start <- "Start (ii)"
plot_df3$start <- "Start (iii)"
plot_df4$start <- "Start (iv)"
plot_df_bigl <- rbind(plot_df1, plot_df2, plot_df3, plot_df4)
plot_df_bigl$l <- "Step size (a)"


plot_df_bigl_asympt <- data.frame(c("CRN", "GCRN", "Reflection"), c(asymptote_crn, 0, 0))
names(plot_df_bigl_asympt) <- c("Coupling", "squaredist.d")
plot_df_bigl_asympt$Linetype <- "ODE asymptote"
plot_df_bigl_asympt$l <- "Step size (a)"


####
# Setup: transient-optimal step size
####
l <- sqrt(2)
asymptote_out <- squareDistAsymptote(l, d, 1, 1) # Asymptote for CRN and reflection
asymptote_crn <- asymptote_out$crn/d

####
# Independent start from stationary
####

# Parameters and starting values ###

rho <- 0

# ODE starting values
x0 <- 1
y0 <- 1
v0 <- sqrt(x0*y0)*rho

z <- rnorm(d)
# MCMC starting values
x <- sqrt(x0) * z
y <- sqrt(y0) * (rho * z + sqrt(1-rho^2) * rnorm(d))

# Compute  ###
plot_df1 <- getPlotDf(x,y,x0,y0,v0,l,d,tmax,length_plot, asymptote_crn)


####
# Positively correlated start, stationary + stationary
####

# Parameters and starting values ###
rho <- 0.9

# ODE starting values
x0 <- 1
y0 <- 1
v0 <- sqrt(x0*y0)*rho

# MCMC starting values
z <- rnorm(d)

x <- sqrt(x0) * z
y <- sqrt(y0) * (rho * z + sqrt(1-rho^2) * rnorm(d))

# Compute ###
plot_df2 <- getPlotDf(x,y,x0,y0,v0,l,d,tmax,length_plot,asymptote_crn)


####
# Independent start, underdispersed + overdispersed
####

rho <- 0

# ODE starting values
x0 <- 1.5
y0 <- 0.5
v0 <- sqrt(x0*y0)*rho

# MCMC starting values
z <- rnorm(d)

x <- sqrt(x0) * z
y <- sqrt(y0) * (rho * z + sqrt(1-rho^2) * rnorm(d))

# Compute ###
plot_df3 <- getPlotDf(x,y,x0,y0,v0,l,d,tmax,length_plot,asymptote_crn)


####
# Negatively correlated start, underdispersed + underdispersed
####

rho <- -0.5

# ODE starting values
x0 <- 0.01
y0 <- 0.25
v0 <- sqrt(x0*y0)*rho

# MCMC starting values
z <- rnorm(d)

x <- sqrt(x0) * z
y <- sqrt(y0) * (rho * z + sqrt(1-rho^2) * rnorm(d))

# Compute ###
plot_df4 <- getPlotDf(x,y,x0,y0,v0,l,d,tmax,length_plot,asymptote_crn)


# Gather data ####
plot_df1$start <- "Start (i)" 
plot_df2$start <- "Start (ii)"
plot_df3$start <- "Start (iii)"
plot_df4$start <- "Start (iv)"
plot_df_smalll <- rbind(plot_df1, plot_df2, plot_df3, plot_df4)
plot_df_smalll$l <- "Step size (b)"


plot_df_smalll_asympt<- data.frame(c("CRN", "GCRN", "Reflection"), c(asymptote_crn, 0, 0))
names(plot_df_smalll_asympt) <- c("Coupling", "squaredist.d")
plot_df_smalll_asympt$Linetype <- "ODE asymptote"
plot_df_smalll_asympt$l <- "Step size (b)"

##########
# Plot ###

plot_df <- rbind(plot_df_bigl, plot_df_smalll)
plot_df_asympt <- rbind(plot_df_bigl_asympt, plot_df_smalll_asympt)

save(plot_df, plot_df_asympt, file = "spherical.RData")

load(file = "spherical.RData")

out_plot <- 
  ggplot(plot_df, aes(x = iter.d, y = squaredist.d)) +
  geom_line(aes(color = Coupling, linetype = Linetype), size = 0.75) +
  scale_linetype_manual(values=c(1,3,2)) +
  ggh4x::facet_grid2(cols = vars(start), rows = vars(l), scales = "free_y", independent = "y") +
  geom_hline(data = plot_df_asympt, aes(yintercept = squaredist.d, 
                                        color = Coupling, 
                                        linetype = Linetype), size = 0.75) +
  theme_bw() +
  theme(legend.position = "bottom", 
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text=element_text(size = 11)) +
  xlab(paste("Iteration", TeX("$t/d$"))) +
  ylab(TeX("$|| X_t - Y_t||^2/d$")) + guides(linetype = guide_legend(override.aes = list(size = 0.75)),
                                             color = guide_legend(override.aes = list(size = 2)))
out_plot 

ggsave(filename = "spherical.pdf", 
       plot = out_plot,
       device = "pdf",  width = 30, height = 16, units = "cm", bg = "transparent")








