library(rwmcouplings)
library(reshape2)
library(ggplot2)
library(latex2exp)

#################
# Elliptical case
#################

eccentricity <- c(1,1.2,2,10)
lambda       <- seq(0.01, 5, length.out = 100) # Natural step size parameter

# # Acceptance rates
# 2 * pnorm( - lambda/2)

# ESJD data frame for plotting
esjd_df <- data.frame("lambda" = lambda, 
                      "ESJD" = (2 * lambda^2)*pnorm(-lambda/2), 
                      "col" = "ESJD")


getAsymptotes <- function(lambda, eccentricity) {
  refl <- matrix(NA, nrow = length(lambda), ncol = length(eccentricity))
  crn <- rep(NA, length(lambda))
  
  for(i in seq_along(lambda)) {
    for(j in seq_along(eccentricity)) {
      out <- scaledSquaredistAsymptote(lambda[i],eccentricity[j])
      
      if(j == 1) {crn[i] <- out$crn}
      refl[i,j] <- out$refl
    }
  }
  
  # Create plotting dataframes
  crn_df <- data.frame(crn) # CRN
  crn_df[,2] <- lambda
  names(crn_df) <- c("asymptote", "lambda")
  crn_df$col <- "CRN"
   
  refl_df <- data.frame(refl) # Reflection
  refl_df$lambda <- lambda
  refl_df <- melt(refl_df, id.vars = "lambda")
  levels(refl_df$variable) <- eccentricity
  refl_df$variable <- as.numeric(levels(refl_df$variable))[refl_df$variable]
  names(refl_df) <- c("lambda", "eccentricity", "asymptote")
  refl_df$col <- "Reflection"

  return(list("crn" = crn_df, "refl" = refl_df))
}

out <- getAsymptotes(lambda, eccentricity)

crn_df_plt  <- out$crn
refl_df_plt <- out$refl

fixe_varyl <-
  ggplot(refl_df_plt, aes(x = lambda, y = asymptote, linetype = as.factor(eccentricity), color = col)) +
  geom_line(data = esjd_df, aes(x = lambda, y = ESJD, linetype = as.factor(1), color = col), lwd = 0.75) +
  geom_line(data = crn_df_plt, aes(x = lambda, y = asymptote, linetype = as.factor(1), color = col), lwd = 0.75) +
  geom_line(data = refl_df_plt, lwd = 0.75) +
  #coord_cartesian(xlim = c(0, 5), ylim = c(0,2)) +
  labs(linetype = "Eccentricity", color = "Coupling") +
  xlab(TeX("$\\lambda$ (step size)")) + ylab(TeX("$s_{coup}^*$ (scaled square distance)")) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title.y = element_text(size = 14),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text=element_text(size = 12)) + 
  guides(linetype = guide_legend(order = 2, override.aes = list(size = 0.75)),
         color = guide_legend(order = 1, override.aes = list(size = 2))) +
  annotate("text", x=2.38, y=1.38, label= "ESJD",colour = "#00BA38") + 
  scale_color_discrete(breaks=c("CRN", "Reflection", NA)) +
  scale_linetype_discrete(breaks = as.character(eccentricity))
fixe_varyl

ggsave(filename = "elliptasymptote.pdf",
       plot = fixe_varyl,
       device = "pdf",  width = 25, height = 14, units = "cm", bg = "transparent")


