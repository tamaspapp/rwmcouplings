library(ggplot2)
library(latex2exp)
library(dplyr)
library(ggpubr)
library(ggh4x)

load("binreg-prelim.RData")

# Label chains that didn't meet by maximum number of iterations run for
chainsDidntMeet <- function(df, iter_max) { df$tau[df$tau == -1] <- iter_max; return(df)}

# Full variance matrix ####
load("binreg-ablation-fullvar.RData")

mala_abl       <- chainsDidntMeet(mala_abl,iter_max)
rwm_abl_gcrefl <- chainsDidntMeet(rwm_abl_gcrefl, iter_max)
rwm_abl_gcrn   <- chainsDidntMeet(rwm_abl_gcrn, iter_max)

rwm_df <- rbind(cbind(rwm_abl_gcrefl, "algorithm" = "RWM + GCRefl"),
                cbind(rwm_abl_gcrn,   "algorithm" = "RWM + GCRN"))

full_df <- rbind(cbind(rwm_abl_gcrefl, "algorithm" = "RWM + GCRefl"),
                 cbind(rwm_abl_gcrn,   "algorithm" = "RWM + GCRN"),
                 cbind(mala_abl,   "algorithm" = "MALA + CRN"))
full_df <- cbind(full_df, "prec" = "Full preconditioner")

# Diagonal variance matrix ####
load("binreg-ablation-diagvar.RData")

mala_abl       <- chainsDidntMeet(mala_abl,iter_max)
rwm_abl_gcrefl <- chainsDidntMeet(rwm_abl_gcrefl, iter_max)
rwm_abl_gcrn   <- chainsDidntMeet(rwm_abl_gcrn, iter_max)

rwm_df <- rbind(cbind(rwm_abl_gcrefl, "algorithm" = "RWM + GCRefl"),
                cbind(rwm_abl_gcrn,   "algorithm" = "RWM + GCRN"))

diag_df <- rbind(cbind(rwm_abl_gcrefl, "algorithm" = "RWM + GCRefl"),
                 cbind(rwm_abl_gcrn,   "algorithm" = "RWM + GCRN"),
                 cbind(mala_abl,   "algorithm" = "MALA + CRN"))
diag_df <- cbind(diag_df, "prec" = "Diagonal preconditioner")
#####

###### Joint plot
joint_df <- rbind(diag_df, full_df)

plot_params <- list(
  # Means as bots with black fill
  stat_summary(fun = mean, position = position_dodge(0.75),
               geom = "point", shape = 21, size = 1.5, show.legend = FALSE, fill = "black"), 
  # Axis and legend, where group is specified by color
  labs(x = TeX("Step size $h$"), y = "Meeting time", color = TeX("Threshold $\\delta$")), 
  coord_trans(y = "log10"), scale_y_continuous(breaks = scales::breaks_log(7), minor_breaks = NULL), 
  theme_bw(),
  # Legend on bottom and one one row
  theme(legend.position = "bottom",
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text=element_text(size = 11)), guides(colour = guide_legend(nrow = 1))
)

joint_plot <-
  ggplot(joint_df, aes(color = factor(thresh), y = tau, x = factor(l))) +
  geom_boxplot() +
  plot_params +
  facet_grid2(prec~algorithm, scales = "free", independent = "x") +
  geom_hline(data = joint_df %>% filter(prec == "Diagonal preconditioner"),
             aes(yintercept = iter_max), col = "black", linetype = "dashed")
joint_plot

ggsave(filename = "binreg-ablation.pdf",
       plot = joint_plot,
       device = "pdf",  width = 24, height = 16, units = "cm", bg = "transparent")
