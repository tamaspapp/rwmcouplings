library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggrepel)
library(ggh4x)

load("binreg-prelim.RData")
load("binreg-targetsample.RData")

get_summaries_for_plot <- function(cpl_asympvar, setting){
  time_rwm <- setting$time_rwm
  time_mala <- setting$time_mala 
  acc_rates <- setting$acc_rates
  
  #####
  # 1. Seconds per <number> effective samples (worst coordinate)
  #####
  
  ##
  # 1.1 Asymptotic variance
  ##
  eff_df <- group_by(cpl_asympvar, l, coord, algorithm) %>% 
    summarise(asympvar.mean = mean(asympvar), asympvar.stderr = sd(asympvar)/sqrt(R))
  
  ##
  # 1.2 Integrated autocorrelation time
  ##
  eff_df <- merge(eff_df, posterior[, c("coord","var","var.stderr")], sort = F, all.x = TRUE, suffixes = NULL)
  eff_df$iact.mean <- eff_df$asympvar.mean / eff_df$var
  eff_df$iact.stderr  <- eff_df$asympvar.stderr / eff_df$var
  
  
  ##
  # 1.3 Summarise: maximum of the average IACT, plus its standard deviation
  ##
  max_iact <- function(datafr){
    group_by(datafr, l, algorithm) %>%
      summarise(iact.max = max(iact.mean), iact.max.stderr = iact.stderr[which.max(iact.mean)])
  }
  min_eff_df <- max_iact(eff_df)

  #####
  # 2. Average meeting time and standard error
  #####
  avg_meet <- function(datafr){
    group_by(datafr, l, algorithm) %>%
      summarise(meet.mean = mean(tau), meet.stderr = sd(tau)/sqrt(R))
  }
  # avg_time <- function(datafr){
  #   group_by(datafr, l, algorithm) %>%
  #     summarise(time.mean = as.double(mean(time)), time.stderr = as.double(sd(time))/sqrt(R))
  # }
  meet_df <- avg_meet(cpl_meet)
  # time_df <- avg_time(cpl_meet)
  
  meet_plot <- merge(merge(min_eff_df, meet_df, suffixes = NULL), acc_rates, suffixes = NULL)
  
  # Now scale by time
  scale_iact <- function(df){
    scale <- ifelse(df$algorithm == "RWM", time_rwm, time_mala)
    df$iact.max <- df$iact.max * scale
    df$iact.max.stderr <- df$iact.max.stderr * scale
    return(df)
  }
  scale_meeting_time <- function(df){
    # Account for two chains, and extra gradient evaluations made by RWM coupling
    scale <- 2 * ifelse(df$algorithm == "RWM", time_rwm + df$acc * (time_mala - time_rwm), time_mala)
    df$meet.mean <- df$meet.mean * scale
    df$meet.stderr <- df$meet.stderr * scale
    return(df)
  }
  time_plot <- meet_plot
  time_plot <- scale_iact(time_plot)
  time_plot <- scale_meeting_time(time_plot)

  return(list("unscaled" = meet_plot,
              "scaled_by_time" = time_plot))
}
print_summaries <- function(plot_dfs, setting) {
  time_rwm <- setting$time_rwm
  time_mala <- setting$time_mala 

  ## Summaries ##
  # How much slower is MALA per iteration?
  print(paste0("MALA is a factor: ", 
               signif(as.numeric(time_mala) / as.numeric(time_rwm),2), 
               " slower than RWM per iteration."))
  
  # How much slower is RWM mixing?
  mixing <- group_by(plot_dfs$scaled_by_time, algorithm) %>% summarise(time = min(iact.max))
  print(paste0("RWM is a factor: ", 
               signif(as.double(mixing$time[mixing$algorithm == "RWM"]) / as.double(mixing$time[mixing$algorithm == "MALA"]),2), 
               " slower to mix than MALA."))
  
  # How much slower is the RWM coalescence?
  mintime <- group_by(plot_dfs$scaled_by_time, algorithm) %>% summarise(time = min(meet.mean))
  print(paste0("RWM is a factor: ", 
               signif(as.double(mintime$time[mintime$algorithm == "RWM"]) / as.double(mintime$time[mintime$algorithm == "MALA"]),2), 
               " slower to meet than MALA."))
}

# Full proposal variance ####
load("binreg-coupling-fullvar.RData")
plot_dfs_fullvar <- get_summaries_for_plot(cpl_asympvar, FullVar)
print_summaries(plot_dfs_fullvar, FullVar)
#####

# Diagonal proposal variance  ####
load("binreg-coupling-diagvar.RData")
plot_dfs_diagvar <- get_summaries_for_plot(cpl_asympvar, DiagVar)
print_summaries(plot_dfs_diagvar, DiagVar)
#####

walltime <- 
  rbind(
    cbind(plot_dfs_fullvar$scaled_by_time, "prec" = "Full preconditioner"),
    cbind(plot_dfs_diagvar$scaled_by_time, "prec" = "Diagonal preconditioner")
  )

unscaled <- 
  rbind(
    cbind(plot_dfs_fullvar$unscaled, "prec" = "Full preconditioner"),
    cbind(plot_dfs_diagvar$unscaled, "prec" = "Diagonal preconditioner")
  )

get_plot_by_iter <- function(df, repel_seed = 1){
  ggplot(df, 
         aes(color = algorithm, x = iact.max, y = meet.mean)) +
    geom_point() +
    geom_path() +
    geom_errorbar(aes(ymax = meet.mean + 2*meet.stderr, 
                      ymin = meet.mean - 2*meet.stderr)) +
    geom_errorbarh(aes(xmax = (iact.max + 2*iact.max.stderr), 
                       xmin = (iact.max - 2*iact.max.stderr)), height = 0.01) +
    geom_text_repel(aes(label = signif(acc,2)), size = 2.5, box.padding = 0.05,
                     max.overlaps = 20, seed = repel_seed, max.time = 5,
                     show.legend = F, color = "black") +
    labs(x = "Integrated autocorrelation time", y = "Average meeting time", color = "Algorithm") + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank(),
                axis.title = element_text(size = 14),
                strip.text.x = element_text(size = 14),
                strip.text.y = element_text(size = 14),
                legend.text=element_text(size = 11)) +
    annotation_logticks(base = 10)
}
get_plot_by_time <- function(df, nsample, repel_seed = 1){
  # X-axis below is "nsample" number of effective samples per second
  
  df$iact.max <- nsample * df$iact.max
  df$iact.max.stderr <- nsample * df$iact.max.stderr
  
  ggplot(df, 
         aes(color = algorithm, x = iact.max, y = meet.mean)) +
    geom_point() +
    geom_path() +
    geom_errorbar(aes(ymax = meet.mean + 2*meet.stderr, 
                      ymin = meet.mean - 2*meet.stderr)) +
    geom_errorbarh(aes(xmax = (iact.max + 2*iact.max.stderr), 
                       xmin = (iact.max - 2*iact.max.stderr)), height = 0.01) +
    geom_text_repel(aes(label = signif(acc,2)), size = 2.5, box.padding = 0.05,
                     max.overlaps = 20, seed = repel_seed, max.time = 5,
                     show.legend = F, color = "black") +
    labs(x = paste("Seconds per", nsample, "effective samples when at equilibrium"),
         y = "Average meeting time (in seconds)", 
         color = "Algorithm") + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank(),
          axis.title = element_text(size = 14),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 14),
          legend.text=element_text(size = 11))
}

plt_walltime <- 
  get_plot_by_time(walltime, 10, 1) + 
  facet_grid2(~prec, scales = "free", independent = "all") +
  scale_x_log10(breaks = scales::breaks_log(6), expand = c(0.2, 0.15)) +
  scale_y_log10(breaks = scales::breaks_log(6), expand = c(0.05, 0.05)) +
  annotation_logticks(base = 10)
plt_walltime

plt_iter <- 
  get_plot_by_iter(unscaled, 1) + 
  facet_grid2(~prec, scales = "free", independent = "all") +
  scale_x_log10(breaks = scales::breaks_log(7), expand = c(0.2, 0.15)) +
  scale_y_log10(breaks = scales::breaks_log(7), expand = c(0.1, 0.05)) +
  annotation_logticks(base = 10)
plt_iter


ggsave(filename = "binreg-time.pdf",
       plot = plt_walltime,
       device = "pdf", width = 24, height = 10, units = "cm", bg = "transparent")
ggsave(filename = "binreg-iter.pdf",
       plot = plt_iter,
       device = "pdf", width = 24, height = 10, units = "cm", bg = "transparent")
