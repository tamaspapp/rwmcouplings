######
# Plots for asymptotically optimal Markovian coupling
######

library(ggplot2)
library(latex2exp)
library(gridExtra)
# Auxiliary functions
f <- function(u,l,x) {
  pmin(1, exp(-l * sqrt(x) * qnorm(u) - 0.5  *l^2))
}
plot_grid <- function(us,l,x) {
  swap_thresh <- pnorm(-l/(2*x^0.5))
  us <- c(c(0, swap_thresh), us[us > swap_thresh ])
  us
}

# Generate data
l <- 2.38
unifs <- seq(0, 1, 0.01)

x1 <- 1
unifs1 <-  plot_grid(unifs,l,x1)

x2 <- 10
unifs2 <-  plot_grid(unifs,l,x2)

df1 <- data.frame ("x" = unifs1,  "y" = f(unifs1, l, x1), "chain" = "X-chain")
df2 <- data.frame ("x" = unifs2,  "y" = f(unifs2, l, x2), "chain" = "Y-chain")
df_plot <- rbind(df1, df2)

# Plots
x_plot <- 
  ggplot(df1, aes(x=x, y=y, group=chain, fill=chain)) +
  geom_line(size=.5) +
  geom_ribbon(aes(x = x, ymax = y),ymin = 0,alpha = 0.3) +
  theme_bw() + 
  facet_grid(cols = vars(chain)) +
  ylab(TeX("$U_x$")) +
  xlab(TeX("$\\tilde{U}_x$")) +
  scale_fill_manual(values = "red") +
  theme(legend.position = "none")
x_plot

y_plot <- 
  ggplot(df2, aes(x=x, y=y, group=chain, fill=chain)) +
  geom_line(size=.5) +
  geom_ribbon(aes(x = x, ymax = y),ymin = 0,alpha = 0.3) +
  theme_bw() + 
  facet_grid(cols = vars(chain)) +
  ylab(TeX("$U_y$")) +
  xlab(TeX("$\\tilde{U}_y$")) + 
  scale_fill_manual(values = "blue") +
  theme(legend.position = "none")
y_plot

out_plot <- grid.arrange(x_plot, y_plot, nrow = 1)

ggsave(filename = "optimalmarkovian.pdf", 
       plot = out_plot,
       device = "pdf",  width = 24, height = 10, units = "cm", bg = "transparent")

# Plot both overlapped
# out_plot <- 
#   ggplot(df_plot, aes(x=x, y=y, group=chain, fill=chain)) +
#   geom_line(size=.5) +
#   facet_grid(cols = vars(chain)) +
#   geom_ribbon(aes(x = x, ymax = y),ymin = 0,alpha = 0.3) +
#   theme(legend.position = "none") +
#   theme_bw()
# out_plot

