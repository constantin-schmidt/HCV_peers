##  Plot of posterior of theta1 ##
##################################
rm(list = ls())

library(ggplot2)
library(bayesplot)
library(rstan)
library(tidyr)
library(dplyr)

################
##  Load data ##
################

# Posterior fit
fit_full_bayes <- readRDS(file = "./03_Data/ZZ_Temp/fit_full_bayes.rds")

####################
##  Posterior     ##
####################

theta2 <- rstan::extract(fit_full_bayes, pars = "theta2")$theta2 |>
  exp()
d <- density(theta2)

# Get quantiles for intervals
q80 <- quantile(theta2, probs = c(0.10, 0.90))
q95 <- quantile(theta2, probs = c(0.025, 0.975))

# Create data frame from density
plot_df <- data.frame(x = d$x, y = d$y)

# Create flags for interval inclusion
plot_df$in_95 <- plot_df$x >= q95[1] & plot_df$x <= q95[2]
plot_df$in_80 <- plot_df$x >= q80[1] & plot_df$x <= q80[2]


####################
##  Posterior plot ##
####################

# Plot
p <- ggplot(plot_df, aes(x = x, y = y)) +
  geom_line(color = "black") +
  
  # 95% interval fill
  geom_ribbon(data = subset(plot_df, in_95), aes(ymin = 0, ymax = y), fill = "#a6cee3", alpha = 0.7) +
  
  # 80% interval fill (darker)
  geom_ribbon(data = subset(plot_df, in_80), aes(ymin = 0, ymax = y), fill = "#1f78b4", alpha = 0.7) +
  
  theme_minimal() +
  labs(
    title = "(c)",
    x = NULL,
    y = "Density"
  ) 

####################
##  Save to file  ##
####################

ggsave(
  filename = "./05_ResultsAndDiagnostics/results_section_graphs/posterior_density_theta2.pdf",
  plot = p, 
  width = 4.5, 
  height = 4, 
  units = "in"
)
