##  Graphs showing the share of the intervention effect during COVID-19  ##
###########################################################################
rm(list =ls())

library(ggplot2)
library(dplyr)
library(tidyverse)

################
##  Load data ##
################

load('./03_Data/ZZ_Temp/panel_het_blue.RData')
load('./03_Data/ZZ_Temp/posterior_predictions_copula_list_full_bayes.RData')
list_assign_model <- posterior_predictions_copula_list

#############################################
##  Define a function to calculate effects ##
#############################################

compute_share_tau_c <- function(posterior_list, panel_data_Y, panel_data_D) {
  treat_array_list <- list()
  dimensions <- dim(posterior_list[[1]])
  
  nIter <- dimensions[1]
  nUnits <- dimensions[2]
  nTimes <- dimensions[3]
  
  Y_array <- array(rep(panel_data_Y, nIter), dim = c(nUnits, nTimes, nIter)) |>
    aperm(c(3, 1, 2))
  D_array <- array(rep(panel_data_D, nIter), dim = c(nUnits, nTimes, nIter)) |>
    aperm(c(3, 1, 2))
  
  for (i in seq_along(posterior_list)) {
    Y_0 <- posterior_list[[i]]
    treat_array_temp <- Y_array - Y_0
    treat_array_list[[i]] <- treat_array_temp * D_array
  }
  
  names(treat_array_list) <- names(posterior_list)
  
  cum_effects_df <- data.frame(matrix(ncol = length(treat_array_list), nrow = dim(treat_array_list[[1]])[1]))
  names(cum_effects_df) <- names(treat_array_list)
  
  for (i in seq_along(treat_array_list)) {
    cum_effects_df[[i]] <- apply(treat_array_list[[i]]
                                 , 1, function(x) sum(x*panel_list$LDtreat)/sum(x))
  }
  
  return(cum_effects_df)
}

#############################################
##  Compute share tau_c in overall effect  ##
#############################################

share_tau_c_in_tau_df <- compute_share_tau_c(list_assign_model, 
                                         panel_list$Y,
                                         panel_list$D)

############
##  Plot  ##
############

plot_long_df <- share_tau_c_in_tau_df |>
  select("rho=1","rho=0") |>
  pivot_longer(cols = everything(),
               names_to = "rho", values_to = "share_tau_c")

# Plot
p <- ggplot(plot_long_df, aes(x = share_tau_c, fill = rho)) +
  geom_density(alpha = 0.4, color = "black") +
  theme_minimal() +
  labs(
    title = "(d)",
    x = NULL,
    y = "Density",
    fill = NULL
  ) +
  xlim(0, .6)+
  scale_fill_manual(
    values = c("rho=1" = "#1f77b4", "rho=0" = "#33a02c"),
    labels = c(expression(rho == 1), expression(rho == 0))
  ) +
  theme(
    legend.position = c(0.95, 0.95),         
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white")
  )


print(p)
ggsave(
  filename = paste0("./05_ResultsAndDiagnostics/results_section_graphs/posterior_distribution_cumulative_effect_during_lockdown_share_in_overall.pdf"), 
  plot = p, 
  width = 4.5, 
  height = 4, 
  units = "in"
)

###########################################
##  Get the credible intervals and mean  ##
###########################################

summary_stats <- rbind(
  "2.5%" = apply(share_tau_c_in_tau_df, 2, quantile, probs = 0.025),
  "Mean" = colMeans(share_tau_c_in_tau_df),
  "97.5%" = apply(share_tau_c_in_tau_df, 2, quantile, probs = 0.975)
)
print(summary_stats)
# mean inf means that one cumualtive intervention effect is estimated as exactly 0
