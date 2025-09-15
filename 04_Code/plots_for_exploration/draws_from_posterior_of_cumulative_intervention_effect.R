################################################################################
##  Caterpillar plot for SATE with different values for rho
################################################################################

library(ggplot2)
library(dplyr)
library(tidyverse)

################
##  Load data ##
################

load('./03_Data/ZZ_Temp/panel_het_blue.RData')
load('./03_Data/ZZ_Temp/posterior_predictions_copula_list_full_bayes.RData')
list_assign_model <- posterior_predictions_copula_list

########################################################
##  Define a function to calculate average effects ##
########################################################

compute_cum_effects <- function(posterior_list, panel_data_Y, panel_data_D) {
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
  
  cum_effects_df <- data.frame(matrix(ncol = length(treat_array_list),
                                      nrow = dim(treat_array_list[[1]])[1]))
  for (i in seq_along(treat_array_list)) {
    cum_effects_df[[i]] <- apply(treat_array_list[[i]], 1, sum)
  }
  names(cum_effects_df) <- names(treat_array_list)
  
  return(cum_effects_df)
}

cum_effects_assign <- compute_cum_effects(list_assign_model, panel_list$Y, panel_list$D)

apply(cum_effects_assign, 2, function(x) quantile(x, .975)-quantile(x, 0.025))
apply(cum_effects_assign, 2, var)

####################################
##  Posterior draws scatter plots ##
####################################

##  Get means
mean_vector_0 <- as.vector(cum_effects_assign[["rho=0"]])
mean_vector_1 <- as.vector(cum_effects_assign[["rho=1"]])

##  Plot vs. imputation model
plot_df <- data.frame(
  mean_0 = mean_vector_0,
  mean_1  = mean_vector_1
)
range_limit <- range(c(plot_df$mean_0, plot_df$mean_1), na.rm = TRUE)
p <- ggplot(plot_df, aes(x = mean_0, y = mean_1)) +
  geom_point(alpha = 0.7) +
  labs(
    title = NULL,
    x = "rho = 0",
    y = "rho = 1"
  ) +
  theme_minimal() +
  xlim(range_limit) +
  ylim(range_limit)
print(p)
