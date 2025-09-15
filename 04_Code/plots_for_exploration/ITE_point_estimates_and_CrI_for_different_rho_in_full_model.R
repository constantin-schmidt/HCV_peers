##  Comparing the ITE point estimates and CrI for different values of rho
##  in full model
###########################################################################
rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(rlang)
library(grid)
library(gtable)

################
##  Load data ##
################

load("./03_Data/ZZ_Temp/panel_het_blue.RData")
load("./03_Data/ZZ_Temp/HEP_dta_short.RData")

################################################################
##  Generate a treatment effect array for full Bayesian model ##
################################################################
load('./03_Data/ZZ_Temp/posterior_predictions_copula_list_full_bayes.RData')
Y_0_0 <- posterior_predictions_copula_list[["rho=0"]]
Y_0_1 <- posterior_predictions_copula_list[["rho=1"]]

nIter <- dim(Y_0_0)[1]
nUnits <- dim(Y_0_0)[2]
nTimes <- dim(Y_0_0)[3]
Y_array <- array(rep(panel_list$Y, nIter), dim = c(nUnits, nTimes, nIter)) |>
  aperm(c(3, 1, 2))

treat_array_0_bayes_0 <- Y_array - Y_0_0
treat_array_0_bayes_1 <- Y_array - Y_0_1

############################################
##  Functions and generate summary stats  ##
############################################

##  Function to calculate the mean, 2.5% and 97.5% quantiles
generate_summary_stats <- function(arr) {
  # Mean along the first dimension (nIter)
  mean_arr <- apply(arr, c(2, 3), mean)
  
  # 2.5% quantile along the first dimension (nIter)
  q2_5_arr <- apply(arr, c(2, 3), function(x) quantile(x, 0.025))
  
  # 97.5% quantile along the first dimension (nIter)
  q97_5_arr <- apply(arr, c(2, 3), function(x) quantile(x, 0.975))
  
  # Combine results into an array
  result_plot_df <- array(NA, dim = c(3, dim(arr)[2], dim(arr)[3]))
  
  result_plot_df[1, , ] <- mean_arr
  result_plot_df[2, , ] <- q2_5_arr
  result_plot_df[3, , ] <- q97_5_arr
  
  # Return the result array
  return(result_plot_df)
}

##  Generate the summary stats
summary_0 <- generate_summary_stats(treat_array_0_bayes_0)
summary_1 <- generate_summary_stats(treat_array_0_bayes_1)

########################
##  ITE scatter plots ##
########################

##  Get means
mean_vector_0 <- as.vector(summary_0[1,,][panel_list$D==1])
mean_vector_1 <- as.vector(summary_1[1,,][panel_list$D==1])

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

######################
##  CrI scatterplot ##
######################

CrI_vector_0 <- as.vector(summary_0[3,,][panel_list$D==1]) -
  as.vector(summary_0[2,,][panel_list$D==1])
CrI_vector_1 <- as.vector(summary_1[3,,][panel_list$D==1]) -
  as.vector(summary_1[2,,][panel_list$D==1])

plot_df <- data.frame(
  CrI_0 = CrI_vector_0,
  CrI_1  = CrI_vector_1
)
range_limit <- range(c(plot_df$CrI_0, plot_df$CrI_1), na.rm = TRUE) # Ensures x- and y-axis have same range
p <- ggplot(plot_df, aes(x = CrI_0, y = CrI_1)) +
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
