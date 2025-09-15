##  Plots comparing point estimates and CrI length for all 3 models ##
######################################################################
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
Y_0 <- posterior_predictions_copula_list[["rho=0"]]

nIter <- dim(Y_0)[1]
nUnits <- dim(Y_0)[2]
nTimes <- dim(Y_0)[3]
Y_array <- array(rep(panel_list$Y, nIter), dim = c(nUnits, nTimes, nIter)) |>
  aperm(c(3, 1, 2))

treat_array_full_bayes_indep <- Y_array - Y_0

################################################################
##  Generate a treatment effect array for independent model ##
################################################################
load("./03_Data/ZZ_Temp/posterior_predictions_copula_list_outcome.RData")
Y_0 <- posterior_predictions_copula_list[["rho=0"]]
nIter <- dim(Y_0)[1]
Y_array <- array(rep(panel_list$Y, nIter), dim = c(nUnits, nTimes, nIter)) |>
  aperm(c(3, 1, 2))

treat_array_cut_indep <- Y_array - Y_0

######################################################################
##  Generate a treatment effect array for Bayesian imputation model ##
######################################################################
load("./03_Data/ZZ_Temp/posterior_predictions_imp.RData")

treat_array_imp <- Y_array - Y0_imp

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
summary_full_bayes <- generate_summary_stats(treat_array_full_bayes_indep)
summary_outcome_indep <- generate_summary_stats(treat_array_cut_indep)
summary_imp_indep <- generate_summary_stats(treat_array_imp)

########################
##  ITE scatter plots ##
########################

##  Get means
mean_vector_full <- as.vector(summary_full_bayes[1,,][panel_list$D==1])
mean_vector_out <- as.vector(summary_outcome_indep[1,,][panel_list$D==1])
mean_vector_imp <- as.vector(summary_imp_indep[1,,][panel_list$D==1])
  # [panel_list$D==1] ensures only observations in intervention are considered

##  Plot vs. imputation model
plot_df <- data.frame(
  mean_full = mean_vector_full,
  mean_imp  = mean_vector_imp
)
range_limit <- range(c(plot_df$mean_full, plot_df$mean_imp), na.rm = TRUE)
p <- ggplot(plot_df, aes(x = mean_full, y = mean_imp)) +
  geom_point(alpha = 0.7) +
  labs(
    title = NULL,
    x = "Full model",
    y = "Pre-intervention outcome model"
  ) +
  theme_minimal() +
  xlim(range_limit) +
  ylim(range_limit)
  ##  Save  ##
  ggsave(
    filename = "./05_ResultsAndDiagnostics/results_section_graphs/compare_point_estimates_CrI_length_all_models/means_vs_imputation_model.pdf", 
    plot = p, 
    width = 4, 
    height = 4, 
    units = "in"
  )

##  Plot vs. outcome model
plot_df <- data.frame(
  mean_full = mean_vector_full,
  mean_out  = mean_vector_out
)
range_limit <- range(c(plot_df$mean_full, plot_df$mean_out), na.rm = TRUE)
p <- ggplot(plot_df, aes(x = mean_full, y = mean_out)) +
  geom_point(alpha = 0.7) +
  labs(
    title = NULL,
    x = "Full model",
    y = "Outcome model"
  ) +
  theme_minimal() +
  xlim(range_limit) +
  ylim(range_limit)
  ##  Save  ##
  ggsave(
    filename = "./05_ResultsAndDiagnostics/results_section_graphs/compare_point_estimates_CrI_length_all_models/means_vs_outcome_model.pdf", 
    plot = p, 
    width = 4, 
    height = 4, 
    units = "in"
  )

######################
##  CrI scatterplot ##
######################

CrI_vector_full <- as.vector(summary_full_bayes[3,,][panel_list$D==1]) -
  as.vector(summary_full_bayes[2,,][panel_list$D==1])
CrI_vector_imp <- as.vector(summary_imp_indep[3,,][panel_list$D==1]) -
  as.vector(summary_imp_indep[2,,][panel_list$D==1])
CrI_vector_out <- as.vector(summary_outcome_indep[3,,][panel_list$D==1]) -
  as.vector(summary_outcome_indep[2,,][panel_list$D==1])

##  Plot vs. imputation model
plot_df <- data.frame(
  CrI_full = CrI_vector_full,
  CrI_imp  = CrI_vector_imp
)
range_limit <- range(c(plot_df$CrI_full, plot_df$CrI_imp), na.rm = TRUE) # Ensures x- and y-axis have same range
p <- ggplot(plot_df, aes(x = CrI_full, y = CrI_imp)) +
  geom_point(alpha = 0.7) +
  labs(
    title = NULL,
    x = "Full model",
    y = "Pre-intervention outcome model"
  ) +
  theme_minimal() +
  xlim(range_limit) +
  ylim(range_limit)
  ##  Save  ##
  ggsave(
    filename = "./05_ResultsAndDiagnostics/results_section_graphs/compare_point_estimates_CrI_length_all_models/CrI_vs_imputation_model.pdf", 
    plot = p, 
    width = 4, 
    height = 4, 
    units = "in"
  )

##  Plot vs. outcome model
plot_df <- data.frame(
  CrI_full = CrI_vector_full,
  CrI_out  = CrI_vector_out
)
range_limit <- range(c(plot_df$CrI_full, plot_df$CrI_out), na.rm = TRUE) # Ensures x- and y-axis have same range
p <- ggplot(plot_df, aes(x = CrI_full, y = CrI_out)) +
  geom_point(alpha = 0.7) +
  labs(
    title = NULL,
    x = "Full model",
    y = "Outcome model"
  ) +
  theme_minimal() +
  xlim(range_limit) +
  ylim(range_limit)
  ##  Save  ##
  ggsave(
    filename = "./05_ResultsAndDiagnostics/results_section_graphs/compare_point_estimates_CrI_length_all_models/CrI_vs_outcome_model.pdf", 
    plot = p, 
    width = 4, 
    height = 4, 
    units = "in"
  )