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
load("./03_Data/ZZ_Temp/posterior_predictions_copula_list_outcome.RData")
list_outcome_only <- posterior_predictions_copula_list
load("./03_Data/ZZ_Temp/posterior_predictions_imp.RData")

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
  
  summary_list <- lapply(names(cum_effects_df), function(method_name) {
    vec <- cum_effects_df[[method_name]]
    data.frame(
      p2.5 = quantile(vec, 0.025),
      mean = mean(vec),
      median = median(vec),
      p97.5 = quantile(vec, 0.975),
      perc_greater_than_0 = mean(vec > 0) * 100,
      method = method_name,
      stringsAsFactors = FALSE
    )
  })
  
  # Combine all summary data frames into one
  cum_effects <- do.call(rbind, summary_list)
  
  # Reorder columns to place 'method' first
  cum_effects <- cum_effects[, c("method", setdiff(names(cum_effects), "method"))]
  
  return(cum_effects)
}

## Imputation model ##
compute_cum_effects_imp <- function(Y0_imp, panel_data_Y, panel_data_D) {
  
  nIter <- dim(Y0_imp)[1]
  nUnits <- dim(Y0_imp)[2]
  nTimes <- dim(Y0_imp)[3]
  
  Y_array <- array(rep(panel_data_Y, nIter), dim = c(nUnits, nTimes, nIter)) |>
    aperm(c(3, 1, 2))
  D_array <- array(rep(panel_data_D, nIter), dim = c(nUnits, nTimes, nIter)) |>
    aperm(c(3, 1, 2))
  
  treat_array_imp <- (Y_array - Y0_imp) * D_array

  # summary_df
    vec <- apply(treat_array_imp, 1, sum)
    cum_effect_imp <- data.frame(
      p2.5 = quantile(vec, 0.025),
      mean = mean(vec),
      median = median(vec),
      p97.5 = quantile(vec, 0.975),
      perc_greater_than_0 = mean(vec > 0) * 100,
      method = "rho=0",
      stringsAsFactors = FALSE
    )
    cum_effect_imp <- cum_effect_imp[, c("method", setdiff(names(cum_effect_imp), "method"))]
  
  return(cum_effect_imp)
}

########################################
##  Average Effects for Both Lists ##
########################################

cum_effects_assign <- compute_cum_effects(list_assign_model, panel_list$Y, panel_list$D)
cum_effects_outcome_only <- compute_cum_effects(list_outcome_only, panel_list$Y, panel_list$D)

# Average effects for Imp model
cum_effects_imp <- compute_cum_effects_imp(Y0_imp, panel_list$Y, panel_list$D)

##############################
## Combine data for plotting ##
##############################

# Add model identifiers with new names
cum_effects_assign <- cum_effects_assign |> mutate(model = "Full model")
cum_effects_outcome_only <- cum_effects_outcome_only |> mutate(model = "Outcome model")
cum_effects_imp <- cum_effects_imp |> mutate(model = "Pre-intervention outcome model")

# Combine all into one dataframe
cum_effects_combined <- bind_rows(cum_effects_assign, cum_effects_outcome_only, cum_effects_imp)

##################
## Create expression-based labels ##
##################

# Define expression-based labels for proper rendering
expression_labels <- c(
  "rho=1"             = expression(rho == 1),
  "rho=0.75"          = expression(rho == 0.75),
  "rho=0.5"           = expression(rho == 0.5),
  "rho=0"             = expression(rho == 0),
  "rho=-1"            = expression(rho == -1),
  "rho=unif[0.75,1]"  = expression(rho %~% Unif(0.75, 1)),
  "rho=unif[0.5,1]"   = expression(rho %~% Unif(0.5, 1)),
  "rho=unif[0,1]"     = expression(rho %~% Unif(0, 1)),
  "rho=unif[-1,1]"    = expression(rho %~% Unif(-1, 1))
)

# Add labels as factor
cum_effects_combined$method_label <- factor(
  cum_effects_combined$method,
  levels = rev(c(
    "rho=1", "rho=0.75", "rho=0.5", "rho=0", "rho=-1",
    "rho=unif[0.75,1]", "rho=unif[0.5,1]", "rho=unif[0,1]", "rho=unif[-1,1]"
  ))
)

##########################
## Caterpillar Plot ##
##########################

# Ensure the 'model' column has the correct factor levels
cum_effects_combined$model <- factor(cum_effects_combined$model, 
                                    levels = c("Full model", "Outcome model", "Pre-intervention outcome model"))

# Integrate label text
cum_effects_combined <- cum_effects_combined |>
  mutate(label_text = paste0(
    round(mean, 1), " [", round(p2.5, 1), ", ", round(p97.5, 1), "]; ",
    round(perc_greater_than_0, 1), "%"
  )) |>
  mutate(label_text = if_else(
    perc_greater_than_0 >= 99.950,
    str_replace(label_text, "100%", "> 99.9%"),
    label_text
  ))


# Calculate space for right-aligned labels
label_position <- max(cum_effects_combined$p97.5) + 100

# Define the plot
p <- ggplot(cum_effects_combined, aes(y = method_label, x = mean, color = model)) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  geom_errorbar(
    aes(xmin = p2.5, xmax = p97.5),
    position = position_dodge(width = 0.7),
    width = 0.2
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_text(
    aes(
      x = label_position,
      label = label_text
    ),
    position = position_dodge(width = 0.7),
    size = 4,
    hjust = 0,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Full model" = "#1f78b4",
    "Outcome model" = "#33a02c",
    "Pre-intervention outcome model" = "#e31a1c"
  )) +
  scale_x_continuous(
    breaks = c(0, 1000, 2000, 3000),
    labels = c("0", "1000", "2000", "3000"),
    expand = expansion(mult = c(0.1, 0.5))
  ) +
  scale_y_discrete(labels = expression_labels, expand = expansion(add = c(0, 1))) +
  labs(
    title = NULL,
    x = NULL,
    y = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank(), 
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5), 
    axis.ticks.x = element_line() 
  )

# Add annotation above the plot
p <- p + 
  annotate(
    "text", 
    x = label_position,
    y = length(unique(cum_effects_combined$method_label)) + .75, 
    label = "Mean [95%-CrI], PPos",
    hjust = 0, 
    size = 4
  )

# Show the plot
print(p)

# Save to file
ggsave(
  filename = "./05_ResultsAndDiagnostics/results_section_graphs/sample_cumulative_TE_with_differing_rho_all_models.pdf", 
  plot = p, 
  width = 10, 
  height = 7, 
  units = "in"
)
