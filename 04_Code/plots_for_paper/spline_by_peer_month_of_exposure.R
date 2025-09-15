##  Plot showing how psi changes with cumulative exposure ##
############################################################
rm(list =ls())

library(ggplot2)
library(rstan)
library(splines)
library(dplyr)

fit_full_bayes <- readRDS(file = "./03_Data/ZZ_Temp/fit_full_bayes.rds")
load('./03_Data/ZZ_Temp/HEP_dta_short.RData')

####################
##  Prepare data  ##
####################

##  Extract parameters  ##
theta1 <- rstan::extract(fit_full_bayes, par='theta1')[[1]]

##  Generate peer-months spline ##
input_vec <- HEP.dta.short$cum_peers 
n_knots <- 3                         
knots <- quantile(input_vec[HEP.dta.short$Dit==1], 
                  probs = seq(0, 1, length.out = n_knots + 2))[1:n_knots+1]
spline_base_long <- bs(input_vec, degree = 3, knots = knots,
                       intercept = F)

colnames(spline_base_long) <- paste0("spline_base_", seq_len(ncol(spline_base_long)))
spline_basis <- cbind(HEP.dta.short, spline_base_long) |>
  ungroup() |>
  select(cum_peers, spline_base_1, spline_base_2,
         spline_base_3, spline_base_4, spline_base_5,
         spline_base_6) |>
  distinct() |>
  arrange(cum_peers)

########################
##  Define functions  ##
########################

generate_spline <- function(spline_basis, theta1) {

  spline_draws <- exp(as.matrix(spline_basis) %*% t(theta1))
  
  # Return mean, 97.5 percentile, and 0.25 percentile
  return(list(
    mean = rowMeans(spline_draws),
    p97_5 = apply(spline_draws, 1, quantile, probs = 0.975),
    p0_25 = apply(spline_draws, 1, quantile, probs = 0.025)
  ))
}

################################
##  Create data set for plot  ##
################################

##  Calculate the psis  ##
spline_summary <- generate_spline(spline_basis[,-1], theta1)

plot_df <- data.frame(
  cum_peers = spline_basis[, 1],
  mean_psi   = unlist(spline_summary[["mean"]]),
  psi_p97_5  = unlist(spline_summary[["p97_5"]]),
  psi_p2_5   = unlist(spline_summary[["p0_25"]])
) |>
  filter(cum_peers <= 42)


############
##  Plot  ##
############

p <- ggplot(plot_df, aes(x = cum_peers, y = mean_psi)) +
  geom_line(color = "#1f78b4", size = 1) +
  geom_ribbon(aes(ymin = psi_p2_5, ymax = psi_p97_5), fill = "#1f78b4", alpha = 0.2) +
  labs(
    x = "Cumulative number of peer-months",
    y = expression(psi(bar(bold(a))[t])),
    title = "(b)"
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#e31a1c") +
  theme_minimal()

## Save ##
ggsave(
  filename = "./05_ResultsAndDiagnostics/results_section_graphs/spline_by_peer_month_of_exposure.pdf",
  plot = p, 
  width = 4.5, 
  height = 4, 
  units = "in"
)

print(p)
