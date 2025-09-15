####################################################
##   Posterior and Prior Densities for phi        ##
####################################################

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

# Prior predictive data
load('./03_Data/ZZ_Temp/HEP_dta_short.RData')
load('./03_Data/ZZ_Temp/panel_het_blue.RData')

####################
##  Posterior     ##
####################

# Convert posterior to array
fit_full_bayes_array <- as.array(fit_full_bayes)
summary_stats <- summary(fit_full_bayes)

# Extract phi parameters
summary_phi <- summary_stats$summary[grep("phi", rownames(summary_stats$summary)), ]
paras <- rownames(summary_phi)

####################
##  Posterior plot ##
####################

# Create the posterior density plot
p <- mcmc_areas(
  fit_full_bayes_array, 
  pars = paras,
  prob = 0.8,
  prob_outer = 0.99,
  point_est = "mean"
) + 
  scale_y_discrete(
    labels = c(
      "phi_a" = expression(varphi^1),
      "phi_0" = expression(varphi^0)
    )
  )

# Show plot
print(p)
