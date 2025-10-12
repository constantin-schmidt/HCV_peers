##  Run all scripts in correct order  ##
########################################
rm(list = ls())

####################
##  Simulate data ##
####################
source('./04_Code/simulate_data.R')

####################
##  Prepare data  ##
####################
# For rstan, the data has to be in matrix format
source('./04_Code/analysis/prepare_data_panel.R')

##########################
##  Full Bayesian model ##
##########################
source('./04_Code/analysis/heterogeneity-model/01_RunModel.R')
source('./04_Code/analysis/heterogeneity-model/02_CounterfactualPredictions_Copula.R')
source('./04_Code/analysis/heterogeneity-model/03_PosteriorPredictions.R')
# Code for model checks is provided in the folder 'ModelChecks_NumberOfFactors'.
# It is very computationally expensive, so not run here.

####################
##  Outcome model ##
####################
source('./04_Code/analysis/only-outcome-models/01_RunModel_Cut.R')
source('./04_Code/analysis/only-outcome-models/02_CounterfactualPredictions_Copula.R')

######################################
##  Pre-intervention outcome model  ##
######################################
source('./04_Code/analysis/imputation-estimator/02_RunModel.R')
source('./04_Code/analysis/imputation-estimator/03_PosteriorPredictions.R')

######################
##  Produce graphs  ##
######################
# Fig. 1
source('./04_Code/plots_for_paper/number_of_peers.R')
source('./04_Code/plots_for_paper/number_of_patients_identified.R')

# Fig. 3
source('./04_Code/plots_for_paper/sample_cumulative_intervention_effect.R')

# Fig. 4
source('./04_Code/plots_for_paper/mean_ITEs_for_full_Bayesian_model_with_rho_0.R')
source('./04_Code/plots_for_paper/spline_by_peer_month_of_exposure.R')
source('./04_Code/plots_for_paper/posterior_density_theta2.R')
source('./04_Code/plots_for_paper/posterior_distribution_cumulative_effect_during_lockdown_share_in_overall.R')

# Fig. 5
source('./04_Code/plots_for_paper/ITE_by_time_for_all_ODNs_rho_1.Rmd.R')
source('./04_Code/plots_for_paper/ITE_by_time_for_all_ODNs_rho_0.Rmd')
