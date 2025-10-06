##  Posterior predictions of the counterfactuals using a copula ##
##################################################################
# Example with overdispersed count outcome Y
# As count outcomes are discrete, any value y corresponds to a range of qunatiles
# of the CDF.
# Thus, we use the following steps to implement our copula approach:
# 1.) Set values for rho, the parameter governing the correlation between
#     exposed and unexposed potential outcomes. rho might be a point prior or
#     uniformly distributed (in which case minimum and maximum are set).
# 2.) Randomly draw a quantile from th possible quantiles of the CDF for which
#     Y is equal the observed value.
# 3.) Use a copula to introduce correlation between potential outcomes and
#     generate a counterfactual
#     3.1)  Using a Gaussian copula, draw from a value from a bivariate standard
#           normal distribution conditional and find corresponding quantile.
#     3.2)  Find the value in the posterior predictive of the counterfactual
#           outcome that corresponds to this qunatile

########################################
##  Run all scripts in correct order  ##
########################################
rm(list = ls())

################################
##  Create data and posterior ##
################################
# Up to 30 minutes
source('./04_Code/simulate_data.R')
source('./04_Code/analysis/pepare_data_panel.R')
source('./04_Code/analysis/heterogeneity-model/01_RunModel.R')
source('./04_Code/analysis/heterogeneity-model/02_CounterfactualPredictions_Copula.R')
source('./04_Code/analysis/heterogeneity-model/03_PosteriorPredictions.R')

##  Load the stan fit ##
fit_full_bayes <- readRDS(file = "./03_Data/ZZ_Temp/fit_full_bayes.rds")

##  Load the original data as matrix list ##
load('./03_Data/ZZ_Temp/panel_het_blue.RData')
# The data needs to be in matrix format with units as rows and time as columns
# See 'prepare_data_panel.R' for code to do that

##################
## Prepare data ##
##################

##  Extract the relevant parameters ##
# Dispersion parameter for outcomes under intervention
phi_a <- rstan::extract(fit_full_bayes, par='phi_a')[[1]]
# Dispersion parameter for outcomes under no intervention
phi_0 <- rstan::extract(fit_full_bayes, par='phi_0')[[1]]
# Mean parameter for observed outcomes
mu <- rstan::extract(fit_full_bayes, par='Mu')[[1]]
# Mean parameter for outcomes under no intervention
count <- rstan::extract(fit_full_bayes, par='Count')[[1]]

##  Load data ##
# Observed outcomes
Y <- panel_list$Y
# Intervention start times
G <- panel_list$G

##  Set dimensions  ##
# Number of iterations
nIter <- dim(mu)[1]
# Number of units
nUnits <- dim(mu)[2]
# Number of time periods
nTimes <- dim(mu)[3]

##################################################
##  1.) Set prior for the correlation parameter ##
##################################################

##  Set values for rho  ##
# Point priors
rho_vals <- c(1, .75, .5, 0, -1)
# Maximum and minimum for rho ~ unif(min, max)
rho_mins <- c(.75, .5, 0, -1)
rho_maxs <- c(1, 1, 1, 1)

##  Create a data frame with rho values ##
# Transform the point priors to data frame
rho_static <- as.data.frame(matrix(NA, nrow = length(rho_vals), ncol = nIter))
for (val in seq_along(rho_vals)) { 
  rho_static[val, ] <- rep(rho_vals[val], nIter)
}

# Draw a value from rho ~ unif(min, max)
rho_prior <- as.data.frame(matrix(NA, nrow = length(rho_maxs), ncol = nIter))
for (val in 1:length(rho_mins)) { 
  rho_prior[val, ] <- runif(nIter,
                            min = rho_mins[val],
                            max = rho_maxs[val])
}

# Bind all the rho values into one data frame
rho_df <- rbind(rho_static, rho_prior)
rho_vector <- c(as.character(rho_vals), paste0("unif[", rho_mins, ",", rho_maxs,"]"))

######################################################
##  2.) Get the quantiles from outcome distribution ##
######################################################

## Get P(y <= Y)  ##
# Probability of observing an outcome less or equal observed outcome
# given the posterior of the parameters
percentile_max <- NA*mu
for (b in 1:nIter) {
  for (i in 1:nUnits) {
    for (t in 1:nTimes) {
      if (t < G[i]) {
        percentile_max[b,i,t] <- pnbinom(Y[i,t], 
                                         mu=mu[b,i,t],
                                         size=phi_0[b])
      } else {
        percentile_max[b,i,t] <- pnbinom(Y[i,t], 
                                         mu=mu[b,i,t],
                                         size=phi_a[b])
      }
    }
  }
}

## Get P(x <= X-1)  ##
# Probability of observing an outcome less or equal observed outcome minus 1
# given the posterior of the parameters
percentile_min <- NA*mu
for (b in 1:nIter) {
  for (i in 1:nUnits) {
    for (t in 1:nTimes) {
      if (t < G[i]) {
        percentile_min[b,i,t] <- pnbinom(Y[i,t]-1,
                                         mu=mu[b,i,t],
                                         size=phi_0[b])
      } else {
        percentile_min[b,i,t] <- pnbinom(Y[i,t]-1,
                                         mu=mu[b,i,t],
                                         size=phi_a[b])
      }
    }
  }
}

##  Draw a quantile ##
# Draw a quantile from a unif(percentile_min, percentile_max) as the 
# quantile for the observed outcome
percentile <- NA*mu
for (b in 1:nIter) {
  for (i in 1:nUnits) {
    for (t in 1:nTimes) {
      percentile[b,i,t] <- runif(1, 
                                 min = percentile_min[b,i,t], 
                                 max = percentile_max[b,i,t])
    }
  }
}

##################################################################
##  3.) Introduce the correlation and generate counterfactuals  ##
##################################################################
posterior_predictions_copula_list <- list()

# Loop over each correlation value
for (row_idx in 1:nrow(rho_df)) {
  
  print(paste("Generating posteriors for: rho=", rho_vector[row_idx]))
  
  # Get the specific rho value for this iteration
  rho <- as.numeric(rho_df[row_idx, ])
  
  ######################################
  ##  3.1)  Introduce the correlation ##
  ######################################
  
  # Predefine
  percentile2 <- NA*mu
  norm_var1 <- NA*mu
  mu2_given_norm_var1 <- NA*mu
  sd2_given_norm_var1 <- NA*nIter
  norm_var2 <- NA*mu
  
  # Loop over iterations, units, and time periods
  for (b in 1:nIter) {
    for (i in 1:nUnits) {
      for (t in 1:nTimes) {
        
        ## Get first normally distributed variable for copula
        # i.e. get the value that corresponds to the drawn quantile on a
        # standard normal distribution
        norm_var1[b,i,t] <- qnorm(percentile[b,i,t], mean = 0, sd = 1)
        
        ## Compute conditional distribution of norm_var2 | norm_var1
        # This is where rho is used as governing the correlation
        mu2_given_norm_var1[b,i,t] <- rho[b] * norm_var1[b,i,t]
        sd2_given_norm_var1[b] <- sqrt(1 - rho[b]^2)
        
        ## Draw norm_var2 from a bivariate normal distribution given norm_var1
        norm_var2[b,i,t] <- rnorm(1, 
                                  mean = mu2_given_norm_var1[b,i,t], 
                                  sd = sd2_given_norm_var1[b])
        
        ## Compute percentile of norm_var2 in the conditional distribution
        percentile2[b,i,t] <- pnorm(norm_var2[b,i,t], mean = 0, sd = 1)
      }
    }
  }
  
  ###################################
  ##  3.2) Generate counterfactual ##
  ###################################
  # Generate the counterfactual using the posteriors of the parameters and 
  # the percentile of the drawn norm_var2
  Y_0 <- NA*mu
  for (b in 1:nIter) {
    for (i in 1:nUnits) {
      for (t in 1:nTimes) {
        Y_0[b,i,t] <- qnbinom(p = percentile2[b,i,t],
                              mu = count[b,i,t],
                              size = phi_0[b])
      }
    }
  }
  
  ## Save the results for this rho value to the list
  posterior_predictions_copula_list[[row_idx]] <- Y_0
  names(posterior_predictions_copula_list)[row_idx] <- paste0("rho=",
                                                              rho_vector[row_idx])
  
}
