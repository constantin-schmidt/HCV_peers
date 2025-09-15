##  Posterior predictions of the counterfactuals using a copula
##############################################################################

fit_het_blue <- readRDS(file = "./03_Data/ZZ_Temp/fit_het_blue.rds")
load('./03_Data/ZZ_Temp/panel_het_blue.RData')

############################
## Extract the parameters ##
############################

phi_a <- rstan::extract(fit_het_blue, par='phi_a')[[1]]
phi_0 <- rstan::extract(fit_het_blue, par='phi_0')[[1]]
mu <- rstan::extract(fit_het_blue, par='Mu')[[1]]
count <- rstan::extract(fit_het_blue, par='Count')[[1]]

Y <- panel_list$Y
G <- panel_list$G

nIter <- dim(mu)[1]
nUnits <- dim(mu)[2]
nTimes <- dim(mu)[3]

##############################################
##  Set prior for the correlation parameter ##
##############################################

rho_vals <- c(1, .75, .5, 0, -1)
rho_mins <- c(.75, .5, 0, -1)
rho_maxs <- c(1, 1, 1, 1)

# one correlation value
rho_static <- as.data.frame(matrix(NA, nrow = length(rho_vals), ncol = nIter))
for (val in seq_along(rho_vals)) { 
  rho_static[val, ] <- rep(rho_vals[val], nIter)
}

# correlation values drawn from prior distribution
rho_prior <- as.data.frame(matrix(NA, nrow = length(rho_maxs), ncol = nIter))
for (val in 1:length(rho_mins)) { 
  rho_prior[val, ] <- runif(nIter,
                            min = rho_mins[val],
                            max = rho_maxs[val])
}

rho_df <- rbind(rho_static, rho_prior)
rho_vector <- c(as.character(rho_vals), paste0("unif[", rho_mins, ",", rho_maxs,"]"))

####################################################
##  Get the percentiles from outcome distribution ##
####################################################

## Get P(x <= X)  ##
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

########################
##  Draw a percentile ##
########################

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

##############################################################
##  Introduce the correlation and generate counterfactuals  ##
##############################################################
posterior_predictions_copula_list <- list()

# Loop over each correlation value
for (row_idx in 1:nrow(rho_df)) {
  
  # Get the specific rho value for this iteration
  rho <- as.numeric(rho_df[row_idx, ])
  
  ################################
  ##  Introduce the correlation ##
  ################################
  
  percentile2 <- NA*mu
  norm_var1 <- NA*mu
  mu2_given_norm_var1 <- NA*mu
  sd2_given_norm_var1 <- NA*nIter
  norm_var2 <- NA*mu
  
  for (b in 1:nIter) {
    for (i in 1:nUnits) {
      for (t in 1:nTimes) {
        ## Get first normally distributed variable
        norm_var1[b,i,t] <- qnorm(percentile[b,i,t], mean = 0, sd = 1)
        
        ## Compute conditional distribution of norm_var2 | norm_var1
        mu2_given_norm_var1[b,i,t] <- rho[b] * norm_var1[b,i,t]
        sd2_given_norm_var1[b] <- sqrt(1 - rho[b]^2)
        
        ## Sample norm_var2 from conditional distribution
        norm_var2[b,i,t] <- rnorm(1, 
                                  mean = mu2_given_norm_var1[b,i,t], 
                                  sd = sd2_given_norm_var1[b])
        
        ## Compute percentile of norm_var2 in the conditional distribution
        percentile2[b,i,t] <- pnorm(norm_var2[b,i,t], mean = 0, sd = 1)
      }
    }
  }
  
  ##############################
  ##  Generate counterfactual ##
  ##############################
  
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

##########
## Save ##
##########

save(posterior_predictions_copula_list, 
     file = "./03_Data/ZZ_Temp/posterior_predictions_copula_list_outcome.RData")
