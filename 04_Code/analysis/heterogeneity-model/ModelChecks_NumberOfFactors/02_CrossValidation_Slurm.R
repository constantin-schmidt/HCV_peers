##  This should be done on a high performance computer  ##

rm(list=ls())
Start = Sys.time()
library(rstan)
library(rstanarm)
library(bayesplot)
library(splines)
library(foreach)
library(doParallel)
library(doRNG)

##################
##  Get task id ##
##################

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", unset=1))

################
##  Set seed  ##
################

starting_seed <- 42
set.seed(starting_seed + task_id)

####################
##  Load the data ##
####################

## load the stan model ##
load(file = 'compiled_MC_full_bayes_CV.rds')
#load('./04_Code/Analysis_BMCM_BlueTeq/heterogeneity-model/Functions/compiled_MC_full_bayes_CV.rds')

## load the full list of panel data ##
#load('./03_Data/ZZ_Temp/panel_het_blue.RData')
load('../Data/panel_het_blue.RData')
#load('./03_Data/ZZ_Temp/HEP_dta_short.RData')
load('../Data/HEP_dta_short.RData')
#load('./03_Data/ZZ_Temp/CV_het_data_blue.RData')
load('../Data/CV_het_data_blue.RData')

###########################
## extract and prep data ##
###########################

##  Prepare number of intervention periods ##
n_treat <- nrow(panel_list$Y) * ncol(panel_list$Y) - 
  sum(pmin(panel_list$G, ncol(panel_list$Y)+1)) + nrow(panel_list$Y)

##  Spline  ##  
##  create the spline base  ##
input_vec <- HEP.dta.short$cum_peers ## Decide the variable to base spline on
n_knots <- 3                         ## Decide number of knots
knots <- quantile(input_vec[HEP.dta.short$Dit==1], 
                  probs = seq(0, 1, length.out = n_knots + 2))[1:n_knots+1]
spline_base_long <- bs(input_vec, degree = 3, knots = knots)

##  transform into matrix ##
n_spline_parameters <- ncol(spline_base_long)
spline_base_array <- array(NA, dim = c(dim(panel_list$Y)[1],
                                       dim(panel_list$Y)[2],
                                       n_spline_parameters))
for (i in 1:n_spline_parameters) {
  spline_base_array[,,i] <- matrix(spline_base_long[,i], 
                                   nrow = dim(panel_list$Y)[1],
                                   ncol = dim(panel_list$Y)[2],
                                   byrow = T)
}

##  Intervention start times ##
  G = panel_list$G

################
##  Functions ##
################

# Function to compute squared error only for exposed ODNs
compute_squared_error <- function(A, B, O) {
  nIter <- dim(B)[1]
  nUnits <- dim(A)[1]
  squared_error <- matrix(NA, nrow = nUnits, ncol = nIter)
  
  # Iterate over ODNs
  for (i in 1:nUnits) {
    # Compute squared errors for exposed ODNs
    if (O[i] < max(O)) {
      for (b in 1:nIter) {
        squared_error[i, b] <- (A[i, O[i]] - B[b, i, O[i]])^2
      }
    } else if (O[i] == max(O)) {
      # Leave NA for unexposed ODNs
      squared_error[i, ] <- NA
    }
  }
  
  return(squared_error)
}

# Function to compute interval score only for exposed ODNs
  ##  set confidence level  ##
  alpha <- 0.05
  ## Function ##
  compute_interval_score <- function(alpha, A, B, O) {
    # Initialize interval score
    interval_score <- numeric(length(O))
    
    # Precompute 2.5th and 97.5th percentiles
    quantile_97_5 <- apply(B, c(2, 3), function(x) quantile(x, 0.975))
    quantile_2_5 <- apply(B, c(2, 3), function(x) quantile(x, 0.025))
    
    # Compute the interval scores for exposed ODNs
    for (i in seq_along(O)) {
      if (O[i] < max(O)) {
        
        # Retrieve appropriate indices
        o_index <- O[i]
        a_val <- A[i, o_index]
        q_2_5 <- quantile_2_5[i, o_index]
        q_97_5 <- quantile_97_5[i, o_index]
        
        # Compute the interval score
        interval_score[i] <- (q_97_5 - q_2_5) +
          (2 / alpha) * (q_2_5 - a_val) * (a_val < q_2_5) +
          (2 / alpha) * (a_val - q_97_5) * (a_val > q_97_5)
      }
      
    # Assign NA for unexposed ODNs
      if (O[i] == max(O)) {
        interval_score[i] <- NA
      }
    }
    
    return(interval_score)
  }
  
###########################
##  model specifications ##
###########################

## Decide which parameters to return ##
PARS = c('Mu', 'phi_0', 'phi_a')

## For CV ##
K <- 5 # Maximum number of factors
cl <- makeCluster(K+1, outfile="/dev/stdout")
registerDoParallel(cl) #cores = K+1, outfile="/dev/stdout")
m <- task_id # Which CV data set is used

## MCMC specification ##
nChains = 2
nThin   = 5
nMCMC   = 2000

##  rstan options
rstan_options(auto_write = TRUE)

########################
## run the stan model ##
########################
  
  ##  Prepare training data
  CV_dat <- CV_het_data_blue[[m]]$modified_mat
  CV_dat[is.na(CV_dat)] <- 99999
  
  ##  Vector of test observations for each unit
  O <- CV_het_data_blue[[m]]$g_col
  O[is.na(O)] <- ncol(panel_list$Y) + 1
  O <- as.vector(O)
  
## Loop through different number of factors
score_mat <- foreach(k = 0:K, .packages=c("rstan")) %dorng% {
    
    ## set up the list of data for stan call ##
  stan_in <- list(k = k, ## k is user-specified number of factors ##
                  m = ncol(panel_list$D),       # m corresponds to T in our notation (number of time points) ##
                  n = nrow(panel_list$D),       # n is number of counties ##
                  n_treat = n_treat, ##  number of treated observations ##
                  p = panel_list$Peers, ## p is the current exposure matrix ##
                  cp = panel_list$cum_Peers, ## cp is the cumulative peer-month exposure matrix ##
                  cps = spline_base_array, ## base for spline of cumulatve peer month ##
                  v = n_spline_parameters, ## number of spline parameters ##
                  g = as.vector(panel_list$G), ## G is a vector with first time periods units are treated ##
                  tt = panel_list$TreatTime, ## Time in treatment
                  d = panel_list$D, ## Treatment indicator matrix
                  ld = panel_list$LDtreat, ## Treatment/ Lockdown interaction
                  y = panel_list$Y, ##  Outcome ##
                  atp = panel_list$avail_tp, ##  First time period peers could be employed ##
                    
                  ## Cross validation inputs
                  y = CV_dat, ## Outcome matrix with missing values for CV ##
                  o = O ## vector indicating the observations left out for CV ##
    )
    
      ##  Estimate model  ##
      fit <- sampling(object = MC_full_bayes_mod,
                      data = stan_in,
                      chains = nChains,
                      iter = nMCMC,
                      cores = nChains, 
                      thin = nThin, 
                      pars = PARS)
      
      ##  Extract parameters  ##
      mu = rstan::extract(fit, par = 'Mu')[[1]]
      phi_0 = rstan::extract(fit, par='phi_0')[[1]]
      phi_a = rstan::extract(fit, par='phi_a')[[1]]
      
      ##  Posterior predictions ##
      PosteriorPredictions = NA*mu
      nIter = dim(mu)[1]
      nUnits = dim(mu)[2]
      nTimes = dim(mu)[3]
      for (b in 1:nIter) {
        for (i in 1:nUnits) {
          for (t in 1:nTimes) {
            if (t < G[i]) {
              PosteriorPredictions[b,i,t] <- rnbinom(1, 
                                               mu=mu[b,i,t], size=phi_0[b])
            } else {
              PosteriorPredictions[b,i,t] <- rnbinom(1, 
                                               mu=mu[b,i,t], size=phi_a[b])
            }
          }
        }
      }
      
      ##  Calculate mean squared error ##
      squared_error <- compute_squared_error(panel_list$Y, PosteriorPredictions, O)
      mse <- mean(squared_error, na.rm = T)
      
      ##  Calculate interval score ##
      interval_score <- compute_interval_score(alpha, panel_list$Y, PosteriorPredictions, O)
      mean_interval_score <- mean(interval_score, na.rm = T)
      
      ## Return ##
      return(list(mse = mse, 
                  mean_interval_score = mean_interval_score, 
                  squared_error = squared_error, 
                  interval_score = interval_score))
}
  
################
## Save data  ##
################
  
  # Combine results into separate outputs
  mses <- sapply(score_mat, function(x) x$mse)
  interval_scores <- sapply(score_mat, function(x) x$mean_interval_score)
  
  # Bind mse and interval scores into a matrix
  score_matrix <- cbind(MSE = mses, Interval_Score = interval_scores)
  
  # Save the score matrix
  save(score_matrix, file = paste0("../Output/scores_", task_id, ".rds"))
  
  # Save detailed squared errors
  squared_errors <- lapply(score_mat, function(x) x$squared_error)
  save(squared_errors, file = paste0("../Output/squared_errors_", task_id, ".rds"))
  
  # Save detailed interval scores
  interval_scores_full <- lapply(score_mat, function(x) x$interval_score)
  save(interval_scores_full, file = paste0("../Output/interval_scores_", task_id, ".rds"))
  
  ## Time the script
  Finish = Sys.time()
  print(Finish-Start)
  