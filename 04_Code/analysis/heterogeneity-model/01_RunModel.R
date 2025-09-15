##  Run full Bayesian factor model 
##  (including models for outcome and intervention assignment)
################################################################################
rm(list=ls())
library(rstan)
library(rstanarm)
library(bayesplot)
library(splines)
library(loo)

# maximum number of latent factors ##
k <- 1

## load the stan model ##
m_mod <- stan_model("./04_Code/analysis/heterogeneity-model/Functions/MC_full_bayes.stan",
                    auto_write=T)

## load data ##
load('./03_Data/ZZ_Temp/panel_het_blue.RData')

###########################
## extract and prep data ##
###########################
  
  print(paste("The outcome variable is:", panel_list$Outcome_var))
  
  ## set up the list of data for stan call ##
  stan_in <- list(k = k, ## k is user-specified number of factors ##
                  m = ncol(panel_list$D),       # m corresponds to T in our notation (number of time points) ##
                  n = nrow(panel_list$D),       # n is number of counties ##
                  p = panel_list$Peers, ## p is the current exposure matrix ##
                  cp = panel_list$cum_Peers, ## cp is the cumulative peer-month exposure matrix ##
                  cps = panel_list$spline_base, ## base for spline of cumulative peer month ##
                  v = dim(panel_list$spline_base)[3], ## number of spline parameters ##
                  g = as.vector(panel_list$G), ## G is a vector with first time periods units are treated ##
                  tt = panel_list$TreatTime, ## Time in treatment
                  d = panel_list$D, ## Treatment indicator matrix
                  ld = panel_list$LDtreat, ## Treatment/ Lockdown interaction
                  y = panel_list$Y, ##  Outcome ##
                  atp = panel_list$avail_tp ##  First time period peers could be employed ##
                  )


###########################
##  model specifications ##
###########################
  
  ## Decide which parameters to return ##
  PARS = c( 'Mu', 'Mu_assgn', 'phi_a', 'phi_0', 'FS', 'L', 'Ups', 
            'theta1', 'theta2',
            'gamma0', 'gamma1', 'gamma2',
            'd0', 'c0', 'Count')
  
  ## MCMC specification ##
  nChains = 2
  nThin   = 5
  nMCMC   = 20000
  nCores  = 30
  
########################
## run the stan model ##
########################
  
  rstan_options(auto_write = TRUE)
  
  Start = Sys.time()
  fit_full_bayes <- sampling(object = m_mod,
                   data = stan_in,
                   chains = nChains,
                   iter = nMCMC,
                   cores = nCores, 
                   seed = 42, 
                   thin = nThin, 
                   pars = PARS)
  Finish = Sys.time()
  Finish-Start
  
##################################################
##  First glance at treatment effect parameters ##
##################################################

  summary_stats <- summary(fit_full_bayes)
  summary_stats$summary[grep("theta", rownames(summary_stats$summary)), ]
  summary_stats$summary[grep("phi",
                             rownames(summary_stats$summary)), ]
  summary_stats$summary[grep("gamma",
                             rownames(summary_stats$summary)), ]
  
#################
##  Save Model ##
#################

  saveRDS(fit_full_bayes, file = "./03_Data/ZZ_Temp/fit_full_bayes.rds")
