##  Model using pre-intervention outcomes only
################################################################################
rm(list=ls())
library(rstan)
library(rstanarm)
library(bayesplot)

# number of latent factors ##
k <- 1

## load the stan model ##
m_mod <- stan_model("./04_Code/analysis/imputation-estimator/Functions/MC_neg_binom.stan",auto_write=T)

## load the full list of panel data ##
load('./03_Data/ZZ_Temp/panel_het_blue.RData')

###########################
## extract and prep data ##
###########################

## stan doesn't allow missing values so we put a numeric "placeholder" in the missings
Y0_obs <- ifelse(panel_list$D==1,
                 99999,
                 panel_list$Y)

## set up the list of data for stan call ##
stan_in <- list(k = k, ## k is user-specified number of factors ##
                m = ncol(Y0_obs), # m corresponds to T in our notation (number of time points) ##
                n = nrow(Y0_obs), # n is number of counties ##
                g = as.vector(panel_list$G), ## G is a vector with first time periods units are treated ##
                y = Y0_obs ## y is the observed Y(0) matrix ##
)

###########################
##  model specifications ##
###########################

## Decide which parameters to return ##
PARS = c('L', 'FS', 'phi', 'Mu', 'd0', 'c0', 'Ups')

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
fit_imp_blue <- sampling(object = m_mod,
                 data = stan_in,
                 chains = nChains,
                 iter = nMCMC,
                 cores=nCores, 
                 seed=42, 
                 thin=nThin, 
                 pars=PARS)
Finish = Sys.time()
Finish-Start

#################
##  Save Model ##
#################

saveRDS(fit_imp_blue, file = "./03_Data/ZZ_Temp/fit_imp.rds")
