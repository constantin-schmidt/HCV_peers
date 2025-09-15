# Extract model paramters and posterior predictions
##############################################################

fit_imp <- readRDS(file = "./03_Data/ZZ_Temp/fit_imp.rds")

############################
## Extract the parameters ##
############################

phi = rstan::extract(fit_imp, par='phi')[[1]]
mu = rstan::extract(fit_imp, par='Mu')[[1]]

##########################################
## Posterior predictive counterfactuals ##
##########################################

Y0_imp = NA*mu
nIter = dim(mu)[1]
nUnits = dim(mu)[2]
nTimes = dim(mu)[3]
for (b in 1:nIter) {
  for (i in 1:nUnits) {
    for (t in 1:nTimes) {
      Y0_imp[b,i,t] = rnbinom(1, mu=mu[b,i,t], size=phi[b])
    }
  }
}

##########
## Save ##
##########

save(Y0_imp, file = "./03_Data/ZZ_Temp/posterior_predictions_imp.RData")
