# Extract model parameters and posterior predictions
##############################################################

fit_het_blue <- readRDS(file = "./03_Data/ZZ_Temp/fit_full_bayes.rds")
load('./03_Data/ZZ_Temp/panel_het_blue.RData')
G <- panel_list$G

############################
## Extract the parameters ##
############################

phi_a <- rstan::extract(fit_full_bayes, par='phi_a')[[1]]
phi_0 <- rstan::extract(fit_full_bayes, par='phi_0')[[1]]
mu <- rstan::extract(fit_full_bayes, par='Mu')[[1]]

##########################################
## Posterior predictive counterfactuals ##
##########################################

PosteriorPredictions = NA*mu
nIter = dim(mu)[1]
nUnits = dim(mu)[2]
nTimes = dim(mu)[3]
for (b in 1:nIter) {
  for (i in 1:nUnits) {
    for (t in 1:nTimes) {
      if (t < G[i]) {
          PosteriorPredictions[b,i,t] = rnbinom( 1, mu=mu[b,i,t], size=phi_0[b])
        } else { 
          PosteriorPredictions[b,i,t] = rnbinom( 1, mu=mu[b,i,t], size=phi_a[b])
      }    
    }
  }
}

##########
## Save ##
##########

save(PosteriorPredictions, 
     file = "./03_Data/ZZ_Temp/posterior_predictions_full_bayes.RData")

