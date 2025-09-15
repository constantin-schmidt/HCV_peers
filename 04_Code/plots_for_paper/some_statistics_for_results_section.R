##  Computes some of the statistics mentioned in the results section  ##
########################################################################

##########################
##  Prob phi_a > phi_0  ##
##########################

fit_full_bayes <- readRDS(file = "./03_Data/ZZ_Temp/fit_full_bayes.rds")

phi_a <- rstan::extract(fit_full_bayes, par='phi_a')[[1]]
phi_0 <- rstan::extract(fit_full_bayes, par='phi_0')[[1]]

print(paste("The probability that phi_a > phi_0 is:", mean(phi_a < phi_0)))

######################################################
##  Number of exposed observations during lockdown  ##
######################################################

load('./03_Data/ZZ_Temp/panel_het_blue.RData')

print(paste("The number of exposed observations overall",
            sum(panel_list$D)))
print(paste("The number of exposed observations during lockdown",
            sum(panel_list$LDtreat)))