##  Some prior predictive checks
##################################################################

load('./03_Data/HEP_dta.RData')
load('./03_Data/ZZ_Temp/panel_het_blue.RData')

##----------------------------
##  phi

## Unexposed observations ##
sqrt_disp_0 = abs(rnorm(100000, 0, .11))
phi_0 = 1/sqrt_disp_0^2
mu = mean(HEP.dta[[panel_list$Outcome_var]]
          [HEP.dta$Dit==0])
hist(1 + mu/phi_0, breaks = 100)
mean(1 + mu/phi_0 > 3)

## Exposed observations ##
sqrt_disp_a = abs(rnorm(100000, 0, .1165))
phi_a = 1/sqrt_disp_a^2
mu = mean(HEP.dta[[panel_list$Outcome_var]]
          [HEP.dta$Dit==1])
hist(1 + mu/phi_a, breaks = 100)
mean(1 + mu/phi_a > 3)
