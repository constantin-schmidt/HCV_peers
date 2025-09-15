## Prepare and format the data for the stan models
#####################################################################################
rm(list = ls())

library(tibble)
library(tidyr)
library(dplyr)
library(rlang)

#--------------------------------
##  Load data
load("./03_Data/HEP_dta.RData")

#--------------------------------
# Generate the vectors & matrices

## Treatment IDs
ODN_id <- unique(HEP.dta$ODN) |>
  as.factor() |>
  as.integer()

## Outcomes

##  received blueteq  ##
received_blueteq <- HEP.dta |>
  dplyr::select(received_blueteq, ODN, Month) |>
  pivot_wider(names_from = Month, values_from = received_blueteq) |>
  column_to_rownames(var="ODN") |>
  as.matrix() 

## Matrix of current exposure (number of peers)
Peers <- HEP.dta |>
  dplyr::select(NumberOfPeers, ODN, Month) |>
  pivot_wider(names_from = Month, values_from = NumberOfPeers) |>
  column_to_rownames(var="ODN") |>
  as.matrix()

## Matrix of cumulative exposure up to t (cumulative peer-month)
  
  ## Create cumulative exposure
  HEP.dta <- HEP.dta |>
    group_by(ODN) |>
    mutate(cum_peers = cumsum(NumberOfPeers))

cum_Peers <- HEP.dta |>
  dplyr::select(cum_peers, ODN, Month) |>
  pivot_wider(names_from = Month, values_from = cum_peers) |>
  column_to_rownames(var="ODN") |>
  as.matrix()

## Matrix of the treatment
D <- HEP.dta |>
  dplyr::select(Dit, ODN, Month) |>
  pivot_wider(names_from = Month, values_from = Dit) |>
  column_to_rownames(var="ODN") |>
  as.matrix()

## Vector of treatment start periods
G <- HEP.dta |>
  group_by(ODN) |>
  mutate(TimePeriod = 1:n()) |>
  slice(which.max(Dit)) |>
  mutate(G = if_else(TimePeriod == 1, ncol(D)+1, TimePeriod)) |>
  column_to_rownames(var="ODN") |>
  select(G) |>
  as.matrix()
  # The start period of ODNs that start treatment after the study period is 
  # set to the maximum time period + 1

##  Covid lockdown / treatment interaction  ##
LD <- HEP.dta |>
  dplyr::select(first_lockdown, ODN, Month) |>
  pivot_wider(names_from = Month, values_from = first_lockdown) |>
  column_to_rownames(var="ODN") |>
  as.matrix()

LDtreat <- LD*D

##  Matrix of time in treatment
TreatTime <- HEP.dta |>
  group_by(ODN) |>
  mutate(TimePeriod = 1:n(),
         G = min(TimePeriod[Dit == 1]),
         TimeInTreat = if_else(Dit == 0,
                               0,
                               TimePeriod - G)) |>
  dplyr::select(TimeInTreat, ODN, Month) |>
  pivot_wider(names_from = Month, values_from = TimeInTreat) |>
  column_to_rownames(var="ODN") |>
  as.matrix()
# Warnings can be ignored, as the are circumvented by the if_else statement

##  Spline  ##  
  ##  create the spline base  ##
  input_vec <- HEP.dta$cum_peers ## Decide the variable to base spline on
  n_knots <- 3                         ## Decide number of knots
  knots <- quantile(input_vec[HEP.dta$Dit==1], 
                    probs = seq(0, 1, length.out = n_knots + 2))[1:n_knots+1]
  spline_base_long <- bs(input_vec, degree = 3, knots = knots,
                         intercept = F)
  
  basis_df <- as.data.frame(spline_base_long)
  names(basis_df) <- paste("Basis", 1:ncol(basis_df), sep = "_")
  
  # Plot each of the basis splines
  matplot(input_vec, basis_df, type = "l", lty = 1, lwd = 2,
          col = rainbow(ncol(basis_df)),
          xlab = "x", ylab = "Spline Basis",
          main = "B-spline Basis Functions")
  legend("topright", legend = names(basis_df), col = rainbow(ncol(basis_df)), lwd = 2)
  
  
  ##  transform into matrix ##
  n_spline_parameters <- ncol(spline_base_long)
  spline_base_array <- array(NA,
                             dim = c(nrow(D),
                                     ncol(D),
                                     n_spline_parameters))
  for (i in 1:n_spline_parameters) {
    spline_base_array[,,i] <- matrix(spline_base_long[,i], 
                                     nrow = nrow(D), 
                                     ncol = ncol(D),
                                     byrow = T)
  }

#--------------------------------
##  First time of peer availability
avail_tp <- min(G)
  
#------------------------------
# Create the list
panel_list <- list("ODN_id" = ODN_id,
                   "received_blueteq" = received_blueteq,
                   "D" = D,
                   "G" = G,
                   "Peers" = Peers,
                   "cum_Peers" = cum_Peers,
                   "TreatTime" = TreatTime,
                   "LDtreat" = LDtreat,
                   "avail_tp" = avail_tp,
                   "spline_base" = spline_base_array)

#------------------------------
## set main outcome ##
Outcome_var <- "received_blueteq"
panel_list$Outcome_var <- Outcome_var
panel_list[["Y"]] <- panel_list[[Outcome_var]]

#------------------------------
##  Save
save(panel_list, file = "./03_Data/ZZ_Temp/panel_het_blue.RData")
