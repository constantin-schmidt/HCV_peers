##  Plots tau_e for rho=0 ##
############################

library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)

######################
##  Load other data ##
######################
load("./03_Data/ZZ_Temp/panel_het_blue.RData")
load("./03_Data/ZZ_Temp/HEP_dta_short.RData")
load("./03_Data/ZZ_Temp/posterior_predictions_copula_list_full_bayes.RData")
Y_0 <- posterior_predictions_copula_list[["rho=0"]]

nIter <- dim(Y_0)[1]
nUnits <- dim(Y_0)[2]
nTimes <- dim(Y_0)[3]

####################
##  Prepare data  ##
####################

Y_array <- array(rep(panel_list$Y, nIter), dim = c(nUnits, nTimes, nIter)) |>
  aperm(c(3, 1, 2))

treat_array <- Y_array - Y_0

cum_peers <- panel_list$cum_Peers
max_cum_peers <- max(HEP.dta.short$cum_peers)

## Summarize results by cumulative exposure ##
result_matrix <- matrix(0, nrow = nIter, ncol = max_cum_peers)
for (k in 1:nIter) {
  current_matrix <- treat_array[k, , ]
  means_for_slice <- numeric(max_cum_peers)
  for (val in 1:max_cum_peers) {
    means_for_slice[val] <- mean(current_matrix[cum_peers==val])
  }
  result_matrix[k,] <- means_for_slice
}

##  Summarize results ##
mean_vector <- colMeans(result_matrix)
lower_vector <- apply(result_matrix, 2, function(x) quantile(x, 0.025, na.rm=T))
upper_vector <- apply(result_matrix, 2, function(x) quantile(x, 0.975, na.rm=T))

## Count the number of untis that are in treatment at least as long ##
n_in_treat <- sapply(1:max_cum_peers, function(n) 
  sum(panel_list$cum_Peers[ , ncol(panel_list$cum_Peers)] >= n))

##  Plot  ##
xaxis_length <- max_cum_peers-3
plot_data <- data.frame(
  x = 0:xaxis_length,
  mean = c(0,mean_vector[1:(xaxis_length)]),
  lower = c(0,lower_vector[1:(xaxis_length)]),
  upper = c(0,upper_vector[1:(xaxis_length)])
)

# Define breaks and labels for the primary and secondary axes
primary_breaks <- c(seq(0, max(plot_data$x), by = 5), xaxis_length)
secondary_breaks <- c(1, which(diff(n_in_treat) != 0) + 1)

primary_labels <- paste0(primary_breaks)
secondary_labels <- n_in_treat[secondary_breaks]

## plot ##
p <- ggplot(plot_data, aes(x = x)) +
  geom_line(aes(y = mean, color = "Intervention effect")) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#1f78b4", alpha = 0.2) +
  labs(
    title = NULL,
    x = "Cumulative number of peer-months",
    y = "Number of additional DDA eligible HCV patients identified"
    ) +
  scale_color_manual(name = "", values = c("Intervention effect" = "#1f78b4")) +
  theme_minimal() +  
  theme(
    panel.grid.major = element_blank(),    
    panel.grid.minor = element_blank(),    
    axis.line = element_line(color = "black"),  
    legend.position = c(0.15, 0.85),  
    legend.direction = "horizontal",
    legend.background = element_blank(),  
    legend.box.background = element_blank()
  ) +
  geom_vline(aes(xintercept = 0, linetype = "Intervention start"), color = "#e31a1c") +
  scale_linetype_manual(name = "", values = c("Intervention start" = "dotted")) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  # Add a secondary axis with different breaks and labels
  scale_x_continuous(
    breaks = primary_breaks, 
    labels = primary_labels, 
    sec.axis = sec_axis(
      ~ .,  
      breaks = secondary_breaks,  
      labels = secondary_labels,  
      name = "Number of Operation Delivery Networks exposed for at least that many peer-months"  
    )
  )

print(p)
ggsave(
  filename = paste0("./05_ResultsAndDiagnostics/results_section_graphs/sample_average_intervention_effects_by_peer_months.pdf"), 
  plot = p, 
  width = 9, 
  height = 4, 
  units = "in"
)
