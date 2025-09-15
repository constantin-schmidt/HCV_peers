##  IIEs (tauits) with rho = 0 for the full Bayesian model  ##
##############################################################

library(ggplot2)

################
##  Load data ##
################

load('./03_Data/ZZ_Temp/posterior_predictions_copula_list_full_bayes.RData')
Y_0 <- posterior_predictions_copula_list[["rho=0"]]

load('./03_Data/ZZ_Temp/panel_het_blue.RData')

############################
## Generate the IEE array ##
############################

nIter <- dim(Y_0)[1]
nUnits <- dim(Y_0)[2]
nTimes <- dim(Y_0)[3]

Y_array <- array(rep(panel_list$Y, nIter), dim = c(nUnits, nTimes, nIter)) |>
  aperm(c(3, 1, 2))

treat_array <- Y_array - Y_0

####################################
##  Calculate mean and quantiles  ##
####################################

plot.dta <- HEP.dta.short

mean_te <- apply(treat_array, c(2, 3), function(x) mean(x)) |>
  t() |>
  as.vector()

quantile_025 <- apply(treat_array, c(2, 3), function(x) quantile(x, 0.025)) |>
  t() |>
  as.vector()

quantile_975 <- apply(treat_array, c(2, 3), function(x) quantile(x, 0.975)) |>
  t() |>
  as.vector()

# Add the computed values to plot.dta
plot.dta$mean_individual_TEs <- mean_te
plot.dta$quantile_025 <- quantile_025
plot.dta$quantile_975 <- quantile_975

# Create a new column indicating whether the conditions for color change are met
plot.dta$color_condition <- ifelse(plot.dta$quantile_025 > 0 |
                                     plot.dta$quantile_975 < 0, 
                                   "Significant", "Not significant")

# Filter ##
filtered_data <- subset(plot.dta, plot.dta$Dit==1)
filtered_data$Month <- factor(filtered_data$Month, 
                              levels = unique(filtered_data$Month[order(filtered_data$Month)]))

############
##  Plot  ##
############
x_breaks <- c(1, 14, 28, 41)

p <- ggplot(filtered_data, aes(x = Month, y = mean_individual_TEs, color = color_condition)) +
  geom_rect(data = filtered_data %>% 
              filter(Month %in% c("2020-03", "2020-04", "2020-05")), 
            aes(xmin = as.numeric(Month) - 0.5, 
                xmax = as.numeric(Month) + 0.5, 
                ymin = -Inf, ymax = Inf, fill = "Covid-19 lockdown"), 
            alpha = 0.03, inherit.aes = FALSE) +
  geom_point(size = 2, show.legend = FALSE) + 
  geom_hline(yintercept = 0, linetype = "dotted", color = "red") + 
  labs(x = "Month",
       y = "Number of additional DDA eligible \nHCV patients identified",
       title = "(a)") +
  scale_color_manual(name = "", 
                     values = c("Significant" = "#e31a1c",
                                "Not significant" = "#1f78b4"),
                     guide = "none") +  
  scale_fill_manual(
    name = "", 
    values = c("Covid-19 lockdown" = "#33a02c"), 
    guide = guide_legend(override.aes = list(alpha = 0.4)) 
  ) +
  scale_x_discrete(breaks = levels(filtered_data$Month)[x_breaks]) +
  theme_minimal() +
  theme(
    legend.position = c(0.3, 0.9),
    axis.text.x = element_text(size = 8, hjust = .75, margin = margin(t = 5, r = 0, b = 0, l = 0)),
  )

## Print the plot ##
print(p)

## Save ##
ggsave(
  filename = paste0("./05_ResultsAndDiagnostics/results_section_graphs/mean_ITEs_for_full_Bayesian_model_with_rho_0.pdf"), 
  plot = p, 
  width = 4.5, 
  height = 4, 
  units = "in"
)
