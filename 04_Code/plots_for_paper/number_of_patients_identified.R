##  Graph showing the number of new DDA eligible patients identified each month
##  in each ODN
################################################################################
rm(list = ls())
library(ggplot2)
library(ggnewscale)

load("./03_Data/HEP_dta.RData")

####################
##  Prepare data  ##
####################

HEP.dta.plot <- HEP.dta |>
  filter(Month < "2021-06")

# Order ODNs
levels <- unique(HEP.dta.plot$ODN[order(HEP.dta.plot$InterventionStart)])
HEP.dta.plot$ODN <- factor(HEP.dta.plot$ODN, levels = levels)

##########################
## Labels and boundary  ##
##########################

# Boundary
boundary <- length(unique(HEP.dta.plot$ODN)) - 
  length(unique(HEP.dta.plot$ODN[HEP.dta.plot$Dit==1]))

# Custom labels
custom_labels <- rep("", length(levels(HEP.dta.plot$ODN))) 
custom_labels[14] <- "Exposed" 
custom_labels[3] <- "Unexposed"
# Label position chosen so it fits well with plot

############
##  Plot  ##
############

p <- ggplot(data = HEP.dta.plot,
            aes(Month, ODN , fill = received_blueteq)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightgreen", high = "darkblue") +
  theme_minimal() +
  scale_x_discrete(breaks = c("2017-01", "2018-01", "2019-01", "2020-01", 
                              "2021-01", "2022-01"),
                   labels = c("Jan 17", "Jan 18", "Jan 19", "Jan 20", 
                              "Jan 21", "Jan 22")) +
  scale_y_discrete(labels = custom_labels, limits = rev(levels(HEP.dta.plot$ODN))) + 
  labs(x = "Month", y = "Operation Delivery Networks", fill = "# patients identified:") +
  geom_hline(yintercept = boundary + 0.5, color = "darkred", size = 0.5) +
  theme(
    axis.text.y = element_text(angle = 90, hjust = 0.5),
    legend.position = "bottom",                
    legend.justification = "center",           
    legend.box = "horizontal",                 
    legend.text = element_text(size = 6),      
    legend.title = element_text(size = 10),    
    plot.title = element_text(hjust = -0.1, vjust = 2, face = "bold", size = 14),
    plot.margin = margin(t = 10, r = 10, b = 0, l = 10) 
  ) +
  ggtitle("(b)")  

##  Save  ##
ggsave(
  filename = "./05_ResultsAndDiagnostics/descriptive_graphs/number_of_patients_identified.pdf", 
       plot = p,
       width = 4,  
       height = 4, 
       units = "in")