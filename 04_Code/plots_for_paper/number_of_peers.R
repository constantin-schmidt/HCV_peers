##  Plotting the start times of each peer in the study period ##
################################################################
rm(list = ls())
library(ggplot2)
library(wesanderson)
library(dplyr)

load("./03_Data/HEPana.RData")

####################
##  Prepare data  ##
####################

HEP.plot <- HEP.ana |>
  filter(ODN != "Prison/Detention centre", Month < "2021-06") |>
  mutate(NumberOfPeers = as.factor(NumberOfPeers))

levels <- unique(HEP.plot$ODN[order(HEP.plot$InterventionStart)])
HEP.plot$ODN <- factor(HEP.plot$ODN, levels = levels)

# Custom labels
custom_labels <- rep("", length(levels(HEP.plot$ODN))) 
custom_labels[14] <- "Exposed" 
custom_labels[3] <- "Unexposed"

# Boundary
boundary <- 5.5

############
##  Plot  ##
############

p <- ggplot(HEP.plot, aes(x = Month, y = ODN, fill = NumberOfPeers)) + 
  geom_tile() + 
  theme(
    legend.position = "bottom", 
    legend.justification = "right",
    legend.box = "horizontal",
    axis.text.y = element_text(angle = 90, hjust = 0.5),  
    axis.ticks.y = element_blank(),                       
    axis.text.x = element_text(),                         
    axis.ticks.x = element_blank(),                       
    axis.title.y = element_text(),                        
    axis.title.x = element_text(),                        
    plot.title = element_text(hjust = -0.1, vjust = 2, face = "bold", size = 14), 
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10) 
  ) + 
  scale_y_discrete(labels = custom_labels, limits = rev(levels(HEP.plot$ODN))) +
  scale_x_discrete(
    breaks = c("2017-01", "2018-01", "2019-01", "2020-01", "2021-01", "2022-01"),
    labels = c("Jan 17", "Jan 18", "Jan 19", "Jan 20", "Jan 21", "Jan 22")
  ) + 
  scale_fill_discrete(name = "# peers:") +
  labs(y = "Operation Delivery Networks") + 
  geom_hline(yintercept = boundary, color = "darkred", size = 0.5) + 
  ggtitle("(a)") +
  guides(fill = guide_legend(ncol = 7, nrow = 1))
print(p)

# Save
ggsave(
  filename = "./05_ResultsAndDiagnostics/descriptive_graphs/number_of_peers.pdf", 
  plot = p, 
  width = 4, 
  height = 4, 
  units = "in"
)
