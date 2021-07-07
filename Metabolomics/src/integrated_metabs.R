Treatment      <- c(rep(c("U5", "U24", "P5", "P24"), each = 2))
Treatment <- factor(Treatment, levels=c("U5", "U24", "P5", "P24"))
Model  <- c(rep(c("Not Integrated", "Integrated")))
Metabolites <- c(10,8,15,14,6,10,7,7)
Data      <- data.frame(Treatment, Model, Metabolites)
library(ggplot2)
library(viridis)
ggplot(Data, aes(fill=Model, y=Metabolites, x=Treatment, label=Metabolites)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity") + 
  geom_text(size = 4, position = position_stack(vjust = 0.5, reverse=TRUE)) +
  theme(text = element_text(size = 18)) + labs(x="Model Condition", y="Number of meatbolites") +
  guides(fill = guide_legend(reverse = TRUE))

