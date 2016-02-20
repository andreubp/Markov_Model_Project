

library(RColorBrewer)
library(ggplot2)
ggplot(df,aes(FPR,TPR,color=GeneSet))+geom_line(size = 1, alpha = 0.9)+
  labs(title= "ROC curve", 
       x = "False Positive Rate (1-Specificity)", 
       y = "True Positive Rate (Sensitivity)")