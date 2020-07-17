library(ggplot2)
library(dplyr)
library(hrbrthemes)
numModels <- c(12, 25, 50, 75, 100, 150, 200, 250, 300, 400)
rfModels <- data.frame(n = numeric(), time.t = numeric())
xgbModels <- data.frame(n = numeric(), time.t = numeric())
# read results
for(i in 1:length(numModels)){
  rfModels[i,] <- readRDS(file = paste0("./test/t", numModels[i], ".rds"))
  xgbModels[i,] <- readRDS(file = paste0("./test/xgboost", numModels[i], ".rds"))
}
# convert to minutes
rfModels[c(8:10),2] <- rfModels[c(8:10),2]*60
xgbModels[c(1:2),2] <- xgbModels[c(1:2),2]/60
# plot rf time
p1 <- ggplot(aes(x = n, y = time.t), data = rfModels)+
        geom_line(color = "purple4", size = 1) +
        geom_point(shape=20, color="grey69", fill="#69b3a2", size=4) +
        xlab("Number of models") +
        ylab("Time (in minutes)") +
        ggtitle("Random Forest Classifier - Computation Time.")
p1 + theme_bw()
# save plot
ggsave(filename = "./test/rfPlot.png", width = 15, height = 15, units = "cm")
# plot xgboost time    
p2 <- ggplot(aes(x = n, y = time.t), data = xgbModels) +
        geom_line(color = "purple4", size = 1) +
        geom_point(shape=20, color="grey69", fill="#69b3a2", size=4) +
        xlab("Number of models") +
        ylab("Time (in minutes)") +
        ggtitle("xGBoost Classifier - Computation Time.")
p2 + theme_bw()
# save plot
ggsave(filename = "./test/xgbPlot.png", width = 15, height = 15, units = "cm")
