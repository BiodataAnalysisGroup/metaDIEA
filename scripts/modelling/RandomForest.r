library(randomForest)
library(ROCR)
library(MLmetrics)
library(caret)
library(hrbrthemes)
library(ggplot2)
# read data
finalData <- read.csv(file = "./FinalData.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
finalData$target_id <- as.character(finalData$target_id)
# drop possible all NA columns
finalData <- finalData[, !apply(is.na(finalData), 2, all)]
# keep way from training data with class = NA
classNaData <- finalData[which(is.na(finalData$DE)),]
# clear data - keep only ones with meaningful classes
finalData <- finalData[which(!is.na(finalData$DE)),]
# replace NAs with 0
finalData[is.na(finalData)] <- 0
# convert character columns to factors
finalData[sapply(finalData, is.character)] <- lapply(finalData[sapply(finalData, is.character)],as.factor)
# remove factors with many levels
finalData <- finalData[, which(lapply(X = finalData, FUN = nlevels) < 50)]
# convert boolean DE column to factor
finalData$DE <- as.factor(finalData$DE)
# class distribution
fig <- ggplot(data = finalData, aes(x = DE)) +
  geom_histogram(stat = "count", fill = "turquoise", color = "purple", binwidth = 5) +
  ggtitle("Class Distribution") +
  xlab(label = "Differential Expression") +
  ylab("Number of Observations") +
  theme_ipsum() +
  theme(plot.title = element_text(hjust = 0.5)) 
fig
ggsave(filename = "./test/rfClassDistribution.png", width = 15, height = 15, units = "cm")
# set seed
set.seed(123)
# shuffle data
finalData <- finalData[sample(nrow(finalData)),]
# split to training and testing datasets
rndsplt <- sample(x = 2, size = nrow(finalData), replace = TRUE, prob = c(0.7,0.3))
trainingData <- finalData[rndsplt == 1, ]
testingData <- finalData[rndsplt == 2, ]
# create search grid
rfGrid <- expand.grid(ntreeArg = (1:5)*500,
                      replaceAge = c(TRUE, FALSE),
                      mtryArg = (round(sqrt(length(trainingData))) + (-2:2)),
                      sampsizeArg  = c(nrow(trainingData), round(nrow(trainingData)*.9), round(nrow(trainingData)*.8), round(nrow(trainingData)*.75)),
                      nodesize = c(1, 3, 5, 7)
                      )
colnames(rfGrid) <- c("ntree", "replace", "mtry", "sampsize", "nodesize")
numModels <- 10; # increase
currParams <- as.data.frame(matrix(data = NA, nrow = 1, ncol = ncol(rfGrid) + 3))
colnames(currParams) <- c(colnames(rfGrid), "Accuracy", "AUC", "F1score")
# loop for numModels
for (i in 1:numModels){
  currParams[i,] <- rfGrid[sample(nrow(rfGrid), size = 1),] # sample random combination and train the forest
  # Train Model
  rfModel <- randomForest(x = trainingData[,-ncol(trainingData)], y = trainingData$DE,
                          ntree = currParams[i,'ntree'],
                          replace = currParams[i,'replace'],
                          mtry = currParams[i,'mtry'],
                          sampsize = currParams[i,'sampsize'],
                          nodesize = currParams[i,'nodesize'])
  # Prediction
  pred <- predict(object = rfModel, testingData[,-ncol(testingData)], type = 'response')
  predprob <- predict(object = rfModel, testingData[,-ncol(testingData)], type = 'prob')
  # Evaluate model
  # confusion Matrix
  # confusionMatrix(data = pred, reference = testingData$DE, positive = 'TRUE')
  # Accuracy
  currParams[i,'Accuracy'] <- Accuracy(y_pred = pred, y_true = testingData[,ncol(testingData)])
  pred_obj = prediction(predprob[,2], testingData[,ncol(testingData)], label.ordering = c("FALSE", "TRUE"))
  # AUC
  currParams[i,'AUC'] <- performance(prediction.obj = pred_obj, measure = "auc")@y.values
  # F1 score
  currParams[i,'F1score'] <- F1_Score(y_true = testingData$DE, y_pred = pred, positive = 'TRUE')
}
# best model
bestModel <- currParams[which(currParams$Accuracy == max(currParams$Accuracy)),]
bestModel <- bestModel[which(bestModel$F1score == max(bestModel$F1score)),]
bestModel <- bestModel[which(bestModel$AUC == max(bestModel$AUC)),]
# Take best model Parameters
bestParams <- bestModel[1,c(1:5)]
# Stability Analysis permuting FinalData
# Create 10 equally size folds
folds <- cut(seq(1,nrow(finalData)),breaks=10,labels=FALSE)
acc <- list()
# Perform 10 fold cross validation
for(i in 1:10){
  # Randomly shuffle the data
  finalData <- finalData[sample(nrow(finalData)),]
  # Segement data by fold using the which() function 
  Indexes <- which(folds==i,arr.ind=TRUE)
  testingData <- finalData[Indexes, ]
  trainingData <- finalData[-Indexes, ]
  pred <- data.frame(matrix(data = NA, nrow = nrow(testingData), ncol = 1))
  # Train Model
  rfModel <- randomForest(x = trainingData[,-ncol(trainingData)], y = trainingData$DE,
                          ntree = bestParams[1,'ntree'],
                          replace = bestParams[1,'replace'],
                          mtry = bestParams[1,'mtry'],
                          sampsize = bestParams[1,'sampsize'],
                          nodesize = bestParams[1,'nodesize'])
  # Predict
  pred <- predict(object = rfModel, testingData[,-ncol(testingData)], type = 'response')
  predprob <- predict(object = rfModel, testingData[,-ncol(testingData)], type = 'prob')
  # Accuracy
  acc[i] <- Accuracy(y_pred = pred, y_true = testingData[,ncol(testingData)])
}
a <- data.frame(acc = unlist(acc), n = c(1:10))
# Plot
p1 <- ggplot(aes(x = n, y = acc), data = a)+
  geom_line(color = "purple4", size = 1) +
  geom_point(shape=20, color="grey69", fill="#69b3a2", size=4) +
  xlab("Permutation i") +
  xlim(1,10) +
  ylab("Accuracy") +
  ylim(0.9,1.1) +
  ggtitle("Stability Analysis") +
  theme_bw()
p1 
ggsave(filename = "./test/rfStability.png", width = 15, height = 15, units = "cm")
