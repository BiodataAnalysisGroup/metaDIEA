library(randomForest)
library(ROCR)
library(MLmetrics)
library(caret)
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
# class imbalance problem
# hist(as.numeric(finalData$DE) - 1, col = 'blue', main = 'Class Imbalance', xlab = 'FALSE = 0 \t TRUE = 1', ylim = c(0,length(finalData$DE)))
# set seed
numModels <- 20 # increase
# search grid
rfGrid <- expand.grid(ntreeArg = (1:5)*500,
                      replaceAge = c(TRUE, FALSE),
                      mtryArg = (round(sqrt(length(finalData))) + (-2:2)),
                      nodesize = c(1, 3, 5, 7)
)
colnames(rfGrid) <- c("ntree", "replace", "mtry", "nodesize")
currParams <- as.data.frame(matrix(data = NA, nrow = 1, ncol = ncol(rfGrid) + 6))
colnames(currParams) <- c(colnames(rfGrid), "sampsizeArg" ,"Accuracy", "AUC", "F1score", "trainingLength", "testingLength")
currParams <- rfGrid[sample(nrow(rfGrid), size = 1),] #sample random combination and train the forest
currParams[c(1:20),c(1:4)] <- currParams
split <- .99
# generate 20 models and measure metrics
for (i in 1:numModels){
  print(paste0("building model", i))
  set.seed(123)
  # shuffle data
  finalData <- finalData[sample(nrow(finalData)),]
  # split to training and testing datasets
  rndsplt <- sample(x = 2, size = nrow(finalData), replace = TRUE, prob = c(split, 1-split))
  trainingData <- finalData[rndsplt == 1, ]
  testingData <- finalData[rndsplt == 2, ]
  # select arguments
  sampsizeArg  = c(nrow(trainingData), round(nrow(trainingData)*.9), round(nrow(trainingData)*.8), round(nrow(trainingData)*.75))
  samp <- sampsizeArg[sample(length(sampsizeArg), 1)]
  currParams[i, "sampsizeArg"] <- samp
  # Train Model
  rfModel <- randomForest(x = trainingData[,-ncol(trainingData)], y = trainingData$DE,
                          ntree = currParams[i, "ntree"],
                          replace = currParams[i, "replace"],
                          mtry = currParams[i, "mtry"],
                          sampsize = currParams[i, "sampsizeArg"],
                          nodesize = currParams[i, "nodesize"])
  # Predict
  pred <- predict(object = rfModel, testingData[,-ncol(testingData)], type = 'response')
  predprob <- predict(object = rfModel, testingData[,-ncol(testingData)], type = 'prob')
  # Evaluate model
  # Confusion Matrix
  # confusionMatrix(data = pred, reference = testingData$DE, positive = 'TRUE')
  # Accuracy
  currParams[i, "Accuracy"] <- Accuracy(y_pred = pred, y_true = testingData[,ncol(testingData)])
  pred_obj = prediction(predprob[,2], testingData[,ncol(testingData)], label.ordering = c("FALSE", "TRUE"))
  # AUC
  currParams[i, "AUC"] <- performance(prediction.obj = pred_obj, measure = "auc")@y.values
  # F1 score
  currParams[i, "F1score"] <- F1_Score(y_true = testingData$DE, y_pred = pred, positive = 'TRUE')
  currParams[i, "trainingLength"] <- nrow(trainingData)
  currParams[i, "testingLength"] <- nrow(testingData)
  split <- split - .05
}
saveRDS(object = currParams[,c(6:10)], file = "./test/rfAccAucF1.rds")
