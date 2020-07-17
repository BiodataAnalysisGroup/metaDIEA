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
# set seed
set.seed(123)
# shuffle data
finalData <- finalData[sample(nrow(finalData)),]
# split to training and testing datasets
rndsplt <- sample(x = 2, size = nrow(finalData), replace = TRUE, prob = c(0.99,0.01))
trainingData <- finalData[rndsplt == 1, ]
testingData <- finalData[rndsplt == 2, ]
numModels <- 20 # increase
time.taken <- data.frame(n = numeric(), time.t = numeric())
n <- nrow(trainingData)
# search grid
rfGrid <- expand.grid(ntreeArg = (1:5)*500,
                      replaceAge = c(TRUE, FALSE),
                      mtryArg = (round(sqrt(length(finalData))) + (-2:2)),
                      nodesize = c(1, 3, 5, 7)
)
colnames(rfGrid) <- c("ntree", "replace", "mtry", "nodesize")
currParams <- rfGrid[sample(nrow(rfGrid), size = 1),] #sample random combination and train the forest
# measure time while decreasing training
for (i in 1:numModels){
  print(paste0("building model", i))
  # indx <- sample(n, floor((1/i)*n))
  trainData <- trainingData[c(1:n),]
  sampsizeArg  = c(nrow(trainData), round(nrow(trainData)*.9), round(nrow(trainData)*.8), round(nrow(trainData)*.75))
  samp <- sampsizeArg[sample(length(sampsizeArg), 1)]
  start.t <- Sys.time()
  # Train Model
  rfModel <- randomForest(x = trainData[,-ncol(trainData)], y = trainData$DE,
                          ntree = currParams$ntree,
                          replace = currParams$replace,
                          mtry = currParams$mtry,
                          sampsize = samp,
                          nodesize = currParams$nodesize)
  end.t <- Sys.time()
  time.taken[i,1] <- nrow(trainData)
  time.taken[i,2] <- end.t - start.t
  # Predict
  pred <- predict(object = rfModel, testingData[,-ncol(testingData)], type = 'response')
  predprob <- predict(object = rfModel, testingData[,-ncol(testingData)], type = 'prob')
  # Evaluate model
  # Confusion Matrix
  # confusionMatrix(data = pred, reference = testingData$DE, positive = 'TRUE')
  # Accuracy
  acc <- Accuracy(y_pred = pred, y_true = testingData[,ncol(testingData)])
  pred_obj = prediction(predprob[,2], testingData[,ncol(testingData)], label.ordering = c("FALSE", "TRUE"))
  # AUC
  areauc <- performance(prediction.obj = pred_obj, measure = "auc")@y.values
  # F1 score
  f1 <- F1_Score(y_true = testingData$DE, y_pred = pred, positive = 'TRUE')
  n <- n - 150
}
saveRDS(object = time.taken, file = "./test/rfDecTrainData.rds")
