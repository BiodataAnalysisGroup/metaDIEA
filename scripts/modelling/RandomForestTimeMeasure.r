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
rndsplt <- sample(x = 2, size = nrow(finalData), replace = TRUE, prob = c(0.95,0.05))
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
numModels <- c(400, 300, 250, 200, 150, 100, 75, 50, 25, 12);
currParams <- as.data.frame(matrix(data = NA, nrow = 1, ncol = ncol(rfGrid) + 3))
colnames(currParams) <- c(colnames(rfGrid), "Accuracy", "AUC", "F1score")
bestGenModel <- currParams
# main loop and measure time
for (j in 1:length(numModels)){
  print(paste0("generating up to ", numModels[j], " models..."))
  start.time <- Sys.time()
  for (i in 1:numModels[j]){
    print(paste0("building model ", i))
    currParams[i,] <- rfGrid[sample(nrow(rfGrid), size = 1),] #sample random combination and train the forest
    # Train Model
    rfModel <- randomForest(x = trainingData[,-ncol(trainingData)], y = trainingData$DE,
                            ntree = currParams[i,'ntree'],
                            replace = currParams[i,'replace'],
                            mtry = currParams[i,'mtry'],
                            sampsize = currParams[i,'sampsize'],
                            nodesize = currParams[i,'nodesize'])
    # Predict
    pred <- predict(object = rfModel, testingData[,-ncol(testingData)], type = 'response')
    predprob <- predict(object = rfModel, testingData[,-ncol(testingData)], type = 'prob')
    # Evaluate model
    # Confusion Matrix
    # confusionMatrix(data = pred, reference = testingData$DE, positive = 'TRUE')
    # Accuracy
    currParams[i,'Accuracy'] <- Accuracy(y_pred = pred, y_true = testingData[,ncol(testingData)])
    pred_obj = prediction(predprob[,2], testingData[,ncol(testingData)], label.ordering = c("FALSE", "TRUE"))
    # AUC
    currParams[i,'AUC'] <- performance(prediction.obj = pred_obj, measure = "auc")@y.values
    # F1 score
    currParams[i,'F1score'] <- F1_Score(y_true = testingData$DE, y_pred = pred, positive = 'TRUE')
  }
  # best generated model among numModels[i] models
  bestModel <- currParams[which(currParams$Accuracy == max(currParams$Accuracy)),]
  bestModel <- bestModel[which(bestModel$AUC == max(bestModel$AUC)),]
  bestModel <- bestModel[which(bestModel$F1score == max(bestModel$F1score)),]
  bestGenModel[j,] <- bestModel[1,] 
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  saveRDS(object = data.frame(n = numModels[j], time.t = time.taken), file = paste0("./test/rf", numModels[j], ".rds"))
}
