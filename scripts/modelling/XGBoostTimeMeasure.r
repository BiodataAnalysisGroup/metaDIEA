library(xgboost)
library(magrittr)
library(dplyr)
library(Matrix)
library(ROCR)
library(MLmetrics)
library(caret)
# read data
finalData <- read.csv(file = "./FinalData.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
finalData$target_id <- as.character(finalData$target_id)
# keep away training data with class = NA
classNaData <- finalData[which(is.na(finalData$DE)),]
# clear data - keep only ones with meaningful classes
finalData <- finalData[which(!is.na(finalData$DE)),]
# replace NAs with 0
finalData[is.na(finalData)] <- 0
# convert character columns to factors
finalData[sapply(finalData, is.character)] <- lapply(finalData[sapply(finalData, is.character)],as.factor)
# remove factors with too many levels
finalData <- finalData[, which(lapply(X = finalData, FUN = nlevels) < 50)]
# convert boolean DE column to factor
finalData$DE <- as.factor(finalData$DE)
# class imbalance problem
hist(as.numeric(finalData$DE) - 1, col = 'blue', main = 'Class Imbalance', xlab = 'FALSE = 0 \t TRUE = 1', ylim = c(0,length(finalData$DE)))
# set seed
set.seed(123)
# shuffle data
finalData <- finalData[sample(nrow(finalData)),]
# split to training and testing datasets
rndsplt <- sample(x = 2, size = nrow(finalData), replace = TRUE, prob = c(0.95,0.05))
trainingData <- finalData[rndsplt == 1, ]
testingData <- finalData[rndsplt == 2, ]
# One-Hot Encoding
label <- ifelse(test = finalData[,'DE'] == TRUE, yes = 1, no = 0)
n = nrow(finalData)
train.index = sample(n,floor(0.70*n))
train.data = as.matrix(sapply(finalData[train.index,], as.numeric))
train.label = label[train.index]
test.data = as.matrix(sapply(finalData[-train.index,], as.numeric))
test.label = label[-train.index]
xgb.train = xgb.DMatrix(data=train.data,label=train.label)
xgb.test = xgb.DMatrix(data=test.data,label=test.label)
test.label <- as.factor(ifelse(test = test.label == 1, yes = TRUE, no = FALSE))
# XGBoost Model Parameters
constant_params <- list(booster = 'gbtree', 
                        silent = 0, 
                        objective = 'binary:logistic', 
                        eval_metric = 'error', 
                        gamma = 0,
                        early.stop.round = 50
)
# tunable parameters
tunable_params_grid <- expand.grid(eta = (1:6)*0.05,
                                   nrounds = (1:6)*150,
                                   max_depth = c(4:9),
                                   subsample = c(0.6:1),
                                   colsample_bytree = c(0.5,0.8,1)
)
colnames(tunable_params_grid) <- c("eta", "nrounds", "max_depth", "subsample", "colsample_bytree")
# watchlist <- list(train = xgb.train, test = xgb.test)
numModels <- c(400, 300, 250, 200, 150, 100, 75, 50, 25, 12);
currParams <- as.data.frame(matrix(data = NA, nrow = 1, ncol = ncol(tunable_params_grid) + 3))
colnames(currParams) <- c(colnames(tunable_params_grid), "Accuracy", "AUC", "F1score")
bestGenModel <- currParams
# main loop and measure time
for (j in 1:length(numModels)){
  print(paste0("generating up to ", numModels[j], " models..."))
  start.time <- Sys.time()
  for (i in 1:numModels[j]){
    print(paste0("building model ", i))
    currParams[i, ] <- tunable_params_grid[sample(nrow(tunable_params_grid), size = 1), ] #sample random combination and train the trees
    # Train Model
    xgbModel <- xgb.train(params = constant_params,
                          eta = tunable_params_grid[i, 'eta'],
                          nrounds = tunable_params_grid[i, 'nrounds'],
                          max_depth = tunable_params_grid[i, 'max_depth'],
                          subsample = tunable_params_grid[i, 'subsample'],
                          colsample_bytree = tunable_params_grid[i,'colsample_bytree'],
                          data = xgb.train)
                          # watchlist = watchlist)
    # Predict
    xgb.pred <- predict(xgbModel, test.data)
    xgb.pred <- as.data.frame(xgb.pred)
    pred.class <- ifelse(test = xgb.pred > 0.5, yes = TRUE, no = FALSE)
    # Evaluate model
    # Confusion Matrix
    # confusionMatrix(data = as.factor(pred.class), reference = test.label, positive = 'TRUE')
    currParams[i, 'Accuracy'] <- Accuracy(y_pred = pred.class, y_true = test.label)
    pred_obj = prediction(xgb.pred, test.label, label.ordering = c("FALSE", "TRUE"))
    currParams[i, 'AUC'] <- performance(prediction.obj = pred_obj, measure = "auc")@y.values
    currParams[i, 'F1score'] <- F1_Score(y_true = test.label, y_pred = pred.class, positive = 'TRUE')
  }
  # best generated model among numModels[i] models
  bestModel <- currParams[which(currParams$Accuracy == max(currParams$Accuracy)),]
  bestModel <- bestModel[which(bestModel$AUC == max(bestModel$AUC)),]
  bestModel <- bestModel[which(bestModel$F1score == max(bestModel$F1score)),]
  bestGenModel[j,] <- bestModel[1,] 
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  saveRDS(object = data.frame(n = numModels[j], time.t = time.taken), file = paste0("./test/xgboost", numModels[j], ".rds"))
}
