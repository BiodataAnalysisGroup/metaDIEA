library(xgboost)
library(magrittr)
library(dplyr)
library(Matrix)
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
# enconding
label <- ifelse(test = finalData[,'DE'] == TRUE, yes = 1, no = 0)
n = nrow(finalData)
#XGBoost Model Parameters
constant_params <- list(booster = 'gbtree', 
                        silent = 0, 
                        objective = 'binary:logistic', 
                        eval_metric = 'error', 
                        gamma = 0,
                        early.stop.round = 50
)
tunable_params_grid <- expand.grid(eta = (1:6)*0.05,
                                   nrounds = (1:6)*150,
                                   max_depth = c(4:9),
                                   subsample = c(0.6:1),
                                   colsample_bytree = c(0.5,0.8,1)
)
colnames(tunable_params_grid) <- c("eta", "nrounds", "max_depth", "subsample", "colsample_bytree")
currParams <- tunable_params_grid[sample(nrow(tunable_params_grid), size = 1), ] #sample random combination and train the trees
numModels <- 20; #increase
time.taken <- data.frame(n = numeric(), time.t = numeric())
# measure time while decreasing training set
for (i in 1:numModels){
  print(paste0("building model", i))
  train.index = sample(n, 0.99*n)
  train.data = as.matrix(sapply(finalData[train.index,-ncol(finalData)], as.numeric))
  train.label = label[train.index]
  test.data = as.matrix(sapply(finalData[-train.index,-ncol(finalData)], as.numeric))
  test.label = label[-train.index]
  trainm = xgb.DMatrix(data=train.data,label=train.label)
  testm = xgb.DMatrix(data=test.data,label=test.label)
  test.label <- as.factor(ifelse(test = test.label == 1, yes = TRUE, no = FALSE))
  start.t <- Sys.time()
  # Train Model
  xgbModel <- xgb.train(params = constant_params,
                        eta = currParams$eta,
                        nrounds = currParams$nrounds,
                        max_depth = currParams$max_depth,
                        subsample = currParams$subsample,
                        colsample_bytree = currParams$colsample_bytree,
                        data = trainm,
                        verbose = 0)
  end.t <- Sys.time()
  time.taken[i,1] <- length(train.index)
  time.taken[i,2] <- end.t - start.t
  # Predict
  xgb.pred <- predict(xgbModel, test.data)
  xgb.pred <- as.data.frame(xgb.pred)
  pred.class <- ifelse(test = xgb.pred > 0.5, yes = TRUE, no = FALSE)
  # Evaluate model
  confusionMatrix(data = as.factor(pred.class), reference = test.label, positive = 'TRUE')
  # accuracy
  acc <- Accuracy(y_pred = pred.class, y_true = test.label)
  pred_obj = prediction(xgb.pred, test.label, label.ordering = c("FALSE", "TRUE"))
  # auc
  areauc <- performance(prediction.obj = pred_obj, measure = "auc")@y.values
  # f1score
  f1 <- F1_Score(y_true = test.label, y_pred = pred.class, positive = 'TRUE')
  n <- n - 150
}
saveRDS(object = time.taken, file = "./test/xgbDecTrainData.rds")
