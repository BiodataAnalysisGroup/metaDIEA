library(xgboost)
library(magrittr)
library(dplyr)
library(Matrix)
library(ROCR)
library(MLmetrics)
library(caret)
library(hrbrthemes)
library(ggplot2)
# read data
finalData <- read.csv(file = "./hg19Final.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
colnames(finalData)[1] <- "target_id"
colnames(finalData)[15] <- "DE"
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
# class distribution
fig <- ggplot(data = finalData, aes(x = DE)) +
  geom_histogram(stat = "count", fill = "turquoise", color = "purple", binwidth = 5) +
  ggtitle("Class Distribution") +
  xlab(label = "Differential Expression") +
  ylab("Number of Observations") +
  theme_ipsum() +
  theme(plot.title = element_text(hjust = 0.5)) 
fig
ggsave(filename = "./xgbClassDistribution.png", width = 15, height = 15, units = "cm")
# set seed
set.seed(123)
# shuffle data
finalData <- finalData[sample(nrow(finalData)),]
# One-Hot Encoding
label <- ifelse(test = finalData[,'DE'] == TRUE, yes = 1, no = 0)
n = nrow(finalData)
train.index = sample(n,floor(0.70*n))
train.data = as.matrix(sapply(finalData[train.index,-ncol(finalData)], as.numeric))
train.label = label[train.index]
test.data = as.matrix(sapply(finalData[-train.index,-ncol(finalData)], as.numeric))
test.label = label[-train.index]
xgb.train = xgb.DMatrix(data=train.data,label=train.label)
xgb.test = xgb.DMatrix(data=test.data,label=test.label)
test.label <- as.factor(ifelse(test = test.label == 1, yes = TRUE, no = FALSE))
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
watchlist <- list(train = xgb.train, test = xgb.test)
numModels <- 10; # could increase
currParams <- as.data.frame(matrix(data = NA, nrow = 1, ncol = ncol(tunable_params_grid) + 3))
colnames(currParams) <- c(colnames(tunable_params_grid), "Accuracy", "AUC", "F1score")

for (i in 1:numModels){
  print(paste0("Building model ", i, "..."))
  currParams[i, ] <- tunable_params_grid[sample(nrow(tunable_params_grid), size = 1), ] #sample random combination and train the trees
  #Train Model
  xgbModel <- xgb.train(params = constant_params,
                        eta = currParams[i, 'eta'],
                        nrounds = currParams[i, 'nrounds'],
                        max_depth = currParams[i, 'max_depth'],
                        subsample = currParams[i, 'subsample'],
                        colsample_bytree = currParams[i,'colsample_bytree'],
                        data = xgb.train)
  #Prediction
  xgb.pred <- predict(xgbModel, test.data)
  xgb.pred <- as.data.frame(xgb.pred)
  pred.class <- ifelse(test = xgb.pred > 0.5, yes = TRUE, no = FALSE)
  #Evaluate model
  confusionMatrix(data = as.factor(pred.class), reference = test.label, positive = 'TRUE')
  currParams[i, 'Accuracy'] <- Accuracy(y_pred = pred.class, y_true = test.label)
  pred_obj = prediction(xgb.pred, test.label, label.ordering = c("FALSE", "TRUE"))
  currParams[i, 'AUC'] <- performance(prediction.obj = pred_obj, measure = "auc")@y.values
  currParams[i, 'F1score'] <- F1_Score(y_true = test.label, y_pred = pred.class, positive = 'TRUE')
}
#taking best model parameters
bestModel <- currParams[which(currParams$Accuracy == max(currParams$Accuracy)), ]
bestModel <- bestModel[which(bestModel$AUC == max(bestModel$AUC)), ]
bestParams <- bestModel[1,c(1:5)]
#final model
xgbModel <- xgb.train(params = constant_params,
                      eta = bestParams$eta,
                      nrounds = bestParams$nrounds,
                      max_depth = bestParams$max_depth,
                      subsample = bestParams$subsample,
                      colsample_bytree = bestParams$colsample_bytree,
                      data = xgb.train)
#Prediction
xgb.pred <- predict(xgbModel, test.data)
xgb.pred <- as.data.frame(xgb.pred)
pred.class <- ifelse(test = xgb.pred > 0.5, yes = TRUE, no = FALSE)
#Evaluate model
confusionMatrix(data = as.factor(pred.class), reference = test.label, positive = 'TRUE')
ROCcurve <- performance(pred_obj, "tpr", "fpr")
plot(ROCcurve, col = "blue")
abline(0, 1, col = "grey")
Accuracy <- Accuracy(y_pred = pred.class, y_true = test.label)
pred_obj = prediction(xgb.pred, test.label, label.ordering = c("FALSE", "TRUE"))
AUC <- performance(prediction.obj = pred_obj, measure = "auc")@y.values