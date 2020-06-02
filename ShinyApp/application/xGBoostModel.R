xGBoostModel <- function(finalData, nModels){
  print("XGBoost....")
  # start time
  start.time <- Sys.time()
  # One-Hot Encoding
  label <- ifelse(test = finalData[,'DE'] == TRUE, yes = 1, no = 0)
  n = nrow(finalData)
  train.index = sample(n,floor(0.70*n))
  # Split to training and testing datasets
  train.data = as.matrix(sapply(finalData[train.index,-ncol(finalData)], as.numeric))
  train.label = label[train.index]
  test.data = as.matrix(sapply(finalData[-train.index,-ncol(finalData)], as.numeric))
  test.label = label[-train.index]
  trainm = xgb.DMatrix(data=train.data,label=train.label)
  testm = xgb.DMatrix(data=test.data,label=test.label)
  test.label <- as.factor(ifelse(test = test.label == 1, yes = TRUE, no = FALSE))
  # XGBoost Model Parameters
  # Untunable Parameters
  constant_params <- list(booster = 'gbtree', 
                          silent = 0, 
                          objective = 'binary:logistic', 
                          eval_metric = 'error', 
                          gamma = 0,
                          early.stop.round = 50
  )
  # Search Grid
  tunable_params_grid <- expand.grid(eta = (1:6)*0.05,
                                     nrounds = (1:6)*150,
                                     max_depth = c(4:9),
                                     subsample = c(0.6:1),
                                     colsample_bytree = c(0.5,0.8,1)
  )
  colnames(tunable_params_grid) <- c("eta", "nrounds", "max_depth", "subsample", "colsample_bytree")
  tunable_params_grid <- tunable_params_grid[sample(nrow(tunable_params_grid)),]
  watchlist <- list(train = trainm, test = testm)
  numModels <- nModels; # specified by the user
  currParams <- as.data.frame(matrix(data = NA, nrow = 1, ncol = ncol(tunable_params_grid) + 3))
  colnames(currParams) <- c(colnames(tunable_params_grid), "Accuracy", "AUC", "F1score")
  withProgress(message = paste0("Generating ", numModels, " models"), detail = "May take a while", value = 0, {
    # Generate nModels 
    for (i in 1:numModels){
      currParams[i, ] <- tunable_params_grid[sample(nrow(tunable_params_grid), size = 1), ] #sample random combination and train the trees
      incProgress(1/numModels, detail = paste0("Model ", i))
      # Train Model
      xgbModel <- xgb.train(params = constant_params,
                            eta = currParams[i, 'eta'],
                            nrounds = currParams[i, 'nrounds'],
                            max_depth = currParams[i, 'max_depth'],
                            subsample = currParams[i, 'subsample'],
                            colsample_bytree = currParams[i,'colsample_bytree'],
                            data = trainm,
                            watchlist = watchlist,
                            verbose = 0)
      #Prediction
      xgb.pred <- predict(xgbModel, test.data)
      xgb.pred <- as.data.frame(xgb.pred)
      pred.class <- ifelse(test = xgb.pred > 0.5, yes = TRUE, no = FALSE)
      # Evaluate model
      # confusionMatrix(data = as.factor(pred.class), reference = test.label, positive = 'TRUE')
      # Accuracy
      currParams[i, 'Accuracy'] <- Accuracy(y_pred = pred.class, y_true = test.label)
      pred_obj = prediction(xgb.pred, test.label, label.ordering = c("FALSE", "TRUE"))
      # Area Under the Curve
      currParams[i, 'AUC'] <- performance(prediction.obj = pred_obj, measure = "auc")@y.values
      # F1 measure
      currParams[i, 'F1score'] <- F1_Score(y_true = test.label, y_pred = pred.class, positive = 'TRUE')
    }  
  })
  # Best Model Rarameters
  bestModel <- currParams[which(currParams$Accuracy == max(currParams$Accuracy)), ]
  bestModel <- bestModel[which(bestModel$F1score == max(bestModel$F1score)), ]
  bestModel <- bestModel[which(bestModel$AUC == max(bestModel$AUC)), ]
  bestParams <- bestModel[1,c(1:5)]
  # Final model
  xgbModel <- xgb.train(params = constant_params,
                        eta = bestParams[1, 'eta'],
                        nrounds = bestParams[1, 'nrounds'],
                        max_depth = bestParams[1, 'max_depth'],
                        subsample = bestParams[1, 'subsample'],
                        colsample_bytree = bestParams[1,'colsample_bytree'],
                        data = trainm,
                        watchlist = watchlist,
                        verbose = 0)
  xgb.pred <- predict(xgbModel, test.data)
  xgb.pred <- as.data.frame(xgb.pred)
  pred.class <- ifelse(test = xgb.pred > 0.5, yes = TRUE, no = FALSE)
  modelAcc <- Accuracy(y_pred = pred.class, y_true = test.label)
  pred_obj = prediction(xgb.pred, test.label, label.ordering = c("FALSE", "TRUE"))
  modelAUC <- performance(prediction.obj = pred_obj, measure = "auc")@y.values
  modelF1score <- F1_Score(y_true = test.label, y_pred = pred.class, positive = 'TRUE')
  # end time
  end.time <- Sys.time()
  # time passed
  timePassed <- end.time - start.time
  Model <- list(ACCURACY = modelAcc, AUC = modelAUC, F1SCORE = modelF1score,  TIME = timePassed, PARAMETERS = bestParams[,c(1:5)], MODEL = xgbModel)
  return(Model)
}
