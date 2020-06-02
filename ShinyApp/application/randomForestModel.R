randomForestModel <- function(finalData, nModels){
  print("Random Forest....")
  # start time
  start.time <- Sys.time()
  # Split data to training(0.7) and testing(0.3)
  rndsplt <- sample(x = 2, size = nrow(finalData), replace = TRUE, prob = c(0.7,0.3))
  trainingData <- finalData[rndsplt == 1, ]
  testingData <- finalData[rndsplt == 2, ]
  # Search grid
  if(length(trainingData) <= 5){
    mtry <- 3
  }else{
    mtry <- ceiling(sqrt(length(trainingData))) + (-1:3)
  }
  rfGrid <- expand.grid(ntreeArg = (1:5)*500,
                        replaceAge = c(TRUE, FALSE),
                        mtryArg = mtry,
                        sampsizeArg  = c(nrow(trainingData), round(nrow(trainingData)*.9), round(nrow(trainingData)*.8), round(nrow(trainingData)*.75)),
                        nodesize = c(1, 3, 5, 7)
  )
  colnames(rfGrid) <- c("ntree", "replace", "mtry", "sampsize", "nodesize")
  rfGrid <- rfGrid[sample(nrow(rfGrid)),]
  numModels <- nModels; # specified by the user
  currParams <- as.data.frame(matrix(data = NA, nrow = 1, ncol = ncol(rfGrid) + 3))
  colnames(currParams) <- c(colnames(rfGrid), "Accuracy", "AUC", "F1score")
  withProgress(message = paste0("Generating ", numModels, " models"), detail = "May take a while", value = 0, {
    # Generate nModels
    for (i in 1:numModels){
      currParams[i,] <- rfGrid[sample(nrow(rfGrid), size = 1),] # sample random combination and train the forest
      incProgress(1/numModels, detail = paste0("Model ", i))
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
      # confusionMatrix(data = pred, reference = testingData$DE, positive = 'TRUE')
      # Accuracy
      currParams[i,'Accuracy'] <- Accuracy(y_pred = pred, y_true = testingData[,ncol(testingData)])
      pred_obj = prediction(predprob[,2], testingData[,ncol(testingData)], label.ordering = c("FALSE", "TRUE"))
      # Area Under the Curve
      currParams[i,'AUC'] <- performance(prediction.obj = pred_obj, measure = "auc")@y.values
      # F1 measure
      currParams[i,'F1score'] <- F1_Score(y_true = testingData$DE, y_pred = pred, positive = 'TRUE')
      
    }
  })
  # Best Model Parameters
  bestModel <- currParams[which(currParams$Accuracy == max(currParams$Accuracy)),]
  bestModel <- bestModel[which(bestModel$F1score == max(bestModel$F1score)),]
  bestModel <- bestModel[which(bestModel$AUC == max(bestModel$AUC)),]
  bestParams <- bestModel[1,c(1:5)]
  # Final Model
  rfModel <- randomForest(x = trainingData[,-ncol(trainingData)], y = trainingData$DE,
                          ntree = bestParams[1,'ntree'],
                          replace = bestParams[1,'replace'],
                          mtry = bestParams[1,'mtry'],
                          sampsize = bestParams[1,'sampsize'],
                          nodesize = bestParams[1,'nodesize'])
  pred <- predict(object = rfModel, testingData[,-ncol(testingData)], type = 'response')
  predprob <- predict(object = rfModel, testingData[,-ncol(testingData)], type = 'prob')
  modelAcc <- Accuracy(y_pred = pred, y_true = testingData[,ncol(testingData)])
  modelAUC <- performance(prediction.obj = pred_obj, measure = "auc")@y.values
  modelF1score <- F1_Score(y_true = testingData$DE, y_pred = pred, positive = 'TRUE')
  # end time
  end.time <- Sys.time()
  # time passed
  timePassed <- end.time - start.time
  Model <- list(ACCURACY = modelAcc, AUC = modelAUC, F1SCORE = modelF1score,  TIME = timePassed, PARAMETERS = bestParams[,c(1:5)], MODEL = rfModel)
  return(Model)
}