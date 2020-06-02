xgbPrediction <- function(indata, model){
  predData <- indata
  # replace possible NAs with 0
  indata[is.na(indata)] <- 0
  # if(names(indata)[1] == "target_id"){
  #   indata <- indata[,-1]  
  # }
  if("target_id" %in% names(indata)){
    indata <- indata[, -which(names(indata) == "target_id")]
  }
  # # convert character columns to factors
  indata[sapply(indata, is.character)] <- lapply(indata[sapply(indata, is.character)],as.factor)
  # remove factors with many levels
  indata <- indata[, which(lapply(X = indata, FUN = nlevels) < 50)]
  # numeric encoding for xgboost
  test.data = as.matrix(sapply(indata, as.numeric))
  # predict classes using the model
  xgb.pred <- predict(object = model, test.data)
  xgb.pred <- ifelse(test = xgb.pred > 0.5, yes = TRUE, no = FALSE)
  predData[ncol(predData)+1] <- as.data.frame(xgb.pred)
  names(predData)[ncol(predData)] <- "DE"
  return(predData)
}