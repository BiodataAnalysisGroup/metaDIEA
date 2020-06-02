rfPrediction <- function(indata, model){
  predData <- indata
  # replace possible NAs with 0
  indata[is.na(indata)] <- 0
  # convert character columns to factors
  indata[sapply(indata, is.character)] <- lapply(indata[sapply(indata, is.character)],as.factor)
  # remove factors with many levels
  indata <- indata[, which(lapply(X = indata, FUN = nlevels) < 50)]
  # predict classes using the model
  rf.pred <- predict(object = model, indata)
  predData[ncol(predData)+1] <- as.data.frame(rf.pred)
  names(predData)[ncol(predData)] <- "DE"
  return(predData)
}