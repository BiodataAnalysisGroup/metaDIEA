simdata_prep <- function(simData, genes){
  #keep only the transcript name
  lst <- strsplit(simData[,1]," ")
  simData[1] <- sapply(lst, '[[', 2)
  names(simData)[1] <- "transcript_id"
  #merge simData to get gene names
  simData <- merge(x = simData, y = genes, by = "transcript_id", all.x = TRUE, all.y = FALSE)[, 2:6]
  simData <- as.data.frame(lapply(simData, unlist))
  simData <- simData[,c(5,1:4)]
  #simData <- simData[order(simData$gene_name),]
  #aggregate DE status
  simData[6] <- ifelse(test = (simData$DEstatus.1 == TRUE | simData$DEstatus.2 == TRUE), yes = TRUE, no = FALSE)
  simData <- simData[, c(1:3,6)]
  names(simData)[4] <- "DE"
  names(simData)[1] <- "Id"
  #convert DE status to numeric, in order to aggregate duplicates
  simData[5] <- ifelse(test = (simData$DE == TRUE), yes = 1,no = 0)
  #aggregate duplicates with respect to DE status
  simData <-aggregate(formula = V5 ~ Id, data = simData, FUN = sum)
  names(simData)[2] <- "DE"
  simData[2] <- ifelse(test = simData$DE == 0, yes = FALSE, no = TRUE)  
  return(simData)
}

