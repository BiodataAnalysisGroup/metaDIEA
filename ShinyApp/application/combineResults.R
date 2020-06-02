compineResults <- function(sleuth_data, tophat_data, hisat2_data, bitseq_data, rsem_data, ebseq_data, simData){
  #merge rsemData and hisatData
  rsemHisatData <- merge(hisat2_data, rsem_data, by <- "Id", all = TRUE)
  #merge rsemHisatData and cuffData
  rsemHisatCuffData <- merge(rsemHisatData, tophat_data, by <- "Id", all = TRUE)
  #merge with ebseqData
  rsemHisatCuffEbseqData <- merge(rsemHisatCuffData, ebseq_data, by <- "Id", all = TRUE)
  #merge with kallistoData
  rsemHisatCuffEbseqKallistoData <- merge(rsemHisatCuffEbseqData, sleuth_data, by <- "Id", all = TRUE)
  #merge with bitseqData 
  finalData <- merge(rsemHisatCuffEbseqKallistoData, bitseq_data, by <- "Id", all = TRUE)  
  finalData <- merge(x = finalData, y = simData, by = intersect(names(finalData), names(simData)), all.x = TRUE, all.y = FALSE)
  # drop columns that were not selected
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
  return(finalData)
}








