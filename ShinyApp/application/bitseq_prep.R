bitseq_prep <- function(bitseq_data, bitseq_transcripts, genes){
  print("BitSeq...")
  colnames(bitseq_data) <- c("PPLR", "Log2FC", "ConfLow", "ConfHigh", "LogMean", "MeanExpr")
  colnames(bitseq_transcripts) <- "transcript_id"
  bitseq_data <- cbind(bitseq_data, bitseq_transcripts)
  #include only up or down-regulated transcripts
  bitseq_data <- bitseq_data[which((bitseq_data$PPLR < 0.2) | (bitseq_data$PPLR > 0.8)),]
  bitseq_data <- merge(x = bitseq_data, y = genes, by = intersect(names(bitseq_data), names(genes)), all.x = TRUE)
  bitseq_data <- as.data.frame(lapply(bitseq_data, unlist))
  bitseq_data <- bitseq_data[, c(8,1:7)]
  names(bitseq_data)[1] <- "Id"
  # convert character variables to factor
  bitseq_data[sapply(bitseq_data, is.character)] <- lapply(bitseq_data[sapply(bitseq_data, is.character)],as.factor)
  # aggregate factors
  d_fac <- bitseq_data[, sapply(X = bitseq_data, FUN = is.factor)]
  d_fac <- aggregate(x = d_fac, by = list(bitseq_data$Id), FUN = Mode)[, c(2,3)]
  # aggregate numerics
  d_num <- bitseq_data[, sapply(X = bitseq_data, FUN = is.numeric)]
  d_num <- aggregate(x = d_num, by = list(bitseq_data$Id), FUN = mean)
  names(d_num)[1] <- "Id"
  # merge aggregated results
  bitseq_data <- merge(x = d_num, y = d_fac, by = intersect(names(d_num), names(d_fac)), all = TRUE)
  return(bitseq_data)
}