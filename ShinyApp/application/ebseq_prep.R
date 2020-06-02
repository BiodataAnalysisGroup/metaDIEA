ebseq_prep <- function(ebseq_data, ebseq_data2, genes){
  print("EBSeq...")
  ebseq_data <- merge(ebseq_data, ebseq_data2, by <- intersect(names(ebseq_data), names(ebseq_data2)), all = FALSE)
  names(ebseq_data)[1] <- "transcript_id"
  ebseq_data <- merge(x = ebseq_data, y = genes, by = intersect(names(ebseq_data), names(genes)), all.x = TRUE, all.y = FALSE)
  ebseq_data <- as.data.frame(lapply(ebseq_data, unlist))
  ebseq_data <- ebseq_data[, c(5,1:4)]
  names(ebseq_data)[1] <- "Id"
  names(ebseq_data)[3] <- "PostFC.IsoEBOut..PostFC"
  names(ebseq_data)[4] <- "PPEE_"
  names(ebseq_data)[5] <- "PPDE_"
  # convert character variables to factor
  ebseq_data[sapply(ebseq_data, is.character)] <- lapply(ebseq_data[sapply(ebseq_data, is.character)],as.factor)
  # aggregate factors
  d_fac <- ebseq_data[, sapply(X = ebseq_data, FUN = is.factor)]
  d_fac <- aggregate(x = d_fac, by = list(ebseq_data$Id), FUN = Mode)[, c(2,3)]
  # aggregate numerics
  d_num <- ebseq_data[, sapply(X = ebseq_data, FUN = is.numeric)]
  d_num <- aggregate(x = d_num, by = list(ebseq_data$Id), FUN = mean)
  names(d_num)[1] <- "Id"
  # merge aggregated results
  ebseq_data <- merge(x = d_num, y = d_fac, by = intersect(names(d_num), names(d_fac)), all = TRUE)
  return(ebseq_data)
}