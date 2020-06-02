hisat_prep <- function(hisat2_data, genes){
  print("Hisat2 - DESeq2....")
  names(hisat2_data)[7] <- "transcript_id"
  #merge with respect to transcript_id
  hisat2_data <- merge(x = hisat2_data, y = genes, by = intersect(names(hisat2_data), names(genes)), all.x = TRUE, all.y = FALSE)
  hisat2_data <- as.data.frame(lapply(hisat2_data, unlist))
  hisat2_data <- hisat2_data[, c(8,1:7)]
  names(hisat2_data)[1] <- "Id"
  # convert character variables to factor
  hisat2_data[sapply(hisat2_data, is.character)] <- lapply(hisat2_data[sapply(hisat2_data, is.character)],as.factor)
  # aggregate factors
  d_fac <- hisat2_data[, sapply(X = hisat2_data, FUN = is.factor)]
  d_fac <- aggregate(x = d_fac, by = list(hisat2_data$Id), FUN = Mode)[, -1]
  # aggregate numerics
  d_num <- hisat2_data[, sapply(X = hisat2_data, FUN = is.numeric)]
  d_num <- aggregate(x = d_num, by = list(hisat2_data$Id), FUN = mean)
  names(d_num)[1] <- "Id"
  # merge aggregated results
  hisat2_data <- merge(x = d_num, y = d_fac, intersect(names(d_num), names(d_fac)), all = TRUE)
  return(hisat2_data)
}