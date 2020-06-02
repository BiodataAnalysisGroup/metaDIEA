sleuth_prep <- function(sleuth_data, genes){
  print("Kallisto - Sleuth...")
  names(sleuth_data)[15] <- "transcript_id"
  #merge with genes.gtf
  sleuth_data <- merge(x = sleuth_data, y = genes, by = intersect(names(sleuth_data), names(genes)), all.x = TRUE)
  sleuth_data <- as.data.frame(lapply(sleuth_data, unlist))
  # convert numeric Ids to factor
  for(i in 1:ncol(sleuth_data)){
    if(grepl(pattern = '*_id', x = colnames(sleuth_data)[i])){
      sleuth_data[,i] <- as.factor(sleuth_data[,i])
    }
  }
  names(sleuth_data)[16] <- "Id"
  #sleuth_data <- sleuth_data[, c(ncol(sleuth_data),1:ncol(sleuth_data)-1)]
  # aggregate factor variables
  d_fac <- sleuth_data[, sapply(X = sleuth_data, FUN = is.factor)]
  d_fac <- aggregate(x = d_fac, by = list(sleuth_data$Id), FUN = Mode)[, c(2:4)]
  # aggregate numeric variables
  d_num <- sleuth_data[, sapply(X = sleuth_data, FUN = is.numeric)]
  d_num <- aggregate(x = d_num, by = list(sleuth_data$Id), FUN = mean)
  names(d_num)[1] <- "Id"
  # merge agragated results
  sleuth_data <- merge(x = d_num, y = d_fac, by = intersect(names(d_num), names(d_fac)), all = TRUE)
  return(sleuth_data)
}