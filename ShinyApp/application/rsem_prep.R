rsem_prep <- function(rsem_data, genes){
  print("RSEM...")
  rsem_data <- merge(x = rsem_data, y = genes, by = intersect(names(rsem_data), names(genes)), all.x = TRUE, all.y = FALSE)
  rsem_data <- as.data.frame(lapply(rsem_data, unlist))
  rsem_data <- rsem_data[, c(8,1:7)]
  colnames(rsem_data)[1] <- "Id"
  # convert character variables to factor
  rsem_data[sapply(rsem_data, is.character)] <- lapply(rsem_data[sapply(rsem_data, is.character)],as.factor)
  # aggregate factors
  d_fac <- rsem_data[, sapply(X = rsem_data, FUN = is.factor)]
  d_fac <- aggregate(x = d_fac, by = list(rsem_data$Id), FUN = Mode)[, c(2,3)]
  # aggregate numerics
  d_num <- rsem_data[, sapply(X = rsem_data, FUN = is.numeric)]
  d_num <- aggregate(x = d_num, by = list(rsem_data$Id), FUN = mean)
  names(d_num)[1] <- "Id"
  # merge aggregated results
  rsem_data <- merge(x = d_num, y = d_fac, by = intersect(names(d_num), names(d_fac)), all = TRUE)
  return(rsem_data)
}