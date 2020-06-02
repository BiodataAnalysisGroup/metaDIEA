cufflinks_prep <- function(cufflinksIsoformsDiff, cufflinksIsoformsFpkm){
  print("Tophat2 - Cufflinks...")
  #combine results from files to single dataframe
  tophat_data <- cbind(cufflinksIsoformsDiff, cufflinksIsoformsFpkm[, c(-1, -4, -5, -7)])
  tophat_data <- tophat_data[which(tophat_data$significant == "yes"), c(3,1:2,4:27)]
  colnames(tophat_data)[1] <- "Id"
  # convert character variables to factor
  tophat_data[sapply(tophat_data, is.character)] <- lapply(tophat_data[sapply(tophat_data, is.character)],as.factor)
  # aggregate factor variables
  d_fac <- tophat_data[, sapply(X = tophat_data, FUN = is.factor)]
  d_fac <- aggregate(x = d_fac, by = list(tophat_data$Id), FUN = Mode)[, c(2:15)]
  # aggregate numeric variables
  d_num <- tophat_data[, sapply(X = tophat_data, FUN = is.numeric)]
  d_num <- aggregate(x = d_num, by = list(tophat_data$Id), FUN = mean)
  names(d_num)[1] <- "Id"
  # merge aggregated results
  tophat_data <- merge(x = d_num, y = d_fac, by = intersect(names(d_num), names(d_fac)), all = TRUE)
  return(tophat_data)
}