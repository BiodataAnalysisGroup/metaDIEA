library(GenomicFeatures)
library(pracma)
#read kallisto-sleuth data
kallistoData <- read.csv(file = "C:/Users/Konstantinos/Desktop/Thesis/Pipeline1/simulated_data.transcript_results.sleuth.csv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
names(kallistoData)[15] <- "transcript_id"
# read genes.gtf file
txdb <- makeTxDbFromGFF(file = "./genes.gtf", format = "gtf")
genes <- as.data.frame(transcripts(txdb, columns = c("GENEID", "TXNAME"))@elementMetadata)
colnames(genes) <- c("gene_id", "transcript_id")
# merge data
kallistoData <- merge(x = kallistoData, y = genes, by = intersect(names(kallistoData), names(genes)), all.x = TRUE)
kallistoData <- as.data.frame(lapply(kallistoData, unlist))
# convert numeric target_id to factor
for(i in 1:ncol(kallistoData)){
  if(grepl(pattern = '*_id', x = colnames(kallistoData)[i])){
    kallistoData[,i] <- as.factor(kallistoData[,i])
  }
}
names(kallistoData)[16] <- "Id"
# convert character to factor
kallistoData[sapply(kallistoData, is.character)] <- lapply(kallistoData[sapply(kallistoData, is.character)],as.factor)
# aggregate factors
d_fac <- kallistoData[, sapply(X = kallistoData, FUN = is.factor)]
d_fac <- aggregate(x = d_fac, by = list(kallistoData$Id), FUN = Mode)[, c(2:4)]
# aggregate numerics
d_num <- kallistoData[, sapply(X = kallistoData, FUN = is.numeric)]
d_num <- aggregate(x = d_num, by = list(kallistoData$Id), FUN = mean)
names(d_num)[1] <- "Id"
# merge aggregated data
kallistoData <- merge(x = d_num, y = d_fac, by = intersect(names(d_num), names(d_fac)), all = TRUE)
remove(genes, d_fac, d_num)
# sava Kallisto-Sleuth data to file
write.csv(x = kallistoData, file = "C:/Users/Konstantinos/Desktop/Thesis/MLData/kallistoData.csv", row.names = FALSE)
