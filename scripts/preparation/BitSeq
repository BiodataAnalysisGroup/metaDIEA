# read bitseq results
bitseqData <- read.table(file = "C:/Users/Konstantinos/Desktop/Thesis/Pipeline4/data.pplr", sep = " ", header = FALSE, stringsAsFactors = FALSE)
colnames(bitseqData) <- c("PPLR", "Log2FC", "ConfLow", "ConfHigh", "LogMean", "MeanExpr")
# read transcripts
transcripts <- read.table(file = "C:/Users/Konstantinos/Desktop/Thesis/Pipeline4/sample_01.bam.data.tr", sep = " ", header = FALSE, stringsAsFactors = FALSE)[2]
colnames(transcripts) <- "transcript_id"
# add trancripts to results
bitseqData <- cbind(bitseqData, transcripts)
# bitseqData <- bitseqData[, c(7,1:6)]
# include only up or down-regulated transcripts
bitseqData <- bitseqData[which((bitseqData$PPLR < 0.05) | (bitseqData$PPLR > 0.95)), ]
remove(transcripts)
# read genes.gtf file
library(GenomicFeatures)
txdb <- makeTxDbFromGFF(file = "./genes.gtf", format = "gtf")
genes <- as.data.frame(transcripts(txdb, columns = c("GENEID", "TXNAME"))@elementMetadata)
colnames(genes) <- c("gene_id", "transcript_id")
# merge with genes
bitseqData <- merge(x = bitseqData, y = genes, by = intersect(names(bitseqData), names(genes)), all.x = TRUE)
bitseqData <- bitseqData[, c(8,1:7)]
names(bitseqData)[1] <- "Id"
# # convert character variables to factor
# bitseqData[sapply(bitseqData, is.character)] <- lapply(bitseqData[sapply(bitseqData, is.character)],as.factor)
# # aggregate factors
# d_fac <- bitseqData[, sapply(X = bitseqData, FUN = is.factor)]
# library(pracma)
# d_fac <- aggregate(x = d_fac, by = list(bitseqData$Id), FUN = Mode)[, c(2,3)]
# # aggregate numerics
# d_num <- bitseqData[, sapply(X = bitseqData, FUN = is.numeric)]
# d_num <- aggregate(x = d_num, by = list(bitseqData$Id), FUN = mean)
# names(d_num)[1] <- "Id"
# merge aggregated results
# bitseqData <- merge(x = d_num, y = d_fac, by = intersect(names(d_num), names(d_fac)), all = TRUE)
# save BitSeq data
write.csv(x = bitseqData, file = "C:/Users/Konstantinos/Desktop/Thesis/MLData/bitseqData.csv", row.names = FALSE)
