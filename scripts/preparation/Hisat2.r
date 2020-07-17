library(GenomicFeatures)
library(pracma)
# read hisat2-Deseq2 data
hisat2Data <- read.table(file = "C:/Users/Konstantinos/Desktop/Thesis/Pipeline3/simulated_data.transcript_results.Deseq2.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
names(hisat2Data)[7] <- "transcript_id"
# read genes.gtf file
txdb <- makeTxDbFromGFF(file = "./genes.gtf", format = "gtf")
genes <- as.data.frame(transcripts(txdb, columns = c("GENEID", "TXNAME"))@elementMetadata)
colnames(genes) <- c("gene_id", "transcript_id")
#merge with respect to transcript_id
hisat2Data <- merge(x = hisat2Data, y = genes, by = intersect(names(hisat2Data), names(genes)), all.x = TRUE, all.y = FALSE)[, c(8,1:7)]
hisat2Data <- as.data.frame(lapply(hisat2Data, unlist))
names(hisat2Data)[1] <- "Id"
# convert character to factor so as to aggregate
hisat2Data[sapply(hisat2Data, is.character)] <- lapply(hisat2Data[sapply(hisat2Data, is.character)],as.factor)
# aggregate factors
d_fac <- hisat2Data[, sapply(X = hisat2Data, FUN = is.factor)]
d_fac <- aggregate(x = d_fac, by = list(hisat2Data$Id), FUN = Mode)[, -1]
# aggregate numerics
d_num <- hisat2Data[, sapply(X = hisat2Data, FUN = is.numeric)]
d_num <- aggregate(x = d_num, by = list(hisat2Data$Id), FUN = mean)
names(d_num)[1] <- "Id"
# merge aggregated data
hisat2Data <- merge(x = d_num, y = d_fac, intersect(names(d_num), names(d_fac)), all = TRUE)
# save Hisat2-DESeq2 data to file
write.csv(hisat2Data, "C:/Users/Konstantinos/Desktop/Thesis/MLData/Hisat2Data.csv", row.names = FALSE)
