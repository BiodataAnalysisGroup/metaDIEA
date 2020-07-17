library(GenomicFeatures)
library(pracma)
# read rsemData
rsemData <- read.table(file = "C:/Users/Konstantinos/Desktop/Thesis/Pipeline5/IsoMat.de.txt", fill = TRUE , sep = "\t", header = TRUE, col.names = c("transcript_id", "PPEE", "PPDE", "PostFC", "RealFC", "C1Mean", "C2Mean"), stringsAsFactors = FALSE)
# read gened.gtf file
txdb <- makeTxDbFromGFF(file = "./genes.gtf", format = "gtf")
genes <- as.data.frame(transcripts(txdb, columns = c("GENEID", "TXNAME"))@elementMetadata)
colnames(genes) <- c("gene_id", "transcript_id")
# merge data
rsemData <- merge(x = rsemData, y = genes, by = "transcript_id", all.x = TRUE, all.y = FALSE)[, c(8,1:7)]
colnames(rsemData)[1] <- "Id"
rsemData <- as.data.frame(lapply(rsemData, unlist))
# convert character variables to factor
rsemData[sapply(rsemData, is.character)] <- lapply(rsemData[sapply(rsemData, is.character)],as.factor)
# aggregate factors
d_fac <- rsemData[, sapply(X = rsemData, FUN = is.factor)]
d_fac <- aggregate(x = d_fac, by = list(rsemData$Id), FUN = Mode)[, c(2,3)]
# aggregate numerics
d_num <- rsemData[, sapply(X = rsemData, FUN = is.numeric)]
d_num <- aggregate(x = d_num, by = list(rsemData$Id), FUN = mean)
names(d_num)[1] <- "Id"
# merge aggregated data
rsemData <- merge(x = d_num, y = d_fac, by = intersect(names(d_num), names(d_fac)), all = TRUE)
# save RSEM data to file
write.csv(rsemData, "C:/Users/Konstantinos/Desktop/Thesis/MLData/RsemData.csv", row.names = FALSE)
