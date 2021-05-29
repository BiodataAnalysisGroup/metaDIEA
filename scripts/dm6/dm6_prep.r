library(GenomicFeatures)
library(pracma)
dataDM6 <- read.table(file = './dm6AllMethods.tsv', sep = '\t', header = TRUE)
# read gened.gtf file
txdb <- makeTxDbFromGFF(file = "./dm6_genes.gtf", format = "gtf")
genes <- as.data.frame(transcripts(txdb, columns = c("GENEID", "TXNAME"))@elementMetadata)
colnames(genes) <- c("gene_id", "transcript_id")
# rename transctipt col
colnames(dataDM6)[1] <- "transcript_id"
# merge for genes
dataDM6 <- merge(x = dataDM6, y = genes, by = "transcript_id", all.x = TRUE, all.y = FALSE)
dataDM6 <- dataDM6[!is.na(dataDM6$gene_id),]
dataDM6 <- dataDM6[,c(length(dataDM6),3:length(dataDM6)-1)]
dataDM6 <- dataDM6[!duplicated(dataDM6$gene_id),]
dataDM6 <- as.data.frame(lapply(dataDM6, unlist))
colnames(dataDM6)[1] <- "Id"
dataDM6$class <- ifelse(test = (dataDM6$class == 0),yes = FALSE,no = TRUE)
#write to csv
write.csv(x = dataDM6, file = "./dm6Final.csv", row.names = FALSE)