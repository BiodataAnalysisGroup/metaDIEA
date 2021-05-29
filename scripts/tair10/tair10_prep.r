library(GenomicFeatures)
library(pracma)
datatair10 <- read.table(file = './tair10AllMethods.tsv', sep = '\t', header = TRUE)
# read gened.gtf file
txdb <- makeTxDbFromGFF(file = "./tair10_genes.gtf", format = "gtf")
genes <- as.data.frame(transcripts(txdb, columns = c("GENEID", "TXNAME"))@elementMetadata)
colnames(genes) <- c("gene_id", "transcript_id")
# rename transctipt col
colnames(datatair10)[1] <- "transcript_id"
# merge for genes
datatair10 <- merge(x = datatair10, y = genes, by = "transcript_id", all.x = TRUE, all.y = FALSE)
datatair10 <- datatair10[!is.na(datatair10$gene_id),]
datatair10 <- datatair10[,c(length(datatair10),3:length(datatair10)-1)]
datatair10 <- datatair10[!duplicated(datatair10$gene_id),]
datatair10 <- as.data.frame(lapply(datatair10, unlist))
colnames(datatair10)[1] <- "Id"
datatair10$class <- ifelse(test = (datatair10$class == 0),yes = FALSE,no = TRUE)
#write to csv
write.csv(x = datatair10, file = "./tair10Final.csv", row.names = FALSE)