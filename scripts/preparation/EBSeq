# read ebseq data
ebseqData1 <- read.table(file = "C:/Users/Konstantinos/Desktop/Thesis/Pipeline6/EBSeq.Simulated.All.Out.FC.txt", fill = TRUE, sep = "\t", header = TRUE, col.names = c("trx", "PostFC"), stringsAsFactors = FALSE)
ebseqData2 <- read.table(file = "C:/Users/Konstantinos/Desktop/Thesis/Pipeline6/EBSeq.Simulated.DE.Out.txt", fill = TRUE , sep = "\t", header = TRUE, col.names = c("trx", "PPEE", "PPDE"), stringsAsFactors = FALSE)
# read genes.gtf file
library(GenomicFeatures)
txdb <- makeTxDbFromGFF(file = "./genes.gtf", format = "gtf")
genes <- as.data.frame(transcripts(txdb, columns = c("GENEID", "TXNAME"))@elementMetadata)
colnames(genes) <- c("gene_id", "transcript_id")
names(ebseqData1)[1] <- "transcript_id"
names(ebseqData2)[1] <- "transcript_id"
# merge PostFC data with genes
ebseqData1 <- merge(ebseqData1, genes, by <- intersect(names(ebseqData1), names(genes)), all.x = TRUE) 
ebseqData1 <- as.data.frame(lapply(ebseqData1, unlist))
ebseqData1[sapply(ebseqData1, is.character)] <- lapply(ebseqData1[sapply(ebseqData1, is.character)],as.factor)
# aggregate factors
library(pracma)
d_fac <- ebseqData1[, sapply(X = ebseqData1, FUN = is.factor)]
d_fac <- aggregate(x = d_fac, by = list(ebseqData1$gene_id), FUN = Mode)[, c(2,3)]
# aggregate numerics
d_num <- ebseqData1[, sapply(X = ebseqData1, FUN = is.numeric)]
d_num <- aggregate(x = d_num, by = list(ebseqData1$gene_id), FUN = mean)
names(d_num)[1] <- "gene_id"
# merge data with aggregated PostFC
ebseqData1 <- merge(x = d_num, y = d_fac, by = intersect(names(d_num), names(d_fac)), all = TRUE)
# merge PPEE and PPDE data with genes
ebseqData2 <- merge(ebseqData2, genes, by <- intersect(names(ebseqData2), names(genes)), all.x = TRUE)
ebseqData2 <- as.data.frame(lapply(ebseqData2, unlist))
ebseqData2[sapply(ebseqData2, is.character)] <- lapply(ebseqData2[sapply(ebseqData2, is.character)],as.factor)
# aggregate factors
d_fac <- ebseqData2[, sapply(X = ebseqData2, FUN = is.factor)]
d_fac <- aggregate(x = d_fac, by = list(ebseqData2$gene_id), FUN = Mode)[, c(2,3)]
# aggregate numerics
d_num <- ebseqData2[, sapply(X = ebseqData2, FUN = is.numeric)]
d_num <- aggregate(x = d_num, by = list(ebseqData2$gene_id), FUN = mean)
names(d_num)[1] <- "gene_id"
ebseqData2 <- merge(x = d_num, y = d_fac, by = intersect(names(d_num), names(d_fac)), all = TRUE)
# merge the two dataframes
ebseqData <- merge(ebseqData2, ebseqData1, by <- "gene_id", all.x = TRUE)[,c(1:5)] 
names(ebseqData) <- c("Id", "PPEE_", "PPDE_", "transcript_id" , "PostFC.IsoEBOut..PostFC")
remove(ebseqData1, ebseqData2, d_fac, d_num)
#write to csv
write.csv(ebseqData, "C:/Users/Konstantinos/Desktop/Thesis/MLData/EBSeqData.csv", row.names = FALSE)
