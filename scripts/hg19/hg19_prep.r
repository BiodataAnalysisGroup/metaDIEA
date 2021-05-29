datahg19 <- read.table(file = './hg19AllMethods.tsv', sep = '\t', header = TRUE)
# read gened.gtf file
genes <- read.table('./hg19_genes.gtf', header = FALSE, sep = '\t')[,c(2,13)]
colnames(genes) <- c("transcript_id", "gene_id")
# rename transctipt col
colnames(datahg19)[1] <- "transcript_id"
# merge for genes
datahg19 <- merge(x = datahg19, y = genes, by = "transcript_id", all.x = TRUE, all.y = FALSE)
datahg19 <- datahg19[!is.na(datahg19$gene_id),]
datahg19 <- datahg19[,c(length(datahg19),3:length(datahg19)-1)]
datahg19 <- datahg19[!duplicated(datahg19$gene_id),]
datahg19 <- as.data.frame(lapply(datahg19, unlist))
colnames(datahg19)[1] <- "Id"
#simData <- read.csv(file = "../simdata.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
#dataDM6 <- merge(x = dataDM6, y = simData, by = "Id", all.x = TRUE, all.y = FALSE)
datahg19$class <- ifelse(test = (datahg19$class == 0),yes = FALSE,no = TRUE)
#write to csv
write.csv(x = datahg19, file = "./hg19Final.csv", row.names = FALSE)