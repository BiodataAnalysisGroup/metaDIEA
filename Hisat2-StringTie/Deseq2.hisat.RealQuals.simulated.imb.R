#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

setwd(args[1])
#library(Rsubread)
library("DESeq2")
library("BiocParallel")
library(IRanges)
register(MulticoreParam(1))

countData <- as.matrix(read.csv("transcript_count_matrix.csv", row.names="transcript_id"))
colData <- read.csv("samples.csv", sep=",", row.names=1)
all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData, design = ~ type)
dds <- DESeq(dds,parallel = T)

res <- results(dds,parallel = T)

res.filt=na.omit(res)
a=as.data.frame(res.filt[res.filt$padj<0.05,])
a$trx=rownames(a)


write.csv(a, "simulated_data.transcript_results.Deseq2.csv",
          row.names=FALSE)
