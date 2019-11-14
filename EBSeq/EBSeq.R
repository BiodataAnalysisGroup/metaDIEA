#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(EBSeq)
library("biomaRt")
# source("https://bioconductor.org/biocLite.R")
# biocLite("EBSeq")

IsoMat=read.table("app/IsoMat.txt",as.is = T)
IsoMat=IsoMat[order(rownames(IsoMat)),]
IsoNames=rownames(IsoMat)
IsoMat=as.matrix(IsoMat)

Genes=read.table("app/Refseq2Gene.txt",as.is = T)
A=as.character(unlist(Genes[1]))
B=as.character(unlist(Genes[2]))
Genes2Trans=data.frame(A,B)

colnames(Genes2Trans)=c("Transcript","Gene")

IsosGeneNames=IsoNames
i=1
for (i in 1:length(IsosGeneNames))
{
  print(paste(i,length(IsosGeneNames),sep="/"))
  ToKeep=which(Genes2Trans[,1]==IsosGeneNames[i])
  if (length(ToKeep)>0)
    IsosGeneNames[i]=as.character(Genes2Trans[ToKeep[1],2])
}

  
Conditions=as.factor(c(rep("0",3),rep("1",3)))
IsoSizes=MedianNorm(IsoMat)

NgList=GetNg(IsoNames, IsosGeneNames)
IsoNgTrun=NgList$IsoformNgTrun


IsoEBOut=EBTest(Data=IsoMat, NgVector=IsoNgTrun,
                Conditions=Conditions,sizeFactors=IsoSizes, maxround=5)

 IsoEBDERes=GetDEResults(IsoEBOut, FDR=0.05)

head(IsoEBDERes$PPMat)
str(IsoEBDERes$Status)
IsoEBOut$PPMat
length(IsoEBDERes$DEfound)


write.table(IsoEBDERes$PPMat,file="app/EBSeq.Simulated.All.Out.txt",sep="\t",quote = F)
write.table(IsoEBDERes$PPMat[rownames(IsoEBDERes$PPMat)%in%IsoEBDERes$DEfound,],file="app/EBSeq.Simulated.DE.Out.txt",sep="\t",quote = F)


write.table(as.data.frame(PostFC(IsoEBOut)$PostFC),file="app/EBSeq.Simulated.All.Out.FC.txt",sep="\t",quote = F)

