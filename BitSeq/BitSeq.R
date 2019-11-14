#!/usr/bin/env Rscript
library(BitSeq)
library(doParallel)

# source("https://bioconductor.org/biocLite.R")
# biocLite("BitSeq")
genome="./mm10/transcript.mm10.fa"

Bams=paste("sample_0",seq(1,6),".bam",sep="")

#run in parallel using 6 processes
sn=Bams[1]

cl <- makeCluster(length(Bams))
registerDoParallel(cl)
clusterEvalQ(cl, {library(BitSeq)})

res=foreach(sn=Bams) %dopar% 
{
  parseAlignment( sn,
                  outFile = paste(sn,".prob",sep=""),
                  trSeqFile = genome,
                  trInfoFile = paste(sn,".data.tr",sep=""),
                  verbose = TRUE )
  
  
}

#fix the trInfoFile, by adding the transcript names
tri.orig=tri.load(paste(genome,".tlst",sep=""))

for (sn in Bams)
{
tri=tri.load(paste(sn,".data.tr",sep=""))
tri$transcript=tri.orig$transcript
tri.save(tri,paste(sn,".data.tr",sep=""))
}

res=foreach(sn=Bams) %dopar% 
{
  
  estimateExpression(paste(sn,".prob",sep=""), outFile = paste(sn,".output",sep=""),
                     outputType = "RPKM", trInfoFile = paste(sn,".data.tr",sep=""),
                     MCMC_burnIn = 1000, MCMC_samplesN = 1000, MCMC_samplesSave = 1000,
                     MCMC_chainsN = 2)
  
}


allConditions = list(c("sample_01.bam.output.rpkm","sample_02.bam.output.rpkm","sample_03.bam.output.rpkm"),
                     c("sample_04.bam.output.rpkm","sample_05.bam.output.rpkm","sample_06.bam.output.rpkm"))
getMeanVariance(allConditions, outFile = "data.means", log = TRUE )

estimateHyperPar( outFile = "data.par",
                  conditions = allConditions,
                  meanFile = "data.means",
                  verbose = TRUE )

estimateDE(allConditions, outFile = "data", parFile = "data.par",samples = TRUE )
stopCluster(cl)