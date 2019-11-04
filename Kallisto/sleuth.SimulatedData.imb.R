#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#simulated.imb
setwd(args[1])

library("sleuth")

trx=read.table("/app/trx.tsv")
colnames(trx)=c("target_id","trx")

sample_id <- list.files(path = "/app/", pattern="sample_")
  kal_dirs=paste("/app",sample_id,sep="/")

s2c=as.data.frame(matrix(ncol = 2,nrow = 6))
colnames(s2c)=c("sample","condition")
s2c[,1]=sample_id
s2c[,2]=c(rep("treated",3),rep("untreated",3))

s2c <- dplyr::mutate(s2c, path = kal_dirs)

print(s2c)

so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
so <- sleuth_wt(so, which_beta = 'conditionuntreated', which_model = 'full')

res_lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
res_wt <- sleuth_results(so, 'conditionuntreated')
res <- merge(res_lrt, res_wt[, c('target_id', 'b', 'se_b', 'mean_obs')], on = 'target_id', sort = FALSE)
sleuth_significant <- dplyr::filter(res, qval <= 0.05)

sleuth_significant=merge(sleuth_significant,trx, by="target_id")
head(sleuth_significant, 20)
write.table(sleuth_significant,quote = F,"/app/simulated_data.transcript_results.sleuth.csv",sep="\t",row.names = F)
