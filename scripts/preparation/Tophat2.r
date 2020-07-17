library(pracma) #used for Mode function
# tophat2-cufflinks Pipeline
cufflinksIsoformsDiff <- read.table(file = "C:/Users/Konstantinos/Desktop/Thesis/Pipeline2/diff_out/isoform_exp.diff", sep = "", header = TRUE, stringsAsFactors = FALSE)
cufflinksIsoformsFpkm <- read.table(file = "C:/Users/Konstantinos/Desktop/Thesis/Pipeline2/diff_out/isoforms.fpkm_tracking", sep = "", header = TRUE)
cufflinksAllData <- cbind(cufflinksIsoformsDiff, cufflinksIsoformsFpkm[, -c(4,7)])
# keep only significant DE genes
cuffData <- cufflinksAllData[which(cufflinksAllData$significant == "yes"),]
cuffData <- cuffData[, c(3,1:2,4:ncol(cuffData))]
colnames(cuffData)[1] <- "Id"
cuffData[sapply(cuffData, is.character)] <- lapply(cuffData[sapply(cuffData, is.character)],as.factor)

# split to numeric and factor data
d_fac <- cuffData[, sapply(X = cuffData, FUN = is.factor)]
d_fac <- aggregate(x = d_fac, by = list(cuffData$Id), FUN = Mode)[, c(2:19)]

d_num <- cuffData[, sapply(X = cuffData, FUN = is.numeric)]
d_num <- aggregate(x = d_num, by = list(cuffData$Id), FUN = mean)
names(d_num)[1] <- "Id"

# merge aggregated data
cuffData <- merge(x = d_num, y = d_fac, by = intersect(names(d_num), names(d_fac)), all = TRUE)

# save Tophat2-Cufflinks data
write.csv(cuffData, "C:/Users/Konstantinos/Desktop/Thesis/MLData/CufflinksData.csv", row.names = FALSE)
