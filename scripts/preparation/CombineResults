#read HisatData
hisatData <- read.csv(file = "./Hisat2Data.csv", sep = ",", header = TRUE)

#read RsemData
rsemData <- read.csv(file = "./RsemData.csv", sep = ",", header = TRUE)

#read CufflinksData
cuffData <- read.csv(file = "./CufflinksData.csv", sep = ",", header = TRUE)

#read EBSeqData
ebseqData <- read.csv(file = "./EBSeqData.csv", sep = ",", header = TRUE)

#read KallistoData
kallistoData <- read.csv(file = "./kallistoData.csv", sep = ",", header = TRUE)

#read BitSeqData
bitseqData <- read.csv(file = "./bitseqData.csv", sep = ",", header = TRUE)

#merge rsemData and hisatData
rsemHisatData <- merge(hisatData, rsemData, by <- "Id", all = TRUE)[, c(1:14)]
names(rsemHisatData)[8] <- "transcript_id"
#merge rsemHisatData and cuffData
rsemHisatCuffData <- merge(rsemHisatData, cuffData, by <- "Id", all = TRUE)
#merge with ebseqData
rsemHisatCuffEbseqData <- merge(rsemHisatCuffData, ebseqData, by <- "Id", all = TRUE)[,-45]
names(rsemHisatCuffEbseqData)[8] <- "transcript_id"
names(rsemHisatCuffEbseqData)[18] <- "test_stat_c"
#merge with kallistoData
rsemHisatCuffEbseqKallistoData <- merge(rsemHisatCuffEbseqData, kallistoData, by <- "Id", all = TRUE)[,-59]
names(rsemHisatCuffEbseqKallistoData)[8] <- "transcript_id"
#merge with bitseqData 
combinedData <- merge(rsemHisatCuffEbseqKallistoData, bitseqData, by <- "Id", all = TRUE)[,-60]
names(combinedData)[8] <- "transcript_id"
#write data to file
write.csv(combinedData, "C:/Users/Konstantinos/Desktop/Thesis/MLData/combinedData.csv", row.names = FALSE)
