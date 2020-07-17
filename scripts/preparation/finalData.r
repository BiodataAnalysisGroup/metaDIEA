#combine pipelines' data with simulated
pipData <- read.csv(file = "./combinedData.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
simData <- read.csv(file = "./simdata.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

#merge results
pipData <- merge(x = pipData, y = simData, by = intersect(names(pipData), names(simData)), all.x = TRUE, all.y = FALSE)

#write to csv
write.csv(x = pipData, file = "./FinalData.csv", row.names = FALSE)
