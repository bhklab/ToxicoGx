#extract from phenodata file
phenoData <- readRDS("data/phenoData.rds")
#subset necessary columns
curationCell <- subset(phenoData,select=c(samplename))
#add new column
curationCell$tggates.cellid <- curationCell$samplename
#change first column name
names(curationCell)[1] <- "unique.cellid"

#save object into data/
saveRDS(curationCell,file="data/curationCell.rds")
