#extract from phenodata file
phenoData <- readRDS("rds/phenoData.rds")
#subset necessary columns
curationCell <- subset(phenoData,select=c(samplename))
#add new column
curationCell$tggates.cellid <- curationCell$samplename
#change first column name
names(curationCell)[1] <- "unique.cellid"
rownames(curationCell) <- curationCell$unique.cellid

#save object into rds/
saveRDS(curationCell,file="rds/curationCell.rds")

