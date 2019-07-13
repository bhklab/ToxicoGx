#extract from phenodata file
phenoData <- readRDS("data/phenoData.rds")
#subset necessary columns
curationTissue <- subset(phenoData,select=c(organ_id))
#change column entries (need to be lowercase to match cell object)
curationTissue$organ_id <- "liver"
#add extra column
curationTissue$tggates.tissueid <- "liver"
#rename first column
names(curationTissue)[1] <- "unique.tissueid"

#save object into data/
saveRDS(curationTissue,file="data/curationTissue.rds")
