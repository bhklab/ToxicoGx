#extract from phenodata file
phenoData <- readRDS("rds/phenoData.rds")
#subset necessary columns
curationTissue <- subset(phenoData,select=c(organ_id))
#change column entries (need to be lowercase to match cell object)
curationTissue$organ_id <- "liver"
#add extra column
curationTissue$tggates.tissueid <- "liver"
#rename first column
names(curationTissue)[1] <- "unique.tissueid"
rownames(curationTissue) <- curationCell$unique.cellid

#save object into rds/
saveRDS(curationTissue,file="rds/curationTissue.rds")

