# extract from phenodata file
phenoData <- readRDS("data/phenoData.rds")
#subset necessary columns
cell <- subset(phenoData,select=c(samplename, organ_id, material_id, species, test_type, batchid))
# rename columns
names(cell)<-c("cellid","tissueid","materialid", "species","testType","batchid")
#entries in tissueid must match entries in tissue curation
cell$tissueid<-"liver"

#save object into data/
saveRDS(cell,file="data/cell.rds")