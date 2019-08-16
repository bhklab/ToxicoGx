# extract from phenodata file
phenoData <- readRDS("rds/phenoData.rds")
#subset necessary columns
cell <- subset(phenoData,select=c(samplename, organ_id, material_id, species, test_type, batchid))
# rename columns
names(cell)<-c("cellid","tissueid","materialid", "species","testType","batchid")
#entries in tissueid must match entries in tissue curation
cell$tissueid<-"liver"
rownames(cell) <- cell$cellid

#save object into rds/
saveRDS(cell,file="rds/cell.rds")

