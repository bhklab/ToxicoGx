# extract from phenodata file
phenoData <- readRDS("data/phenoData.rds")
#subset necessary columns
drug <- unique(subset(phenoData,select=c(drugid, drugid_abbr, drugid_no)))
#change rownames
rownames(drug) <- drug$drugid

#save object into data/
saveRDS(drug,file="data/drug.rds")