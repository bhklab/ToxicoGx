# extract from phenodata file
phenoData <- readRDS("data/phenoData.rds")
#subset necessary columns
drug <- unique(subset(phenoData,select=c(drugid, drugid_abbr, drugid_no)))
#change rownames
rownames(drug) <- drug$drugid
#subset out missing drugs
drug <- subset(drug,drugid != "rosiglitazone maleate")
drug <- subset(drug,drugid != "propranolol")

#save object into data/
saveRDS(drug,file="data/drug.rds")