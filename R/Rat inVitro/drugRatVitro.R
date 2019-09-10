# extract from phenodata file
phenoData <- readRDS("rds/phenoData.rds")
#subset necessary columns
drug <- unique(subset(phenoData,select=c(drugid, drugid_abbr, drugid_no)))
#change rownames
rownames(drug) <- drug$drugid
names(drug) <- c("drugid", 'drugid_abbr', 'drugid_no')
rownames(drug) <- drug$drugid

#save object into data/
saveRDS(drug,file="rds/drug.rds")

