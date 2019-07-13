#extract from phenodata file
phenoData <- readRDS("data/phenoData.rds")
#read from drug mapping file
drug_csv <- read.csv("data/curationDrug_mapping.csv")
#subset necessary columns
curationDrug <- subset(drug_csv, select=-c(X, Type.of.Mapping))
#rename columns
names(curationDrug) <- c("unique.drugid", "tggates.drugid")
#change rownames
rownames(curationDrug) <- curationDrug$unique.drugid
#subset out missing drugs
curationDrug <- subset(curationDrug,unique.drugid != "rosiglitazone maleate")
curationDrug <- subset(curationDrug,unique.drugid != "propranolol")

#save object into data/
saveRDS(curationDrug,file="data/curationDrug.rds")