#read from drug mapping file
drug_csv <- read.csv("data/TG-GATES Drug ID (RatInVitro)_mapped.csv",stringsAsFactors = FALSE)
#subset necessary columns
curationDrug <- subset(drug_csv, select=-c(X, Type.of.Mapping))
#rename columns
names(curationDrug) <- c("unique.drugid", "tggates.drugid")
#change rownames
rownames(curationDrug) <- curationDrug$unique.drugid

#save object into rds/
saveRDS(curationDrug,file="rds/curationDrug.rds")
