sensitivityInfo <- sensInfo
sensitivityProfiles <- sensProf
sensitivityRaw <- sensRaw

# rm(sensInfo)
# rm(sensProf)
# rm(sensRaw)

orderedcellid<-sensitivityInfo$cellid
names(orderedcellid)<-c("cellid")

########################################################################

## TISSUE CURATION ##
curationTissue<-orderedcellid
curationTissue<-as.data.frame(curationTissue)
curationTissue$unique.tissueid <- "liver"
curationTissue$tggates.tissueid <- "liver"

rownames(curationTissue)<-curationTissue$curationTissue
curationTissue<-subset(curationTissue,select=-c(curationTissue))

# write.csv(curationTissue, file="curationTissue.csv")
# ?data.frame

#########################################################################

## CELL CURATION ##
curationCell <- orderedcellid
curationCell <- as.data.frame(curationCell)
curationCell$unique.cellid <- curationCell$curationCell
curationCell$tggates.cellid <- curationCell$curationCell

rownames(curationCell) <- curationCell$curationCell
curationCell<-subset(curationCell, select=-c(curationCell))

#########################################################################

## DRUG CURATION ##
drug_csv <- read.csv("data/TG-GATES Drug ID (HumanInVitro)_mapped.csv")
View(curationDrug)
curationDrug <- subset(drug_csv, select=c(TG.GATES.Compound.HumanInvitro, Lab.s.unique.drugid))
names(curationDrug) <- c("unique.drugid", "tggates.drugid")
rownames(curationDrug) <- curationDrug$unique.drugid
rm(drug_csv)
#write.csv(curationDrug, file="curationDrugObject.csv")

#########################################################################

## CELL OBJECT ##
# extract from phenodata file
phenodata <- read.csv("data/phenowLot.csv",sep=",",header=TRUE, stringsAsFactors=TRUE)
# subset necessary columns
cell <- subset(phenodata, phenodata$DOSE_UNIT != "f?g/kg",select=c(X, ORGAN_ID, MATERIAL_ID, SPECIES, TEST_TYPE, CELL_LOT_TYPE))
cell <- arrange(cell, X)
# assign row names as barcode 
rownames(cell)<-cell$X
# rename column names
names(cell)<-c("cellid","tissueid","materialid", "species","testType","batchid")
cell$tissueid<-"liver"

#########################################################################

## DRUG OBJECT ##
drug_csv <- read.csv("data/drugObject.csv")
#head(drug_csv)

drug <- drug_csv
rownames(drug) <- drug_csv$drugid
drug <- drug[,-1]

#########################################################################

## subsetting from assayData ##
#DO NOT RUN THE FOLLOWING 2 LINES IN CASE I MESS UP
# oldAssay <- eset@assayData$exprs
# oldAssaySE <- eset@assayData$se.exprs
###################################################
missingCEL <- as.data.frame(subset(phenodata,phenodata$DOSE_UNIT != "f?g/kg", select=c(X)))
missingCEL <- arrange(missingCEL, X)
missingCEL$X <- paste0("00",missingCEL$X,".CEL")

eset <- readRDS('data/eset.rds')
storageMode(eset)<-"environment"
eset@assayData$exprs<-subset(eset@assayData$exprs,select=missingCEL$X)
storageMode(eset)<-"lockedEnvironment"
#pdata
path.sc=file.path("data")
sampleinfo <- read.table(file.path(path.sc, "TGGATEsfeatureInfo.txt"), sep="\t")
sampleinfo[sampleinfo == "" | sampleinfo == " "] <- NA
annot <- read.csv("data/phenoDataCOLNAMES_DONE.csv", stringsAsFactors=FALSE, check.names=FALSE, header=TRUE, row.names=1)
annot <- subset(annot, annot$concentration_old_units != "f?g/kg")
pData(eset) <- annot
#identify probe names
probes <- eset@assayData$exprs
#map ENSG probes in eset to feature data file of ENSG probes only
Featuredata <- read.csv("data/finalFeatureNAs.csv")
rownames(Featuredata) <- Featuredata$X
Featuredata <- subset(Featuredata, select = -c(X.1, X))
Featuredata$BEST <- NA
names(Featuredata)[3] <- "Symbol"
fData(eset) <- Featuredata
##subset out AFFX probes in eset
neweset <- subset(eset, substr(rownames(probes), 0, 4) != "AFFX")
#map feature data file to no-AFFX eset
rownames(neweset) <- Featuredata$gene_id
#created new eset with feature info and no AFFX probes
featureEset <- neweset

annotation(featureEset) <- "rna"

rna <- eset


a <- pData(eset) #2605 entries
cell <- arrange(cell, cellid)
all(cell$cellid == a$cellid)
pData(rna) <- a
rm(a)



#test <- eset@assayData$exprs
#test<-subset(test,select=missingCEL$X)
