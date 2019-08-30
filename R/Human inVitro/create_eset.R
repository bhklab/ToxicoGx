# To normalize gene expression data and export as expression set in .RDS file form.

# Install Affymetrix chip (cdf) package to manage data arrays of human hepatocytes. 
### Name of package: Human Genome U133 Plus 2.0 Array, found at http://www.affymetrix.com/support/technical/byproduct.affx?product=hg-u133-plus
install.packages("data/hgu133plus2hsensgcdf_22.0.0.tar.gz", repos = NULL, type = "source")

# Install Affymetrix package from BioConductor.
# Install PharmacoGx package from BioConductor.

# Load libraries installed.
setwd("data")
library("affy")
library("hgu133plus2hsensgcdf")
library("PharmacoGx")

# Group all .CEL files in one variable.
celfn <- list.celfiles("data/CELfiles", full.names = TRUE)

# Perfrom RMA Normalization.
eset <- just.rma(filenames = celfn, verbose=TRUE, cdfname = "hgu133plus2hsensgcdf")

# Save expression set (eset) as .rds file.
saveRDS(eset, file='eset.rds')

# Biobase
# pre-requisites
### info about each experiment
message("Read sample information")
library(Biobase)

path.sc=file.path("data")

# Read .txt file of feature info from metadata.
sampleinfo <- read.table(file.path(path.sc, "TGGATEsfeatureInfo.txt"), sep="\t")
sampleinfo[sampleinfo == "" | sampleinfo == " "] <- NA

# Read annotations file and feature info .csv files.
annot <- read.csv("data/phenowLot.csv", stringsAsFactors=FALSE, check.names=FALSE, header=TRUE, row.names=1)
Featuredata <- read.csv("data/finalFeatureNAs.csv", stringAsFactors=FALSE, check.names=FALSE, header=TRUE, row.names=1)

# Map phenodata to rows and columns of eset.
pData(eset) <- as.data.frame(sampleinfo[match(gsub("[.]CEL[.]gz$", "", rownames(pData(eset))), rownames(sampleinfo)), , drop=FALSE])
colnames(exprs(eset)) <- rownames(pData(eset)) <- gsub("[.]CEL[.]gz$", "", colnames(exprs(eset)))
controls <- rownames(exprs(eset))[grep("AFFX", rownames(exprs(eset)))]

# Map feature data to eset.
fData(eset) <- Featuredata

# Identify probe names in eset.
probes <- eset@assayData$exprs

# Subset eset without AFFX control probes.
neweset <- subset(eset, substr(rownames(probes), 0, 4) != "AFFX")

# Map probe names from feature data file to no-AFFX eset.
rownames(neweset) <- Featuredata$gene_id

# Map new eset (no AFFX probes) to feature info .
featureEset <- neweset

# Save the eset, read it.
saveRDS(featureEset, file='feature_eset.rds')
final <- readRDS('feature_eset.rds')
final

####### 20/06/2019 ########
#load library and eset
library(Biobase)
generate_eset <- just.rma(filenames = celfn, verbose=TRUE, cdfname = "hgu133plus2hsensgcdf")
eset <- readRDS('C:/Users/amytx/OneDrive - UHN/tggates/eset.rds')

#pdata
path.sc=file.path("C:/Users/amytx/OneDrive - UHN/tggates")
sampleinfo <- read.table(file.path(path.sc, "TGGATEsfeatureInfo.txt"), sep="\t")
sampleinfo[sampleinfo == "" | sampleinfo == " "] <- NA
annot <- read.csv("C:/Users/amytx/OneDrive - UHN/tggates/phenoDataCOLNAMES_DONE.csv", stringsAsFactors=FALSE, check.names=FALSE, header=TRUE, row.names=1)
pData(eset) <- annot
pData(eset)[,"batchid"] <- NA


##fdata
#identify probe names
probes <- eset@assayData$exprs

#map ENSG probes in eset to feature data file of ENSG probes only
Featuredata <- read.csv("C:/Users/amytx/OneDrive - UHN/tggates/finalFeatureNAs.csv")
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

#save the eset, read it
saveRDS(featureEset, file='feature_eset.rds')
final <- readRDS('feature_eset.rds')
final

length(intersect(pData(eset)[,"cellid"] , cell$cellid))
length(intersect(pData(eset)[,"cellid"] , cell$cellid))
length(intersect(pData(eset)[,"cellid"] , cell$cellid))
length(intersect(unique(sensitivityInfo$cellid) , cell$cellid))

rna <- eset

