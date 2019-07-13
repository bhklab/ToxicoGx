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
library(Biobase)

# Group all .CEL files in one variable.
celfn <- list.celfiles("data/CELfiles", full.names = TRUE)

#generate eset 
generated_eset <- just.rma(filenames = celfn, verbose=TRUE, cdfname = "hgu133plus2hsensgcdf")
##OR
#read existing eset
eset <- readRDS("data/eset.rds")

#pdata
path.sc=file.path("data")
sampleinfo <- read.table(file.path(path.sc, "TGGATEsfeatureInfo.txt"), sep="\t")
sampleinfo[sampleinfo == "" | sampleinfo == " "] <- NA
annot <- read.csv("data/phenoDataCOLNAMES_DONE.csv", stringsAsFactors=FALSE, check.names=FALSE, header=TRUE, row.names=1)
pData(eset) <- annot
pData(eset)[,"batchid"] <- NA


##fdata
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

#save the eset, read it
saveRDS(featureEset, file='feature_eset.rds')
final <- readRDS('feature_eset.rds')
final

#tests
length(intersect(pData(eset)[,"cellid"] , cell$cellid))
length(intersect(pData(eset)[,"cellid"] , cell$cellid))
length(intersect(pData(eset)[,"cellid"] , cell$cellid))
length(intersect(unique(sensitivityInfo$cellid) , cell$cellid))

rna <- eset
