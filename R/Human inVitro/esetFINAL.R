# To normalize gene expression data and export as expression set in .RDS file form.

# Install Affymetrix chip (cdf) package to manage data arrays of human hepatocytes. 
### Name of package: Human Genome U133 Plus 2.0 Array, found at http://www.affymetrix.com/support/technical/byproduct.affx?product=hg-u133-plus
install.packages("data/hgu133plus2hsensgcdf_22.0.0.tar.gz", repos = NULL, type = "source")

# Install Affymetrix package from BioConductor.
# Install PharmacoGx package from BioConductor.

# Load libraries installed.
library("affy")
library("hgu133plus2hsensgcdf")
library("PharmacoGx")

# Group all .CEL files in one variable.
celfn <- list.celfiles("data/CELfiles", full.names = TRUE)

# Perform RMA Normalization.
eset <- just.rma(filenames = celfn, verbose=TRUE, cdfname = "hgu133plus2hsensgcdf")

###########################################################################################################

#eset up to this point (only normalization done):
eset<-readRDS("data/esetNORMALIZED_ONLY.rds")

# Load phenoData
phenoData <- readRDS("data/phenoData.rds")
# Load featureData
featureData <- readRDS("data/featureData.rds")

#unlock eset@assayData
storageMode(eset)<-"environment"
#rename eset@assayData column names
colnames(eset@assayData$exprs)<-substr(colnames(eset@assayData$exprs),3,12)
colnames(eset@assayData$se.exprs)<-substr(colnames(eset@assayData$se.exprs),3,12)
rownames(eset@protocolData@data)<-substr(rownames(eset@protocolData@data),3,12)
#subset eset@assayData : 2605 columns/samples -> 2573
#missingCEL is a data.frame containing the info for 2573 samples
missingCEL <- as.data.frame(subset(phenoData,phenoData$concentration_old_units != "ƒÊg/kg", select=c(samplename)))
#missingCEL <- apply(missingCEL, MARGIN=2, str)
missingCEL["samplename",]<-toString(missingCEL["samplename",])
missingCEL <- head(missingCEL, 2573)
#subsetting samples
eset@assayData$exprs<-subset(eset@assayData$exprs,select=missingCEL$samplename)
eset@assayData$se.exprs<-subset(eset@assayData$se.exprs,select=missingCEL$samplename)
eset@protocolData@data<-subset(eset@protocolData@data,rownames(eset@protocolData@data) %in% missingCEL$samplename)
#subsetting probes
eset<-subset(eset, substr(rownames(eset@assayData$exprs), 0, 4) != "AFFX")
# Add phenoData
pData(eset)<-phenoData
# Add featureData
fData(eset)<-featureData
#replace _at
rownames(eset)<-gsub("_at","",rownames(eset))
#lock eset@assayData environment again
storageMode(eset)<-"lockedEnvironment"



annotation(eset)<-"rna"

# TESTS #
#length(intersect(pData(eset)[,"cellid"] , cell$cellid))
#length(intersect(pData(eset)[,"cellid"] , cell$cellid))
#length(intersect(pData(eset)[,"cellid"] , cell$cellid))
#length(intersect(unique(sensitivityInfo$cellid) , cell$cellid))

# Save expression set (eset) as .rds file.
saveRDS(eset, file="data/eset.rds")
