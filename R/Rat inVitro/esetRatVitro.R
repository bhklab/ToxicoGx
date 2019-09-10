#' OPEN TG-GATES Database - Rat (in vitro)
#' 
#' Program generates object ExpressionSet that is further used to build ToxicoSet.
#' Inputs include gene expression files (.CEL), phenodata, and featuredata.
#' 
#' Data is pre-processed using RMA Normalization Method.

# Install libraries "affy", "Biobase", "Pharmacogx" from Bioconductor
# Load libraries 
library(Biobase)
library(affy)
library(PharmacoGx)


#----------------------------------------------------------------------#
## START FROM SCRATCH. SKIP IF NORMALIZED eSET OBJECT ALREADY IN FILE ##

# Install CDF Normalization Package - Rat 230 2.0 Array Version 23.0 (Item 122 Source C) 
# Download link: http://bit.ly/2Jzu57r
install.packages("data/rat2302rnensgcdf_23.0.0.tar.gz", repos = NULL, source = 'source')
# Load cdf package
library(rat2302rnensgcdf)  
# Call names of CEL files from folder
celFiles <- list.celfiles("data/RatInVitroCEL", full.names = TRUE)
# Normalize CEL files (RMA method)
esetNorm <- just.rma(filenames = celFiles, verbose = TRUE, cdfname = "rat2302rnensgcdf")
# Save object
saveRDS(esetNorm, file = "rds/esetNorm.rds")  

eset <- readRDS("rds/esetNorm.rds")


#-------------------------#
## CONTINUE HERE INSTEAD ##

# Generate eset so far, ONLY includes normalized expression
eset <- readRDS("rds/esetNorm.rds")
# Read phenoinfo
phenoData <- readRDS("rds/phenoData.rds")
# Read featuredata
featureData <- readRDS("rds/featureData.rds")

# Unlock eset@assyadata to modify for desired info 
storageMode(eset) <- "environment"

#delete 4 samples/columns from AssayData, protocolData
eset@assayData$exprs <- subset(eset@assayData$exprs, select=!(colnames(eset@assayData$exprs) %in% c('003017922030_2.CEL', '003017923001.CEL', '003017923004.CEL', '003017923005.CEL')))
eset@assayData$se.exprs <- subset(eset@assayData$se.exprs, select=!(colnames(eset@assayData$se.exprs) %in% c('003017922030_2.CEL', '003017923001.CEL', '003017923004.CEL', '003017923005.CEL')))
eset@protocolData@data <- subset(eset@protocolData@data, !(rownames(eset@protocolData@data) %in% c('003017922030_2.CEL', '003017923001.CEL', '003017923004.CEL', '003017923005.CEL')))


# Rename rows/columns by substring of current names
rownames(eset@protocolData@data) <- substr(rownames(eset@protocolData@data),3,12)
colnames(eset@assayData$exprs) <- substr(colnames(eset@assayData$exprs),3,12)
colnames(eset@assayData$se.exprs) <- substr(colnames(eset@assayData$se.exprs),3,12)
# Subset probes
eset <- subset(eset, substr(rownames(eset@assayData$exprs), 0, 4) != "AFFX")
saveRDS(eset, file = "rds/eset.rds")

# Add phenoData
pData(eset) <- phenoData
# Add featureData
fData(eset) <- featureData
# Replace '_at' with nothing
rownames(eset) <- gsub("_at","",rownames(eset))
# Lock eset@assayData environment when done
storageMode(eset)  <- "lockedEnvironment"

# Set annotation info
annotation(eset) <- "rna"



# TESTS #
#length(intersect(pData(eset)[,"cellid"] , cell$cellid))
#length(intersect(pData(eset)[,"cellid"] , cell$cellid))
#length(intersect(pData(eset)[,"cellid"] , cell$cellid))
#length(intersect(unique(sensitivityInfo$cellid) , cell$cellid))

# Save eset as .rds file
saveRDS(eset, file="rds/eset.rds")

