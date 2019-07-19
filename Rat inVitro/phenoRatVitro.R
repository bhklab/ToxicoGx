#' OPEN TG-GATES Database - Rat (in vitro)
#' 
#' Program generates object PhenoInfo that is used to map phenodata during creation of eSet
#' Inputs include attributes data (.tsv) from LSDB Database

library(dplyr)

# REMOVE IRRELEVANT DATA #

# Load Phenodata attributes from LSDB Arhive
allAttributes <- read.delim("data/Open-tggates_AllAttribute.tsv", stringsAsFactors = F)
# Subset dataframe for SPECIES 'Rat' with TEST_TYPE 'in vitro'. Eliminate BARCODE 'No ChipData'.
allAttributes <- subset(allAttributes,allAttributes$SPECIES == "Rat" & allAttributes$TEST_TYPE == "in vitro"  &  allAttributes$BARCODE != "No ChipData")
# Remove columns from dataframe where entire columns contain unspecified and 'NA' values
allAttributes[allAttributes == 'NA' | allAttributes == 'not specified'] <- NA  # convert characters to NULL
ratAttributes <- allAttributes[,colSums(is.na(allAttributes)) < nrow(allAttributes)]  # remove NULL columns
# Eliminate 2 extra zeroes beginning BARCODE values
ratAttributes$BARCODE <- substr(ratAttributes$BARCODE,3,12)

# ADD COLUMNS FOR ADDITIONAL DATA #

# Column Batch ID #
# Read mapping values
batchId <- read.csv("data/phenoRatVitro.csv",stringsAsFactors = FALSE)
#Rename columns for convenient mapping
colnames(batchId) <-as.character(c('drugid', 'drugid_abbr', 'drugid_no', 'medium_type' ,'cell','low_dose','med_dose','high_dose','batch_id'))
colnames(batchId) <- as.factor(colnames(batchId))
# Eliminate unnecessary rows and columns
batchId <- batchId[-c(1:3),] 
batchId <- batchId[,-c(10:10)]
# Merge batchId with ratAttributes to create new column in dataframe
ratAttributes <- merge(ratAttributes, batchId, by.x = "COMPOUND_NAME", by.y = "drugid")
colnames(ratAttributes) <- as.factor(colnames(ratAttributes))
# Fix dose units 
ratAttributes$DOSE_UNIT <- "μM"

# Column UID #
ratAttributes$UID <- NA

# MODIFY ROW & COLUMNS #

# Copy BARCODE values to cellid
ratAttributes$cellid <- ratAttributes$BARCODE
# Create UID template
ratAttributes$UID <- paste0("drugid_",ratAttributes$COMPOUND.Abbr.,"_",ratAttributes$COMPOUND_NAME,"_",ratAttributes$BARCODE,"_",ratAttributes$DOSE,"_",gsub(" hr","",ratAttributes$SACRI_PERIOD),"hr_rep",ratAttributes$INDIVIDUAL_ID)
# Copy BARCODE values to CEL file names
ratAttributes$celfilename <- paste0("00",as.character(ratAttributes$BARCODE),".CEL")
# Both sensitivity and perturbation info included
ratAttributes$xptype <- "both"

# Rearrange order of columns
ratAttributes <- ratAttributes[,c(2,3,4,5,6,27,15,16,17,14,29,1,20,21,18,19,28,11,12,13,7,8,30,31)]
# Rename columns
colnames(ratAttributes) <- c('samplename','chiptype','exp_id','group_id','individual_id','batchid','concentration','concentration_units','dose_level','duration','cellid','drugid','drugid_abbr','drugid_no','DNA','LDH','UID','species','test_type','sex_type','organ_id','material_id','celfilename','xptype')

# Arrange rows in ascending order of barcode
ratAttributes <- arrange(ratAttributes, ratAttributes$samplename)
# Rename rows to barcode name 
rownames(ratAttributes) <- ratAttributes$samplename

# RETURN #

# Save as object
saveRDS(ratAttributes, file = "rds/ratAttributes.rds")
