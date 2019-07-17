#' OPEN TG-GATES Database - Rat (in vitro)
#' 
#' Program generates object PhenoInfo that is used to map phenodata during creation of eSet
#' Inputs include attributes data (.tsv) from LSDB Database


# Load Phenodata attributes from LSDB Arhive
allAttributes <- read.delim("data/Open-tggates_AllAttribute.tsv", stringsAsFactors = F)
# Subset for SPECIES 'Rat' with DOSE_UNIT 'ƒÊg/kg'. Subset out BARCODE 'No ChipData'.
allAttributes <- subset(allAttributes,allAttributes$SPECIES == "Rat" & allAttributes$DOSE_UNIT == "ƒÊg/kg" & allAttributes$BARCODE != "No ChipData")
