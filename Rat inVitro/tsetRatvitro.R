#' Tset creation TG-GATES rat in vitro 
#' Inputs include curation objects.
#' Program generates ToxicoSet for data analysis

cell <- readRDS("rds/cell.rds")
drug <- readRDS("rds/drug.rds")

curationCell <- readRDS("rds/curationCell.rds")
curationDrug <- readRDS("rds/curationDrug.rds")
curationTissue <- readRDS("rds/curationTissue.rds")

rna <- readRDS("rds/eset.rds")

#######################################################################################################################

#DNA
sensitivityInfo <- readRDS("rds/sensitivityInfo_DNA.rds")
sensitivityProfiles <- readRDS("rds/sensitivityProfiles_DNA.rds")
sensitivityRaw <- readRDS("rds/sensitivityRaw_DNA.rds")

TGGATES_ratDNA <- PharmacoSet("TGGATES_ratDNA",
                                molecularProfiles=list("rna"=rna),
                                cell=cell,
                                drug=drug,
                                sensitivityInfo=sensitivityInfo,
                                sensitivityRaw=sensitivityRaw,
                                sensitivityProfiles=sensitivityProfiles,
                                curationDrug=curationDrug,
                                curationCell=curationCell,
                                curationTissue=curationTissue,
                                datasetType=c("both"),
                                verify = TRUE)

saveRDS(TGGATES_ratDNA,file="rds/TGGATES_ratDNA.rds")

#######################################################################################################################

#LDH
sensitivityInfo <- readRDS("rds/sensitivityInfo_LDH.rds")
sensitivityProfiles <- readRDS("rds/sensitivityProfiles_LDH.rds")
sensitivityRaw <- readRDS("rds/sensitivityRaw_LDH.rds")

TGGATES_ratLDH <- PharmacoSet("TGGATES_ratLDH",
                                molecularProfiles=list("rna"=rna),
                                cell=cell,
                                drug=drug,
                                sensitivityInfo=sensitivityInfo,
                                sensitivityRaw=sensitivityRaw,
                                sensitivityProfiles=sensitivityProfiles,
                                curationDrug=curationDrug,
                                curationCell=curationCell,
                                curationTissue=curationTissue,
                                datasetType=c("both"),
                                verify = TRUE)

saveRDS(TGGATES_ratLDH,file="rds/TGGATES_ratLDH.rds")
