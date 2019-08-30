cell <- readRDS("data/cell.rds")
drug <- readRDS("data/drug.rds")

curationCell <- readRDS("data/curationCell.rds")
curationDrug <- readRDS("data/curationDrug.rds")
curationTissue <- readRDS("data/curationTissue.rds")

rna <- readRDS("data/eset.rds")

#######################################################################################################################

#DNA
sensitivityInfo <- readRDS("data/sensitivityInfo_DNA.rds")
sensitivityProfiles <- readRDS("data/sensitivityProfiles_DNA.rds")
sensitivityRaw <- readRDS("data/sensitivityRaw_DNA.rds")

TGGATES_humanDNA <- PharmacoSet("TGGATES_humanDNA",
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

saveRDS(TGGATES_humanDNA,"TGGATES_humanDNA.rds")

#######################################################################################################################

#LDH
sensitivityInfo <- readRDS("data/sensitivityInfo_LDH.rds")
sensitivityProfiles <- readRDS("data/sensitivityProfiles_LDH.rds")
sensitivityRaw <- readRDS("data/sensitivityRaw_LDH.rds")

TGGATES_humanLDH <- PharmacoSet("TGGATES_humanLDH",
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

saveRDS(TGGATES_humanLDH,"TGGATES_humanLDH.rds")
