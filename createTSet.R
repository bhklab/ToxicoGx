sensitivityInfo <- readRDS("data/sensitivityInfo.rds")
sensitivityProfiles_biorep <- readRDS("data/sensitivityProfiles_biorep.rds")
sensitivityRaw <- readRDS("data/sensitivityRaw.rds")

cell <- readRDS("data/cell.rds")
drug <- readRDS("data/drug.rds")

curationCell <- readRDS("data/curationCell.rds")
curationDrug <- readRDS("data/curationDrug.rds")
curationTissue <- readRDS("data/curationTissue.rds")

rna <- readRDS("data/eset.rds")

TGGATES <- PharmacoSet("TGGATES",
                       molecularProfiles=list("rna"=rna),
                       cell=cell,
                       drug=drug,
                       sensitivityInfo=sensitivityInfo,
                       sensitivityRaw=sensitivityRaw,
                       sensitivityProfiles=sensitivityProfiles_biorep,
                       curationDrug=curationDrug,
                       curationCell=curationCell,
                       curationTissue=curationTissue,
                       datasetType=c("both"),
                       verify = TRUE)

saveRDS(TGGATES,"TGGATES.rds")
