####### 20/06/2019 #######

PharmacoSet("TGGATES",
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

length(intersect(curationDrug$unique.drugid, rownames(drug))) != nrow(drug)
length(intersect(curationCell$unique.cellid, rownames(cell))) != nrow(cell)

length(intersect(pData(eset)[,"cellid"] , cell$cellid))
length(intersect(pData(eset)[,"cellid"] , cell$cellid))
length(intersect(pData(eset)[,"cellid"] , cell$cellid))
length(intersect(sensitivityInfo$cellid , cell$cellid))

all(rownames(sensitivityInfo) == rownames(sensitivityRaw))
all(rownames(sensitivityInfo) == rownames(sensitivityProfiles))
all(rownames(sensitivityInfo) == dimnames(sensitivityRaw)[[1]])
