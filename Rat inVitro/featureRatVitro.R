#' OPEN-TG-GATES Database - Rat (in vitro)
#' 
#'Program generates object featureInfo that maps featuredata to eSet to build ToxicoSet.
#'Input include lab annotations (.csv) and probe names from Rat in vitro eSet

library(affy)
library(biomaRt)

# Load Lab Annotations file for genes
labAnnot <- read.csv("data/annot_ensembl_all_genes.csv", header = T)

# List ensembl options
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)

# Choose dataset to use
ensembl = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl", host="uswest.ensembl.org",ensemblRedirect = FALSE)

# Get values of CEL gene probes
eset <- readRDS("rds/eset.rds")
affxrows <- rownames(eset@assayData$exprs)
saveRDS(affxrows, file = "rds/probesRatVitro.rds")
CELgenes <- readRDS("rds/probesRatVitro.rds")

results2 <- getBM(attributes=c("external_gene_name","ensembl_gene_id","gene_biotype","entrezgene_id","external_transcript_name","ensembl_transcript_id"), filters = "affy_rat230_2",values=CELgenes,mart=ensembl) 
