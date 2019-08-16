#' OPEN-TG-GATES Database - Rat (in vitro)
#' 
#'Program generates object featureInfo that maps featuredata to eSet to build ToxicoSet.
#'Input include Lab annotations

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

storageMode(eset) <- "environment"
affxrows <- rownames(eset@assayData$exprs)
rownames(eset@assayData$exprs) <- substr(rownames(eset@assayData$exprs), 1, nchar(affxrows)-3)

saveRDS(affxrows, file = "rds/probesRatVitro.rds")
CELgenes <- readRDS("rds/probesRatVitro.rds")

results <- getBM(attributes=c("external_gene_name","ensembl_gene_id","gene_biotype","entrezgene_id","external_transcript_name","ensembl_transcript_id"), filters = "ensembl_gene_id",values=CELgenes,mart=ensembl) 

# Account for duplicated genes
uniqueB <- results[!duplicated(results$ensembl_gene_id),]
# Add column for merging
#uniqueB$X <- unique$ensembl_gene_id

CELnotB <- unique(CELgenes)[!unique(CELgenes) %in% uniqueB$ensembl_gene_id]  #44 not in uniqueB

newLabAnnot<-subset(labAnnot,labAnnot$X %in% CELnotB,select=c(X,gene_biotype, gene_id, gene_name,transcript_id, transcript_name,EntrezGene.ID))

names(uniqueB) <- c("gene_name", "gene_id", "gene_biotype", "EntrezGene.ID", "transcript_name", "transcript_id")

# Merge time
finalFeature <- uniqueB

#remainGenes<-CELgenes[!(CELgenes %in% finalFeature$X)]
#leftoverGenes<-as.data.frame(remainGenes,row.names=NULL)

# Create empty columns
#leftoverGenes$A<-NA
#leftoverGenes$B<-NA
#leftoverGenes$C<-leftoverGenes$remainGenes
#leftoverGenes$D<-NA
#leftoverGenes$E<-NA
#leftoverGenes$G<-NA
# Rename columns
#names(leftoverGenes)<-c("X","gene_biotype","gene_id","gene_name","transcript_id","transcript_name", "EntrezGene.ID")
# Merge
#finalFeature<-rbind(finalFeature,leftoverGenes)
# Reformat finalFeature
#finalFeature<-subset(finalFeature,select=-c(X))
#names(finalFeature)[3] <- "Symbol"
finalFeature$BEST <- NA
names(finalFeature) <- c("Symbol", "gene_id", "gene_biotype", "EntrezGene.ID", "transcript_name", "transcript_id", "BEST")
rownames(finalFeature) <- finalFeature$gene_id

#Write into file
saveRDS(finalFeature,"rds/featureData.rds")
