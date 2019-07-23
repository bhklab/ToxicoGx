#load packages
library("xml2")
library(biomaRt)

#load files
#LAB ANNOTATION FILE
labAnnot<-read.csv("data/annot_ensembl_all_genes.csv",header=TRUE)
#PROBE NAMES FROM GENE EXPRESSION FILE
CELgenes<-readRDS("data/probeNames.rds")
#Use biomaRt to do preliminary gene annotations
ensembl<-useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="uswest.ensembl.org",ensemblRedirect = FALSE)
results <- getBM(attributes=c("external_gene_name","ensembl_gene_id","gene_biotype","entrezgene_id","external_transcript_name","ensembl_transcript_id"), filters = "ensembl_gene_id",values=CELgenes,mart=ensembl)

#########################################################################################################################################################################################################################

##BIOMART RESULT##
results <- readRDS("data/bioMart_output.rds")

#Account for missing genes (59)
#48 are accounted for by lab annotation file; 11 were left NA
#Select for the first of any duplicate gene names
uniqueBiomaRt<-results[!duplicated(results$ensembl_gene_id),]
#Add another column for merging
uniqueBiomaRt$X<-uniqueBiomaRt$ensembl_gene_id
#From current lab annotations, select for the 48 that are currently unaccounted for after biomaRt
CELnotbiomaRt<-unique(CELgenes)[!(unique(CELgenes) %in% uniqueBiomaRt$ensembl_gene_id)] #in CELgenes but not in uniqueBiomaRt (59)
newLabAnnot<-subset(labAnnot,labAnnot$X %in% CELnotbiomaRt,select=c(X,gene_biotype, gene_id, gene_name,transcript_id, transcript_name,EntrezGene.ID))
#Rename uniqueBiomaRt columns to prepare for merging of data sets
names(uniqueBiomaRt)<-c("gene_name", "gene_id", "gene_biotype", "EntrezGene.ID", "transcript_name", "transcript_id", "X")
#Merging the data sets
finalFeature<-rbind(uniqueBiomaRt,newLabAnnot)
#Making the rest of the 11 genes NA
#Select for the 11 missing genes
remainGenes<-CELgenes[!(CELgenes %in% finalFeature$X)]
#Turn the matrix into a data.set
leftoverGenes<-as.data.frame(remainGenes,row.names=NULL)
#create empty columns
leftoverGenes$A<-NA
leftoverGenes$B<-NA
leftoverGenes$C<-leftoverGenes$remainGenes
leftoverGenes$D<-NA
leftoverGenes$E<-NA
leftoverGenes$G<-NA
#Rename columns
names(leftoverGenes)<-c("X","gene_biotype","gene_name","gene_id","transcript_id","transcript_name", "EntrezGene.ID")
#Merge
finalFeature<-rbind(finalFeature,leftoverGenes)
#reformat finalFeature
finalFeature<-subset(finalFeature,select=-c(X))
names(finalFeature)[3] <- "Symbol"
finalFeature$BEST <- NA
finalFeature<-arrange(finalFeature,finalFeature$gene_id)
rownames(finalFeature)<-finalFeature$gene_id

#Write into file
<<<<<<< HEAD
saveRDS(finalFeature,"data/featureData.rds")
=======
saveRDS(finalFeature,"data/featureData.rds")
>>>>>>> d3e1a831c5b5c97e7362c1600cdd4e9fb40b1923
