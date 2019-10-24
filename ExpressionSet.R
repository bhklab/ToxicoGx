library(ToxicoGx)
library(Biobase)
library(affy)
library(affyio)
library(BiocManager)
library(rat2302rnensgcdf)
library(dplyr)

#phenodata

phenodata <- read.csv(file = "s_Hepatocyte.csv", head = TRUE, sep="\t")
phenodata
#phenodata contains 939 relevant CEL files from 5590 total file set - the microarray data for all tissue types 

is.data.frame(phenodata)

phenodata[phenodata == 'NA' | phenodata =='not specified'] <- NA 
phenodataC <- phenodata[,colSums(is.na(phenodata)) < nrow(phenodata)]
ncol(phenodata)
ncol(phenodataC)


  
phenodataC

names(phenodataC) <- gsub("Characteristics.", "", names(phenodataC))
phenodataC

phenodataC$DoseUnit. <- "μM"
phenodataC

names(phenodataC) <- gsub("Factor.Value", "", names(phenodataC))
phenodataC

phenodataC$Organism. <- "R.norvegicus"
phenodataC

names(phenodataC) <- gsub("Term.", "", names(phenodataC))
phenodataC

phenodataC <- phenodataC[,-32]
phenodataC
ncol(phenodataC)

colnames(phenodataC)[colnames(phenodataC) == "CompoDosetDoseReplicate#"] <- "Compound.DoseDuration.Dose.BioReplicate#"
names(phenodataC) <- gsub("Strain.", "Strain", names(phenodataC))
names(phenodataC) <- gsub("Organism.", "Organism", names(phenodataC))
names(phenodataC) <- gsub("Subject.ID.", "Subject.ID", names(phenodataC))
names(phenodataC) <- gsub("Sex.", "Sex", names(phenodataC))
names(phenodataC) <- gsub("Cell.", "Cell", names(phenodataC))
names(phenodataC) <- gsub("Assay.Type.", "Assay.Type", names(phenodataC))
names(phenodataC) <- gsub("Biological.Replicate.", "Biological.Replicate", names(phenodataC))
names(phenodataC) <- gsub("Technical.Replicate.", "Technical.Replicate", names(phenodataC))
names(phenodataC) <- gsub(".Compound.", "Compound.Name", names(phenodataC))
names(phenodataC) <- gsub("Control.", "Control", names(phenodataC))
names(phenodataC) <- gsub("Sample.Match.", "Sample.Match", names(phenodataC))
names(phenodataC) <- gsub(".Dose.", "Dose", names(phenodataC))
names(phenodataC) <- gsub("StdInChIKey.", "StdInChIKey", names(phenodataC))
names(phenodataC) <- gsub("Comment.chEMBL.ID.", "Comment.chEMBL.ID", names(phenodataC))
names(phenodataC) <- gsub("DoseUnit.", "Dose.Unit", names(phenodataC))
names(phenodataC) <- gsub("DoseDuration.", "Dose.Duration", names(phenodataC))
names(phenodataC) <- gsub("Dose.DurationUnit.", "Dose.Duration.Unit", names(phenodataC))
names(phenodataC) <- gsub("Vehicle.", "Vehicle", names(phenodataC))
phenodataC

COLphenodataC <- colnames(phenodataC)
as.vector(COLphenodataC)
cat(COLphenodataC, sep="\n")

numbers <- c(1:31)
numbers
paste0(numbers, COLphenodataC)

hep.filename <- sort(phenodata$Characteristics.Subject.ID., decreasing=FALSE,)
hep.filename

phenodataC


saveRDS(phenodataC, file = "/Users/parwaiznijrabi/Desktop/phenoData.rds")

is.vector(phenodata$Characteristics.Subject.ID.)

phenodata$Factor.Value.Compound.

unique(phenodata$Factor.Value.Compound.)
#all drugs used in the study 


hepfilelist <- paste0(hep.filename, ".CEL")
hepfilelist
#hepfilelist, refer to this for further 



#assaydata 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("affy")

browseVignettes("affy")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationDbi")

install.packages("~/Desktop/rat2302rnensgcdf", repos= NULL, type="source")
library(rat2302rnensgcdf)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("affyio")

##eset <- just.rma(filenames = celfn, verbose = TRUE, cdfname = cdf)

dir <- "/Users/parwaiznijrabi/Desktop/hepatocytemicroarray/hepatocyte.data/liver_drug_matrix"
all.files <- list.files(dir)
#listing all files, storing in all.files 

head(all.files)

celFiles <- list.celfiles("/Users/parwaiznijrabi/Desktop/hepatocytemicroarray/hepatocyte.data/liver_drug_matrix", full.names = TRUE)
#all celfiles, 4540 

all.data.dir <- data.frame(all.files, celFiles)
all.data.dir
#dataframe with .cel extensions and full path 

sub.all.data.dir <- all.data.dir[all.data.dir$all.files %in% hepfilelist, ]
#subset
sub.all.data.dir
is.vector(sub.all.data.dir)

dim(sub.all.data.dir)

setdiff(hepfilelist,sub.all.data.dir$all.files)

hepfilelistasc <- sort(hepfilelist, decreasing = FALSE,)
hepfilelistasc
?sort
head(sub.all.data.dir)

list.celfiles("liver_drug_matrix", full.names = TRUE)

#> sub.all.data.dir[1:4]
#> path.id <- as.character(sub.all.data.dir$celFiles)
#> path.id[1:4]
#> esetNorm <- just.rma(filenames = path.id, verbose = TRUE, cdfname = "rat2302rnensgcdf")
#> args(just.rma)
#> path.id[1]
#> names(esetNorm) NULL>   > class(esetNorm)[1] "ExpressionSet"attr(,"package")[1] "Biobase aa <- exprs(esetNorm) dim(aa) head(sub.all.data.dir)
#> exprs(esetNorm)
#> dir <- "/Users/parwaiznijrabi/Desktop/Hep939"
#> id1 <- list.files(dir)
#> length(id1)
#setdiif()

id1 <- list.files(dir)
dir <- "/Users/parwaiznijrabi/Desktop/Hep939"
#id1 and id are the same in length, id1 is just converted to a list 

length(id1) #length 893
length(id)  #length 893
length(hepfilelist) #length 939


eSet.miss <- setdiff(hepfilelist, id1)
sort(eSet.miss, decreasing = FALSE,)
length(eSet.miss) #<- 30 files excluded from normalization analysis either corrupt or could not be read 


int <- intersect(hepfilelist, id1)
int

setdiff(id1,int)
#893 files in HEP939 set, 46 files missing
#of the 46, 93081, 93180, 93368, 93394, 93639 are ERROR: CORRUPT
# of the 46, 93015, 93074, 93079, 93091, 93113, 93188, 93223, 93283,93310,93321,93363,93366, 93421, 93444, 93479, 93504,93616,93672,93707,93713,93741,93783,93876,93947,93969,94212,94213,94219, 96015,98382 = celfile?
# of the 46 93253,4,5,6,7,8,9,93260, 94226, 95508, neet to be transferred to hep939

#corrupt <- c(93081, 93180, 93368, 93394, 93639,100992) <- NO LONGER CORRUPT, BACK IN DATASET
#iscelfileQ <- c(93015,93074,93079,93091,93113,93188, 93223, 93283,93310,93321,93363,93366, 93421, 93444, 93479, 93504,93616,93672,93707,93713,93741,93783,93876,93947,93969,94212,94213,94219, 96015,98382) <- redownloaded and ran function, worked 
transfer <-  c(93253,93254,93255,93256,93257,93258,93259,93260, 94226, 95508, 100992)
length(corrupt)
length(iscelfileQ)
length(transfer)
#=46

id
id <- paste(dir, list.files(dir), sep="/")
#paste files, why?
esetNorm <- just.rma(filenames = id, verbose = TRUE, cdfname = "rat2302rnensgcdf")
saveRDS(esetNorm, file = "/Users/parwaiznijrabi/Desktop/esetNorm.rds")

class(exprs)
colnames(exprs)
dim(exprs)

summary(phenodataC)
#all(rownames(phenodataC)==colnames(exprs))
#________________________________feature data

library(affy)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
phenodataC
library(biomaRt)
esetNorm
exprs(esetNorm)
browseVignettes("biomaRt")
listMarts()

ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)

ensembl = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl", host="uswest.ensembl.org",ensemblRedirect = FALSE)
eset <- readRDS("/Users/parwaiznijrabi/Desktop/esetNorm.rds")

storageMode(eset) <- "environment"
affxrows <- rownames(eset@assayData$exprs)
rownames(eset@assayData$exprs) <- substr(rownames(eset@assayData$exprs), 1, nchar(affxrows)-3)


saveRDS(affxrows, file = "/Users/parwaiznijrabi/Desktop/probesRatVitro.rds")
CELgenes <- readRDS("/Users/parwaiznijrabi/Desktop/probesRatVitro.rds")
results <- getBM(attributes=c("external_gene_name","ensembl_gene_id","gene_biotype","entrezgene_id","external_transcript_name","ensembl_transcript_id"), filters = "ensembl_gene_id",values=CELgenes,mart=ensembl) 
uniqueB <- results[!duplicated(results$ensembl_gene_id),]


CELnotB <- unique(CELgenes) [!unique(CELgenes) %in% uniqueB$ensembl_gene_id]

LABANNOT <- read.csv("annot_ensembl_all_genes.csv")
CELnotB <- unique(CELgenes) [!unique(CELgenes) %in% uniqueB$ensembl_gene_id]

newLabAnnot <- subset(LABANNOT, LABANNOT$X %in% CELnotB, select=c(X,gene_biotype, gene_id, gene_name, transcript_id, transcript_name,EntrezGene.ID))
names(uniqueB) <- c("gene_name", "gene_id", "gene_biotype", "EntrezGene.ID", "transcript_name", "transcript_id")

finalFeature <- uniqueB


remainGenes<-CELgenes[!(CELgenes %in% finalFeature$X)]
leftoverGenes<-as.data.frame(remainGenes,row.names=NULL)

# Create empty columns
leftoverGenes$A<-NA
leftoverGenes$B<-NA
leftoverGenes$C<-leftoverGenes$remainGenes
leftoverGenes$D<-NA
leftoverGenes$E<-NA
leftoverGenes$G<-NA
# Rename columns
names(leftoverGenes)<-c("X","gene_biotype","gene_id","gene_name","transcript_id","transcript_name", "EntrezGene.ID")
# Merge
finalFeature<-rbind(finalFeature,leftoverGenes)
#Reformat finalFeature
finalFeature<-subset(finalFeature,select=-c(X))
names(finalFeature)[3] <- "Symbol"

finalFeature$BEST <- NA
names(finalFeature) <- c("Symbol", "gene_id", "gene_biotype", "EntrezGene.ID", "transcript_name", "transcript_id", "BEST")
rownames(finalFeature) <- finalFeature$gene_id

finalFeature

saveRDS(finalFeature, "/Users/parwaiznijrabi/Desktop/featureData.rds")





#eset up to this point (only normalization done):
eset1<-readRDS("/Users/parwaiznijrabi/Desktop/esetNorm.rds")

# Load phenoData
phenoData1 <- readRDS("/Users/parwaiznijrabi/Desktop/phenoData.rds")
# Load featureData
featureData1 <- readRDS("/Users/parwaiznijrabi/Desktop/featureData.rds")

#unlock eset@assayData
storageMode(eset1)<-"environment"
#rename eset@assayData column names
#colnames(eset1@assayData$exprs)<-substr(colnames(eset@assayData$exprs),3,12)
#colnames(eset1@assayData$se.exprs)<-substr(colnames(eset@assayData$se.exprs),3,12)
#rownames(eset1@protocolData@data)<-substr(rownames(eset@protocolData@data),3,12)
#subset eset@assayData : 2605 columns/samples -> 2573
#missingCEL is a data.frame containing the info for 2573 samples
#missingCEL <- as.data.frame(subset(phenoData,phenoData$concentration_old_units != "ƒÊg/kg", select=c(samplename)))
#missingCEL <- apply(missingCEL, MARGIN=2, str)
#missingCEL["samplename",]<-toString(missingCEL["samplename",])
#missingCEL <- head(missingCEL, 939)
#subsetting samples
#eset@assayData$exprs<-subset(eset@assayData$exprs,select=missingCEL$samplename)
#eset@assayData$se.exprs<-subset(eset@assayData$se.exprs,select=missingCEL$samplename)
#eset@protocolData@data<-subset(eset@protocolData@data,rownames(eset@protocolData@data) %in% missingCEL$samplename)
#subsetting probes
eset1<-subset(eset1, substr(rownames(eset1@assayData$exprs), 0, 4) != "AFFX")
# Add phenoData
pData(eset1)<-phenoData1
# Add featureData
fData(eset1)<-featureData1
#replace _at
rownames(eset)<-gsub("_at","",rownames(eset))
#lock eset@assayData environment again
storageMode(eset)<-"lockedEnvironment"


annotation(eset1)<-"rna"

# TESTS #
#length(intersect(pData(eset)[,"cellid"] , cell$cellid))
#length(intersect(pData(eset)[,"cellid"] , cell$cellid))
#length(intersect(pData(eset)[,"cellid"] , cell$cellid))
#length(intersect(unique(sensitivityInfo$cellid) , cell$cellid))

# Save expression set (eset) as .rds file.
saveRDS(eset1, file="/Users/parwaiznijrabi/Desktop/ExpressionSetRough.rds")


