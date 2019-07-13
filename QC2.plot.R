#BiocManager::install("biomaRt")
#install.packages("rlang")
library(Biobase)
library(ggplot2)
library(ggfortify)
library(biomaRt)
library(car)

#Load tSet
tset <- readRDS("data/TGGATES.rds")

#Quality Check 2 : Principle component analysis of gene expression data (Control Vs. treated)
##########################################################################################

#Map 100 probes (from paper supp) to EnsemblID using bioMart
CELgenes<-read.csv("data/Probe_ID.csv")
ensembl<-useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="http://uswest.ensembl.org:80/biomart/martservice?type=registry&requestid=biomaRt",ensemblRedirect = FALSE)
results <- getBM(attributes=c("affy_hg_u133_plus_2","external_gene_name","ensembl_gene_id","gene_biotype","entrezgene_id"), filters = "affy_hg_u133_plus_2",values=CELgenes,mart=ensembl)
probes <- unique(results$ensembl_gene_id) #99 entries

#READING THE OUTPUT SEPARATELY FOR CODEOCEAN
#############################################################################################
#Genes of interest
probes <- readRDS("data/PCA_probes.rds")

#extracting data from tset
#featureData
f<-fData(tset@molecularProfiles$rna)
#metaData
p<-pData(tset@molecularProfiles$rna)
#expression values
assay <- exprs(tset@molecularProfiles$rna)

#subsetting the samples - Control, high dose at 24 hr time point
samples <- subset(p$samplename, (p$dose_level == "Control" | p$dose_level == "High") & p$duration == "24 hr")
#cytokines used in the paper
cytokines = c("interferon alpha, human","hepatocyte growth factor, human","interleukin 1 beta, human",
              "interleukin 6, human","transforming growth factor beta 1","TNFƒ¿")

control <- subset(p,p$dose_level=="Control" & p$duration == "24 hr",select=c(samplename))

high <- subset(p,p$dose_level=="High" & p$duration == "24 hr" & !(p$drugid %in% cytokines),select=c(samplename))

cyto <- subset(p,p$dose_level=="High" & p$duration == "24 hr" & p$drugid %in% cytokines,select=c(samplename))

expr <- as.data.frame(t(subset(assay,rownames(assay) %in% probes,select=as.character(samples))))
# expr <- t(expr)
# expr <- as.data.frame(expr)
expr$control<-NA
expr$control[rownames(expr) %in% as.character(control$samplename)] <- 1 #Control
expr$control[rownames(expr) %in% as.character(high$samplename)] <- 2 #High
expr$control[rownames(expr) %in% as.character(cyto$samplename)] <- 3 #Cytokine

#expr2 <- as.data.frame(expr)

#PCA & plotting
tset.pca <- prcomp(as.matrix(expr),scale. = TRUE)

png("results/pca.png", width = 800, height = 600)
scatterplot(x = tset.pca$x[,1], 
            y = tset.pca$x[,2], 
            regLine = FALSE, 
            smooth = FALSE, 
            boxplots = FALSE, 
            groups = expr$control, 
            col = c('dark green','chartreuse1','red'),
            cex = 3, 
            pch = c(20,20,20),
            legend=FALSE,
            xlab="PC1",
            ylab="PC2",
            main="Principal Component Analysis : High Dose, 24h")
legend(5.4,-5.5,legend=c("Control","Experiment","Cytokine"),col=c('dark green','light green','red'),pch=c(20,20,20),cex=0.75)
dev.off()
