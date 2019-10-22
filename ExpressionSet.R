library(ToxicoGx)



#phenodata

phenodata <- read.csv(file = "s_Hepatocyte.csv", head = TRUE, sep="\t")
phenodata
hep.filename <- sort(phenodata$Characteristics.Subject.ID., decreasing=FALSE,)

is.vector(phenodata$Characteristics.Subject.ID.)

phenodata$Factor.Value.Compound.

unique(phenodata$Factor.Value.Compound.)

paste0(hep.filename,".CEL")


cdf <- read.csv(file = "i_Investigation-hepatocyte.csv", head = TRUE, sep="\t")
cdf
#assaydata 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("affy")

browseVignettes("affy")

#import only specific files listed from metadata
#cdf download --? 
#normalize 
#4 - To create featureData, you have to map the Ensembl IDs from your RMA output(eset) above to gene names, Entrez IDs etc.. 
#You can use biomaRt for that.   Let me know if you could not map all features
