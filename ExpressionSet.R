library(ToxicoGx)
setwd("~Desktop/practicedata/")
setwd("~Desktop/practicedata")


#phenodata
phenodata <- read.csv(file = "s_Hepatocyte.csv", head = TRUE, sep="\t")
aphenodata
hep.filename <- sort(Hep.Metadata$Characteristics.Subject.ID., decreasing=FALSE,)

is.vector(data$Characteristics.Subject.ID.)

phenodata$Factor.Value.Compound.

unique(phenodata$Factor.Value.Compound.)

paste0(hep.filename,".CEL")

#assaydata 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("affy")

browseVignettes("affy")
