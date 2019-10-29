getwd()
library(g.data)
library(readxl)

install.packages("g.data")
install.packages("readxl")
drugann <- read_xls("drugmatrix.xls", sheet= NULL, col_names = TRUE)

drugann

is.data.frame(drugann)

ncol(drugann)

nrow(drugann)


drugann$`Hybridization Date [C]` <- drugann$`Array Design REF` <- drugann$Organism <- drugann$Strain <- drugann$Sex <- drugann$Organ <- drugann$`Organ or Cell [C]` <- drugann$`Term Accession Number` <- drugann$StdInChIKey <- drugann$`Comment[ChEMBL ID` <- NULL
drugann

ncol(drugann)
str(drugann)

dim(drugann)

drugann
