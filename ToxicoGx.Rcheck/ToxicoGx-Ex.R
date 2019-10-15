pkgname <- "ToxicoGx"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "ToxicoGx-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('ToxicoGx')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("ToxicoSet-class")
### * ToxicoSet-class

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ToxicoSet-class
### Title: Class to contain Toxicogenomic Data
### Aliases: ToxicoSet-class .ToxicoSet cellInfo,ToxicoSet-method
###   cellInfo<-,ToxicoSet,data.frame-method drugInfo,ToxicoSet-method
###   drugInfo<-,ToxicoSet,data.frame-method
###   phenoInfo,ToxicoSet,character-method
###   phenoInfo<-,ToxicoSet,character,data.frame-method
###   molecularProfiles,ToxicoSet,character-method
###   molecularProfiles<-,ToxicoSet,character,matrix-method
###   featureInfo,ToxicoSet,character-method
###   featureInfo<-,ToxicoSet,character,data.frame-method
###   sensitivityInfo,ToxicoSet-method
###   sensitivityInfo<-,ToxicoSet,data.frame-method
###   sensitivityProfiles,ToxicoSet-method
###   sensitivityProfiles<-,ToxicoSet,data.frame-method
###   sensitivityProfiles<-,ToxicoSet,matrix-method
###   sensitivityMeasures,ToxicoSet-method drugNames,ToxicoSet-method
###   drugNames<-,ToxicoSet,character-method cellNames,ToxicoSet-method
###   cellNames<-,ToxicoSet,character-method
###   fNames,ToxicoSet,character-method dateCreated,ToxicoSet-method
###   cSetName,ToxicoSet-method pertNumber,ToxicoSet-method
###   sensNumber,ToxicoSet-method pertNumber<-,ToxicoSet,array-method
###   sensNumber<-,ToxicoSet,matrix-method

### ** Examples

data(TGGATESsmall)
cellInfo <- cellInfo(TGGATESsmall)

data(TGGATESsmall)
cellInfo(TGGATESsmall) <- cellInfo(TGGATESsmall)

data(TGGATESsmall)
phenoInfo <- phenoInfo(TGGATESsmall, mDataType="rna")

data(TGGATESsmall)
phenoInfo(TGGATESsmall, mDataType="rna") <- phenoInfo(TGGATESsmall, mDataType="rna")

data(TGGATESsmall)
TGGATES_mProf <- molecularProfiles(TGGATESsmall, "rna")[1:10,]

molecularProfiles(TGGATESsmall, "rna") <- molecularProfiles(TGGATESsmall, "rna")

data(TGGATESsmall)
featureInfo <- featureInfo(TGGATESsmall, "rna")[1:10,]

data(TGGATESsmall)
featureInfo(TGGATESsmall, "rna") <- featureInfo(TGGATESsmall, "rna")

sensInf<- sensitivityInfo(TGGATESsmall)[1:10,]

data(TGGATESsmall)
sensitivityInfo(TGGATESsmall) <- sensitivityInfo(TGGATESsmall)

data(TGGATESsmall)
sensProf <- sensitivityProfiles(TGGATESsmall)

sensitivityProfiles(TGGATESsmall) <- sensitivityProfiles(TGGATESsmall)

sensitivityMeasures(TGGATESsmall)

cellNames(TGGATESsmall)

data(TGGATESsmall)
cellNames(TGGATESsmall) <- cellNames(TGGATESsmall)

fNames(TGGATESsmall, "rna")[1:10]

dateCreated(TGGATESsmall)

tSetName <- cSetName
tSetName(TGGATESsmall)

pertNumber(TGGATESsmall)

sensNumber(TGGATESsmall)

pertNumber(TGGATESsmall) <- pertNumber(TGGATESsmall)

sensNumber(TGGATESsmall) <- sensNumber(TGGATESsmall)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ToxicoSet-class", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("checkTSetStructure")
### * checkTSetStructure

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: checkTSetStructure
### Title: A function to verify the structure of a ToxicoSet
### Aliases: checkTSetStructure

### ** Examples


checkTSetStructure(TGGATESsmall)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("checkTSetStructure", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("computeAUC")
### * computeAUC

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: computeAUC
### Title: Computes the AUC for a Drug Dose Viability Curve
### Aliases: computeAUC

### ** Examples

dose <- c("0.0025","0.008","0.025","0.08","0.25","0.8","2.53","8")
viability <- c("108.67","111","102.16","100.27","90","87","74","57")
computeAUC(dose, viability)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("computeAUC", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("computeICn")
### * computeICn

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: computeIC50
### Title: Computes the ICn for any n in 0-100 for a Drug Dose Viability
###   Curve
### Aliases: computeIC50 computeICn

### ** Examples

dose <- c("0.0025","0.008","0.025","0.08","0.25","0.8","2.53","8")
viability <- c("108.67","111","102.16","100.27","90","87","74","57")
computeIC50(dose, viability)
computeICn(dose, viability, n=10)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("computeICn", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dim-ToxicoSet-method")
### * dim-ToxicoSet-method

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dim,ToxicoSet-method
### Title: Get the dimensions of a ToxicoSet
### Aliases: dim,ToxicoSet-method

### ** Examples

data(TGGATESsmall)
dim(TGGATESsmall)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dim-ToxicoSet-method", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("drugDoseResponseCurve")
### * drugDoseResponseCurve

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: drugDoseResponseCurve
### Title: Plot drug response curve of a given drug and a given cell for a
###   list of tSets (objects of the ToxicoSet class).
### Aliases: drugDoseResponseCurve

### ** Examples

if (interactive()) {
drugDoseResponseCurve(concentrations=list("Experiment 1"=c(.008, .04, .2, 1)),
 viabilities=list(c(100,50,30,1)), plot.type="Both")
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("drugDoseResponseCurve", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("drugGeneResponseCurve")
### * drugGeneResponseCurve

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: drugGeneResponseCurve
### Title: Compares gene expression for a specificed set of features over
###   specific drug dosages vs time
### Aliases: drugGeneResponseCurve

### ** Examples


if (interactive()) {
drugGeneResponseCurve(TGGATESsmall, dose = c("Control", "Low", "Middle"),
  mDataTypes="rna", drug = "naphthyl isothiocyanate",
  duration = c("2", "8", "24"), features = "ENSG00000000003_at")
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("drugGeneResponseCurve", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("drugInfo-set")
### * drugInfo-set

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: drugInfo<-
### Title: drugInfo<- Generic
### Aliases: drugInfo<-

### ** Examples

data(TGGATESsmall)
drugInfo(TGGATESsmall) <- drugInfo(TGGATESsmall)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("drugInfo-set", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("drugInfo")
### * drugInfo

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: drugInfo
### Title: drugInfo Generic
### Aliases: drugInfo

### ** Examples

data(TGGATESsmall)
drugInfo <- drugInfo(TGGATESsmall)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("drugInfo", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("drugNames-set")
### * drugNames-set

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: drugNames<-
### Title: drugNames<- Generic
### Aliases: drugNames<-

### ** Examples

data(TGGATESsmall)
drugNames(TGGATESsmall) <- drugNames(TGGATESsmall)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("drugNames-set", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("drugNames")
### * drugNames

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: drugNames
### Title: drugNames Generic
### Aliases: drugNames

### ** Examples

data(TGGATESsmall)
drugName <- drugNames(TGGATESsmall)[1:10]




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("drugNames", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("drugPerturbationSig")
### * drugPerturbationSig

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: drugPerturbationSig
### Title: Creates a signature representing gene expression (or other
###   molecular profile) change induced by administrating a drug, for use
###   in drug effect analysis.
### Aliases: drugPerturbationSig

### ** Examples

#data(TGGATES_small)
#drug.perturbation <- drugPerturbationSig(TGGATES_small, mDataType="rna", nthread=1)
#print(drug.perturbation)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("drugPerturbationSig", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("drugResponseCurve")
### * drugResponseCurve

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: drugResponseCurve
### Title: Plot drug response curve of a given drug and a given cell for a
###   list of tSets (objects of the ToxicoSet class).
### Aliases: drugResponseCurve

### ** Examples

if (interactive()) {
drugResponseCurve(concentrations=list("Experiment 1"=c(.008, .04, .2, 1)),
 viabilities=list(c(100,50,30,1)), plot.type="Both")
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("drugResponseCurve", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("drugTimeResponseCurve")
### * drugTimeResponseCurve

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: drugTimeResponseCurve
### Title: Compares viabilities at a given dose over different experimental
###   duration
### Aliases: drugTimeResponseCurve

### ** Examples

if (interactive()) {
  ToxicoGx::drugTimeResponseCurve(TGGATESsmall, cellline = "Hepatocyte",
    dose = c("Control", "Low", "Middle"),
    drug = "naphthyl isothiocyanate", duration = c("2", "8", "24"))
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("drugTimeResponseCurve", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("logLogisticRegression")
### * logLogisticRegression

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: logLogisticRegression
### Title: Fits curves of the form E = E_inf + (1 - E_inf)/(1 +
###   (c/EC50)^HS) to dose-response data points (c, E) given by the user
###   and returns a vector containing estimates for HS, E_inf, and EC50.
### Aliases: logLogisticRegression

### ** Examples

dose <- c("0.0025","0.008","0.025","0.08","0.25","0.8","2.53","8")
viability <- c("108.67","111","102.16","100.27","90","87","74","57")
computeAUC(dose, viability)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("logLogisticRegression", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mDataNames-ToxicoSet-method")
### * mDataNames-ToxicoSet-method

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mDataNames,ToxicoSet-method
### Title: mDataNames
### Aliases: mDataNames,ToxicoSet-method

### ** Examples

mDataNames(TGGATESsmall)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mDataNames-ToxicoSet-method", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("names-ToxicoSet-method")
### * names-ToxicoSet-method

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: names,ToxicoSet-method
### Title: tSet Name
### Aliases: names,ToxicoSet-method

### ** Examples

names(TGGATESsmall)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("names-ToxicoSet-method", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("show-ToxicoSet-method")
### * show-ToxicoSet-method

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: show,ToxicoSet-method
### Title: Show a ToxicoSet
### Aliases: show,ToxicoSet-method

### ** Examples

TGGATESsmall




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("show-ToxicoSet-method", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("show-ToxicoSig-method")
### * show-ToxicoSig-method

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: show,ToxicoSig-method
### Title: Show ToxicoGx Signatures
### Aliases: show,ToxicoSig-method

### ** Examples

data(TGGATESsmall)
drug.perturbation <- drugPerturbationSig(TGGATESsmall, mDataType="rna", nthread=1, duration = "2",
     drugs = head(drugNames(TGGATESsmall)), features = fNames(TGGATESsmall, "rna")[seq_len(2)])
drug.perturbation




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("show-ToxicoSig-method", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("showSigAnnot")
### * showSigAnnot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: showSigAnnot
### Title: Show the Annotations of a signature object
### Aliases: showSigAnnot

### ** Examples

data(TGGATESsmall)
drug.perturbation <- drugPerturbationSig(TGGATESsmall, mDataType="rna", nthread=1, duration = "2",
     drugs = head(drugNames(TGGATESsmall)), features = fNames(TGGATESsmall, "rna")[seq_len(2)])
showSigAnnot(drug.perturbation)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("showSigAnnot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sub-ToxicoSet-ANY-ANY-ANY-method")
### * sub-ToxicoSet-ANY-ANY-ANY-method

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: [,ToxicoSet,ANY,ANY,ANY-method
### Title: '['
### Aliases: [,ToxicoSet,ANY,ANY,ANY-method

### ** Examples

tSet <- TGGATESsmall[cellNames(TGGATESsmall), drugNames(TGGATESsmall)[1:3]]




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sub-ToxicoSet-ANY-ANY-ANY-method", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("subsetTo")
### * subsetTo

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: subsetTo
### Title: A function to subset a ToxicoSet to data containing only
###   specified drugs, cells and genes
### Aliases: subsetTo

### ** Examples

TGGATESDrugNames  <- drugNames(TGGATESsmall)
TGGATESCells <- cellNames(TGGATESsmall)
tSet <- subsetTo(TGGATESsmall,drugs = TGGATESDrugNames[1],
  cells = TGGATESCells[1], duration = "2")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("subsetTo", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summarizeMolecularProfiles")
### * summarizeMolecularProfiles

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summarizeMolecularProfiles
### Title: Takes molecular data from a ToxicoSet, and summarises them into
###   one entry per drug and experimental condition.
### Aliases: summarizeMolecularProfiles

### ** Examples

data(TGGATESsmall)
summMP <- ToxicoGx::summarizeMolecularProfiles(
  tSet = TGGATESsmall, mDataType = "rna",
  cell.lines=cellNames(TGGATESsmall), drugs = head(drugNames(TGGATESsmall)),
  features = fNames(TGGATESsmall,"rna"), duration = "8",
  dose = c("Control", "High"), summary.stat = "median",
  fill.missing = TRUE, verbose=TRUE
  )

#subset into expression matrix for a requested drug
assays <- SummarizedExperiment::assays(summMP)[[drugNames(TGGATESsmall)[1]]]
#summarization of phenoData for requested experiments
phenoData <- SummarizedExperiment::colData(summMP)
#summarization of phenoData for requested experiments
featureData <- SummarizedExperiment::rowData(summMP) #featureData for requested experiments




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summarizeMolecularProfiles", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summarizeSensitivityProfiles")
### * summarizeSensitivityProfiles

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summarizeSensitivityProfiles
### Title: Takes the sensitivity data from a ToxicoSet, and summarises them
###   into a drug vs cell line table
### Aliases: summarizeSensitivityProfiles

### ** Examples

data(TGGATESsmall)
TGGATESauc <- summarizeSensitivityProfiles(TGGATESsmall, sensitivity.measure='auc_recomputed')




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summarizeSensitivityProfiles", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("updateCellId")
### * updateCellId

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: updateCellId
### Title: A function to update cell ids
### Aliases: updateCellId
### Keywords: internal

### ** Examples

data(TGGATESsmall)
updateCellId(TGGATESsmall, new.ids = cellNames(TGGATESsmall))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("updateCellId", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("updateDrugId")
### * updateDrugId

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: updateDrugId
### Title: A function to update drug ids
### Aliases: updateDrugId
### Keywords: internal

### ** Examples

data(TGGATESsmall)
updateDrugId(TGGATESsmall, new.ids = drugNames(TGGATESsmall))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("updateDrugId", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
