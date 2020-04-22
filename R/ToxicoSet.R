#' Class to contain Toxico-genomic Data
#'
#' The ToxicoSet (tSet) class was development to contain and organise large
#' ToxicGenomic datasets as well as provide useful tools for interacting with
#' this data. Functions are included for exploring the relationship between
#' survival fraction and gene expression in cultured human and rat tissues
#' during exposure to a wide ranges of compounds. Features include plotting dose
#' and exposure time curves, calculating AUC, fitting linear models and
#' computing sensitivity signatures.
#'
#' @slot annotation A \code{list} of annotation data about the ToxicoSet,
#'    including the \code{$name} and the session information for how the object
#'    was creating, detailing the exact versions of R and all the packages used
#' @slot molecularProfiles A \code{list} containing \code{SummarizedExperiment}
#'   type object for holding data for RNA, DNA, SNP and CNV
#'   measurements, with associated \code{fData} and \code{pData}
#'   containing the row and column metadata
#' @slot cell A \code{data.frame} containing the annotations for all the cell
#'   lines profiled in the data set, across all data types
#' @slot drug A \code{data.frame} containg the annotations for all the drugs
#'   profiled in the data set, across all data types
#' @slot sensitivity A \code{list} containing all the data for the sensitivity
#'   experiments, including \code{$info}, a \code{data.frame} containing the
#'   experimental info,\code{$raw} a 3D \code{array} containing raw data,
#'   \code{$profiles}, a \code{data.frame} containing sensitivity profiles
#'   statistics, and \code{$n}, a \code{data.frame} detailing the number of
#'   experiments for each cell-drug pair
#' @slot perturbation A \code{list} containting \code{$n}, a \code{data.frame}
#'   summarizing the available perturbation data,
#' @slot curation A \code{list} containing mappings for \code{$drug},
#'   \code{cell}, \code{tissue} names  used in the data set to universal
#'   identifiers used between different ToxicoSet objects
#' @slot datasetType A \code{character} string of 'sensitivity',
#'   'perturbation', or both detailing what type of data can be found in the
#'   ToxicoSet, for proper processing of the data
#'
#' @return An object of the ToxicoSet class
#'
#' @importClassesFrom CoreGx CoreSet
.ToxicoSet <- setClass("ToxicoSet", slots = list(drug="data.frame"),
                       contains="CoreSet")

# The default constructor above does a poor job of explaining the required
# structure of a ToxicoSet. The constructor function defined below guides the
# user into providing the required components of the curation and senstivity
# lists and hides the annotation slot which the user does not need to manually
# fill. This also follows the design of the Expression Set class.

#' ToxicoSet constructor
#'
#' A constructor that simplifies the process of creating ToxicoSets, as well
#' as creates empty objects for data not provided to the constructor. Only
#' objects returned by this constructor are expected to work with the ToxicoSet
#' methods. For a much more detailed instruction on creating ToxicoSets, please
#' see the "CreatingToxicoSet" vignette.
#'
#' @param name A \code{character} string detailing the name of the dataset
#' @param molecularProfiles A \code{list} of ExpressionSet objects containing
#'   molecular profiles
#' @param cell A \code{data.frame} containing the annotations for all the cell
#'   lines profiled in the data set, across all data types
#' @param drug A \code{data.frame} containing the annotations for all the drugs
#'   profiled in the data set, across all data types
#' @param sensitivityInfo A \code{data.frame} containing the information for the
#'   sensitivity experiments
#' @param sensitivityRaw A 3 Dimensional \code{array} contaning the raw drug
#'   dose – response data for the sensitivity experiments
#' @param sensitivityProfiles \code{data.frame} containing drug sensitivity profile
#'   statistics such as IC50 and AUC
#' @param sensitivityN,perturbationN A \code{data.frame} summarizing the
#'   available sensitivity/perturbation data
#' @param curationCell,curationDrug,curationTissue A \code{data.frame} mapping
#'   the names for cells, drugs, and tissues used in the data set to universal
#'   identifiers used between different ToxicoSet objects
#' @param datasetType A \code{character} string of "sensitivity",
#'   "perturbation", or both detailing what type of data can be found in the
#'   ToxicoSet, for proper processing of the data
#' @param verify \code{boolean} Should the function verify the ToxicoSet and
#'   print out any errors it finds after construction?
#'
#' @return An object of class \code{ToxicoSet}
#'
#' @import methods
#' @importFrom utils sessionInfo
#' @importFrom stats na.omit
#' @importFrom SummarizedExperiment rowData colData assay assays
#'   assayNames Assays rowData<- colData<-
#' @importFrom S4Vectors DataFrame SimpleList metadata
#' @importFrom CoreGx CoreSet
#' @export
ToxicoSet <-  function(name,
                       molecularProfiles=list(),
                       cell=data.frame(),
                       drug=data.frame(),
                       sensitivityInfo=data.frame(),
                       sensitivityRaw=array(dim = c(0,0,0)),
                       sensitivityProfiles=matrix(),
                       sensitivityN=matrix(nrow = 0, ncol=0),
                       perturbationN=array(NA, dim = c(0,0,0)),
                       curationDrug=data.frame(),
                       curationCell = data.frame(),
                       curationTissue = data.frame(),
                       datasetType=c("sensitivity", "perturbation", "both"),
                       verify = TRUE)
{
  ##TOOD:: Abstract as much of this contstructor as possible to CoreGx!
  datasetType <- match.arg(datasetType)

  annotation <- list()
  annotation$name <- as.character(name)
  annotation$dateCreated <- date()
  annotation$sessionInfo <- sessionInfo()
  annotation$call <- match.call()

    for (i in seq_len(length(molecularProfiles))){
        if (!is(molecularProfiles[[i]], "SummarizedExperiment")) {
            stop(sprintf("Please provide the %s data as a SummarizedExperiment", names(molecularProfiles[i])))
        } else {
          rowData(molecularProfiles[[i]]) <-
            rowData(molecularProfiles[[i]])[rownames(assays(molecularProfiles[[i]])[[1]]), , drop=FALSE]
          colData(molecularProfiles[[i]]) <-
            colData(molecularProfiles[[i]])[colnames(assays(molecularProfiles[[i]])[[1]]), , drop=FALSE]
        }
    }

  sensitivity <- list()

  if (!all(rownames(sensitivityInfo) == rownames(sensitivityProfiles) &
           rownames(sensitivityInfo) == dimnames(sensitivityRaw)[[1]])) {
    stop("Please ensure all the row names match between the sensitivity data.")
  }

  sensitivity$info <- as.data.frame(sensitivityInfo, stringsAsFactors=FALSE)
  sensitivity$raw <- sensitivityRaw
  sensitivity$profiles <- as.data.frame(sensitivityProfiles, stringsAsFactors=FALSE)
  sensitivity$n <- sensitivityN

  ### TODO:: Make sure to fix the curation to check for matching row names to
  ### the radiation and cell line matrices
  curation <- list()
  curation$cell <- as.data.frame(curationCell, stringsAsFactors = FALSE)
  curation$drug <- as.data.frame(curationDrug, stringsAsFactors = FALSE)
  curation$tissue <- as.data.frame(curationTissue, stringsAsFactors = FALSE)

  perturbation <- list()
  perturbation$n <- perturbationN
  if (datasetType == "perturbation" || datasetType == "both") {
    perturbation$info <- "The metadata for the perturbation experiments is
      available for each molecular type by calling the appropriate info function.
      \n For example, for RNA transcriptome perturbations, the metadata can be
      accessed using rnaInfo(tSet)."
  } else {
    perturbation$info <- "Not a perturbation dataset."
  }

  tSet  <- .ToxicoSet(annotation=annotation,
                      molecularProfiles=molecularProfiles,
                      cell=as.data.frame(cell),
                      drug=as.data.frame(drug),
                      datasetType=datasetType,
                      sensitivity=sensitivity,
                      perturbation=perturbation,
                      curation=curation)
  if (verify) { checkTSetStructure(tSet)}
  if(length(sensitivityN) == 0 & datasetType %in% c("sensitivity", "both")) {
    sensNumber(tSet) <- .summarizeSensitivityNumbers(tSet)
  }
  if(length(perturbationN) == 0  & datasetType %in% c("perturbation", "both")) {
    pertNumber(tSet) <- .summarizePerturbationNumbers(tSet)
  }
  return(tSet)
}

#' name Getter method
#'
#' Retrieves the name of a tSet
#'
#' @examples
#' name(TGGATESsmall)
#'
#' @param object \code{ToxicoSet} A ToxicoSet object
#'
#' @return \code{character} A string of the tSet's name
#'
#' @importFrom CoreGx name
#' @importFrom methods callNextMethod
#' @export
setMethod(name, "ToxicoSet", function(object) {
  callNextMethod(object)
})

#' cellInfo Getter
#'
#' Get the cell line annotations in a ToxicoSet
#'
#' @examples
#' data(TGGATESsmall)
#' cellInfo <- cellInfo(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#'
#' @return a \code{data.frame} with the cell annotations
#'
#' @describeIn ToxicoSet Returns the annotations for all the cell lines tested
#'   on in the ToxicoSet
#'
#' @importFrom CoreGx cellInfo
#' @importFrom methods callNextMethod
#' @export
setMethod(cellInfo, "ToxicoSet", function(object){
  callNextMethod(object)
})

#' cellInfo Replace Method
#'
#' Set cell line annotations for a ToxicoSet object
#'
#' @examples
#' data(TGGATESsmall)
#' cellInfo(TGGATESsmall) <- cellInfo(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#' @param value A \code{data.frame} of replacement values
#'
#' @return Updated \code{ToxicoSet}
#'
#' @describeIn ToxicoSet Returns the annotations for all the cell lines tested
#'   on in the ToxicoSet
#'
#' @importFrom CoreGx cellInfo<-
#' @importFrom methods callNextMethod
#' @export
setReplaceMethod("cellInfo",
                 signature = signature(object = "ToxicoSet",
                                       value = "data.frame"),
                 function(object, value)
{
  callNextMethod(object, value)
})

##TODO:: Abstract this method to CoreGx
#' drugInfo Getter
#'
#' Get the drug annotations in a ToxicoSet object
#'
#' @examples
#' data(TGGATESsmall)
#' drugInfo <- drugInfo(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#'
#' @return a \code{data.frame} with the drug annotations
#'
setGeneric("drugInfo", function(object) standardGeneric("drugInfo"))
#' @describeIn ToxicoSet Returns the annotations for all the drugs tested in
#'   the ToxicoSet
#' @export
setMethod("drugInfo", signature("ToxicoSet"), function(object) {
  object@drug
})


#' drugInfo<- Setter method
#'
#' Set the drug annotations in a ToxicoSet object
#'
#' @examples
#' data(TGGATESsmall)
#' drugInfo(TGGATESsmall) <- drugInfo(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object.
#' @param value A \code{data.frame} of replacement values.
#'
#' @return Updated \code{ToxicoSet}
#'
setGeneric("drugInfo<-", function(object, value) standardGeneric("drugInfo<-"))
#' @describeIn ToxicoSet Update the drug annotations
#' @export
setReplaceMethod("drugInfo", signature = signature(object = "ToxicoSet",value = "data.frame"), function(object, value){
  object@drug <- value
  object
})

#' phenoInfo Getter
#'
#' Get the phenotype annotations for cell lines with the specificed molecular
#'   data type
#'
#' @examples
#' data(TGGATESsmall)
#' phenoInfo <- phenoInfo(TGGATESsmall, mDataType="rna")
#'
#' @param object A \code{ToxicoSet} object
#' @param mDataType \code{character} A string specifying the type of molecular
#'   data to list the phenotype information for.
#'
#' @return a \code{Dframe} with the experiment info
#'
#' @describeIn ToxicoSet Return the experiment info from the given type of
#'   molecular data in ToxicoSet
#'
#' @importFrom CoreGx phenoInfo
#'
#' @export
setMethod("phenoInfo",
          signature("ToxicoSet", "character"),
          function(object, mDataType)
{
  callNextMethod(object, mDataType)
})

#' phenoInfo<- Setter
#'
#' Set the phenotype annotations for cell lines with the selected molecular
#'   data type.
#'
#' @examples
#' data(TGGATESsmall)
#' phenoInfo(TGGATESsmall, mDataType="rna") <-
#'   phenoInfo(TGGATESsmall, mDataType="rna")
#'
#' @param object A \code{ToxicoSet} object.
#' @param mDataType A \code{character} with the type of molecular data to return/update
#' @param value A \code{data.frame}, \code{DataFrame} or \code{DFrame} of
#'   replacement values.
#'
#' @return The updated \code{ToxicoSet}
#'
#' @describeIn ToxicoSet Update the the given type of molecular data experiment
#'   info in the ToxicoSet
#'
#' @importMethodsFrom CoreGx phenoInfo<-
#' @importFrom methods callNextMethod
#' @export
setReplaceMethod("phenoInfo",
                 signature = signature(object="ToxicoSet",
                                       mDataType ="character",
                                       value="data.frame"),
                 function(object, mDataType, value)
  {
  callNextMethod(object, mDataType, value)
})

#' molecularProfiles Getter
#'
#' Get the molecular profile data associated with the specific molecular data
#'
#' @examples
#' data(TGGATESsmall)
#' TGGATES_mProf <- molecularProfiles(TGGATESsmall, "rna")[seq_len(10),]
#'
#' @param object A \code{ToxicoSet} object.
#' @param mDataType \code{character} A string specifying the type of molecular
#'    data to list the phenotype information for.
#' @param assay A \code{character} Name of the desired assay; if excluded
#' defaults to first assay in the SummarizedExperiment for the given mDataType
#'
#' @describeIn ToxicoSet Return the given type of molecular data from the
#'   ToxicoSet
#'
#' @importFrom CoreGx molecularProfiles
#' @importFrom methods callNextMethod
#' @export
setMethod("molecularProfiles", signature("ToxicoSet"), function(object, mDataType, assay){
  callNextMethod(object, mDataType, assay)
})

#' molecularProfiles<- Setter
#'
#' Set the molecular profile data assciated with the specificied molecular
#'   data type
#'
#' @examples
#' molecularProfiles(TGGATESsmall, "rna") <-
#'   molecularProfiles(TGGATESsmall, "rna")
#'
#' @param object A \code{ToxicoSet} object.
#' @param value A \code{matrix} of replacement values.
#' @param assay \code{character} Name or index of the assay data to return
#' @param mDataType \code{character} A string specifying the type of molecular
#'   data to list the phenotype information for.
#'
#' @return Updated \code{ToxicoSet}
#'
#' @importMethodsFrom CoreGx molecularProfiles<-
#'
#' @describeIn ToxicoSet Update the given type of molecular data from the
#'    ToxicoSet
#'
#' @importFrom methods callNextMethod
#' @importFrom CoreGx molecularProfiles<-
#' @export
setReplaceMethod("molecularProfiles",
                 signature = signature(object="ToxicoSet",
                                       mDataType ="character",
                                       assay="character",
                                       value="matrix"),
                 function(object, mDataType, assay, value){
                   callNextMethod(object, mDataType, assay, value)
                 })
#' @describeIn ToxicoSet Update the given type of molecular data from the
#'   ToxicoSet
#' @export
setReplaceMethod("molecularProfiles",
                 signature = signature(object="ToxicoSet",
                                       mDataType ="character",
                                       assay="missing",
                                       value="matrix"),
                 function(object, mDataType, assay, value){
  callNextMethod(object, mDataType, assay, value)
})

#' molecularProfilesSlot Getter
#'
#' @describeIn ToxicoSet Get contents of molecularProfiles slot
#'
#' @examples
#' data(TGGATESsmall)
#' molecularProfilesSlot(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} from which to return a list of all availble
#'   SummarizedExperiment objects
#'
#' @return A \code{list} containing the molecularProfiles from a tSet
#'
#' @importFrom CoreGx molecularProfilesSlot
#' @importFrom methods callNextMethod
#' @export
setMethod("molecularProfilesSlot", signature("ToxicoSet"), function(object) {
  callNextMethod(object)
})

#' molecularProfilesSlot<- Setter
#'
#' @describeIn ToxicoSet Update the molecular profiles slot of a ToxicoSet and
#'    returns the updated copy
#'
#' @examples
#' data(TGGATESsmall)
#' molecularProfilesSlot(TGGATESsmall) <- molecularProfilesSlot(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object for which values will be replaced
#' @param value A \code{list} containing molecular profiles as SummarizedExperiments
#'
#' @return A copy of the \code{ToxicoSet} with the molecularProfiles slot updated
#'
#' @importFrom CoreGx molecularProfilesSlot<-
#' @importFrom methods callNextMethod
#' @export
setReplaceMethod("molecularProfilesSlot", signature("ToxicoSet"),
                 function(object, value) {
                   callNextMethod(object, value)
                 })


#' Getter for featureInfo method
#'
#' @examples
#' data(TGGATESsmall)
#' featureInfo <- featureInfo(TGGATESsmall, "rna")[seq_len(10),]
#'
#' @param object A \code{ToxicoSet} object
#' @param mDataType \code{character} A string specifying the type of molecular
#'   data to list the phenotype information for.
#'
#' @describeIn ToxicoSet Return the feature info for the given molecular data
#'
#' @importFrom CoreGx featureInfo
#' @importFrom methods callNextMethod
#' @export
setMethod("featureInfo",
          signature("ToxicoSet", "character"),
          function(object, mDataType)
{
  callNextMethod(object, mDataType)
})

#' featureInfo<- Setter
#'
#' Set the feature annotations for a specficied molecular data type
#'
#' @examples
#' data(TGGATESsmall)
#' featureInfo(TGGATESsmall, "rna") <- featureInfo(TGGATESsmall, "rna")
#'
#' @param object A \code{ToxicoSet} object
#' @param value A \code{data.frame} of replacement values
#' @param mDataType \code{character} A string specifying the type of molecular
#'   data
#'
#' @return Updated \code{ToxicoSet}
#'
#' @describeIn ToxicoSet Replace the gene info for the molecular data
#'
#' @importMethodsFrom CoreGx featureInfo<-
#' @importFrom methods callNextMethod
#' @export
setReplaceMethod("featureInfo",
                 signature = signature(object="ToxicoSet",
                                       mDataType ="character",
                                       value="data.frame"),
                 function(object, mDataType, value)
{
  callNextMethod(object, mDataType, value)
})

##TODO:: Migrate this to CoreGx
#' sensitivityRaw Generic
#'
#' @examples
#' data(TGGATESsmall)
#' sensitivityRaw(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} to extract the raw sensitivity data from
#' @param ... A \code{list} to allow new parameters in specific methods
#'
#' @return A \code{array} containing the raw sensitivity data
#'
#' @export
setGeneric("sensitivityRaw", function(object, ...) standardGeneric("sensitivityRaw"))
#' @describeIn ToxicoSet Retrive the raw dose and viability data from an tSet
#' @inheritParams sensitivityRaw
#' @export
setMethod("sensitivityRaw", signature("ToxicoSet"), function(object) {
  object@sensitivity$raw
})

##TODO:: Migrate this to CoreGx
#' sensitivityRaw<- Replacement Generic
#'
#' @examples
#' data(TGATESsmall)
#' sensitivityRaw(TGGATESsmall) <- sensitivityRaw(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} to extract the raw sensitivity data from
#' @param ... A \code{list} to allow new parameters in specific methods
#' @param value A \code{array} containing the raw dose and viability data for
#'   the tSet
#'
#' @return A copy of the \code{ToxicoSet} containing the updated sensitivty data
#'
#' @export
setGeneric("sensitivityRaw<-", function(object, ..., value) standardGeneric("sensitivityRaw<-"))
#' @describeIn ToxicoSet Set the raw dose and viability data for a tSet and return
#'   and updated copty
#' @inheritParams sensitivityRaw<-
#' @export
setReplaceMethod("sensitivityRaw", signature("ToxicoSet", "array"),
                 function(object, value) {
                   object@sensitivity$raw <- value
                   object
                 })

##TODO:: Migrate this to CoreGx
#' sensitivitySlot Generic
#'
#' @param object A \code{ToxicoSet} to extract the raw sensitivity data from
#' @param ... Allow new parameters to be defined for this generic
#'
#' @return A \code{list} of the sensitivity slot contents
#'
#' @export
setGeneric("sensitivitySlot", function(object, ...) standardGeneric("sensitivitySlot"))

#' sensitivitySlot Getter
#'
#' @describeIn ToxicoSet Retrieves the contents of the sensitivity slot
#'
#' @examples
#' data(TGGATESsmall)
#' sensitivitySlot(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} to extract the raw sensitivity data from
#'
#'
#' @export
setMethod("sensitivitySlot", signature("ToxicoSet"), function(object) {
  object@sensitivity
})

##TODO:: Migrate this to CoreGx
#' sensitivitySlot<- Replacement Generic
#'
#' @param object A \code{ToxicoSet} to extract the raw sensitivity data from
#' @param ... Allow new parameters to be defined for this generic
#' @param value A \code{list} of new sensitivity slot data for the tSet
#'
#' @return A copy of the \code{ToxicoSet} containing the updated sensitivty slot
#'
#' @export
setGeneric("sensitivitySlot<-", function(object, ..., value) standardGeneric("sensitivitySlot<-"))

#' sensitivity Slot Setter
#'
#' @describeIn ToxicoSet Set the raw dose and viability data for an tSet and return
#'   and updated copy
#'
#' @examples
#' data(TGGATESsmall)
#' sensitivitySlot(TGGATESsmall) <- sensitivitySlot(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} to extract the raw sensitivity data from
#' @param ... Allow new parameters to be defined for this generic
#' @param value A \code{list} of new sensitivity slot data for the tSet
#'
#' @export
setReplaceMethod("sensitivitySlot", signature(object="ToxicoSet", value="list"),
                 function(object, ..., value) {
                   ##TODO:: Implement error handinlg for this slot
                   object@sensitivity <- value
                   object
                 })

#' sensitivityInfo Getter
#'
#' Get the annotations for the sensitivity experiments in the ToxicoSet
#'
#' @examples
#' data(TGGATESsmall)
#' sensInf<- sensitivityInfo(TGGATESsmall)[seq_len(10),]
#'
#' @param object A \code{ToxicoSet} object
#'
#' @return a \code{data.frame} with the experiment info
#'
#' @describeIn ToxicoSet Return the drug dose sensitivity experiment info
#'
#' @importFrom CoreGx sensitivityInfo
#' @importFrom methods callNextMethod
#' @export
setMethod("sensitivityInfo",
          "ToxicoSet",
          function(object)
{
  callNextMethod(object)
})

#' sensitivityInfo<- Setter
#'
#' Set the annotations for sensitivity experiments in this ToxicSet
#'
#' @examples
#' data(TGGATESsmall)
#' sensitivityInfo(TGGATESsmall) <- sensitivityInfo(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#' @param value A \code{data.frame} of replacement values
#'
#' @describeIn ToxicoSet Update the sensitivity experiment info
#'
#' @return Updated \code{ToxicoSet}
#'
#' @importMethodsFrom CoreGx sensitivityInfo<-
#' @export
setReplaceMethod("sensitivityInfo",
                 signature = signature(object="ToxicoSet",
                                       value="data.frame"),
                 function(object, value)
{
  callNextMethod(object, value)
})


#' sensitivityProfiles Getter
#'
#' Get the data for sensitivty experiments on cell lines in a ToxicoSet
#'
#' @examples
#' data(TGGATESsmall)
#' sensProf <- sensitivityProfiles(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#'
#' @return a \code{data.frame} with the experiment info
#'
#' @describeIn ToxicoSet Return the phenotypic data for the drug dose
#'   sensitivity
#'
#' @importFrom CoreGx sensitivityProfiles
#' @export
setMethod("sensitivityProfiles",
          "ToxicoSet",
          function(object)
{
  callNextMethod(object)
})


#' sensitivityProfiles<-
#'
#' Set the data for sensitivity experiments on cell lines in a ToxicoSet
#'
#' @examples
#' sensitivityProfiles(TGGATESsmall) <- sensitivityProfiles(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#' @param value A \code{data.frame} of replacement values
#'
#' @return Updated \code{ToxicoSet}
#'
#' @importFrom CoreGx sensitivityProfiles<-
#'
#' @describeIn ToxicoSet Update the phenotypic data for the drug dose
#'   sensitivity
#' @export
setReplaceMethod("sensitivityProfiles",
                 signature = signature(object="ToxicoSet",
                                       value="data.frame"),
                 function(object, value)
{
  callNextMethod(object, value)
})
#' @describeIn ToxicoSet Update the phenotypic data for the drug dose
#'   sensitivity
#'
## TODO:: Find out how to document overloaded methods (to include multiple parameter types)
#' @export
setReplaceMethod("sensitivityProfiles",
                 signature = signature(object="ToxicoSet",
                                       value="matrix"),
                 function(object, value)
{
  callNextMethod(object, value)
})

#' sensitivityMeasures Getter
#'
#' Get the avilable measurments for sensitivity experiments in a ToxicoSet
#'
#' @examples
#' sensitivityMeasures(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#'
#' @return A \code{character} vector of all the available sensitivity measures
#'
#' @describeIn ToxicoSet Returns the available sensitivity profile
#'   summaries, for example, whether there are IC50 values available
#'
#' @importFrom CoreGx sensitivityMeasures
#' @importFrom methods callNextMethod
#' @export
setMethod("sensitivityMeasures",
          "ToxicoSet",
          function(object)
{
  callNextMethod(object)
})

##TODO:: Abstract this method to CoreGx
#' drugNames Generic
#'
#' A generic for the drugNames method
#'
#' @examples
#' data(TGGATESsmall)
#' drugName <- drugNames(TGGATESsmall)[seq_len(10)]
#'
#' @param object A \code{ToxicoSet} object from which to retrieve the included
#'   drug names
#'
#' @return A vector of the drug names used in the ToxicoSet
setGeneric("drugNames", function(object) standardGeneric("drugNames"))
#' @describeIn ToxicoSet Return the names of the drugs used in the ToxicoSet
#' @export
setMethod(drugNames,
          "ToxicoSet",
          function(object)
{
  rownames(drugInfo(object))
})

##TODO:: Abstract this method to CoreGx
#' drugNames<- Generic
#'
#' A generic for the drugNames replacement method
#'
#' @examples
#' data(TGGATESsmall)
#' drugNames(TGGATESsmall) <- drugNames(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object to modify
#' @param value A \code{character} vector of replacement drug names
#'
#' @return Updated \code{ToxicoSet}
setGeneric("drugNames<-", function(object, value) standardGeneric("drugNames<-"))
#' @describeIn ToxicoSet Update the drug names used in the dataset
#' @export
setReplaceMethod("drugNames",
                 signature = signature(object="ToxicoSet",
                                       value="character"),
                 function(object, value)
{
  object <- updateDrugId(object, value)
  return(object)
})

#' cellNames Getter
#'
#' Get names of cell lines in a ToxicoSet object
#'
#' @examples
#' cellNames(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#'
#' @return A vector of the cell names used in the ToxicoSet
#'
#' @describeIn ToxicoSet Return the cell names used in the dataset
#'
#' @importFrom CoreGx cellNames
#' @importFrom methods callNextMethod
#' @export
setMethod("cellNames",
          signature("ToxicoSet"),
          function(object)
{
  callNextMethod(object)
})

#' cellNames<- Setter
#'
#' Set the cell line names in a TocicoSet object
#'
#' @examples
#' data(TGGATESsmall)
#' cellNames(TGGATESsmall) <- cellNames(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object to modify
#' @param value A \code{character} of replacement cell line names
#'
#' @return Updated \code{ToxicoSet}
#'
#' @describeIn ToxicoSet Update the cell names used in the dataset
#'
#' @importMethodsFrom CoreGx cellNames<-
#' @export
setReplaceMethod("cellNames",
                 signature = signature(object="ToxicoSet",
                                       value="character"),
                 function(object, value)
{
  callNextMethod(object, value)
})

#' fNames Getter
#'
#' Get the feature names in a ToxicoSet object for the specified molecular data
#'   type
#'
#' @examples
#' fNames(TGGATESsmall, "rna")[seq_len(10)]
#'
#' @param object A \code{ToxicoSet} object
#' @param mDataType \code{character} A string specifying the type of molecular
#'   data to list the phenotype information for
#'
#' @return A \code{character} vector of the feature names
#'
#' @describeIn ToxicoSet Return the feature names used in the dataset
#'
#' @importFrom CoreGx fNames
#' @importFrom methods callNextMethod
#' @export
setMethod("fNames",
          signature("ToxicoSet", "character"),
          function(object, mDataType)
{
  callNextMethod(object, mDataType)
})

#' fNames<- Setter
#'
#' Set the feature names in a ToxicoSet object for the specified molecular
#'   data type
#'
#'@examples
#' data(TGGATESsmall)
#' cellNames(TGGATESsmall) <- cellNames(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#' @param mDataType \code{character} A string specifying the type of molecular
#'   data to list the phenotype information for.
#' @param value A \code{character} vector of the new feature names
#'
#' @return Updated \code{ToxicoSet}
#'
#' @describeIn ToxicoSet Update the feature names used in the dataset
#'
#' @importFrom CoreGx fNames<-
#' @importFrom methods callNextMethod
#' @export
setReplaceMethod("fNames",
                 signature("ToxicoSet", "character"),
                 function(object, mDataType, value)
{
  callNextMethod(object, mDataType, value)
})

#' dateCreated Getter
#'
#' Get the date a ToxicoSet object was created
#'
#' @examples
#' dateCreated(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#'
#' @return The date the ToxicoSet was created
#'
#' @describeIn ToxicoSet Return the date the ToxicoSet was created
#'
#' @importFrom CoreGx dateCreated
#' @importFrom methods callNextMethod
#' @export
setMethod("dateCreated",
          signature("ToxicoSet"),
          function(object)
{
  callNextMethod(object)
})

##TODO:: Export this to CoreGx

#' datasetType Generic
#'
#' A generic for retrieving the dataset type of an tSet object
#'
#' @param object A \code{ToxicoSet} from which to retrieve the dataset type
#' @param ... A \code{list} containing fall through arguments; this allows
#'   addition of new parameters to methods for this generic
#'
#' @return A \code{character} vector containing the dataset type
#'
#' @export
setGeneric("datasetType", function(object, ...) standardGeneric("datasetType"))

#' datasetType Getter
#'
#' @examples
#' data(TGGATESsmall)
#' datasetType(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} from which to retrieve the dataset type
#' @param ... A \code{list} containing fall through arguments; this allows
#'   addition of new parameters to methods for this generic
#'
#' @describeIn ToxicoSet Update the dataset type of an tSet and return a copy of
#'     the updated object
#' @export
setMethod("datasetType", signature("ToxicoSet"), function(object) {
  ##TODO:: Add error handling to this function
  object@datasetType
})


##TODO:: Export this to CoreGx
#' datasetType<- Replacement Generic
#'
#' A generic for updating the dataset type of a ToxicoSet object
#'
#' @examples
#' data(TGGATESsmall)
#' datasetType(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} from which to retrieve the dataset type
#' @param value A \code{character} vector containing the dataset type
#'
#' @return A \code{ToxicoSet} with the datasetType slot updated
#'
#' @export
setGeneric("datasetType<-",  function(object, value) standardGeneric("datasetType<-"))
#' @inheritParams datasetType<-
#' @describeIn ToxicoSet Update the dataset type of an tSet and return a copy of
#'     the updated object
#' @export
setReplaceMethod("datasetType", signature("ToxicoSet"), function(object, value) {
  ##TODO:: Add error handling to this function
  object@datasetType <- value
  object
})


#' pertNumber Getter
#'
#' Get an array of the number of pertubration experiments per drug and cell
#'   line in a ToxicoSet object
#'
#' @examples
#' pertNumber(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#'
#' @return A 3D \code{array} with the number of perturbation experiments per
#'   radiation type and cell line, and data type
#'
#' @describeIn ToxicoSet Return the summary of available perturbation
#'   experiments
#'
#' @importMethodsFrom CoreGx pertNumber
#' @importFrom methods callNextMethod
#' @export
setMethod("pertNumber",
          signature("ToxicoSet"),
          function(object)
{
  callNextMethod(object)
})


#' sensNumber Getter
#'
#' Get the number of sensitivity experiments per drug and cell line in a
#'   ToxicoSet
#'
#' @examples
#' sensNumber(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#'
#' @return A \code{data.frame} with the number of sensitivity experiments per
#'   drug and cell line
#'
#' @describeIn ToxicoSet Return the summary of available sensitivity
#'   experiments
#'
#' @importFrom CoreGx sensNumber
#' @importFrom methods callNextMethod
#' @export
setMethod("sensNumber",
          "ToxicoSet",
          function(object)
{
  callNextMethod(object)
})

#' pertNumber<- Setter
#'
#' Set the number of perturbation experiments per drug and cell line and
#'   molecular data type in a ToxicoSet object
#'
#' @examples
#' pertNumber(TGGATESsmall) <- pertNumber(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object to modify
#' @param value An \code{array} of replacement values
#'
#' @return The updated \code{ToxicoSet}
#'
#' @describeIn ToxicoSet Update the summary of available perturbation
#'   experiments
#'
#' @importMethodsFrom CoreGx pertNumber<-
#' @export
setReplaceMethod('pertNumber',
                 signature = signature(object="ToxicoSet",
                                       value="array"),
                 function(object, value)
{
  callNextMethod(object, value)
})

#' sensNumber<- Setter
#'
#' Set the number of sensitivity experiments per drug and cell line in a
#'   ToxicoSet object
#'
#' @examples
#' sensNumber(TGGATESsmall) <- sensNumber(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object to modify
#' @param value A \code{matrix} of replacement values
#'
#' @return The updated \code{ToxicoSet}
#'
#' @describeIn ToxicoSet Update the summary of available sensitivity
#'   experiments
#'
#' @importMethodsFrom CoreGx sensNumber<-
#' @importFrom methods callNextMethod
#' @export
setReplaceMethod('sensNumber',
                 signature = signature(object="ToxicoSet",
                                       value="matrix"),
                 function(object, value)
{
  callNextMethod(object, value)
})

##TODO:: Implement a limited show method in CoreGx which can be extended
##  here
#' Show a ToxicoSet
#'
#' @param object A \code{ToxicoSet} object to print a summary for
#'
#' @examples
#' TGGATESsmall
#'
#' @return Prints the ToxicoSet object to the output stream, and returns
#'   invisible NULL.
#'
#' @export
setMethod("show", signature=signature(object="ToxicoSet"),
          function(object) {
            cat("Name: ", name(object), "\n")
            cat("Date Created: ", dateCreated(object), "\n")
            cat("Number of cell lines: ", nrow(cellInfo(object)), "\n")
            cat("Number of drugs: ", nrow(drugInfo(object)), "\n")
            if("dna" %in% names(molecularProfilesSlot(object))){cat("DNA: \n");
              cat("\tDim: ", dim(molecularProfiles(object, mDataType="dna")),
                "\n")}
            if("rna" %in% names(molecularProfilesSlot(object))){cat("RNA: \n");
              cat("\tDim: ", dim(molecularProfiles(object, mDataType="rna")),
                "\n")}
            if("rnaseq" %in% names(molecularProfilesSlot(object))){cat("RNASeq: \n");
              cat("\tDim: ", dim(molecularProfiles(object, mDataType="rnaseq")),
                "\n")}
            if("snp" %in% names(molecularProfilesSlot(object))){cat("SNP: \n");
              cat("\tDim: ", dim(molecularProfiles(object, mDataType="snp")),
                "\n")}
            if("cnv" %in% names(molecularProfilesSlot(object))){cat("CNV: \n");
              cat("\tDim: ", dim(molecularProfiles(object, mDataType="cnv")),
                "\n")}
            cat("Drug pertubation: \n")
            cat("\tPlease look at pertNumber(tSet) to determine number of
              experiments for each drug-cell combination.\n")
            cat("Drug sensitivity: \n")
            cat("\tNumber of Experiments: ",nrow(sensitivityInfo(object)),"\n")
            cat("\tPlease look at sensNumber(tSet) to determine number of
              experiments for each drug-cell combination.\n")
          })

#' mDataNames
#'
#' Returns the names of the molecular data types available in a ToxicoSet
#'   object
#'
#' @examples
#' mDataNames(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#'
#' @return Vector of names of the molecular data types
#'
#' @importFrom CoreGx mDataNames cellInfo<-
#' @export
setMethod(
  "mDataNames",
  signature("ToxicoSet"),
  function(object)
{
    callNextMethod(object)
})

##TODO:: Export to CoreGx
##FIXME:: How do I import generics from BiocGenerics?
#' annotation Slot Getter
#'
#' @param object A \code{ToxicoSet}
#' @param ... A \code{list} to allow definition of new parameters on this generic
#'
#' @return A \code{list} of named annotaiton
#'
#' @examples
#' data(TGGATESsmall)
#' annotation(TGGATESsmall)
#'
#' @export
setGeneric("annotation", function(object, ...) standardGeneric("annotation"))
#' @describeIn ToxicoSet Retrieve the annotations slot form an tSet
#' @inheritParams annotation<-
#' @export
setMethod('annotation', signature("ToxicoSet"), function(object) {
  object@annotation
})

##TODO:: Export to CoreGx
##FIXME:: How do I import generics from BiocGenerics?
#' annotation<- Slot Setter
#'
#' @param object A \code{ToxicoSet}
#' @param ... A \code{list} to allow definition of new parameters on this generic
#' @param value A \code{list} of annotations to add to the annotatiosn slot of
#'   an tSet
#'
#' @return A copy of the \code{ToxicoSet} with the updated annotation slot
#'
#' @examples
#' data(TGGATESsmall)
#' annotation(TGGATESsmall) <- annotation(TGGATESsmall)
#'
#' @export
setGeneric("annotation<-", function(object, ..., value) standardGeneric("annotation<-"))
#' @describeIn ToxicoSet Update the annotation slot of a tSet
#' @inheritParams annotation<-
#' @export
setReplaceMethod("annotation", signature("ToxicoSet", "list"), function(object, value) {
  object@annotation <- value
  object
})

##TODO:: Export to CoreGx
#' curation Slot Getter
#'
#' @param object A \code{ToxicoSet}
#' @param ... A \code{list} to allow definition of new parameters on this generic
#'
#' @return A \code{list} of unique cell and tissue identifiers to check validity
#'   of an tSet
#'
#' @examples
#' data(TGGATESsmall)
#' curation(TGGATESsmall)
#'
#' @export
setGeneric("curation", function(object, ...) standardGeneric("curation"))
#' @describeIn ToxicoSet Retrieve the curation slot form an tSet
#' @inheritParams curation
#' @export
setMethod('curation', signature("ToxicoSet"), function(object) {
  object@curation
})

##TODO:: Export to CoreGx
##FIXME:: How do I import generics from BiocGenerics?
#' curation<- Slot Setter
#'
#' @param object A \code{ToxicoSet}
#' @param ... A \code{list} to allow definition of new parameters on this generic
#' @param value A \code{list} of curations for the cell and tissues types in the
#'   tSet object
#'
#' @return A copy of the \code{ToxicoSet} with the updated curation slot
#'
#' @examples
#' data(TGGATESsmall)
#' curation(TGGATESsmall) <- curation(TGGATESsmall)
#'
#' @export
setGeneric("curation<-", function(object, ..., value) standardGeneric("curation<-"))
#' @describeIn ToxicoSet Update the annotation slot of a tSet
#' @inheritParams annotation<-
#' @export
setReplaceMethod("curation", signature("ToxicoSet", "list"), function(object, value) {
  object@curation <- value
  object
})

#'`[`
#'
#' @examples
#' tSet <- TGGATESsmall[cellNames(TGGATESsmall), drugNames(TGGATESsmall)[seq_len(3)]]
#'
#'@param x tSet
#'@param i Cell lines to keep in tSet
#'@param j Drugs to keep in tSet
#'@param ... further arguments
#'@param drop A boolean flag of whether to drop single dimensions or not
#'@return Returns the subsetted tSet
#'@export
setMethod(`[`, "ToxicoSet", function(x, i, j, ..., drop = FALSE){
  if(is.character(i) && is.character(j)) {
    return(subsetTo(x, cells=i, drugs=j,  molecular.data.cells=i))
  }
  else if(is.numeric(i) && is.numeric(j) &&
           (as.integer(i)==i) && (as.integer(j)==j)) {
    return(subsetTo(x, cells=cellNames(x)[i], drugs=drugNames(x)[j],
                    molecular.data.cells=cellNames(x)[i]))
  }
})

#' Get the dimensions of a ToxicoSet
#'
#' @examples
#' data(TGGATESsmall)
#' dim(TGGATESsmall)
#'
#' @param x ToxicoSet
#' @return A named vector with the number of Cells and Drugs in the ToxicoSet
#' @export
setMethod(
  "dim",
  signature("ToxicoSet"),
  function(x)
{
  return(c(Cells=length(cellNames(x)), Drugs=length(drugNames(x))))
})

#### subsetTo ####

## FIXED? TODO:: Subset function breaks if it doesnt find cell line in sensitivity info
#' A function to subset a ToxicoSet to data containing only specified drugs, cells and genes
#'
#' This is the prefered method of subsetting a ToxicoSet. This function allows
#' abstraction of the data to the level of biologically relevant objects: drugs
#' and cells. The function will automatically go through all of the
#' combined data in the ToxicoSet and ensure only the requested radiations
#' and cell lines are found in any of the slots. This allows quickly picking out
#' all the experiments for a radiation or cell of interest, as well removes the need
#' to keep track of all the metadata conventions between different datasets.
#'
#' @examples
#' TGGATESDrugNames  <- drugNames(TGGATESsmall)
#' TGGATESCells <- cellNames(TGGATESsmall)
#' tSet <- subsetTo(TGGATESsmall,drugs = TGGATESDrugNames[1],
#'   cells = TGGATESCells[1], duration = "2")
#'
#' @param object A \code{ToxicoSet} to be subsetted
#' @param cell_lines A list or vector of cell names as used in the dataset to which
#'   the object will be subsetted. If left blank, then all cells will be left in
#'   the dataset.
#' @param drugs A list or vector of drug names as used in the dataset to which
#'   the object will be subsetted. If left blank, then all drugs will be left in
#'   the dataset.
#' @param features A list or vector of feature names as used in the dataset from
#'   which the object will be subsetted. If left blank that all features will
#'   be left in.
#' @param molecular.data.cells A list or vector of cell names to keep in the
#'   molecular data
#' @param duration A \code{list} or \code{vector} of the experimental durations
#'   to include in the subset as strings. Defaults to all durations if parameter
#'   is not specified.
#' @param ... Other arguments passed to other functions within the package
#'
#' @return A ToxicoSet with only the selected drugs and cells
#'
#' @importFrom CoreGx .unionList
#' @export
## TODO:: Include dose parmater to subset on
subsetTo <- function(object, cell_lines = NULL,
                     drugs=NULL,
                     molecular.data.cells=NULL,
                     duration=NULL,
                     features=NULL,
                     ...
                     ) {
  ## TODO:: Remove this or add it to the function parameters?
  drop = FALSE

  ####
  # PARSING ARGUMENTS
  ####
  adArgs = list(...)
  if ("exps" %in% names(adArgs)) {
    exps <- adArgs[["exps"]]
    if(is(exps, "data.frame")) {
      exps2 <- exps[[name(object)]]
      names(exps2) <- rownames(exps)
      exps <- exps2
    } else{
      exps <- exps[[name(object)]]
    }
  }else {
    exps <- NULL
  }

  if ("dose" %in% names(adArgs)) {
    ## TODO:: Add subsetting on dose
    stop("Due to the structure of tSets, subsetting on dose can only be done on
      specific slots - not on the entire tSet")
  }

  ## MISSING VALUE HANDLING FOR PARAMETERS
  # Get named list of default values for missing parameters
  argDefaultList <-
    paramMissingHandler(funName = "subsetTo", tSet = object,
                        drugs = drugs, cell_lines = cell_lines, features = features,
                        duration = duration)
  # Assign any missing parameter default values to function environment
  if (length(argDefaultList) > 0) {
    for (idx in seq_along(argDefaultList)) {
      assign(names(argDefaultList)[idx], argDefaultList[[idx]])
    }
  }

  # ERROR HANDLING FOR PARAMETERS
  paramErrorChecker(funName = "subsetTo", tSet = object,
                    cell_lines = cell_lines,
                    drugs = drugs, features = features,
                    duration = duration)

  ##TODO:: Add a value to tSet which indicates the experimental design!
  ##FIXME:: Don't hard code object names!
  if (name(object) == "drugMatrix") {
    if (!('DMSO' %in% drugs)) {
      drugs <- c(drugs, 'DMSO')
    }
  }

  ######
  # SUBSETTING MOLECULAR PROFILES SLOT
  ######
  ### TODO:: implement strict subsetting at this level!!!!

  ### the function missing does not work as expected in the context below, because the arguments are passed to the anonymous
  ### function in lapply, so it does not recognize them as missing
  molecularProfilesSlot(object) <-
    lapply(molecularProfilesSlot(object),
      function(SE, cell_lines, drugs, molecular.data.cells, duration, features){

    if (!is.null(features)) {
      SE <- SE[which(rownames(SummarizedExperiment::rowData(SE)) %in% features), ]
    }

    ##FIXME:: Why is are all these if conditions being checked against length? Just use grepl?
    molecular.data.type <-
      ifelse(
        length(grep("rna", S4Vectors::metadata(SE)$annotation) > 0),
        "rna",
        S4Vectors::metadata(SE)$annotation
        )

    if (length(grep(molecular.data.type, names(molecular.data.cells))) > 0) {
      cell_lines <- molecular.data.cells[[molecular.data.type]]
    }
    column_indices <- NULL

    if (length(cell_lines) == 0 && length(drugs) == 0) {
      column_indices <- seq_len(ncol(SE))
    }
    if (length(cell_lines) == 0 && datasetType(object) == "sensitivity") {
      column_indices <- seq_len(ncol(SE))
    }

    # Selecting indices which match the cells argument
    cell_line_index <- NULL
    if (length(cell_lines) != 0) {
      if (!all(cell_lines %in% cellNames(object))) {
        stop("Some of the cell names passed to function did not match to names
          in the ToxicoSet. Please ensure you are using cell names as
          returned by the cellNames function")
      }
      cell_line_index <- which(SummarizedExperiment::colData(SE)[["cellid"]] %in% cell_lines)
    }

    # Selecting indexes which match drugs arguement
    drugs_index <- NULL
    if (datasetType(object) == "perturbation" || datasetType(object) == "both") {
      if (length(drugs) != 0) {
        if (!all(drugs %in% drugNames(object))){
          stop("Some of the drug names passed to function did not match to names in the ToxicoSet Please ensure you are using drug names as returned by the drugNames function")
        }
        drugs_index <- which(SummarizedExperiment::colData(SE)[["drugid"]] %in% drugs)
      }
    }

    if (length(drugs_index) != 0 && length(cell_line_index) != 0) {
      if (length(intersect(drugs_index, cell_line_index)) == 0) {
        stop("This Drug - Cell Line combination was not tested together.")
      }
      column_indices <- intersect(drugs_index, cell_line_index)
    } else {
      if (length(drugs_index) != 0) {
        column_indices <- drugs_index
      }
      if (length(cell_line_index) != 0) {
        column_indices <- cell_line_index
      }
    }

    # LOGIC TO SUBSET BASED ON DURATION
    ## TODO:: Determine if this works for other SE data types
    if (!is.null(duration)){
      if (all(!(duration %in% unique(SummarizedExperiment::colData(SE[, column_indices])$duration)))) {
        # Error when other parameters are passed in
        if ( !is.null(cell_lines) | !is.null(drugs) | !is.null(molecular.data.cells)) {
          stop(paste0(
            "There are no molecular profiles with duration of ",
            duration, " in the tSet with the selected parameters."
          ))
        } else { # Error when no other parameters are passed in
          stop(paste0(
            "There are no molecular profiles with duration of ",
            duration, " in the tSet."
          ))
        }
      }
      duration_indices <- which(SummarizedExperiment::colData(SE)$duration %in% duration)
      column_indices <- intersect(column_indices, duration_indices)
    }

    row_indices <- seq_len(nrow(SummarizedExperiment::assay(SE, 1)))

    # Final SE
    SE <- SE[row_indices, column_indices]
    return(SE)

  }, cell_lines = cell_lines,
      drugs = drugs,
      molecular.data.cells = molecular.data.cells,
      duration = duration,
      features = features)


  ######
  # SUBSET SENSITIVITY SLOT
  ######
  # Logic if any "..." arguments are passed to subsetTo
  if ((datasetType(object) == "sensitivity" | datasetType(object) == "both") & length(exps) != 0) {
    sensitivityInfo(object) <- sensitivityInfo(object)[exps, , drop=drop]
    rownames(sensitivityInfo(object)) <- names(exps)
    if (length(sensitivityRaw(object)) > 0) {
      sensitivityRaw(object) <- sensitivityRaw(object)[exps, , , drop=drop]
      dimnames(sensitivityRaw(object))[[1]] <- names(exps)
    }
    sensitivityProfiles(object) <- sensitivityProfiles(object)[exps, , drop=drop]
    rownames(sensitivityProfiles(object)) <- names(exps)

    sensNumber(object) <- .summarizeSensitivityNumbers(object)
  }
  # Logic if drug or cell parameters are passed to subsetTo
  else if (
    (datasetType(object) == "sensitivity" | datasetType(object) == "both") &
    (length(drugs) != 0 | length(cell_lines) != 0 | !is.null(duration) )
  ) {

    drugs_index <- which(sensitivityInfo(object)[, "drugid"] %in% drugs)
    cell_line_index <- which(sensitivityInfo(object)[,"cellid"] %in% cell_lines)
    if (length(drugs_index) !=0 & length(cell_line_index) !=0 ) {
      if (length(intersect(drugs_index, cell_line_index)) == 0) {
        stop("This Drug - Cell Line combination was not tested together.")
      }
      row_indices <- intersect(drugs_index, cell_line_index)
    } else {
      if(length(drugs_index)!=0 & length(cell_lines)==0) {
        row_indices <- drugs_index
      } else {
        if(length(cell_line_index)!=0 & length(drugs)==0){
          row_indices <- cell_line_index
        } else {
          # Includes all rows if cell or drug arguments are absent
          row_indices <- seq_len(nrow(sensitivityInfo(object)))
        }
      }
    }
    # LOGIC TO SUBSET BASED ON DURATION
    if(!is.null(duration)){
      if(all(!(duration %in% unique(sensitivityInfo(object)[row_indices,]$duration_h)))) {
        # Error when other parameters are passed in
        if(!is.null(cell_lines) | !is.null(drugs) | !is.null(molecular.data.cells)) {
          stop(paste0(
            ## TODO:: Is sample the correct way to refer to one treatment/duration combination in TGx experiments?
            "There are no samples with duration of ",
            duration, " in the tSet with the selected parameters."
          ))
        } else { # Error when no other parameters are passed in
          stop(paste0(
            "There are no samples with duration of ",
            duration, " in the tSet"
          ))
        }
      }
      duration_indices <- which(sensitivityInfo(object)$duration_h %in% duration)
      row_indices <- intersect(row_indices, duration_indices)
    }
    sensItemNames <- names(sensitivitySlot(object))
    sensitivityVals <-
      lapply(sensItemNames, function(sensItemName, drop){
        if (sensItemName == "n") {
          sensItem <- sensitivitySlot(object)[[sensItemName]]
          if (!is.null(cell_lines)) {
            sensItem[which(rownames(sensItem) %in% cell_lines),
                     which(colnames(sensItem) %in% drugs), drop = drop]
          } else {
            sensItem[ , which(colnames(sensItem) %in% drugs), drop = drop]
          }
        } else {
          sensItem <- sensitivitySlot(object)[[sensItemName]]
          if (length(dim(sensItem)) == 3) {
            sensItem[row_indices, , , drop = drop]
          } else {
            sensItem[row_indices, , drop = drop]
          }
        }
      }, drop = drop)
    names(sensitivityVals) <- sensItemNames
    sensitivitySlot(object) <- sensitivityVals
  }

  #####
  # SUBSET DRUG SLOT
  #####
  if (length(drugs) == 0) {
    if (datasetType(object) == "sensitivity" | datasetType(object) == "both"){
      drugs <- unique(sensitivityInfo(object)[["drugid"]])
    }
    if(datasetType(object) == "perturbation" | datasetType(object) == "both"){
      drugs <- union(drugs, na.omit(.unionList(lapply(molecularProfilesSlot(object), function(SE){unique(SummarizedExperiment::colData(SE)[["drugid"]])}))))
    }
  }

  #####
  # SUBSET CELLS SLOT
  #####
  if (length(cell_lines) == 0) {
    cell_lines <- union(cell_lines, na.omit(.unionList(lapply(molecularProfilesSlot(object), function(SE){unique(SummarizedExperiment::colData(SE)[["cellid"]])}))))
    if (datasetType(object) == "sensitivity" | datasetType(object) == "both"){
      cell_lines <- union(cell_lines, sensitivityInfo(object)[["cellid"]])
    }
  }
  #####
  # ASSIGN SUBSETS BACK TO TOXICOSET OBJECT
  #####
  drugInfo(object) <- drugInfo(object)[drugs , , drop=drop]
  cellInfo(object) <- cellInfo(object)[cell_lines , , drop=drop]
  object@curation$drug <- object@curation$drug[drugs , , drop=drop]
  object@curation$cell <- object@curation$cell[cell_lines , , drop=drop]
  object@curation$tissue <- object@curation$tissue[cell_lines , , drop=drop]
  return(object)
}

#
# END SUBSET TO FUNCTION
#



### TODO:: Add updating of sensitivity Number tables
#' A function to update cell ids
#' @examples
#' data(TGGATESsmall)
#' updateCellId(TGGATESsmall, new.ids = cellNames(TGGATESsmall))
#'
#' @param tSet A \code{ToxicoSet} object to be updates
#' @param new.ids A \code{character} vector of ids to update with
#'
#' @return \code{none} Updates the cell ids in the ToxicoSet
#'
#' @importFrom CoreGx updateCellId
#' @keywords internal
#' @export
updateCellId <- function(tSet, new.ids = vector("character")){
  CoreGx::updateCellId(tSet, new.ids)
}

### TODO:: Add updating of sensitivity Number tables
#' A function to update drug ids
#' @examples
#' data(TGGATESsmall)
#' updateDrugId(TGGATESsmall, new.ids = drugNames(TGGATESsmall))
#'
#' @param tSet [object] A ToxicoSet object to be updates
#' @param new.ids [character] A character vector of ids to update with
#'
#' @return \code{none} Updates the drug ids in the ToxicoSet
#'
#' @keywords internal
#' @export
updateDrugId <- function(tSet, new.ids = vector("character")){

  if (length(new.ids)!= nrow(drugInfo(tSet))){
    stop("Wrong number of drug identifiers")
  }

  if(tSet@datasetType=="sensitivity" | tSet@datasetType=="both"){
    myx <- match(sensitivityInfo(tSet)[,"drugid"],rownames(drugInfo(tSet)))
    sensitivityInfo(tSet)[,"drugid"] <- new.ids[myx]

  }
  if(tSet@datasetType=="perturbation"|tSet@datasetType=="both"){
    tSet@molecularProfiles <- lapply(tSet@molecularProfiles, function(SE){

      myx <- match(SummarizedExperiment::colData(SE)[["drugid"]],rownames(drugInfo(tSet)))
      SummarizedExperiment::colData(SE)[["drugid"]]  <- new.ids[myx]
      return(SE)
    })
  }


  if(any(duplicated(new.ids))){
    warning("Duplicated ids passed to updateDrugId. Merging old ids into the
            same identifier")

    if(ncol(sensNumber(tSet)) > 0){
      sensMatch <- match(colnames(sensNumber(tSet)), rownames(drugInfo(tSet)))
    }
    if(dim(pertNumber(tSet))[[2]] > 0){
      pertMatch <- match(dimnames(pertNumber(tSet))[[2]], rownames(drugInfo(tSet)))
    }
    curMatch <- match(rownames(tSet@curation$drug),rownames(drugInfo(tSet)))

    duplId <- unique(new.ids[duplicated(new.ids)])
    for(id in duplId){

      if (ncol(sensNumber(tSet))>0){
        myx <- which(new.ids[sensMatch] == id)
        sensNumber(tSet)[,myx[1]] <- apply(sensNumber(tSet)[,myx], 1, sum)
        sensNumber(tSet) <- sensNumber(tSet)[,-myx[-1]]
        # sensMatch <- sensMatch[-myx[-1]]
      }
      if (dim(pertNumber(tSet))[[2]]>0){
        myx <- which(new.ids[pertMatch] == id)
        pertNumber(tSet)[,myx[1],] <- apply(pertNumber(tSet)[,myx,], c(1,3), sum)
        pertNumber(tSet) <- pertNumber(tSet)[,-myx[-1],]
        # pertMatch <- pertMatch[-myx[-1]]
      }

      myx <- which(new.ids[curMatch] == id)
      tSet@curation$drug[myx[1],] <- apply(tSet@curation$drug[myx,], 2, paste, collapse="///")
      tSet@curation$drug <- tSet@curation$drug[-myx[-1],]
      # curMatch <- curMatch[-myx[-1]]

      myx <- which(new.ids == id)
      drugInfo(tSet)[myx[1],] <- apply(drugInfo(tSet)[myx,], 2, paste, collapse="///")
      drugInfo(tSet) <- drugInfo(tSet)[-myx[-1],]
      new.ids <- new.ids[-myx[-1]]
      if(ncol(sensNumber(tSet))>0){
        sensMatch <- match(colnames(sensNumber(tSet)), rownames(drugInfo(tSet)))
      }
      if(dim(pertNumber(tSet))[[2]]>0){
        pertMatch <- match(dimnames(pertNumber(tSet))[[2]], rownames(drugInfo(tSet)))
      }
      curMatch <- match(rownames(tSet@curation$drug),rownames(drugInfo(tSet)))
    }
  } else {
    if (dim(pertNumber(tSet))[[2]]>0){
      pertMatch <- match(dimnames(pertNumber(tSet))[[2]], rownames(drugInfo(tSet)))
    }
    if (ncol(sensNumber(tSet))>0){
      sensMatch <- match(colnames(sensNumber(tSet)), rownames(drugInfo(tSet)))
    }
    curMatch <- match(rownames(tSet@curation$drug),rownames(drugInfo(tSet)))
  }

  if (dim(pertNumber(tSet))[[2]]>0){
    dimnames(pertNumber(tSet))[[2]] <- new.ids[pertMatch]
  }
  if (ncol(sensNumber(tSet))>0){
    colnames(sensNumber(tSet)) <- new.ids[sensMatch]
  }
  rownames(tSet@curation$drug) <- new.ids[curMatch]
  rownames(drugInfo(tSet)) <- new.ids


  return(tSet)
}

.summarizeSensitivityNumbers <- function(tSet) {

  if (tSet@datasetType != "sensitivity" && tSet@datasetType != "both") {
    stop ("Data type must be either sensitivity or both")
  }

  ## consider all drugs
  drugn <- rownames(tSet@drug)

  ## consider all cell lines
  celln <- rownames(tSet@cell)

  sensitivity.info <- matrix(0, nrow=length(celln), ncol=length(drugn), dimnames=list(celln, drugn))
  drugids <- tSet@sensitivity$info[ , "drugid"]
  cellids <- tSet@sensitivity$info[ , "cellid"]
  cellids <- cellids[grep("///", drugids, invert=TRUE)]
  drugids <- drugids[grep("///", drugids, invert=TRUE)]


  tt <- table(cellids, drugids)
  sensitivity.info[rownames(tt), colnames(tt)] <- tt

  return(sensitivity.info)
}


.summarizePerturbationNumbers <- function(tSet) {

  if (tSet@datasetType != "perturbation" && tSet@datasetType != "both") {
    stop ("Data type must be either perturbation or both")
  }

  ## consider all drugs
  drugn <- rownames(tSet@drug)

  ## consider all cell lines
  celln <- rownames(tSet@cell)

  perturbation.info <- array(0, dim=c(length(celln), length(drugn), length(tSet@molecularProfiles)), dimnames=list(celln, drugn, names((tSet@molecularProfiles))))

  for (i in seq_along(tSet@molecularProfiles)) {
    if (nrow(SummarizedExperiment::colData(tSet@molecularProfiles[[i]])) > 0 &&
        all(
          is.element(c("cellid", "drugid"),
                     colnames(SummarizedExperiment::colData(tSet@molecularProfiles[[i]]))))) {
      tt <- table(SummarizedExperiment::colData(tSet@molecularProfiles[[i]])[ , "cellid"], SummarizedExperiment::colData(tSet@molecularProfiles[[i]])[ , "drugid"])
      perturbation.info[rownames(tt), colnames(tt), names(tSet@molecularProfiles)[i]] <- tt
    }
  }

  return(perturbation.info)
}

#' A function to verify the structure of a ToxicoSet
#'
#' This function checks the structure of a ToxicoSet, ensuring that the
#' correct annotations are in place and all the required slots are filled so
#' that matching of cells and drugs can be properly done across different types
#' of data and with other studies.
#'
#' @examples
#'
#' checkTSetStructure(TGGATESsmall)
#'
#' @param tSet A \code{ToxicoSet} object
#' @param plotDist Should the function also plot the distribution of molecular data?
#' @param result.dir The path to the directory for saving the plots as a string, defaults to `tempdir()`
#' @return Prints out messages whenever describing the errors found in the structure of the pset object passed in.
#' @export
#' @importFrom graphics hist
#' @importFrom grDevices dev.off pdf
checkTSetStructure <-
  function(tSet, plotDist=FALSE, result.dir=".") {
    if(!file.exists(result.dir) & plotDist) { dir.create(result.dir, showWarnings=FALSE, recursive=TRUE) }
    for( i in seq_along(tSet@molecularProfiles)) {
      profile <- tSet@molecularProfiles[[i]]
      nn <- names(tSet@molecularProfiles)[i]
      if((S4Vectors::metadata(profile)$annotation == "rna" | S4Vectors::metadata(profile)$annotation == "rnaseq") & plotDist)
      {
        pdf(file=file.path(result.dir, sprintf("%s.pdf", nn)))
        hist(SummarizedExperiment::assay(profile, 1), breaks = 100)
        dev.off()
      }
      warning(ifelse(nrow(SummarizedExperiment::rowData(profile)) != nrow(SummarizedExperiment::assay(profile, 1)), sprintf("%s: number of features in fData is different from expression slots", nn), sprintf("%s: fData dimension is OK", nn)))
      warning(ifelse(nrow(SummarizedExperiment::colData(profile)) != ncol(SummarizedExperiment::assay(profile, 1)), sprintf("%s: number of cell lines in pData is different from expression slots", nn), sprintf("%s: pData dimension is OK", nn)))
      warning(ifelse("cellid" %in% colnames(SummarizedExperiment::colData(profile)), "", sprintf("%s: cellid does not exist in pData columns", nn)))
      warning(ifelse("batchid" %in% colnames(SummarizedExperiment::colData(profile)), "", sprintf("%s: batchid does not exist in pData columns", nn)))
      if(S4Vectors::metadata(profile)$annotation == "rna" | S4Vectors::metadata(profile)$annotation == "rnaseq")
      {
        warning(ifelse("BEST" %in% colnames(SummarizedExperiment::rowData(profile)), "BEST is OK", sprintf("%s: BEST does not exist in fData columns", nn)))
        warning(ifelse("Symbol" %in% colnames(SummarizedExperiment::rowData(profile)), "Symbol is OK", sprintf("%s: Symbol does not exist in fData columns", nn)))
      }
      if("cellid" %in% colnames(SummarizedExperiment::colData(profile))) {
        if(!all(SummarizedExperiment::colData(profile)[,"cellid"] %in% rownames(tSet@cell))) {
          warning(sprintf("%s: not all the cell lines in this profile are in cell lines slot", nn))
        }
      }else {
        warning(sprintf("%s: cellid does not exist in pData", nn))
      }
    }
    if("tissueid" %in% colnames(tSet@cell)) {
      if("unique.tissueid" %in% colnames(tSet@curation$tissue))
      {
        if(length(intersect(rownames(tSet@curation$tissue), rownames(tSet@cell))) != nrow(tSet@cell)) {
          message("rownames of curation tissue slot should be the same as cell slot (curated cell ids)")
        } else{
          if(length(intersect(tSet@cell$tissueid, tSet@curation$tissue$unique.tissueid)) != length(table(tSet@cell$tissueid))){
            message("tissueid should be the same as unique tissue id from tissue curation slot")
          }
        }
      } else {
        message("unique.tissueid which is curated tissue id across data set should be a column of tissue curation slot")
      }
      if(any(is.na(tSet@cell[,"tissueid"]) | tSet@cell[,"tissueid"] == "", na.rm = TRUE)) {
        message(sprintf("There is no tissue type for this cell line(s): %s", paste(rownames(tSet@cell)[which(is.na(tSet@cell[,"tissueid"]) | tSet@cell[,"tissueid"] == "")], collapse = " ")))
      }
    } else {
      warning("tissueid does not exist in cell slot")
    }

    if("unique.cellid" %in% colnames(tSet@curation$cell)) {
      if(length(intersect(tSet@curation$cell$unique.cellid, rownames(tSet@cell))) != nrow(tSet@cell)) {
        message("rownames of cell slot should be curated cell ids")
      }
    } else {
      message("unique.cellid which is curated cell id across data set should be a column of cell curation slot")
    }

    if(length(intersect(rownames(tSet@curation$cell), rownames(tSet@cell))) != nrow(tSet@cell)) {
      message("rownames of curation cell slot should be the same as cell slot (curated cell ids)")
    }

    if("unique.drugid" %in% colnames(tSet@curation$drug)) {
      if(length(intersect(tSet@curation$drug$unique.drugid, rownames(tSet@drug))) != nrow(tSet@drug)) {
        message("rownames of drug slot should be curated drug ids")
      }
    } else {
      message("unique.drugid which is curated drug id across data set should be a column of drug curation slot")
    }

    if(length(intersect(rownames(tSet@curation$cell), rownames(tSet@cell))) != nrow(tSet@cell)) {
      message("rownames of curation drug slot should be the same as drug slot (curated drug ids)")
    }

    if(!is(tSet@cell, "data.frame")) {
      warning("cell slot class type should be dataframe")
    }
    if(!is(tSet@drug, "data.frame")) {
      warning("drug slot class type should be dataframe")
    }
    if(tSet@datasetType %in% c("sensitivity", "both"))
    {
      if(!is(tSet@sensitivity$info, "data.frame")) {
        warning("sensitivity info slot class type should be dataframe")
      }
      if("cellid" %in% colnames(tSet@sensitivity$info)) {
        if(!all(tSet@sensitivity$info[,"cellid"] %in% rownames(tSet@cell))) {
          warning("not all the cell lines in sensitivity data are in cell slot")
        }
      }else {
        warning("cellid does not exist in sensitivity info")
      }
      if("drugid" %in% colnames(tSet@sensitivity$info)) {
        drug.ids <- unique(tSet@sensitivity$info[,"drugid"])
        drug.ids <- drug.ids[grep("///",drug.ids, invert=TRUE)]
        if(!all(drug.ids %in% rownames(tSet@drug))) {
          message("not all the drugs in sensitivity data are in drug slot")
        }
      }else {
        warning("drugid does not exist in sensitivity info")
      }

      if(any(!is.na(tSet@sensitivity$raw))) {
        if(!all(dimnames(tSet@sensitivity$raw)[[1]] %in% rownames(tSet@sensitivity$info))) {
          warning("For some experiments there is raw sensitivity data but no experimet information in sensitivity info")
        }
      }
      if(!all(rownames(tSet@sensitivity$profiles) %in% rownames(tSet@sensitivity$info))) {
        warning("For some experiments there is sensitivity profiles but no experimet information in sensitivity info")
      }
    }
  }
