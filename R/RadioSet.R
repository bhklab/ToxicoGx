#' A Class to Contain RadioGenomic datasets together with their curations
#'
#' The RadioSet (RSet) class was developed to contain and organise large
#' RadioGenomic datasets, and aid in their metanalysis. It was designed
#' primarily to allow bioinformaticians and biologists to work with data at the
#' level of genes and cell lines, providing a more naturally intuitive
#' interface and simplifying analyses between several datasets. As such, it was
#' designed to be flexible enough to hold datasets of two different natures
#' while providing a common interface. The class can accomidate datasets
#' containing both radiation dose response data, as well as datasets contaning
#' genetic profiles of cell lines pre and post treatement with compounds, known
#' respecitively as sensitivity and perturbation datasets.
#'
#'
#' @slot annotation A \code{list} of annotation data about the RadioSet,
#'    including the \code{$name} and the session information for how the object
#'    was creating, detailing the exact versions of R and all the packages used
#' @slot molecularProfiles A \code{list} containing 4 \code{Biobase::ExpressionSet}
#'   type object for holding data for RNA, DNA, SNP and Copy Number Variation
#'   measurements respectively, with associated \code{fData} and \code{pData}
#'   containing the row and column metadata
#' @slot cell A \code{data.frame} containg the annotations for all the cell
#'   lines profiled in the data set, across all data types
#' @slot radiation A \code{data.frame} containg the annotations for all the radiation treatment types
#'   used in the in the dataset, across all data types
#' @slot sensitivity A \code{list} containing all the data for the sensitivity
#'   experiments, including \code{$info}, a \code{data.frame} containing the
#'   experimental info,\code{$raw} a 3D \code{array} containing raw data,
#'   \code{$profiles}, a \code{data.frame} containing sensitivity profiles
#'   statistics, and \code{$n}, a \code{data.frame} detailing the number of
#'   experiments for each cell-radiation type pair
#' @slot perturbation A \code{list} containting \code{$n}, a \code{data.frame}
#'   summarizing the available perturbation data,
#' @slot curation A \code{list} containing mappings for
#'   \code{cell} and \code{tissue} names used in the data set to universal
#'   identifiers used between different RadioSet objects
#' @slot datasetType A \code{character} string of 'sensitivity',
#'   'perturbation', or both detailing what type of data can be found in the
#'   RadioSet, for proper processing of the data
#' @return An object of the RadioSet class
#' @importClassesFrom CoreGx CoreSet
.RadioSet <- setClass("RadioSet", slots = list(radiation="data.frame"
                                                     # tables="array",
                                                     # table.summary="list",
                                                     # dateCreated="character",
                                                     ),
                                                     contains = "CoreSet")


# The default constructor above does a poor job of explaining the required structure of a RadioSet.
# The constructor function defined below guides the user into providing the required components of the curation and senstivity lists
# and hides the annotation slot which the user does not need to manually fill.
# This also follows the design of the Expression Set class.

#' RadioSet constructor
#'
#' A constructor that simplifies the process of creating RadioSets, as well
#' as creates empty objects for data not provided to the constructor. Only
#' objects returned by this constructor are expected to work with the RadioSet
#' methods. For a much more detailed instruction on creating RadioSets, please
#' see the "CreatingRadioSet" vignette.
#'
#' @examples
#' ## For help creating a RadioSet object, please see the following vignette:
#' browseVignettes("PharmacoGx")
#'
#' @param name A \code{character} string detailing the name of the dataset
#' @param molecularProfiles A \code{list} of ExpressionSet objects containing
#'   molecular profiles
#' @param cell A \code{data.frame} containg the annotations for all the cell
#'   lines profiled in the data set, across all data types
#' @param radiation A \code{data.frame} containg the annotations for all the radiations
#'   profiled in the data set, across all data types
#' @param sensitivityInfo A \code{data.frame} containing the information for the
#'   sensitivity experiments
#' @param sensitivityRaw A 3 Dimensional \code{array} contaning the raw radiation
#'   dose – response data for the sensitivity experiments
#' @param sensitivityProfiles \code{data.frame} containing radiation sensitivity profile
#'   statistics such as IC50 and AUC
#' @param sensitivityN,perturbationN A \code{data.frame} summarizing the
#'   available sensitivity/perturbation data
#' @param curationCell,curationTissue A \code{data.frame} mapping
#'   the names for radiations, cells and tissues used in the data set to universal
#'   identifiers used between different RadioSet objects
#' @param datasetType A \code{character} string of 'sensitivity',
#'   'preturbation', or both detailing what type of data can be found in the
#'   RadioSet, for proper processing of the data
#' @param verify \code{boolean} Should the function verify the RadioSet and
#'   print out any errors it finds after construction?
#' @return An object of class RadioSet
#' @export
#' @import methods
#' @importFrom utils sessionInfo
#' @importFrom stats na.omit
RadioSet <-  function(name,
                          molecularProfiles=list(),
                          cell=data.frame(),
                          radiation=data.frame(),
                          sensitivityInfo=data.frame(),
                          sensitivityRaw=array(dim=c(0,0,0)),
                          sensitivityProfiles=matrix(),
                          sensitivityN=matrix(nrow=0, ncol=0),
                          perturbationN=array(NA, dim=c(0,0,0)),
                          # curationDrug=data.frame(),
                          curationCell = data.frame(),
                          curationTissue = data.frame(),
                          datasetType=c("sensitivity", "perturbation", "both"),
                          verify = TRUE)
{
    datasetType <- match.arg(datasetType)

    annotation <- list()
    annotation$name <- as.character(name)
    annotation$dateCreated <- date()
    annotation$sessionInfo <- sessionInfo()
    annotation$call <- match.call()

    #molecularProfiles <- list("dna"=dna, "rna"=rna, "snp"=snp, "cnv"=cnv)
    for (i in seq_along(molecularProfiles)){
        if (class(molecularProfiles[[i]]) != "ExpressionSet"){
            stop(sprintf("Please provide the %s data as an ExpressionSet", names(molecularProfiles[i])))
        }else{
      Biobase::fData(molecularProfiles[[i]]) <- Biobase::fData(molecularProfiles[[i]])[rownames(Biobase::exprs(molecularProfiles[[i]])), , drop=FALSE]
      Biobase::pData(molecularProfiles[[i]]) <- Biobase::pData(molecularProfiles[[i]])[colnames(Biobase::exprs(molecularProfiles[[i]])), , drop=FALSE]
        }

    }
    #if (class(cell)!="data.frame"){
    #    stop("Please provide the cell line annotations as a data frame.")
    #}
    #if (class(radiation)!="data.frame"){
    #    stop("Please provide the radiation annotations as a data frame.")
    #}

    sensitivity <- list()

    if (!all(rownames(sensitivityInfo) == rownames(sensitivityProfiles) & rownames(sensitivityInfo) == dimnames(sensitivityRaw)[[1]])){
        stop("Please ensure all the row names match between the sensitivity data.")
    }

    sensitivity$info <- as.data.frame(sensitivityInfo, stringsAsFactors = FALSE)
    sensitivity$raw <- sensitivityRaw
    sensitivity$profiles <- as.data.frame(sensitivityProfiles, stringsAsFactors = FALSE)
    sensitivity$n <- sensitivityN

    curation <- list()
    # curation$radiation <- as.data.frame(curationDrug, stringsAsFactors = FALSE)
    curation$cell <- as.data.frame(curationCell, stringsAsFactors = FALSE)
    curation$tissue <- as.data.frame(curationTissue, stringsAsFactors = FALSE)
    ### TODO:: Make sure to fix the curation to check for matching row names to the radiation and cell line matrices!!!!!!


    perturbation <- list()
    perturbation$n <- perturbationN
    if (datasetType == "perturbation" || datasetType == "both") {
        perturbation$info <- "The metadata for the perturbation experiments is available for each molecular type by calling the appropriate info function. \n For example, for RNA transcriptome perturbations, the metadata can be accessed using rnaInfo(rSet)."
    } else {
        perturbation$info <- "Not a perturbation dataset."
    }

    rSet  <- .RadioSet(annotation=annotation, molecularProfiles=molecularProfiles, cell=as.data.frame(cell), radiation=as.data.frame(radiation), datasetType=datasetType, sensitivity=sensitivity, perturbation=perturbation, curation=curation)
    if (verify) { checkRSetStructure(rSet)}
  if(length(sensitivityN) == 0 & datasetType %in% c("sensitivity", "both")) {
    rSet@sensitivity$n <- .summarizeSensitivityNumbers(rSet)
  }
    if(length(perturbationN) == 0  & datasetType %in% c("perturbation", "both")) {
      rSet@perturbation$n <- .summarizePerturbationNumbers(rSet)
    }
  return(rSet)
}

# Initialize global variable
rSet <- NULL

#' cellInfo Generic
#'
#' @param rSet A \code{RadioSet} object
#' @param cSet Parameter name for parent method inherited from CoreGx
#'
#' @return a \code{data.frame} with the cell annotations
#' @importFrom CoreGx cellInfo
#' @describeIn RadioSet Returns the annotations for all the cell lines tested on in the RadioSet
#' @export
setMethod(cellInfo,
          "RadioSet",
          function(cSet=rSet){
   callNextMethod(cSet)
})

#' cellInfo Replace Method
#'
#' @param object A \code{RadioSet} object
#' @param value A replacement value
#'
#' @return Updated \code{RadioSet}
#' @importFrom CoreGx cellInfo<-
#' @describeIn RadioSet Update the cell line annotations
#' @export
setReplaceMethod("cellInfo", signature = signature(object="RadioSet",value="data.frame"), function(object, value){
  if(is.null(rownames(value))){
    stop("Please provide the cell_id as rownames for the cell line annotations")
  }
  object <- callNextMethod(object, value)
  object
})

#' radiationInfo Generic
#'
#' Generic for radiationInfo method
#'
#' @examples

#' radiationInfo(Cleveland_small)
#'
#' @inheritParams cellInfo
#' @param rSet A \code{RadioSet} object
#'
#' @return a \code{data.frame} with the radiation annotations
setGeneric("radiationInfo", function(rSet) standardGeneric("radiationInfo"))
#' @describeIn RadioSet Returns the annotations for all the radiations tested in the RadioSet
#' @export
setMethod(radiationInfo, "RadioSet", function(rSet){
  rSet@radiation
})

#' radiationInfo<- Generic
#'
#' Generic for radiationInfo replace method
#'
#' @examples

#' radiationInfo(Cleveland_small) <- radiationInfo(Cleveland_small)
#'
#' @inheritParams cellInfo
#' @inheritParams cellInfo<-
#'
#' @param object The \code{RadioSet} to replace radiation info in
#' @param value A \code{data.frame} with the new radiation annotations
#'
#' @return Updated \code{RadioSet}
setGeneric("radiationInfo<-", function(object, value) standardGeneric("radiationInfo<-"))
#' @describeIn RadioSet Update the radiation annotations
#' @export
setReplaceMethod("radiationInfo", signature = signature(object="RadioSet",value="data.frame"), function(object, value){
  object@radiation <- value
  object
})

#' phenoInfo Generic
#'
#' Generic for phenoInfo method
#'
#' @examples

#' phenoInfo(Cleveland_small, mDataType="rna")
#'
#' @inheritParams cellInfo
#' @param mDataType A \code{character} with the type of molecular data to return/update
#' @return a \code{data.frame} with the experiment info
# setGeneric("phenoInfo", function(rSet, mDataType) standardGeneric("phenoInfo"))
#' @importFrom CoreGx phenoInfo
#' @describeIn RadioSet Return the experiment info from the given type of molecular data in RadioSet
#' @export
setMethod("phenoInfo",
          signature = c("RadioSet", "character"),
          function(cSet=rSet, mDataType){
  callNextMethod(cSet, mDataType)
})

#' phenoInfo<- Generic
#'
#' Generic for phenoInfo replace method
#'
#' @examples
#'

#' phenoInfo(Cleveland_small, mDataType="rna") <- phenoInfo(Cleveland_small, mDataType="rna")
#'
#' @inheritParams cellInfo<-
#' @inheritParams phenoInfo
#'
#' @return The updated \code{RadioSet}
# setGeneric("phenoInfo<-", function(object, mDataType, value) standardGeneric("phenoInfo<-"))
#' @importMethodsFrom CoreGx phenoInfo<-
#' @describeIn RadioSet Update the the given type of molecular data experiment info in the RadioSet
#' @export
setReplaceMethod("phenoInfo", signature = signature(object="RadioSet", mDataType ="character",value="data.frame"), function(object, mDataType, value){

  object <- callNextMethod(object, mDataType, value)
  object
})

#' molecularProfiles Generic
#'
#' Generic for molecularProfiles method
#'
#' @examples

#' Cleveland_mProf <- molecularProfiles(Cleveland_small, "rna")
#' Cleveland_mProf[1:10,]
#'
#' @inheritParams phenoInfo
#'
#' @return a \code{data.frame} with the experiment info
# setGeneric("molecularProfiles", function(rSet, mDataType) standardGeneric("molecularProfiles"))
#' @describeIn RadioSet Return the given type of molecular data from the RadioSet
#' @importMethodsFrom CoreGx molecularProfiles
#' @export
setMethod(molecularProfiles,
          signature("RadioSet", "character"),
          function(cSet=rSet, mDataType){
   callNextMethod(cSet, mDataType)
})

#' molecularProfiles<- Generic
#'
#' Generic for molecularProfiles replace method
#'
#' @examples

#' molecularProfiles(Cleveland_small, "rna") <- molecularProfiles(Cleveland_small, "rna")
#'
#' @inheritParams phenoInfo<-
#' @return Updated \code{RadioSet}
#' @importMethodsFrom CoreGx molecularProfiles<-
#' @describeIn RadioSet Update the given type of molecular data from the RadioSet
#' @export
setReplaceMethod("molecularProfiles", signature = signature(object="RadioSet", mDataType ="character",value="matrix"), function(object, mDataType, value){

  object <- callNextMethod(object, mDataType, value)
  object
})

#'
#' Generic for featureInfo method
#'
#' @examples

#' featureInfo(Cleveland_small, "rna")[1:10,]
#'
#' @inheritParams phenoInfo
#'
#'
#' @describeIn RadioSet Return the feature info for the given molecular data
#' @importFrom CoreGx featureInfo
#' @export
setMethod("featureInfo",
          signature("RadioSet", "character"),
          function(cSet=rSet, mDataType){
  callNextMethod(cSet, mDataType)
})

#' featureInfo<- Generic
#'
#' Generic for featureInfo replace method
#'
#' @examples

#' featureInfo(Cleveland_small, "rna") <- featureInfo(Cleveland_small, "rna")
#'
#' @inheritParams phenoInfo<-
#'
# @param object The \code{RadioSet} to replace gene annotations in
# @param mDataType The type of molecular data to be updated
# @param value A \code{data.frame} with the new feature annotations
#'
#' @return Updated \code{RadioSet}
#' @describeIn RadioSet Replace the gene info for the molecular data
#' @importMethodsFrom CoreGx featureInfo<-
#' @export
setReplaceMethod("featureInfo", signature = signature(object="RadioSet", mDataType ="character",value="data.frame"), function(object, mDataType, value){

  # if(mDataType %in% names(object@molecularProfiles)){Biobase::fData(object@molecularProfiles[[mDataType]]) <- value}
  object <- callNextMethod(object, mDataType, value)
  object

})

#' sensitivityInfo Generic
#'
#' Generic for sensitivityInfo method
#'
#' @examples

#' sensInf<- sensitivityInfo(Cleveland_small)
#' sensInf[1:10,]
#'
#' @inheritParams cellInfo
#'
#' @return a \code{data.frame} with the experiment info
#' @describeIn RadioSet Return the radiation dose sensitivity experiment info
#' @importMethodsFrom CoreGx sensitivityInfo
#' @export
setMethod(sensitivityInfo,
          "RadioSet",
          function(cSet=rSet){
    callNextMethod(cSet)
})

#' sensitivityInfo<- Generic
#'
#' A generic for the sensitivityInfo replacement method
#'
#'
#' @examples

#' sensitivityInfo(Cleveland_small) <- sensitivityInfo(Cleveland_small)
#'
#' @inheritParams cellInfo<-
# @param object The \code{RadioSet} to update
# @param value A \code{data.frame} with the new sensitivity annotations
#'
#'
#' @return Updated \code{RadioSet}
# setGeneric("sensitivityInfo<-", function(object, value) standardGeneric("sensitivityInfo<-"))
#' @importMethodsFrom CoreGx sensitivityInfo<-
#' @describeIn RadioSet Update the sensitivity experiment info
#' @export
setReplaceMethod("sensitivityInfo", signature = signature(object="RadioSet",value="data.frame"), function(object, value){
    # object@sensitivity$info <- value
    object <- callNextMethod(object, value)
    object
})


#' sensitivityProfiles Generic
#'
#' Generic for sensitivityProfiles method
#'
#' @examples

#' sensitivityProfiles(Cleveland_small)
#'
#' @inheritParams cellInfo
#'
#' @return a \code{data.frame} with the experiment info
# setGeneric("sensitivityProfiles", function(rSet) standardGeneric("sensitivityProfiles"))
#' @describeIn RadioSet Return the phenotypic data for the radiation dose sensitivity
#' @importFrom CoreGx sensitivityProfiles
#' @export
setMethod(sensitivityProfiles,
          "RadioSet",
          function(cSet=rSet){

  callNextMethod(cSet)

})

#' sensitivityProfiles<- Generic
#'
#' A generic for the sensitivityProfiles replacement method
#'
#' @examples

#' sensitivityProfiles(Cleveland_small) <- sensitivityProfiles(Cleveland_small)
#'
#' @inheritParams cellInfo<-
#'
# @param object The \code{RadioSet} to update
# @param value A \code{data.frame} with the new sensitivity profiles. If a matrix object is passed in, converted to data.frame before assignment
#'
#' @return Updated \code{RadioSet}
# setGeneric("sensitivityProfiles<-", function(object, value) standardGeneric("sensitivityProfiles<-"))
#' @importFrom CoreGx sensitivityProfiles<-
#' @describeIn RadioSet Update the phenotypic data for the radiation dose
#'   sensitivity
#' @export
setReplaceMethod("sensitivityProfiles", signature = signature(object="RadioSet",value="data.frame"), function(object, value){

    object <- callNextMethod(object, value)
    object
})
#' @describeIn RadioSet Update the phenotypic data for the radiation dose
#'   sensitivity
#' @export
setReplaceMethod("sensitivityProfiles", signature = signature(object="RadioSet",value="matrix"), function(object, value){

    object <- callNextMethod(object, value)
    object
})

#' sensitivityMeasures Generic
#'
#' A generic for the sensitivityMeasures  method
#'
#' @examples

#' sensitivityMeasures(Cleveland_small)
#'
#'
#' @inheritParams cellInfo
#' @return A \code{character} vector of all the available sensitivity measures
#' @describeIn RadioSet Returns the available sensitivity profile
#'   summaries, for example, whether there are IC50 values available
#' @importFrom CoreGx sensitivityMeasures
#' @export
setMethod(sensitivityMeasures,
          "RadioSet",
          function(cSet=rSet){
  callNextMethod(cSet)

})

#' radiationTypes Generic
#'
#' A generic for the radiationTypes method
#'
#' @examples

#' radType <- radiationTypes(Cleveland_small)
#' radType[1:10]
#'
#' @inheritParams cellInfo
#' @param rSet A \code{RadioSet}
#'
#' @return A vector of the radiation names used in the RadioSet
setGeneric("radiationTypes", function(rSet) standardGeneric("radiationTypes"))
#' @describeIn RadioSet Return the names of the radiations used in the RadioSet
#' @export
setMethod(radiationTypes,
          "RadioSet",
          function(rSet){

  rownames(radiationInfo(rSet))

})

#' radiationTypes<- Generic
#'
#' A generic for the radiationTypes replacement method
#'
#'
#' @examples

#' radiationTypes(Cleveland_small) <- radiationTypes(Cleveland_small)
#'
#' @inheritParams cellInfo<-
#' @param object A \code{RadioSet} object to update
#' @param value A \code{character} vector of the new radiation names
#'
#' @return Updated \code{RadioSet}
setGeneric("radiationTypes<-", function(object, value) standardGeneric("radiationTypes<-"))
#' @describeIn RadioSet Update the radiation names used in the dataset
#' @export
setReplaceMethod("radiationTypes", signature = signature(object="RadioSet",value="character"), function(object, value){

    object <- updateRadId(object, value)
    return(object)
})

#' cellNames Generic
#'
#' A generic for the cellNames method
#'
#' @examples

#' cellNames(Cleveland_small)
#'
#' @inheritParams cellInfo
#'
#' @return A vector of the cell names used in the RadioSet
# setGeneric("cellNames", function(rSet) standardGeneric("cellNames"))
#' @describeIn RadioSet Return the cell names used in the dataset
#' @importFrom CoreGx cellNames
#' @export
setMethod("cellNames",
          "RadioSet",
          function(cSet=rSet){
  callNextMethod(cSet)
})

#' cellNames<- Generic
#'
#' A generic for the cellNames replacement method
#'
#' @examples

#' cellNames(Cleveland_small) <- cellNames(Cleveland_small)
#'
#' @inheritParams cellInfo<-
# @param object The \code{RadioSet} to update
# @param value A \code{character} vector of the new cell names
#'
#'
#' @return Updated \code{RadioSet}
#' @importMethodsFrom CoreGx cellNames<-
#' @describeIn RadioSet Update the cell names used in the dataset
#' @export
setReplaceMethod("cellNames", signature = signature(object="RadioSet",value="character"), function(object, value){
    object <- callNextMethod(object, value)
    return(object)
})

#' fNames Generic
#'
#' A generic for the fNames method
#'
#' @examples

#' fNames(Cleveland_small, "rna")[1:10]
#'
#' @inheritParams cellInfo
#'
#' @return A \code{character} vector of the feature names
#' @describeIn RadioSet Return the feature names used in the dataset
#' @importFrom CoreGx fNames
#' @export
setMethod("fNames",
          signature("RadioSet", "character"),
          function(cSet=rSet, mDataType){
  callNextMethod(cSet, mDataType)
})

###TODO:: Define this method in CoreGx and import it; doesn't work here because no generic defined
# fNames<- Generic
#
# A generic for the feature name replacement method
#
#@examples
# data(Cleveland_small)
# cellNames(Cleveland_small) <- cellNames(Cleveland_small)
#
# @inheritParams phenoInfo<-
# @param value A \code{character} vector of the new feature names
# @return Updated \code{RadioSet}
# @describeIn RadioSet Update the feature names used in the dataset
# @export
#setReplaceMethod("fNames", signature = signature(object="RadioSet",value="character",mDataType="character"), function(object, value){
#      rownames(featureInfo(object, mDataType)) <- value
#})

#' dateCreated Generic
#'
#' A generic for the dateCreated method
#'
#' @examples

#' dateCreated(Cleveland_small)
#'
#' @inheritParams cellInfo
#'
#' @return The date the RadioSet was created
#' @describeIn RadioSet Return the date the RadioSet was created
#' @importFrom CoreGx dateCreated
#' @export
setMethod(dateCreated,
          signature = c("RadioSet"),
          function(cSet=rSet) {
  callNextMethod(cSet)
})


#' rSetName Generic
#'
#' A generic for the rSetName method
#'
#' @examples

#' rSetName <- cSetName
#' rSetName(Cleveland_small)
#'
#' @inheritParams cellInfo
#'
#' @return The name of the RadioSet
#' @describeIn RadioSet Return the name of the RadioSet
#' @importFrom CoreGx cSetName
#' @export
setMethod("cSetName",
          signature = c("RadioSet"),
          function(cSet=rSet){
    callNextMethod(cSet)
})
rSetName <- cSetName
###TODO:: Figure out if there is a better way to rename imported generics/methods

#' pertNumber Generic
#'
#' A generic for the pertNumber method
#'
#' @examples

#' pertNumber(Cleveland_small)
#'
#' @inheritParams cellInfo
#'
#' @return A 3D \code{array} with the number of perturbation experiments per radiation type and cell line, and data type
# setGeneric("pertNumber", function(rSet) standardGeneric("pertNumber"))
#' @describeIn RadioSet Return the summary of available perturbation
#'   experiments
#' @importMethodsFrom CoreGx pertNumber
#' @export
setMethod(pertNumber,
          "RadioSet",
          function(cSet=rSet){
  callNextMethod(cSet)
})


#' sensNumber Generic
#'
#' A generic for the sensNumber method
#'
#' @examples

#' sensNumber(Cleveland_small)
#'
#' @inheritParams cellInfo
#'
#' @return A \code{data.frame} with the number of sensitivity experiments per radiation type and cell line
#' @describeIn RadioSet Return the summary of available sensitivity
#'   experiments
#' @importFrom CoreGx sensNumber
#' @export
setMethod(sensNumber,
          "RadioSet",
          function(cSet=rSet){
    callNextMethod(cSet)
})

#' pertNumber<- Generic
#'
#' A generic for the pertNumber method
#'
#' @examples

#' pertNumber(Cleveland_small) <- pertNumber(Cleveland_small)
#'
#' @inheritParams cellInfo<-
# @param object A \code{RadioSet}
# @param value A new 3D \code{array} with the number of perturbation experiments per radiation type and cell line, and data type
#'
#'
#' @return The updated \code{RadioSet}
# setGeneric("pertNumber<-", function(object, value) standardGeneric("pertNumber<-"))
#' @importMethodsFrom CoreGx pertNumber<-
#' @describeIn RadioSet Update the summary of available perturbation
#'   experiments
#' @export
setReplaceMethod('pertNumber', signature = signature(object="RadioSet",value="array"), function(object, value){

  object <- callNextMethod(object, value)
  object

})

#' sensNumber<- Generic
#'
#' A generic for the sensNumber method
#'
#'
#' @examples

#' sensNumber(Cleveland_small) <- sensNumber(Cleveland_small)
#'
#' @inheritParams cellInfo<-
# @param value A new \code{data.frame} with the number of sensitivity experiments per radiation and cell line
# @param object A \code{RadioSet}
#'
#' @return The updated \code{RadioSet}
# setGeneric("sensNumber<-", function(object, value) standardGeneric("sensNumber<-"))
#' @importMethodsFrom CoreGx sensNumber<-
#' @describeIn RadioSet Update the summary of available sensitivity
#'   experiments
#' @export
setReplaceMethod('sensNumber', signature = signature(object="RadioSet",value="matrix"), function(object, value){

  object <- callNextMethod(object, value)
  object

})

#' Show a RadioSet
#'
#' @param object A \code{RadioSet} object
#'
#' @examples

#' Cleveland_small
#'
#' @return Prints the RadioSet object to the output stream, and returns invisible NULL.
#' @export
setMethod("show", signature=signature(object="RadioSet"),
    function(object) {
        cat("Name: ", rSetName(object), "\n")
        cat("Date Created: ", dateCreated(object), "\n")
    cat("Number of cell lines: ", nrow(cellInfo(object)), "\n")
    cat("Number of radiation types: ", nrow(radiationInfo(object)), "\n")
        if("dna" %in% names(object@molecularProfiles)){cat("DNA: \n");cat("\tDim: ", dim(molecularProfiles(object, mDataType="dna")), "\n")}
      if("rna" %in% names(object@molecularProfiles)){cat("RNA: \n");cat("\tDim: ", dim(molecularProfiles(object, mDataType="rna")), "\n")}
      if("rnaseq" %in% names(object@molecularProfiles)){cat("RNASeq: \n");cat("\tDim: ", dim(molecularProfiles(object, mDataType="rnaseq")), "\n")}
      if("snp" %in% names(object@molecularProfiles)){cat("SNP: \n");cat("\tDim: ", dim(molecularProfiles(object, mDataType="snp")), "\n")}
      if("cnv" %in% names(object@molecularProfiles)){cat("CNV: \n");cat("\tDim: ", dim(molecularProfiles(object, mDataType="cnv")), "\n")}
        cat("Drug pertubation: \n")
        cat("\tPlease look at pertNumber(rSet) to determine number of experiments for each radiation-cell combination.\n")
        cat("Drug sensitivity: \n")
        cat("\tNumber of Experiments: ",nrow(sensitivityInfo(object)),"\n")
        cat("\tPlease look at sensNumber(rSet) to determine number of experiments for each radiation-cell combination.\n")
    })

#' mDataNames
#'
#' Returns the molecular data names for the RadioSet.
#'
#' @examples

#' mDataNames(Cleveland_small)
#'
#' @inheritParams cellInfo
#' @param cSet The parameter
#'
# @param cSet Function parameter name from CoreSet parent class
#' @return Vector of names of the molecular data types
# Imports generic
#' @importFrom CoreGx mDataNames
#' @export
setMethod(
  "mDataNames",
  signature = c("RadioSet"),
  definition = function(cSet=rSet){
    callNextMethod(cSet)
  }
)

#'`[`
#'
#'@param x RSet
#'@param i Cell lines to keep in RSet
#'@param j Drugs to keep in RSet
#'@param ... further arguments
#'@param drop A boolean flag of whether to drop single dimensions or not
#'@return Returns the subsetted RSet
#'@export
setMethod(`[`, "RadioSet", function(x, i, j, ..., drop = FALSE){
  if(is.character(i)&&is.character(j)){
    return(subsetTo(x, cells=i, radiations=j,  molecular.data.cells=i))
  }
  else if(is.numeric(i) && is.numeric(j) && (as.integer(i)==i) && (as.integer(j)==j)){
    return(subsetTo(x, cells=cellNames(x)[i], radiations=radiationTypes(x)[j],  molecular.data.cells=cellNames(x)[i]))
  }
})

#' Get the dimensions of a RadioSet
#'
#' @param x RadioSet
#' @return A named vector with the number of Cells and Drugs in the RadioSet
#' @export
setMethod("dim", signature=signature(x="RadioSet"), function(x){

  return(c(Cells=length(cellNames(x)), Radiation=length(radiationTypes(x))))

})


## FIXED? TODO:: Subset function breaks if it doesnt find cell line in sensitivity info
#' A function to subset a RadioSet to data containing only specified radiations, cells and genes
#'
#' This is the prefered method of subsetting a RadioSet. This function allows
#' abstraction of the data to the level of biologically relevant objects: radiations
#' and cells. The function will automatically go through all of the
#' combined data in the RadioSet and ensure only the requested radiations
#' and cell lines are found in any of the slots. This allows quickly picking out
#' all the experiments for a radiation or cell of interest, as well removes the need
#' to keep track of all the metadata conventions between different datasets.
#'
#' @examples

#' clevelandRadiationTypes  <- radiationTypes(Cleveland_small)
#' clevelandCells <- cellNames(Cleveland_small)
#' RSet <- subsetTo(Cleveland_small,radiationTypes = clevelandRadiationTypes[1],
#'   cells = clevelandCells[1])
#' RSet
#'
#' @param rSet A \code{RadioSet} to be subsetted
#' @param cells A list or vector of cell names as used in the dataset to which
#'   the object will be subsetted. If left blank, then all cells will be left in
#'   the dataset.
#' @param radiationTypes A list or vector of radiation names as used in the dataset to which
#'   the object will be subsetted. If left blank, then all radiationTypes will be left in
#'   the dataset.
#' @param molecular.data.cells A list or vector of cell names to keep in the
#'   molecular data
#' @param keep.controls If the dataset has perturbation type experiments, should
#'   the controls be kept in the dataset? Defaults to true.
#' @param ... Other arguments passed by other function within the package
#' @return A RadioSet with only the selected radiation types and cells
#' @importFrom CoreGx unionList
#' @export
# subsetTo <- function(rSet, cells=NULL, radiationTypes=NULL, exps=NULL, molecular.data.cells=NULL, keep.controls=TRUE) {
subsetTo <- function(rSet, cells=NULL, radiationTypes=NULL, molecular.data.cells=NULL, keep.controls=TRUE, ...) {
  drop=FALSE

  adArgs = list(...)
  if ("exps" %in% names(adArgs)) {
  	exps <- adArgs[["exps"]]
  	if(class(exps)=="data.frame"){
  		exps2 <- exps[[rSetName(rSet)]]
  		names(exps2) <- rownames(exps)
  		exps <- exps2
  	} else{
  		exps <- exps[[rSetName(rSet)]]
  	}
  }else {
    exps <- NULL
  }
  if(!missing(cells)){
    cells <- unique(cells)
  }

  if(!missing(radiationTypes)){
    radiationTypes <- unique(radiationTypes)
  }

  if(!missing(molecular.data.cells)){
    molecular.data.cells <- unique(molecular.data.cells)
  }

    ### TODO:: implement strict subsetting at this level!!!!

    ### the function missing does not work as expected in the context below, because the arguments are passed to the anonymous
    ### function in lapply, so it does not recognize them as missing

  rSet@molecularProfiles <- lapply(rSet@molecularProfiles, function(eset, cells, radiationTypes, molecular.data.cells){

    molecular.data.type <- ifelse(length(grep("rna", Biobase::annotation(eset)) > 0), "rna", Biobase::annotation(eset))
    if (length(grep(molecular.data.type, names(molecular.data.cells))) > 0) {
      cells <- molecular.data.cells[[molecular.data.type]]
    }
      column_indices <- NULL

      if (length(cells)==0 && length(radiationTypes)==0) {
          column_indices <- 0:ncol(eset)
      }
      if(length(cells)==0 && rSet@datasetType=="sensitivity") {
        column_indices <- 0:ncol(eset)
      }

      cell_line_index <- NULL
      if(length(cells)!=0) {
        if (!all(cells %in% cellNames(rSet))) {
              stop("Some of the cell names passed to function did not match to names in the RadioSet. Please ensure you are using cell names as returned by the cellNames function")
        }
          cell_line_index <- which(Biobase::pData(eset)[["cellid"]] %in% cells)
        # if (length(na.omit(cell_line_index))==0){
    #       stop("No cell lines matched")
    #     }
      }
      radiationTypes_index <- NULL
      if(rSet@datasetType=="perturbation" || rSet@datasetType=="both"){
        if(length(radiationTypes) != 0) {
            if (!all(radiationTypes %in% radiationTypes(rSet))){
                  stop("Some of the radiation names passed to function did not match to names in the RadioSet. Please ensure you are using radiation names as returned by the radiationTypes function")
            }
          radiationTypes_index <- which(Biobase::pData(eset)[["radiation.type"]] %in% radiationTypes)
          # if (length(radiationTypes_index)==0){
    #         stop("No radiationTypes matched")
    #       }
          if(keep.controls) {
            control_indices <- which(Biobase::pData(eset)[["xptype"]]=="control")
            radiationTypes_index <- c(radiationTypes_index, control_indices)
          }
        }
      }

      if(length(radiationTypes_index) != 0 && length(cell_line_index) != 0) {
        if(length(intersect(radiationTypes_index, cell_line_index)) == 0) {
          stop("This Drug - Cell Line combination was not tested together.")
        }
        column_indices <- intersect(radiationTypes_index, cell_line_index)
      } else {
        if(length(radiationTypes_index) !=0) {
        column_indices <- radiationTypes_index
      }
        if(length(cell_line_index) !=0) {
        column_indices <- cell_line_index
      }
      }

      row_indices <- 0:nrow(Biobase::exprs(eset))

      eset <- eset[row_indices,column_indices]
      return(eset)

  }, cells=cells, radiationTypes=radiationTypes, molecular.data.cells=molecular.data.cells)

  if ((rSet@datasetType == "sensitivity" | rSet@datasetType == "both") & length(exps) != 0) {
      rSet@sensitivity$info <- rSet@sensitivity$info[exps, , drop=drop]
      rownames(rSet@sensitivity$info) <- names(exps)
      if(length(rSet@sensitivity$raw) > 0) {
        rSet@sensitivity$raw <- rSet@sensitivity$raw[exps, , , drop=drop]
        dimnames(rSet@sensitivity$raw)[[1]] <- names(exps)
      }
      rSet@sensitivity$profiles <- rSet@sensitivity$profiles[exps, , drop=drop]
      rownames(rSet@sensitivity$profiles) <- names(exps)

      rSet@sensitivity$n <- .summarizeSensitivityNumbers(rSet)
  }
  else if ((rSet@datasetType == "sensitivity" | rSet@datasetType == "both") & (length(radiationTypes) != 0 | length(cells) != 0)) {

        radiationTypes_index <- which (sensitivityInfo(rSet)[, "radiation.type"] %in% radiationTypes)
        cell_line_index <- which (sensitivityInfo(rSet)[,"cellid"] %in% cells)
        if (length(radiationTypes_index) !=0 & length(cell_line_index) !=0 ) {
          if (length(intersect(radiationTypes_index, cell_line_index)) == 0) {
            stop("This Drug - Cell Line combination was not tested together.")
          }
          row_indices <- intersect(radiationTypes_index, cell_line_index)
        } else {
          if(length(radiationTypes_index)!=0 & length(cells)==0) {
                row_indices <- radiationTypes_index
          } else {
              if(length(cell_line_index)!=0 & length(radiationTypes)==0){
                  row_indices <- cell_line_index
              } else {
              row_indices <- vector()
              }
          }
       }
        rSet@sensitivity[names(rSet@sensitivity)[names(rSet@sensitivity)!="n"]] <- lapply(rSet@sensitivity[names(rSet@sensitivity)[names(rSet@sensitivity)!="n"]], function(x,i, drop){
            #browser()
          if (length(dim(x))==2){
            return(x[i,,drop=drop])
          }
          if (length(dim(x))==3){
            return(x[i,,,drop=drop])
          }
          }, i=row_indices, drop=drop)
  }

	if (length(radiationTypes)==0) {
		if(rSet@datasetType == "sensitivity" | rSet@datasetType == "both"){
			radiationTypes <- unique(sensitivityInfo(rSet)[["radiation.type"]])
		}
		if(rSet@datasetType == "perturbation" | rSet@datasetType == "both"){
			radiationTypes <- union(radiationTypes, na.omit(unionList(lapply(rSet@molecularProfiles, function(eSet){unique(Biobase::pData(eSet)[["radiation.type"]])}))))
		}
	}
	if (length(cells)==0) {
		cells <- union(cells, na.omit(unionList(lapply(rSet@molecularProfiles, function(eSet){unique(Biobase::pData(eSet)[["cellid"]])}))))
        if (rSet@datasetType =="sensitivity" | rSet@datasetType == "both"){
            cells <- union(cells, sensitivityInfo(rSet)[["cellid"]])
        }
	}
	radiationInfo(rSet) <- radiationInfo(rSet)[radiationTypes , , drop=drop]
	cellInfo(rSet) <- cellInfo(rSet)[cells , , drop=drop]
	rSet@curation$radiation <- rSet@curation$radiation[radiationTypes , , drop=drop]
	rSet@curation$cell <- rSet@curation$cell[cells , , drop=drop]
	rSet@curation$tissue <- rSet@curation$tissue[cells , , drop=drop]
	if (rSet@datasetType == "sensitivity" | rSet@datasetType == "both"  & length(exps) == 0) {
	  rSet@sensitivity$n <- rSet@sensitivity$n[cells, radiationTypes , drop=drop]
	}
	if (rSet@datasetType == "perturbation" | rSet@datasetType == "both") {
	    rSet@perturbation$n <- rSet@perturbation$n[cells,radiationTypes, , drop=drop]
    }
      return(rSet)
}
### TODO:: Add updating of sensitivity Number tables
updateCellId <- function(rSet, new.ids = vector("character")){

  if (length(new.ids)!=nrow(cellInfo(rSet))){
    stop("Wrong number of cell identifiers")
  }

  if(rSet@datasetType=="sensitivity"|rSet@datasetType=="both"){
    myx <- match(sensitivityInfo(rSet)[,"cellid"],rownames(cellInfo(rSet)))
    sensitivityInfo(rSet)[,"cellid"] <- new.ids[myx]

  }


  rSet@molecularProfiles <- lapply(rSet@molecularProfiles, function(eset){

      myx <- match(Biobase::pData(eset)[["cellid"]],rownames(cellInfo(rSet)))
      Biobase::pData(eset)[["cellid"]]  <- new.ids[myx]
      return(eset)
        })





  if(any(duplicated(new.ids))){
    warning("Duplicated ids passed to updateCellId. Merging old ids into the same identifier")

    if(ncol(sensNumber(rSet))>0){
      sensMatch <- match(rownames(sensNumber(rSet)), rownames(cellInfo(rSet)))
    }
    if(dim(pertNumber(rSet))[[2]]>0){
      pertMatch <- match(dimnames(pertNumber(rSet))[[1]], rownames(cellInfo(rSet)))
    }
    curMatch <- match(rownames(rSet@curation$cell),rownames(cellInfo(rSet)))

    duplId <- unique(new.ids[duplicated(new.ids)])
    for(id in duplId){

      if (ncol(sensNumber(rSet))>0){
        myx <- which(new.ids[sensMatch] == id)
        sensNumber(rSet)[myx[1],] <- apply(sensNumber(rSet)[myx,], 2, sum)
        sensNumber(rSet) <- sensNumber(rSet)[-myx[-1],]
        # sensMatch <- sensMatch[-myx[-1]]
      }
      if (dim(pertNumber(rSet))[[1]]>0){
        myx <- which(new.ids[pertMatch] == id)
        pertNumber(rSet)[myx[1],,] <- apply(pertNumber(rSet)[myx,,], c(1,3), sum)
        pertNumber(rSet) <- pertNumber(rSet)[-myx[-1],,]
        # pertMatch <- pertMatch[-myx[-1]]
      }

      myx <- which(new.ids[curMatch] == id)
      rSet@curation$cell[myx[1],] <- apply(rSet@curation$cell[myx,], 2, paste, collapse="///")
      rSet@curation$cell <- rSet@curation$cell[-myx[-1],]
      rSet@curation$tissue[myx[1],] <- apply(rSet@curation$tissue[myx,], 2, paste, collapse="///")
      rSet@curation$tissue <- rSet@curation$tissue[-myx[-1],]
      # curMatch <- curMatch[-myx[-1]]

      myx <- which(new.ids == id)
      cellInfo(rSet)[myx[1],] <- apply(cellInfo(rSet)[myx,], 2, paste, collapse="///")
      cellInfo(rSet) <- cellInfo(rSet)[-myx[-1],]
      new.ids <- new.ids[-myx[-1]]
      if(ncol(sensNumber(rSet))>0){
        sensMatch <- match(rownames(sensNumber(rSet)), rownames(cellInfo(rSet)))
      }
      if(dim(pertNumber(rSet))[[1]]>0){
        pertMatch <- match(dimnames(pertNumber(rSet))[[1]], rownames(cellInfo(rSet)))
      }
      curMatch <- match(rownames(rSet@curation$cell),rownames(cellInfo(rSet)))
    }
  } else {
    if (dim(pertNumber(rSet))[[1]]>0){
      pertMatch <- match(dimnames(pertNumber(rSet))[[1]], rownames(cellInfo(rSet)))
    }
    if (ncol(sensNumber(rSet))>0){
      sensMatch <- match(rownames(sensNumber(rSet)), rownames(cellInfo(rSet)))
    }
    curMatch <- match(rownames(rSet@curation$cell),rownames(cellInfo(rSet)))
  }

  if (dim(pertNumber(rSet))[[1]]>0){
    dimnames(pertNumber(rSet))[[1]] <- new.ids[pertMatch]
  }
  if (ncol(sensNumber(rSet))>0){
    rownames(sensNumber(rSet)) <- new.ids[sensMatch]
  }
  rownames(rSet@curation$cell) <- new.ids[curMatch]
  rownames(rSet@curation$tissue) <- new.ids[curMatch]
  rownames(cellInfo(rSet)) <- new.ids





  # myx <- match(rownames(rSet@curation$cell),rownames(cellInfo(rSet)))
  # rownames(rSet@curation$cell) <- new.ids[myx]
  # rownames(rSet@curation$tissue) <- new.ids[myx]
  # if (dim(pertNumber(rSet))[[1]]>0){
  #   myx <- match(dimnames(pertNumber(rSet))[[1]], rownames(cellInfo(rSet)))
  #   dimnames(pertNumber(rSet))[[1]] <- new.ids[myx]
  # }
  # if (nrow(sensNumber(rSet))>0){
  #   myx <- match(rownames(sensNumber(rSet)), rownames(cellInfo(rSet)))
  #   rownames(sensNumber(rSet)) <- new.ids[myx]
  # }
  # rownames(cellInfo(rSet)) <- new.ids
  return(rSet)

}

# updateFeatureNames <- function(rSet, new.ids = vector("character")){
#
#   if (length(new.ids)!=nrow(cellInfo(rSet))){
#     stop("Wrong number of cell identifiers")
#   }
#
#   if(rSet@datasetType=="sensitivity"|rSet@datasetType=="both"){
#     myx <- match(sensitivityInfo(rSet)[,"cellid"],rownames(cellInfo(rSet)))
#     sensitivityInfo(rSet)[,"cellid"] <- new.ids[myx]
#
#   }
#
#   rSet@molecularProfiles <- lapply(rSet@molecularProfiles, function(eset){
#
#     myx <- match(pData(eset)[["cellid"]],rownames(cellInfo(rSet)))
#     pData(eset)[["cellid"]]  <- new.ids[myx]
#     return(eset)
#       })
#   myx <- match(rownames(rSet@curation$cell),rownames(cellInfo(rSet)))
#   rownames(rSet@curation$cell) <- new.ids[myx]
#   rownames(rSet@curation$tissue) <- new.ids[myx]
#   if (dim(pertNumber(rSet))[[1]]>0){
#     myx <- match(dimnames(pertNumber(rSet))[[1]], rownames(cellInfo(rSet)))
#     dimnames(pertNumber(rSet))[[1]] <- new.ids[myx]
#   }
#   if (nrow(sensNumber(rSet))>0){
#     myx <- match(rownames(sensNumber(rSet)), rownames(cellInfo(rSet)))
#     rownames(sensNumber(rSet)) <- new.ids[myx]
#   }
#   rownames(cellInfo(rSet)) <- new.ids
#   return(rSet)
#
# }

### TODO:: Add updating of sensitivity Number tables
updateRadId <- function(rSet, new.ids = vector("character")){

  if (length(new.ids)!=nrow(radiationInfo(rSet))){
     stop("Wrong number of radiation identifiers")
  }

   if(rSet@datasetType=="sensitivity"|rSet@datasetType=="both"){
     myx <- match(sensitivityInfo(rSet)[,"radiation.type"],rownames(radiationInfo(rSet)))
     sensitivityInfo(rSet)[,"radiation.type"] <- new.ids[myx]

   }
   if(rSet@datasetType=="perturbation"|rSet@datasetType=="both"){
     rSet@molecularProfiles <- lapply(rSet@molecularProfiles, function(eset){

       myx <- match(Biobase::pData(eset)[["radiation.type"]],rownames(radiationInfo(rSet)))
       Biobase::pData(eset)[["radiation.type"]]  <- new.ids[myx]
       return(eset)
     })
   }


   if(any(duplicated(new.ids))){
     warning("Duplicated ids passed to updateDrugId. Merging old ids into the same identifier")

     if(ncol(sensNumber(rSet))>0){
       sensMatch <- match(colnames(sensNumber(rSet)), rownames(radiationInfo(rSet)))
     }
     if(dim(pertNumber(rSet))[[2]]>0){
       pertMatch <- match(dimnames(pertNumber(rSet))[[2]], rownames(radiationInfo(rSet)))
     }
     curMatch <- match(rownames(rSet@curation$radiation),rownames(radiationInfo(rSet)))

     duplId <- unique(new.ids[duplicated(new.ids)])
     for(id in duplId){

       if (ncol(sensNumber(rSet))>0){
         myx <- which(new.ids[sensMatch] == id)
         sensNumber(rSet)[,myx[1]] <- apply(sensNumber(rSet)[,myx], 1, sum)
         sensNumber(rSet) <- sensNumber(rSet)[,-myx[-1]]
         sensMatch <- sensMatch[-myx[-1]]
       }
       if (dim(pertNumber(rSet))[[2]]>0){
         myx <- which(new.ids[pertMatch] == id)
         pertNumber(rSet)[,myx[1],] <- apply(pertNumber(rSet)[,myx,], c(1,3), sum)
         pertNumber(rSet) <- pertNumber(rSet)[,-myx[-1],]
         pertMatch <- pertMatch[-myx[-1]]
       }

       myx <- which(new.ids[curMatch] == id)
       rSet@curation$radiation[myx[1],] <- apply(rSet@curation$radiation[myx,], 2, paste, collapse="///")
       rSet@curation$radiation <- rSet@curation$radiation[-myx[-1],]
       curMatch <- curMatch[-myx[-1]]

       myx <- which(new.ids == id)
       radiationInfo(rSet)[myx[1],] <- apply(radiationInfo(rSet)[myx,], 2, paste, collapse="///")
       radiationInfo(rSet) <- radiationInfo(rSet)[-myx[-1],]
       new.ids <- new.ids[-myx[-1]]
       if(ncol(sensNumber(rSet))>0){
         sensMatch <- match(colnames(sensNumber(rSet)), rownames(radiationInfo(rSet)))
       }
       if(dim(pertNumber(rSet))[[2]]>0){
         pertMatch <- match(dimnames(pertNumber(rSet))[[2]], rownames(radiationInfo(rSet)))
       }
       curMatch <- match(rownames(rSet@curation$radiation),rownames(radiationInfo(rSet)))
     }
   } else {
     if (dim(pertNumber(rSet))[[2]]>0){
       pertMatch <- match(dimnames(pertNumber(rSet))[[2]], rownames(radiationInfo(rSet)))
     }
     if (ncol(sensNumber(rSet))>0){
       sensMatch <- match(colnames(sensNumber(rSet)), rownames(radiationInfo(rSet)))
     }
     curMatch <- match(rownames(rSet@curation$radiation),rownames(radiationInfo(rSet)))
   }

   if (dim(pertNumber(rSet))[[2]]>0){
     dimnames(pertNumber(rSet))[[2]] <- new.ids[pertMatch]
   }
   if (ncol(sensNumber(rSet))>0){
     colnames(sensNumber(rSet)) <- new.ids[sensMatch]
   }
   #rownames(rSet@curation$radiation) <- new.ids[curMatch]
   rownames(radiationInfo(rSet)) <- new.ids


   return(rSet)
 }

.summarizeSensitivityNumbers <- function(rSet) {

  if (rSet@datasetType != "sensitivity" && rSet@datasetType != "both") {
    stop ("Data type must be either sensitivity or both")
  }

  ## unique radiation identifiers
  # radiationn <- sort(unique(rSet@sensitivity$info[ , "radiation.type"]))

  ## consider all radiations
  radiationn <- rownames(rSet@radiation)

  ## unique radiation identifiers
  # celln <- sort(unique(rSet@sensitivity$info[ , "cellid"]))

  ## consider all cell lines
  celln <- rownames(rSet@cell)

  sensitivity.info <- matrix(0, nrow=length(celln), ncol=length(radiationn), dimnames=list(celln, radiationn))
  radiation.types <- rSet@sensitivity$info[ , "radiation.type"]
  cellids <- rSet@sensitivity$info[ , "cellid"]
  cellids <- cellids[grep("///", radiation.types, invert=TRUE)]
  radiation.types <- radiation.types[grep("///", radiation.types, invert=TRUE)]


  tt <- table(cellids, radiation.types)
  sensitivity.info[rownames(tt), colnames(tt)] <- tt

    return(sensitivity.info)
}

.summarizeMolecularNumbers <- function(rSet) {

  ## consider all molecular types
  mDT <- mDataNames(rSet)

  ## consider all cell lines
  celln <- rownames(rSet@cell)

  molecular.info <- matrix(0, nrow=length(celln), ncol=length(mDT), dimnames=list(celln, mDT))

  for(mDataType in mDT) {
    tt <- table(phenoInfo(rSet, mDataType)$cellid)
    molecular.info[names(tt), mDataType] <- tt

  }
  return(molecular.info)
}


.summarizePerturbationNumbers <- function(rSet) {

  if (rSet@datasetType != "perturbation" && rSet@datasetType != "both") {
    stop ("Data type must be either perturbation or both")
  }

  ## unique radiation identifiers
  # radiationn <- sort(unique(unlist(lapply(rSet@molecularProfiles, function (x) {
  #   res <- NULL
  #   if (nrow(pData(x)) > 0 & "radiation.type" %in% colnames(pData(x))) {
  #     res <- pData(x)[ , "radiation.type"]
  #   }
  #   return (res)
  # }))))

  ## consider all radiations
  radiationn <- rownames(rSet@radiation)

  ## unique cell line identifiers
  # celln <- sort(unique(unlist(lapply(rSet@molecularProfiles, function (x) {
  #   res <- NULL
  #   if (nrow(pData(x)) > 0 & "cellid" %in% colnames(pData(x))) {
  #     res <- pData(x)[ , "cellid"]
  #   }
  #   return (res)
  # }))))

  ## consider all cell lines
  celln <- rownames(rSet@cell)

  perturbation.info <- array(0, dim=c(length(celln), length(radiationn), length(rSet@molecularProfiles)), dimnames=list(celln, radiationn, names((rSet@molecularProfiles))))

    for (i in 1:length(rSet@molecularProfiles)) {
      if (nrow(Biobase::pData(rSet@molecularProfiles[[i]])) > 0 && all(is.element(c("cellid", "radiation.type"), colnames(Biobase::pData(rSet@molecularProfiles[[i]]))))) {
      tt <- table(Biobase::pData(rSet@molecularProfiles[[i]])[ , "cellid"], Biobase::pData(rSet@molecularProfiles[[i]])[ , "radiation.type"])
        perturbation.info[rownames(tt), colnames(tt), names(rSet@molecularProfiles)[i]] <- tt
      }
    }

    return(perturbation.info)
}

#' A function to verify the structure of a RadioSet
#'
#' This function checks the structure of a PharamcoSet, ensuring that the
#' correct annotations are in place and all the required slots are filled so
#' that matching of cells and radiations can be properly done across different types
#' of data and with other studies.
#'
#' @examples

#'
#' checkRSetStructure(Cleveland_small)
#'
#' @param rSet A \code{RadiOSet} object
#' @param plotDist Should the function also plot the distribution of molecular data?
#' @param result.dir The path to the directory for saving the plots as a string, defaults to `tempdir()``
#' @return Prints out messages whenever describing the errors found in the structure of the pset object passed in.
#' @export
#' @importFrom graphics hist
#' @importFrom grDevices dev.off pdf

checkRSetStructure <-
  function(rSet, plotDist=FALSE, result.dir=tempdir()) {
    if(!file.exists(result.dir) & plotDist) { dir.create(result.dir, showWarnings=FALSE, recursive=TRUE) }
    for( i in seq_along(rSet@molecularProfiles)) {
      profile <- rSet@molecularProfiles[[i]]
      nn <- names(rSet@molecularProfiles)[i]
      if((Biobase::annotation(profile) == "rna" | Biobase::annotation(profile) == "rnaseq") & plotDist)
      {
        pdf(file=file.path(result.dir, sprintf("%s.pdf", nn)))
        hist(Biobase::exprs(profile), breaks = 100)
        dev.off()
      }
      warning(ifelse(nrow(Biobase::fData(profile)) != nrow(Biobase::exprs(profile)), sprintf("%s: number of features in fData is different from expression slots", nn), sprintf("%s: fData dimension is OK", nn)))
      warning(ifelse(nrow(Biobase::pData(profile)) != ncol(Biobase::exprs(profile)), sprintf("%s: number of cell lines in pData is different from expression slots", nn), sprintf("%s: pData dimension is OK", nn)))
      warning(ifelse("cellid" %in% colnames(Biobase::pData(profile)), "", sprintf("%s: cellid does not exist in pData columns", nn)))
      warning(ifelse("batchid" %in% colnames(Biobase::pData(profile)), "", sprintf("%s: batchid does not exist in pData columns", nn)))
      if(Biobase::annotation(profile) == "rna" | Biobase::annotation(profile) == "rnaseq")
      {
        warning(ifelse("BEST" %in% colnames(Biobase::fData(profile)), "BEST is OK", sprintf("%s: BEST does not exist in fData columns", nn)))
        warning(ifelse("Symbol" %in% colnames(Biobase::fData(profile)), "Symbol is OK", sprintf("%s: Symbol does not exist in fData columns", nn)))
      }
      if("cellid" %in% colnames(Biobase::pData(profile))) {
        if(!all(Biobase::pData(profile)[,"cellid"] %in% rownames(rSet@cell))) {
          warning(sprintf("%s: not all the cell lines in this profile are in cell lines slot", nn))
        }
      }else {
        warning(sprintf("%s: cellid does not exist in pData", nn))
      }
    }
    if("tissueid" %in% colnames(rSet@cell)) {
      if("unique.tissueid" %in% colnames(rSet@curation$tissue))
      {
        if(length(intersect(rownames(rSet@curation$tissue), rownames(rSet@cell))) != nrow(rSet@cell)) {
          message("rownames of curation tissue slot should be the same as cell slot (curated cell ids)")
        } else{
          if(length(intersect(rSet@cell$tissueid, rSet@curation$tissue$unique.tissueid)) != length(table(rSet@cell$tissueid))){
            message("tissueid should be the same as unique tissue id from tissue curation slot")
          }
        }
      } else {
        message("unique.tissueid which is curated tissue id across data set should be a column of tissue curation slot")
      }
      if(any(is.na(rSet@cell[,"tissueid"]) | rSet@cell[,"tissueid"]=="", na.rm=TRUE)){
        message(sprintf("There is no tissue type for this cell line(s): %s", paste(rownames(rSet@cell)[which(is.na(rSet@cell[,"tissueid"]) | rSet@cell[,"tissueid"]=="")], collapse=" ")))
      }
    } else {
      warning("tissueid does not exist in cell slot")
    }

    if("unique.cellid" %in% colnames(rSet@curation$cell)) {
      if(length(intersect(rSet@curation$cell$unique.cellid, rownames(rSet@cell))) != nrow(rSet@cell)) {
        print("rownames of cell slot should be curated cell ids")
      }
    } else {
      print("unique.cellid which is curated cell id across data set should be a column of cell curation slot")
    }
#     if("cellid" %in% colnames(rSet@cell)) {
#       if(length(intersect(rSet@curation$cell$cellid, rownames(rSet@cell))) != nrow(rSet@cell)) {
#         print("values of cellid column should be curated cell line ids")
#       }
#     } else {
#       print("cellid which is curated cell id across data set should be a column of cell slot")
#     }

    if(length(intersect(rownames(rSet@curation$cell), rownames(rSet@cell))) != nrow(rSet@cell)) {
      print("rownames of curation cell slot should be the same as cell slot (curated cell ids)")
    }

    ##TODO:: Determine if rSet@curation$radiation is intended to be defined in rSet objects
    ## It is not currently defined in rSet class defintion of this package
    #if("unique.radiation.type" %in% colnames(rSet@curation$radiation)) {
    #  if(length(intersect(rSet@curation$radiation$unique.radiation.type, rownames(rSet@radiation))) != nrow(rSet@radiation)) {
    #    print("rownames of radiation slot should be curated radiation ids")
    #  }
    #} else {
    #  print("unique.radiation.type which is curated radiation id across data set should be a column of radiation curation slot")
    #}

#     if("radiation.type" %in% colnames(rSet@radiation)) {
#       if(length(intersect(rSet@curation$radiation$radiation.type, rownames(rSet@radiation))) != nrow(rSet@radiation)) {
#         print("values of radiation.type column should be curated radiation ids")
#       }
#     } else {
#       print("radiation.type which is curated radiation id across data set should be a column of radiation slot")
#     }

    if(length(intersect(rownames(rSet@curation$cell), rownames(rSet@cell))) != nrow(rSet@cell)) {
      print("rownames of curation radiation slot should be the same as radiation slot (curated radiation ids)")
    }

    if(class(rSet@cell) != "data.frame") {
      warning("cell slot class type should be dataframe")
    }
    if(class(rSet@radiation) != "data.frame") {
      warning("radiation slot class type should be dataframe")
    }
    if(rSet@datasetType %in% c("sensitivity", "both"))
    {
      if(class(rSet@sensitivity$info) != "data.frame") {
        warning("sensitivity info slot class type should be dataframe")
      }
      if("cellid" %in% colnames(rSet@sensitivity$info)) {
        if(!all(rSet@sensitivity$info[,"cellid"] %in% rownames(rSet@cell))) {
          warning("not all the cell lines in sensitivity data are in cell slot")
        }
      }else {
        warning("cellid does not exist in sensitivity info")
      }
      if("radiation.type" %in% colnames(rSet@sensitivity$info)) {
        radiation.ids <- unique(rSet@sensitivity$info[,"radiation.type"])
        radiation.ids <- radiation.ids[grep("///",radiation.ids, invert=TRUE)]
        if(!all(radiation.ids %in% rownames(rSet@radiation))) {
          print("not all the radiations in sensitivity data are in radiation slot")
        }
      }else {
        warning("radiation.type does not exist in sensitivity info")
      }

      if(any(!is.na(rSet@sensitivity$raw))) {
        if(!all(dimnames(rSet@sensitivity$raw)[[1]] %in% rownames(rSet@sensitivity$info))) {
          warning("For some experiments there is raw sensitivity data but no experimet information in sensitivity info")
        }
      }
      if(!all(rownames(rSet@sensitivity$profiles) %in% rownames(rSet@sensitivity$info))) {
        warning("For some experiments there is sensitivity profiles but no experimet information in sensitivity info")
      }
    }
  }

