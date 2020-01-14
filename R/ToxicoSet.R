#' Class to contain Toxico-genomic Data
#'
#' A description which has yet to be added to this class. This is just a place
#' holder.
#'
#' @inheritParams CoreSet
#' @slot drug A \code{data.frame} containg the annotations for all the drugs
#'   profiled in the in the dataset, across all data types
#'
#' @return An object of the ToxicoSet class
#'
#' @importClassesFrom CoreGx CoreSet
#' @export
.ToxicoSet <- setClass("ToxicoSet", slots = list(drug = "data.frame"),
                       contains = "CoreSet")


# The default constructor above does a poor job of explaining the required structure of a ToxicoSet.
# The constructor function defined below guides the user into providing the required components of
# the curation and senstivity lists and hides the annotation slot which the user does not need to
# manually fill. This also follows the design of the Expression Set class.


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
#' @return An object of class ToxicoSet
#'
#' @import methods
#' @importFrom utils sessionInfo
#' @importFrom stats na.omit
#'
#' @export
#'
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
  datasetType <- match.arg(datasetType)

  annotation <- list()
  annotation$name <- as.character(name)
  annotation$dateCreated <- date()
  annotation$sessionInfo <- sessionInfo()
  annotation$call <- match.call()

  #molecularProfiles <- list("dna"=dna, "rna"=rna, "snp"=snp, "cnv"=cnv)
  for (i in seq_along(molecularProfiles)) {
    if (!is(molecularProfiles[[i]], "ExpressionSet")) {
      stop(sprintf("Please provide the %s data as an ExpressionSet", names(molecularProfiles[i])))
    }else{
      Biobase::fData(molecularProfiles[[i]]) <- Biobase::fData(molecularProfiles[[i]])[rownames(Biobase::exprs(molecularProfiles[[i]])), , drop = FALSE]
      Biobase::pData(molecularProfiles[[i]]) <- Biobase::pData(molecularProfiles[[i]])[colnames(Biobase::exprs(molecularProfiles[[i]])), , drop = FALSE]
    }
  }

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
  curation$drug <- as.data.frame(curationDrug, stringsAsFactors = FALSE)
  curation$tissue <- as.data.frame(curationTissue, stringsAsFactors = FALSE)
  ### TODO:: Make sure to fix the curation to check for matching row names to the radiation and cell line matrices!!!!!!


  perturbation <- list()
  perturbation$n <- perturbationN
  if (datasetType == "perturbation" || datasetType == "both") {
    perturbation$info <- "The metadata for the perturbation experiments is available for each molecular type by calling the appropriate info function. \n For example, for RNA transcriptome perturbations, the metadata can be accessed using rnaInfo(tSet)."
  } else {
    perturbation$info <- "Not a perturbation dataset."
  }

  tSet  <- .ToxicoSet(annotation=annotation, molecularProfiles=molecularProfiles, cell=as.data.frame(cell), drug=as.data.frame(drug), datasetType=datasetType, sensitivity=sensitivity, perturbation=perturbation, curation=curation)
  if (verify) { checkTSetStructure(tSet)}
  if(length(sensitivityN) == 0 & datasetType %in% c("sensitivity", "both")) {
    tSet@sensitivity$n <- .summarizeSensitivityNumbers(tSet)
  }
  if(length(perturbationN) == 0  & datasetType %in% c("perturbation", "both")) {
    tSet@perturbation$n <- .summarizePerturbationNumbers(tSet)
  }
  return(tSet)
}

## TODO:: Implement this in CoreGx
#' tSet Name
#'
#' Retrieves the name of a tSet
#'
#' @examples
#' names(TGGATESsmall)
#'
#' @param tSet [ToxicoSet] A ToxcioSet object
#' @param x [param] The named parameter from the base R names function. For
#'   internal use only.
#'
#' @return [character] A string of the tSet's name
#'
#' @export
setMethod(names,
          "ToxicoSet",
          function(x=tSet) {
            x@annotation$name
          })


## TODO:: Find more elegant solution to 'No visibile global function definition'
tSet <- NULL

#' cellInfo Generic
#'
#' A generic for cellInfo method
#'
#' @examples
#' data(TGGATESsmall)
#' cellInfo <- cellInfo(TGGATESsmall)
#'
#' @param tSet A \code{ToxicoSet} object
#' @param cSet Parameter name for parent method inherited from CoreGx
#'
#' @return a \code{data.frame} with the cell annotations
#'
#' @importFrom CoreGx cellInfo
#'
#' @describeIn ToxicoSet Returns the annotations for all the cell lines tested on in the ToxicoSet
#'
#' @export
setMethod(cellInfo,
          "ToxicoSet",
          function(cSet=tSet){
            callNextMethod(cSet)
          })

#' cellInfo Replace Method
#'
#' Generic for cellInfo replace method
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
#' @importFrom CoreGx cellInfo<-
#'
#' @describeIn ToxicoSet Update the cell line annotations
#'
#' @export
setReplaceMethod("cellInfo", signature = signature(object = "ToxicoSet", value = "data.frame"), function(object, value){
  if (is.null(rownames(value))){
    stop("Please provide the cell_id as rownames for the cell line annotations")
  }
  object <- callNextMethod(object, value)
  object
})

#' drugInfo Generic
#'
#' The generic for drugInfo method
#'
#' @examples
#' data(TGGATESsmall)
#' drugInfo <- drugInfo(TGGATESsmall)
#'
#' @param tSet A \code{ToxicoSet} object
#'
#' @return a \code{data.frame} with the drug annotations
setGeneric("drugInfo", function(tSet) standardGeneric("drugInfo"))
#' @describeIn ToxicoSet Returns the annotations for all the drugs tested in the ToxicoSet
#' @export
setMethod(drugInfo, "ToxicoSet", function(tSet){
  tSet@drug
})

#' drugInfo<- Generic
#'
#' Generic for drugInfo replace method
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

#' phenoInfo Generic
#'
#' Generic for phenoInfo method
#'
#' @examples
#' data(TGGATESsmall)
#' phenoInfo <- phenoInfo(TGGATESsmall, mDataType="rna")
#'
#' @param tSet A \code{ToxicoSet} object #'
#' @param cSet Parameter name for parent method inherited from CoreGx
#' @param mDataType \code{character} A string specifying the type of molecular
#'   data to list the phenotype information for.
#'
#' @return a \code{data.frame} with the experiment info
#'
#' @importFrom CoreGx phenoInfo
#'
#' @describeIn ToxicoSet Return the experiment info from the given type of molecular data in ToxicoSet
#'
#' @export
#'
setMethod("phenoInfo",
          signature = c("ToxicoSet", "character"),
          function(cSet=tSet, mDataType){
            callNextMethod(cSet, mDataType)
          })

#' phenoInfo<- Generic
#'
#' Generic for phenoInfo replace method
#'
#' @examples
#' data(TGGATESsmall)
#' phenoInfo(TGGATESsmall, mDataType="rna") <- phenoInfo(TGGATESsmall, mDataType="rna")
#'
#' @param object A \code{ToxicoSet} object.
#' @param value A \code{data.frame} of replacement values.
#' @param mDataType A \code{character} with the type of molecular data to return/update
#'
#' @return The updated \code{ToxicoSet}
#'
#' @importMethodsFrom CoreGx phenoInfo<-
#'
#' @describeIn ToxicoSet Update the the given type of molecular data experiment info in the ToxicoSet
#'
#' @export
#'
setReplaceMethod("phenoInfo", signature = signature(object="ToxicoSet", mDataType ="character", value="data.frame"), function(object, mDataType, value) {
  object <- callNextMethod(object, mDataType, value)
  object
})

#' molecularProfiles Generic
#'
#' Generic for molecularProfiles method
#'
#' @examples
#' data(TGGATESsmall)
#' TGGATES_mProf <- molecularProfiles(TGGATESsmall, "rna")[1:10,]
#'
#' @param object A \code{ToxicoSet} object.
#' @param value A \code{character} vector of replacement values.
#' @param mDataType \code{character} A string specifying the type of molecular #'   data to list the phenotype information for.
#'
#' @describeIn ToxicoSet Return the given type of molecular data from the ToxicoSet
#'
#' @importMethodsFrom CoreGx molecularProfiles
#'
#' @export
setMethod("molecularProfiles",
          signature("ToxicoSet", "character"),
          function(cSet=tSet, mDataType){
            callNextMethod(cSet, mDataType)
          })

#' molecularProfiles<- Generic
#'
#' Generic for molecularProfiles replace method
#'
#' @examples
#' molecularProfiles(TGGATESsmall, "rna") <- molecularProfiles(TGGATESsmall, "rna")
#'
#' @param object A \code{ToxicoSet} object.
#' @param value A \code{matrix} of replacement values.
#' @param mDataType \code{character} A string specifying the type of molecular
#'   data to list the phenotype information for.
#'
#' @return Updated \code{ToxicoSet}
#'
#' @importMethodsFrom CoreGx molecularProfiles<-
#'
#' @describeIn ToxicoSet Update the given type of molecular data from the ToxicoSet
#'
#' @export
setReplaceMethod("molecularProfiles", signature = signature(object = "ToxicoSet", mDataType = "character", value = "matrix"), function(object, mDataType, value) {
  callNextMethod(object, mDataType, value)
})

#' Generic for featureInfo method
#'
#' @examples
#' data(TGGATESsmall)
#' featureInfo <- featureInfo(TGGATESsmall, "rna")[1:10,]
#'
#' @param tSet A \code{ToxicoSet} object #'
#' @param cSet Parameter name for parent method inherited from CoreGx
#' @param mDataType \code{character} A string specifying the type of molecular
#'   data to list the phenotype information for.
#'
#' @describeIn ToxicoSet Return the feature info for the given molecular data
#'
#' @importFrom CoreGx featureInfo
#'
#' @export
setMethod("featureInfo",
          signature("ToxicoSet", "character"),
          function(cSet=tSet, mDataType){
            callNextMethod(cSet, mDataType)
          })

#' featureInfo<- Generic
#'
#' Generic for featureInfo replace method
#'
#' @examples
#' data(TGGATESsmall)
#' featureInfo(TGGATESsmall, "rna") <- featureInfo(TGGATESsmall, "rna")
#'
#'
#' @param object A \code{ToxicoSet} object
#' @param value A \code{data.frame} of replacement values
#' @param mDataType \code{character} A string specifying the type of molecular
#'   data
#'
#'
#' @return Updated \code{ToxicoSet}
#'
#' @describeIn ToxicoSet Replace the gene info for the molecular data
#'
#' @importMethodsFrom CoreGx featureInfo<-
#'
#' @export
setReplaceMethod("featureInfo", signature = signature(object="ToxicoSet", mDataType ="character",value="data.frame"), function(object, mDataType, value){

  # if(mDataType %in% names(object@molecularProfiles)){Biobase::fData(object@molecularProfiles[[mDataType]]) <- value}
  object <- callNextMethod(object, mDataType, value)
  object
})

#' sensitivityInfo Generic
#'
#' Generic for sensitivityInfo method
#'
#' @examples
#' sensInf<- sensitivityInfo(TGGATESsmall)[1:10,]
#'
#' @param tSet A \code{ToxicoSet} object
#' @param cSet Parameter name for parent method inherited from CoreGx
#'
#' @return a \code{data.frame} with the experiment info
#'
#' @describeIn ToxicoSet Return the drug dose sensitivity experiment info
#'
#' @importMethodsFrom CoreGx sensitivityInfo
#'
#' @export
setMethod(sensitivityInfo,
          "ToxicoSet",
          function(cSet=tSet){
            callNextMethod(cSet)
          })

#' sensitivityInfo<- Generic
#'
#' A generic for the sensitivityInfo replacement method
#'
#'
#' @examples
#' data(TGGATESsmall)
#' sensitivityInfo(TGGATESsmall) <- sensitivityInfo(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#' @param value A \code{data.frame} of replacement values
#'
#' @return Updated \code{ToxicoSet}
#'setGeneric("sensitivityInfo<-", function(object, value) standardGeneric("sensitivityInfo<-"))
#' @importMethodsFrom CoreGx sensitivityInfo<-
#'
#' @describeIn ToxicoSet Update the sensitivity experiment info
#'
#' @export
setReplaceMethod("sensitivityInfo", signature = signature(object="ToxicoSet",value="data.frame"), function(object, value){
  # object@sensitivity$info <- value
  object <- callNextMethod(object, value)
  object
})


#' sensitivityProfiles Generic
#'
#' Generic for sensitivityProfiles method
#'
#' @examples
#' data(TGGATESsmall)
#' sensProf <- sensitivityProfiles(TGGATESsmall)
#'
#' @param tSet A \code{ToxicoSet} object
#' @param cSet Parameter name for parent method inherited from CoreGx
#'
#' @return a \code{data.frame} with the experiment info
#'setGeneric("sensitivityProfiles", function(tSet) standardGeneric("sensitivityProfiles"))
#'
#' @describeIn ToxicoSet Return the phenotypic data for the drug dose sensitivity
#'
#' @importFrom CoreGx sensitivityProfiles
#'
#' @export
setMethod(sensitivityProfiles,
          "ToxicoSet",
          function(cSet=tSet){
            callNextMethod(cSet)
          })

#' sensitivityProfiles<- Generic
#'
#' A generic for the sensitivityProfiles replacement method
#'
#' @examples
#' sensitivityProfiles(TGGATESsmall) <- sensitivityProfiles(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#' @param value A \code{data.frame} of replacement values
#'
#' @return Updated \code{ToxicoSet}
#'setGeneric("sensitivityProfiles<-", function(object, value) standardGeneric("sensitivityProfiles<-"))
#' @importFrom CoreGx sensitivityProfiles<-
#' @describeIn ToxicoSet Update the phenotypic data for the drug dose
#'   sensitivity
#' @export
setReplaceMethod("sensitivityProfiles", signature = signature(object="ToxicoSet",value="data.frame"), function(object, value){
  object <- callNextMethod(object, value)
  object
})
#' @describeIn ToxicoSet Update the phenotypic data for the drug dose
#'   sensitivity
#'
## TODO:: Find out how to document overloaded methods (to include multiple parameter types)
#' @export
setReplaceMethod("sensitivityProfiles", signature = signature(object="ToxicoSet",value="matrix"), function(object, value){
  object <- callNextMethod(object, value)
  object
})

#' sensitivityMeasures Generic
#'
#' A generic for the sensitivityMeasures  method
#'
#' @examples
#' sensitivityMeasures(TGGATESsmall)
#'
#' @param tSet A \code{ToxicoSet} object
#' @param cSet Parameter name for parent method inherited from CoreGx
#'
#' @return A \code{character} vector of all the available sensitivity measures
#'
#' @describeIn ToxicoSet Returns the available sensitivity profile
#'   summaries, for example, whether there are IC50 values available
#'
#' @importFrom CoreGx sensitivityMeasures
#'
#' @export
setMethod(sensitivityMeasures,
          "ToxicoSet",
          function(cSet=tSet){
            callNextMethod(cSet)
          })

#' drugNames Generic
#'
#' A generic for the drugNames method
#'
#' @examples
#' data(TGGATESsmall)
#' drugName <- drugNames(TGGATESsmall)[1:10]
#'
#' @param tSet A \code{ToxicoSet} object from which to retrieve the included
#'   drug names
#'
#' @return A vector of the drug names used in the ToxicoSet
setGeneric("drugNames", function(tSet) standardGeneric("drugNames"))
#' @describeIn ToxicoSet Return the names of the drugs used in the ToxicoSet
#' @export
setMethod(drugNames,
          "ToxicoSet",
          function(tSet){
            rownames(drugInfo(tSet))
          })

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
setReplaceMethod("drugNames", signature = signature(object="ToxicoSet",value="character"), function(object, value){
  object <- updateDrugId(object, value)
  return(object)
})

#' cellNames Generic
#'
#' A generic for the cellNames method
#'
#' @examples
#' cellNames(TGGATESsmall)
#'
#' @param tSet A \code{ToxicoSet} object
#' @param cSet Parameter name for parent method inherited from CoreGx
#'
#' @return A vector of the cell names used in the ToxicoSet
# setGeneric("cellNames", function(tSet) standardGeneric("cellNames"))
#'
#' @describeIn ToxicoSet Return the cell names used in the dataset
#'
#' @importFrom CoreGx cellNames
#'
#' @export
setMethod("cellNames",
          "ToxicoSet",
          function(cSet=tSet){
            callNextMethod(cSet)
          })

#' cellNames<- Generic
#'
#' A generic for the cellNames replacement method
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
#' @importMethodsFrom CoreGx cellNames<-
#'
#' @describeIn ToxicoSet Update the cell names used in the dataset
#'
#' @export
setReplaceMethod("cellNames", signature = signature(object="ToxicoSet",value="character"), function(object, value){
  object <- callNextMethod(object, value)
  return(object)
})

#' fNames Generic
#'
#' A generic for the fNames method
#'
#' @examples
#' fNames(TGGATESsmall, "rna")[1:10]
#'
#' @param tSet A \code{ToxicoSet} object
#' @param cSet Parameter name for parent method inherited from CoreGx
#' @param mDataType \code{character} A string specifying the type of molecular
#'   data to list the phenotype information for.
#'
#' @return A \code{character} vector of the feature names
#'
#' @describeIn ToxicoSet Return the feature names used in the dataset
#'
#' @importFrom CoreGx fNames
#'
#' @export
setMethod("fNames",
          signature("ToxicoSet", "character"),
          function(cSet=tSet, mDataType){
            callNextMethod(cSet, mDataType)
          })

###TODO:: Define this method in CoreGx and import it; doesn't work here because no generic defined
# fNames<- Generic
#
# A generic for the feature name replacement method
#
#@examples
# data(TGGATESsmall)
# cellNames(TGGATESsmall) <- cellNames(TGGATESsmall)
#
# @param mDataType \code{character} A string specifying the type of molecular #'   data to list the phenotype information for.
# @param value A \code{character} vector of the new feature names
# @return Updated \code{ToxicoSet}
# @describeIn ToxicoSet Update the feature names used in the dataset
# @export
#setReplaceMethod("fNames", signature = signature(object="ToxicoSet",value="character",mDataType="character"), function(object, value){
#      rownames(featureInfo(object, mDataType)) <- value
#})

#' dateCreated Generic
#'
#' A generic for the dateCreated method
#'
#' @examples
#' dateCreated(TGGATESsmall)
#'
#' @param tSet A \code{ToxicoSet} object
#' @param cSet Parameter name for parent method inherited from CoreGx
#'
#' @return The date the ToxicoSet was created
#'
#' @describeIn ToxicoSet Return the date the ToxicoSet was created
#'
#' @importFrom CoreGx dateCreated
#'
#' @export
setMethod(dateCreated,
          signature = c("ToxicoSet"),
          function(cSet=tSet) {
            callNextMethod(cSet)
          })


#' tSetName Generic
#'
#' A generic for the tSetName method
#'
#' @examples
#' tSetName <- cSetName
#' tSetName(TGGATESsmall)
#'
#' @param tSet A \code{ToxicoSet} object
#' @param cSet Parameter name for parent method inherited from CoreGx
#'
#' @return The name of the ToxicoSet
#'
#' @describeIn ToxicoSet Return the name of the ToxicoSet
#'
#' @importFrom CoreGx cSetName
#'
#' @export
setMethod("cSetName",
          signature = c("ToxicoSet"),
          function(cSet=tSet){
            callNextMethod(cSet)
          })
tSetName <- cSetName

#' pertNumber Generic
#'
#' A generic for the pertNumber method
#'
#' @examples
#' pertNumber(TGGATESsmall)
#'
#' @param tSet A \code{ToxicoSet} object
#' @param cSet Parameter name for parent method inherited from CoreGx
#'
#' @return A 3D \code{array} with the number of perturbation experiments per radiation type and cell line, and data type
# setGeneric("pertNumber", function(tSet) standardGeneric("pertNumber"))
#'
#' @describeIn ToxicoSet Return the summary of available perturbation
#'   experiments
#'
#' @importMethodsFrom CoreGx pertNumber
#'
#' @export
setMethod(pertNumber,
          "ToxicoSet",
          function(cSet=tSet){
            callNextMethod(cSet)
          })


#' sensNumber Generic
#'
#' A generic for the sensNumber method
#'
#' @examples
#' sensNumber(TGGATESsmall)
#'
#' @param tSet A \code{ToxicoSet} object
#' @param cSet Parameter name for parent method inherited from CoreGx
#'
#' @return A \code{data.frame} with the number of sensitivity experiments per drug and cell line
#'
#' @describeIn ToxicoSet Return the summary of available sensitivity
#'   experiments
#'
#' @importFrom CoreGx sensNumber
#'
#' @export
setMethod(sensNumber,
          "ToxicoSet",
          function(cSet=tSet){
            callNextMethod(cSet)
          })

#' pertNumber<- Generic
#'
#' A generic for the pertNumber method
#'
#' @examples
#' pertNumber(TGGATESsmall) <- pertNumber(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object to modify
#' @param value An \code{array} of replacement values
## TODO:: Remove duplicate param names from CoreGx; this documentation is incorrect as is
# @param value A new 3D \code{array} with the number of perturbation experiments per drug and cell line, and data type
#'
#' @return The updated \code{ToxicoSet}
#'
#' @importMethodsFrom CoreGx pertNumber<-
#'
#' @describeIn ToxicoSet Update the summary of available perturbation
#'   experiments
#'
#' @export
setReplaceMethod('pertNumber', signature = signature(object="ToxicoSet",value="array"), function(object, value){
  object <- callNextMethod(object, value)
  object

})

#' sensNumber<- Generic
#'
#' A generic for the sensNumber method
#'
#' @examples
#' sensNumber(TGGATESsmall) <- sensNumber(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object to modify
#' @param value A \code{matrix} of replacement values
#'
#' @return The updated \code{ToxicoSet}
#'
#' @importMethodsFrom CoreGx sensNumber<-
#'
#' @describeIn ToxicoSet Update the summary of available sensitivity
#'   experiments
#'
#' @export
setReplaceMethod('sensNumber', signature = signature(object="ToxicoSet",value="matrix"), function(object, value){
  object <- callNextMethod(object, value)
  object

})

#' Show a ToxicoSet
#'
#' @param object A \code{ToxicoSet} object to print a summary for
#'
#' @examples
#' TGGATESsmall
#'
#' @return Prints the ToxicoSet object to the output stream, and returns invisible NULL.
#' @export
setMethod("show", signature=signature(object="ToxicoSet"),
          function(object) {
            cat("Name: ", tSetName(object), "\n")
            cat("Date Created: ", dateCreated(object), "\n")
            cat("Number of cell lines: ", nrow(cellInfo(object)), "\n")
            cat("Number of drugs: ", nrow(drugInfo(object)), "\n")
            if("dna" %in% names(object@molecularProfiles)){cat("DNA: \n");cat("\tDim: ", dim(molecularProfiles(object, mDataType="dna")), "\n")}
            if("rna" %in% names(object@molecularProfiles)){cat("RNA: \n");cat("\tDim: ", dim(molecularProfiles(object, mDataType="rna")), "\n")}
            if("rnaseq" %in% names(object@molecularProfiles)){cat("RNASeq: \n");cat("\tDim: ", dim(molecularProfiles(object, mDataType="rnaseq")), "\n")}
            if("snp" %in% names(object@molecularProfiles)){cat("SNP: \n");cat("\tDim: ", dim(molecularProfiles(object, mDataType="snp")), "\n")}
            if("cnv" %in% names(object@molecularProfiles)){cat("CNV: \n");cat("\tDim: ", dim(molecularProfiles(object, mDataType="cnv")), "\n")}
            cat("Drug pertubation: \n")
            cat("\tPlease look at pertNumber(tSet) to determine number of experiments for each drug-cell combination.\n")
            cat("Drug sensitivity: \n")
            cat("\tNumber of Experiments: ",nrow(sensitivityInfo(object)),"\n")
            cat("\tPlease look at sensNumber(tSet) to determine number of experiments for each drug-cell combination.\n")
          })

#' mDataNames
#'
#' Returns the molecular data names for the ToxicoSet.
#'
#' @examples
#' mDataNames(TGGATESsmall)
#'
#' @param tSet A \code{ToxicoSet} object
#' @param cSet Parameter name for parent method inherited from CoreGx
#'
#' @return Vector of names of the molecular data types
# Imports generic
#' @importFrom CoreGx mDataNames cellInfo<-
#' @export
setMethod(
  "mDataNames",
  signature = c("ToxicoSet"),
  definition = function(cSet=tSet){
    callNextMethod(cSet)
  }
)

#'`[`
#'
#' @examples
#' tSet <- TGGATESsmall[cellNames(TGGATESsmall), drugNames(TGGATESsmall)[1:3]]
#'
#'@param x tSet
#'@param i Cell lines to keep in tSet
#'@param j Drugs to keep in tSet
#'@param ... further arguments
#'@param drop A boolean flag of whether to drop single dimensions or not
#'@return Returns the subsetted tSet
#'@export
setMethod(`[`, "ToxicoSet", function(x, i, j, ..., drop = FALSE){
  if(is.character(i)&&is.character(j)){
    return(subsetTo(x, cells=i, drugs=j,  molecular.data.cells=i))
  }
  else if(is.numeric(i) && is.numeric(j) && (as.integer(i)==i) && (as.integer(j)==j)){
    return(subsetTo(x, cells=cellNames(x)[i], drugs=drugNames(x)[j],  molecular.data.cells=cellNames(x)[i]))
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
setMethod("dim", signature=signature(x="ToxicoSet"), function(x){

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
#' @param tSet A \code{ToxicoSet} to be subsetted
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
#' @param ... Other arguments passed by other function within the package
#' @return A ToxicoSet with only the selected drugs and cells
#' @importFrom CoreGx unionList
#' @export
## TODO:: Include dose parmater to subset on
subsetTo <- function(tSet, cell_lines = NULL,
                     drugs=NULL,
                     molecular.data.cells=NULL,
                     duration=NULL, features=NULL, ...
                     ) {
  ## TODO:: Remove this or add it to the function parameters?
  drop = FALSE

  ####
  # PARSING ARGUMENTS
  ####
  adArgs = list(...)
  if ("exps" %in% names(adArgs)) {
    exps <- adArgs[["exps"]]
    if(class(exps) == "data.frame"){
      exps2 <- exps[[cSetName(tSet)]]
      names(exps2) <- rownames(exps)
      exps <- exps2
    } else{
      exps <- exps[[cSetName(tSet)]]
    }
  }else {
    exps <- NULL
  }

  if ("dose" %in% names(adArgs)) {
    ## TODO:: Add subsetting on dose
    stop("Due to the structure of tSets, subsetting on dose can only be done on specific slots - not on the entire tSet")
  }

  ## MISSING VALUE HANDLING FOR PARAMETERS
  # Get named list of defualt values for missing parameters
  argDefaultList <-
    paramMissingHandler(funName = "subsetTo", tSet = tSet,
                        drugs = drugs, cell_lines = cell_lines, features = features,
                        duration = duration)
  # Assign any missing parameter default values to function environment
  if (length(argDefaultList) > 0) {
    for (idx in seq_along(argDefaultList)) {
      assign(names(argDefaultList)[idx], argDefaultList[[idx]])
    }
  }

  # ERROR HANDLING FOR PARAMETERS
  paramErrorChecker(funName = "subsetTo", tSet = tSet,
                    cell_lines = cell_lines,
                    drugs = drugs, features = features,
                    duration = duration)

  ######
  # SUBSETTING MOLECULAR PROFILES SLOT
  ######
  ### TODO:: implement strict subsetting at this level!!!!

  ### the function missing does not work as expected in the context below, because the arguments are passed to the anonymous
  ### function in lapply, so it does not recognize them as missing
  tSet@molecularProfiles <- lapply(tSet@molecularProfiles, function(eset, cell_lines, drugs, molecular.data.cells, duration, features){

    if (!is.null(features)) {
      eset <- eset[which(Biobase::featureNames(eset) %in% features), ]
    }

    molecular.data.type <- ifelse(length(grep("rna", Biobase::annotation(eset)) > 0), "rna", Biobase::annotation(eset))
    if (length(grep(molecular.data.type, names(molecular.data.cells))) > 0) {
      cell_lines <- molecular.data.cells[[molecular.data.type]]
    }
    column_indices <- NULL

    if (length(cell_lines) == 0 && length(drugs) == 0) {
      column_indices <- seq_len(ncol(eset))
    }
    if (length(cell_lines) == 0 && tSet@datasetType == "sensitivity") {
      column_indices <- seq_len(ncol(eset))
    }

    # Selecting indices which match the cells argument
    cell_line_index <- NULL
    if (length(cell_lines) != 0) {
      if (!all(cell_lines %in% cellNames(tSet))) {
        stop("Some of the cell names passed to function did not match to names in the PharmacoSet. Please ensure you are using cell names as returned by the cellNames function")
      }
      cell_line_index <- which(Biobase::pData(eset)[["cellid"]] %in% cell_lines)
    }

    # Selecting indexes which match drugs arguement
    drugs_index <- NULL
    if (tSet@datasetType == "perturbation" || tSet@datasetType == "both") {
      if (length(drugs) != 0) {
        if (!all(drugs %in% drugNames(tSet))){
          stop("Some of the drug names passed to function did not match to names in the ToxicoSet Please ensure you are using drug names as returned by the drugNames function")
        }
        drugs_index <- which(Biobase::pData(eset)[["drugid"]] %in% drugs)
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
    ## TODO:: Determine if this works for other eSet data types
    if (!is.null(duration)){
      if (all(!(duration %in% unique(Biobase::pData(eset[,column_indices])$duration)))) {
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
      duration_indices <- which(Biobase::pData(eset)$duration %in% duration)
      column_indices <- intersect(column_indices, duration_indices)
    }

    row_indices <- seq_len(nrow(Biobase::exprs(eset)))

    # Final eSet
    eset <- eset[row_indices, column_indices]
    return(eset)

  }, cell_lines = cell_lines, drugs = drugs, molecular.data.cells = molecular.data.cells, duration = duration, features = features)


  ######
  # SUBSET SENSITIVITY SLOT
  ######
  # Logic if any "..." arguments are passed to subsetTo
  if ((tSet@datasetType == "sensitivity" | tSet@datasetType == "both") & length(exps) != 0) {
    tSet@sensitivity$info <- tSet@sensitivity$info[exps, , drop=drop]
    rownames(tSet@sensitivity$info) <- names(exps)
    if (length(tSet@sensitivity$raw) > 0) {
      tSet@sensitivity$raw <- tSet@sensitivity$raw[exps, , , drop=drop]
      dimnames(tSet@sensitivity$raw)[[1]] <- names(exps)
    }
    tSet@sensitivity$profiles <- tSet@sensitivity$profiles[exps, , drop=drop]
    rownames(tSet@sensitivity$profiles) <- names(exps)

    tSet@sensitivity$n <- .summarizeSensitivityNumbers(tSet)
  }
  # Logic if drug or cell parameters are passed to subsetTo
  else if (
    (tSet@datasetType == "sensitivity" | tSet@datasetType == "both") &
    (length(drugs) != 0 | length(cell_lines) != 0 | !is.null(duration) )
  ) {

    drugs_index <- which(sensitivityInfo(tSet)[, "drugid"] %in% drugs)
    cell_line_index <- which(sensitivityInfo(tSet)[,"cellid"] %in% cell_lines)
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
          row_indices <- seq_len(nrow(sensitivityInfo(tSet)))
        }
      }
    }
    # LOGIC TO SUBSET BASED ON DURATION
    if(!is.null(duration)){
      if(all(!(duration %in% unique(sensitivityInfo(tSet)[row_indices,]$duration_h)))) {
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
            duration, " in the tSet."
          ))
        }
      }
      duration_indices <- which(sensitivityInfo(tSet)$duration_h %in% duration)
      row_indices <- intersect(row_indices, duration_indices)
    }
    sensItemNames <- names(tSet@sensitivity)
    sensitivityVals <-
      lapply(sensItemNames, function(sensItemName, drop){
        if (sensItemName == "n") {
          sensItem <- tSet@sensitivity[[sensItemName]]
          if (!is.null(cell_lines)) {
            sensItem[which(rownames(sensItem) %in% cell_lines), which(colnames(sensItem) %in% drugs), drop = drop]
          } else {
            sensItem[ , which(colnames(sensItem) %in% drugs), drop = drop]
          }
        } else {
          sensItem <- tSet@sensitivity[[sensItemName]]
          if (length(dim(sensItem)) == 3) {
            sensItem[row_indices, , , drop = drop]
          } else {
            sensItem[row_indices, , drop = drop]
          }
        }
      }, drop = drop)
    names(sensitivityVals) <- sensItemNames
    tSet@sensitivity <- sensitivityVals
  }

  #####
  # SUBSET DRUG SLOT
  #####
  if (length(drugs) == 0) {
    if (tSet@datasetType == "sensitivity" | tSet@datasetType == "both"){
      drugs <- unique(sensitivityInfo(tSet)[["drugid"]])
    }
    if(tSet@datasetType == "perturbation" | tSet@datasetType == "both"){
      drugs <- union(drugs, na.omit(unionList(lapply(tSet@molecularProfiles, function(eSet){unique(Biobase::pData(eSet)[["drugid"]])}))))
    }
  }
  #####
  # SUBSET CELLS SLOT
  #####
  if (length(cell_lines) == 0) {
    celll.lines <- union(cell_lines, na.omit(unionList(lapply(tSet@molecularProfiles, function(eSet){unique(Biobase::pData(eSet)[["cellid"]])}))))
    if (tSet@datasetType == "sensitivity" | tSet@datasetType == "both"){
      cell_lines <- union(cell_lines, sensitivityInfo(tSet)[["cellid"]])
    }
  }
  #####
  # ASSIGN SUBSETS BACK TO TOXICOSET OBJECT
  #####
  drugInfo(tSet) <- drugInfo(tSet)[drugs , , drop=drop]
  cellInfo(tSet) <- cellInfo(tSet)[cell_lines , , drop=drop]
  tSet@curation$drug <- tSet@curation$drug[drugs , , drop=drop]
  tSet@curation$cell <- tSet@curation$cell[cell_lines , , drop=drop]
  tSet@curation$tissue <- tSet@curation$tissue[cell_lines , , drop=drop]
  return(tSet)
}

#
# END SUBSET TO FUNCTION
#







### TODO:: Add updating of sensitivity Number tables
#' A function to update cell ids
#' @examples
#' data(TGGATESsmall)
#' updateCellId(TGGATESsmall, new.ids = cellNames(TGGATESsmall))
#' @param tSet [object] A ToxicoSet object to be updates
#' @param new.ids [character] A character vector of ids to update with
#' @keywords internal
#' @export
updateCellId <- function(tSet, new.ids = vector("character")){

  if (length(new.ids)!=nrow(cellInfo(tSet))){
    stop("Wrong number of cell identifiers")
  }

  if(tSet@datasetType=="sensitivity"|tSet@datasetType=="both"){
    myx <- match(sensitivityInfo(tSet)[,"cellid"],rownames(cellInfo(tSet)))
    sensitivityInfo(tSet)[,"cellid"] <- new.ids[myx]

  }


  tSet@molecularProfiles <- lapply(tSet@molecularProfiles, function(eset){

    myx <- match(Biobase::pData(eset)[["cellid"]], rownames(cellInfo(tSet)))
    Biobase::pData(eset)[["cellid"]]  <- new.ids[myx]
    return(eset)
  })

  if(any(duplicated(new.ids))){
    warning("Duplicated ids passed to updateCellId. Merging old ids into the same identifier")

    if(ncol(sensNumber(tSet))>0){
      sensMatch <- match(rownames(sensNumber(tSet)), rownames(cellInfo(tSet)))
    }
    if(dim(pertNumber(tSet))[[2]]>0){
      pertMatch <- match(dimnames(pertNumber(tSet))[[1]], rownames(cellInfo(tSet)))
    }
    curMatch <- match(rownames(tSet@curation$cell),rownames(cellInfo(tSet)))

    duplId <- unique(new.ids[duplicated(new.ids)])
    for(id in duplId){

      if (ncol(sensNumber(tSet))>0){
        myx <- which(new.ids[sensMatch] == id)
        sensNumber(tSet)[myx[1],] <- apply(sensNumber(tSet)[myx,], 2, sum)
        sensNumber(tSet) <- sensNumber(tSet)[-myx[-1],]
        # sensMatch <- sensMatch[-myx[-1]]
      }
      if (dim(pertNumber(tSet))[[1]]>0){
        myx <- which(new.ids[pertMatch] == id)
        pertNumber(tSet)[myx[1],,] <- apply(pertNumber(tSet)[myx,,], c(1,3), sum)
        pertNumber(tSet) <- pertNumber(tSet)[-myx[-1],,]
        # pertMatch <- pertMatch[-myx[-1]]
      }

      myx <- which(new.ids[curMatch] == id)
      tSet@curation$cell[myx[1],] <- apply(tSet@curation$cell[myx,], 2, paste, collapse="///")
      tSet@curation$cell <- tSet@curation$cell[-myx[-1],]
      tSet@curation$tissue[myx[1],] <- apply(tSet@curation$tissue[myx,], 2, paste, collapse="///")
      tSet@curation$tissue <- tSet@curation$tissue[-myx[-1],]
      # curMatch <- curMatch[-myx[-1]]

      myx <- which(new.ids == id)
      cellInfo(tSet)[myx[1],] <- apply(cellInfo(tSet)[myx,], 2, paste, collapse="///")
      cellInfo(tSet) <- cellInfo(tSet)[-myx[-1],]
      new.ids <- new.ids[-myx[-1]]
      if(ncol(sensNumber(tSet))>0){
        sensMatch <- match(rownames(sensNumber(tSet)), rownames(cellInfo(tSet)))
      }
      if(dim(pertNumber(tSet))[[1]]>0){
        pertMatch <- match(dimnames(pertNumber(tSet))[[1]], rownames(cellInfo(tSet)))
      }
      curMatch <- match(rownames(tSet@curation$cell),rownames(cellInfo(tSet)))
    }
  } else {
    if (dim(pertNumber(tSet))[[1]]>0){
      pertMatch <- match(dimnames(pertNumber(tSet))[[1]], rownames(cellInfo(tSet)))
    }
    if (ncol(sensNumber(tSet))>0){
      sensMatch <- match(rownames(sensNumber(tSet)), rownames(cellInfo(tSet)))
    }
    curMatch <- match(rownames(tSet@curation$cell),rownames(cellInfo(tSet)))
  }

  if (dim(pertNumber(tSet))[[1]]>0){
    dimnames(pertNumber(tSet))[[1]] <- new.ids[pertMatch]
  }
  if (ncol(sensNumber(tSet))>0){
    rownames(sensNumber(tSet)) <- new.ids[sensMatch]
  }
  rownames(tSet@curation$cell) <- new.ids[curMatch]
  rownames(tSet@curation$tissue) <- new.ids[curMatch]
  rownames(cellInfo(tSet)) <- new.ids

  return(tSet)
}





# updateFeatureNames <- function(tSet, new.ids = vector("character")){
#
#   if (length(new.ids)!=nrow(cellInfo(tSet))){
#     stop("Wrong number of cell identifiers")
#   }
#
#   if(tSet@datasetType=="sensitivity"|tSet@datasetType=="both"){
#     myx <- match(sensitivityInfo(tSet)[,"cellid"],rownames(cellInfo(tSet)))
#     sensitivityInfo(tSet)[,"cellid"] <- new.ids[myx]
#
#   }
#
#   tSet@molecularProfiles <- lapply(tSet@molecularProfiles, function(eset){
#
#     myx <- match(pData(eset)[["cellid"]],rownames(cellInfo(tSet)))
#     pData(eset)[["cellid"]]  <- new.ids[myx]
#     return(eset)
#       })
#   myx <- match(rownames(tSet@curation$cell),rownames(cellInfo(tSet)))
#   rownames(tSet@curation$cell) <- new.ids[myx]
#   rownames(tSet@curation$tissue) <- new.ids[myx]
#   if (dim(pertNumber(tSet))[[1]]>0){
#     myx <- match(dimnames(pertNumber(tSet))[[1]], rownames(cellInfo(tSet)))
#     dimnames(pertNumber(tSet))[[1]] <- new.ids[myx]
#   }
#   if (nrow(sensNumber(tSet))>0){
#     myx <- match(rownames(sensNumber(tSet)), rownames(cellInfo(tSet)))
#     rownames(sensNumber(tSet)) <- new.ids[myx]
#   }
#   rownames(cellInfo(tSet)) <- new.ids
#   return(tSet)
#
# }

### TODO:: Add updating of sensitivity Number tables
#' A function to update drug ids
#' @examples
#' data(TGGATESsmall)
#' updateDrugId(TGGATESsmall, new.ids = drugNames(TGGATESsmall))
#' @param tSet [object] A ToxicoSet object to be updates
#' @param new.ids [character] A character vector of ids to update with
#' @keywords internal
#' @export
updateDrugId <- function(tSet, new.ids = vector("character")){

  if (length(new.ids)!= nrow(drugInfo(tSet))){
    stop("Wrong number of drug identifiers")
  }

  if(tSet@datasetType=="sensitivity"|tSet@datasetType=="both"){
    myx <- match(sensitivityInfo(tSet)[,"drugid"],rownames(drugInfo(tSet)))
    sensitivityInfo(tSet)[,"drugid"] <- new.ids[myx]

  }
  if(tSet@datasetType=="perturbation"|tSet@datasetType=="both"){
    tSet@molecularProfiles <- lapply(tSet@molecularProfiles, function(eset){

      myx <- match(Biobase::pData(eset)[["drugid"]],rownames(drugInfo(tSet)))
      Biobase::pData(eset)[["drugid"]]  <- new.ids[myx]
      return(eset)
    })
  }


  if(any(duplicated(new.ids))){
    warning("Duplicated ids passed to updateDrugId. Merging old ids into the same identifier")

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

  ## unique drug identifiers
  # drugn <- sort(unique(tSet@sensitivity$info[ , "drugid"]))

  ## consider all drugs
  drugn <- rownames(tSet@drug)

  ## unique drug identifiers
  # celln <- sort(unique(tSet@sensitivity$info[ , "cellid"]))

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
    if (nrow(Biobase::pData(tSet@molecularProfiles[[i]])) > 0 && all(is.element(c("cellid", "drugid"), colnames(Biobase::pData(tSet@molecularProfiles[[i]]))))) {
      tt <- table(Biobase::pData(tSet@molecularProfiles[[i]])[ , "cellid"], Biobase::pData(tSet@molecularProfiles[[i]])[ , "drugid"])
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
    for( i in 1:length(tSet@molecularProfiles)) {
      profile <- tSet@molecularProfiles[[i]]
      nn <- names(tSet@molecularProfiles)[i]
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
        if(!all(Biobase::pData(profile)[,"cellid"] %in% rownames(tSet@cell))) {
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
        print("rownames of cell slot should be curated cell ids")
      }
    } else {
      print("unique.cellid which is curated cell id across data set should be a column of cell curation slot")
    }

    if(length(intersect(rownames(tSet@curation$cell), rownames(tSet@cell))) != nrow(tSet@cell)) {
      print("rownames of curation cell slot should be the same as cell slot (curated cell ids)")
    }

    if("unique.drugid" %in% colnames(tSet@curation$drug)) {
      if(length(intersect(tSet@curation$drug$unique.drugid, rownames(tSet@drug))) != nrow(tSet@drug)) {
        print("rownames of drug slot should be curated drug ids")
      }
    } else {
      print("unique.drugid which is curated drug id across data set should be a column of drug curation slot")
    }

    if(length(intersect(rownames(tSet@curation$cell), rownames(tSet@cell))) != nrow(tSet@cell)) {
      print("rownames of curation drug slot should be the same as drug slot (curated drug ids)")
    }

    if(class(tSet@cell) != "data.frame") {
      warning("cell slot class type should be dataframe")
    }
    if(class(tSet@drug) != "data.frame") {
      warning("drug slot class type should be dataframe")
    }
    if(tSet@datasetType %in% c("sensitivity", "both"))
    {
      if(class(tSet@sensitivity$info) != "data.frame") {
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
          print("not all the drugs in sensitivity data are in drug slot")
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
