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
.ToxicoSet <- setClass("ToxicoSet",
                       slots = list(drug="data.frame"),
                       contains="CoreSet")

### -------------------------------------------------------------------------
### Constructor -------------------------------------------------------------
### -------------------------------------------------------------------------

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

# Helper Functions --------------------------------------------------------

.summarizeSensitivityNumbers <- function(tSet) {

  if (datasetType(tSet) != "sensitivity" && datasetType(tSet) != "both") {
    stop ("Data type must be either sensitivity or both")
  }

  ## consider all drugs
  drugn <- rownames(tSet@drug)

  ## consider all cell lines
  celln <- rownames(tSet@cell)

  sensitivity.info <- matrix(0, nrow=length(celln), ncol=length(drugn),
                             dimnames=list(celln, drugn))
  drugids <- sensitivityInfo(tSet)[ , "drugid"]
  cellids <- sensitivityInfo(tSet)[ , "cellid"]
  cellids <- cellids[grep("///", drugids, invert=TRUE)]
  drugids <- drugids[grep("///", drugids, invert=TRUE)]


  tt <- table(cellids, drugids)
  sensitivity.info[rownames(tt), colnames(tt)] <- tt

  return(sensitivity.info)
}

.summarizePerturbationNumbers <- function(tSet) {

  if (datasetType(tSet) != "perturbation" && datasetType(tSet) != "both") {
    stop ("Data type must be either perturbation or both")
  }

  ## consider all drugs
  drugn <- rownames(tSet@drug)

  ## consider all cell lines
  celln <- rownames(tSet@cell)

  perturbation.info <- array(0, dim=c(length(celln), length(drugn), length(molecularProfilesSlot(tSet))), dimnames=list(celln, drugn, names((molecularProfilesSlot(tSet)))))

  for (i in seq_along(molecularProfilesSlot(tSet))) {
    if (nrow(SummarizedExperiment::colData(molecularProfilesSlot(tSet)[[i]])) > 0 &&
        all(
          is.element(c("cellid", "drugid"),
                     colnames(SummarizedExperiment::colData(molecularProfilesSlot(tSet)[[i]]))))) {
      tt <- table(SummarizedExperiment::colData(molecularProfilesSlot(tSet)[[i]])[ , "cellid"], SummarizedExperiment::colData(molecularProfilesSlot(tSet)[[i]])[ , "drugid"])
      perturbation.info[rownames(tt), colnames(tt), names(molecularProfilesSlot(tSet))[i]] <- tt
    }
  }

  return(perturbation.info)
}


### -------------------------------------------------------------------------
### Class Validity ----------------------------------------------------------
### -------------------------------------------------------------------------


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
#' @importFrom S4Vectors metadata
checkTSetStructure <-
  function(tSet, plotDist=FALSE, result.dir=".") {
    if(!file.exists(result.dir) && plotDist) { dir.create(result.dir, showWarnings=FALSE, recursive=TRUE) }
    for( i in seq_along(molecularProfilesSlot(tSet))) {
      profile <- molecularProfilesSlot(tSet)[[i]]
      if (is.null(names(metadata(profile))))
        stop(paste0("Please ensure all items in the metadata slot of
             SummarizedExperiments are named. Item ", i, " of molecualrProfiles
             does not have metadata names."))
      if (!("annotation" %in% names(metadata(profile))))
        stop(paste0("At minimum the SummarizedExperiments in molecularProfiles must contain
                a metadata item names 'annotation' specifying the molecular datatype
               the SummarizedExperiment contains! Item ", i, " of
               molecularProfilesis missing annotation metadata."))
      nn <- names(molecularProfilesSlot(tSet))[i]
      if(plotDist) {
        if (S4Vectors::metadata(profile)$annotation == "rna" || S4Vectors::metadata(profile)$annotation == "rnaseq") {
          pdf(file=file.path(result.dir, sprintf("%s.pdf", nn)))
          hist(SummarizedExperiment::assay(profile, 1), breaks = 100)
          dev.off()
        }
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
    if(datasetType(tSet) %in% c("sensitivity", "both"))
    {
      if(!is(sensitivityInfo(tSet), "data.frame")) {
        warning("sensitivity info slot class type should be dataframe")
      }
      if("cellid" %in% colnames(sensitivityInfo(tSet))) {
        if(!all(sensitivityInfo(tSet)[,"cellid"] %in% rownames(tSet@cell))) {
          warning("not all the cell lines in sensitivity data are in cell slot")
        }
      }else {
        warning("cellid does not exist in sensitivity info")
      }
      if("drugid" %in% colnames(sensitivityInfo(tSet))) {
        drug.ids <- unique(sensitivityInfo(tSet)[,"drugid"])
        drug.ids <- drug.ids[grep("///",drug.ids, invert=TRUE)]
        if(!all(drug.ids %in% rownames(tSet@drug))) {
          message("not all the drugs in sensitivity data are in drug slot")
        }
      }else {
        warning("drugid does not exist in sensitivity info")
      }

      if(any(!is.na(tSet@sensitivity$raw))) {
        if(!all(dimnames(tSet@sensitivity$raw)[[1]] %in% rownames(sensitivityInfo(tSet)))) {
          warning("For some experiments there is raw sensitivity data but no experimet information in sensitivity info")
        }
      }
      if(!all(rownames(tSet@sensitivity$profiles) %in% rownames(sensitivityInfo(tSet)))) {
        warning("For some experiments there is sensitivity profiles but no experimet information in sensitivity info")
      }
    }
  }



# -------------------------------------------------------------------------
# Method Definitions ------------------------------------------------------
# -------------------------------------------------------------------------


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
