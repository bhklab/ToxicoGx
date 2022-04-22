#' @include allGenerics.R
NULL

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
#' @slot sample A \code{data.frame} containing the annotations for all the cell
#'   lines profiled in the data set, across all data types
#' @slot treatment A \code{data.frame} containg the annotations for all the drugs
#'   profiled in the data set, across all data types
#' @slot treatmentResponse A \code{list} containing all the data for the sensitivity
#'   experiments, including \code{$info}, a \code{data.frame} containing the
#'   experimental info,\code{$raw} a 3D \code{array} containing raw data,
#'   \code{$profiles}, a \code{data.frame} containing sensitivity profiles
#'   statistics, and \code{$n}, a \code{data.frame} detailing the number of
#'   experiments for each cell-drug pair
#' @slot perturbation A \code{list} containting \code{$n}, a \code{data.frame}
#'   summarizing the available perturbation data,
#' @slot curation A \code{list} containing mappings for \code{$treatment},
#'   \code{sample}, \code{tissue} names  used in the data set to universal
#'   identifiers used between different ToxicoSet objects
#' @slot datasetType A \code{character} string of 'sensitivity',
#'   'perturbation', or both detailing what type of data can be found in the
#'   ToxicoSet, for proper processing of the data
#'
#' @return An object of the ToxicoSet class
#'
#' @importClassesFrom CoreGx CoreSet
.ToxicoSet <- setClass("ToxicoSet", contains="CoreSet")

## TODO:: implement .intern slot to hold arbitrary metadata about a tSet

## Variables for dynamic inheritted roxygen2 docs

.local_class="ToxicoSet"
.local_data="TGGATESsmall"
.local_sample="cell"

#### CoreGx inherited methods
####
#### Note: The raw documentation lives in CoreGx, see the functions called
#### in @eval tags for the content of the metaprogrammed roxygen2 docs.
####
#### See .parseToRoxygen method in utils-messages.R file of CoreGx to
#### create similar metaprogrammed docs.
####
#### Warning: for dynamic docs to work, you must set
#### Roxygen: list(markdown = TRUE, r6=FALSE)
#### in the DESCRPTION file!


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
#' @inheritParams CoreGx::CoreSet
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
                       sample=data.frame(),
                       treatment=data.frame(),
                       sensitivityInfo=data.frame(),
                       sensitivityRaw=array(dim = c(0,0,0)),
                       sensitivityProfiles=matrix(),
                       sensitivityN=matrix(nrow = 0, ncol=0),
                       perturbationN=array(NA, dim = c(0,0,0)),
                       curationTreatment=data.frame(),
                       curationSample = data.frame(),
                       curationTissue = data.frame(),
                       datasetType=c("sensitivity", "perturbation", "both"),
                       #sharedControls=FALSE,
                       verify = TRUE)
{
    # .Deprecated("ToxicoSet2", package=packageName(), msg="The ToxicoSet class is
    #     being redesigned. Please use the new constructor to ensure forwards
    #     compatibility with future releases! Old objects can be updated with
    #     the updateObject method.", old="ToxicoSet")

    cSet <- CoreGx::CoreSet(
        name=name,
        sample=sample,
        treatment=treatment,
        sensitivityInfo=sensitivityInfo,
        sensitivityRaw=sensitivityRaw,
        sensitivityProfiles=sensitivityProfiles,
        sensitivityN=sensitivityN,
        perturbationN=perturbationN,
        curationTreatment=curationTreatment,
        curationSample=curationSample,
        curationTissue=curationTissue,
        datasetType=datasetType,
        verify=verify
    )

    tSet  <- .ToxicoSet(
        annotation=cSet@annotation,
        molecularProfiles=cSet@olecularProfiles,
        sample=cSet@sample,
        treatment=cSet@treatment,
        datasetType=cSet@datasetTypes,
        treatmentResponse=cSet@treatmentResponse,
        perturbation=cSet@perturbation,
        curation=cSet@curation
    )
    if (verify) { checkTSetStructure(tSet)}
    if (length(sensitivityN) == 0 & datasetType %in% c("sensitivity", "both")) {
        sensNumber(tSet) <- .summarizeSensitivityNumbers(tSet)
    }
    if (length(perturbationN) == 0  & datasetType %in% c("perturbation", "both")) {
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
  drugn <- treatmentNames(tSet)

  ## consider all cell lines
  celln <- rownames(sampleInfo(tSet))

  sensitivity.info <- matrix(0, nrow=length(celln), ncol=length(drugn),
                             dimnames=list(celln, drugn))
  drugids <- sensitivityInfo(tSet)[, "treatmentid"]
  cellids <- sensitivityInfo(tSet)[, "sampleid"]
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
  drugn <- treatmentNames(tSet)

  ## consider all cell lines
  celln <- rownames(sampleInfo(tSet))

  perturbation.info <- array(0, dim=c(length(celln), length(drugn), length(molecularProfilesSlot(tSet))), dimnames=list(celln, drugn, names((molecularProfilesSlot(tSet)))))

  for (i in seq_along(molecularProfilesSlot(tSet))) {
    if (nrow(SummarizedExperiment::colData(molecularProfilesSlot(tSet)[[i]])) > 0 &&
        all(
          is.element(c("sampleid", "treatmentid"),
                     colnames(SummarizedExperiment::colData(molecularProfilesSlot(tSet)[[i]]))))) {
      tt <- table(SummarizedExperiment::colData(molecularProfilesSlot(tSet)[[i]])[ , "sampleid"], SummarizedExperiment::colData(molecularProfilesSlot(tSet)[[i]])[ , "treatmentid"])
      perturbation.info[rownames(tt), colnames(tt), names(molecularProfilesSlot(tSet))[i]] <- tt
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
#' checkTSetStructure(TGGATESsmall)
#'
#' @param tSet A \code{ToxicoSet} object
#' @param plotDist Should the function also plot the distribution of molecular data?
#' @param result.dir The path to the directory for saving the plots as a string, defaults to `tempdir()`
#'
#' @return Prints out messages whenever describing the errors found in the structure of the pset object passed in.
#'
#' @importFrom graphics hist
#' @importFrom grDevices dev.off pdf
#' @importFrom S4Vectors metadata
#' @importFrom CoreGx .message .warning .error
#'
#' @export
checkTSetStructure <- function(tSet, plotDist=FALSE, result.dir=".") {

    if(!file.exists(result.dir) && plotDist)
        dir.create(result.dir, showWarnings=FALSE, recursive=TRUE)

    for( i in seq_along(molecularProfilesSlot(tSet))) {
        profile <- molecularProfilesSlot(tSet)[[i]]
        if (is.null(names(metadata(profile))))
            .error(paste0("Please ensure all items in the metadata slot of
                 SummarizedExperiments are named. Item ", i, " of molecualrProfiles
                 does not have metadata names."))
        if (!("annotation" %in% names(metadata(profile))))
            .error(paste0("At minimum the SummarizedExperiments in molecularProfiles must contain
                    a metadata item names 'annotation' specifying the molecular datatype
                   the SummarizedExperiment contains! Item ", i, " of
                   molecularProfilesis missing annotation metadata."))
        nn <- names(molecularProfilesSlot(tSet))[i]

        if(plotDist) {
            if (S4Vectors::metadata(profile)$annotation == "rna" ||
                S4Vectors::metadata(profile)$annotation == "rnaseq")
            {
                pdf(file=file.path(result.dir, sprintf("%s.pdf", nn)))
                hist(SummarizedExperiment::assay(profile, 1), breaks = 100)
                dev.off()
            }
        }
        if (nrow(SummarizedExperiment::rowData(profile)) !=
            nrow(SummarizedExperiment::assay(profile, 1)))
        {
            .warning(sprintf("%s: number of features in fData is different from expression slots", nn))
        } else {
            .message(sprintf("%s: fData dimension is OK", nn))
        }
        if (nrow(SummarizedExperiment::colData(profile)) != ncol(SummarizedExperiment::assay(profile, 1)))
        {
            .warning(sprintf("%s: number of cell lines in pData is different from expression slots", nn))
        } else {
            .message(sprintf("%s: pData dimension is OK", nn))
        }

        if ("sampleid" %in% colnames(SummarizedExperiment::colData(profile))) {
            .message("sampleid OK!")
        } else {
            .warning(sprintf("%s: sampleid does not exist in pData columns", nn))
        }
        if ("batchid" %in% colnames(SummarizedExperiment::colData(profile))) {
            .message("batchid OK!")
        } else {
            .warning(sprintf("%s: batchid does not exist in pData columns", nn))
        }
        if (S4Vectors::metadata(profile)$annotation == "rna" ||
            S4Vectors::metadata(profile)$annotation == "rnaseq")
        {
            if ("BEST" %in% colnames(SummarizedExperiment::rowData(profile))) {
                .message("BEST is OK")
            } else {
                .warning(sprintf("%s: BEST does not exist in fData columns", nn))
            }

            if ("Symbol" %in% colnames(SummarizedExperiment::rowData(profile))) {
                .message("Symbol is OK")
            } else {
                .warning(sprintf("%s: Symbol does not exist in fData columns", nn))
            }
        }
        if ("sampleid" %in% colnames(SummarizedExperiment::colData(profile))) {
            if (!all(SummarizedExperiment::colData(profile)[, "sampleid"] %in% rownames(sampleInfo(tSet)))) {
                .warning(sprintf("%s: not all the cell lines in this profile are in cell lines slot", nn))
            }
        } else {
            .warning(sprintf("%s: sampleid does not exist in pData", nn))
        }
    }
    if ("tissueid" %in% colnames(sampleInfo(tSet))) {
        if ("unique.tissueid" %in% colnames(curation(tSet)$tissue)) {
        if (length(intersect(rownames(curation(tSet)$tissue), rownames(sampleInfo(tSet)))) != nrow(sampleInfo(tSet))) {
            .message("rownames of curation tissue slot should be the same as cell slot (curated cell ids)")
        } else {
            if(length(intersect(sampleInfo(tSet)$tissueid, curation(tSet)$tissue$unique.tissueid)) !=
                length(table(sampleInfo(tSet)$tissueid)))
            {
                .message("tissueid should be the same as unique tissue id from tissue curation slot")
            }
        }
        } else {
            .message("unique.tissueid which is curated tissue id across data set should be a column of tissue curation slot")
        }
        if(any(is.na(sampleInfo(tSet)[,"tissueid"]) | sampleInfo(tSet)[,"tissueid"] == "", na.rm = TRUE)) {
            .message(sprintf("There is no tissue type for this cell line(s): %s", paste(rownames(sampleInfo(tSet))[which(is.na(sampleInfo(tSet)[,"tissueid"]) | sampleInfo(tSet)[,"tissueid"] == "")], collapse = " ")))
        }
    } else {
        .warning("tissueid does not exist in cell slot")
    }

    if ("unique.sampleid" %in% colnames(curation(tSet)$cell)) {
        if(length(intersect(curation(tSet)$cell$unique.sampleid, rownames(sampleInfo(tSet)))) != nrow(sampleInfo(tSet))) {
            .message("rownames of cell slot should be curated cell ids")
        }
    } else {
        .message("unique.sampleid which is curated cell id across data set should be a column of cell curation slot")
    }

    if (length(intersect(rownames(curation(tSet)$cell), rownames(sampleInfo(tSet)))) != nrow(sampleInfo(tSet))) {
        .message("rownames of curation cell slot should be the same as cell slot (curated cell ids)")
    }

    if ("unique.treatmentid" %in% colnames(curation(tSet)$treatment)) {
        if(length(intersect(curation(tSet)$treatment$unique.treatmentid, treatmentNames(tSet))) != nrow(treatmentInfo(tSet))) {
            .message("rownames of drug slot should be curated drug ids")
        }
    } else {
        .message("unique.treatmentid which is curated drug id across data set should be a column of drug curation slot")
    }

    if (length(intersect(rownames(curation(tSet)$cell), rownames(sampleInfo(tSet)))) != nrow(sampleInfo(tSet))) {
        .message("rownames of curation drug slot should be the same as drug slot (curated drug ids)")
    }

    if (!is(sampleInfo(tSet), "data.frame")) {
        .warning("cell slot class type should be dataframe")
    }
    if (!is(treatmentInfo(tSet), "data.frame")) {
        .warning("drug slot class type should be dataframe")
    }
    if (datasetType(tSet) %in% c("sensitivity", "both"))
    {
        if(!is(sensitivityInfo(tSet), "data.frame")) {
            .warning("sensitivity info slot class type should be dataframe")
        }
        if("sampleid" %in% colnames(sensitivityInfo(tSet))) {
            if(!all(sensitivityInfo(tSet)[,"sampleid"] %in% rownames(sampleInfo(tSet)))) {
                .warning("not all the cell lines in sensitivity data are in cell slot")
            }
        } else {
            .warning("sampleid does not exist in sensitivity info")
        }
        if ("treatmentid" %in% colnames(sensitivityInfo(tSet))) {
            drug.ids <- unique(sensitivityInfo(tSet)[, "treatmentid"])
            drug.ids <- drug.ids[grep("///",drug.ids, invert=TRUE)]
            if (!all(drug.ids %in% treatmentNames(tSet))) {
                .message("not all the drugs in sensitivity data are in drug slot")
            }
        } else {
            .warning("treatmentid does not exist in sensitivity info")
        }

        if (any(!is.na(sensitivityRaw(tSet)))) {
            if(!all(dimnames(sensitivityRaw(tSet))[[1]] %in% rownames(sensitivityInfo(tSet)))) {
                .warning("For some experiments there is raw sensitivity data but no experimet information in sensitivity info")
            }
        }
        if (!all(rownames(sensitivityProfiles(tSet)) %in% rownames(sensitivityInfo(tSet)))) {
            .warning("For some experiments there is sensitivity profiles but no experimet information in sensitivity info")
        }
    }
}



# -------------------------------------------------------------------------
# Method Definitions ------------------------------------------------------
# -------------------------------------------------------------------------


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
#' @importMethodsFrom CoreGx show
#' @export
setMethod("show", signature=signature(object="ToxicoSet"), function(object) {
    callNextMethod(object)
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
setMethod("dim", signature("ToxicoSet"), function(x) {
    return(c(Cells=length(sampleNames(x)), Drugs=length(treatmentNames(x))))
})


#' @importFrom CoreGx updateSampleId
#' @aliases updateCellId
updateSampleId <- updateCellId <- function(object, new.ids=vector('character')) {
    CoreGx::updateSampleId(object, new.ids)
}

#' @importFrom CoreGx updateTreatmentId
#' @aliases updateDrugId
updateTreatmentId <- updateDrugId <- function(object, new.ids=vector('character')) {
    CoreGx:::updateTreamentId(object, new.ids)
}