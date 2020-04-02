#' tSet molecularProfiles from eSets to SEs
#'
#' Converts all ExpressionSet objects within the molecularProfiles slot of a
#'   ToxicoSet to SummarizedExperiments
#'
#' @param tSet \code{S4} A ToxicoSet containing molecular data in ExpressionSets
#'
#' @return \code{S4} A ToxicoSet containing molecular data in a SummarizedExperiments
#'
#' @importFrom parallel mclapply
#' @importFrom SummarizedExperiment assay assays assayNames
#' @importClassesFrom SummarizedExperiment SummarizedExperiment Assays
#' @importFrom Biobase exprs fData pData annotation protocolData
#' @importFrom S4Vectors SimpleList DataFrame
#' @importFrom stats setNames
#'
#' @export
#' @keywords internal
.convertTsetMolecularProfilesToSE <- function(tSet) {

  eSets <- tSet@molecularProfiles # Extract eSet data

  tSet@molecularProfiles <-
    lapply(eSets,
           function(eSet){

             # Build summarized experiment from eSet
             SE <- SummarizedExperiment::SummarizedExperiment(
               ## TODO:: Do we want to pass an environment for better memory efficiency?
               assays=SimpleList(as.list(Biobase::assayData(eSet))
               ),
               # Switch rearrange columns so that IDs are first, probes second
               rowData=S4Vectors::DataFrame(Biobase::fData(eSet),
                                            rownames=rownames(Biobase::fData(eSet))
               ),
               colData=S4Vectors::DataFrame(Biobase::pData(eSet),
                                            rownames=rownames(Biobase::pData(eSet))
               ),
               metadata=list("experimentData" = eSet@experimentData,
                             "annotation" = Biobase::annotation(eSet),
                             "protocolData" = Biobase::protocolData(eSet)
               )
             )
             ## TODO:: Determine if this can be done in the SE constructor?
             # Extract names from expression set
             SummarizedExperiment::assayNames(SE) <- Biobase::assayDataElementNames(eSet)
             # Assign SE to tSet
             mDataType <- Biobase::annotation(eSet)
             tSet@molecularProfiles[[mDataType]] <- SE
           })
  setNames(tSet@molecularProfiles, names(eSets))
  tSet
}

##' Validate tSet molecularProfiles Conversion
##'
##' Checks that all the information contained in an ExpressionSet molecularProfile
##'   was successfully tranferred to the SummarizedExperiment molecularProfile
##'
##' @param tSet_new \code{S4} a tSet containing molecularProfiles as SummarizedExperiments
##' @param tSet_old \code{S4} a tSet containing molecularProfiles as ExpressionSets
##'
##' @return \code{message} Any slots which are not the same
##'
##' @importFrom testthat expect_equal test_that
##' @import SummarizedExperiment
##' @import Biobase
##'
##' @export
##' @keywords internal
.validateTsetMolecularProfilesToSEConversion <- function(tSet_old, tSet_new) {

  # Testing that tSets are in correct order
  print("Checking is tSet structures are correct")

  testthat::expect_true(
    all(vapply(tSet_old@molecularProfiles, function(x) { is(x, "ExpressionSet") }, FUN.VALUE = logical(1))),
    info = "Old tSet doesn't contain ExpressionSet objects, maybe argument order is wrong?"
  )

  testthat::expect_true(
    all(vapply(tSet_new@molecularProfiles, function(x) { is(x, "SummarizedExperiment") }, FUN.VALUE = logical(1))),
    info = "New tSet doesn't contain SummarizedExperiment objects, maybe argument order is wrong?"
  )

  # Comparing molecularProfiles slot data
  print("Checking molecularProfiles slots hold equivalent data.")

  for (i in seq_len(length(tSet_old@molecularProfiles))) {
    for (j in seq_along(assays(tSet_new@molecularProfiles[[i]]))) {
      testthat::expect_true(
        all(
          as.list(assayData(tSet_old@molecularProfiles[[i]]))[[j]] ==
            assay(tSet_new@molecularProfiles[[i]], j),
          na.rm = TRUE
        ),
        info = "The assay data is not equivalent"
      )
    }
  }
  ## TODO:: Rewrite this as an apply statement
  for (i in seq_len(length(tSet_old@molecularProfiles))) { # Have to compare like this due to NAs in data
    # Checking phenoData
    testthat::expect_true(
      if (nrow(pData(tSet_old@molecularProfiles[[i]])) > 0) {
        all(
          as(tSet_old@molecularProfiles[[i]]@phenoData, "data.frame") ==
            as.data.frame(tSet_new@molecularProfiles[[i]]@colData[
              seq_len(length(tSet_new@molecularProfiles[[i]]@colData) -1)]),
          na.rm = TRUE)
      } else { TRUE },
      info = "The phenoData is not equivalent",
    )
    # Checking featureData
    testthat::expect_true(
      if (nrow(fData(tSet_old@molecularProfiles[[i]])) > 0) {
        all(
          as(tSet_old@molecularProfiles[[i]]@featureData, "data.frame") ==
            as.data.frame(tSet_new@molecularProfiles[[i]]@elementMetadata[
              seq_len(length(tSet_new@molecularProfiles[[i]]@elementMetadata) -1)]),
          na.rm = TRUE)
      } else { TRUE },
      info = "The featureData is not equivalent",
    )
    # Checking protocolData
    testthat::expect_true(
      all(
        as(tSet_old@molecularProfiles[[i]]@protocolData, "data.frame") ==
          as(tSet_new@molecularProfiles[[i]]@metadata$protocolData, "data.frame"),
        na.rm = TRUE),
      info = "The protocolData is not equivalent"
    )
  }

  testthat::expect_equal(
    lapply(tSet_old@molecularProfiles, function(x) { x@annotation }),
    lapply(tSet_new@molecularProfiles, function(x) { x@metadata$annotation }),
    info = "The annotation is not equivalent"
  )

  testthat::expect_equal(
    lapply(tSet_old@molecularProfiles, function(x) { x@experimentData }),
    lapply(tSet_new@molecularProfiles, function(x) { x@metadata$experimentData }),
    info = "The experimentData is not equivalent"
  )

  ##TODO:: Removed .__classVersion__ from SE as it is a property specific to eSet
  # testthat::expect_equal(
  #   lapply(tSet_old@molecularProfiles, function(x) { x@.__classVersion__ }),
  #   lapply(tSet_new@molecularProfiles, function(x) { x@metadata$.__classVersion__}),
  #   info = "The .__classVersion__ is not equivalent"
  # )

  # Comparing remainder of tSet slots; should not be affect by conversion
  print("Comparing remainder of tSet slots")

  testthat::test_that("Checking tSet@annotation slot is unchanged.", {
    testthat::expect_equal(tSet_old@annotation, tSet_new@annotation)
  })

  testthat::test_that("Checking tSet@cell slot is unchanged.", {
    testthat::expect_equal(tSet_old@cell, tSet_new@cell)
  })

  testthat::test_that("Checking tSet@drug slot is unchanged.", {
    testthat::expect_equal(tSet_old@drug, tSet_new@drug)
  })

  testthat::test_that("Checking tSet@sensitivity slot is unchanged.", {
    testthat::expect_equal(tSet_old@sensitivity, tSet_new@sensitivity)
  })

  testthat::test_that("Checking tSet@datasetType slot is unchanged.", {
    testthat::expect_equal(tSet_old@datasetType, tSet_new@datasetType)
  })

  testthat::test_that("Checking tSet@perturbation slot is unchanged.", {
    testthat::expect_equal(tSet_old@perturbation, tSet_new@perturbation)
  })

  testthat::test_that("Checking tSet@curation slot is unchanged.", {
    testthat::expect_equal(tSet_old@curation, tSet_new@curation)
  })
  message("Tests pass!")
}

##TODO:: Determine why CCLEsmall is 3x larger in memory after conversion?
#' Utility function to resave all datasets after modifying converttSetMolecularProfiles
#'
#' Converts all example dastasets specificed as an argument from
#'   molecularProfiles as ExpressionSet to molecularProfiles as
#'   SummarizedExperiment and saves them in the data folder
#'
#' @param datasets \code{character} A list of the example datasets to update
#'
#' @return \code{none} Works by side effects alone to resave all example
#'   datasets in a package to have SummarizedExperiments for molecularProfiles
#'
#' @export
#' @keywords internal
.resaveAllExampleDatasets <- function(datasets) {
  for (dataset in datasets) {
    dataDir <- paste0(grep('data', list.dirs(), value=TRUE))
    load(paste0(dataDir, '/', dataset, '_old.rda'))
    assign(dataset, .convertTsetMolecularProfilesToSE(get(dataset)))
    save(list=dataset, file=paste0(dataDir, '/', dataset, '.rda'), compress='xz')
  }
}
