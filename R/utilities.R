# tSet molecularProfiles from eSets to SEs
#
# Converts all ExpressionSet objects within the molecularProfiles slot of a
#   ToxicoSet to SummarizedExperiments
#
# @param tSet \code{S4} A ToxicoSet containing molecular data in ExpressionSets
#
# @return \code{S4} A ToxicoSet containing molecular data in a SummarizedExperiments
#
#' @importFrom SummarizedExperiment assay assays assayNames
#' @importClassesFrom SummarizedExperiment SummarizedExperiment Assays
#' @importFrom Biobase exprs fData pData annotation protocolData assayData experimentData
#' @importFrom S4Vectors SimpleList DataFrame
#' @importFrom stats setNames
#' @export
#' @keywords internal
.convertTsetMolecularProfilesToSE <- function(tSet) {

  eSets <- molecularProfilesSlot(tSet) # Extract eSet data

  molecularProfilesSlot(tSet) <-
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
             molecularProfilesSlot(tSet)[[mDataType]] <- SE
           })
  setNames(molecularProfilesSlot(tSet), names(eSets))
  tSet
}

# Validate tSet molecularProfiles Conversion
#
# Checks that all the information contained in an ExpressionSet molecularProfile
#   was successfully tranferred to the SummarizedExperiment molecularProfile
#
# @param tSet_new \code{S4} a tSet containing molecularProfiles as SummarizedExperiments
# @param tSet_old \code{S4} a tSet containing molecularProfiles as ExpressionSets
#
# @return \code{message} Any slots which are not the same
#
#' @importFrom assertthat are_equal
#' @importFrom SummarizedExperiment SummarizedExperiment Assays assay
#'   assayNames assayNames<-
#' @importFrom Biobase exprs fData pData annotation protocolData
#'   assayDataElementNames experimentData assayData
#' @keywords internal
.validateTsetMolecularProfilesToSEConversion <- function(tSet_old, tSet_new) {

  # Testing that tSets are in correct order
  message("Checking is tSet structures are correct")

  if(!all(vapply(tSet_old@molecularProfiles,
                 function(x) { is(x, "ExpressionSet") },
                 FUN.VALUE = logical(1)))
  ) message("Old tSet doesn't contain ExpressionSet objects, maybe argument
            order is wrong?")

  if(!all(vapply(molecularProfilesSlot(tSet_new),
                 function(x) { is(x, "SummarizedExperiment") },
                 FUN.VALUE = logical(1)))
  ) message("New tSet doesn't contain SummarizedExperiment objects, maybe
            argument order is wrong?")

  # Comparing molecularProfiles slot data
  message("Checking molecularProfiles slots hold equivalent data.")

  for (i in seq_len(length(tSet_old@molecularProfiles))) {
    for (j in seq_along(assays(molecularProfilesSlot(tSet_new)[[i]]))) {
      if(!all(as.list(assayData(tSet_old@molecularProfiles[[i]]))[[j]] ==
              as.list(assays(molecularProfilesSlot(tSet_new)[[i]]))[[j]],
          na.rm = TRUE)
        ) message("The assay data is not equivalent")
    }
  }
  ## TODO:: Rewrite this as an apply statement
  for (i in seq_len(length(tSet_old@molecularProfiles))) { # Have to compare like this due to NAs in data
    # Checking phenoData
    if(
      if (nrow(pData(tSet_old@molecularProfiles[[i]])) > 0) {
        !all(
          as(tSet_old@molecularProfiles[[i]]@phenoData, "data.frame") ==
            as.data.frame(molecularProfilesSlot(tSet_new)[[i]]@colData[
              seq_len(length(molecularProfilesSlot(tSet_new)[[i]]@colData) - 1)]),
          na.rm = TRUE)
      } else { FALSE }
    ) message("The phenoData is not equivalent")
    # Checking featureData
    if(
      if (nrow(fData(tSet_old@molecularProfiles[[i]])) > 0) {
        !all(
          as(tSet_old@molecularProfiles[[i]]@featureData, "data.frame") ==
            as.data.frame(molecularProfilesSlot(tSet_new)[[i]]@elementMetadata[
              seq_len(length(molecularProfilesSlot(tSet_new)[[i]]@elementMetadata) - 1)]),
          na.rm=TRUE)
      } else { FALSE }
    ) message("The featureData is not equivalent")
    # Checking protocolData
    if(
      !all(
        as(tSet_old@molecularProfiles[[i]]@protocolData, "data.frame") ==
          as(molecularProfilesSlot(tSet_new)[[i]]@metadata$protocolData, "data.frame"),
        na.rm = TRUE)
    ) message("The protocolData is not equivalent")
  }

  if(!assertthat::are_equal(
    lapply(tSet_old@molecularProfiles, function(x) { annotation(x) }),
    lapply(molecularProfilesSlot(tSet_new), function(x) { metadata(x)$annotation }))
  ) message("The annotation is not equivalent")

  if(!assertthat::are_equal(
    lapply(tSet_old@molecularProfiles, function(x) { experimentData(x) }),
    lapply(molecularProfilesSlot(tSet_new), function(x) { metadata(x)$experimentData }))
  ) message("The experimentData is not equivalent")

  # Comparing remainder of tSet slots; should not be affect by conversion
  message("Comparing remainder of tSet slots")

  if(!assertthat::are_equal(annotation(tSet_old), annotation(tSet_new)))
    message("The annotation slots are not equivalent!")

  if(!assertthat::are_equal(cellInfo(tSet_old), cellInfo(tSet_new)))
    message("The cell slots are not equivalent!")

  if(!assertthat::are_equal(drugInfo(tSet_old), drugInfo(tSet_new)))
    message("The drug slots are not equivalent!")

  if(!assertthat::are_equal(sensitivitySlot(tSet_old), sensitivitySlot(tSet_new)))
    message("The sensitivity slots are not equivalent!")

  if(!assertthat::are_equal(datasetType(tSet_old), datasetType(tSet_new)))
    message("The datasetType slots are not equivalent!")

  if(!assertthat::are_equal(tSet_old@perturbation, tSet_new@perturbation))
    message("The perturbation slots are not equivalent!")

 if(!assertthat::are_equal(curation(tSet_old), curation(tSet_new)))
   message("The curation slots are not equivalent!")
}

# Utility function to resave all datasets after modifying converttSetMolecularProfiles
#
# Converts all example dastasets specificed as an argument from
#   molecularProfiles as ExpressionSet to molecularProfiles as
#   SummarizedExperiment and saves them in the data folder
#
# @param datasets \code{character} A list of the example datasets to update
#
# @return \code{none} Works by side effects alone to resave all example
#   datasets in a package to have SummarizedExperiments for molecularProfiles
.resaveAllExampleDatasets <- function(datasets) {
  for (dataset in datasets) {
    dataDir <- paste0(grep('data', list.dirs(), value=TRUE))
    load(paste0(dataDir, '/', dataset, '_old.rda'))
    assign(dataset, .convertTsetMolecularProfilesToSE(get(dataset)))
    save(list=dataset, file=paste0(dataDir, '/', dataset, '.rda'), compress='xz')
  }
}


.eSetToSE <- function(eSet) {
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
    return(SE)
}

.validateESetToSEConversions <- function(eSet, SE) {
        for (j in seq_along(assays(SE))) {
            if(!all(as.list(assayData(eSet))[[j]] ==
                    as.list(assays(SE))[[j]],
                    na.rm = TRUE)
            ) message("The assay data is not equivalent")
        }
    ## TODO:: Rewrite this as an apply statement
        # Checking phenoData
        if(
            if (nrow(pData(eSet)) > 0) {
                !all(
                    as(eSet@phenoData, value="data.frame") ==
                    as.data.frame(SE@colData[
                        seq_len(length(SE@colData) - 1)]),
                    na.rm = TRUE)
            } else { FALSE }
        ) message("The phenoData is not equivalent")
        # Checking featureData
        if(
            if (nrow(fData(eSet)) > 0) {
                !all(
                    as(eSet@featureData, value="data.frame") ==
                    as.data.frame(SE@elementMetadata[
                        seq_len(length(SE@elementMetadata) - 1)]),
                    na.rm=TRUE)
            } else { FALSE }
        ) message("The featureData is not equivalent")
        # Checking protocolData
        if(
            !all(
                as(eSet@protocolData, value="data.frame") ==
                as(SE@metadata$protocolData, value="data.frame"),
                na.rm = TRUE)
        ) message("The protocolData is not equivalent")

    if(!assertthat::are_equal(
        annotation(eSet),
        metadata(SE)$annotation
    )) message("The annotation is not equivalent")

    if(!assertthat::are_equal(
        experimentData(eSet),
        metadata(SE)$experimentData
    )) message("The experimentData is not equivalent")
}