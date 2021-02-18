#' Return a table of ToxicoSets available for download
#'
#' The function fetches a table of all ToxicoSets available for download from
#' the ToxicoGx server. The table includes the names of the ToxicoSet, the
#' types of data available in the object, and the date of last update.
#' 
#' Much more information on the processing of the data and data provenance can be found at:
#' www.orcestra.ca 
#'
#' @examples
#' if (interactive()){
#' availablePSets()
#' }
#'
#' @param canonical [`logical`] Should available TSets show only official TSets, or should
#'   user generated TSets be included?
#' 
#' @return A data.frame with details about the available ToxicoSet objects
#' @export
#' @import jsonlite
availableTSets <- function(canonical=TRUE){
  require(jsonlite)
  if (canonical) {
    avail.tsets <- fromJSON("https://www.orcestra.ca/api/toxicosets/canonical")
  } else {
    return("Only canonical TSets are available at the moment")
  }
  
  tSetTable <- data.frame("ToxicoSet.Name" = avail.tsets$dataset$name,
                          "Data.Source" = avail.tsets$dataset$versionInfo$data$rawMicroarrayData,
                          "Date.Created" = avail.tsets$dateCreated,
                          "URL" = avail.tsets$downloadLink, stringsAsFactors = FALSE, check.names = FALSE)
  return(tSetTable)
}

#' Download a ToxicoSet object
#'
#' This function allows you to download a \code{ToxicoSet} object for use with this
#' package. The \code{ToxicoSets} have been extensively curated and organised within
#' a ToxicoSet class, enabling use with all the analysis tools provided in
#' \code{ToxicoGx}.
#'
#' @examples
#' if (interactive()) {
#' drugMatrix_rat <- downloadtSet("drugMatrix_rat")
#' }
#'
#' @param name \code{Character} string, the name of the PhamracoSet to download.
#' @param saveDir \code{Character} string with the folder path where the
#'     ToxicoSet should be saved. Defaults to \code{'./tSets/'}. Will create
#'     directory if it does not exist.
#' @param tSetFileName \code{character} string, the file name to save the dataset under
#' @param verbose \code{bool} Should status messages be printed during download.
#'   Defaults to TRUE.
#' @return A tSet object with the dataset, downloaded from our server
#' @import downloader
#' @export
downloadTSet <- function(name, saveDir = tempdir(), tSetFileName = NULL, verbose = TRUE) {

  if (missing(saveDir)) {message("Downloading to temporary folder... Use saveDir parameter to save to a specific path")}
  tSetTable <- availableTSets(canonical=T)

  whichx <- match(name, tSetTable[, 1])
  if (is.na(whichx)) {
    stop('Unknown Dataset. Please use the availabletSets() function for the
         table of available ToxicoSets.')
  }

  if (!file.exists(saveDir)) {
    dir.create(saveDir, recursive = TRUE)
  }

  if (is.null(tSetFileName)) {
    tSetFileName <- paste0(tSetTable[whichx,"ToxicoSet.Name"], ".rds")
  }
  if (!file.exists(file.path(saveDir, tSetFileName))) {
    downloader::download(url = as.character(tSetTable[whichx,"URL"]),
                         destfile = file.path(saveDir, tSetFileName),
                         quiet = !verbose, mode='wb')
  }

  print(file.path(saveDir, tSetFileName))
  tSet <- readRDS(file.path(saveDir, tSetFileName))

  return(tSet)
}

#' @importFrom utils read.table write.table
.createtSetEntry <- function(tSet, outfn) {

  if (file.exists(outfn)) {
    tSetTable <- read.table(outfn, as.is = TRUE)
    newrow <- c(name(tSet), datasetType(tSet), paste(names(molecularProfilesSlot(tSet)), collapse = "/"), annotation(tSet)$dateCreated, NA)
    tSetTable <- rbind(tSetTable, newrow)
    rownames(tSetTable) <- tSetTable[, 1]
    write.table(tSetTable, file = outfn)
  } else {
    newrow <- c(name(tSet), datasetType(tSet), paste(names(molecularProfilesSlot(tSet)), collapse = "/"), annotation(tSet)$dateCreated, NA)
    tSetTable <- t(matrix(newrow))
    colnames(tSetTable) <- c("ToxicoSet.Name","Data.Source","Date.Updated","URL")
    rownames(tSetTable) <- tSetTable[,1]
    write.table(tSetTable, file = outfn)
  }
}
