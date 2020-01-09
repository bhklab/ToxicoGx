#' Return a table of ToxicoSets available for download
#'
#' The function fetches a table of all ToxicoSets available for download from
#' the PharmacoGx server. The table includes the names of the PharamcoSet, the
#' types of data available in the object, and the date of last update.
#'
#' @examples
#' if (interactive()){
#' availabletSets()
#' }
#'
#' @param saveDir \code{character} Directory to save the table of tSets
#' @param myfn \code{character} The filename for the table of tSets
#' @param verbose \code{bool} Should status messages be printed during download.
#' @return A data.frame with details about the available ToxicoSet objects
#' @export
#' @import downloader
#' @importFrom utils read.table write.table
availableTSets <- function(saveDir=tempdir(), myfn="availableToxicoSets.csv", verbose=TRUE){

  if (!file.exists(saveDir)) {
    dir.create(saveDir, recursive = TRUE)
  }

  downloader::download("https://ndownloader.figshare.com/files/19347956?private_link=d286d7386d5f5e778585",
                       destfile = file.path(saveDir, myfn),
                       quiet = !verbose)

  tSetTable <- read.csv(file.path(saveDir, myfn), header = TRUE, stringsAsFactors = FALSE)
  return(tSetTable)
}

#' Download a ToxicoSet object
#'
#' This function allows you to download a \code{ToxicoSet} object for use with this
#' package. The \code{ToxicoSets} have been extensively curated and organised within
#' a PharacoSet class, enabling use with all the analysis tools provided in
#' \code{PharmacoGx}.
#'
#' @examples
#' if (interactive()){
#' downloadtSet("TGGATESvignette")
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
#' @export
#' @import downloader
downloadTSet <- function(name, saveDir = tempdir(), tSetFileName = NULL, verbose = TRUE) {

  if (missing(saveDir)) {message("Downloading to temporary folder... Use saveDir parameter to save to a specific path")}
  tSetTable <- availableTSets(saveDir = saveDir)

  whichx <- match(name, tSetTable[, 1])
  if (is.na(whichx)) {
    stop('Unknown Dataset. Please use the availabletSets() function for the table of available PharamcoSets.')
  }

  if (!file.exists(saveDir)) {
    dir.create(saveDir, recursive = TRUE)
  }

  if (is.null(tSetFileName)) {
    tSetFileName <- paste0(tSetTable[whichx,"ToxicoSet.Name"], ".rda")
  }
  if (!file.exists(file.path(saveDir, tSetFileName))) {
    downloader::download(url = as.character(tSetTable[whichx,"URL"]), destfile = file.path(saveDir, tSetFileName), quiet = !verbose)
  }

  load(file.path(saveDir, tSetFileName), envir = globalenv())
}

#' @importFrom utils read.table write.table
.createtSetEntry <- function(tSet, outfn) {

  if (file.exists(outfn)) {
    tSetTable <- read.table(outfn, as.is = TRUE)
    newrow <- c(tSetName(tSet), tSet@datasetType, paste(names(tSet@molecularProfiles), collapse = "/"), tSet@annotation$dateCreated, NA)
    tSetTable <- rbind(tSetTable, newrow)
    rownames(tSetTable) <- tSetTable[, 1]
    write.table(tSetTable, file = outfn)
  } else {
    newrow <- c(tSetName(tSet), tSet@datasetType, paste(names(tSet@molecularProfiles), collapse = "/"), tSet@annotation$dateCreated, NA)
    tSetTable <- t(matrix(newrow))
    colnames(tSetTable) <- c("ToxicoSet.Name","Description", "Available.Molecular.Profiles","Date.Updated","URL")
    rownames(tSetTable) <- tSetTable[,1]
    write.table(tSetTable, file = outfn)
  }
}
