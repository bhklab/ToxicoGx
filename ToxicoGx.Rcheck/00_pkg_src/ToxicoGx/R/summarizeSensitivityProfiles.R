#' Takes the sensitivity data from a ToxicoSet, and summarises them into a
#' drug vs cell line table
#'
#' This function creates a table with drug as rows and cell lines as columns,
#' summarising the drug senstitivity data of a ToxicoSet into drug-cell line
#' pairs for a specified experiment duration.
#'
#' @examples
#' data(TGGATESsmall)
#' TGGATESauc <- summarizeSensitivityProfiles(TGGATESsmall, sensitivity.measure='auc_recomputed')
#'
#' @param tSet [ToxicoSet] The ToxicoSet from which to extract the data
#' @param sensitivity.measure [character] which sensitivity sensitivity.measure to use? Use the
#'   sensitivityMeasures function to find out what measures are available for each TSet.
#' @param cell.lines \code{character} The cell lines to be summarized.
#'    If any cell lines has no data, it will be filled with
#'   missing values
#' @param drugs \code{character} The drugs to be summarized.
#'   If any drugs has no data, it will be filled with
#'   missing values
#' @param duration \code{numeric} The duration at which to summarize
#'   the drug-cell combo.
#' @param summary.stat \code{character} which summary method to use if there are repeated
#'   cell line-drug experiments? Choices are "mean", "median", "first", "last", "max", or "min"
#' @param fill.missing \code{boolean} should the missing cell lines not in the
#'   molecular data object be filled in with missing values?
#' @param verbose Should the function print progress messages?
#'
#' @return [matrix] A matrix with drugs going down the rows, cell lines across
#'   the columns, with the selected sensitivity statistic for each pair.
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats median
#' @importFrom reshape2 acast
#'
#' @export
#'
summarizeSensitivityProfiles <- function(tSet,
                                         sensitivity.measure="auc_recomputed",
                                         cell.lines,
                                         drugs,
                                         duration,
                                         summary.stat=c("mean", "median", "first", "last", "max", "min"),
                                         fill.missing=TRUE, verbose=TRUE) {

  summary.stat <- match.arg(summary.stat)
  # sensitivity.measure <- match.arg(sensitivity.measure)
  if (!(sensitivity.measure %in% c(colnames(sensitivityProfiles(tSet)),"max.conc"))) {
    #if the sensitivity.measure specified is not available for the tSet
    stop (sprintf("Invalid sensitivity measure for %s, choose among: %s", tSet@annotation$name, paste(colnames(sensitivityProfiles(tSet)), collapse=", ")))
  }

  if (missing(cell.lines)) {
    #if cell.lines was not specified
    cell.lines <- cellNames(tSet)
  }
  if (missing(drugs)) {
    #if drugs was not specified
    if (sensitivity.measure != "Synergy_score")
    {
      #if the sensitivity.measure specified was not "Synergy_score"
      drugs <- drugNames(tSet)
    }else{
      #wtf is this
      drugs <- sensitivityInfo(tSet)[grep("///", sensitivityInfo(tSet)$drugid), "drugid"]
    }
  }
  if (missing(duration)) { # Selects the first row's duration if no duration is specified in argument
    duration <- sensitivityInfo(tSet)$duration_h[1]
  }
  if (length(duration) > 1 ) {
    stop("Please enter only one duration value to be summarized.")
  }

  pp <- sensitivityInfo(tSet)
  ppRows <- which(pp$cellid %in% cell.lines & pp$drugid %in% drugs & pp$duration_h %in% duration) ### NEEDED to deal with duplicated rownames!!!!!!!
  #ppRows <- which()
  if(sensitivity.measure != "max.conc") {
    #if the sensitivity.measure specified is not "max.conc"
    dd <- sensitivityProfiles(tSet)
  } else {
    #if the sensitivity.measure specified is "max.conc"
    if(!"max.conc"%in% colnames(sensitivityInfo(tSet))){
      # if max.conc is not a column in sensitivityInfo:
      # call updateMaxConc, which finds the maximum dosage within sensitivity raw, puts
      # the value in a new column of sensitivity info called max.conc, and returns the tSet
      tSet <- updateMaxConc(tSet)

    }
    ##dd contains the sensitivity Info of the tSet
    dd <- sensitivityInfo(tSet)

  }

  #result is a matrix of NA's where # of rows, # columns is as specified:
  result <- matrix(NA_real_, nrow=length(drugs), ncol=length(cell.lines))
  #specify the row, column names of the result matrix
  rownames(result) <- drugs
  colnames(result) <- cell.lines

  # if(verbose){

  #   message(sprintf("Summarizing %s sensitivity data for:\t%s", sensitivity.measure, tSet@annotation$name))
  #   total <- length(drugs)*length(cell.lines)
  #   # create progress bar
  #   pb <- utils::txtProgressBar(min=0, max=total, style=3)
  #   i <- 1


  # }

  pp_dd <- cbind(pp[,c("cellid", "drugid","duration_h")], "sensitivity.measure"=dd[, sensitivity.measure])


  summary.function <- function(x) {
    if(all(is.na(x))){
      return(NA_real_)
    }
    switch(summary.stat,
           "mean" = {
             return(mean(as.numeric(x), na.rm=TRUE))
           },
           "median" = {
             return(median(as.numeric(x), na.rm=TRUE))
           },
           "first" = {
             return(as.numeric(x)[[1]])
           },
           "last" = {
             return(as.numeric(x)[[length(x)]])
           },
           "max"= {
             return(max(as.numeric(x), na.rm=TRUE))
           },
           "min" = {
             return(min(as.numeric(x), na.rm=TRUE))
           })

  }

  pp_dd <- pp_dd[pp_dd[,"cellid"]%in%cell.lines & pp_dd[,"drugid"]%in%drugs & pp_dd[,"duration_h"]%in%duration,]

  tt <- reshape2::acast(pp_dd, drugid~cellid, fun.aggregate=summary.function, value.var="sensitivity.measure")
  # tt <- tt[drugs, cell.lines]



  result[rownames(tt), colnames(tt)] <- tt

  if (!fill.missing) {

    myRows <- apply(result, 1, function(x) !all(is.na(x)))
    myCols <- apply(result, 2, function(x) !all(is.na(x)))
    result <- result[myRows, myCols]
  }
  return(result)
}
