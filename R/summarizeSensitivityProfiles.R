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
#' @param tSet \code{ToxicoSet} The ToxicoSet from which to extract the data
#' @param sensitivity.measure \code{character} which sensitivity sensitivity.measure to use? Use the
#'   sensitivityMeasures function to find out what measures are available for each TSet.
#' @param cell_lines \code{character} The cell lines to be summarized.
#'    If any cell lines has no data, it will be filled with missing values
#' @param drugs \code{character} The drugs to be summarized.
#'   If any drugs has no data, it will be filled with
#'   missing values. Defaults to include all drugs in the given tSet.
#' @param duration \code{numeric} The duration at which to summarize
#'   the drug-cell combo. This is a required parameter.
#' @param summary.stat \code{character} which summary method to use if there are repeated
#'   cell line-drug experiments? Choices are "mean", "median", "first", "last", "max", or "min"
#' @param fill.missing \code{boolean} should the missing cell lines not in the
#'   molecular data object be filled in with missing values?
#' @param verbose Should the function print progress messages?
#'
#' @return \code{matrix} A matrix with drugs going down the rows, cell lines across
#'   the columns, with the selected sensitivity statistic for each pair.
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats median
#' @importFrom reshape2 acast
#'
#' @export
summarizeSensitivityProfiles <- function(tSet,
                                         duration = NULL,
                                         cell_lines = NULL,
                                         drugs = NULL,
                                         sensitivity.measure="auc_recomputed",
                                         summary.stat = c("mean",
                                                          "median", "first",
                                                          "last", "max", "min"),
                                         fill.missing=TRUE, verbose=TRUE)
  {

  ## MISSING VALUE HANDLING FOR PARAMETERS
  # Get named list of defualt values for missing parameters
  argDefaultList <-
    paramMissingHandler(
      funName = "summarizeSensitivityProfiles", tSet = tSet,
      cell_lines = cell_lines, drugs = drugs, duration = duration
    )
  # Assign any missing parameter default values to function environment
  ## TODO:: I think we can do a for loop over index names?
  ## TODO:: Refactor to lapply
  if (length(argDefaultList) > 0) {
    for (idx in seq_along(argDefaultList)) {
      assign(names(argDefaultList)[idx], argDefaultList[[idx]])
    }
  }

  ## ERROR HANDLING FOR FUNCTION PARAMETERS
  paramErrorChecker(
    "summarizeSensitivtyProfiles", tSet = tSet, drugs = drugs,
    sensivity.measure = sensitivity.measure, duration = duration,
    summary.stat = summary.stat
    )

  summary.stat <- match.arg(summary.stat)

  pp <- ToxicoGx::sensitivityInfo(tSet)
  ## TODO:: Determine what this supposed to do?
  #ppRows <- which(pp$cellid %in% cell_lines & pp$drugid %in% drugs & pp$duration_h %in% duration) ### NEEDED to deal with duplicated rownames!!!!!!!

  if (sensitivity.measure != "max.conc") {
    #if the sensitivity.measure specified is not "max.conc"
    dd <- sensitivityProfiles(tSet)
  } else {
    #if the sensitivity.measure specified is "max.conc"
    if (!"max.conc" %in% colnames(ToxicoGx::sensitivityInfo(tSet))) {
      # if max.conc is not a column in sensitivityInfo:
      # call updateMaxConc, which finds the maximum dosage within sensitivity raw, puts
      # the value in a new column of sensitivity info called max.conc, and returns the tSet
      tSet <- updateMaxConc(tSet)

    }
    ##dd contains the sensitivity Info of the tSet
    dd <- ToxicoGx::sensitivityInfo(tSet)

  }

  #result is a matrix of NA's where # of rows, # columns is as specified:
  result <- matrix(NA_real_, nrow = length(drugs), ncol = length(cell_lines))
  #specify the row, column names of the result matrix
  rownames(result) <- drugs
  colnames(result) <- cell_lines


  ## TODO:: Finish progress bar
  # if(verbose){

  #   message(sprintf("Summarizing %s sensitivity data for:\t%s", sensitivity.measure, tSet@annotation$name))
  #   total <- length(drugs)*length(cell_lines)
  #   # create progress bar
  #   pb <- utils::txtProgressBar(min=0, max=total, style=3)
  #   i <- 1


  # }

  pp_dd <- cbind(pp[,c("cellid", "drugid","duration_h")], "sensitivity.measure" = dd[, sensitivity.measure])

  summary.function <- function(x) {
    if (all(is.na(x))) {
      return(NA_real_)
    }
    switch(summary.stat,
           "mean" = {
             return(mean(as.numeric(x), na.rm = TRUE))
           },
           "median" = {
             return(median(as.numeric(x), na.rm = TRUE))
           },
           "first" = {
             return(as.numeric(x)[[1]])
           },
           "last" = {
             return(as.numeric(x)[[length(x)]])
           },
           "max" = {
             return(max(as.numeric(x), na.rm = TRUE))
           },
           "min" = {
             return(min(as.numeric(x), na.rm = TRUE))
           })

  }

  pp_dd <- pp_dd[
    pp_dd[,"cellid"] %in% cell_lines &
    pp_dd[,"drugid"] %in% drugs &
    pp_dd[,"duration_h"] %in% duration,

    ]

  tt <- reshape2::acast(pp_dd, drugid~cellid, fun.aggregate = summary.function,
                        value.var = "sensitivity.measure")
  # tt <- tt[drugs, cell_lines]



  result[rownames(tt), colnames(tt)] <- tt

  if (!fill.missing) {

    myRows <- apply(result, 1, function(x) !all(is.na(x)))
    myCols <- apply(result, 2, function(x) !all(is.na(x)))
    result <- result[myRows, myCols]
  }
  return(result)
}
