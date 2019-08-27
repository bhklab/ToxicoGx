#' Takes the sensitivity data from a RadioSet, and summarises them into a
#' drug vs cell line table
#'
#' This function creates a table with cell lines as rows and radiation types as columns,
#' summarising the drug senstitivity data of a RadioSet into drug-cell line
#' pairs
#'
#' @examples
#' data(Cleveland_small)
#' GDSCauc <- summarizeSensitivityProfiles(Cleveland_small, sensitivity.measure='AUC_published')
#'
#' @param rSet [RadioSet] The RadioSet from which to extract the data
#' @param sensitivity.measure [character] which sensitivity sensitivity.measure to use? Use the
#'   sensitivityMeasures function to find out what measures are available for each PSet.
#' @param cell.lines \code{character} The cell lines to be summarized.
#'    If any cell lines has no data, it will be filled with
#'   missing values
#' @param radiation.types \code{character} The radiation types to be summarized.
#'   If any radiation type has no data, it will be filled with
#'   missing values
#' @param summary.stat \code{character} which summary method to use if there are repeated
#'   cell line-drug experiments? Choices are "mean", "median", "first", or "last"
#' @param fill.missing \code{boolean} should the missing cell lines not in the
#'   molecular data object be filled in with missing values?
#' @param verbose Should the function print progress messages?
#' @return [matrix] A matrix with cell lines going down the rows, radiation types across
#'   the columns, with the selected sensitivity statistic for each pair.
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats median
#' @importFrom reshape2 acast
#' @export

summarizeSensitivityProfiles <- function(rSet, sensitivity.measure="AUC_recomputed", cell.lines, radiation.types, summary.stat=c("mean", "median", "first", "last", "max", "min"), fill.missing=TRUE, verbose=TRUE){
	summary.stat <- match.arg(summary.stat)
  #sensitivity.measure <- match.arg(sensitivity.measure)
  if (!(sensitivity.measure %in% c(colnames(sensitivityProfiles(rSet)),"max.conc"))) {
    stop (sprintf("Invalid sensitivity measure for %s, choose among: %s", rSet@annotation$name, paste(colnames(sensitivityProfiles(rSet)), collapse=", ")))
  }
  if (missing(cell.lines)) {
    cell.lines <- cellNames(rSet)
  }
  if (missing(radiation.types)) {
    if (sensitivity.measure != "Synergy_score")
    {
      radTypes <- radiationTypes(rSet)
    }else{
      radTypes <- sensitivityInfo(rSet)[grep("///", sensitivityInfo(rSet)$radiation.type), "radiation.type"]
    }
  }

  pp <- sensitivityInfo(rSet)
  ppRows <- which(pp$cellid %in% cell.lines & pp$radiation.type %in% radTypes) ### NEEDED to deal with duplicated rownames!!!!!!!
  if(sensitivity.measure != "max.conc") {
    dd <- sensitivityProfiles(rSet)
  } else {

    if(!"max.conc"%in% colnames(sensitivityInfo(rSet))){

      rSet <- updateMaxConc(rSet)

    }
    dd <- sensitivityInfo(rSet)

  }

  result <- matrix(NA_real_, nrow=length(radTypes), ncol=length(cell.lines))
  rownames(result) <- radTypes
  colnames(result) <- cell.lines

  # if(verbose){

  #   message(sprintf("Summarizing %s sensitivity data for:\t%s", sensitivity.measure, rSet@annotation$name))
  #   total <- length(radTypes)*length(cell.lines)
  #   # create progress bar
  #   pb <- utils::txtProgressBar(min=0, max=total, style=3)
  #   i <- 1


  # }

  pp_dd <- cbind(pp[,c("cellid", "radiation.type")], "sensitivity.measure"=dd[, sensitivity.measure])


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

  pp_dd <- pp_dd[pp_dd[,"cellid"]%in%cell.lines & pp_dd[,"radiation.type"]%in%radTypes,]

  tt <- reshape2::acast(pp_dd, radiation.type~cellid, fun.aggregate=summary.function, value.var="sensitivity.measure")
 # tt <- tt[radTypes, cell.lines]



  result[rownames(tt), colnames(tt)] <- tt

	if (!fill.missing) {

    myRows <- apply(result, 1, function(x) !all(is.na(x)))
    myCols <- apply(result, 2, function(x) !all(is.na(x)))
    result <- result[myRows, myCols]
	}
  return(result)
}
