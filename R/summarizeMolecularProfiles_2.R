#' Takes molecular data from a ToxicoSet, and summarises them
#' into one entry per drug for the requested duration
#'
#' Given a ToxicoSet with molecular data, this function will summarize
#' the data into one profile per cell line, using the chosen summary.stat.
#'
#' @examples
#' data(TGGATES_small)
#' TGGATES_small <- summarizeMolecularProfiles2(TGGATES_small,
#'                     mDataType = "rna", cell.lines=cellNames(TGGATES_small),
#'                     summary.stat = 'median', fill.missing = TRUE, verbose=TRUE)
#' TGGATES_small
#'
#' @param tSet \code{ToxicoSet} The ToxicoSet to summarize
#' @param mDataType \code{character} which one of the molecular data types
#' to use in the analysis, out of all the molecular data types available for the tSet
#' for example: rna
#' @param cell.lines \code{character} The cell lines to be summarized.
#'   If any cell.line has no data, missing values will be created
#' @param drugs \code{character} The drugs to be summarized
#' @param features \code{character} A vector of the feature names to include in the summary
#' @param duration \code{numeric} A vector of durations to summarize across
#' @param dose \code{character} The dose level to summarize replicates across
#' @param summary.stat \code{character} which summary method to use if there are repeated
#'   cell.lines? Choices are "mean", "median", "first", or "last"
#' @param fill.missing \code{boolean} should the missing cell lines not in the
#'   molecular data object be filled in with missing values?
#' @param summarize A flag which when set to FALSE (defaults to TRUE) disables summarizing and
#'   returns the data unchanged as a ExpressionSet
#' @param verbose \code{boolean} should messages be printed
#' @return \code{matrix} An ExpressionSet with the molecular data summarized
#'   per cell line.
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom Biobase ExpressionSet exprs pData AnnotatedDataFrame assayDataElement assayDataElement<- fData<-
#' @export

summarizeToxicoMolecularProfiles2 <- function(tSet, 
                                              mDataType, 
                                              cell.lines, # eg. "Hepatocyte"
                                              drugs, # eg. "Valproic.acid"
                                              features,  # eg. "ENSG00000000003"
                                              duration = c(2,8,24), # eg. c(2,8,24)
                                              dose = c("Control","Low","Middle","High"), # eg. "Control"
                                              summary.stat = c("mean", "median", "first", "last"), 
                                              fill.missing = TRUE, 
                                              summarize = TRUE, 
                                              verbose = TRUE
) {
  
  ############## ERRORTRAPPING NOT DONE ##############
  
  dd <- molecularProfiles(tSet, mDataType)[features, , drop = F] #expression matrix of the tSet
  pp <- phenoInfo(tSet, mDataType) #phenoData of the tSet
  
  unique.cells <- unique(cell.lines) #unique cell types (row names of the result)
  pp2 <- pp[(pp[,"cellid"] %in% unique.cells & pp[,"drugid"] %in% drugs 
             & pp[,"duration"] %in% duration & pp[,"dose_level"] %in% dose), , drop = F] #only the phenoData that is relevant to the request input
  dd2 <- dd[features,rownames(pp2), drop = F] #only the gene expression data that is relevant to the request input
  
  ddt <- dd2[,NA][,1, drop = F]
  ppt <- pp2[FALSE,]
  
  for (i in unique.cells){
    pp3 <- pp2[pp2[,"cellid"] == i, , drop=F] #only the phenoData for the specific cell type
    dd3 <- dd2[features, rownames(pp3), drop = F] #only the gene expression data for the specific cell type
    
    #if there are replicates
    if (ncol(dd3) > 1){
      switch(summary.stat, #ddr, ppr contains gene expression data, phenoData, for replicates
             "mean" = {
               ddr <- apply(dd3, 1, mean)
             },
             "median"={
               ddr <- apply(dd3, 1, median)
             }, 
             "first"={
               ddr <- dd3[ ,1 , drop=FALSE]
             },
             "last" = {
               ddr <- dd3[ , ncol(dd3), drop=FALSE]
             },
      )
      ppr <- apply(pp3[, , drop=FALSE], 2, function (x) {
        x <- paste(unique(as.character(x[!is.na(x)])), collapse="///")
        return (x)
      })
      ppr <- as.data.frame(t(ppr))
      ppr[!is.na(ppr) & ppr == ""] <- NA
      
      # ddt, ppt contains final result
      ddt <- cbind(ddt,ddr)
      ppt <- rbind(ppt,ppr)
    } else { #if no replicates
      ddt <- cbind(ddt,dd3)
      ppt <- cbind(ppt,pp3)
    }
  }
  ddt <- ddt[,-1, drop = F]
  
  colnames(ddt) <- rownames(ppt) <- unique.cells
  ddt <- ddt[, unique.cells, drop = F]
  ppt <- ppt[unique.cells, , drop = F]
  ppt[,"cellid"] <- unique.cells
  
  res <- ExpressionSet(ddt)
  ppt <- as.data.frame(ppt, stringsAsFactors=FALSE)
  ppt$tissueid <- cellInfo(tSet)[ppt$cellid, "tissueid"]
  
  Biobase::pData(res) <- ppt
  Biobase::fData(res) <- featureInfo(tSet, mDataType)
  
  res <- res[features,]
  Biobase::protocolData(res) <- Biobase::AnnotatedDataFrame()
  if(!is.null(assayDataElement(res, "se.exprs"))) assayDataElement(res,"se.exprs") <- NULL
  return(res)
}
