#' Takes molecular data from a ToxicoSet, and summarises them
#' into one entry per drug and experimental condition.
#'
#' Given a ToxicoSet with molecular data, this function will summarize
#' the data into one profile per experimental condition (duration, dose level) using the chosen summary.stat
#' and return a SummarizedExperiment object, with one Assay corresponding to a requested drug.
#'
#' @examples
#' data(TGGATES_small)
#' summMP <- summarizeMolecularProfiles1(TGGATES_small, mDataType = "rna", cell.lines=cellNames(TGGATES_small),
#'                     drugs = head(drugNames(TGGATES_small)), features = fNames(TGGATES_small,"rna"),
#'                     duration = 8, dose = c("Control", "High"), summary.stat = 'median', fill.missing = TRUE, verbose=TRUE)
#' assays(summMP)[[drugNames(TGGATES_small)[1]]] #subset into expression matrix for a requested drug
#' colData(summMP) #summarization of phenoData for requested experiments
#' rowData(summMP) #featureData for requested experiments
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
#' @return \code{SummarizedExperiment} A SummarizedExperiment object with the molecular data summarized
#'   per cell line.
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom Biobase ExpressionSet exprs pData AnnotatedDataFrame assayDataElement assayDataElement<- fData<-
#' @import SummarizedExperiment
#' @import PharmacoGx
#' @export

summarizeToxicoMolecularProfiles1 <- function(tSet,
                                              mDataType,
                                              cell.lines,
                                              drugs,
                                              features,
                                              duration,
                                              dose,
                                              summary.stat = c("mean", "median", "first", "last"),
                                              fill.missing = TRUE, 
                                              summarize = TRUE, 
                                              verbose = TRUE
) {
  
  ############## ERRORTRAPPING NOT DONE ##############
  
  dd <- PharmacoGx::molecularProfiles(tSet, mDataType)[features, , drop = F] #expression matrix of the tSet
  pp <- PharmacoGx::phenoInfo(tSet, mDataType) #phenoData of the tSet
  ff <- PharmacoGx::featureInfo(tSet, mDataType)[features,,drop = F]
  
  unique.cells <- unique(cell.lines) #unique cell types (row names of the result)
  #subset phenoData to include only the experiments requested
  pp2 <- pp[(pp[,"cellid"] %in% unique.cells & pp[,"drugid"] %in% drugs 
             & pp[,"duration"] %in% duration & pp[,"dose_level"] %in% dose), , drop = F] #only the phenoData that is relevant to the request input
  dd2 <- dd[features,rownames(pp2), drop = F] #only the gene expression data that is relevant to the request input
  
  #vector of experimental conditions requested for each drug
  a <- paste(expand.grid(dose,duration)[,1], expand.grid(dose,duration)[,2], sep = ";")
  #b <- expand.grid(dose,duration)
  
  ddt <- dd2[,NA][,c(1:length(a)), drop = F]
  ppt <- pp2[FALSE,]
  
  exp.list <- list()
  cnt <- 0
  
  for (drug in drugs){
    cnt <- cnt + 1
    for (i in a){
      curr_dose <- sub(';.*$','', i)
      curr_dur <- sub('.*;','',i)
      
      pp3 <- pp2[(pp2[,"dose_level"] == curr_dose 
                  & pp2[,"duration"] == curr_dur
                  & pp2[,"drugid"] == drug), , drop = F]
      dd3 <- dd2[features,rownames(pp3), drop = F]
      
      if (ncol(dd3) > 1){
        switch(summary.stat, #ddr, ppr contains gene expression data, phenoData, for replicates
               "mean" = { ddr <- apply(dd3, 1, mean) },
               "median"={ ddr <- apply(dd3, 1, median) }, 
               "first"={ ddr <- dd3[ ,1 , drop=FALSE] },
               "last" = { ddr <- dd3[ , ncol(dd3), drop=FALSE] },
        )
        ppr <- apply(pp3[, , drop=FALSE], 2, function (x) {
          x <- paste(unique(as.character(x[!is.na(x)])), collapse="/")
          return (x)
        })
        ppr <- as.data.frame(t(ppr))
        ppr[!is.na(ppr) & ppr == ""] <- NA
        
        ddt <- cbind(ddt,ddr)
        ppt <- rbind(ppt,ppr)
      } else {
        ddt <- cbind(ddt,dd3)
        ppt <- rbind(ppt,pp3)
      }
    }
    ddt <- ddt[,-(1:length(a)), drop = F] #ddt contains the final expression matrix for a single drug
    colnames(ddt) <- a
    
    exp.list[[cnt]] <- ddt
  }
  names(exp.list) <- drugs
  ppf <- pp2[FALSE,]
  for (i in unique(ppt[,"dose_level"])){
    for (j in unique(ppt[,"duration"])){
      pp4 <- apply(ppt[ppt[,"dose_level"] == i & ppt[,"duration"] == j,,drop=F], 2, function (x) {
        x <- paste(unique(as.character(x[!is.na(x)])), collapse="///")
        return (x)
      })
      pp4 <- as.data.frame(t(pp4))
      pp4[!is.na(pp4) & pp4 == ""] <- NA
      ppf <- rbind(ppf,pp4)
    }
  }
  ppf <- as.data.frame(ppf,stringsAsFactors=FALSE)
  rownames(ppf) <- paste(ppf[,"dose_level"],";",ppf[,"duration"], sep = "")
  ppf <- ppf[match(rownames(ppf), colnames(exp.list[[1]])),]
  
  res <- SummarizedExperiment(assays = exp.list, rowData = ff, colData = ppf)
  
  return (res)
}
