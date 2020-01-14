#' Takes molecular data from a ToxicoSet, and summarises them
#' into one entry per drug and experimental condition.
#'
#' Given a ToxicoSet with molecular data, this function will summarize
#' the data into one profile per experimental condition (duration, dose level)
#' using the chosen summary.stat and return a SummarizedExperiment object, with
#' one Assay corresponding to a requested drug.
#'
#' @examples
#' data(TGGATESsmall)
#' summMP <- ToxicoGx::summarizeMolecularProfiles(
#'   tSet = TGGATESsmall, mDataType = "rna",
#'   cell_lines=cellNames(TGGATESsmall), drugs = head(drugNames(TGGATESsmall)),
#'   features = fNames(TGGATESsmall,"rna")[1:100], duration = "8",
#'   dose = c("Control", "High"), summary.stat = "median",
#'   fill.missing = TRUE, verbose=TRUE
#'   )
#'
#' #subset into expression matrix for a requested drug
#' assays <- SummarizedExperiment::assays(summMP)[[drugNames(TGGATESsmall)[1]]]
#' #summarization of phenoData for requested experiments
#' phenoData <- SummarizedExperiment::colData(summMP)
#' #summarization of phenoData for requested experiments
#' featureData <- SummarizedExperiment::rowData(summMP) #featureData for requested experiments
#'
#' @param tSet \code{ToxicoSet} The ToxicoSet to summarize
#' @param mDataType \code{character} which one of the molecular data types
#' to use in the analysis, out of all the molecular data types available for the tSet
#' for example: rna
#' @param cell_lines \code{character} The cell lines to be summarized.
#'   If any cell.line has no data, missing values will be created
#' @param drugs \code{character} The drugs to be summarized
#' @param features \code{character} A vector of the feature names to include in the summary
#' @param duration \code{character} A vector of durations to summarize across
#' @param dose \code{character} The dose level to summarize replicates across
#' @param summary.stat \code{character} which summary method to use if there are repeated
#'   cell_lines? Choices are "mean", "median", "first", or "last"
#' @param fill.missing \code{boolean} should the missing cell lines not in the
#'   molecular data object be filled in with missing values?
#' @param summarize A flag which when set to FALSE (defaults to TRUE) disables summarizing and
#'   returns the data unchanged as a ExpressionSet
#' @param verbose \code{boolean} should messages be printed
#' @return \code{SummarizedExperiment} A SummarizedExperiment object with the molecular data summarized
#'   per cell line.
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom Biobase ExpressionSet exprs pData AnnotatedDataFrame assayDataElement assayDataElement<- fData<-
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @export
#'
# Output is of SummarizedExperiment class
## TODO:: Rewrite this using apply functions instead of for loops
summarizeMolecularProfiles <-
  function(tSet,
           mDataType,
           cell_lines = NULL, # Defaults get set in paramMissingHandler call
           drugs = NULL,
           features = NULL,
           duration = NULL,
           dose = c("Control", "Low", "Middle", "High"),
           summary.stat = c("mean", "median", "first", "last"),
           fill.missing = TRUE,
           summarize = TRUE,
           verbose = TRUE
  ) {

    ##### CHECKING INPUT VALIDITY #####

    ## MISSING VALUE HANDLING FOR PARAMETERS
    # Get named list of defualt values for missing parameters
    argDefaultList <-
      paramMissingHandler(
        funName = "summarizeMolecularProfiles", tSet = tSet,
        mDataType = mDataType, cell_lines = cell_lines, drugs = drugs,
        features = features, duration = duration
      )

    # Assign any missing parameter default values to function environment
    if (length(argDefaultList) > 0) {
      for (idx in seq_along(argDefaultList)) {
        assign(names(argDefaultList)[idx], argDefaultList[[idx]])
      }
    }

    ## TODO:: Standardized parameter names across all function
    ## ERROR HANDLING FOR PARAMETERS
    paramErrorChecker(
      "summarizeMolecularProfiles", tSet = tSet,
      mDataType = mDataType, cell_lines = cell_lines, drugs = drugs,
      features = features, duration = duration, dose = dose,
      summary.stat = summary.stat
    )

    ##### FUNCTION LOGIC BEGINS #####

    dd <- ToxicoGx::molecularProfiles(tSet, mDataType)[features, , drop = F] #expression matrix of the tSet
    pp <- ToxicoGx::phenoInfo(tSet, mDataType) #phenoData of the tSet
    ff <- ToxicoGx::featureInfo(tSet, mDataType)[features,,drop = F]

    unique.cells <- unique(cell_lines) #unique cell types (row names of the result)
    #subset phenoData to include only the experiments requested
    pp2 <- pp[(pp[,"cellid"] %in% unique.cells & pp[,"drugid"] %in% drugs
               & pp[,"duration"] %in% duration & pp[,"dose_level"] %in% dose), , drop = F] #only the phenoData that is relevant to the request input
    dd2 <- dd[features,rownames(pp2), drop = F] #only the gene expression data that is relevant to the request input

    #vector of experimental conditions requested for each drug
    a <- paste(expand.grid(dose,duration)[,1], expand.grid(dose, duration)[,2], sep = ";")
    #b <- expand.grid(dose,duration)

    ddt <- dd[,NA][,c(1:length(a)), drop = F]
    ppt <- pp[FALSE,]

    exp.list <- list()
    cnt <- 0
    blank <- ddt[,1,drop = F]

    #lapply(drugs, function(drug){
    for (drug in drugs) {
      cnt <- cnt + 1
      #lapply(a, function(i){
      for (i in a) {
        if (verbose == TRUE) {
          message(i)
        }
        ## TODO:: Is the print error occuring here?
        curr_dose <- sub(';.*$','', i)
        curr_dur <- sub('.*;','', i)

        pp3 <- pp2[(pp2[,"dose_level"] == curr_dose
                    & pp2[,"duration"] == curr_dur
                    & pp2[,"drugid"] == drug), , drop = F]
        dd3 <- dd2[features,rownames(pp3), drop = F]

        if (ncol(dd3) > 1){ #if there are replicates
          switch(summary.stat, #ddr, ppr contains gene expression data, phenoData, for replicates
                 "mean" = { ddr <- apply(dd3, 1, mean) },
                 "median"= { ddr <- apply(dd3, 1, median) },
                 "first"= { ddr <- dd3[ ,1 , drop=FALSE] },
                 "last" = { ddr <- dd3[ , ncol(dd3), drop=FALSE] },
          )
          ppr <- apply(pp3[, , drop=FALSE], 2, function (x) {
            x <- paste(unique(as.character(x[!is.na(x)])), collapse="/")
            return(x)
          })
          ppr <- as.data.frame(t(ppr))
          ppr[!is.na(ppr) & ppr == ""] <- NA

          ddt <- cbind(ddt,ddr)
          ppt <- rbind(ppt,ppr)
        } else if (ncol(dd3) == 0){ #experiment does not exist
          ddt <- cbind(ddt,blank)
          # ppt <- rbind(ppt,pp3[NA,])
        }
        else{#no replicates
          ddt <- cbind(ddt,dd3)
          ppt <- rbind(ppt,pp3)
        }
      }#)
      ddt <- ddt[,-(seq_len(length(a))), drop = F] #ddt contains the final expression matrix for a single drug
      colnames(ddt) <- a

      exp.list[[cnt]] <- ddt
    }#)
    names(exp.list) <- drugs
    ppf <- pp2[FALSE,]
    #lapply(unique(ppt[,"dose_level"]), function(i){
    for (i in unique(ppt[,"dose_level"])) {
      if (verbose == TRUE ) {
        message(i)
      }
      #lapply(unique(ppt[,"duration"]), function(j){
      for (j in unique(ppt[,"duration"])) {
        if (verbose == TRUE) {
          message(j)
        }
        pp4 <- apply(ppt[ppt[,"dose_level"] == i & ppt[,"duration"] == j,,drop = F], 2, function(x) {
          x <- paste(unique(as.character(x[!is.na(x)])), collapse = "///")
          #if (is.na(x)){x <- paste("Exp ",,sep="")}
          return(x)
        })
        pp4 <- as.data.frame(t(pp4))
        pp4[!is.na(pp4) & pp4 == ""] <- NA
        ppf <- rbind(ppf,pp4)
      }#)
    }#)
    ppf <- as.data.frame(ppf,stringsAsFactors=FALSE)
    rownames(ppf) <- paste(ppf[,"dose_level"],";",ppf[,"duration"], sep = "")
    vec <- as.vector(colnames(exp.list[[1]]))
    ppf <- ppf[vec,]

    res <- SummarizedExperiment(assays = exp.list, rowData = ff, colData = ppf)

    return(res)
  }
