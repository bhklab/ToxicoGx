#' Drug perturbation analysis
#'
#' Creates a signature representing gene expression (or other molecular profile)
#' change induced by administrating a drug, for use in drug effect analysis.
#'
#' Given a Toxicoset of the perturbation experiment type, and a character vector
#' of drugs, the function will compute a signature for the effect of
#' drug concentration on the molecular profile of a cell. The algorithm uses a
#' regression model which corrects for experimental batch effects, cell specific
#' differences, and duration of experiment to isolate the effect of the
#' concentration of the drug applied. The function returns the estimated
#' coefficient for concentration, the t-stat, the p-value and the false
#' discovery rate associated with that coefficient, in a 3 dimensional array,
#' with genes in the first direction, drugs in the second, and the selected
#' return values in the third.
#'
#' @examples
#' if (interactive()) {
#' data(TGGATESsmall)
#' drug.perturbation <- drugPerturbationSig(TGGATESsmall, mDataType="rna",
#'  features = head(fNames(TGGATESsmall, "rna")), nthread=1)
#' }
#'
#' @param tSet \code{ToxicoSet} a ToxicoSet of the perturbation experiment type
#' @param mDataType \code{character} which one of the molecular data types to use
#'   in the analysis, out of dna, rna, rnaseq, snp, cnv (only rna currently supported)
#' @param drugs \code{character} a vector of drug names for which to compute the
#'   signatures. Should match the names used in the ToxicoSet.
#' @param cell_lines \code{character} a vector of cell names to use in computing the
#'   signatures. Should match the names used in the ToxicoSet.
#' @param features \code{character} a vector of features for which to compute the
#'   signatures. Should match the names used in correspondant molecular data in ToxicoSet.
#' @param duration \code{character} a vector of experiment durations for which to inlcude in the
#'   computed the signatures.
#' @param dose \code{character} a vector of dose levels to include in the results
#' @param nthread \code{numeric} if multiple cores are available, how many cores
#'   should the computation be parallelized over?
#' @param returnValues \code{character} Which of estimate, t-stat, p-value and fdr
#'   should the function return for each gene drug pair
#' @param verbose \code{bool} Should diagnostive messages be printed? (default false)
#'
#' @return \code{ToxicoSig} An object composed of a 3D array with genes in the
#'   first dimension, drugs in the second, and return values in the third.
#'
#' @importFrom tibble as_tibble
#' @import ggplot2
#'
#' @export
#'
drugPerturbationSig <- function(
  tSet, mDataType,
  drugs = NULL,
  cell_lines = NULL,
  features = NULL,
  duration = NULL,
  dose = NULL,
  nthread=1,
  returnValues=c("estimate","tstat", "pvalue", "fdr"),
  verbose=FALSE
){

  # ALLOCATE CORES FOR PARALLEL PROCESSING
  availcore <- parallel::detectCores()
  if ( nthread > availcore) {
    nthread <- availcore
  }
  options("mc.cores" = nthread)

  ## MISSING VALUE HANDLING FOR PARAMETERS
  # Get named list of defualt values for missing parameters
  argDefaultList <-
    paramMissingHandler("drugPerturbationSig", tSet = tSet, mDataType = mDataType,
                        drugs = drugs, cell_lines = cell_lines, features = features,
                        duration = duration, dose = dose)

  # Assign any missing parameter default values to function environment
  if (length(argDefaultList) > 0) {
    for (idx in seq_along(argDefaultList)) {
      assign(names(argDefaultList)[idx], argDefaultList[[idx]])
    }
  }

  # ERROR HANDLING FOR PARAMETERS
  paramErrorChecker("drugPerturbationSig", tSet = tSet,
                    mDataType = mDataType, cell_lines = cell_lines,
                    drugs = drugs, features = features,
                    duration = duration, dose = dose)

  returnValues <- match.arg(returnValues, several.ok = TRUE)

  # SUBSET tSET BASED ON PARAMETERS
  tSetSubsetOnParams <-
    suppressWarnings(subsetTo(tSet, mDataType = mDataType, cells = cell_lines, drugs = drugs,
             features = features, duration = duration))

  # SUBSET SAMPLES BASED ON DOSE
  samples <- rownames(phenoInfo(tSetSubsetOnParams, mDataType)[which(phenoInfo(tSetSubsetOnParams, mDataType)$dose %in% dose),])

  # LOOP OVER DRUGS TO CALCULATE PER DRUG SUMMARY STATISTICS
  mcres <- lapply(drugs, function(x, exprs, sampleinfo) {

    # Subset to correct drugs
    exprs <- exprs[which(sampleinfo[ , "drugid"] %in% x),]
    sampleinfo <- sampleinfo[which(sampleinfo[ , "drugid"] %in% x),]

    # Warning that rankGeneDrugPerturbation will return a matrix of NAs for this drug
    if (length(unique(as.character(sampleinfo[, "xptype"]))) < 2) {
      warning(paste0("There are only controls available at dose levels ",
                     paste(dose, collapse = " ") ," for ", x, ",
                     summary statistics for this drug will be excluded for the results.
                     Adding another dose level will likely generate results."))
    }
    res <- NULL
    i <- x

    ## using a linear model (x ~ concentration + cell + batch + duration)
    res <- ToxicoGx::rankGeneDrugPerturbation(
      data = exprs, drug = x, drug.id = as.character(sampleinfo[ , "drugid"]),
      drug.concentration = as.numeric(sampleinfo[ , "concentration"]),
      type = as.character(sampleinfo[ , "cellid"]),
      xp = as.character(sampleinfo[ , "xptype"]),
      batch = as.character(sampleinfo[ , "batchid"]),
      duration = as.character(sampleinfo[ , "duration"]) ,
      single.type = FALSE, nthread = nthread,
      verbose = FALSE)$all[ , returnValues, drop = FALSE]
    res <- list(res)
    names(res) <- i
    return(res)
  }, exprs = t(molecularProfiles(tSetSubsetOnParams, mDataType)[features, samples, drop=FALSE]),
     sampleinfo = phenoInfo(tSetSubsetOnParams, mDataType)[which(phenoInfo(tSetSubsetOnParams, mDataType)$samplename %in% samples), ]
  )

  # ASSEMBLE RESULTS TO BE INCLUDED IN TOXICOSIG OBJECT
  res <- do.call(c, mcres)
  res <- res[!sapply(res, is.null)]
  drug.perturbation <- array(NA, dim = c(nrow(featureInfo(tSet, mDataType)[features,, drop = FALSE]), length(res), ncol(res[[1]])), dimnames = list(rownames(featureInfo(tSet, mDataType)[features,,drop = FALSE]), names(res), colnames(res[[1]])))
  for (j in seq_len(ncol(res[[1]]))) {
    ttt <- sapply(res, function(x, j, k) {
      xx <- array(NA, dim = length(k), dimnames = list(k))
      xx[rownames(x)] <- x[ , j, drop = FALSE]
      return(xx)
    }, j = j, k = rownames(featureInfo(tSet, mDataType)[features,, drop = FALSE]))
    drug.perturbation[rownames(featureInfo(tSet, mDataType)[features,, drop = FALSE]), names(res), j] <- ttt
  }

  # CREATE TOXICOSIG OBJECT
  drug.perturbation <- ToxicoGx::ToxicoSig(drug.perturbation, tSetName = cSetName(tSet), Call = as.character(match.call()), SigType = 'Perturbation')

  # RETURN TOXICOSIG OBJECT
  return(drug.perturbation)
}
