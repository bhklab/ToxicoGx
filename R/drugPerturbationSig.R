#	Get Drug Perturbation Signatures from a ToxicoSet
###############################################################################
## Drug perturbation analysis
## create profiles before vs after drug for each drug 
###############################################################################

#' Creates a signature representing gene expression (or other molecular profile)
#' change induced by administrating a drug, for use in drug effect analysis.
#' 
#' Given a Toxicoset of the perturbation experiment type, and a list of drugs,
#' the function will compute a signature for the effect of drug concentration on
#' the molecular profile of a cell. The algorithm uses a regression model which
#' corrects for experimental batch effects, cell specific differences, and
#' duration of experiment to isolate the effect of the concentration of the drug
#' applied. The function returns the estimated coefficient for concentration,
#' the t-stat, the p-value and the false discovery rate associated with that
#' coefficient, in a 3 dimensional array, with genes in the first direction,
#' drugs in the second, and the selected return values in the third.
#' 
#' @examples
#' data(TGGATES_small)
#' drug.perturbation <- drugPerturbationSig(TGGATES_small, mDataType="rna", nthread=1)
#' print(drug.perturbation)
#' 
#' @param tSet [ToxicoSet] a ToxicoSet of the perturbation experiment type
#' @param mDataType [character] which one of the molecular data types to use
#'   in the analysis, out of dna, rna, rnaseq, snp, cnv (only rna currently supported)
#' @param drugs [character] a vector of drug names for which to compute the
#'   signatures. Should match the names used in the ToxicoSet.
#' @param cells [character] a vector of cell names to use in computing the 
#'   signatures. Should match the names used in the ToxicoSet.
#' @param features [character] a vector of features for which to compute the
#'   signatures. Should match the names used in correspondant molecular data in ToxicoSet.
#' @param duration [numeric] a vector of experiment durations for which to
#'   compute the signatures.
#' @param nthread [numeric] if multiple cores are available, how many cores
#'   should the computation be parallelized over?
#' @param returnValues [character] Which of estimate, t-stat, p-value and fdr
#'   should the function return for each gene drug pair?
#' @param verbose [bool] Should diagnostive messages be printed? (default false)
#' @return [list] a 3D array with genes in the first dimension, drugs in the
#'   second, and return values in the third.
#' @export

drugPerturbationSig <- function(tSet, mDataType, drugs, cells, features, duration, nthread=1, returnValues=c("estimate","tstat", "pvalue", "fdr"), verbose=FALSE){
  availcore <- parallel::detectCores()
  if ( nthread > availcore) {
    nthread <- availcore
  }
  options("mc.cores"=nthread)
  if(!missing(cells)){
    if(!all(cells%in%cellNames(tSet))){
      stop("The cell names should match to the names used in cellNames(tSet)")
    }
    tSet <- subsetTo(tSet, cells=cells)
  }
  if (mDataType %in% names(tSet@molecularProfiles)) {
    #eset <- tSet@molecularProfiles[[mDataType]]
    if(Biobase::annotation(tSet@molecularProfiles[[mDataType]])!="rna"){
      stop(sprintf("Only rna data type perturbations are currently implemented"))
    }
  } else {
    stop (sprintf("This tSet does not have any molecular data of type %s, choose among: %s", mDataType), paste(names(tSet@molecularProfiles), collapse=", "))
  }
  
  
  if (missing(drugs)) {
    drugn <- drugNames(tSet)
  } else {
    drugn <- drugs
  }
  dix <- is.element(drugn, phenoInfo(tSet, mDataType)[ , "drugid"])
  if (verbose && !all(dix)) {
    warning (sprintf("%i/%i drugs can be found", sum(dix), length(drugn)))
  }
  if (!any(dix)) {
    stop("None of the drugs were found in the dataset")
  }
  drugn <- drugn[dix]
  
  if (missing(features)) {
    features <- rownames(featureInfo(tSet, mDataType))
  } else {
    fix <- is.element(features, rownames(featureInfo(tSet, mDataType)))
    if (verbose && !all(fix)) {
      warning (sprintf("%i/%i features can be found", sum(fix), length(features)))
    }
    features <- features[fix]
  }
  
  samples <- rownames(phenoInfo(tSet, mDataType)[phenoInfo(tSet, mDataType)[,"duration"] %in% duration & phenoInfo(tSet, mDataType)[,"drugid"] %in% drugs,])
  # splitix <- parallel::splitIndices(nx=length(drugn), ncl=nthread)
  # splitix <- splitix[sapply(splitix, length) > 0]
  mcres <- lapply(drugn, function(x, exprs, sampleinfo) {
    res <- NULL
    i = x 
    ## using a linear model (x ~ concentration + cell + batch + duration)
    res <- rankGeneDrugPerturbation(data=exprs, drug=i, drug.id=as.character(sampleinfo[ , "drugid"]), drug.concentration=as.numeric(sampleinfo[ , "concentration"]), type=as.character(sampleinfo[ , "cellid"]), xp=as.character(sampleinfo[ , "xptype"]), batch=as.character(sampleinfo[ , "batchid"]), duration=as.character(sampleinfo[ , "duration"]) ,single.type=FALSE, nthread=nthread, verbose=FALSE)$all[ , returnValues, drop=FALSE]
    res <- list(res)
    names(res) <- i
    return(res)
  }, exprs=t(molecularProfiles(tSet, mDataType)[features, samples, drop=FALSE]), sampleinfo=phenoInfo(tSet, mDataType)[phenoInfo(tSet, mDataType)[,"duration"] == duration & phenoInfo(tSet, mDataType)[,"drugid"] %in% drugs,])
  res <- do.call(c, mcres)
  res <- res[!sapply(res, is.null)]
  drug.perturbation <- array(NA, dim=c(nrow(featureInfo(tSet, mDataType)[features,, drop=FALSE]), length(res), ncol(res[[1]])), dimnames=list(rownames(featureInfo(tSet, mDataType)[features,,drop=FALSE]), names(res), colnames(res[[1]])))
  for(j in 1:ncol(res[[1]])) {
    ttt <- sapply(res, function(x, j, k) {
      xx <- array(NA, dim=length(k), dimnames=list(k))
      xx[rownames(x)] <- x[ , j, drop=FALSE]
      return (xx)
    }, j=j, k=rownames(featureInfo(tSet, mDataType)[features,, drop=FALSE]))
    drug.perturbation[rownames(featureInfo(tSet, mDataType)[features,, drop=FALSE]), names(res), j] <- ttt
  }
  
  drug.perturbation <- ToxicoSig(drug.perturbation, tSetName = tSetName(tSet), Call = as.character(match.call()), SigType='Perturbation')
  
  return(drug.perturbation)
}
