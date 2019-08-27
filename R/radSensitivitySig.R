#' Creates a signature representing the association between gene expression (or
#' other molecular profile) and radiation dose response, for use in radiation sensitivity
#' analysis.
#'
#' Given a RadioSet of the sensitivity experiment type, and a list of drugs,
#' the function will compute a signature for the effect gene expression on the
#' molecular profile of a cell. The function returns the estimated coefficient,
#' the t-stat, the p-value and the false discovery rate associated with that
#' coefficient, in a 3 dimensional array, with genes in the first direction,
#' drugs in the second, and the selected return values in the third.
#'
#' @examples
#' data(Cleveland_small)
#' rad.sensitivity <- radSensitivitySig(Cleveland_small, mDataType="rna",
#'              nthread=1, features = fNames(Cleveland_small, "rna")[1],
#'              radiation.types=radiationTypes(Cleveland_small))
#' print(rad.sensitivity)
#'
#' @param rSet [PharmacoSet] a PharmacoSet of the perturbation experiment type
#' @param mDataType [character] which one of the molecular data types to use
#'   in the analysis, out of dna, rna, rnaseq, snp, cnv
#' @param radiation.types [character] a vector of radiation.types for which to compute the
#'   signatures. Should match the names used in the PharmacoSet.
#' @param features [character] a vector of features for which to compute the
#'   signatures. Should match the names used in correspondant molecular data in PharmacoSet.
#' @param nthread [numeric] if multiple cores are available, how many cores
#'   should the computation be parallelized over?
#' @param returnValues [character] Which of estimate, t-stat, p-value and fdr
#'   should the function return for each gene?
#' @param sensitivity.measure [character] which measure of the radiation
#'   sensitivity should the function use for its computations? Use the
#'   sensitivityMeasures function to find out what measures are available for each PSet.
#' @param molecular.summary.stat What summary statistic should be used to
#'   summarize duplicates for cell line molecular profile measurements?
#' @param sensitivity.summary.stat What summary statistic should be used to
#'   summarize duplicates for cell line sensitivity measurements?
#' @param sensitivity.cutoff Allows to provide upper and lower bounds to
#'   sensitivity measures in the cases where the values exceed physical values
#'   due to numerical or other errors.
#' @param standardize [character] One of "SD", "rescale", or "none", for the form of standardization of
#'   the data to use. If "SD", the the data is scaled so that SD = 1. If rescale, then the data is scaled so that the 95%
#'   interquantile range lies in [0,1]. If none no rescaling is done.
#' @param verbose [boolean] 'TRUE' if the warnings and other infomrative message shoud be displayed
#' @param ... additional arguments not currently fully supported by the function
#' @return [list] a 3D array with genes in the first dimension, radiation.types in the
#'   second, and return values in the third.
#' @export
#' @import parallel

radSensitivitySig <- function(rSet,
 mDataType,
 radiation.types,
 features,
 sensitivity.measure = "AUC_recomputed",
 molecular.summary.stat = c("mean", "median", "first", "last", "or", "and"),
 sensitivity.summary.stat = c("mean", "median", "first", "last"),
 returnValues = c("estimate", "pvalue", "fdr"),
 sensitivity.cutoff, standardize = c("SD", "rescale", "none"),
 nthread = 1,
 verbose=TRUE, ...) {

  ### This function needs to: Get a table of AUC values per cell line / drug
  ### Be able to recompute those values on the fly from raw data if needed to change concentration
  ### Be able to choose different summary methods on fly if needed (need to add annotation to table to tell what summary method previously used)
  ### Be able to extract genomic data
  ### Run rankGeneDrugSens in parallel at the drug level
  ### Return matrix as we had before

  #sensitivity.measure <- match.arg(sensitivity.measure)
  molecular.summary.stat <- match.arg(molecular.summary.stat)
  sensitivity.summary.stat <- match.arg(sensitivity.summary.stat)
  standardize <- match.arg(standardize)

  dots <- list(...)
  ndots <- length(dots)


  if (!all(sensitivity.measure %in% colnames(sensitivityProfiles(rSet)))) {
    stop (sprintf("Invalid sensitivity measure for %s, choose among: %s", rSet@annotation$name, paste(colnames(sensitivityProfiles(rSet)), collapse=", ")))
  }

  if (!(mDataType %in% names(rSet@molecularProfiles))) {
    stop (sprintf("Invalid mDataType for %s, choose among: %s", rSet@annotation$name, paste(names(rSet@molecularProfiles), collapse=", ")))
  }
  switch (Biobase::annotation(rSet@molecularProfiles[[mDataType]]),
    "mutation" = {
      if (!is.element(molecular.summary.stat, c("or", "and"))) {
        stop ("Molecular summary statistic for mutation must be either 'or' or 'and'")
      }
    },
    "fusion" = {
      if (!is.element(molecular.summary.stat, c("or", "and"))) {
        stop ("Molecular summary statistic for fusion must be either 'or' or 'and'")
      }
    },
    "rna" = {
      if (!is.element(molecular.summary.stat, c("mean", "median", "first", "last"))) {
        stop ("Molecular summary statistic for rna must be either 'mean', 'median', 'first' or 'last'")
      }
    },
    "cnv" = {
      if (!is.element(molecular.summary.stat, c("mean", "median", "first", "last"))) {
        stop ("Molecular summary statistic for cnv must be either 'mean', 'median', 'first' or 'last'")
      }
    },
    "rnaseq" = {
      if (!is.element(molecular.summary.stat, c("mean", "median", "first", "last"))) {
        stop ("Molecular summary statistic for rna must be either 'mean', 'median', 'first' or 'last'")
      }},
      stop (sprintf("No summary statistic for %s has been implemented yet", Biobase::annotation(rSet@molecularProfiles[[mDataType]])))
      )

  if (!is.element(sensitivity.summary.stat, c("mean", "median", "first", "last"))) {
    stop ("Sensitivity summary statistic for sensitivity must be either 'mean', 'median', 'first' or 'last'")
  }

  if (missing(sensitivity.cutoff)) {
    sensitivity.cutoff <- NA
  }
  if (missing(radiation.types)){
    drugn <- radiationTypes(rSet)
  } else {
    drugn <- radiation.types
  }
  availcore <- parallel::detectCores()
  if ( nthread > availcore) {
    nthread <- availcore
  }

  if (missing(features)) {
    features <- rownames(featureInfo(rSet, mDataType))
  } else {
    fix <- is.element(features, rownames(featureInfo(rSet, mDataType)))
    if (verbose && !all(fix)) {
      warning (sprintf("%i/%i features can be found", sum(fix), length(features)))
    }
    features <- features[fix]
  }

  if(is.null(dots[["sProfiles"]])){
    drugpheno.all <- lapply(sensitivity.measure, function(sensitivity.measure) {

      return(t(summarizeSensitivityProfiles(rSet,
        sensitivity.measure = sensitivity.measure,
        summary.stat = sensitivity.summary.stat,
        verbose = verbose)))

    })} else {
      sProfiles <- dots[["sProfiles"]]
      drugpheno.all <- list(t(sProfiles))
    }

    dix <- is.element(drugn, do.call(colnames, drugpheno.all))
    if (verbose && !all(dix)) {
      warning (sprintf("Only %i/%i radiation types can be found", sum(dix), length(drugn)))
    }
    if (!any(dix)) {
      stop("None of the chosen radiation types were found in the dataset")
    }
    drugn <- drugn[dix]

    rSet@molecularProfiles[[mDataType]] <- summarizeMolecularProfiles(rSet = rSet,
      mDataType = mDataType,
      summary.stat = molecular.summary.stat,
      verbose = verbose)[features, ]

    if(!is.null(dots[["mProfiles"]])){
      mProfiles <- dots[["mProfiles"]]
      Biobase::exprs(rSet@molecularProfiles[[mDataType]]) <- mProfiles[features, colnames(rSet@molecularProfiles[[mDataType]]), drop = FALSE]

    }

    drugpheno.all <- lapply(drugpheno.all, function(x) {x[phenoInfo(rSet, mDataType)[ ,"cellid"], , drop = FALSE]})

    type <- as.factor(cellInfo(rSet)[phenoInfo(rSet, mDataType)[ ,"cellid"], "tissueid"])
    batch <- phenoInfo(rSet, mDataType)[, "batchid"]
    batch[!is.na(batch) & batch == "NA"] <- NA
    batch <- as.factor(batch)
    names(batch) <- phenoInfo(rSet, mDataType)[ , "cellid"]
    batch <- batch[rownames(drugpheno.all[[1]])]
    if (verbose) {
      message("Computing radiation sensitivity signatures...")
    }

    splitix <- parallel::splitIndices(nx = length(drugn), ncl = 1)
    splitix <- splitix[sapply(splitix, length) > 0]
    mcres <-  parallel::mclapply(splitix, function(x, drugn, expr, drugpheno, type, batch, standardize, nthread) {
      res <- NULL
      for(i in drugn[x]) {
        ## using a linear model (x ~ concentration + cell + batch)
        dd <- lapply(drugpheno, function(x) x[,i])
        dd <- do.call(cbind, dd)
        colnames(dd) <- seq_len(ncol(dd))
        if(!is.na(sensitivity.cutoff)) {
          dd <- factor(ifelse(dd > sensitivity.cutoff, 1, 0), levels=c(0, 1))
        }
        rr <- rankGeneDrugSensitivity(data=expr, drugpheno=dd, type=type, batch=batch, single.type=FALSE, standardize=standardize, nthread=nthread, verbose=verbose)
        res <- c(res, list(rr$all))
      }
      names(res) <- drugn[x]
      return(res)
    }, drugn=drugn, expr=t(molecularProfiles(rSet, mDataType)[features, , drop=FALSE]), drugpheno=drugpheno.all, type=type, batch=batch, nthread=nthread, standardize=standardize)

    res <- do.call(c, mcres)
    res <- res[!sapply(res, is.null)]
    drug.sensitivity <- array(NA,
      dim = c(nrow(featureInfo(rSet, mDataType)[features,, drop=FALSE]),
        length(res), ncol(res[[1]])),
      dimnames = list(rownames(featureInfo(rSet, mDataType)[features,]), names(res), colnames(res[[1]])))
    for(j in 1:ncol(res[[1]])) {
      ttt <- sapply(res, function(x, j, k) {
        xx <- array(NA, dim = length(k), dimnames = list(k))
        xx[rownames(x)] <- x[ , j, drop=FALSE]
        return (xx)
      },
      j = j,
      k = rownames(featureInfo(rSet, mDataType)[features,, drop = FALSE]))
      drug.sensitivity[rownames(featureInfo(rSet, mDataType)[features,, drop = FALSE]), names(res), j] <- ttt
    }

    drug.sensitivity <- RadioSig(drug.sensitivity, PSetName = rSetName(rSet), Call ="as.character(match.call())", SigType='Sensitivity')

    return(drug.sensitivity)
  }
