#' Conduct differential expression analysis using the limma R pacakge
#'
#' WARNING: This function can take a very long time to compute!
#'
#' @examples
#' if (interactive()) {
#' data(TGGATESsmall)
#' analysis <- computeLimmaDiffExpr(TGGATESsmall)
#' }
#'
#' @param object A [`ToxicoSet`] object with a molecular profile named 'rna'
#' @param buildTable [`logical`] Should the result of the eBayes function
#'  from limma be assembled into a data.table containing the result along
#'  with the gene, compound and durations names. Default it TRUE, otherwise
#'  this function with return the object produced by eBayes.
#'
#' @return A [`data.table`] containing the results the limma differential
#'  expression analysis comparing control vs each dose level for each compound
#'  within each duration.
#'
#' @import data.table
#' @import Biobase
#' @import limma
#' @importFrom SummarizedExperiment SummarizedExperiment coerce
#' @importFrom stats model.matrix model.frame
#' @importFrom BiocParallel bplapply
#'
#' @export
setMethod('computeLimmaDiffExpr', signature(object='ToxicoSet'), function(object, buildTable=TRUE) {

    ## TODO:: Add messaages to keep track of where the function execution is at

    # ---- 1. Get the required input data
    SE <- molecularProfilesSlot(object)$rna  # Extract the rna expression SummarizedExperiment
    eset <- as(SE, 'ExpressionSet')  # Coerce to an ExpressionSet
    eset$drugid <- make.names(eset$drugid)
    eset$cellid <- make.names(eset$cellid)

    # ---- 2. Extract the metadata needed to build the design matrix

    # Get the sample name, drug, dose and duration and cell type from the
    #   experiments phenotypic data
    targets <- as.data.frame(pData(eset)[, c("samplename", "cellid", "drugid",
                                             "dose_level", "duration")])
    colnames(targets) <- c('sample', 'cell', 'compound', 'dose', 'duration')
    # to prevent dropping numeric columns when converting to factors
    targets <- data.frame(lapply(targets, as.character))

    # ---- 3. Create the design matrix

    # Construct the design matrix with an intercept at 0
    # This follows the make the simplest design matrix possible strategy outlined
    # on page 37 of the limma user guide.
    hasMultipleCells <- length(unique(targets$cell)) > 1
    if (hasMultipleCells) {
        design <- model.matrix(~0 + compound:dose:duration:cell,
            data=model.frame(targets))
    } else {
        design <- model.matrix(~0 + compound:dose:duration,
            data=model.frame(targets))
    }


    # Make the names valid factor names and match with contrasts
    colnames(design) <- gsub(':', '_', colnames(design))

    # ---- 4. Fit a linear model based on design matrix

    # Every possible comparison between compound:dose:duration will
    # be fit in the linear model - this is less computationally efficient but
    # allows arbitrary comparisons to be extracted after the model is fit.
    # We chose this strategy to avoid the need to prespecify an experimental
    # design (less statistical planning required).
    fit <- lmFit(eset, design)

    # ---- 5. Setup the contrasts to performs our desired statistical tests
    setDT(targets)  # Convert to data.table

    columns <- 'duration'
    hasSharedControls <- name(object) %in% c('drugMatrix_rat', 'EMEXP2458')
    if (!hasSharedControls)
        columns <- c('compound', columns)
    if (hasMultipleCells)
        columns <- c('cell', columns)

    controls <- unique(targets[dose == 'Control',
                        .(controlLevels=paste0('compound', compound, '_dose',
                                               dose, '_duration', duration,
                                               '_cell', cell)),
                        by=c('cell', 'compound', 'duration')])
    levels <- unique(targets[dose != 'Control',
                      .(treatmentLevels=paste0('compound', compound, '_dose',
                                               dose, '_duration', duration,
                                               '_cell', cell)),
                      by=c('cell', 'compound', 'duration')])

    # fix levels if there is only one cell type
    if (!hasMultipleCells) {
        controls$controlLevels <- gsub('_cell[^_]*$', '',
            controls$controlLevels)
        levels$treatmentLevels <- gsub('_cell[^_]*$', '',
            levels$treatmentLevels)
    }

    contrastStrings <- unique(controls[levels,
                              .(contrasts=paste0(treatmentLevels, '-', controlLevels)),
                              on=columns, all=!hasSharedControls]$contrasts)


    # Make contrasts matrix for fitting statistical testing
    contrasts <- makeContrasts(contrasts=contrastStrings, levels=design)

    # Perfrom statistical tests for the specified contrasts
    contrFit <- contrasts.fit(fit, contrasts)

    # ---- 6. Predict coefficients using emperical Bayes moderation of SE
    # Generate a t-stat, moderated F-stat and log-odds of differential expression
    stats <- eBayes(contrFit)

    if (!buildTable) return(stats)

    # ---- 7. Assemble the results into a data.table object
    compoundNames <- make.names(drugInfo(object)$drugid)
    compounds <- drugInfo(object)$drugid
    ## TODO:: refactor this into muliple lapply statements!
    resultList <- lapply(contrastStrings, function(comparison) {
      # Disassmble contrasts into annotations for this statistical test
      comparisonNoLabels <- gsub('cell|compound|dose|duration',
        '', comparison)
      annotations <- unlist(strsplit(comparisonNoLabels, '-'))
      annotations <- strsplit(annotations, '_')

      if (!(annotations[[1]][3] == annotations[[2]][3]))
        stop("Duration mismatch between treatment and control!")

      # Get the compound_id based on the contrast name then get the compound name
      compound_id <- which(compoundNames %in% annotations[[1]][1])
      compound <- compounds[compound_id]

      # Get the tests per gene for the given comparison
      statCols <- c("logFC", "B", "P.Value", "adj.P.Val", "AveExpr")
      results <- topTreat(stats, coef=comparison, sort.by="none",
                          number='all',adjust.method="BH")[, statCols]
      statCols <- c("fold_change", "log_odds", "p_value", "fdr", "avg_expr")
      colnames(results) <- statCols

      DT <- data.table(
        "gene" = rownames(results),
        "compound" = rep(compound, nrow(results)),
        "dose" = rep(annotations[[1]][2], nrow(results)),
        "duration" = rep(annotations[[1]][3], nrow(results)),
        results
      )
      if (length(annotations[[1]]) > 3) {
          DT$cell <- rep(annotations[[1]][4], nrow(results))
          nonStatCols <- setdiff(colnames(DT), statCols)
          setcolorder(DT, c(nonStatCols, statCols))
      }
      return(DT)
    })
    analysis <- rbindlist(resultList)

    ## TODO:: Do we want to only return significant results? Will be a lot less rows

    # --- 8. Annotate and return the results
    return(analysis)
})

