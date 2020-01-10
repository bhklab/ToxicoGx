#' Compares gene expression for a specificed set of features over specific
#'   drug dosages vs time
#'
#' This function generates a plot visualizing the relationship between gene
#'   expression, time and dose level for the selected tSet. The plot is generated
#'   with ggplot2 and can be customized using ggplot plot + function() syntax.
#'
#' @examples
#'
#' if (interactive()) {
#'   drugGeneResponseCurve(TGGATESsmall, dose = c("Control", "Low", "Middle"),
#'   mDataTypes="rna", drug = drugNames(TGGATESsmall)[1],
#'   duration = c("2", "8", "24"), features = "ENSG00000000003_at")
#' }
#'
#' @param tSet \code{ToxicoSet} A ToxicoSet to be plotted in this graph. Currently
#'   only a single tSet is supported.
#' @param dose \code{character} A vector of dose levels to be included in the
#'   plot. Default to include all dose levels available for a drug. If you specify
#'   more than two features you may only pass in up to two dose levels.
#' @param mDataTypes \code{vector} A vector specifying the molecular data types to
#'   include in this plot. Defaults to the first mDataType if not specified.ex
#'   This release version only accepts one mDataType, more to be added in
#'   forthcoming releases.
#' @param features \code{character} A vector of feature names to include in the plot.
#'   If you specify more than two dose levels, you may only pass in up to two features.
#' @param drug \code{character} A drug name to include in this plot.
#'   See drugNames(tSet) for a list of options.
#' @param duration \code{character} A vector of durations to include in the plot.
#' @param summarize_replicates \code{logical} If TRUE will average viability
#'   across replicates for each unique drug-dose-duration combination.
#' @param cell_lines \code{character} A vector of cell lines to include in the plot.
#' @param line_width \code{numeric} A number specifying the thickness of lines
#'   in the plot, as passed to size in geom_line(). Defaults to 1.
#' @param point_size \code{numeric} A number specifying how large points should
#'   be in the plot, as passed to size in geom_point(). Defaults to 2.5.
#' @param ggplot_args \code{list} A list of ggplot2 functions which can be called using the plot + function()
#'   syntax. This allows arbitrary customization of the plot including changing the title, axis labels, colours, etc.
#'   Please see the included examples for basic usage or ggplot2 documentation for advanced customization.
#' @param verbose \code{boolean} Should warning messages about the data passed
#'   in be printed?
#'
#' @return Plot of the viabilities for each drug vs time of exposure
#'
#' @import data.table
#' @import ggplot2
#' @importFrom tibble as_tibble

#' @export
drugGeneResponseCurve <- function(
  tSet,
  duration = NULL,
  cell_lines = NULL,
  mDataTypes = NULL,
  features = NULL,
  dose = NULL,
  drug = NULL,
  summarize_replicates = TRUE,
  line_width = 1,
  point_size = 2.5,
  ggplot_args = NULL,
  verbose=TRUE
) {

  # Place tSet in a list if not already
  if (!is(tSet, "list")) { tSet <- list(tSet) }

  ## Tempary warnings until function is finished
  if (length(tSet) > 1) { stop("This function currently only supports one tSet per plot...")}
  if (length(drug) > 1) { stop("This function currently only supports one drug per plot...")}
  if (length(mDataTypes) > 1) {stop("This function currently only supports one molecular data type per plot...")}
  if (length(features) > 2) { if (length(dose) > 2) { stop("To plot more than one feature, please specify only up to two dose levels...")}}
  if (length(dose) > 2) { if (length(features) > 2) { stop("To plot more than one dose level, please specify up to two molecular feature...")}}

  # Deal with controls (i.e., treated with DMSO)
  if (any(vapply(tSet, function(tSet) { names(tSet) == "drugMatrix"}, FUN.VALUE = logical(1)))) {
    drug <- c("DMSO", drug)
  }

  ## TODO:: Generalize this to work with multiple data types
  if (missing(mDataTypes)) { mDataTypes <- names(tSet[[1]]@molecularProfiles) }

  if (is.null(features)) {
    features <- lapply(tSet, function(tSet) {
      rownames(featureInfo(tSet, "rna"))[seq_len(5)]
    })
  }

  if (missing(cell_lines)) {cell_lines <- unique(phenoInfo(tSet[[1]], mDataTypes[1])$cellid)}
  if (length(cell_lines) > 1) { stop("Only one cell type per plot is currently supported...")}

  # Places features in list if not already
  if (!is(features, "list")) {
    features <- list(features)
  }
  names(features) <- vapply(tSet, function(x) names(x), FUN.VALUE = character(1))

  # Subsetting the tSet based on parameter arguments
  tSet <- lapply(tSet, function(tSet) {
    suppressWarnings({ToxicoGx::subsetTo(tSet, mDataType = mDataTypes, drugs = drug,
                       duration = duration, features = unique(unlist(features)))})
  })

  # Get only the dose levels available for that drug
  dose <- intersect(dose, as.character(unique(phenoInfo(tSet[[1]], "rna")$dose_level)))

  # Gather the plot data
  plotData <- lapply(tSet, function(tSet) {
    m <- lapply(mDataTypes, function(mDataType) {
      mProf <- molecularProfiles(tSet, mDataType)
      list(
        "data" = data.table(
          mProf,
          keep.rownames = TRUE
        ),
        "pInfo" = data.table(phenoInfo(tSet, mDataType))
      )
    })
    names(m) <- mDataTypes; m
  })
  names(plotData) <-  vapply(tSet, function(x) names(x), FUN.VALUE = character(1))

  #### Assembling the plot data ####
  d <- plotData[[1]][[1]]$data
  pInfo <- plotData[[1]][[1]]$pInfo
  data <- melt.data.table(d, id.vars = 1, variable.factor = FALSE)
  colnames(data) <- c("feature", "samplename", "expression")
  data[, samplename := as.numeric(samplename)]
  pInfo[, samplename := as.numeric(samplename)]
  fInfo <- data.table(featureInfo(tSet[[1]], "rna"))
  colnames(fInfo)[2] <- "feature"

  plotData <- merge(data, pInfo[, .(samplename, individual_id,
                                    drugid, dose_level, duration)],
                    by = "samplename")
  plotData <- merge(plotData, fInfo[, .(Symbol, feature)], by = "feature")
  plotData[, dose_level := as.factor(dose_level)]
  plotData[dose_level == "Control",
                       expression := mean(expression),
                       by = .(dose_level, duration, Symbol)]
  max_rep <- max(plotData[dose_level != 'Control', unique(individual_id)])
  plotData <- plotData[individual_id %in% seq_len(max_rep), .SD, by = .(dose_level, duration)]
  plotData <- plotData[dose_level %in% dose, .SD]

  #### Rendering the plot ####
  if (summarize_replicates == FALSE) {
    plot <- ggplot(as_tibble(plotData),
           aes(x = as.numeric(duration),
               y = expression,
               color = dose_level,
               linetype = as.factor(individual_id),
               shape = Symbol,
               group = interaction(dose_level, individual_id, Symbol))) +
      geom_line(size = line_width) +
      geom_point(size = point_size)
  } else {
    plotData <- plotData[, expression := mean(expression), by = .(dose_level, duration, Symbol)][individual_id == 1]
    plot <- ggplot(plotData, aes(as.numeric(duration), expression, color = dose_level,)) +
      geom_line(aes(linetype = Symbol), size = line_width) +
      geom_point(size = point_size)
  }
  plot <- plot +
    labs(
      title = paste0("Drug Gene Response Curve for ", paste(drug, collapse = " & "), " in ", paste(cell_lines, collapse = " & "), collapse = " & "),
      color = "Dose Level",
      linetype = "Replicate",
      shape = "Feature"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14)
    ) + xlab("Duration (hrs)") +
    ylab("Expression") +
    scale_x_continuous(breaks=as.numeric(duration), labels = duration)

  if (!is.null(ggplot_args)) {
    plot <- plot + ggplot_args
  }
  plot
}
