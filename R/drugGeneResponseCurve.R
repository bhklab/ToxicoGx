#' Compares gene expression for a specificed set of features over specific
#'   drug dosages vs time
#'
#' This function generates a plot visualizing the relationship between gene
#'   expression, time and dose level for the selected tSet.
#'
#' @examples
#'
#' if (interactive()) {
#' drugGeneResponseCurve(TGGATESsmall, dose = c("Control", "Low", "Middle"), mDataTypes="rna", drug = drugNames(TGGATESsmall)[1], duration = c("2", "8", "24"), features = "ENSG00000000003_at")
#' }
#'
#' @param tSet \code{ToxicoSet} A ToxicoSet to be plotted in this graph. Currently
#'   only a single tSet is supported, passing more may results in errors.
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
#' @param cell.lines \code{character} A vector of cell lines to include in the plot.
#'   Currently limited to one cell lines per plot with plans to add support for
#'   more in upcoming releases.
# @param xlim \code{numeric} A vector of minimum and maximum values for the x-axis
#'   of the returned plot.
# @param ylim \code{numeric} A vector of minimum and miximum values for the y-axis
#   of the returned plot.
# @param title \code{character} A string containing the desired plot name. If excluded
#'   a title wil be generated automatically.
# @param legend.loc \code{character} The location of the legend as passed to the plot()
#'   function from base graphics. Suggested values are either "topleft" or
#'   "topright", with the latter as the default.
# @param mycol `vector`` A vector of length equal to the product of the
#'   number of drugs, features and doses passed to the function. Takes colour
#'   arguments as passed to `col` parameter in the `plot()` function.
#'   Default palette is used when unspecified.
#' @param summarize.replicates \code{logical} If true will take the average of all
#'  replicates at each time point per gene and duration. This release has not
#'  yet implemented this feature.
# @param lwd \code{numeric} The line width to plot width
# @param cex \code{numeric} The cex parameter passed to plot. Controls the size of
#'   plot points and the font size of the legend and defaults to 0.7.
# @param cex.main \code{numeric} The cex.main parameter passed to plot,
#'   controls the size of the titles and defaults to 1.0.
# @param trunc \code{bool} Should the viability values be truncated to lie in
#   \code{0-100} before doing the fitting
# @param verbose \code{boolean} Should warning messages about the data passed
#'   in be printed?
#'
#' @return Plot of the viabilities for each drug vs time of exposure
#'
#' @import RColorBrewer
#' @importFrom graphics plot rect points lines legend
#' @importFrom grDevices rgb
#' @importFrom magicaxis magaxis
#' @import dplyr
#' @import data.table
#' @importFrom magrittr %<>%
#' @importFrom rlist list.map list.select list.filter list.search
#'
#' @export
drugGeneResponseCurve <- function(
  tSet,
  duration,
  cell.lines,
  mDataTypes,
  features = NULL,
  dose,
  drug,
  summarize.replicates = TRUE,
  xlim=c(0, 24),
  ylim=c(0, 15),
  mycol,
  title,
  lwd = 1.5,
  cex = 1,
  cex.main = 0.9,
  legend.loc = "topright",
  verbose=TRUE
) {

  # Place tSet in a list if not already
  if (!is(tSet, "list")) { tSet <- list(tSet) }

  ## Tempary warnings until function is finished
  if (length(tSet) > 1) { warning("Multiple tSet plotting has not been tested in this release...")}
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

  if (missing(cell.lines)) {cell.lines <- unique(phenoInfo(tSet[[1]], mDataTypes[1])$cellid)}
  if (length(cell.lines) > 1) { stop("Only one cell type per plot is currently supported...")}

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
  dose <- intersect(dose, unique(phenoInfo(tSet[[1]], "rna")$dose_level))

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

  plotData <- merge(data, pInfo[, .(samplename, individual_id, drugid, dose_level, duration)], by = "samplename")
  plotData <- merge(plotData, fInfo[, .(Symbol, feature)], by = "feature")
  plotData[, dose_level := as.factor(dose_level)]
  plotData[dose_level == "Control",
                       expression := mean(expression),
                       by = .(dose_level, duration, Symbol)]
  max_rep <- max(plotData[dose_level != 'Control', unique(individual_id)])
  plotData <- plotData[individual_id %in% seq_len(max_rep), .SD, by = .(dose_level, duration)]

  #### Rendering the plot ####
  if (summarize.replicates == FALSE) {
    ggplot(as_tibble(plotData),
           aes(x = sort(as.numeric(duration)),
               y = expression,
               color = dose_level,
               linetype = as.factor(individual_id),
               shape = Symbol,
               group = interaction(dose_level, individual_id, Symbol))) +
      geom_line() +
      geom_point() +
      labs(
        title = paste0("Drug Gene Response Curve for ", paste(drug, collapse = " & "), " in ", paste(cell.lines, collapse = " & "), collapse = " & "),
        color = "Dose Level",
        linetype = "Replicate",
        shape = "Feature"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14)
      ) + xlab("Duration (hrs)") +
      ylab("Expression")
  } else {
    plotData <- plotData[, expression := mean(expression), by = .(dose_level, duration, Symbol)][individual_id == 1]
    ggplot(plotData, aes(as.numeric(duration), expression)) +
      geom_line(aes(color = dose_level, linetype = Symbol), size = 1) +
      geom_point(aes(color = dose_level, shape = interaction(dose_level, Symbol)), size = 2) +
      labs(
          title = paste0("Drug Gene Response Curve for ", paste(drug, collapse = " & "), " in ", paste(cell.lines, collapse = " & "), collapse = " & "),
          color = "Dose Level",
          linetype = "Feature"
        ) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 16)
        ) +
      xlab("Duration (hrs)") +
      ylab("Expression")
  }
}
