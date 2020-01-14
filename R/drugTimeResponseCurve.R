#' Compares viabilities at a given dose over different experimental durations
#'
#' This function generates a plot visualizing the relationship between gene
#'   expression, time and dose level for the selected tSet. The plot is generated
#'   with ggplot2 and can be customized using ggplot plot + function() syntax.
#'
#' @examples
#'   library(ggplot2)
#'
#'   # Default settings
#'   plot <- drugTimeResponseCurve(TGGATESsmall, cell_lines = "Hepatocyte",
#'   dose = c("Control", "Low", "Middle"), drugs = drugNames(TGGATESsmall)[6],
#'   duration = c("2", "8", "24"))
#'
#'   # Customize title, x/y labels, x/y limits, colour palette and define
#'   # custom ticks for x axis using the function argument ggplot2_args
#'   customizations <- list(labs(title= 'My Custom Title', ylab = 'The y-axis'),
#'                          xlim(c(2, 24)), ylim(c(99,105)),
#'                          scale_color_brewer(palette="Set1"),
#'                          scale_x_continuous(breaks=c(2, 8, 24),
#'                            labels = c("Two", "Eight", "Twenty-Four"))
#'                          )
#'
#'    if(interactive()) {
#'       drugTimeResponseCurve(TGGATESsmall, cell_lines = "Hepatocyte",
#'         dose = c("Control", "Low", "Middle"),
#'         drugs = drugNames(TGGATESsmall)[6], duration = c("2", "8", "24"),
#'         ggplot_args = customizations)
#'    }
#'
#'    # Customize the plot using standard ggplot2 syntax
#'    if(interactive()) {
#'       plot + labs(title= 'My Custom Title', ylab = 'The y-axis') +
#'         xlim(c(2, 24)) + ylim(c(99,105)) + scale_color_brewer(palette="Set1")
#'    }
#'
#' @param tSet \code{ToxicoSet} A ToxicoSet to be plotted in
#'   this figure
#' @param dose \code{character} A vector of dose levels to be included in the
#'   plot. Default to include all dose levels available for a drug. Must include
#'   at minimum two dose levels, one of witch is "Control".
#' @param drugs \code{character} A drugs or pair of drugs to be plotted.
#' @param duration \code{character} A vector of durations to include in the plot.
#' @param summarize_replicates \code{logical} If TRUE will average viability
#'   across replicates for each unique drug-dose-duration combination.
#' @param cell_lines \code{character} A vector of cell lines to include in the
#'   plot.
#' @param line_width \code{numeric} A number specifying the thickness of lines
#'   in the plot, as passed to size in geom_line(). Defaults to 1.
#' @param point_size \code{numeric} A number specifying how large points should
#'   be in the plot, as passed to size in geom_point(). Defaults to 2.5.
#' @param verbose \code{boolean} Should warning messages about the data passed
#'   in be printed?
#' @param ggplot_args \code{list} A list of ggplot2 functions which can be
#'   called using the plot + function() syntax. This allows arbitrary
#'   customization of the plot including changing the title, axis labels,
#'   colours, etc. Please see the included examples for basic usage or ggplot2
#'   documentation for advanced customization. Alternatively, you could assign
#'   the return value to a variable and add the customization yourself using
#'   plot + function().
#'
#' @return Plot of the viabilities for each drugs vs time of exposure
#'
#' @import ggplot2
#' @importFrom magrittr %<>%
#' @importFrom dplyr %>% filter group_by mutate
#' @importFrom tidyr gather
#'
#' @export
drugTimeResponseCurve <- function(
  tSet,
  duration = NULL,
  cell_lines = NULL,
  dose = NULL,
  drugs = NULL,
  summarize_replicates = TRUE,
  line_width = 1,
  point_size = 2.5,
  verbose=TRUE,
  ggplot_args=NULL
) {
  # Place tSet in a list if not already
  if (!is(tSet, "list")) {
    tSet <- list(tSet)
  }


  paramErrorChecker("drugTimeResponseCurve",
                    tSets = tSet, drugs = drugs, duration = duration,
                    cell_lines = cell_lines, dose = dose)

  ## TODO:: Throw warning if a dose level or time point is not available for
    # a specific drug

  # Subsetting the tSet based on parameter arguments
  tSet <- lapply(tSet, function(tSet) {
    suppressWarnings({subsetTo(tSet, mDataType = "rna", drugs = drugs,
                               duration = duration, cells = cell_lines)})
  })

  # Gather data for the plot
  plotData <- lapply(tSet, function(tSet) {
    sInfo <- sensitivityInfo(tSet)[, seq_len(4)]
    sValues <- tSet@sensitivity$raw[,,2]
    plotData <- cbind(sInfo, sValues)
    cols <- c('Low', 'Middle', 'High')
    colnames(plotData)[which(colnames(plotData) %in% c('doses1', 'doses2', 'doses3'))] <-
      cols[which(c('doses1', 'doses2', 'doses3') %in% colnames(plotData))]
    plotData %<>% gather('dose_level', 'viability', Control, Low, Middle, High)
  })

  for (data in plotData) {
    if (summarize_replicates) {
      data %<>% group_by(dose_level, duration_h) %>% mutate(viability = mean(viability))
      plot <- ggplot(as_tibble(data) %>% filter(replicate == 1),
                     aes(as.numeric(duration_h), viability, color = dose_level)) +
        geom_point(size = point_size) +
        geom_line(size = line_width)
    } else {
      plot <- ggplot(as_tibble(data), aes(as.numeric(duration_h), viability,
                                          color = dose_level,
                                          shape = as.factor(replicate),
                                          linetype = as.factor(replicate))) +
        geom_point(size = point_size) +
        geom_line(size = line_width)
    }
  }

  plot <- plot + labs(
    title = paste0("Drug Response Curve for ",
                   paste(drugs, collapse = " & "), " in ",
                   paste(cell_lines, collapse = " & "), collapse = " & "),
    color = "Dose Level",
    shape = "Replicate"
  ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14)
    ) +
    xlab("Duration (hrs)") +
    ylab("Viability (%)") +
    scale_x_continuous(breaks=as.numeric(duration), labels = duration)

  # Pass in any additional ggplot2 customizations
  if (!(is.null(ggplot_args))) {
    plot <- plot + ggplot_args
  }
  plot
}
