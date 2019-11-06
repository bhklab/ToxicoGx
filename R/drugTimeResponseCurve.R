#' Compares viabilities at a given dose over different experimental duration
#'
#' Description of this function
#'
#' @examples
#' if (interactive()) {
#'   ToxicoGx::drugTimeResponseCurve(TGGATESsmall, cellline = "Hepatocyte", dose = c("Control", "Low", "Middle"), drug = drugNames(TGGATESsmall)[1], duration = c("2", "8", "24"))
#' }
#'
#' @param tSet \code{ToxicoSet} A ToxicoSet to be plotted in
#'   this graph.
#' @param dose \code{character} A vector of dose levels to be included in the
#'   plot. Default to include all dose levels available for a drug. Must include
#'   at minimum two dose levels, one of witch is "Control".
#' @param drug \code{character} A vector of drugs to be included in this plot.
#' @param duration \code{character} A vector of durations to include in the plot.
#' @param cell.lines \code{character} A vector of cell lines to include in the plot.
#' @param viability_as_pct \code{logical} A vector specifying if viabilities should
#'   be plotted as a percentage. Defaults to TRUE.
#' @param xlim \code{numeric} A vector of minimum and maximum values for the x-axis
#'   of the returned plot.
#' @param ylim \code{numeric} A vector of minimum and miximum values for the y-axis
#'   of the returned plot.
#' @param title \code{character} A string containing the desired plot name. If excluded
#'   a title wil be generated automatically.
#' @param legend.loc \code{character} The location of the legend as passed to the
#'   legends function from base graphics. Recommended values are 'topright' or
#'   'topleft'. Default is 'topleft'.
#' @param mycol \code{vector} A vector of length equal to the lenth of the tSet
#'   argument specifying which RColorBrewer colour to use per tSet. Default
#'   colours will be used if this parameter is excluded.
#' @param summarize.replicates \code{logical} If true will take the average of all
#'  replicates at each time point per gene and duration. This release has not
#'  yet implemented this feature.
#' @param x.custom.ticks \code{vector} A numeric vector of the distance between major
#'   and minor ticks on the x-axis. If excluded ticks appear only where duration
#'   values are specified.
#' @param lwd \code{numeric} The line width to plot width
#' @param cex \code{numeric} The cex parameter passed to plot
#' @param cex.main \code{numeric} The cex.main parameter passed to plot, controls
#' the size of the titles
# @param trunc \code{bool} Should the viability values be truncated to lie in \code{0-100}
# before doing the fitting
#' @param verbose \code{boolean} Should warning messages about the data passed in be printed?
#'
#' @return Plot of the viabilities for each drug vs time of exposure
#'
#' @import RColorBrewer
#' @importFrom graphics plot rect points lines legend axis
#' @importFrom grDevices rgb
#' @importFrom magicaxis magaxis
#'
#' @export
drugTimeResponseCurve <- function(
  tSet,
  duration = NULL,
  cell.lines = NULL,
  dose = NULL,
  drug = NULL,
  summarize.replicates = TRUE,
  viability_as_pct = TRUE,
  xlim = NULL,
  ylim = NULL,
  mycol,
  x.custom.ticks = NULL,
  title,
  lwd = 1,
  cex = 0.5,
  cex.main = 0.9,
  legend.loc = "topright",
  verbose=TRUE
) {

  # Place tSet in a list if not already
  if (!is(tSet, "list")) {
    tSet <- list(tSet)
  }

  ## TODO:: REMOVE ::: BEFORE CRAN SUBMISSION!
  ToxicoGx:::paramErrorChecker("drugTimeResponseCurve",
                    tSets = tSet, duration = duration, cell.lines = cell.lines,
                    dose = dose)

  ## TODO:: Make this function work with multiple tSet
  ## TODO:: Make this function work with multiple drugs
  ## TODO:: Throw warning if a dose level or time point is not available for a specific drug
  ## TODO:: Add logic to handle viability_as_pct = FALSE

  # Subsetting the tSet based on parameter arguments
  tSet <- lapply(tSet, function(tSet) {
    subsetTo(tSet, mDataType = "rna", drugs = drug, duration = duration)
  })

  # Extracting the data required for plotting into a list of data.frames
  plotData <- lapply(tSet, function(tSet) {
    cbind(
      tSet@sensitivity$raw[,,2][, which(c("Control", "Low", "Middle", "High") %in% dose)],
      ToxicoGx::sensitivityInfo(tSet)[, c("duration_h", "replicate")]
    )
  })

  # Extractiing the duration values for each row of plotData
  times <- lapply(plotData, function(data) {
    lapply(seq_along(unique(data$replicate)), function(idx) {
      as.numeric(data$duration_h)[which(data$replicate == idx)]
    })
  })

  # Assembling the legend names for each line to be plotted
  legendValues <- lapply(seq_along(plotData), function(d_idx){
    lapply(dose, function(level) {
      lapply(seq_along(unique(plotData[[d_idx]]$replicate)), function(r_idx) {
        paste(drug[d_idx], level, r_idx, sep = '_')
      })
    })
  })

  # Extracting the viability values for each row of plotData
  responses <- lapply(plotData, function(data) {
    lapply(cell.lines, function(cell){
      lapply(drug, function(drug){
        lapply(dose, function(level) {
          lapply(seq_along(unique(data$replicate)), function(idx) {
            as.vector(data[ , which(c("Control", "Low", "Middle", "High") %in% level)])[which(data$replicate == idx)]
          })
        })
      })
    })
  })

  # Summarizing replicate values
  if (summarize.replicates == TRUE) {
    responses <- lapply(seq_along(tSet), function(t_idx) {
      lapply(seq_along(dose), function(d_idx) {
        responseVals <- NULL
        responseVect <- unlist(responses[[t_idx]][[d_idx]])
        for (time in seq_along(unique(unlist(times[[t_idx]])))) {
          responseVals <- c(responseVals, mean(responseVect[time], responseVect[time + length(unique(unlist(times[[t_idx]])))]))
        }
        responseVals
      })
    })
    # Take unique values of all time replicates and place into a list
    times <- list(unique(unlist(times)))
    legendValues <- lapply(legendValues, function(legendLevel) {
      lapply(legendLevel, function(legendName) {unique(gsub("_[^_]*$", "", unlist(legendName)))})
    })

  }

  # Set x and y axis ranges based on time and viability values
  time.range <- c(min(unlist(unlist(times))), max(unlist(unlist(times))))
  viability.range <- c(min(unlist(responses, recursive = TRUE)), max(unlist(responses, recursive = TRUE)))
  for (i in seq_along(tSet)) {
    ## TODO:: Generalize this to n replicates
    time.range <- c(min(time.range[1], min(unlist(times[[i]], recursive = TRUE), na.rm = TRUE), na.rm = TRUE), max(time.range[2], max(unlist(times[[i]], recursive = TRUE), na.rm = TRUE), na.rm = TRUE))
    viability.range <- c(0, max(viability.range[2], max(unlist(responses[[i]], recursive = TRUE), na.rm = TRUE), na.rm = TRUE))
  }
  x1 <- 24; x2 <- 0

  ## FINDS INTERSECTION OF RANGES IF MORE THAN ONE tSet PLOTTED
  if (length(times) > 1) {
    common.ranges <- .getCommonConcentrationRange(times)

    for (i in seq_along(times)) {
      x1 <- min(x1, min(common.ranges[[i]]))
      x2 <- max(x2, max(common.ranges[[i]]))
    }
  }

  # SETS CUSTOM RANGE FOR X-AXIS IF PASSED AS ARGUEMENT
  if (!is.null(xlim)) {
    time.range <- xlim
  }

  ## SETS CUSTOM RANGE FOR Y-AXIS IF PASSED AS ARGUEMENT
  if (!is.null(ylim)) {
    print("ylim not null")
    viability.range <- ylim
  }

  ## SETS PLOT TITLE
  if(missing(title)){
    if(!missing(drug) && !missing(cell.lines)){
      title <- sprintf("%s:%s", drug, cell.lines)
    } else {
      title <- "Drug Time Response Curve"
    }
  }

  ## SETS DEFAULT COLOUR PALETTE
  if (missing(mycol)) {
    mycol <- RColorBrewer::brewer.pal(n = 9, name = "Set1")
    legends.col <- mycol
  }

  #### DRAWING THE PLOT ####
  plot(NA, xlab = "Time (hr)", ylab = "% Viability", axes = FALSE, main = title, ylim = viability.range, xlim = time.range, cex = cex, cex.main = cex.main)
  # Adds plot axes
  if (!is.null(x.custom.ticks)) {
    magicaxis::magaxis(side = 1:2, frame.plot = TRUE, tcl = -.3, majorn = c(x.custom.ticks[1], 3), minorn = c(x.custom.ticks[2], 2))
  } else {
    magicaxis::magaxis(side = 2, frame.plot = TRUE, tcl = -.3, majorn = c(3), minorn = c(2))
    axis(1, labels = as.numeric(duration), at = as.numeric(duration))
  }

  # Initialize legends variables
  legends <- NULL
  pch.val <- NULL
  legends.col <- NULL
  # TBD what this dose
  if (length(times) > 1) {
    rect(xleft = x1, xright = x2, ybottom = viability.range[1] , ytop = viability.range[2] , col = rgb(240, 240, 240, maxColorValue = 255), border = FALSE)
  }
  if (summarize.replicates == FALSE) {
    # Loop over tSet
    for (i in seq_along(times)) {
      j <- 1
      # Loop over dose level
      for (level in seq_along(dose)) {
        # Loop over replicates per dose level
        ## TODO:: Generalize this for n replicates
        for (replicate in seq_along(unique(responses[[i]][[level]]))) {
          # Plot per tSet, per dose level, per replicate points
          points(times[[i]][[replicate]], responses[[i]][[level]][[replicate]], pch = replicate, col = mycol[j], cex = cex)
          lines(times[[i]][[replicate]], responses[[i]][[level]][[replicate]], lty = replicate, lwd = lwd, col = mycol[j])

          legends <- c(legends, legendValues[[i]][[level]][[replicate]])
          legends.col <- c(legends.col, mycol[j])
          pch.val <- c(pch.val, replicate)
        }
        j <- j + 1
      }
    }
  } else {
    # Loop over tSet
    for (i in seq_along(times)) {
      j <- 1
      # Loop over dose level
      for (level in seq_along(dose)) {

        # Plot per tSet, per dose level, per replicate points
        points(times[[i]], responses[[i]][[level]], pch = 1, col = mycol[j], cex = cex + 0.2)
        lines(times[[i]], responses[[i]][[level]], lty = 1, lwd = lwd, col = mycol[j])

        legends <- c(legends, legendValues[[i]][[level]])
        legends.col <- c(legends.col, mycol[j])
        pch.val <- c(pch.val, 1)
        j <- j + 1
      }
    }
  }
  legend(legend.loc, legend = legends, col = legends.col, bty = "n", inset = c(0,2), xpd = TRUE, cex = cex, pch = pch.val)
}
