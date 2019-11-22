#' Compares viabilities at a given dose over different experimental duration
#'
#' Description of this function
#'
#' @examples
#' if (interactive()) {
#'   drugTimeResponseCurve(TGGATESsmall, cell.lines = "Hepatocyte", dose = c("Control", "Low", "Middle"), drugs = drugNames(TGGATESsmall)[6], duration = c("2", "8", "24"))
#' }
#'
#' @param tSet \code{ToxicoSet} A ToxicoSet to be plotted in
#'   this figure
#' @param dose \code{character} A vector of dose levels to be included in the
#'   plot. Default to include all dose levels available for a drug. Must include
#'   at minimum two dose levels, one of witch is "Control".
#' @param drugs \code{character} A drugs or pair of drugss to be plotted.
#' @param duration \code{character} A vector of durations to include in the plot.
#' @param cell.lines \code{character} A vector of cell lines to include in the plot.
#' @param xlim \code{numeric} A vector of minimum and maximum values for the x-axis
#'   of the returned plot.
#' @param ylim \code{numeric} A vector of minimum and miximum values for the y-axis
#'   of the returned plot.
#' @param title \code{character} A string containing the desired plot name. If excluded
#'   a title wil be generated automatically.
#' @param mycol `vector`` A vector of length equal to the product of the
#'   number of drugs, features and doses passed to the function. Takes colour
#'   arguments as passed to `col` parameter in the `plot()` function.
#'   Default palette is used when unspecified.
#' @param summarize.replicates \code{logical} If true will take the average of all
#'  replicates at each time point per gene and duration. This release has not
#'  yet implemented this feature.
#' @param x.custom.ticks \code{vector} A numeric vector of the distance between major
#'   and minor ticks on the x-axis. If excluded ticks appear only where duration
#'   values are specified.
#' @param legend.loc \code{character} A string specifying the legend location, passed
#    to the base R \code{legend} function. Defaults to "bottomright".
#' @param lwd \code{numeric} The line width to plot width
#' @param cex \code{numeric} The cex parameter passed to plot
#' @param cex.main \code{numeric} The cex.main parameter passed to plot, controls
#' the size of the titles
# @param trunc \code{bool} Should the viability values be truncated to lie in \code{0-100}
# before doing the fitting
#' @param verbose \code{boolean} Should warning messages about the data passed in be printed?
#'
#' @return Plot of the viabilities for each drugs vs time of exposure
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
  drugs = NULL,
  summarize.replicates = TRUE,
  xlim = NULL,
  ylim = NULL,
  mycol,
  x.custom.ticks = NULL,
  title,
  lwd = 1.5,
  cex = 1,
  cex.main = 1.0,
  legend.loc = "bottomright",
  verbose=TRUE
) {
  # Place tSet in a list if not already
  if (!is(tSet, "list")) {
    tSet <- list(tSet)
  }

  ## TODO:: REMOVE ::: BEFORE CRAN SUBMISSION!
  ToxicoGx:::paramErrorChecker("drugTimeResponseCurve",
                    tSets = tSet, drugs = drugs, duration = duration,
                    cell.lines = cell.lines, dose = dose)

  ## TODO:: Make this function work with multiple drugs
  ## TODO:: Throw warning if a dose level or time point is not available for a specific drug
  ## TODO:: Add logic to handle viability_as_pct = FALSE?

  # Subsetting the tSet based on parameter arguments
  tSet <- lapply(tSet, function(tSet) {
    suppressWarnings({subsetTo(tSet, mDataType = "rna", drugs = drugs, duration = duration, cells = cell.lines)})
  })

  # Extracting the data required for plotting into a list of data.frames
  plotData <- lapply(tSet, function(tSet) {
    c <- lapply(cell.lines, function(cell) {
      d <- lapply(drugs, function(drg) {
        cbind(
          tSet@sensitivity$raw[,,2][
            which(rownames(sensitivityInfo(tSet)) %in%
                    rownames(subset(sensitivityInfo(tSet), cellid == cell & drugid == drg))),
            which(dose %in% c("Control", "Low", "Middle", "High"))],
          sensitivityInfo(tSet)[
            which(rownames(sensitivityInfo(tSet)) %in%
                    rownames(subset(sensitivityInfo(tSet), cellid == cell & drugid == drg))),
            c("duration_h", "replicate")
            ])
      })
      names(d) <- drugs; d
    })
    names(c) <- cell.lines; c
  })
  names(plotData) <- vapply(tSet, names, FUN.VALUE = character(1))

  # Extractiing the duration values for each replicate of plot
  times <- lapply(plotData, function(data) {
    c <- lapply(cell.lines, function(cell) {
      d <- lapply(drugs, function(drg) {
        lapply(unique(data[[cell]][[drg]]$replicate), function(rep) {
          as.numeric(subset(data[[cell]][[drg]], replicate == rep)$duration_h)
        })
      })
      names(d) <- drugs; d
    })
    names(c) <- cell.lines; c
  })
  names(times) <- vapply(tSet, names, FUN.VALUE = character(1))

  # Assembling the legend names for each line to be plotted
  legendValues <- lapply(plotData, function(data){
    c <- lapply(cell.lines, function(cell){
      d <- lapply(drugs, function(drg){
        ds <- lapply(dose, function(level) {
          lapply(unique(data[[cell]][[drg]]$replicate), function(rep) {
            paste(cell, drg, level, rep, sep = '_')
          })
        })
        names(ds) <- dose; ds
      })
      names(d) <- drugs; d
    })
    names(c) <- cell.lines; c
  })
  names(legendValues) <- vapply(tSet, names, FUN.VALUE = character(1))

  # Extracting the viability values for each row of plotData
  responses <- lapply(plotData, function(data) {
    c <- lapply(cell.lines, function(cell){
      d <- lapply(drugs, function(drg){
        ds <- lapply(dose, function(level) {
          lapply(unique(data[[cell]][[drg]]$replicate), function(rep) {
            as.vector(data[[cell]][[drg]][ , which(c("Control", "Low", "Middle", "High") %in% level)])[which(data[[cell]][[drg]]$replicate == rep)]
          })
        })
        names(ds) <- dose; ds
      })
      names(d) <- drugs; d
    })
    names(c) <- cell.lines; c
  })
  names(responses) <- vapply(tSet, names, FUN.VALUE = character(1))

  # Summarizing replicate values
  if (summarize.replicates == TRUE) {
    ## TODO:: Make this work in one lapply statement for all values
    responses <- lapply(responses, function(data) {
      c <- lapply(cell.lines, function(cell) {
        d <- lapply(drugs, function(drg) {
          ds <- lapply(dose, function(level) {
            ## TODO:: Generalize to n replicates
            v <- vapply(seq_along(unique(duration)), function(t) {
              mean(c(data[[cell]][[drg]][[level]][[1]][t],
                   data[[cell]][[drg]][[level]][[2]][t]))
            }, FUN.VALUE = numeric(1))
            v[!is.na(v)] # Remove NA values if mesaurements are missing for a t
          })
          names(ds) <- dose; ds
        })
        names(d) <- drugs; d
      })
      names(c) <- cell.lines; c
    })
    names(responses) <- vapply(tSet, names, FUN.VALUE = character(1))

    # Match time values to summarized replicates
    times <- lapply(times, function(data){
      c <- lapply(cell.lines, function(cell) {
        d <- lapply(drugs, function(drg){
          unique(unlist(data[[cell]][[drg]]))
        })
        names(d) <- drugs; d
      })
      names(c) <- cell.lines; c
    })
    names(times) <- vapply(tSet, names, FUN.VALUE = character(1))

    # Truncate replicate from the legend labels
    legendValues <- lapply(legendValues, function(data){
      c <- lapply(cell.lines, function(cell) {
        d <- lapply(drugs, function(drg) {
          ds <- lapply(dose, function(level) {
            unique(gsub("_[^_]*$", "", unlist(data[[cell]][[drg]][[level]])))
          })
          names(ds) <- dose; ds
        })
        names(d) <- drugs; d
      })
      names(c) <- cell.lines; c
    })
    names(legendValues) <- vapply(tSet, names, FUN.VALUE = character(1))
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

  # SETS CUSTOM RANGE FOR X-AXIS IF PASSED AS ARGUMENT
  if (!is.null(xlim)) {
    time.range <- xlim
  }

  ## SETS CUSTOM RANGE FOR Y-AXIS IF PASSED AS ARGUMENT
  if (!is.null(ylim)) {
    viability.range <- ylim
  }

  ## SETS PLOT TITLE
  if (missing(title)) {
    if (!missing(drugs) && !missing(cell.lines)){
      title <- sprintf("%s\n%s:%s", "Drug Time Response Curve", paste(drugs, collapse = " & "), paste(cell.lines, collapse = " & "))
    } else {
      title <- "Drug Time Response Curve"
    }
  }

  ## SETS DEFAULT COLOUR PALETTE
  if (missing(mycol)) {
    mycol <- c(RColorBrewer::brewer.pal(n = 9, name = "Set1"), RColorBrewer::brewer.pal(n = 12, name = "Set3"))
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
  if (!summarize.replicates) {
    # Loop over tSet
    for (data in seq_along(tSet)) {
      for (cell in cell.lines) {
        if (length(drugs) > 1) {j <- 1}
        for (drg in drugs) {
            if (length(drugs) == 1) {j <- 1}
            # Loop over dose level
            for (level in dose) {
              # Loop over replicates per dose level
              ## TODO:: Generalize this for n replicates
              for (replicate in seq_along(responses[[data]][[cell]][[drg]][[level]])) {
                # Plot per tSet, per dose level, per replicate points
                points(times[[data]][[cell]][[drg]][[replicate]],
                       responses[[data]][[cell]][[drg]][[level]][[replicate]],
                       pch = j, col = mycol[j], cex = cex)
                lines(times[[data]][[cell]][[drg]][[replicate]],
                      responses[[data]][[cell]][[drg]][[level]][[replicate]],
                      lty = replicate, lwd = lwd, col = mycol[j])
                legends <- c(legends, legendValues[[data]][[cell]][[drg]][[level]][[replicate]])
                legends.col <- c(legends.col, mycol[j])
                pch.val <- c(pch.val, j)
              }
              j <- j + 1
            }
        j <- j + 1
        }
      }
    }
  } else {
    # Loop over tSet
    for (data in seq_along(tSet)) {
      for (cell in cell.lines) {
        if (length(drugs) > 1) { j <- 1 }
        for (drg in drugs) {
          if (length(drugs) == 1 ) { j <- 1 }
          # Loop over dose level
          for (level in dose) {
            points(times[[data]][[cell]][[drg]],
                   responses[[data]][[cell]][[drg]][[level]],
                   pch = j, col = mycol[j], cex = cex)
            lines(times[[data]][[cell]][[drg]],
                  responses[[data]][[cell]][[drg]][[level]],
                  lty = which(cell.lines %in% cell), lwd = lwd, col = mycol[j])
            legends <- c(legends, legendValues[[data]][[cell]][[drg]][[level]])
            legends.col <- c(legends.col, mycol[j])
            pch.val <- c(pch.val, which(dose %in% level))
            j <- j + 1
          }
          j <- j + 1
        }
      }
    }
  }
  legend(legend.loc, legend = legends, col = legends.col, bty = "n", cex = cex, pch = pch.val)
}

