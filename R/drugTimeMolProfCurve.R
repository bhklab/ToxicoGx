#' Compares gene expression for a specificed set of features over specific
#'   drug dosages vs time
#'
#' This function generates a plot visualizing the relationship between gene
#'   expression, time and dose level for the selected tSet.
#'
#' @examples
#'
#' if (interactive()) {
#' drugTimeMolProfCurve(TGGATESsmall, dose = c("Control", "Low", "Middle"), mDataTypes="rna", drug = drugNames(TGGATESsmall)[1], duration = c("2", "8", "24"), features = "ENSG00000000003_at")
#' }
#'
#' @param tSet \code{ToxicoSet} A ToxicoSet to be plotted in this graph. Currently
#'   only a single tSet is supported, passing more may results in errors.
#' @param dose \code{character} A vector of dose levels to be included in the
#'   plot. Default to include all dose levels available for a drug. Must include
#'   at minimum two dose levels, one of which must be "Control".
#' @param mDataTypes \code{vector} A vector specifying the molecular data types to
#'   include in this plot. Defaults to the first mDataType if not specified.
#'   This release version only accepts one mDataType, more to be added in
#'   forthcoming releases.
#' @param features \code{character} A vector of feature names to include in the plot.
#'   Please note that using too many features will have a significant computational
#'   cost and will likely result in a over crowded plot.
#' @param drug \code{character} A vector of drugs to be included in this plot. In
#'   this release, only one drug is supported.
#' @param duration \code{character} A vector of durations to include in the plot.
#' @param cellline \code{character} A vector of cell lines to include in the plot.
#' @param xlim \code{numeric} A vector of minimum and maximum values for the x-axis
#'   of the returned plot.
#' @param ylim \code{numeric} A vector of minimum and miximum values for the y-axis
#'   of the returned plot.
#' @param title \code{character} A string containing the desired plot name. If excluded
#'   a title wil be generated automatically.
#' @param legend.loc \code{character} The location of the legend as passed to the plot()
#'   function from base graphics. Suggested values are either "topleft" or
#'   "topright", with the latter as the default.
#' @param mycol \code{vector} A vector of length equal to the number of features
#'   argument specifying which RColorBrewer colour to use per feature. Default
#'   colours will be used if this parameter is excluded.
#' @param summarize.replicates \code{logical} If true will take the average of all
#'  replicates at each time point per gene and duration. This release has not
#'  yet implemented this feature.
#' @param x.custom.ticks \code{vector} A numeric vector of the distance between major
#'   and minor ticks on the x-axis. If excluded ticks appear only where duration
#'   values are specified.
#' @param lwd \code{numeric} The line width to plot width
#' @param cex \code{numeric} The cex parameter passed to plot. Controls the size of
#'   plot points and the font size of the legend and defaults to 0.7.
#' @param cex.main \code{numeric} The cex.main parameter passed to plot,
#'   controls the size of the titles and defaults to 1.0.
# @param trunc \code{bool} Should the viability values be truncated to lie in
#   \code{0-100} before doing the fitting
#' @param verbose \code{boolean} Should warning messages about the data passed
#'   in be printed?
#'
#' @return Plot of the viabilities for each drug vs time of exposure
#'
#' @import RColorBrewer
#' @importFrom graphics plot rect points lines legend
#' @importFrom grDevices rgb
#' @importFrom magicaxis magaxis
#' @importFrom foreach foreach
#'
#' @export
drugTimeMolProfCurve <- function(
  tSet,
  duration,
  cellline,
  mDataTypes,
  features = NULL,
  dose,
  drug,
  summarize.replicates = FALSE,
  xlim=c(0, 24),
  ylim=c(0, 15),
  mycol,
  x.custom.ticks = NULL,
  title,
  lwd = 1.5,
  cex = 0.7,
  cex.main = 1.0,
  legend.loc = "topright",
  verbose=TRUE
) {

  # Place tSet in a list if not already
  if (!is(tSet, "list")) { tSet <- list(tSet) }

  ## Tempary warnings until function is finished
  if (length(tSet) > 1) { warning("Multiple tSet plotting has not been tested in this release...")}
  if (length(drug) > 1) { stop("This function currently only supports one drug per plot...")}
  if (length(mDataTypes) > 1) {stop("This function currently only supports one molecular data type per plot...")}
  if (length(features) > 1) {warning("Plot scaling for multiple features has not yet been implemented for this release...")}

  ## TODO:: Generalize this to work with multiple data types
  if (missing(mDataTypes)) { mDataTypes <- names(tSet[[1]]@molecularProfiles) }

  if (is.null(features)) {
    features <- lapply(tSet, function(tSet) {
      rownames(featureInfo(tSet, "rna"))[seq_len(5)]
    })
  }

  # Places features in list if not already
  if (!is(features, "list")) { features <- list(features) }

  # Subsetting the tSet based on parameter arguments
  tSet <- lapply(tSet, function(tSet) {
    ToxicoGx::subsetTo(tSet, mDataType = mDataTypes, drugs = drug,
                       duration = duration, features = unique(unlist(features)))
  })

  # Extracting the data required for plotting into a list of data.frames
  # list of tSet < list of mDataTypes <df of plotData
  plotData <- lapply(tSet, function(tSet) {
    mDataTypesData <- lapply(mDataTypes, function(mDataType) {
      profileMatrix <- molecularProfiles(tSet, mDataType)
      relevantFeatureInfo <- featureInfo(tSet, mDataType)[, c("gene_id", "transcript_name") ]
      relevantPhenoInfo <- phenoInfo(tSet, mDataType)[, c("samplename", "cellid", "drugid", "concentration", "dose_level", "duration", "species", "individual_id")]
      relevantSensitivityInfo <- ToxicoGx::sensitivityInfo(tSet)[, c("drugid", "duration_h", "replicate", "Control", "Low", "Middle", "High") ]
      data <- list(profileMatrix, relevantFeatureInfo, relevantPhenoInfo, relevantSensitivityInfo)
      names(data) <- c("data", "featureInfo", "phenoInfo", "sensitivityInfo")
      data
    })
    names(mDataTypesData) <- mDataTypes # Name list items for easy to understand subsetting
    mDataTypesData
  })
  names(plotData) <- vapply(tSet, names, FUN.VALUE = character(1))

  # Get a list of times per tSet per mDataType per dose level
  # This will also need to be per drug if we extend the function to multiple drugs
  times <- lapply(plotData, function(tSetData) {
    mDataTimes <- lapply(mDataTypes, function(mDataType) {
      doseLevels <- lapply(unique(tSetData[[mDataType]]$phenoInfo$dose_level), function(doseLvl){ ## TODO:: Fix this to only include listed doses
        as.numeric(unique(tSetData[[mDataType]][["phenoInfo"]][["duration"]][ which(tSetData[[mDataType]]$phenoInfo$dose_level %in% doseLvl)]))
      })
      names(doseLevels) <- dose; doseLevels
    })
    names(mDataTimes) <- mDataTypes; mDataTimes
  })
  names(times) <- vapply(tSet, function(x) names(x), FUN.VALUE = character(1)) # Get the names for each tSet

  # Assembling the legend names for each line to be plotted
  legendValues <- lapply(plotData, function(tSetData) {
    mDataLegends <- lapply(mDataTypes, function(mDataType) {
      doseLegends <- lapply(dose, function(doseLvl) {
        replicateLegends <- lapply(unique(tSetData[[mDataType]]$sensitivityInfo$replicate), function(rep){
          legendValues <- vapply(tSetData[[mDataType]]$featureInfo$gene_id, function(feature) {
            paste(
              doseLvl,
              paste(gsub("_at", "", tSetData[[mDataType]]$featureInfo[feature, "gene_id"])),
              paste(gsub("-.*", "", tSetData[[mDataType]]$featureInfo[feature, "transcript_name"])),
              rep,
              sep = "_" )
          }, FUN.VALUE = character(1))
          names(legendValues) <- unique(unlist(features)); legendValues
        })
      })
      names(doseLegends) <- dose; doseLegends
    })
    names(mDataLegends) <- mDataTypes; mDataLegends # Note: FUN.VALUE refers to the type and length of EACH function call, not of the returned vector
  })
  names(legendValues) <-  vapply(tSet, function(x) names(x), FUN.VALUE = character(1))

  # Expression
  expression <- lapply(plotData, function(tSetData) {
    mDataExpr <- lapply(mDataTypes, function(mDataType) {
      doseExpr <- lapply(dose, function(doseLvl) {
        replicateExpr <- lapply(unique(tSetData[[mDataType]]$sensitivityInfo$replicate), function(rep) {
          featureExpr <- lapply(tSetData[[mDataType]]$featureInfo$gene_id, function(feature) {
            sensInf <- tSetData[[mDataType]]$sensitivityInfo
            expressionVals <- tSetData[[mDataType]]$data[
              feature, # The feature of interest
              c(sensInf[ which(sensInf$replicate == rep) , doseLvl]) # Gets the sample ids for the given dose level and replicate
              ]
            #names(expressionVals) <- sensInf[which(sensInf$replicate == rep) , "duration_h"] # Name based on included durations
            #expressionVals
          })
          names(featureExpr) <- unique(unlist(features)); featureExpr
        })
      })
      names(doseExpr) <- dose; doseExpr
    })
    names(mDataExpr) <- mDataTypes; mDataExpr
  })
  names(expression) <-  vapply(tSet, function(x) names(x), FUN.VALUE = character(1))

  #### SUMMARIZATION ####

  # Summarizing replicate values
  if (summarize.replicates) {
   expression <- lapply(seq_along(expression), function(t_idx) {
     lapply(seq_along(mDataTypes), function(m_idx) {
       lapply(seq_along(dose), function(d_idx) {
         lapply(seq_along(features), function(f_idx) {
            vapply(seq_along(unique(duration)), function(idx) {
              ## TODO:: Generalize this to n replicates
                mean(expression[[t_idx]][[m_idx]][[d_idx]][[1]][[f_idx]][[idx]],
                     expression[[t_idx]][[m_idx]][[d_idx]][[1]][[f_idx]][[idx]]
                     )
            }, FUN.VALUE = numeric(1))
           })
         })
       })
     })
   # Take unique values of all time replicates and place into a list
   times <- lapply(times, function(tSet) {
     lapply(tSet, function(mDataType) {
      lapply(mDataType, function(dose) {
        unique(dose)
      })
     })
   })
   legendValues <- lapply(legendValues, function(legendLevel) {
      lapply(legendLevel, function(legendName) {unique(gsub("_[^_]*$", "", unlist(legendName)))})
   })
  }

  #### AXIS RANGES ####

  # Set x and y axis ranges based on time and viability values
  time.range <- as.numeric(c(min(unlist(times)), max(unlist(times))))
  expression.range <- c(floor(min(unlist(expression, recursive = TRUE))), ceiling(max(unlist(expression, recursive = TRUE))))
  for (i in seq_along(tSet)) {
    ## TODO:: Generalize this to n replicates
    time.range <- c(min(time.range[1], min(unlist(times[[i]], recursive = TRUE), na.rm = TRUE), na.rm = TRUE), max(time.range[2], max(unlist(times[[i]], recursive = TRUE), na.rm = TRUE), na.rm = TRUE))
    expression.range <- c(0, max(expression.range[2], max(unlist(expression[[i]], recursive = TRUE), na.rm = TRUE), na.rm = TRUE))
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
  if (!missing(xlim)) {
    time.range <- xlim
  }

  ## SETS CUSTOM RANGE FOR Y-AXIS IF PASSED AS ARGUEMENT
  if (!missing(ylim)) {
    expression.range <- ylim
  }

  ## SETS PLOT TITLE
  if (missing(title)) {
    if (!missing(drug) && !missing(cellline)){
      title <- sprintf("%s:%s", drug, cellline)
    } else {
      title <- "Expression Time Response Curve"
    }
  }

  ## SETS DEFAULT COLOUR PALETTE
  if (missing(mycol)) {
    mycol <- RColorBrewer::brewer.pal(n = 9, name = "Set1")
    legends.col <- mycol
  }

  #### DRAWING THE PLOT ####
  plot(NA, xlab = "Time (hr)", ylab = "Expression", axes = FALSE, main = title, ylim = expression.range, xlim = time.range, cex = cex, cex.main = cex.main)
  # Adds plot axes
  if (!is.null(x.custom.ticks)) {
    magicaxis::magaxis(side = 1:2, frame.plot = TRUE, tcl = -.3, majorn = c(x.custom.ticks[1], 3), minorn = c(x.custom.ticks[2], 2))
  } else {
    magicaxis::magaxis(2, frame.plot = TRUE, tcl = -.3, majorn = expression.range[2], minorn=2)
    #axis(2, labels = seq_len(expression.range[2] + 1), at = seq_len(expression.range[2] + 1))
    axis(1, labels = as.numeric(duration), at = as.numeric(duration))
  }

  # Initialize legends variables
  legends <- NULL
  pch.val <- NULL
  legends.col <- NULL
  # TBD what this does?
  if (length(times) > 1) {
    rect(xleft = x1, xright = x2, ybottom = expression.range[1] , ytop = expression.range[2] , col=rgb(240, 240, 240, maxColorValue = 255), border = FALSE)
  }

  if (summarize.replicates == FALSE) {
    # Loop over tSet
    for (i in seq_along(times)) { # tSet subset
      for (mDataType in seq_along(mDataTypes)) {
        # Loop over dose level
        j <- 1
        for (level in seq_along(dose)) {
          # Loop over replicates per dose level
          ## TODO:: Generalize this for n replicates
          for (replicate in seq_along(expression[[i]][[mDataType]][[level]])) {
            # Plot per tSet, per dose level, per replicate points
            for (feature in seq_along(features[[mDataType]])) {
              points(times[[i]][[mDataType]][[level]], expression[[i]][[mDataType]][[level]][[replicate]][[feature]], pch = replicate, col = mycol[j], cex = cex)
              # Select plot type
              lines(times[[i]][[mDataType]][[level]], expression[[i]][[mDataType]][[level]][[replicate]][[feature]], lty = replicate, lwd = lwd, col = mycol[j])
              legends <- c(legends, legendValues[[i]][[mDataType]][[level]][[replicate]][feature])
              legends.col <- c(legends.col, mycol[j])
              pch.val <- c(pch.val, feature)
            }
          }
          j <- j + 1
        }
      }
    }
  } else {
    # Loop over tSet
    for (i in seq_along(times)) { # tSet subset
      for (mDataType in seq_along(mDataTypes)) {
        # Loop over dose level
        j <- 1
        for (level in seq_along(dose)) {
            # Plot per tSet, per dose level, per replicate points
            for (feature in seq_along(features[[mDataType]])) {
              points(times[[i]][[mDataType]][[level]], expression[[i]][[mDataType]][[level]][[feature]], col = mycol[j], cex = cex)
              # Select plot type
              lines(times[[i]][[mDataType]][[level]], expression[[i]][[mDataType]][[level]][[feature]], lwd = lwd, col = mycol[j])
              legends <- c(legends, legendValues[[i]][[mDataType]][[level]][[feature]])
              legends.col <- c(legends.col, mycol[j])
              pch.val <- c(pch.val, feature)
            }
          j <- j + 1
        }
      }
    }
  }
  legend(legend.loc, legend = legends, col = legends.col, bty = "n", cex = cex, pch = pch.val)
}
