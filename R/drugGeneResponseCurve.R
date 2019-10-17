#' Compares gene expression for a specificed set of features over specific
#'   drug dosages vs time
#'
#' This function generates a plot visualizing the relationship between gene
#'   expression, time and dose level for the selected tSet.
#'
#' @examples
#'
#' if (interactive()) {
#' drugGeneResponseCurve(TGGATESsmall, dose = c("Control", "Low", "Middle"),
#'   mDataTypes="rna", drug = "naphthyl isothiocyanate",
#'   duration = c("2", "8", "24"), features = "ENSG00000000003_at")
#' }
#'
#' @param tSets [ToxicoSet] A ToxicoSet to be plotted in this graph. Currently
#'   only a single tSet is supported, passing more may results in errors.
#' @param dose [character] A vector of dose levels to be included in the
#'   plot. Default to include all dose levels available for a drug. Must include
#'   at minimum two dose levels, one of which must be "Control".
#' @param mDataTypes [vector] A vector specifying the molecular data types to
#'   include in this plot. Defaults to the first mDataType if not specified.
#'   This release version only accepts one mDataType, more to be added in
#'   forthcoming releases.
#' @param features [character] A vector of feature names to include in the plot.
#'   Please note that using too many features will have a significant computational
#'   cost and will likely result in a over crowded plot.
#' @param drug [character] A vector of drugs to be included in this plot. In
#'   this release, only one drug is supported.
#' @param duration [character] A vector of durations to include in the plot.
#' @param cellline [character] A vector of cell lines to include in the plot.
#' @param xlim [numeric] A vector of minimum and maximum values for the x-axis
#'   of the returned plot.
#' @param ylim [numeric] A vector of minimum and miximum values for the y-axis
#'   of the returned plot.
#' @param title [character] A string containing the desired plot name. If excluded
#'   a title wil be generated automatically.
#' @param legend.loc [character] The location of the legend as passed to the plot()
#'   function from base graphics. Suggested values are either "topleft" or
#'   "topright", with the latter as the default.
#' @param mycol [vector] A vector of length equal to the number of features
#'   argument specifying which RColorBrewer colour to use per feature. Default
#'   colours will be used if this parameter is excluded.
#' @param plot.type [character] The type of plot which you would like returned. Options
#'   are 'Actual' for unfitted curve, 'Fitted' for the fitted curve and 'Both'
#'   to display 'Actual and 'Fitted' in the sample plot. Currently this function
#'   only supports 'Actual'.
#' @param summarize.replicates [logical] If true will take the average of all
#'  replicates at each time point per gene and duration. This release has not
#'  yet implemented this feature.
#' @param x.custom.ticks [vector] A numeric vector of the distance between major
#'   and minor ticks on the x-axis. If excluded ticks appear only where duration
#'   values are specified.
#' @param lwd [numeric] The line width to plot width
#' @param cex [numeric] The cex parameter passed to plot. Controls the size of
#'   plot points and the font size of the legend and defaults to 0.7.
#' @param cex.main [numeric] The cex.main parameter passed to plot,
#'   controls the size of the titles and defaults to 1.0.
# @param trunc [bool] Should the viability values be truncated to lie in
#   [0-100] before doing the fitting
#' @param verbose [boolean] Should warning messages about the data passed
#'   in be printed?
#'
#' @return Plot of the viabilities for each drug vs time of exposure
#'
#' @import RColorBrewer
#' @importFrom graphics plot rect points lines legend
#' @importFrom grDevices rgb
#' @importFrom magicaxis magaxis
#'
#' @export
drugGeneResponseCurve <- function(
  tSets,
  duration,
  cellline,
  mDataTypes,
  features = NULL,
  dose,
  drug,
  plot.type="Actual",
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

  # Place tSets in a list if not already
  if (!is(tSets, "list")) { tSets <- list(tSets) }

  if (length(tSets) > 1) { warning("Multiple tSet plotting has not been tested in this release...")}
  if (length(drug) > 1) { stop("This function currently only supports one drug per plot...")}
  if (length(mDataTypes) > 1) {stop("This function currently only supports one molecular data type per plot...")}
  if (length(features) > 1) {warning("Plot scaling for multiple features has not yet been implemented for this release...")}

  ## TODO:: Generalize this to work with multiple data types
  if (missing(mDataTypes)) { mDataTypes <- names(tSets[[1]]@molecularProfiles) }

  if (is.null(features)) {
    features <- lapply(tSets, function(tSet) {
      rownames(featureInfo(tSet, "rna"))[seq_len(5)]
    })
  }

  # Places features in list if not already
  if (!is(features, "list")) { features <- list(features) }

  # Subsetting the tSets based on parameter arguments
  tSets <- lapply(tSets, function(tSet) {
    #if (is.null(features)) { features <- lapply(mDataTypes, function(mDataType) { featureInfo(tSet, mDataType)$gene_id }) }
    ToxicoGx::subsetTo(tSet, mDataType = mDataTypes, drugs = drug,
                       duration = duration, features = unique(unlist(features)))
  })

  # Extracting the data required for plotting into a list of data.frames
  # list of tSets < list of mDataTypes <df of plotData
  plotData <- lapply(tSets, function(tSet) {
    mDataTypesData <- lapply(mDataTypes, function(mDataType) {
      profileMatrix <- molecularProfiles(tSet, mDataType) # Sensitivity
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
  names(plotData) <- vapply(tSets, names, FUN.VALUE = character(1))

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
  names(times) <- vapply(tSets, function(x) names(x), FUN.VALUE = character(1)) # Get the names for each tSet

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
  names(legendValues) <-  vapply(tSets, function(x) names(x), FUN.VALUE = character(1))

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
  names(expression) <-  vapply(tSets, function(x) names(x), FUN.VALUE = character(1))

  #### SUMMARIZATION ####

  # Summarizing replicate values
  #if (summarize.replicates == TRUE) {
  #  expression <- lapply(seq_along(tSets), function(t_idx) {
  #    lapply(seq_along(dose), function(d_idx) {
  #      expressionVals <- NULL
  #      responseVect <- unlist(expression[[t_idx]][[d_idx]])
  #      for (time in seq_along(unique(duration))) {
  #        expressionVals <- c(expressionVals, mean(responseVect[time], responseVect[time + length(unique(duration))]))
  #      }
  #      expressionVals
  #    })
  #  })
  #  # Take unique values of all time replicates and place into a list
  #  times <- list(unique(unlist(times)))
  #  legendValues <- lapply(legendValues, function(legendLevel) {
  #    lapply(legendLevel, function(legendName) {unique(gsub("_[^_]*$", "", unlist(legendName)))})
  #  })
  #
  #}

  #### AXIS RANGES ####

  # Set x and y axis ranges based on time and viability values
  time.range <- as.numeric(c(min(unlist(unlist(times))), max(unlist(unlist(times)))))
  expression.range <- c(floor(min(unlist(expression, recursive = TRUE))), ceiling(max(unlist(expression, recursive=TRUE))))
  for (i in seq_along(tSets)) {
    ## TODO:: Generalize this to n replicates
    time.range <- c(min(time.range[1], min(unlist(times[[i]], recursive = TRUE), na.rm = TRUE), na.rm = TRUE), max(time.range[2], max(unlist(times[[i]], recursive = TRUE), na.rm = TRUE), na.rm = TRUE))
    expression.range <- c(0, max(expression.range[2], max(unlist(expression[[i]], recursive=TRUE), na.rm = TRUE), na.rm = TRUE))
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
  if(missing(title)){
    if(!missing(drug) && !missing(cellline)){
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
    magicaxis::magaxis(side = 1:2, frame.plot = TRUE, tcl= -.3, majorn= c(x.custom.ticks[1], 3), minorn = c(x.custom.ticks[2], 2))
  } else {
    magicaxis::magaxis(2, frame.plot = TRUE, tcl = -.3, majorn=expression.range[2], minorn=2)
    #axis(2, labels = seq_len(expression.range[2] + 1), at = seq_len(expression.range[2] + 1))
    axis(1, labels = as.numeric(duration), at = as.numeric(duration))
  }

  # Initialize legends variables
  legends <- NULL
  pch.val <- NULL
  legends.col <- NULL
  # TBD what this dose
  if (length(times) > 1) {
    rect(xleft = x1, xright = x2, ybottom = expression.range[1] , ytop = expression.range[2] , col=rgb(240, 240, 240, maxColorValue = 255), border = FALSE)
  }
  if (summarize.replicates == FALSE) {
    # Loop over tSets
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
              switch(plot.type,
                     "Actual" = {
                       lines(times[[i]][[mDataType]][[level]], expression[[i]][[mDataType]][[level]][[replicate]][[feature]], lty = replicate, lwd = lwd, col = mycol[j])
                     },
                     "Fitted" = {
                       #log_logistic_params <- logLogisticRegression(conc=times[[i]][[replicate]], expression = expression[[i]][[level]][[replicate]])
                       #x_vals <- .GetSupportVec(times[[i]][[replicate]])
                       #lines(10 ^ x_vals, .Hill(x_vals, pars=c(log_logistic_params$HS, log_logistic_params$E_inf/100, log10(log_logistic_params$EC50))) * 100 ,lty = replicate, lwd = lwd, col = mycol[j])
                     },
                     "Both" = {
                       lines(times[[i]][[mDataType]][[level]], expression[[i]][[mDataType]][[level]][[replicate]][[feature]], lty = replicate, lwd = lwd, col = mycol[j])
                       #log_logistic_params <- logLogisticRegression(conc = times[[i]][[replicate]], expression = expression[[i]][[level]][[replicate]])
                       #x_vals <- .GetSupportVec(times[[i]][[replicate]])
                       #lines(10 ^ x_vals, .Hill(x_vals, pars = c(log_logistic_params$HS, log_logistic_params$E_inf/100, log10(log_logistic_params$EC50))) * 100, lty=replicate, lwd=lwd, col=mycol[j])
                     })
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
    # Loop over tSets
    stop("Summarized replicates has net been implemented in this package yet.")
    for (i in seq_along(times)) {
      j <- 1
      for (mDataType in seq_along(mDataTypes)) {
        # Loop over dose level
        for (level in seq_along(dose)) {
          # Loop over replicates per dose level
          ## TODO:: Generalize this for n replicate
          # Plot per tSet, per dose level, per replicate points
          points(times[[i]][[mDataType]][[level]], expression[[i]][[mDataType]][[level]], pch = 1, col = mycol[j], cex = cex + 0.2)
          # Select plot type
          switch(plot.type,
                 "Actual" = {
                   lines(times[[i]][[mDataType]], expression[[i]][[level]], lty = 1, lwd = lwd, col = mycol[j])
                 },
                 "Fitted" = {
                   #log_logistic_params <- logLogisticRegression(conc = times[[i]], expression=expression[[i]][[level]])
                   stop("Curve fitting has not been implemented in this function yet. Feature coming soon!")
                   #x_vals <- .GetSupportVec(times[[i]])
                   #lines(10 ^ x_vals, .Hill(x_vals, pars = c(log_logistic_params$HS, log_logistic_params$E_inf/100, log10(log_logistic_params$EC50))) * 100 ,lty=1, lwd=lwd, col=mycol[j])
                 },
                 "Both" = {
                   warning("Curve fitting has not been implemented in this function yet. Feature coming soon!")
                   lines(times[[i]],expression[[i]][[mDataType]][[level]],lty=1, lwd = lwd, col = mycol[j])
                   #log_logistic_params <- logLogisticRegression(conc = times[[i]], expression = expression[[i]][[level]])
                   #x_vals <- .GetSupportVec(times[[i]])
                   #lines(10 ^ x_vals, .Hill(x_vals, pars = c(log_logistic_params$HS, log_logistic_params$E_inf/100, log10(log_logistic_params$EC50))) * 100, lty = 1, lwd=lwd, col=mycol[j])
                 })
          legends <- c(legends, legendValues[[i]][[level]])
          legends.col <- c(legends.col, mycol[j])
          pch.val <- c(pch.val, 1)
          j <- j + 1
        }
      }
    }
  }
  legend(legend.loc, legend = legends, col = legends.col, bty="n", cex = cex, pch = pch.val)
}
