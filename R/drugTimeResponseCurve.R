#' Compares viabilities at a given dose over different time points
#'
#' Description of this function
#'
#' @examples
#' ToxicoGx:::drugTimeResponseCurve(TGGATESsmall, dose = c("Control", "Low", "Middle"), drug = "naphthyl isothiocyanate", duration = c("2", "8", "24))
#'
#' @param tSets [ToxicoSet] A ToxicoSet or list of ToxicoSets to be plotted in
#'   this graph.
#' @param dose [character] A vector of dose levels to be included in the
#'   plot. Default to include all dose levels available for a drug. Must include
#'   at minimum two dose levels, one of witch is "Control".
#' @param drug [character] A vector of drugs to be included in this plot.
#' @param duration [character] A vector of durations to include in the plot.
#' @param cellLine [character] A vector of cell lines to include in the plot.
#' @param viabilityAsPct [logical] A vector specifying if viabilities should
#'   be plotted as a percentage. Defaults to TRUE.
#' @param xlim [numeric] A vector of minimum and maximum values for the x-axis
#'   of the returned plot.
#' @param ylim [numeric] A vector of minimum and miximum values for the y-axis
#'   of the returned plot.
#' @param title [character] A string containing the desired plot name. If excluded
#'   a title wil be generated automatically.
#' @param legendLoc [character] The location of the legend as passed to the plot()
#'   function from base graphics.
#' @param mycol [vector] A vector of length equal to the lenth of the tSets
#'   argument specifying which RColorBrewer colour to use per tSet. Default
#'   colours will be used if this parameter is excluded.
#' @param plotType [character] The type of plot which you would like returned. Options
#'   are 'Actual' for unfitted curve, 'Fitted' for the fitted curve and 'Both'
#'   to display 'Actual and 'Fitted' in the sample plot.
#' @param lwd [numeric] The line width to plot with
#' @param cex [numeric] The cex parameter passed to plot
#' @param cex.main [numeric] The cex.main parameter passed to plot, controls the size of the titles
#' @param legend.loc And argument passable to xy.coords for the position to place the legend.
#' @param trunc [bool] Should the viability values be truncated to lie in [0-100] before doing the fitting
#' @param verbose [boolean] Should warning messages about the data passed in be printed?
#'
#' @return Plot of the viabilities for each drug vs time of exposure
#'
#' @import RColorBrewer
#' @importFrom graphics plot rect points lines legend
#' @importFrom grDevices rgb
#' @importFrom magicaxis magaxis
#'
#' @export
drugTimeResponseCurve <- function(
  tSets,
  duration,
  cellline,
  dose,
  drug,
  plotType="Actual",
  viabilityAsPct=TRUE,
  xlim=c(0, 24),
  ylim=c(0, 100),
  mycol,
  title,
  lwd = 0.5,
  cex = 0.5,
  cex.main = 0.9,
  legend.loc = "topright",
  verbose=TRUE
  ) {

  # Place tSets in a list if not already
  if (!is(tSets, "list")) {
    tSets <- list(tSets)
  }

  # paramErrorChecker()

  # Subsetting the tSets based on parameter arguments
  tSets <- lapply(tSets, function(tSet) {
    subsetTo(tSet, mDataType = "rna", drugs=drug, duration=duration)
  })

  # Extracting the data required for plotting into a list of data.frames
  plotData <- lapply(tSets, function(tSet) {
    cbind(
      tSet@sensitivity$raw[,,2][, which(c("Control", "Low", "Middle", "High") %in% dose)],
      sensitivityInfo(tSet)[, c("duration_h", "replicate")]
    )
  })

  # Getting the names for each
  legendValues <- lapply(plotData, function(data){
    lapply(seq_along(unique(data$replicate)), function(idx){
      rownames(data)[which(data$replicate == idx)]
    })
  })

  # Extractiing the duration values for each row of plotData
  times <- lapply(plotData, function(data) {
    lapply(seq_along(unique(data$replicate)), function(idx) {
      as.numeric(data$duration_h)[which(data$replicate == idx)]
    })
  })

  # Extracting the viability values for each row of plotData
  responses <- lapply(plotData, function (data) {
    responseVect <- lapply(dose, function(level) {
      lapply(seq_along(unique(data$replicate)), function(idx) {
        as.vector(data[ , which(c("Control", "Low", "Middle", "High") %in% level)])[which(data$replicate == idx)]
      })
    })
  })

  # Summarizing replicate values


  # Set x and y axis ranges based on time and viability values
  time.range <- c(min(unlist(unlist(times))), max(unlist(unlist(times))))
  viability.range <- c(min(unlist(responses, recursive = TRUE)), max(unlist(responses, recursive=TRUE)))
  for(i in seq_along(times)) {
    ## TODO:: Generalize this to n replicates
    time.range <- c(min(time.range[1], min(unlist(times[[i]], recursive=TRUE), na.rm=TRUE), na.rm=TRUE), max(time.range[2], max(unlist(times[[i]], recursive=TRUE), na.rm=TRUE), na.rm=TRUE))
    viability.range <- c(0, max(viability.range[2], max(unlist(responses[[i]], recursive=TRUE), na.rm=TRUE), na.rm=TRUE))
  }
  x1 <- 24; x2 <- 0

  ## FINDS INTERSECTION OF RANGES IF MORE THAN ONE tSet PLOTTED
  if(length(times) > 1) {
    common.ranges <- .getCommonConcentrationRange(times)

    for(i in seq_along(times)) {
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
    viability.range <- ylim
  }

  ## SETS PLOT TITLE
  if(missing(title)){
    if(!missing(drug) && !missing(cellline)){
      title <- sprintf("%s:%s", drug, cellline)
    } else {
      title <- "Drug Time Response Curve"
    }
  }

  ## SETS DEFAULT COLOUR PALETTE
  if (missing(mycol)) {
    mycol <- RColorBrewer::brewer.pal(n=9, name="Set1")
    legends.col <- mycol
  }

  print(viability.range)
  print(time.range)

  #### DRAWING THE PLOT ####
  plot(NA, xlab="Time (hr)", ylab="% Viability", axes =FALSE, main=title, ylim=viability.range, xlim=time.range, cex=cex, cex.main=cex.main)
  # Adds plot axes
  magicaxis::magaxis(side=1:2, frame.plot=TRUE, tcl=-.3, majorn=c(5,3), minorn=c(5,2))
  # Initialize legends variables
  legends <- NULL
  legends.col <- NULL
  # TBD what this dose
  if (length(times) > 1) {
    rect(xleft=x1, xright=x2, ybottom=viability.range[1] , ytop=viability.range[2] , col=rgb(240, 240, 240, maxColorValue = 255), border=FALSE)
  }
  # Loop over tSets
  for (i in seq_along(times)) {
    j <- 1
    # Loop over dose level
    for (level in seq_along(dose)) {
      # Loop over replicates per dose level
      ## TODO:: Generalize this for n replicates
      for(idx in seq_len(2)) {
        # Plot per tSet, per dose level, per replicate points
        points(times[[i]][[idx]], responses[[i]][[level]][[idx]], pch=20, col = mycol[j], cex=cex)
        # Select plot type
        switch(plotType,
             "Actual"={
               lines(times[[i]][[idx]], responses[[i]][[level]][[idx]], lty=idx, lwd=lwd, col=mycol[j])
             },
             "Fitted"={
               log_logistic_params <- logLogisticRegression(conc=times[[i]][[idx]], viability=responses[[i]][[level]][[idx]])
               x_vals <- .GetSupportVec(times[[i]][[idx]])
               lines(10 ^ x_vals, ToxicoGx:::.Hill(x_vals, pars=c(log_logistic_params$HS, log_logistic_params$E_inf/100, log10(log_logistic_params$EC50))) * 100 ,lty=idx, lwd=lwd, col=mycol[j])
             },
             "Both"={
               lines(times[[i]][[idx]],responses[[i]][[level]][[idx]],lty=idx, lwd=lwd, col = mycol[j])
               log_logistic_params <- logLogisticRegression(conc = times[[i]][[idx]], viability = responses[[i]][[level]][[idx]])
               x_vals <- .GetSupportVec(times[[i]][[idx]])
               lines(10 ^ x_vals, ToxicoGx:::.Hill(x_vals, pars=c(log_logistic_params$HS, log_logistic_params$E_inf/100, log10(log_logistic_params$EC50))) * 100 ,lty=idx, lwd=lwd, col=mycol[j])
             })
        legends<- c(legends, legendValues[[i]][[idx]][level])
        print(legends)
        legends.col <- c(legends.col, mycol[i])
      }
      j <- j + 1

    }
  }
  legend(legend.loc, legend=legends, col=legends.col, bty="n", cex=cex, pch=c(15,15))
}
