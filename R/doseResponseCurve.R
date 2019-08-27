#' Plot drug response curve of a given drug and a given cell for a list of rSets (objects of the RadioSet class).
#'
#' Given a list of RadioSets, the function will plot the drug_response curve,
#' for a given drug/cell pair. The y axis of the plot is the viability percentage
#' and x axis is the log transformed Ds. If more than one rSet is
#' provided, a light gray area would show the common concentration range between rSets.
#' User can ask for type of sensitivity measurment to be shown in the plot legend.
#' The user can also provide a list of their own Ds and viability values,
#' as in the examples below, and it will be treated as experiments equivalent to values coming
#' from a pset. The names of the concentration list determine the legend labels.
#'
#' @examples
#' doseResponseCurve(Ds=list("Experiment 1" = c(0, 2, 4, 6)),
#' SFs=list("Experiment 1" = c(1,.6,.4,.2)), plot.type="Both")
#'
#' @param rad.type [string] The type(s) of radiation dosage to be
#' plotted. If the plot is desirable for more than one radioset, A unique drug id
#' should be provided.
#' @param cellline [string] A cell line name for which the radiation response curve should be
#' plotted. If the plot is desirable for more than one radioset, a unique cell id
#' should be provided.
#' @param rSets [list] a list of RadioSet objects, for which the function
#' should plot the curves.
#' @param Ds,SFs [list] A list of Doses and SFs to plot, the function assumes that
#' Ds[[i]] is plotted against SFs[[i]]. The names of the D list are used to create the legend labels
#' @param legends.label [vector] A vector of sensitivity measurment types which could
#' be any combination of  ic50_published, auc_published, auc_recomputed and auc_recomputed_star.
#' A legend will be displayed on the top right of the plot which each line of the legend is
#' the values of requested sensitivity measerments for one of the requested rSets.
#' If this parameter is missed no legend would be provided for the plot.
#' @param ylim [vector] A vector of two numerical values to be used as ylim of the plot.
#' If this parameter would be missed c(0,100) would be used as the ylim of the plot.
#' @param xlim [vector] A vector of two numerical values to be used as xlim of the plot.
#' If this parameter would be missed the minimum and maximum comncentrations between all
#' the rSets would be used as plot xlim.
#' @param mycol [vector] A vector with the same lenght of the rSets parameter which
#' will determine the color of the curve for the pharmaco sets. If this parameter is
#' missed default colors from Rcolorbrewer package will be used as curves color.
#' @param plot.type [character] Plot type which can be the actual one ("Actual") or
#' the one fitted by logl logistic regression ("Fitted") or both of them ("Both").
#' If this parameter is missed by default actual curve is plotted.
#' @param summarize.replicates [character] If this parameter is set to true replicates
#' are summarized and replicates are plotted individually otherwise
#' @param title [character] The title of the graph. If no title is provided, then it defaults to
#' 'Drug':'Cell Line'.
#' @param lwd [numeric] The line width to plot with
#' @param cex [numeric] The cex parameter passed to plot
#' @param cex.main [numeric] The cex.main parameter passed to plot, controls the size of the titles
#' @param legend.loc And argument passable to xy.coords for the position to place the legend.
#' @param trunc [bool] Should the viability values be truncated to lie in [0-100] before doing the fitting
#' @param verbose [boolean] Should warning messages about the data passed in be printed?
#' @return Plots to the active graphics device and returns and invisible NULL.
#' @export
#' @import RColorBrewer
#' @importFrom graphics plot rect axis
#' @importFrom grDevices rgb
#' @importFrom graphics plot
#' @importFrom graphics rect
#' @importFrom grDevices rgb
#' @importFrom graphics points
#' @importFrom graphics lines
#' @importFrom graphics legend
#' @importFrom magicaxis magaxis



doseResponseCurve <-
function(rad.type = "radiation",
 cellline,
 rSets=list(),
 Ds=list(),
 SFs=list(),
         # conc_as_log = FALSE,
         # viability_as_pct = TRUE,
 trunc=TRUE,
 legends.label = c("alpha", "beta","rsquared"),
 ylim=c(0,100),
 xlim, mycol,
 title,
 plot.type=c("Fitted","Actual", "Both"),
 summarize.replicates=TRUE,
 lwd = 1,
 cex = 0.7,
 cex.main = 0.9,
 legend.loc = "topright",
 verbose=TRUE) {
  if(!missing(rSets)){
    if (class(rSets) != "list") {
      if (class(rSets) == "RadioSet") {
        temp <- rSetName(rSets)
        rSets <- list(rSets)
        names(rSets) <- temp
      } else {
        stop("Type of rSets parameter should be either a rSet or a list of rSets.")
      }
    }
  }
  if(!all(legends.label %in% c("alpha", "beta","rsquared"))){
    stop(paste("Only", paste(c("'alpha'", "'beta'","'rsquared'"), collapse = ", "), "implemented for legend labels.", split = " "))
  }
  # if(!missing(rSets) && (missing(rad.type) || missing(cellline))){
  ## XXX:: HACK
  if(!missing(rSets) && (missing(cellline))){
    stop("If you pass in a rSet then drug and cellline must be set") }
  # } else {
  #   if(missing(drug)){
  #   drug <- "Drug"}
  #   if(missing(cellline))
  #   cellline <- "Cell Line"
  # }
    if(!missing(Ds)){
      if(missing(SFs)){

        stop("Please pass in the Survival Fractions to Plot with the Doses.")

      }
      if (class(Ds) != "list") {
        if (mode(Ds) == "numeric") {
          if(mode(SFs)!="numeric"){
            stop("Passed in 1 vector of Doses but the Survival Fractions are not numeric!")
          }
        # cleanData <- sanitizeInput(Ds,
        #   SFs,
        #   conc_as_log = conc_as_log,
        #   viability_as_pct = viability_as_pct,
        #   trunc = trunc,
        #   verbose = verbose)
        # Ds <- 10^cleanData[["log_conc"]]
          Ds <- list(Ds)
        # SFs <- 100*cleanData[["viability"]]
          SFs <- list(SFs)
          names(Ds) <- "Exp1"
          names(SFs) <- "Exp1"
        } else {
          stop("Mode of Doses parameter should be either numeric or a list of numeric vectors")
        }
      } else{
        if(length(SFs)!= length(Ds)){
          stop("The number of D and SF vectors passed in differs")
        }
        if(is.null(names(Ds))){
          names(Ds) <- paste("Exp", 1:length(Ds))
        }
        for(i in 1:length(Ds)){

          if (mode(Ds[[i]]) == "numeric") {
            if(mode(SFs[[i]])!="numeric"){
              stop(sprintf("Ds[[%d]] are numeric but the SFs[[%d]] are not numeric!",i,i))
            }
          # cleanData <- sanitizeInput(Ds[[i]],
          #   SFs[[i]],
          #   conc_as_log = conc_as_log,
          #   viability_as_pct = viability_as_pct,
          #   trunc = trunc,
          #   verbose = verbose)
          # Ds[[i]] <- 10^cleanData[["log_conc"]]
          # SFs[[i]] <- 100*cleanData[["viability"]]
          } else {
            stop(sprintf("Mode of Ds[[%d]] parameter should be numeric",i))
          }

        }

      }
    }

    common.range.star <- FALSE

    if (missing(plot.type)) {
      plot.type <- "Actual"
    }

    doses <- list(); responses <- list(); legend.values <- list(); j <- 0; rSetNames <- list()
    if(!missing(rSets)){
      for(i in 1:length(rSets)) {
        exp_i <- which(sensitivityInfo(rSets[[i]])[ ,"cellid"] == cellline & sensitivityInfo(rSets[[i]])[ ,"radiation.type"] == rad.type)
        if(length(exp_i) > 0) {
          if (summarize.replicates) {
            rSetNames[[i]] <- rSetName(rSets[[i]])
            if (length(exp_i) == 1) {
              drug.responses <- as.data.frame(cbind("Dose"=as.numeric(as.vector(rSets[[i]]@sensitivity$raw[exp_i, , "Dose"])),
                "Viability"=as.numeric(as.vector(rSets[[i]]@sensitivity$raw[exp_i, , "Viability"])), stringsAsFactors=FALSE))
              drug.responses <- drug.responses[complete.cases(drug.responses), ]
            }else{
              drug.responses <- as.data.frame(cbind("Dose"=apply(rSets[[i]]@sensitivity$raw[exp_i, , "Dose"], 2, function(x){median(as.numeric(x), na.rm=TRUE)}),
                "Viability"=apply(rSets[[i]]@sensitivity$raw[exp_i, , "Viability"], 2, function(x){median(as.numeric(x), na.rm=TRUE)}), stringsAsFactors=FALSE))
              drug.responses <- drug.responses[complete.cases(drug.responses), ]
            }
            doses[[i]] <- drug.responses$Dose
            responses[[i]] <- drug.responses$Viability
            names(doses[[i]]) <- names(responses[[i]]) <- 1:length(doses[[i]])
            if (!missing(legends.label)) {
              if(length(legends.label)>0) {
                linQuad_params <- linearQuadraticModel(D = doses[[i]], SF = responses[[i]])
                if(any(grepl("alpha", x=legends.label))){
                  legend.values[[i]] <- paste(legend.values[i][[1]],sprintf("%s = %s", "alpha", round(linQuad_params[1], digits=2)), sep=", ")
                }
                if(any(grepl("beta", x=legends.label))){
                  legend.values[[i]] <- paste(legend.values[i][[1]],sprintf("%s = %s", "beta", round(linQuad_params[2], digits=2)), sep=", ")
                }
                if(any(grepl("rsquared", x=legends.label))){
                  legend.values[[i]] <- paste(legend.values[i][[1]],sprintf("%s = %s", "R^2", round(CoreGx::examineGOF(linQuad_params)[1], digits=2)), sep=", ")
                }
              } else {
                legend.values[[i]] <- ""
              }
            }
          } else {
            for (exp in exp_i) {
              j <- j + 1
              rSetNames[[j]] <- rSetName(rSets[[i]])

              drug.responses <- as.data.frame(cbind("Dose"=as.numeric(as.vector(rSets[[i]]@sensitivity$raw[exp, , "Dose"])),
                "Viability"=as.numeric(as.vector(rSets[[i]]@sensitivity$raw[exp, , "Viability"])), stringsAsFactors=FALSE))
              drug.responses <- drug.responses[complete.cases(drug.responses), ]
              doses[[j]] <- drug.responses$Dose
              responses[[j]] <- drug.responses$Viability
              names(doses[[j]]) <- names(responses[[j]]) <- 1:length(doses[[j]])
              if (!missing(legends.label)) {
                if(length(legends.label)>0){
                  linQuad_params <- linearQuadraticModel(D = doses2[[i]], SF = responses2[[i]])
                  if(any(grepl("alpha", x=legends.label))){
                    legend.values2[[i]] <- paste(legend.values2[i][[1]],sprintf("%s = %s", "alpha", round(linQuad_params[1], digits=2)), sep=", ")
                  }
                  if(any(grepl("beta", x=legends.label))){
                    legend.values2[[i]] <- paste(legend.values2[i][[1]],sprintf("%s = %s", "beta", round(linQuad_params[2], digits=2)), sep=", ")
                  }
                  if(any(grepl("rsquared", x=legends.label))){
                    legend.values2[[i]] <- paste(legend.values2[i][[1]],sprintf("%s = %s", "R^2", round(CoreGx::examineGOF(linQuad_params)[1], digits=2)), sep=", ")
                  }
                }
              } else {
                tt <- unlist(strsplit(rownames(rSets[[i]]@sensitivity$info)[exp], split="_"))
                if (tt[1] == "radiation.type") {
                  legend.values[[j]] <- tt[2]
                }else{
                  legend.values[[j]] <- rownames(rSets[[i]]@sensitivity$info)[exp]
                }
              }
            }
          }
        } else {
          warning("The cell line and drug combo were not tested together. Aborting function.")
          return()
        }
      }
    }
    if(!missing(Ds)){
      doses2 <- list(); responses2 <- list(); legend.values2 <- list(); j <- 0; rSetNames2 <- list();
      for (i in 1:length(Ds)){
        doses2[[i]] <- Ds[[i]]
        responses2[[i]] <- SFs[[i]]
        if(length(legends.label)>0){
          linQuad_params <- linearQuadraticModel(D = doses2[[i]], SF = responses2[[i]])
          if(any(grepl("alpha", x=legends.label))){
            legend.values2[[i]] <- paste(legend.values2[i][[1]],sprintf("%s = %s", "alpha", round(linQuad_params[1], digits=2)), sep=", ")
          }
          if(any(grepl("beta", x=legends.label))){
            legend.values2[[i]] <- paste(legend.values2[i][[1]],sprintf("%s = %s", "beta", round(linQuad_params[2], digits=2)), sep=", ")
          }
          if(any(grepl("rsquared", x=legends.label))){
            legend.values2[[i]] <- paste(legend.values2[i][[1]],sprintf("%s = %s", "R^2", round(CoreGx::examineGOF(linQuad_params)[1], digits=2)), sep=", ")
          }
        } else{ legend.values2[[i]] <- ""}

        rSetNames2[[i]] <- names(Ds)[[i]]
      }
      doses <- c(doses, doses2)
      responses <- c(responses, responses2)
      legend.values <- c(legend.values, legend.values2)
      rSetNames <- c(rSetNames, rSetNames2)
    }

    if (missing(mycol)) {
    # require(RColorBrewer) || stop("Library RColorBrewer is not available!")
      mycol <- RColorBrewer::brewer.pal(n=7, name="Set1")
    }

    dose.range <- c(10^100 , 0)
    viability.range <- c(1 , 1)
    for(i in 1:length(doses)) {
      # dose.range <- c(min(dose.range[1], min(doses[[i]], na.rm=TRUE), na.rm=TRUE), max(dose.range[2], max(doses[[i]], na.rm=TRUE), na.rm=TRUE))
      dose.range <- c(0, max(dose.range[2], max(doses[[i]], na.rm=TRUE), na.rm=TRUE))
      viability.range <- c(min(viability.range[1], min(responses[[i]], na.rm=TRUE), na.rm=TRUE), 1)
    }
    x1 <- 10 ^ 10; x2 <- 0

  # if(length(doses) > 1) {
  #   common.ranges <- PharmacoGx::.getCommonConcentrationRange(doses)

  #   for(i in 1:length(doses)) {
  #     x1 <- min(x1, min(common.ranges[[i]]))
  #     x2 <- max(x2, max(common.ranges[[i]]))
  #   }
  # }
    if (!missing(xlim)) {
      dose.range <- xlim
    }
    if (!missing(ylim)) {
      viability.range <- ylim
    }
    if(missing(title)){
    ## FIXME:: HACK
    # if(!missing(drug)&&!missing(cellline)){
    #   title <- sprintf("%s:%s", drug, cellline)
    # } else {
      if (length(rSets)){
        title <- sprintf("Radiation Response Curve for: %s", cellline)
      } else {
        title <- "Radiation Response Curve"
      }
    # }
    }
    plot(NA, xlab="Dose (Gy)", ylab="Survival Fraction", axes =FALSE, main=title, log="y", ylim=viability.range, xlim=dose.range, cex=cex, cex.main=cex.main)
    magicaxis::magaxis(side=1:2, frame.plot=TRUE, tcl=-.3, majorn=c(5,5), minorn=c(5,3), label=c(TRUE,FALSE))
    if(max(viability.range)/min(viability.range)<50){
      ticks <- magicaxis::maglab(viability.range, exptext = TRUE)
    } else {
      ticks <- magicaxis::maglab(viability.range, exptext = TRUE, log=TRUE)
    }
    ticks$exp <- sapply(ticks$exp, function(x) return(as.expression(bquote(10^ .(round(log10(eval(x)), 2))))))
    axis(2, at=ticks$labat,labels=ticks$exp)
    legends <- NULL
    legends.col <- NULL
  # if (length(doses) > 1) {
  #   rect(xleft=x1, xright=x2, ybottom=viability.range[1] , ytop=viability.range[2] , col=rgb(240, 240, 240, maxColorValue = 255), border=FALSE)
  # }

    for (i in 1:length(doses)) {
      points(doses[[i]],responses[[i]],pch=20,col = mycol[i], cex=cex)

      switch(plot.type , "Actual"={
        lines(doses[[i]], responses[[i]], lty=1, lwd=lwd, col=mycol[i])
      }, "Fitted"={
        linQuad_params <- linearQuadraticModel(D = doses[[i]], SF = responses[[i]])
        x_vals <- CoreGx::.GetSupportVec(c(0,doses[[i]]))
        lines(x_vals, (.linearQuadratic(x_vals, pars=linQuad_params, SF_as_log=FALSE)),lty=1, lwd=lwd, col=mycol[i])
      },"Both"={
      # lines(doses[[i]],responses[[i]],lty=1,lwd=lwd,col = mycol[i])
        linQuad_params <- linearQuadraticModel(D = doses[[i]], SF = responses[[i]])
        x_vals <- CoreGx::.GetSupportVec(c(0,doses[[i]]))
        lines(x_vals, (.linearQuadratic(x_vals, pars=linQuad_params, SF_as_log=FALSE)),lty=1, lwd=lwd, col=mycol[i])
      })
      if (length(legend.values)){
              legends<- c(legends, sprintf("%s%s", rSetNames[[i]], legend.values[[i]]))
      } else {
                      legends<- c(legends, sprintf("%s", rSetNames[[i]]))
      }
      legends.col <-  c(legends.col, mycol[i])
    }
  # if (common.range.star) {
  #   if (length(doses) > 1) {
  #     for (i in 1:length(doses)) {
  #       points(common.ranges[[i]], responses[[i]][names(common.ranges[[i]])], pch=8, col=mycol[i])
  #     }
  #   }
  # }
    legend(legend.loc, legend=legends, col=legends.col, bty="n", cex=cex, pch=c(15,15))
    return(invisible(NULL))
  }

