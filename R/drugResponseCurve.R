#' Plot drug response curve of a given drug and a given cell for a list of tSets (objects of the ToxicoSet class).
#'
#' Given a list of ToxicoSets, the function will plot the drug_response curve,
#' for a given drug/cell pair. The y axis of the plot is the viability percentage
#' and x axis is the log transformed concentrations. If more than one tSet is
#' provided, a light gray area would show the common concentration range between tSets.
#' User can ask for type of sensitivity measurment to be shown in the plot legend.
#' The user can also provide a list of their own concentrations and viability values,
#' as in the examples below, and it will be treated as experiments equivalent to values coming
#' from a tSet. The names of the concentration list determine the legend labels.
#'
#' @examples
#' if (interactive()) {
#' drugResponseCurve(concentrations=list("Experiment 1"=c(.008, .04, .2, 1)),
#'  viabilities=list(c(100,50,30,1)), plot.type="Both")
#' }
#'
#' @param drug [string] A drug name for which the drug response curve should be
#' plotted. If the plot is desirable for more than one toxico set, A unique drug id
#' should be provided.
#' @param cellline [string] A cell line name for which the drug response curve should be
#' plotted. If the plot is desirable for more than one toxico set, A unique cell id
#' should be provided.
#' @param durations [numeric] A duration for which the drug response curve should be plotted.
#' @param tSets [list] a list of ToxicoSet objects, for which the function
#' should plot the curves.
#' @param concentrations,viabilities [list] A list of concentrations and viabilities to plot, the function assumes that
#' concentrations[[i]] is plotted against viabilities[[i]]. The names of the concentration list are used to create the legend labels
#' @param conc_as_log [logical], if true, assumes that log10-concentration data has been given rather than concentration data,
#' and that log10(ICn) should be returned instead of ICn. Applies only to the concentrations parameter.
#' @param viability_as_pct [logical], if false, assumes that viability is given as a decimal rather
#' than a percentage, and that E_inf passed in as decimal. Applies only to the viabilities parameter.
#' @param legends.label [vector] A vector of sensitivity measurment types.
#' A legend will be displayed on the top right of the plot which each line of the legend is
#' the values of requested sensitivity measurements for one of the requested tSets.
#' If this parameter is missed no legend would be provided for the plot.
#' @param ylim [vector] A vector of two numerical values to be used as ylim of the plot.
#' If this parameter would be missed c(0,100) would be used as the ylim of the plot.
#' @param xlim [vector] A vector of two numerical values to be used as xlim of the plot.
#' If this parameter would be missed the minimum and maximum comncentrations between all
#' the tSets would be used as plot xlim.
#' @param mycol [vector] A vector with the same length of the tSets parameter which
#' will determine the color of the curve for the toxico sets. If this parameter is
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
#'
#' @return Plots to the active graphics device and returns and invisible NULL.
#'
#' @import RColorBrewer
#' @importFrom graphics plot rect points lines legend
#' @importFrom grDevices rgb
#' @importFrom magicaxis magaxis
#'
#' @export
#'
drugResponseCurve <-
  ## TODO:: Note that all arguments must default to null for the paramErrorChecker to work
  ##   The logic for the user remains unchanged.
  function(drug=NULL,
           cellline=NULL,
           durations=NULL,
           tSets=NULL,
           concentrations=NULL,
           viabilities=NULL,
           conc_as_log = FALSE,
           viability_as_pct = TRUE,
           trunc=TRUE,
           legends.label = c("ic50_published", "gi50_published","auc_published","auc_recomputed","ic50_recomputed"),
           ylim=c(0,100),
           xlim, mycol,
           title,
           plot.type=c("Fitted","Actual", "Both"),
           summarize.replicates=TRUE,
           lwd = 0.5,
           cex = 0.5,
           cex.main = 0.9,
           legend.loc = "topright",
           verbose=TRUE) {

    ## HANDLE INCORRECT ARGUMENT TYPES
    if (!is.null(tSets)) { if(!is(tSets, "list")) { tSets <- list(tSets) }}

    #### CHECK PARAMETERS MEET FUNCTION REQUIREMENTS ####
    paramErrorChecker("drugDoseResponseCurve",
                      tSets=tSets, concentrations=concentrations,
                      viabilities=viabilities, cell.lines=cellline, drugs=drug,
                      duration=durations)


    #### SANITIZE CONCENTRATION INPUTS ####
    if(!is.null(concentrations)){ #if a concentrations argument has been passed in
      ## IF CONCENTRATIONS IS A VECTOR

      ## TODO:: I think much of this logic is unnecessary due to concentration being discrete in ToxicoGx
      if (!is(concentrations, "list")) {
          #sanitizeInput returns a list of length 2, where [[1]] is the conc, [[2]] is the viabilities as the user requested
          cleanData <- sanitizeInput(concentrations,
                                     viabilities,
                                     conc_as_log = conc_as_log,
                                     viability_as_pct = viability_as_pct,
                                     trunc = trunc,
                                     verbose = verbose)
          #well idk what it's doing here but ok, the concentrations and viabilities are turned into lists and given names
          concentrations <- 10^cleanData[["log_conc"]]
          concentrations <- list(concentrations)
          viabilities <- 100*cleanData[["viability"]]
          viabilities <- list(viabilities)
          names(concentrations) <- "Exp1"
          names(viabilities) <- "Exp1"
      ## IF CONCENTRATIONS IS LIST
      } else {
        if(is.null(names(concentrations))){ # Assign numbered names if no names provided
          names(concentrations) <- paste("Exp", seq_along(concentrations))
        }
        # Loop over each concentration vector in the list
        ## TODO:: Benchmark against apply statement
        for(i in seq_along(concentrations)){
          # if multiple concentrations were passed in
            # sanitize input for each of the multiple concentrations
            cleanData <- sanitizeInput(concentrations[[i]],
                                       viabilities[[i]],
                                       conc_as_log = conc_as_log,
                                       viability_as_pct = viability_as_pct,
                                       trunc = trunc,
                                       verbose = verbose)
            concentrations[[i]] <- 10^cleanData[["log_conc"]]
            viabilities[[i]] <- 100*cleanData[["viability"]]
        }
      }
    }

    ## TODO:: This appear to do nothing? It is only assigned once as false
    common.range.star <- FALSE

    # SETS DEFAULT PLOT TYPE
    if (missing(plot.type)) { plot.type <- "Actual" }

    ## GETS CONCENTRATIONS AND VIABILITIES FROM tSET BASED ON INTERSECTION
    doses <- list(); responses <- list(); legend.values <- list(); j <- 0; tSetNames <- list(); tSetNames_temp <- list(); doses_temp <- list(); responses_temp <- list(); legend.values_temp <- list();
    if(!is.null(tSets)){
      for(i in seq_along(tSets)) {
        # exp_i contains the indices of sensitivity object rows that correspond to the requested cell line, drug
        exp_i <- which(sensitivityInfo(tSets[[i]])[ ,"cellid"] == cellline & sensitivityInfo(tSets[[i]])[ ,"drugid"] == drug)
        if(length(exp_i) > 0) { #if there is a UID that corresponds to the cell line - drug combo
          if (summarize.replicates) { #if replicates should be summarized
            if (length(exp_i) == 1) { #if there is only 1 UID that corresponds to the cell line - drug combo
              tSetNames[[i]] <- tSetName(tSets[[i]])
              drug.responses <- as.data.frame(cbind("Dose"=as.numeric(as.vector(tSets[[i]]@sensitivity$raw[exp_i, colnames(tSets[[i]]@sensitivity$raw[,,"Dose"]) != "Control", "Dose"])),
                                                    "Viability"=as.numeric(as.vector(tSets[[i]]@sensitivity$raw[exp_i, colnames(tSets[[i]]@sensitivity$raw[,,"Viability"]) != "Control", "Viability"])), stringsAsFactors=FALSE))
              drug.responses <- drug.responses[complete.cases(drug.responses), ]
            }else{ #if there are multiple UIDs that correspond to the cell line - drug combo
              drug.responses <- data.frame()
              k = 0
              for (j in durations){
                k = k + 1
                tSetNames_temp[[k]] <- tSetName(tSets[[i]])
                exp_j <- which(sensitivityInfo(tSets[[i]])[ ,"cellid"] == cellline & sensitivityInfo(tSets[[i]])[ ,"drugid"] == drug &
                                 sensitivityInfo(tSets[[i]])[ ,"duration_h"] == j)

                drug.responses <- as.data.frame(cbind("Dose"=apply(tSets[[i]]@sensitivity$raw[exp_j, colnames(tSets[[i]]@sensitivity$raw[,,"Dose"]) != "Control", "Dose"], 2, function(x){median(as.numeric(x), na.rm=TRUE)}),
                                                      "Viability"=apply(tSets[[i]]@sensitivity$raw[exp_j, colnames(tSets[[i]]@sensitivity$raw[,,"Viability"]) != "Control", "Viability"], 2, function(x){median(as.numeric(x), na.rm=TRUE)}), stringsAsFactors=FALSE))
                rownames(drug.responses) <- paste(drug.responses,"hr_",rownames(drug.responses), sep = "")
                drug.responses <- drug.responses[complete.cases(drug.responses), ]
                doses_temp[[k]] <- drug.responses$Dose
                responses_temp[[k]] <- drug.responses$Viability
                names(doses_temp[[k]]) <- names(responses_temp[[k]]) <- 1:length(doses_temp[[k]])

                if (!missing(legends.label)) { #if the user specified a legend label
                  # legend.values_temp contains the values for all sensitivity measurements corresponding to the experiment UIDs
                  if (length(legends.label) > 1) {
                    legend.values_temp[[i]] <- paste(unlist(lapply(legends.label, function(x){
                      sprintf("%s = %s", x, round(as.numeric(tSets[[i]]@sensitivity$profiles[exp_j,x]), digits=2))
                    })), collapse = ", ")
                  } else {
                    legend.values_temp[[i]] <- sprintf("%s = %s", legends.label, round(as.numeric(tSets[[i]]@sensitivity$profiles[exp_j, legends.label]), digits=2))
                  }
                } else { #if the user did not specify a legend label
                  legend.values_temp[[i]] <- paste("_",j,"hr", sep = "")
                }
              }
              tSetNames <- c(tSetNames, tSetNames_temp)
              doses <- c(doses, doses_temp)
              responses <- c(responses, responses_temp)
              legend.values <- c(legend.values, legend.values_temp)
            }
          }else { #if replicates should not be summarized
            for (exp in exp_i) {
              j <- j + 1
              tSetNames[[j]] <- tSetName(tSets[[i]])

              drug.responses <- as.data.frame(cbind("Dose"=as.numeric(as.vector(tSets[[i]]@sensitivity$raw[exp, colnames(tSets[[i]]@sensitivity$raw[,,"Dose"]) != "Control", "Dose"])),
                                                    "Viability"=as.numeric(as.vector(tSets[[i]]@sensitivity$raw[exp, colnames(tSets[[i]]@sensitivity$raw[,,"Viability"]) != "Control", "Viability"])), stringsAsFactors=FALSE))
              drug.responses <- drug.responses[complete.cases(drug.responses), ]
              doses[[j]] <- drug.responses$Dose
              responses[[j]] <- drug.responses$Viability
              names(doses[[j]]) <- names(responses[[j]]) <- 1:length(doses[[j]])
              if (!missing(legends.label)) {
                if (length(legends.label) > 1) {
                  legend.values[[j]] <- paste(unlist(lapply(legends.label, function(x){
                    sprintf("%s = %s", x, round(as.numeric(tSets[[i]]@sensitivity$profiles[exp, x]), digits=2))
                  })), collapse = ", ")
                } else {
                  legend.values[[j]] <- sprintf(" Exp %s %s = %s", rownames(tSets[[i]]@sensitivity$info)[exp], legends.label, round(as.numeric(tSets[[i]]@sensitivity$profiles[exp, legends.label]), digits=2))
                }
              } else {
                tt <- unlist(strsplit(rownames(tSets[[i]]@sensitivity$info)[exp], split="_"))
                if (tt[1] == "drugid") {
                  legend.values[[j]] <- paste(tt[2],"_", tt[5], sep = "")
                }else{
                  legend.values[[j]] <- rownames(tSets[[i]]@sensitivity$info)[exp]
                }
              }
            }
          }

          for (i in length(doses)){
            names(doses[[i]]) <- names(responses[[i]]) <- i
          }
        } else {
          warning("The cell line and drug combo were not tested together. Aborting function.")
          return()
        }
      }
    }

    if(!is.null(concentrations)){
      doses2 <- list(); responses2 <- list(); legend.values2 <- list(); j <- 0; tSetNames2 <- list();
      for (i in seq_along(concentrations)){
        doses2[[i]] <- concentrations[[i]]
        responses2[[i]] <- viabilities[[i]]
        if(length(legends.label)>0){
          if(any(grepl("AUC", x=toupper(legends.label)))){
            legend.values2[[i]] <- paste(legend.values2[i][[1]],sprintf("%s = %s", "AUC", round(computeAUC(concentrations[[i]],viabilities[[i]], conc_as_log=FALSE, viability_as_pct=TRUE)/100, digits=2)), sep=", ")
          }
          if(any(grepl("IC50", x=toupper(legends.label)))){
            legend.values2[[i]] <- paste(legend.values2[i][[1]],sprintf("%s = %s", "IC50", round(computeIC50(concentrations[[i]],viabilities[[i]], conc_as_log=FALSE, viability_as_pct=TRUE), digits=2)), sep=", ")
          }

        } else{ legend.values2[[i]] <- ""}

        tSetNames2[[i]] <- names(concentrations)[[i]]
      }
      doses <- c(doses, doses2)
      responses <- c(responses, responses2)
      legend.values <- c(legend.values, legend.values2)
      tSetNames <- c(tSetNames, tSetNames2)
    }

    ## SETS DEFAULT COLOUR PALETTE
    if (missing(mycol)) {
      mycol <- RColorBrewer::brewer.pal(n=9, name="Set1")
    }

    #
    dose.range <- c(10^100 , 0)
    viability.range <- c(0 , 10)
    for(i in seq_along(doses)) {
      dose.range <- c(min(dose.range[1], min(doses[[i]], na.rm=TRUE), na.rm=TRUE), max(dose.range[2], max(doses[[i]], na.rm=TRUE), na.rm=TRUE))
      viability.range <- c(0, max(viability.range[2], max(responses[[i]], na.rm=TRUE), na.rm=TRUE))
    }
    x1 <- 10 ^ 10; x2 <- 0

    ## FINDS INTERSECTION OF RANGES IF MORE THAN ONE EXPERIMENT PLOTTED
    if(length(doses) > 1) {
      common.ranges <- .getCommonConcentrationRange(doses)

      for(i in seq_along(doses)) {
        x1 <- min(x1, min(common.ranges[[i]]))
        x2 <- max(x2, max(common.ranges[[i]]))
      }
    }
    # SETS CUSTOM RANGE FOR X-AXIS IF PASSED AS ARGUEMENT
    if (!missing(xlim)) {
      dose.range <- xlim
    }
    ## SETS CUSTOM RANGE FOR Y-AXIS IF PASSED AS ARGUEMENT
    if (!missing(ylim)) {
      viability.range <- ylim
    }

    ## SETS PLOT TITLE
    if(missing(title)){
      if(!missing(drug)&&!missing(cellline)){
        title <- sprintf("%s:%s", drug, cellline)
      } else {
        title <- "Drug Dose Response Curve"
      }
    }

    #### DRAWING THE PLOT ####
    plot(NA, xlab="Concentration (uM)", ylab="% Viability", axes =FALSE, main=title, log="x", ylim=viability.range, xlim=dose.range, cex=cex, cex.main=cex.main)
    magicaxis::magaxis(side=1:2, frame.plot=TRUE, tcl=-.3, majorn=c(5,3), minorn=c(5,2))
    legends <- NULL
    legends.col <- NULL
    if (length(doses) > 1) {
      rect(xleft=x1, xright=x2, ybottom=viability.range[1] , ytop=viability.range[2] , col=rgb(240, 240, 240, maxColorValue = 255), border=FALSE)
    }
    for (i in seq_along(doses)) {
      points(doses[[i]],responses[[i]],pch=20,col = mycol[i], cex=cex)
      switch(plot.type,
        "Actual"={
          lines(doses[[i]], responses[[i]], lty=1, lwd=lwd, col=mycol[i])
          },
        "Fitted"={
          log_logistic_params <- logLogisticRegression(conc=doses[[i]], viability=responses[[i]])
          log10_x_vals <- .GetSupportVec(log10(doses[[i]]))
          lines(10 ^ log10_x_vals, .Hill(log10_x_vals, pars=c(log_logistic_params$HS, log_logistic_params$E_inf/100, log10(log_logistic_params$EC50))) * 100 ,lty=1, lwd=lwd, col=mycol[i])
          },
        "Both"={
          lines(doses[[i]],responses[[i]],lty=1,lwd=lwd,col = mycol[i])
          log_logistic_params <- logLogisticRegression(conc = doses[[i]], viability = responses[[i]])
          log10_x_vals <- .GetSupportVec(log10(doses[[i]]))
          lines(10 ^ log10_x_vals, .Hill(log10_x_vals, pars=c(log_logistic_params$HS, log_logistic_params$E_inf/100, log10(log_logistic_params$EC50))) * 100 ,lty=1, lwd=lwd, col=mycol[i])
          })
      legends<- c(legends, sprintf("%s%s", tSetNames[[i]], legend.values[[i]]))
      legends.col <- c(legends.col, mycol[i])
    }
    if (common.range.star) {
      if (length(doses) > 1) {
        for (i in seq_along(doses)) {
          points(common.ranges[[i]], responses[[i]][names(common.ranges[[i]])], pch=8, col=mycol[i])
        }
      }
    }
    legend(legend.loc, legend=legends, col=legends.col, bty="n", cex=cex, pch=c(15,15))
    ## TODO:: Don't think we need a return statement
    #return(invisible(NULL))
}
