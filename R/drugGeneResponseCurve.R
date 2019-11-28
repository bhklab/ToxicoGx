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
#' @param xlim \code{numeric} A vector of minimum and maximum values for the x-axis
#'   of the returned plot.
#' @param ylim \code{numeric} A vector of minimum and miximum values for the y-axis
#'   of the returned plot.
#' @param title \code{character} A string containing the desired plot name. If excluded
#'   a title wil be generated automatically.
#' @param legend.loc \code{character} The location of the legend as passed to the plot()
#'   function from base graphics. Suggested values are either "topleft" or
#'   "topright", with the latter as the default.
#' @param mycol `vector`` A vector of length equal to the product of the
#'   number of drugs, features and doses passed to the function. Takes colour
#'   arguments as passed to `col` parameter in the `plot()` function.
#'   Default palette is used when unspecified.
#' @param summarize.replicates \code{logical} If true will take the average of all
#'  replicates at each time point per gene and duration. This release has not
#'  yet implemented this feature.
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

  plotData <- lapply(tSet, function(tSet) {
    lapply(mDataTypes, function(mDataType) {
      mProf <- molecularProfiles(tSet, mDataType)
      list(
        "data" = data.table(
          mProf,
          keep.rownames = TRUE
        ),
        "fInfo" = data.table(featureInfo(tSet, mDataType)),
        "pInfo" = data.table(phenoInfo(tSet, mDataType))
      )
    })
  })

  # Get a list of times per tSet per mDataType per dose level
  # This will also need to be per drug if we extend the function to multiple drugs
  times <- lapply(plotData, function(tData) {
    m <- lapply(tData, function(mData) {
        lapply(split(mData$pInfo[, unique(duration), by = dose_level], by = "dose_level"),
               function(x) { x$V1 })
               })
    names(m) <- mDataTypes; m
  })
  names(times) <- vapply(tSet, function(x) names(x), FUN.VALUE = character(1)) # Get the names for each tSet

  ## TODO:: Rewrite legendValues to generate after summarization; will remove a bunch of code
  # Assembling the legend names for each line to be plotted
  legendValues <- lapply(plotData, function(tData) {
    m <- lapply(tData, function(mData) {
      dose_rep <- split(mData$pInfo[, paste(dose_level,unique(individual_id), sep = "_"), by = dose_level],
        by = "dose_level")
      lapply(dose_rep, function(dLevel) {
          paste(gsub("-[^-]*$", "", mData$fInfo[, unique(transcript_name)]), dLevel$V1, sep = "_")
        })
    })
    names(m) <- mDataTypes; m # Note: FUN.VALUE refers to the type and length of EACH function call, not of the returned vector
  })
  names(legendValues) <-  vapply(tSet, function(x) names(x), FUN.VALUE = character(1))

  # Get expression values per tSet, mDataType, dose, time, gene and replicate
  expression <-
    list.map(plotData, f(t) ~
      list.map(t, f(m) ~
        lapply(split(m$pInfo[, as.character(samplename), by = .(dose_level, duration)], by = "dose_level"), function(i) {
          lapply((unname(split(i, by = "duration", )) %>% list.map(f(x) ~ x$V1)),
                 function(s) {
                   split(m$data[, c("rn", s), by = rn, with = F], by = "rn") %>%
                     list.map(f(x) ~  as.numeric(unlist(x), use.names = F)[-1])
             })
          })
       )
    )

  # Collapse the number of Control replicates in the plot to match the does level with the highest number of replicates
  if (any(vapply(tSet, function(tSet) { names(tSet) == "drugMatrix"}, FUN.VALUE = logical(1)))) {
    max <- list.flatten(expression) %>% list.map(i ~ if (!grepl('Control.*', .name)) length(i) else 0) %>% unlist() %>% max()
    expression <-
      list.map(expression, f(t) ~
        list.map(t, f(m) ~
          list.map(m, f(d) ~
            if (.name == 'Control') list.map(d, f(i) ~ list.map(i, f(s) ~ rep(mean(s), max))) else d
          )
        )
      )
    legendValues <-
      list.map(legendValues, t ~
        list.map(t, m ~
          list.map(m, d ~
            if (.name == 'Control') d[seq_len(max)] else d
            )
          )
        )
  }

  #### SUMMARIZATION ####
  if (summarize.replicates) {

    expression <- list.map(expression, f(t) ~
                    list.map(t, f(m) ~
                      list.map(m, f(d) ~
                        list.map(d, f(t) ~
                          list.map(t, f(f) ~
                            mean(f)
                                )
                              )
                            )
                          )
                        )

    legendValues <- list.map(legendValues, f(x) ~
                      list.map(x, f(x) ~
                        list.map(x, f(x) ~
                          unique(gsub("_[^_]*$", '', x)
                            )
                          )
                        )
                      )
  }

  #### AXIS RANGES ####

  # Set x and y axis ranges based on time and viability values
  time.range <- as.numeric(c(min(unlist(times)), max(unlist(times))))
  expression.range <- c(floor(min(unlist(expression))), ceiling(max(unlist(expression))))
  for (i in seq_along(tSet)) {
    ## TODO:: Generalize this to n replicates
    time.range <- as.numeric(c(min(time.range[1], min(unlist(times[[i]]), na.rm = TRUE), na.rm = TRUE), max(time.range[2], max(unlist(times[[i]]), na.rm = TRUE), na.rm = TRUE)))
    expression.range <- c(0, max(expression.range[2], max(unlist(expression[[i]]), na.rm = TRUE), na.rm = TRUE))
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
      title <- sprintf("%s\n%s:%s", "Expression Time Response Curve", paste(drug, collapse = " & "), paste(cell.lines, collapse = " & "))
  }

  ## SETS DEFAULT COLOUR PALETTE
  if (missing(mycol)) {
    mycol <- c(RColorBrewer::brewer.pal(n = 9, name = "Set1"), RColorBrewer::brewer.pal(n = 12, name = "Set3"))
    legends.col <- mycol
  }

  #### DRAWING THE PLOT ####

  # Reset par after function compeltes
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))

  # Modify par in functions scope

  b <- l <- t <- r <- 0
  if (grepl("^bottom$", legend.loc)) {b <- 10; inset = c(0, -0.7)}
  if (grepl(".*left.*", legend.loc)) {l <- 10; inset = c(-0.4, 0)}
  if (grepl("^top$", legend.loc)) {t <- 10; inset <- c(0, -0.7)}
  if (grepl(".*right,*", legend.loc)) {r <- 10; inset <- c(-0.3, 0)}
  par(xpd = T, mar = par()$mar + c(b,l,t,r))

  # Create the plot frame
  plot(NA, xlab = "Time (hr)", ylab = "Expression", axes = FALSE, main = title,
       ylim = expression.range, xlim = time.range, cex = cex, cex.main = cex.main)
  box()

  # Adds plot axes
  axis(1, at = as.numeric(duration), tcl = -0.3)
  axis(2, at = seq.int(expression.range[1], expression.range[2]), tcl = -0.3)

  # Initialize legends variables
  legends <- NULL
  pch.val <- NULL
  legends.col <- NULL

  list.map(expression)

}


# if (!summarize.replicates) { # i.e., if summarize replicates is false
#   # Loop over tSet
#   k <- j <- 1
#   for (i in seq_along(times)) { # tSet subset
#     for (mDataType in seq_along(mDataTypes)) {
#       # Loop over dose level
#       for (level in seq_along(dose)) {
#         # Loop over replicates per dose level
#         ## TODO:: Generalize this for n replicates
#           if (length(seq_along(features[[mDataType]])) > 2 ) { j <- 1 }
#           # Plot per tSet, per dose level, per replicate points
#           for (feature in seq_along(features[[mDataType]])) {
#             for (replicate in seq_along(expression[[i]][[mDataType]][[level]])) {
#             points(times[[i]][[mDataType]][[level]],
#                    expression[[i]][[mDataType]][[level]][[replicate]][[feature]],
#                    pch = k, col = mycol[j], cex = cex)
#             # Select plot type
#             lines(times[[i]][[mDataType]][[level]],
#                   expression[[i]][[mDataType]][[level]][[replicate]][[feature]],
#                   lty = replicate, lwd = lwd, col = mycol[j])
#             legends <- c(legends,
#                          legendValues[[i]][[mDataType]][[level]][[replicate]][feature])
#             legends.col <- c(legends.col, mycol[j])
#             pch.val <- c(pch.val, k)
#             k <- k + 1
#           }
#           if (length(seq_along(features[[mDataType]])) > 1) { j <- j + 1 }
#         }
#         j <- j + 1
#       }
#     }
#   }
# } else {
#   # Loop over tSet
#   k <- j <- 1
#   for (i in seq_along(times)) { # tSet subset
#     for (mDataType in seq_along(mDataTypes)) {
#       # Loop over dose level
#       if (length(features[[i]][[mDataType]]) > 1) {j <- 1 }
#       for (level in seq_along(dose)) {
#         # Loop over replicates per dose level
#         ## TODO:: Generalize this for n replicates
#         # Plot per tSet, per dose level, per replicate points
#         for (feature in seq_along(features[[mDataType]])) {
#           points(times[[i]][[mDataType]][[level]],
#                  expression[[i]][[mDataType]][[level]][[feature]],
#                  pch = k, col = mycol[j], cex = cex)
#           # Select plot type
#           lines(times[[i]][[mDataType]][[level]],
#                 expression[[i]][[mDataType]][[level]][[feature]],
#                 lty = 1, lwd = lwd, col = mycol[j])
#           legends <- c(legends,
#                        legendValues[[i]][[mDataType]][[level]][feature])
#           legends.col <- c(legends.col, mycol[j])
#           pch.val <- c(pch.val, k)
#           k <- k + 1
#           j <- j + 1
#         }
#       if (length(features[[i]][[mDataType]]) > 1) {j <- j + 1 }
#       }
#     }
#   }
# }
# legend(legend.loc, legend = legends, col = legends.col, bty = "L", cex = cex, pch = pch.val, inset = inset)
