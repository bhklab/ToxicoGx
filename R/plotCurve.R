#' Plot radiation dose-response curve
#'
#' This function plots doses of radiation against the cancer cell survival fractions thereby observed.
#'
#' @param D vector of radiation doses
#' @param SF vector of survival fractions corresponding to the doses
#' @param pars parameters (alpha, beta) in the equation SF = exp(-alpha * D - beta * D ^ 2)
#' @param filename name of PDF which will be created by the function
#' @param fit_curve should the graph include a linear-quadratic curve of best fit? Defaults to TRUE
#' @param SF_as_log should SF be expressed in log10 on the graph? Defaults to TRUE
#' @examples plotCurve(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
#'     c(1.1, 0.8, 0.7, 0.45, 0.15, -0.1, -0.1, -0.4, -0.65, -0.75, -1.1),
#'         filename = NULL)
#' @importFrom graphics lines plot points axis
#' @importFrom grDevices dev.off pdf
#' @export
plotCurve <- function(D, SF, pars, filename = "dose_response_plot.pdf", fit_curve = TRUE, SF_as_log = TRUE) {
  CoreGx::.sanitizeInput(x = D,
                          y = SF,
                          x_as_log = FALSE,
                          y_as_log = FALSE,
                          y_as_pct = FALSE,
                          trunc = FALSE,
                          verbose = FALSE)

  padding <- 1.1 # whitespace on graph around function range

  if (fit_curve) {
    if (missing(pars)) {
      pars <- unlist(linearQuadraticModel(D, SF))
    } else {
      CoreGx::.sanitizeInput(pars = pars,
                              x_as_log = FALSE,
                              y_as_log = FALSE,
                              y_as_pct = FALSE,
                              trunc = FALSE,
                              verbose = FALSE)
    }
    print(paste0("A linear-quadratic curve was fit to the data with parameters alpha = ", pars[[1]], " and beta = ", pars[[2]], "."))
    trendlineDs <- CoreGx::.GetSupportVec(D)
    trendlineSFs <- .linearQuadratic(trendlineDs, pars = pars, SF_as_log = TRUE)
  }

  xlim <- range(D)
  xlim <- mean(xlim) + padding * c((xlim[1] - mean(xlim)), xlim[2] - mean(xlim))
  if (!missing(SF)) {
    DSF <- CoreGx::.reformatData(x = D,
                                 y = SF,
                                 x_to_log = FALSE,
                                 y_to_log = TRUE,
                                 y_to_frac = FALSE,
                                 trunc = FALSE)
    D <- DSF[["x"]]
    SF <- DSF[["y"]]
  }


  if (TRUE) {
    if (!missing(SF)) {
      if (fit_curve) {
        ylim <- padding * c(min(c(SF, trendlineSFs[length(trendlineSFs)])), 0)
      } else {
        ylim <- padding * c(min(SF), 0)
      }
    } else {
      ylim <- padding * c(trendlineSFs[length(trendlineSFs)], 0)
    }
  } else {
    ylim <- c(0, 1)
  }

  pdf(file = filename)

  plot(NULL,
       xlim = xlim,
       ylim = ylim,
       xlab = "Dose (Gy)",
       ylab = "Survival Fraction",
       col = "red",
       yaxt="n")

  if (!missing(SF)) {
    points(D, SF, col = "red", pch = 19)
  }

  if (missing(SF) || fit_curve) {
    lines(trendlineDs, trendlineSFs, col = "blue", pch = 19)
  }

  ticks <- CoreGx::.GetSupportVec(x=signif(ylim,1), 10)
  labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(round(i, 2)))))
  # labels <- sapply(ticks, function(i) return(sprintf("%0.1e", 10^i)))
  axis(2, at=ticks, labels=labels)
  dev.off()

  return(invisible(0))
}
