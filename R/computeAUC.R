#' computeAUC: computes AUC
#'
#' @description This function computes the area under a dose-response curve of the form survival fraction SF = exp(-alpha * D - beta * D ^ 2).
#'
#' @examples computeAUC(D=c(0.1, 0.5, 0.7, 0.9), pars=c(0.2, 0.1), lower = 0, upper = 1) # Returns 0.7039296
#'
#' @param D vector of dosages
#' @param SF vector of survival fractions
#' @param pars parameters (alpha, beta) in equation y = exp(-alpha * x - beta * x ^ 2)
#' @param lower lower bound of dose region to compute AUC over
#' @param upper upper bound of dose region to compute AUC over
#' @param trunc should survival fractions be truncated downward to 1 if they exceed 1?
#' @param SF_as_log A boolean indicating whether survival fraction is displayed on a log axis. Defaults to FALSE
#' @param area.type should the AUC of the raw (D, SF) points be returned, or should the AUC of a curve fit to said points be returned instead?
#' @param verbose how detailed should error and warning messages be? See details.
#'
#' @details If lower and/or upper are missing, the function assumes their values
#'   to be the minimum and maximum D-values, respectively. For all warnings to
#'   be silent, set trunc = FALSE. For warnings to be output, set trunc = TRUE.
#'   For warnings to be output along with the arguments that triggered them,
#'   set trunc = 2.
#' @export
#' @importFrom stats pnorm
#' @importFrom caTools trapz

# Added SF_as_log arguement with default as false to match condition on line 93
computeAUC <- function(D, SF, pars, lower, upper, trunc = TRUE, SF_as_log = FALSE, area.type = c("Fitted", "Actual"), verbose = TRUE) {
  area.type <- match.arg(area.type)

  if (!missing(SF)) {
    CoreGx::.sanitizeInput(x = D,
                            y = SF,
                            x_as_log = FALSE,
                            y_as_log = FALSE,
                            y_as_pct = FALSE,
                            trunc = trunc,
                            verbose = FALSE)

    DSF <- CoreGx::.reformatData(x = D,
                                 y = SF,
                                 x_to_log = FALSE,
                                 y_to_log = FALSE,
                                 y_to_frac = FALSE,
                                 trunc = trunc)
    D <- DSF[["x"]]
    SF <- DSF[["y"]]
  } else if (!missing(pars)) {
    CoreGx::.sanitizeInput(pars = pars,
                            x_as_log = FALSE,
                            y_as_log = FALSE,
                            y_as_pct = FALSE,
                            trunc = trunc,
                            verbose = FALSE)
    Dpars <- CoreGx::.reformatData(x = D,
                                    pars = pars,
                                    x_to_log = FALSE,
                                    y_to_log = FALSE,
                                    y_to_frac = FALSE,
                                    trunc = trunc)
    D <- Dpars[["x"]]
    pars <- Dpars[["pars"]]
  } else {
    stop("SF and pars can't both be missing.")
  }

  if (!missing(lower) && !missing(upper)) {
    ###TODO:: Check if this function still works correctly
    CoreGx::.sanitizeInput(pars = pars, # Added this line to resolve error returned from CoreGx
                            lower = lower,
                            upper = upper,
                            x_as_log = FALSE,
                            y_as_log = FALSE,
                            y_as_pct = FALSE,
                            trunc = trunc,
                            verbose = verbose)
  }

  if (area.type == "Fitted") {
    if (missing(pars)) {
      pars <- unlist(linearQuadraticModel(D = D,
                                          SF = SF,
                                          trunc = trunc,
                                          verbose = verbose))
    }
    if (missing(lower)) {
      lower <- min(D)
    }
    if (missing(upper)) {
      upper <- max(D)
    }

    if (SF_as_log == TRUE) { # Modified condition to correct error
      return(pars[[1]] / 2 * (lower ^ 2 - upper ^ 2) + pars[[2]] / 3 * (lower ^ 3 - upper ^ 3))
    } else {
      if (pars[[2]] == 0) {
        if (pars[[1]] == 0) {
          return(upper - lower)
        } else {
          return((exp(-pars[[1]] * lower) - exp(-pars[[1]] * upper)) / pars[[1]])
        }
      } else {
        # return(exp(pars[[1]] ^ 2 / 4 / pars[[2]]) *
        #        sqrt(pi / pars[[2]]) *
        #        (pnorm(sqrt(2 * pars[[2]]) * (upper + pars[[1]] / 2 / pars[[2]])) -
        #         pnorm(sqrt(2 * pars[[2]]) * (lower + pars[[1]] / 2 / pars[[2]]))))
        # return(sqrt(pi / pars[[2]]) *
        #       (exp(pars[[1]] ^ 2 / 4 / pars[[2]] + pnorm(sqrt(2 * pars[[2]]) * (upper + pars[[1]] / 2 / pars[[2]]), log.p = TRUE))
        #        -
        #        exp(pars[[1]] ^ 2 / 4 / pars[[2]] + pnorm(sqrt(2 * pars[[2]]) * (lower + pars[[1]] / 2 / pars[[2]]), log.p = TRUE))))
        x <- CoreGx::.GetSupportVec(x=D, output_length = 1000)
        y <- .linearQuadratic(D=x, pars=pars, SF_as_log=FALSE)
        return(caTools::trapz(x, y))

      }
    }

  } else if (area.type == "Actual") {
    #print("Actual")
    if (missing(SF)) {
      stop("Please pass in SF-values.")
    } else {
      return(caTools::trapz(x = D, y = SF))
    }
  }
}
