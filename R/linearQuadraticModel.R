#' Fit linear-quadratic curves to dose-response data
#'
#' @description This function fits a linear-quadratic curve to dose-response data.
#' @param D vector of radiation doses
#' @param SF vector of survival fractions corresponding to the doses
#' @param lower_bounds vector of length 2 containing minimum allowed values of fitted alpha and beta, respectively
#' @param upper_bounds vector of length 2 containing maximum allowed values of fitted alpha and beta, respectively
#' @param scale parameter of the assumed error distribution of the data; see details
#' @param family family of distributions of the error terms in the data; currently supported options are "normal" and "cauchy"
#' @param median_n see details
#' @param trunc should survival fractions be truncated downward to 1? Defaults to FALSE.
#' @param verbose see details
#' @details 'verbose' outputs warnings that are otherwised suppressed when the function sanity-checks user inputs. 'median_n' denotes the number of distributions from family 'family' that are medianned. (Note that setting n = 1 (the default) is equivalent to using a simple normal or cauchy distribution without taking any medians.)
#' @examples linearQuadraticModel(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10), 
#'  c(1.1, 0.8, 0.7, 0.45, 0.15, -0.1, -0.1, -0.4, -0.65, -0.75, -1.1))
#' @export

linearQuadraticModel <- function (D,
                                  SF,
                                  lower_bounds = c(0, 0),
                                  upper_bounds = c(1, 1),
                                  scale = 5,
                                  family = c("normal", "Cauchy"),
                                  median_n = 1,
                                  # SF_as_log = TRUE,
                                  trunc = FALSE,
                                  verbose = FALSE) {
  match.arg(family)

  CoreGx::.sanitizeInput(x = D,
                          y = SF,
                          x_as_log = FALSE,
                          y_as_log = FALSE,
                          y_as_pct = FALSE,
                          trunc = trunc,
                          verbose = verbose)

  DSF <- CoreGx::.reformatData(x = D,
                               y = SF,
                               x_to_log = FALSE,
                               y_to_log = TRUE,
                               y_to_frac = FALSE,
                               trunc = trunc)
  D <- DSF[["x"]]
  SF <- DSF[["y"]]

  if (!(all(lower_bounds < upper_bounds))) {
    if (verbose == 2) {
      print("lower_bounds:")
      print(lower_bounds)
      print("upper_bounds:")
      print(upper_bounds)
    }
    stop ("All lower bounds must be less than the corresponding upper_bounds.")
  }

  if(!((0 %in% D) || SF[D==0] == 0)){
    D <- c(0,D)
    SF <- c(0,SF)
  }

  gritty_guess <- .makeGrittyGuess(lower_bounds = lower_bounds,
                                   upper_bounds = upper_bounds,
                                   D = D,
                                   SF = SF)

  guess <- CoreGx::.fitCurve(x = D,
                              y = SF,
                              f = .linearQuadratic,
                              density = c(100, 100),
                              step = c(0.005, 0.005),
                              precision = 0.005,
                              lower_bounds = lower_bounds,
                              upper_bounds = upper_bounds,
                              scale = scale,
                              family = family,
                              median_n = median_n,
                              trunc = FALSE,
                              verbose = verbose,
                              gritty_guess = gritty_guess,
                              span = 0.1)

  names(guess) <- c("alpha", "beta")

  return(guess)
}
