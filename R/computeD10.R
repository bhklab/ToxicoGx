#' Compute D10
#'
#' @description This function computes the radiation dose at which only 10% of cancer cells survive under the exponential model SF = exp(-alpha * D - beta * D ^ 2) given alpha and beta, where D is the radiation dose given and SF is the fraction of cells surviving
#' @param pars parameters (alpha, beta) in equation y = exp(-alpha * x - beta * x ^ 2)
#' @details The units of the returned dose are the inverses of the units of the alpha and beta passed in.
#' @examples computeD10(c(0.2, 0.1))
#' @export

computeD10 <- function(pars) {
  CoreGx::.sanitizeInput(pars = pars,
                          x_as_log = FALSE,
                          y_as_log = FALSE,
                          y_as_pct = FALSE,
                          trunc = FALSE,
                          verbose = FALSE)

  return(.linearQuadraticInv(SF = 0.1,
                             pars = pars,
                             SF_as_log = FALSE))
}
