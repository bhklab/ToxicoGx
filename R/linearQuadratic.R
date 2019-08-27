#' Linear quadratic
#'
#' @param D A vector of drug concentrations
#' @param pars Parameters (a, b) of the linear model
#' @param SF_as_log Boolen indicating whether survival fraction is logged
#'
#' @export
.linearQuadratic <- function(D, pars, SF_as_log = TRUE) {

  SF <- -(pars[[1]] * D + pars[[2]] * D ^ 2)

  if (!SF_as_log) {
    SF <- exp(SF)
  }

  return(SF)
}
