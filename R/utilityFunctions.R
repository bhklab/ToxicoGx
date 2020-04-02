# UTILITY FUNCTIONS

#' Warning if verbose
#'
#' @param warning \code{character} A warning message to display if verbose is True
#'
#' @keywords internal
#' @export
#'
warn_verbose <- function(warning) {
  if (verbose == TRUE) {
    warning(warning)
  }
}

#' Message if verbose
#'
#' @param message \code{character} A message to display if verbose is True
#'
#' @keywords internal
#' @export
#'
msg_verbose <- function(message) {
  if (verbose == TRUE) {
    message(message)
  }
}

