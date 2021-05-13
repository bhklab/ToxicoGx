#' Generic method for performing differential expression analysis on an S4 object
#'   using the limma package
#'
#' @param object [`S4`] An S4 object to conduct differential expression analysis
#'   on.
#' @param ... Allow new parameters to be added to this generic.
#'
#' @return To be defined by the method implementation.
#'
#' @export
setGeneric('computeLimmaDiffExpr',
    function(object, ...) setGeneric('methods-computeLimmaDiffExpr.R'))