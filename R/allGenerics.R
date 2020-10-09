### FIXME:: Move this into CoreGx
##' Private getter for .intern environment slot inside an object.
##'
##' @param object What to object to retrieve a value from .intern for
##' @parma ... Allow definition of new parmaters on this generic
##'
##' @keywords internal
##' @exprot
##' @noRd
#setGeneric('.getIntern', function(object, ...)
#    standardGeneric('.getIntern'))

## TODO:: Move this to CoreGx
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