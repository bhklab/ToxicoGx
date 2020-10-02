## FIXME:: Remove this and define as methods on coerce generic instead
#' Define an S4 Generic for the methods::as function
#'
#' This will allow creation of new definitons for object conversions
#'
#' @param object An object to coerce to a different class
#' @param ... To allow new parameters for this generc
#' @param value The type to coerce the object to
#'
#' @return The object updated to the new type
#'
#' @export
setGeneric('as', function(object, ..., value) methods::as(object, ..., value))

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