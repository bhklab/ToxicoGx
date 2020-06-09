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