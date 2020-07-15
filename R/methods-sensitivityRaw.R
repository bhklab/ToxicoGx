##TODO:: Migrate this to CoreGx
#' sensitivityRaw Generic
#'
#' @examples
#' data(TGGATESsmall)
#' sensitivityRaw(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} to extract the raw sensitivity data from
#' @param ... A \code{list} to allow new parameters in specific methods
#'
#' @return A \code{array} containing the raw sensitivity data
#'
#' @export
setGeneric("sensitivityRaw", function(object, ...) standardGeneric("sensitivityRaw"))
#' @describeIn ToxicoSet Retrive the raw dose and viability data from an tSet
#' @inheritParams sensitivityRaw
#' @export
setMethod("sensitivityRaw", signature("ToxicoSet"), function(object) {
    object@sensitivity$raw
})

##TODO:: Migrate this to CoreGx
#' sensitivityRaw<- Replacement Generic
#'
#' @examples
#' data(TGATESsmall)
#' sensitivityRaw(TGGATESsmall) <- sensitivityRaw(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} to extract the raw sensitivity data from
#' @param ... A \code{list} to allow new parameters in specific methods
#' @param value A \code{array} containing the raw dose and viability data for
#'   the tSet
#'
#' @return A copy of the \code{ToxicoSet} containing the updated sensitivty data
#'
#' @export
setGeneric("sensitivityRaw<-", function(object, ..., value) standardGeneric("sensitivityRaw<-"))
#' @describeIn ToxicoSet Set the raw dose and viability data for a tSet and return
#'   and updated copty
#' @inheritParams sensitivityRaw<-
#' @export
setReplaceMethod("sensitivityRaw", signature("ToxicoSet", "array"),
                 function(object, value) {
                     object@sensitivity$raw <- value
                     object
                 })