##TODO:: Export to CoreGx
#' curation Slot Getter
#'
#' @param object A \code{ToxicoSet}
#' @param ... A \code{list} to allow definition of new parameters on this generic
#'
#' @return A \code{list} of unique cell and tissue identifiers to check validity
#'   of an tSet
#'
#' @examples
#' data(TGGATESsmall)
#' curation(TGGATESsmall)
#'
#' @export
setGeneric("curation", function(object, ...) standardGeneric("curation"))
#' @describeIn ToxicoSet Retrieve the curation slot form an tSet
#' @inheritParams curation
#' @export
setMethod('curation', signature("ToxicoSet"), function(object) {
    object@curation
})

##TODO:: Export to CoreGx
##FIXME:: How do I import generics from BiocGenerics?
#' curation<- Slot Setter
#'
#' @param object A \code{ToxicoSet}
#' @param ... A \code{list} to allow definition of new parameters on this generic
#' @param value A \code{list} of curations for the cell and tissues types in the
#'   tSet object
#'
#' @return A copy of the \code{ToxicoSet} with the updated curation slot
#'
#' @examples
#' data(TGGATESsmall)
#' curation(TGGATESsmall) <- curation(TGGATESsmall)
#'
#' @export
setGeneric("curation<-", function(object, ..., value) standardGeneric("curation<-"))
#' @describeIn ToxicoSet Update the annotation slot of a tSet
#' @inheritParams annotation<-
#' @export
setReplaceMethod("curation", signature("ToxicoSet", "list"), function(object, value) {
    object@curation <- value
    object
})