##TODO:: Export to CoreGx
##FIXME:: How do I import generics from BiocGenerics?
#' annotation Slot Getter
#'
#' @param object A \code{ToxicoSet}
#' @param ... A \code{list} to allow definition of new parameters on this generic
#'
#' @return A \code{list} of named annotaiton
#'
#' @examples
#' data(TGGATESsmall)
#' annotation(TGGATESsmall)
#'
#' @export
setGeneric("annotation", function(object, ...) standardGeneric("annotation"))
#' @describeIn ToxicoSet Retrieve the annotations slot form an tSet
#' @inheritParams annotation<-
#' @export
setMethod('annotation', signature("ToxicoSet"), function(object) {
    object@annotation
})

##TODO:: Export to CoreGx
##FIXME:: How do I import generics from BiocGenerics?
#' annotation<- Slot Setter
#'
#' @param object A \code{ToxicoSet}
#' @param ... A \code{list} to allow definition of new parameters on this generic
#' @param value A \code{list} of annotations to add to the annotatiosn slot of
#'   an tSet
#'
#' @return A copy of the \code{ToxicoSet} with the updated annotation slot
#'
#' @examples
#' data(TGGATESsmall)
#' annotation(TGGATESsmall) <- annotation(TGGATESsmall)
#'
#' @export
setGeneric("annotation<-", function(object, ..., value) standardGeneric("annotation<-"))
#' @describeIn ToxicoSet Update the annotation slot of a tSet
#' @inheritParams annotation<-
#' @export
setReplaceMethod("annotation", signature("ToxicoSet", "list"), function(object, value) {
    object@annotation <- value
    object
})
