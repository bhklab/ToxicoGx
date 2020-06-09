##TODO:: Export this to CoreGx
#' datasetType<- Replacement Generic
#'
#' A generic for updating the dataset type of a ToxicoSet object
#'
#' @examples
#' data(TGGATESsmall)
#' datasetType(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} from which to retrieve the dataset type
#' @param value A \code{character} vector containing the dataset type
#'
#' @return A \code{ToxicoSet} with the datasetType slot updated
#'
#' @export
setGeneric("datasetType<-",  function(object, value) standardGeneric("datasetType<-"))
#' @inheritParams datasetType<-
#' @describeIn ToxicoSet Update the dataset type of an tSet and return a copy of
#'     the updated object
#' @export
setReplaceMethod("datasetType", signature("ToxicoSet"), function(object, value) {
    ##TODO:: Add error handling to this function
    object@datasetType <- value
    object
})