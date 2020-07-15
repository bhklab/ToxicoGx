#' dateCreated Getter
#'
#' Get the date a ToxicoSet object was created
#'
#' @examples
#' dateCreated(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#'
#' @return The date the ToxicoSet was created
#'
#' @describeIn ToxicoSet Return the date the ToxicoSet was created
#'
#' @importFrom CoreGx dateCreated
#' @importFrom methods callNextMethod
#' @export
setMethod("dateCreated",
          signature("ToxicoSet"),
          function(object)
          {
              callNextMethod(object)
          })

##TODO:: Export this to CoreGx

#' datasetType Generic
#'
#' A generic for retrieving the dataset type of an tSet object
#'
#' @param object A \code{ToxicoSet} from which to retrieve the dataset type
#' @param ... A \code{list} containing fall through arguments; this allows
#'   addition of new parameters to methods for this generic
#'
#' @return A \code{character} vector containing the dataset type
#'
#' @export
setGeneric("datasetType", function(object, ...) standardGeneric("datasetType"))

#' datasetType Getter
#'
#' @examples
#' data(TGGATESsmall)
#' datasetType(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} from which to retrieve the dataset type
#' @param ... A \code{list} containing fall through arguments; this allows
#'   addition of new parameters to methods for this generic
#'
#' @describeIn ToxicoSet Update the dataset type of an tSet and return a copy of
#'     the updated object
#' @export
setMethod("datasetType", signature("ToxicoSet"), function(object) {
    ##TODO:: Add error handling to this function
    object@datasetType
})
