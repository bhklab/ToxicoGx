##TODO:: Migrate this to CoreGx
#' sensitivitySlot Generic
#'
#' @param object A \code{ToxicoSet} to extract the raw sensitivity data from
#' @param ... Allow new parameters to be defined for this generic
#'
#' @return A \code{list} of the sensitivity slot contents
#'
#' @export
setGeneric("sensitivitySlot", function(object, ...) standardGeneric("sensitivitySlot"))

#' sensitivitySlot Getter
#'
#' @describeIn ToxicoSet Retrieves the contents of the sensitivity slot
#'
#' @examples
#' data(TGGATESsmall)
#' sensitivitySlot(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} to extract the raw sensitivity data from
#'
#'
#' @export
setMethod("sensitivitySlot", signature("ToxicoSet"), function(object) {
    object@sensitivity
})

##TODO:: Migrate this to CoreGx
#' sensitivitySlot<- Replacement Generic
#'
#' @param object A \code{ToxicoSet} to extract the raw sensitivity data from
#' @param ... Allow new parameters to be defined for this generic
#' @param value A \code{list} of new sensitivity slot data for the tSet
#'
#' @return A copy of the \code{ToxicoSet} containing the updated sensitivty slot
#'
#' @export
setGeneric("sensitivitySlot<-", function(object, ..., value) standardGeneric("sensitivitySlot<-"))

#' sensitivity Slot Setter
#'
#' @describeIn ToxicoSet Set the raw dose and viability data for an tSet and return
#'   and updated copy
#'
#' @examples
#' data(TGGATESsmall)
#' sensitivitySlot(TGGATESsmall) <- sensitivitySlot(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} to extract the raw sensitivity data from
#' @param ... Allow new parameters to be defined for this generic
#' @param value A \code{list} of new sensitivity slot data for the tSet
#'
#' @export
setReplaceMethod("sensitivitySlot", signature(object="ToxicoSet", value="list"),
                 function(object, ..., value) {
                     ##TODO:: Implement error handinlg for this slot
                     object@sensitivity <- value
                     object
                 })
