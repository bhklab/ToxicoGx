#' molecularProfilesSlot Getter
#'
#' @describeIn ToxicoSet Get contents of molecularProfiles slot
#'
#' @examples
#' data(TGGATESsmall)
#' molecularProfilesSlot(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} from which to return a list of all availble
#'   SummarizedExperiment objects
#'
#' @return A \code{list} containing the molecularProfiles from a tSet
#'
#' @importFrom CoreGx molecularProfilesSlot
#' @importFrom methods callNextMethod
#' @export
setMethod("molecularProfilesSlot", signature("ToxicoSet"), function(object) {
    callNextMethod(object)
})

#' molecularProfilesSlot<- Setter
#'
#' @describeIn ToxicoSet Update the molecular profiles slot of a ToxicoSet and
#'    returns the updated copy
#'
#' @examples
#' data(TGGATESsmall)
#' molecularProfilesSlot(TGGATESsmall) <- molecularProfilesSlot(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object for which values will be replaced
#' @param value A \code{list} containing molecular profiles as SummarizedExperiments
#'
#' @return A copy of the \code{ToxicoSet} with the molecularProfiles slot updated
#'
#' @importFrom CoreGx molecularProfilesSlot<-
#' @importFrom methods callNextMethod
#' @export
setReplaceMethod("molecularProfilesSlot", signature("ToxicoSet"),
                 function(object, value) {
                     callNextMethod(object, value)
                 })