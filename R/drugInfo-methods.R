##TODO:: Abstract this method to CoreGx
#' drugInfo Getter
#'
#' Get the drug annotations in a ToxicoSet object
#'
#' @examples
#' data(TGGATESsmall)
#' drugInfo <- drugInfo(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#'
#' @return a \code{data.frame} with the drug annotations
#'
setGeneric("drugInfo", function(object) standardGeneric("drugInfo"))
#' @describeIn ToxicoSet Returns the annotations for all the drugs tested in
#'   the ToxicoSet
#' @export
setMethod("drugInfo", signature("ToxicoSet"), function(object) {
    object@drug
})


#' drugInfo<- Setter method
#'
#' Set the drug annotations in a ToxicoSet object
#'
#' @examples
#' data(TGGATESsmall)
#' drugInfo(TGGATESsmall) <- drugInfo(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object.
#' @param value A \code{data.frame} of replacement values.
#'
#' @return Updated \code{ToxicoSet}
#'
setGeneric("drugInfo<-", function(object, value) standardGeneric("drugInfo<-"))
#' @describeIn ToxicoSet Update the drug annotations
#' @export
setReplaceMethod("drugInfo", signature = signature(object = "ToxicoSet",value = "data.frame"), function(object, value){
    object@drug <- value
    object
})
