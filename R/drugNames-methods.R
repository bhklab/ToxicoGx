##TODO:: Abstract this method to CoreGx
#' drugNames Generic
#'
#' A generic for the drugNames method
#'
#' @examples
#' data(TGGATESsmall)
#' drugName <- drugNames(TGGATESsmall)[seq_len(10)]
#'
#' @param object A \code{ToxicoSet} object from which to retrieve the included
#'   drug names
#'
#' @return A vector of the drug names used in the ToxicoSet
setGeneric("drugNames", function(object) standardGeneric("drugNames"))
#' @describeIn ToxicoSet Return the names of the drugs used in the ToxicoSet
#' @export
setMethod(drugNames,
          "ToxicoSet",
          function(object)
          {
              rownames(drugInfo(object))
          })

##TODO:: Abstract this method to CoreGx
#' drugNames<- Generic
#'
#' A generic for the drugNames replacement method
#'
#' @examples
#' data(TGGATESsmall)
#' drugNames(TGGATESsmall) <- drugNames(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object to modify
#' @param value A \code{character} vector of replacement drug names
#'
#' @return Updated \code{ToxicoSet}
setGeneric("drugNames<-", function(object, value) standardGeneric("drugNames<-"))
#' @describeIn ToxicoSet Update the drug names used in the dataset
#' @export
setReplaceMethod("drugNames",
                 signature = signature(object="ToxicoSet",
                                       value="character"),
                 function(object, value)
                 {
                     object <- updateDrugId(object, value)
                     return(object)
                 })
