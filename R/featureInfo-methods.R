#' Getter for featureInfo method
#'
#' @examples
#' data(TGGATESsmall)
#' featureInfo <- featureInfo(TGGATESsmall, "rna")[seq_len(10),]
#'
#' @param object A \code{ToxicoSet} object
#' @param mDataType \code{character} A string specifying the type of molecular
#'   data to list the phenotype information for.
#'
#' @describeIn ToxicoSet Return the feature info for the given molecular data
#'
#' @importFrom CoreGx featureInfo
#' @importFrom methods callNextMethod
#' @export
setMethod("featureInfo",
          signature("ToxicoSet", "character"),
          function(object, mDataType)
          {
              callNextMethod(object, mDataType)
          })

#' featureInfo<- Setter
#'
#' Set the feature annotations for a specficied molecular data type
#'
#' @examples
#' data(TGGATESsmall)
#' featureInfo(TGGATESsmall, "rna") <- featureInfo(TGGATESsmall, "rna")
#'
#' @param object A \code{ToxicoSet} object
#' @param value A \code{data.frame} of replacement values
#' @param mDataType \code{character} A string specifying the type of molecular
#'   data
#'
#' @return Updated \code{ToxicoSet}
#'
#' @describeIn ToxicoSet Replace the gene info for the molecular data
#'
#' @importMethodsFrom CoreGx featureInfo<-
#' @importFrom methods callNextMethod
#' @export
setReplaceMethod("featureInfo",
                 signature = signature(object="ToxicoSet",
                                       mDataType ="character",
                                       value="data.frame"),
                 function(object, mDataType, value)
                 {
                     callNextMethod(object, mDataType, value)
                 })
