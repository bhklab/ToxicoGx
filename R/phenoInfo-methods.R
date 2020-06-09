#' phenoInfo Getter
#'
#' Get the phenotype annotations for cell lines with the specificed molecular
#'   data type
#'
#' @examples
#' data(TGGATESsmall)
#' phenoInfo <- phenoInfo(TGGATESsmall, mDataType="rna")
#'
#' @param object A \code{ToxicoSet} object
#' @param mDataType \code{character} A string specifying the type of molecular
#'   data to list the phenotype information for.
#'
#' @return a \code{Dframe} with the experiment info
#'
#' @describeIn ToxicoSet Return the experiment info from the given type of
#'   molecular data in ToxicoSet
#'
#' @importFrom CoreGx phenoInfo
#'
#' @export
setMethod("phenoInfo",
          signature("ToxicoSet", "character"),
          function(object, mDataType)
          {
              callNextMethod(object, mDataType)
          })

#' phenoInfo<- Setter
#'
#' Set the phenotype annotations for cell lines with the selected molecular
#'   data type.
#'
#' @examples
#' data(TGGATESsmall)
#' phenoInfo(TGGATESsmall, mDataType="rna") <-
#'   phenoInfo(TGGATESsmall, mDataType="rna")
#'
#' @param object A \code{ToxicoSet} object.
#' @param mDataType A \code{character} with the type of molecular data to return/update
#' @param value A \code{data.frame}, \code{DataFrame} or \code{DFrame} of
#'   replacement values.
#'
#' @return The updated \code{ToxicoSet}
#'
#' @describeIn ToxicoSet Update the the given type of molecular data experiment
#'   info in the ToxicoSet
#'
#' @importMethodsFrom CoreGx phenoInfo<-
#' @importFrom methods callNextMethod
#' @export
setReplaceMethod("phenoInfo",
                 signature = signature(object="ToxicoSet",
                                       mDataType ="character",
                                       value="data.frame"),
                 function(object, mDataType, value)
                 {
                     callNextMethod(object, mDataType, value)
                 })
