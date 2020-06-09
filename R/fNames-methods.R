#' fNames Getter
#'
#' Get the feature names in a ToxicoSet object for the specified molecular data
#'   type
#'
#' @examples
#' fNames(TGGATESsmall, "rna")[seq_len(10)]
#'
#' @param object A \code{ToxicoSet} object
#' @param mDataType \code{character} A string specifying the type of molecular
#'   data to list the phenotype information for
#'
#' @return A \code{character} vector of the feature names
#'
#' @describeIn ToxicoSet Return the feature names used in the dataset
#'
#' @importFrom CoreGx fNames
#' @importFrom methods callNextMethod
#' @export
setMethod("fNames",
          signature("ToxicoSet", "character"),
          function(object, mDataType)
          {
              callNextMethod(object, mDataType)
          })

#' fNames<- Setter
#'
#' Set the feature names in a ToxicoSet object for the specified molecular
#'   data type
#'
#'@examples
#' data(TGGATESsmall)
#' cellNames(TGGATESsmall) <- cellNames(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#' @param mDataType \code{character} A string specifying the type of molecular
#'   data to list the phenotype information for.
#' @param value A \code{character} vector of the new feature names
#'
#' @return Updated \code{ToxicoSet}
#'
#' @describeIn ToxicoSet Update the feature names used in the dataset
#'
#' @importFrom CoreGx fNames<-
#' @importFrom methods callNextMethod
#' @export
setReplaceMethod("fNames",
                 signature("ToxicoSet", "character"),
                 function(object, mDataType, value)
                 {
                     callNextMethod(object, mDataType, value)
                 })