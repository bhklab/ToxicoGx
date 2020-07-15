#' cellNames Getter
#'
#' Get names of cell lines in a ToxicoSet object
#'
#' @examples
#' cellNames(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#'
#' @return A vector of the cell names used in the ToxicoSet
#'
#' @describeIn ToxicoSet Return the cell names used in the dataset
#'
#' @importFrom CoreGx cellNames
#' @importFrom methods callNextMethod
#' @export
setMethod("cellNames",
          signature("ToxicoSet"),
          function(object)
          {
              callNextMethod(object)
          })

#' cellNames<- Setter
#'
#' Set the cell line names in a TocicoSet object
#'
#' @examples
#' data(TGGATESsmall)
#' cellNames(TGGATESsmall) <- cellNames(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object to modify
#' @param value A \code{character} of replacement cell line names
#'
#' @return Updated \code{ToxicoSet}
#'
#' @describeIn ToxicoSet Update the cell names used in the dataset
#'
#' @importMethodsFrom CoreGx cellNames<-
#' @export
setReplaceMethod("cellNames",
                 signature = signature(object="ToxicoSet",
                                       value="character"),
                 function(object, value)
                 {
                     callNextMethod(object=object, value=value)
                 })