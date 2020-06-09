#' pertNumber Getter
#'
#' Get an array of the number of pertubration experiments per drug and cell
#'   line in a ToxicoSet object
#'
#' @examples
#' pertNumber(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#'
#' @return A 3D \code{array} with the number of perturbation experiments per
#'   radiation type and cell line, and data type
#'
#' @describeIn ToxicoSet Return the summary of available perturbation
#'   experiments
#'
#' @importMethodsFrom CoreGx pertNumber
#' @importFrom methods callNextMethod
#' @export
setMethod("pertNumber",
          signature("ToxicoSet"),
          function(object)
          {
              callNextMethod(object)
          })

#' pertNumber<- Setter
#'
#' Set the number of perturbation experiments per drug and cell line and
#'   molecular data type in a ToxicoSet object
#'
#' @examples
#' pertNumber(TGGATESsmall) <- pertNumber(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object to modify
#' @param value An \code{array} of replacement values
#'
#' @return The updated \code{ToxicoSet}
#'
#' @describeIn ToxicoSet Update the summary of available perturbation
#'   experiments
#'
#' @importMethodsFrom CoreGx pertNumber<-
#' @export
setReplaceMethod('pertNumber',
                 signature = signature(object="ToxicoSet",
                                       value="array"),
                 function(object, value)
                 {
                     callNextMethod(object, value)
                 })