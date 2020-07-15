#' sensNumber Getter
#'
#' Get the number of sensitivity experiments per drug and cell line in a
#'   ToxicoSet
#'
#' @examples
#' sensNumber(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#'
#' @return A \code{data.frame} with the number of sensitivity experiments per
#'   drug and cell line
#'
#' @describeIn ToxicoSet Return the summary of available sensitivity
#'   experiments
#'
#' @importFrom CoreGx sensNumber
#' @importFrom methods callNextMethod
#' @export
setMethod("sensNumber",
          "ToxicoSet",
          function(object)
          {
              callNextMethod(object)
          })

#' sensNumber<- Setter
#'
#' Set the number of sensitivity experiments per drug and cell line in a
#'   ToxicoSet object
#'
#' @examples
#' sensNumber(TGGATESsmall) <- sensNumber(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object to modify
#' @param value A \code{matrix} of replacement values
#'
#' @return The updated \code{ToxicoSet}
#'
#' @describeIn ToxicoSet Update the summary of available sensitivity
#'   experiments
#'
#' @importMethodsFrom CoreGx sensNumber<-
#' @importFrom methods callNextMethod
#' @export
setReplaceMethod('sensNumber',
                 signature = signature(object="ToxicoSet",
                                       value="matrix"),
                 function(object, value)
                 {
                     callNextMethod(object=object, value=value)
                 })