#' sensitivityMeasures Getter
#'
#' Get the avilable measurments for sensitivity experiments in a ToxicoSet
#'
#' @examples
#' sensitivityMeasures(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#'
#' @return A \code{character} vector of all the available sensitivity measures
#'
#' @describeIn ToxicoSet Returns the available sensitivity profile
#'   summaries, for example, whether there are IC50 values available
#'
#' @importFrom CoreGx sensitivityMeasures
#' @importFrom methods callNextMethod
#' @export
setMethod("sensitivityMeasures",
          "ToxicoSet",
          function(object)
          {
              callNextMethod(object)
          })