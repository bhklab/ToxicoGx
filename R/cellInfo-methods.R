#' cellInfo Getter
#'
#' Get the cell line annotations in a ToxicoSet
#'
#' @examples
#' data(TGGATESsmall)
#' cellInfo <- cellInfo(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#'
#' @return a \code{data.frame} with the cell annotations
#'
#' @describeIn ToxicoSet Returns the annotations for all the cell lines tested
#'   on in the ToxicoSet
#'
#' @importFrom CoreGx cellInfo
#' @importFrom methods callNextMethod
#' @export
setMethod(cellInfo, "ToxicoSet", function(object){
    callNextMethod(object)
})

#' cellInfo Replace Method
#'
#' Set cell line annotations for a ToxicoSet object
#'
#' @examples
#' data(TGGATESsmall)
#' cellInfo(TGGATESsmall) <- cellInfo(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#' @param value A \code{data.frame} of replacement values
#'
#' @return Updated \code{ToxicoSet}
#'
#' @describeIn ToxicoSet Returns the annotations for all the cell lines tested
#'   on in the ToxicoSet
#'
#' @importFrom CoreGx cellInfo<-
#' @importFrom methods callNextMethod
#' @export
setReplaceMethod("cellInfo",
                 signature = signature(object = "ToxicoSet",
                                       value = "data.frame"),
                 function(object, value)
                 {
                     callNextMethod(object, value)
                 })
