#' mDataNames
#'
#' Returns the names of the molecular data types available in a ToxicoSet
#'   object
#'
#' @examples
#' mDataNames(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#'
#' @return Vector of names of the molecular data types
#'
#' @importFrom CoreGx mDataNames cellInfo<-
#' @export
setMethod(
    "mDataNames",
    signature("ToxicoSet"),
    function(object)
    {
        callNextMethod(object)
    })