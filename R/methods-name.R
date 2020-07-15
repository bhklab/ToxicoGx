#' name Getter method
#'
#' Retrieves the name of a tSet
#'
#' @examples
#' name(TGGATESsmall)
#'
#' @param object \code{ToxicoSet} A ToxicoSet object
#'
#' @return \code{character} A string of the tSet's name
#'
#' @importFrom CoreGx name
#' @importFrom methods callNextMethod
#' @export
setMethod(name, "ToxicoSet", function(object) {
    callNextMethod(object)
})

##TODO:: Implement name Setter method and add the generic to CoreGx