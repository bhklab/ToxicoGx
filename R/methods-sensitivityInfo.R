#' sensitivityInfo Getter
#'
#' Get the annotations for the sensitivity experiments in the ToxicoSet
#'
#' @examples
#' data(TGGATESsmall)
#' sensInf<- sensitivityInfo(TGGATESsmall)[seq_len(10),]
#'
#' @param object A \code{ToxicoSet} object
#'
#' @return a \code{data.frame} with the experiment info
#'
#' @describeIn ToxicoSet Return the drug dose sensitivity experiment info
#'
#' @importFrom CoreGx sensitivityInfo
#' @importFrom methods callNextMethod
#' @export
setMethod("sensitivityInfo",
          "ToxicoSet",
          function(object)
          {
              callNextMethod(object)
          })

#' sensitivityInfo<- Setter
#'
#' Set the annotations for sensitivity experiments in this ToxicSet
#'
#' @examples
#' data(TGGATESsmall)
#' sensitivityInfo(TGGATESsmall) <- sensitivityInfo(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#' @param value A \code{data.frame} of replacement values
#'
#' @describeIn ToxicoSet Update the sensitivity experiment info
#'
#' @return Updated \code{ToxicoSet}
#'
#' @importMethodsFrom CoreGx sensitivityInfo<-
#' @export
setReplaceMethod("sensitivityInfo",
                 signature = signature(object="ToxicoSet",
                                       value="data.frame"),
                 function(object, value)
                 {
                     callNextMethod(object=object, value=value)
                 })