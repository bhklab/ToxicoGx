#' sensitivityProfiles Getter
#'
#' Get the data for sensitivty experiments on cell lines in a ToxicoSet
#'
#' @examples
#' data(TGGATESsmall)
#' sensProf <- sensitivityProfiles(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#'
#' @return a \code{data.frame} with the experiment info
#'
#' @describeIn ToxicoSet Return the phenotypic data for the drug dose
#'   sensitivity
#'
#' @importFrom CoreGx sensitivityProfiles
#' @export
setMethod("sensitivityProfiles",
          "ToxicoSet",
          function(object)
          {
              callNextMethod(object)
          })


#' sensitivityProfiles<-
#'
#' Set the data for sensitivity experiments on cell lines in a ToxicoSet
#'
#' @examples
#' sensitivityProfiles(TGGATESsmall) <- sensitivityProfiles(TGGATESsmall)
#'
#' @param object A \code{ToxicoSet} object
#' @param value A \code{data.frame} of replacement values
#'
#' @return Updated \code{ToxicoSet}
#'
#' @importFrom CoreGx sensitivityProfiles<-
#'
#' @describeIn ToxicoSet Update the phenotypic data for the drug dose
#'   sensitivity
#' @export
setReplaceMethod("sensitivityProfiles",
                 signature = signature(object="ToxicoSet",
                                       value="data.frame"),
                 function(object, value)
                 {
                     callNextMethod(object=object, value=value)
                 })
#' @describeIn ToxicoSet Update the phenotypic data for the drug dose
#'   sensitivity
#'
## TODO:: Find out how to document overloaded methods (to include multiple parameter types)
#' @export
setReplaceMethod("sensitivityProfiles",
                 signature = signature(object="ToxicoSet",
                                       value="matrix"),
                 function(object, value)
                 {
                     callNextMethod(object=object, value=value)
                 })