setOldClass('sessionInfo', sessionInfo)

#' @importFrom utils sessionInfo
.RadioSig <- setClass('RadioSig', slots=list(

            RSetName='character',
            DateCreated = 'character',
            SigType = 'character',
            SessionInfo = 'sessionInfo',
            Call = 'character'), contains='array')

RadioSig <- function(Data=array(NA, dim=c(0,0,0)), PSetName='', DateCreated=date(), SigType='sensitivity', SessionInfo=sessionInfo(), Call='No Call Recorded'){

#attr(SessionInfo, 'class') <- NULL

return(.RadioSig(Data, PSetName=PSetName, DateCreated=DateCreated, SigType=SigType, SessionInfo=SessionInfo, Call=Call))}


#' Show RadioGx Signatures
#'
#' @examples
#' data(Cleveland_small)
#' rad.sensitivity <- radSensitivitySig(Cleveland_small, mDataType="rna",
#'              nthread=1, features = fNames(Cleveland_small, "rna")[1])
#' rad.sensitivity
#'
#' @param object \code{RadioSig}
#' @return Prints the RadioGx Signatures object to the output stream, and returns invisible NULL.
#' @export
setMethod("show", signature=signature(object='RadioSig'),
        function(object) {
        cat('RadioSet Name: ', attr(object, 'RSetName'), "\n")
        cat('Signature Type: ', attr(object, 'SigType'), "\n")
        cat("Date Created: ", attr(object, 'DateCreated'), "\n")
        cat("Number of Radiation Types: ", dim(object)[[2]], "\n")
        cat("Number of Genes/Probes: ", dim(object)[[1]], "\n")
           })

#' Show the Annotations of a signature object
#'
#' This funtion prints out the information about the call used to compute the rad signatures, and the session info
#' for the session in which the computation was done. Useful for determining the exact conditions used to generate signatures.
#'
#' @examples
#' data(Cleveland_small)
#' rad.sensitivity <- radSensitivitySig(Cleveland_small, mDataType="rna",
#'              nthread=1, features = fNames(Cleveland_small, "rna")[1])
#' showSigAnnot(rad.sensitivity)
#'
#' @param Sigs An object of the \code{RadioSig} Class, as
#' returned by \code{radPerturbationSig} or \code{radSensitivitySig}
#' @return Prints the RadioGx Signatures annotations to the output stream, and returns invisible NULL.
#' @export
showSigAnnot <- function(Sigs){

  print(Sigs@Call)
  print(Sigs@SessionInfo)
  return(invisible(NULL))
}





