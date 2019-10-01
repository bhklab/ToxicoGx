setOldClass('sessionInfo', sessionInfo)

#' @importFrom utils sessionInfo
.ToxicoSig <- setClass('ToxicoSig', slots=list(
            Arguments = "list",
            tSetName='character',
            DateCreated = 'character',
            SigType = 'character',
            SessionInfo = 'sessionInfo',
            Call = 'character'), contains='array')

ToxicoSig <- function(Data=array(NA, dim=c(0,0,0)), tSetName='', DateCreated=date(), SigType='sensitivity', SessionInfo=sessionInfo(), Call='No Call Recorded', Arguments = list()){
  return(.ToxicoSig(Data, Arguments = Arguments, tSetName=tSetName, DateCreated=DateCreated, SigType=SigType, SessionInfo=SessionInfo, Call=Call))}

#' Show ToxicoGx Signatures
#'
#' @examples
#' data(TGGATESsmall)
#' drug.perturbation <- drugPerturbationSig(TGGATESsmall, mDataType="rna", nthread=1, duration = "2",
#'      drugs = head(drugNames(TGGATESsmall)), features = fNames(TGGATESsmall, "rna")[1:2])
#' drug.perturbation
#'
#' @param object \code{ToxicoSig}
#'
#' @return Prints the ToxicoGx Signatures object to the output stream, and returns invisible NULL.
#'
#' @export
#'
setMethod("show", signature=signature(object='ToxicoSig'),
          function(object) {
            cat('ToxicoSet Name: ', attr(object, 'PSetName'), "\n")
            cat('Signature Type: ', attr(object, 'SigType'), "\n")
            cat("Date Created: ", attr(object, 'DateCreated'), "\n")
            cat("Number of Drugs: ", dim(object)[[2]], "\n")
            cat("Number of Genes/Probes: ", dim(object)[[1]], "\n")
          })

#' Show the Annotations of a signature object
#'
#' This funtion prints out the information about the call used to compute the drug signatures, and the session info
#' for the session in which the computation was done. Useful for determining the exact conditions used to generate signatures.
#'
#' @examples
#' data(TGGATESsmall)
#' drug.perturbation <- drugPerturbationSig(TGGATESsmall, mDataType="rna", nthread=1, duration = "2",
#'      drugs = head(drugNames(TGGATESsmall)), features = fNames(TGGATESsmall, "rna")[1:2])
#' showSigAnnot(drug.perturbation)
#'
#' @param Sigs An object of the \code{ToxicoSig} Class, as returned by \code{drugPerturbationSig}
#'
#' @return Prints the ToxicoGx Signatures annotations to the output stream, and returns invisible NULL.
#'
#' @export
#'
showSigAnnot <- function(Sigs){

  print(Sigs@Call)
  print(Sigs@SessionInfo)
  return(invisible(NULL))
}
