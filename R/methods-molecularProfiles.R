#' molecularProfiles Getter
#'
#' Get the molecular profile data associated with the specific molecular data
#'
#' @examples
#' data(TGGATESsmall)
#' TGGATES_mProf <- molecularProfiles(TGGATESsmall, "rna")[seq_len(10),]
#'
#' @param object A \code{ToxicoSet} object.
#' @param mDataType \code{character} A string specifying the type of molecular
#'    data to list the phenotype information for.
#' @param assay A \code{character} Name of the desired assay; if excluded
#' defaults to first assay in the SummarizedExperiment for the given mDataType
#'
#' @describeIn ToxicoSet Return the given type of molecular data from the
#'   ToxicoSet
#'
#' @importFrom CoreGx molecularProfiles
#' @importFrom methods callNextMethod
#' @export
setMethod("molecularProfiles", signature("ToxicoSet"), function(object, mDataType, assay) {
    callNextMethod(object=object, mDataType=mDataType, assay=assay)
})

#' molecularProfiles<- Setter
#'
#' Set the molecular profile data assciated with the specificied molecular
#'   data type
#'
#' @examples
#' molecularProfiles(TGGATESsmall, "rna") <-
#'   molecularProfiles(TGGATESsmall, "rna")
#'
#' @param object A \code{ToxicoSet} object.
#' @param value A \code{matrix} of replacement values.
#' @param assay \code{character} Name or index of the assay data to return
#' @param mDataType \code{character} A string specifying the type of molecular
#'   data to list the phenotype information for.
#'
#' @return Updated \code{ToxicoSet}
#'
#' @importMethodsFrom CoreGx molecularProfiles<-
#'
#' @describeIn ToxicoSet Update the given type of molecular data from the
#'    ToxicoSet
#'
#' @importFrom methods callNextMethod
#' @importFrom CoreGx molecularProfiles<-
#' @export
setReplaceMethod("molecularProfiles",
                 signature = signature(object="ToxicoSet",
                                       mDataType ="character",
                                       assay="character",
                                       value="matrix"),
                 function(object, mDataType, assay, value){
                     callNextMethod(object=object, mDataType=mDataType, assay=assay, value=value)
                 })
#' @describeIn ToxicoSet Update the given type of molecular data from the
#'   ToxicoSet
#' @export
setReplaceMethod("molecularProfiles",
                 signature = signature(object="ToxicoSet",
                                       mDataType ="character",
                                       assay="missing",
                                       value="matrix"),
                 function(object, mDataType, assay, value){
                     callNextMethod(object=object, mDataType=mDataType, assay=assay, value=value)
                 })