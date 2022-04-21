#' @include ToxicoSet-accessors.R
NULL

#' Update the ToxicoSet class after changes in it struture or API
#'
#' @param object A `ToxicoSet` object to update the class structure for.
#'
#' @return `ToxicoSet` with update class structure.
#'
#' @md
#' @importMethodsFrom CoreGx updateObject
#' @export
setMethod("updateObject", signature("ToxicoSet"), function(object) {
    cSet <- callNextMethod(object)
    tSet <- as(cSet, "ToxicoSet")
    names(curation(tSet)) <- gsub("drug", "treatment", names(curation(tSet)))
    if ("treatment" %in% names(curation(tSet))) {
        colnames(curation(tSet)$treatment) <- gsub("treatmentid", "treatmentid",
            colnames(curation(tSet)$treatment))
    }
    validObject(tSet)
    return(tSet)
})