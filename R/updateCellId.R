### TODO:: Add updating of sensitivity Number tables
#' A function to update cell ids
#' @examples
#' data(TGGATESsmall)
#' updateCellId(TGGATESsmall, new.ids = cellNames(TGGATESsmall))
#'
#' @param tSet A \code{ToxicoSet} object to be updates
#' @param new.ids A \code{character} vector of ids to update with
#'
#' @return \code{none} Updates the cell ids in the ToxicoSet
#'
#' @importFrom CoreGx updateCellId
#' @keywords internal
#' @export
updateCellId <- function(tSet, new.ids = vector("character")){
    CoreGx::updateCellId(tSet, new.ids)
}