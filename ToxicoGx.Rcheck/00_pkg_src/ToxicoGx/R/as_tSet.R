#' Coerce pSet to tSet
#'
#' Forces a \code{PharmacoSet} class objects to be a \code{ToxicoSet} object.
#'
#' @param pSet A \code{PharmacoSet} class object.
#'
#' @return A \code{ToxicoSet} class object containing all data from the pSet.
#'
#'
#' @export
pharmacoSettoToxicoSet <-
          function(pSet) {
            tSet <- ToxicoSet(
              name = pSet@annotation$name,
              molecularProfiles = pSet@molecularProfiles,
              cell = pSet@cell,
              drug = pSet@drug,
              sensitivityInfo = pSet@sensitivity$info,
              sensitivityRaw = pSet@sensitivity$raw,
              sensitivityProfiles = pSet@sensitivity$profiles,
              sensitivityN = pSet@sensitivity$n,
              perturbationN = pSet@perturbation$n,
              curationDrug = pSet@curation$drug,
              curationCell = pSet@curation$cell,
              curationTissue = pSet@curation$tissue,
              datasetType = pSet@datasetType,
              verify = TRUE)
            return(tSet)
          }
