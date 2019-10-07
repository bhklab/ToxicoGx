#' Generate data to test a tSet object and class accessor methods
#'
#' Unit tests in this package will assume that this data is the ground truth
#'   and any deviation from these files represents a breakdown in the curation
#'   and/or accessor methods of the tSet object
#'
#' @param tSet A \code{ToxicoSet} object to be used as ground truth for tests
#'
.saveTSetAsGroundTruth <- function(tSet) {
  # Get tSet name for file names
  name <- tSet@annotation$name

  ## TODO:: Find the path dynamically
  # Set location of tests/testthat folder
  path <- "tests/testthat/"

  # To test drugInfo() against and validate the slot data is unchanged
  saveRDS(tSet@drug , file = paste0(path, "drugInfo.", name, ".rds"))
  # To test drugNames() against
  saveRDS(tSet@drug$drugid , file = paste0(path, "drugNames.", name, ".rds"))

  # Validate annotation slot is unchanged
  saveRDS(tSet@annotation , file = paste0(path, "annotation.", name, ".rds"))
  # To test tSetName() against
  saveRDS(tSet@annotation$name, file = paste0(path, "name.", name, ".rds"))

  # To test mDataNames() against
  saveRDS(names(tSet@molecularProfiles), file = paste0(path, "mDataNames.", name, ".rds"))
  # Save appropriate data for ALL molecularProfiles in tSet
  parallel::mclapply(names(tSet@molecularProfiles),
                     function(dataType) {
                       # To test molecularProfiles() against
                       saveRDS(Biobase::exprs(tSet@molecularProfiles[[dataType]]) , file = paste0(path, dataType, ".molecularProfiles.", name, ".rds"))
                       # To test phenoInfo() against
                       saveRDS(Biobase::pData(tSet@molecularProfiles[[dataType]]) , file = paste0(path, dataType, ".phenoInfo.", name, ".rds"))
                       # To test fNames() against
                       saveRDS(rownames(Biobase::fData(tSet@molecularProfiles[[dataType]])), file = paste0(path, dataType, ".fNames.", name, ".rds"))
                       # To test featureInfo() against
                       saveRDS(Biobase::fData(tSet@molecularProfiles[[dataType]]), file = paste0(path, dataType, ".featureInfo.", name, ".rds"))
                     })

  # To test cellInfo() against
  saveRDS(tSet@cell , file = paste0(path, "cellInfo.", name, ".rds"))
  # To test cellNames against
  saveRDS(tSet@cell$cellid , file = paste0(path, "cellNames.", name, ".rds"))

  # Validate sensitivity slot data is unchanged
  saveRDS(tSet@sensitivity , file = paste0(path, "sensitivity.", name, ".rds"))
  # To test sesnstivityInfo() against
  saveRDS(tSet@sensitivity$info , file = paste0(path, "sensitivityInfo.", name, ".rds"))
  # To test sensitivityProfiles() against
  saveRDS(tSet@sensitivity$profiles , file = paste0(path, "sensitivityProfiles.", name, ".rds"))
  # To test sensitivityMeasures() against
  saveRDS(names(tSet@sensitivity$profiles) , file = paste0(path, "sensitivityMeasures.", name, ".rds"))
  # To test sensNumber() against
  saveRDS(tSet@sensitivity$n , file = paste0(path, "sensNumber.", name, ".rds"))

  # Validate perturbation slot data is unchanged
  saveRDS(tSet@perturbation , file = paste0(path, "perturbation.", name, ".rds"))
  # To test pertNumber() against
  saveRDS(tSet@perturbation$n , file = paste0(path, "pertNumber.", name, ".rds"))

  # Validate curation slot data is unchanged
  saveRDS(tSet@curation , file = paste0(path, "curation.", name, ".rds"))

  # To test subsetTo against
  saveRDS(subsetTo(tSet, drugs=drugNames(tSet)[1]), file = paste0(path, "subsetTo.", name, ".rds"))

}
