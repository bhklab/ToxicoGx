#' @keywords internal
.saveTSetAsGroundTruth <- function(tSet) {
  # Get tSet name for file names
  name <- tSet@annotation$name
  # Save accessor method outputs for use in unit tests
  saveRDS(tSet@cell , file = paste("cellInfo", name, "rds", sep = "."))
  saveRDS(tSet@drug , file = paste("drugInfo", name, "rds", sep = "."))

  ## TODO:: Improve this so that it saves ALL moleclularProfiles
  saveRDS(tSet@molecularProfiles$rna , file = paste("molecularProfiles", name, "rds", sep = "."))
  saveRDS(Biobase::pData(tSet@molecularProfiles$rna) , file = paste("phenoInfo", name, "rds", sep = "."))
  saveRDS(Biobase::fData(tSet@molecularProfiles$rna), file = paste("featureInfo", name, "rds", sep = "."))

  saveRDS(tSet , file = paste("drugInfo", name, "rds", sep = "."))
  saveRDS(tSet , file = paste("drugInfo", name, "rds", sep = "."))
  saveRDS(tSet , file = paste("drugInfo", name, "rds", sep = "."))
  saveRDS(tSet , file = paste("drugInfo", name, "rds", sep = "."))
  saveRDS(tSet , file = paste("drugInfo", name, "rds", sep = "."))
  saveRDS(tSet , file = paste("drugInfo", name, "rds", sep = "."))
  saveRDS(tSet , file = paste("drugInfo", name, "rds", sep = "."))
  saveRDS(tSet , file = paste("drugInfo", name, "rds", sep = "."))
}
