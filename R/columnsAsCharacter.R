# Change all factor columns of tSet data.frames to character
#
#
# @importFrom maggitr %<>%
# @importFrom dplyr mutate_if
columnsAsCharacter <- function(tSet) {
  for (mDataType in names(tSet@molecularProfiles)) {
    phenoInfo(tSet, mDataType) %<>% dplyr::mutate_if(is.factor, as.character)
    featureInfo(tSet, mDataType) %<>% dplyr::mutate_if(is.factor, as.character)
  }
  if (length(sensitivityInfo(tSet)) > 0) {
    sensitivityInfo(tSet) %<>% dplyr::mutate_if(is.factor, as.character)
  }
  cellInfo(tSet) %<>% dplyr::mutate_if(is.factor, as.character)
  drugInfo(tSet) %<>% dplyrmutate_if(is.factor, as.character)
  return(tSet)
}

## TODO:: Determine how to use list.filter to recursively return all data.frames in a nested lsit
## TODO:: Write a recursive functions to return all non-list elements in a nested list
