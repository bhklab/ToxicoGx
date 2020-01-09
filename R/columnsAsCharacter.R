# Change all factor columns of tSet data.frames to character
#
#
# @importFrom maggitr %<>%
# @importFrom dplyr mutate_if
columnsAsCharacter <- function(tSet) {
  for (mDataType in names(tSet@molecularProfiles)) {
    phenoInfo(tSet, mDataType) %<>%
      tibble::rownames_to_column() %>%
      dplyr::mutate_if(is.factor, as.character) %>%
      tibble::column_to_rownames()

    featureInfo(tSet, mDataType) %<>%
      tibble::rownames_to_column() %>%
      dplyr::mutate_if(is.factor, as.character) %>%
      tibble::column_to_rownames()
  }
  if (length(sensitivityInfo(tSet)) > 0) {
    sensitivityInfo(tSet) %<>%
      tibble::rownames_to_column() %>%
      dplyr::mutate_if(is.factor, as.character) %>%
      tibble::column_to_rownames()
  }
  cellInfo(tSet) %<>%
    tibble::rownames_to_column() %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    tibble::column_to_rownames()

  drugInfo(tSet) %<>%
    tibble::rownames_to_column() %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    tibble::column_to_rownames()

  return(tSet)
}

## TODO:: Determine how to use list.filter to recursively return all data.frames in a nested lsit
## TODO:: Write a recursive functions to return all non-list elements in a nested list
