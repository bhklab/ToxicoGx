# A Handler to Assign Default Values for Missing Parameters
#
# This function will take in the params of a function as well as its name.
#   Missing values will then be assigned to a list, which will be used to
#   populate the parent functions scope with the correct default argument for
#   each missing parameter
#
# @param funName \code{character} A string of the function name. This argument is
#   used to match the correct parameter checking conditions with each functionss
# @param ... \code{pairlist} A list of all parameters passed as arguements to the
#   function "funName".
#
# @return \code{list} A list of all missing parameter argument values, named
#    with the respective missing parameters,
#
#' @keywords internal
paramMissingHandlerTests <- function(funName, tSet, mDataType, ...) {

  ## Errors if tSet parameter not passed an argument
  if (missing(tSet)) {
    stop(paste0(funName, " requires a tSet argument!"))
  } else if (missing(mDataType)) {
    stop(paste0(funName, " requires an mDataType argument!"))
  }

  ## Interection of missing values for similar functions
  intersectMissingChecks <- c(
    "cell_lines", "drugs", "features", "duration"
  )

  missingChecks <-
    switch(funName,
           "summarizeMolecularProfiles" =
             #c(
             intersectMissingChecks#,
           #),
    )

  # Assigns values for missing parameters and throws messages
  .checkParamsForMissingTests(
    funName = funName, tSet = tSet, mDataType = mDataType,
    missingChecks = missingChecks, ...
    )
}

#' @keywords internal
.checkParamsForMissingTests <- function(
  funName = funName, tSet = tSet, missingChecks, mDataType, ...) {
  # Initialize variable names in the local environment
  cell_lines <- duration <- drugs <- features <- NULL
  # Extract named arguments into local environment
  argList <- list(...)
  for (idx in seq_len(length(argList))) { ## TODO:: Make this work with seq_along()
    assign(names(argList)[idx], argList[[idx]])
  }

  testthat::context(
    paste("Testing paramMissingHandler returns correct messages for ")
          )

  for (missing in missingChecks) {
    switch(
      missing,
      "cell_lines" = {
      message(paste0(missing, " parameter not specified, defaults to all cell lines in the given tSet!"))
      },
      "drugs" = {if (is.null(drugs)) { missingParamValues[[missing]] <- unique(drugNames(tSet));
      message(paste0(missing, " parameter not specified, defaults to all drugs in the given tSet!"))}
      },
      "features" = {if (is.null(features)) {missingParamValues[[missing]] <- unique(fNames(tSet, mDataType));
      message(paste0(missing, " parameter not specified, defaults to all features in the given tSet for the specified mDataType!"))}
      },
      "duration" = {if (is.null(duration)) {missingParamValues[[missing]] <- unique(as.character(ToxicoGx::sensitivityInfo(tSet)$duration_h));
      message(paste0(missing, " parameter not specified, defaults to all experimental durations in given tSet!"))}
      }
    )
  }
}
