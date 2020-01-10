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
paramMissingHandler <- function(funName, tSet, mDataType, ...) {

  ## Errors if tSet parameter not passed an argument
  if (missing(tSet)) {
    stop(paste0(funName, " requires a tSet argument!"))
  }

  ## Errors if mDataType parameter not passed an argument
  if (missing(mDataType)) {
    if (funName == "summarizeSensitivityProfiles" | funName == "subsetTo") {
      mDataType <- names(tSet@molecularProfiles)
    } else {
      stop(paste0(funName, " requires an mDataType argument!"))
    }
  }

  ## Interection of missing values for similar functions
  intersectMissingChecks <- c(
    "cell_lines", "drugs"
  )

  missingChecks <-
    switch(funName,
           "summarizeMolecularProfiles" =
             c(intersectMissingChecks,
               "features", "durations"
             ),
           "summarizeSensitivityProfiles" =
             c(intersectMissingChecks,
              "duration"
              ),
           "drugPerturbationSig" =
             c(intersectMissingChecks,
               "features", "durations", "dose"
               ),
           "subsetTo" =
             c(intersectMissingChecks,
              "features", "durations", "dose"
             ),
    )

  # Assigns values for missing parameters and throws warnings
  return(
    .checkParamsForMissing(funName = funName, tSet = tSet, mDataType = mDataType,
                         missingChecks = missingChecks, ...)
  )
}

#' @keywords internal
.checkParamsForMissing <- function(
  funName = funName, tSet = tSet, missingChecks, mDataType, ...) {
  # Initialize variable names in the local environment
  cell_lines <- duration <- drugs <- features <- dose <- NULL
  # Extract named arguments into local environment
  argList <- list(...)
  for (idx in seq_len(length(argList))) { ## TODO:: Make this work with seq_along()
    assign(names(argList)[idx], argList[[idx]])
  }

  # Prealocate memory for returned list, this speeds up execution as R is bad
  # at dynamic memory allocation
  missingParamValues <- list()

  for (missing in missingChecks) {
    switch(
      missing,
      "cell_lines" = {if (is.null(cell_lines)) { missingParamValues[[missing]] <- unique(cellNames(tSet));
        warning(paste0(missing, " parameter not specified, defaults to all cell lines in the given tSet!"))}
        },
      "drugs" = {if (is.null(drugs)) { missingParamValues[[missing]] <- unique(drugNames(tSet));
        warning(paste0(missing, " parameter not specified, defaults to all drugs in the given tSet!"))}
        },
      "features" = {if (is.null(features)) {missingParamValues[[missing]] <- unique(fNames(tSet, mDataType[1]));
        warning(paste0(missing, " parameter not specified, defaults to all features in the given tSet for the specified mDataType!"))}
        },
      "durations" = {if (is.null(duration)) {missingParamValues[[missing]] <- unique(as.character(ToxicoGx::sensitivityInfo(tSet)$duration_h));
      warning(paste0(missing, " parameter not specified, defaults to all experimental durations in given tSet!"))}
        },
      "duration" = {if (is.null(duration)) {missingParamValues[[missing]] <- unique(as.character(ToxicoGx::sensitivityInfo(tSet)$duration_h))[1];
      warning(paste0(missing, " parameter not specified, defaults to ", missingParamValues[[missing]]))}
      },
      "dose" = {if (is.null(dose)) {missingParamValues[[missing]] <- unique(phenoInfo(tSet, mDataType)$dose_level);
      warning(paste0(missing, " parameter not specified, defaults to all dose levels in the given tSet for the specified mDataType!"))}}
    )
  }
  return(missingParamValues)
}

