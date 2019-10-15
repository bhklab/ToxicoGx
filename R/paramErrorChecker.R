
# An Error Checker for Function Parameters
#
# This function will take in the params of a function as well as its name.
#   Error checking will then be conducted on each of it's parameter arguments
#   to ensure they meet the input requirements for that function. Descriptive
#   errors are returned if the the arguements do not meet the criteria for that
#   function.
#
# @param funName [character] A string of the function name. This argument is
#   used to match the correct parameter checking conditions with each function.
# @param ... [pairlist] A list of all parameters passed as arguements to the
#   function "funName".
#
# @return Returns nothing, this function works by side effects only
#
#' @keywords internal
paramErrorChecker <- function(funName, tSet, ...) {

  # Intersection of all function's parameter tests
  intersectParamChecks <- c("tSetNotIs", "tSetGt1",
                            "mDataTypeNotChar", "mDataTypeNotIn",
                            "cell.linesNotChar","cell.linesNotIn",
                            "drugsNotChar",
                            "durationNotChar"
  )
  intersectViabPlotParamChecks <- c("tSetNotIs",
                                    "vaibilitiesNotMissing", "viabilitiesNotNum"
  )

  # Matches the correct parameter constraints to each function name
  paramChecks <-
    switch(funName,
           "drugPerturbationSig" =
             c(intersectParamChecks,
               "mDataTypeGt1",
               "featuresLt2",
               "doseLt2", "doseNotCtl"
             ),
           "summarizeMolecularProfiles" =
             c(intersectParamChecks,
               "summary.statNotChar", "summary.statNotIn", "summary.statGt1"
             ),
           "drugDoseResponseCurve" =
             c(intersectViabPlotParamChecks, "viabilitiesDiffLenConc",
               "concentrationsNotNum"
             ),
           "drugTimeResponseCurve" =
             c(intersectViabPlotParamChecks, "viabilitiesDiffLenDur",
               "durationNotChar")
    )

  ## Handle missing values
  if (missing(tSet)) {
    tSet <- NULL
  }

  ## DEBUG STATEMENTS
  #print(paste0("Function: ", funName))
  #print(paste0("Checks: ",paste(paramChecks, collapse=",")))

  # Conducts the parameter checks based on function matched paramChecks
  .checkParamsForErrors(tSet = tSet, paramChecks = paramChecks, ...)
}

#' @keywords internal
.checkParamsForErrors <- function(tSet, paramChecks, ...) {

  cell.lines = concentrations = dose = drugs = duration = features = mDataType =
    summary.stat = tSets = viabilities <- NULL

  # Extract named arguments into local environment
  argList <- list(...)
  for (idx in seq_len(length(argList))) { ## TODO:: Make this work with seq_along()
    assign(names(argList)[idx], argList[[idx]])
  }

  ## TODO:: Write a cases function that lazily evaluates LHS to replace this switch
  ## TODO:: Benmark for loop vs apply statement for this code
  # Runs the parameter checks specific for the given funName
  for (check in paramChecks) {
    invisible(
      switch(
        check,
        # tSet checks
        "tSetNotIs" = {if (!is(tSet, "ToxicoSet")) { stop(paste0(tSet@annotations$name, " is a ", class(tSet), ", not a ToxicoSet.")) }},
        "tSetGt1" = {if (length(tSet) > 1) { stop("You may only pass in one tSet.") }},
        # tSets checks
        "tSetsNotIs" = {if (!is.null(tSets)) { if (!all(vapply(tSets, function (tSet) { is(tSet, "ToxicoSet") }, FUN.VALUE = logical(1) ))) { stop("One or more arguments to tSets parameter is not a 'ToxicoSet'.")}}},
        # mDataType checks
        "mDataTypeGt1" = {if (length(unlist(mDataType)) > 1) { stop("Please only pass in one molecular data type.") }},
        "mDataTypeNotChar" = {if (!is.character(mDataType)) { stop("mDataType must be a string.") }},
        "mDataTypeNotIn" = {if (!(mDataType %in% mDataNames(tSet))) { stop(paste0("The molecular data type(s) ", paste(mDataType[which(!(mDataType %in% mDataNames(tSet)))], collapse = ", " ), " is/are not present in ", tSet@annotation$name, ".")) }},
        # cell.lines checks
        "cell.linesNotChar" = {if (!is.character(unlist(cell.lines))) { stop("cell.lines parameter must contain strings.") }},
        "cell.linesNotIn" = {if (all(!(cell.lines %in% cellNames(tSet)))) { stop(paste0("The cell line(s) ", paste(cell.lines[which(!(cell.lines %in% cellNames(tSet)))], collapse = ", "), " is/are not present in ", tSet@annotation$name, "with the specified parameters.")) }},
        # drugs checks
        "drugsNotChar" = {if (!is.character(unlist(drugs))) { stop("drugs parameter must contain strings.") }},
        "drugsNotIn" = {if (all(!(drugs %in% drugNames(tSet)))) { stop(paste0("The drug(s) ", paste(drugs[which(!(drugs %in% drugNames(tSet)))], collapse = ", "), " is/are not present in ", tSet@annotation$name, ".")) }},
        # features checks
        "featuresLt2" = {if (length(fNames(tSet, mDataType)) < 2) { stop("Must include at least 2 features to calculate summary statistics") }},
        "featuresNotChar" = {if (!is.character(unlist(features))) { stop("features parameter contain strings.") }},
        "featuresNotIn" = {if (all(!(fNames(tSet, mDataType[1]) %in% features))) { stop(paste0("The feature(s) ", paste(features[which(!(features %in% fNames(tSet, mDataType[1])))], collapse = ", "), " is/are not present in ", tSet@annotation$name, ".")) }},
        # duration checks
        "durationNotChar" = {if (!is.character(unlist(duration))) { stop("duration parameter must contain strings.") }},
        "durationNotIn" = {if (all(!(duration %in% sensitivityInfo(tSet)$duration_h))) { stop(paste0("The duration(s) ", paste(duration[which(!(duration %in% sensitivityInfo(tSet)$duration_h))]), collapse = ", ", "is/are not present in ", tSet@annotation$name, ".")) }},
        # dose checks
        "doseLt2" = {if (length(dose) < 2) { stop("To fit a linear model we need at least two dose levels, please add anothor to the dose argument in the function call.") }},
        "doseNoCtl" = {if (!("Control" %in% dose)) { stop("You should not calculate summary statistics without including a control! Please add 'Control' to the dose argument vector.") }},
        "doseNotChar" = {if (!is.character(dose)) { stop("Dose must be a string or character vector.") }},
        "doseNotIn" = {if (all(!(dose %in% phenoInfo(tSet, mDataType)$dose_level))) { stop(paste0("The dose level(s) ", dose, " is/are not present in ", tSet@annotation$name, " with the specified parameters.")) }},
        # summary.stat
        "summary.statNotChar" = {if (!is.character(summary.stat)) { stop("The parameter summary.stat must be a string or character vector.") }},
        "summary.statGt1" = {if (length(summary.stat) > 1)  {stop("Please pick only one summary statistic") }},
        "summary.statNotIn" = {if (!(summary.stat %in% c("mean", "median", "first", "last"))) { stop(paste0("The the statistic ", summary.stat, " is not implemented in this package")) }},
        # viabilties checks
        "viabilitiesNotNum" = {if (!is.null(viabilities)) { if (!all(vapply(viabilities, function (viability) { is(viability, "numeric") }, FUN.VALUE = logical(1) ))) { stop("Viability values must be numeric.") }}},
        "viabilitiesNotMissing" = {if (!is.null(concentrations)) { if (is.null(viabilities)) { stop("If you pass in an argument for concentrations, you must also pass in an argument for viabilities.")}}},
        "viabilitiesDiffLenConc" = {
          if (!is.null(concentrations) && !is.null(viabilities)) {
            if (length(viabilities) != length(concentrations)) {
              stop(paste0(ifelse(is(viabilities, "list"), "List", "Vector"), " of viabilities is ", length(viabilities), "long, but concentrations is ", length(concentrations), "long.")) }}},
        "viabilitiesDiffLenDur" = {
          if (!is.null(duration) && !is.null(viabilities)) {
            if (length(viabilities) != length(duration)) {
              stop(paste0(ifelse(is(viabilities, "list"), "List", "Vector"), " of viabilities is ", length(viabilities), "long, but concentrations is ", length(duration), "long.")) }}},
        # concentrations checks
        "concentrationsNotNum" = {
          if (!is.null(concentrations)) {
            if (!all(vapply(concentrations, function (concentration) {
              is(concentration, "numeric") }, FUN.VALUE = logical(1) ))) {
              stop("Concentration values must be numeric.") }
          }}
      )
    )
  }
}
