#' An Error Checker for Function Parameters
#'
#' This function will take in the params of a function as well as its name.
#'   Error checking will then be conducted on each of it's parameter arguments
#'   to ensure they meet the input requirements for that function. Descriptive
#'   errosr are returned if the the parameters do not meet the criteria for that
#'   function.
#'
#' @param funName [character] A string of the function name. This argument is
#'   used to match the correct conditions with each function.
#' @param ... [list] A list of all parameters passed into the function whos
#'   name is supplied in funName parameter.
#'
#' @return Returns nothing, this function works by side effects only
#'
#' @export
#' @keywords internal
paramErrorChecker <- function(funName, tSet, ...) {

  # Extract named arguments into local environment
  argList <- list(...)
  print(argList)
  for (idx in seq_len(length(argList))) {
    assign(
      gsub("^.", "", names(argList)[idx]), # Removes $ at start of var names
      argList[[idx]])
  }

  # Matches the correct parameter constraints to each function name
  paramChecks <-
    switch(funName,
           "drugPerturbationSig" =
               c("tSetGt1",
                 "mDataTypeGt1", "mDataTypeNotStr", "mDataTypeNotIn",
                 "cell.linesNotChar",  "cell.linesNotIn",
                 "drugsNotChar", "drugsNotIn",
                 "featuresLt2", "featuresNotChar", "featuresNotIn",
                 "drugsNotChar", "durationNotIn",
                 "doseNotChar", "doseNotIn", , "doseLt2", "doseNoCtl", 
                  ),
           "summarizeMolecularProfiles" = 
              c("tSetGt1",
                "mDataTypeNotStr","mDataTypeNotIn",
                "cell.linesNotChar", "cell.linesNotIn",
                "drugsNotChar", "drugsNotIn",
                "featuresNotChar", "featuresNotIn",
                "durationNotChar", "durationNotIn",
                "doseNotChar", "doseNotIn",
                "summary.statNotChar", "summary.statNotIn", "summary.statGt1"
               )
           )

  # Runs the parameter checks specific to the given funName
  invisible(
    dplyr::case_when(
      # Checks for tSet
      "tSetGt1" %in% paramChecks ~ ifelse(length(unlist(tSet)) > 1, stop("You may only pass in one tSet."), "" ),
      # mDataType checks
      "mDataTypeGt1" %in% paramChecks ~ ifelse(length(unlist(mDataType)) > 1, stop("Please only pass in one molecular data type."), "" ),
      "mDataTypeNotChar" %in% paramChecks ~ ifelse(!is.character(mDataType), stop("mDataType must be a string."), "" ),
      "mDataTypeNotIn" %in% paramChecks ~ ifelse(all(!(mDataType %in% mDataNames(unlist(tSet))), stop(paste0("The molecular data type(s) ", paste(mDataType[which(!(mDataType %in% mDataNames(tSet)))], collapse = ", " ), " is/are not present in ", tSet@annotation$name, ".")), "")),
      # cell.lines checks
      "cell.linesNotChar" %in% paramChecks ~ ifelse(!is.character(unlist(cell.lines)), stop("cell.lines parameter must contain strings."), ""),
      "cell.linesNotIn" %in% paramChecks ~ ifelse(all(!(cell.lines %in% cellNames(tSet))), stop(paste0("The cell line(s) ", paste(cell.lines[which(!(cell.lines %in% cellNames(tSet)))], collapse = ", "), " is/are not present in ", tSet@annotation$name, "with the specified parameters.")), ""),
      # drugs checks
      "drugsNotChar" %in% paramChecks ~ ifelse(!is.character(unlist(drugs)), stop("drugs parameter must contain strings."), ""),
      "drugsNotIn" %in% paramChecks ~ ifelse(all(!(drugs %in% drugNames(tSet))), stop(paste0("The drug(s) ", paste(drugs[which(!(drugs %in% drugNames(tSet)))], collapse = ", "), " is/are not present in ", tSet@annotation$name, ".")), ""),
      # features checks
      ## TODO:: Do we want to implement this function with 1 feature?
      "featuresLt2" %in% paramChecks ~ ifelse(length(fNames(tSet, mDataType)) < 2, stop("Must include at least 2 features to calculate summary statistics"), ""),
      "featuresNotChar" %in% paramChecks ~ ifelse(!is.character(unlist(features)), stop("features parameter contain strings."), ""),
      "featuresNotIn" %in% paramChecks ~ ifelse(all(!(fNames(tSet, mDataType[1]) %in% features)), stop(paste0("The feature(s) ", paste(features[which(!(features %in% fNames(tSet, mDataType[1])))], collapse = ", "), " is/are not present in ", tSet@annotation$name, ".")), ""),
      # duration checks
      "durationNotChar" %in% paramChecks ~ ifelse(!is.character(unlist(duration)), stop("duration parameter must contain strings."), ""),
      "durationNotIn" %in% paramChecks ~ ifelse(all(!(duration %in% sensitivityInfo(tSet)$duration_h)), stop(paste0("The duration(s) ", paste(duration[which(!(duration %in% sensitivityInfo(tSet)$duration_h))]), collapse = ", ", "is/are not present in ", tSet@annotation$name, ".")), ""),
      # dose checks
      "doseLt2" %in% paramChecks ~ ifelse(length(dose) < 2, stop("To fit a linear model we need at least two dose levels, please add anothor to the dose argument in the function call."), ""),
      "doseNoCtl" %in% paramChecks ~ ifelse(!("Control" %in% dose), stop("You should not calculate summary statistics without including a control! Please add 'Control' to the dose argument vector."),""),
      "doseNotChar" %in% paramChecks ~ ifelse(!is.character(dose), stop("Dose must be a string or character vector."), ""),
      "doseNotIn" %in% paramChecks ~ ifelse(all(!(dose %in% phenoInfo(tSet, mDataType)$dose_level)), stop(paste0("The dose level(s) ", dose, " is/are not present in ", tSet@annotation$name, " with the specified parameters.")), ""),
      # summary.stat
      "summary.statNotChar" %in% paramChecks ~ ifelse(!is.character(dose), stop("Dose must be a string or character vector."), ""),
      "summary.statGt1" %in% paramChecks ~ ifelse(length(unlist(summary.stat)) > 1, stop("Please pick only one summary statistic"), "" )),
      "summary.statNotIn" %in% paramChecks ~ ifelse(!(summary.stat %in% c("mean", "median", "first", "last")), stop(paste0("The the statistic ", summary.stat, " is not implemented in this package")), "")
  )
}





