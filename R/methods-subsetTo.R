#'`[`
#'
#' @examples
#' tSet <- TGGATESsmall[cellNames(TGGATESsmall), drugNames(TGGATESsmall)[seq_len(3)]]
#'
#'@param x tSet
#'@param i Cell lines to keep in tSet
#'@param j Drugs to keep in tSet
#'@param ... further arguments
#'@param drop A boolean flag of whether to drop single dimensions or not
#'@return Returns the subsetted tSet
#'@export
setMethod(`[`, "ToxicoSet", function(x, i, j, ..., drop = FALSE){
    if(is.character(i) && is.character(j)) {
        return(subsetTo(x, cells=i, drugs=j,  molecular.data.cells=i))
    }
    else if(is.numeric(i) && is.numeric(j) &&
            (as.integer(i)==i) && (as.integer(j)==j)) {
        return(subsetTo(x, cells=cellNames(x)[i], drugs=drugNames(x)[j],
                        molecular.data.cells=cellNames(x)[i]))
    }
})

#### subsetTo ####

## FIXED? TODO:: Subset function breaks if it doesnt find cell line in sensitivity info
#' A function to subset a ToxicoSet to data containing only specified drugs, cells and genes
#'
#' This is the prefered method of subsetting a ToxicoSet. This function allows
#' abstraction of the data to the level of biologically relevant objects: drugs
#' and cells. The function will automatically go through all of the
#' combined data in the ToxicoSet and ensure only the requested radiations
#' and cell lines are found in any of the slots. This allows quickly picking out
#' all the experiments for a radiation or cell of interest, as well removes the need
#' to keep track of all the metadata conventions between different datasets.
#'
#' @examples
#' TGGATESDrugNames  <- drugNames(TGGATESsmall)
#' TGGATESCells <- cellNames(TGGATESsmall)
#' tSet <- subsetTo(TGGATESsmall,drugs = TGGATESDrugNames[1],
#'   cells = TGGATESCells[1], duration = "2")
#'
#' @param object A \code{ToxicoSet} to be subsetted
#' @param cell_lines A list or vector of cell names as used in the dataset to which
#'   the object will be subsetted. If left blank, then all cells will be left in
#'   the dataset.
#' @param drugs A list or vector of drug names as used in the dataset to which
#'   the object will be subsetted. If left blank, then all drugs will be left in
#'   the dataset.
#' @param features A list or vector of feature names as used in the dataset from
#'   which the object will be subsetted. If left blank that all features will
#'   be left in.
#' @param molecular.data.cells A list or vector of cell names to keep in the
#'   molecular data
#' @param duration A \code{list} or \code{vector} of the experimental durations
#'   to include in the subset as strings. Defaults to all durations if parameter
#'   is not specified.
#' @param ... Other arguments passed to other functions within the package
#'
#' @return A ToxicoSet with only the selected drugs and cells
#'
#' @importFrom CoreGx .unionList
#' @export
## TODO:: Include dose parmater to subset on
subsetTo <- function(object, cell_lines = NULL,
                     drugs=NULL,
                     molecular.data.cells=NULL,
                     duration=NULL,
                     features=NULL,
                     ...
) {
    ## TODO:: Remove this or add it to the function parameters?
    drop = FALSE

    ####
    # PARSING ARGUMENTS
    ####
    adArgs = list(...)
    if ("exps" %in% names(adArgs)) {
        exps <- adArgs[["exps"]]
        if(is(exps, "data.frame")) {
            exps2 <- exps[[name(object)]]
            names(exps2) <- rownames(exps)
            exps <- exps2
        } else{
            exps <- exps[[name(object)]]
        }
    }else {
        exps <- NULL
    }

    if ("dose" %in% names(adArgs)) {
        ## TODO:: Add subsetting on dose
        stop("Due to the structure of tSets, subsetting on dose can only be done on
      specific slots - not on the entire tSet")
    }

    ## MISSING VALUE HANDLING FOR PARAMETERS
    # Get named list of default values for missing parameters
    argDefaultList <-
        paramMissingHandler(funName = "subsetTo", tSet = object,
                            drugs = drugs, cell_lines = cell_lines, features = features,
                            duration = duration)
    # Assign any missing parameter default values to function environment
    if (length(argDefaultList) > 0) {
        for (idx in seq_along(argDefaultList)) {
            assign(names(argDefaultList)[idx], argDefaultList[[idx]])
        }
    }

    # ERROR HANDLING FOR PARAMETERS
    paramErrorChecker(funName = "subsetTo", tSet = object,
                      cell_lines = cell_lines,
                      drugs = drugs, features = features,
                      duration = duration)

    ##TODO:: Add a value to tSet which indicates the experimental design!
    ##FIXME:: Don't hard code object names!
    if (name(object) %in% c("drugMatrix_rat", "EMEXP2458")) {
        if (!('DMSO' %in% drugs)) {
            drugs <- c(drugs, 'DMSO')
        }
    }

    ######
    # SUBSETTING MOLECULAR PROFILES SLOT
    ######
    ### TODO:: implement strict subsetting at this level!!!!

    ### the function missing does not work as expected in the context below, because the arguments are passed to the anonymous
    ### function in lapply, so it does not recognize them as missing
    molecularProfilesSlot(object) <-
        lapply(molecularProfilesSlot(object),
               function(SE, cell_lines, drugs, molecular.data.cells, duration, features){

                   if (!is.null(features)) {
                       SE <- SE[which(rownames(SummarizedExperiment::rowData(SE)) %in% features), ]
                   }

                   ##FIXME:: Why is are all these if conditions being checked against length? Just use grepl?
                   molecular.data.type <-
                       ifelse(
                           length(grep("rna", S4Vectors::metadata(SE)$annotation) > 0),
                           "rna",
                           S4Vectors::metadata(SE)$annotation
                       )

                   if (length(grep(molecular.data.type, names(molecular.data.cells))) > 0) {
                       cell_lines <- molecular.data.cells[[molecular.data.type]]
                   }
                   column_indices <- NULL

                   if (length(cell_lines) == 0 && length(drugs) == 0) {
                       column_indices <- seq_len(ncol(SE))
                   }
                   if (length(cell_lines) == 0 && datasetType(object) == "sensitivity") {
                       column_indices <- seq_len(ncol(SE))
                   }

                   # Selecting indices which match the cells argument
                   cell_line_index <- NULL
                   if (length(cell_lines) != 0) {
                       if (!all(cell_lines %in% cellNames(object))) {
                           stop("Some of the cell names passed to function did not match to names
          in the ToxicoSet. Please ensure you are using cell names as
          returned by the cellNames function")
                       }
                       cell_line_index <- which(SummarizedExperiment::colData(SE)[["cellid"]] %in% cell_lines)
                   }

                   # Selecting indexes which match drugs arguement
                   drugs_index <- NULL
                   if (datasetType(object) == "perturbation" || datasetType(object) == "both") {
                       if (length(drugs) != 0) {
                           if (!all(drugs %in% drugNames(object))){
                               stop("Some of the drug names passed to function did not match to names in the ToxicoSet Please ensure you are using drug names as returned by the drugNames function")
                           }
                           drugs_index <- which(SummarizedExperiment::colData(SE)[["drugid"]] %in% drugs)
                       }
                   }

                   if (length(drugs_index) != 0 && length(cell_line_index) != 0) {
                       if (length(intersect(drugs_index, cell_line_index)) == 0) {
                           stop("This Drug - Cell Line combination was not tested together.")
                       }
                       column_indices <- intersect(drugs_index, cell_line_index)
                   } else {
                       if (length(drugs_index) != 0) {
                           column_indices <- drugs_index
                       }
                       if (length(cell_line_index) != 0) {
                           column_indices <- cell_line_index
                       }
                   }

                   # LOGIC TO SUBSET BASED ON DURATION
                   ## TODO:: Determine if this works for other SE data types
                   if (!is.null(duration)){
                       if (all(!(duration %in% unique(SummarizedExperiment::colData(SE[, column_indices])$duration)))) {
                           # Error when other parameters are passed in
                           if ( !is.null(cell_lines) | !is.null(drugs) | !is.null(molecular.data.cells)) {
                               stop(paste0(
                                   "There are no molecular profiles with duration of ",
                                   duration, " in the tSet with the selected parameters."
                               ))
                           } else { # Error when no other parameters are passed in
                               stop(paste0(
                                   "There are no molecular profiles with duration of ",
                                   duration, " in the tSet."
                               ))
                           }
                       }
                       duration_indices <- which(SummarizedExperiment::colData(SE)$duration %in% duration)
                       column_indices <- intersect(column_indices, duration_indices)
                   }

                   row_indices <- seq_len(nrow(SummarizedExperiment::assay(SE, 1)))

                   # Final SE
                   SE <- SE[row_indices, column_indices]
                   return(SE)

               }, cell_lines = cell_lines,
               drugs = drugs,
               molecular.data.cells = molecular.data.cells,
               duration = duration,
               features = features)


    ######
    # SUBSET SENSITIVITY SLOT
    ######
    # Logic if any "..." arguments are passed to subsetTo
    if ((datasetType(object) == "sensitivity" | datasetType(object) == "both") & length(exps) != 0) {
        sensitivityInfo(object) <- sensitivityInfo(object)[exps, , drop=drop]
        rownames(sensitivityInfo(object)) <- names(exps)
        if (length(sensitivityRaw(object)) > 0) {
            sensitivityRaw(object) <- sensitivityRaw(object)[exps, , , drop=drop]
            dimnames(sensitivityRaw(object))[[1]] <- names(exps)
        }
        sensitivityProfiles(object) <- sensitivityProfiles(object)[exps, , drop=drop]
        rownames(sensitivityProfiles(object)) <- names(exps)

        sensNumber(object) <- .summarizeSensitivityNumbers(object)
    }
    # Logic if drug or cell parameters are passed to subsetTo
    else if (
        (datasetType(object) == "sensitivity" | datasetType(object) == "both") &
        (length(drugs) != 0 | length(cell_lines) != 0 | !is.null(duration) )
    ) {

        drugs_index <- which(sensitivityInfo(object)[, "drugid"] %in% drugs)
        cell_line_index <- which(sensitivityInfo(object)[,"cellid"] %in% cell_lines)
        if (length(drugs_index) !=0 & length(cell_line_index) !=0 ) {
            if (length(intersect(drugs_index, cell_line_index)) == 0) {
                stop("This Drug - Cell Line combination was not tested together.")
            }
            row_indices <- intersect(drugs_index, cell_line_index)
        } else {
            if(length(drugs_index)!=0 & length(cell_lines)==0) {
                row_indices <- drugs_index
            } else {
                if(length(cell_line_index)!=0 & length(drugs)==0){
                    row_indices <- cell_line_index
                } else {
                    # Includes all rows if cell or drug arguments are absent
                    row_indices <- seq_len(nrow(sensitivityInfo(object)))
                }
            }
        }
        # LOGIC TO SUBSET BASED ON DURATION
        if(!is.null(duration)){
            if(all(!(duration %in% unique(sensitivityInfo(object)[row_indices,]$duration_h)))) {
                # Error when other parameters are passed in
                if(!is.null(cell_lines) | !is.null(drugs) | !is.null(molecular.data.cells)) {
                    stop(paste0(
                        ## TODO:: Is sample the correct way to refer to one treatment/duration combination in TGx experiments?
                        "There are no samples with duration of ",
                        duration, " in the tSet with the selected parameters."
                    ))
                } else { # Error when no other parameters are passed in
                    stop(paste0(
                        "There are no samples with duration of ",
                        duration, " in the tSet"
                    ))
                }
            }
            duration_indices <- which(sensitivityInfo(object)$duration_h %in% duration)
            row_indices <- intersect(row_indices, duration_indices)
        }
        sensItemNames <- names(sensitivitySlot(object))
        sensitivityVals <-
            lapply(sensItemNames, function(sensItemName, drop){
                if (sensItemName == "n") {
                    sensItem <- sensitivitySlot(object)[[sensItemName]]
                    if (!is.null(cell_lines)) {
                        sensItem[which(rownames(sensItem) %in% cell_lines),
                                 which(colnames(sensItem) %in% drugs), drop = drop]
                    } else {
                        sensItem[ , which(colnames(sensItem) %in% drugs), drop = drop]
                    }
                } else {
                    sensItem <- sensitivitySlot(object)[[sensItemName]]
                    if (length(dim(sensItem)) == 3) {
                        sensItem[row_indices, , , drop = drop]
                    } else {
                        sensItem[row_indices, , drop = drop]
                    }
                }
            }, drop = drop)
        names(sensitivityVals) <- sensItemNames
        sensitivitySlot(object) <- sensitivityVals
    }

    #####
    # SUBSET DRUG SLOT
    #####
    if (length(drugs) == 0) {
        if (datasetType(object) == "sensitivity" | datasetType(object) == "both"){
            drugs <- unique(sensitivityInfo(object)[["drugid"]])
        }
        if(datasetType(object) == "perturbation" | datasetType(object) == "both"){
            drugs <- union(drugs, na.omit(.unionList(lapply(molecularProfilesSlot(object), function(SE){unique(SummarizedExperiment::colData(SE)[["drugid"]])}))))
        }
    }

    #####
    # SUBSET CELLS SLOT
    #####
    if (length(cell_lines) == 0) {
        cell_lines <- union(cell_lines, na.omit(.unionList(lapply(molecularProfilesSlot(object), function(SE){unique(SummarizedExperiment::colData(SE)[["cellid"]])}))))
        if (datasetType(object) == "sensitivity" | datasetType(object) == "both"){
            cell_lines <- union(cell_lines, sensitivityInfo(object)[["cellid"]])
        }
    }
    #####
    # ASSIGN SUBSETS BACK TO TOXICOSET OBJECT
    #####
    drugInfo(object) <- drugInfo(object)[drugs , , drop=drop]
    cellInfo(object) <- cellInfo(object)[cell_lines , , drop=drop]
    curation(object)$drug <- curation(object)$drug[drugs , , drop=drop]
    curation(object)$cell <- curation(object)$cell[cell_lines , , drop=drop]
    curation(object)$tissue <- curation(object)$tissue[cell_lines , , drop=drop]
    return(object)
}

#
# END SUBSET TO FUNCTION
#
