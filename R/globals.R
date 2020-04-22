# Declaring global variables for dplyr and data.table column names
utils::globalVariables(c('cellid','drugid','read.csv','samplename','.',
                         'Symbol', 'feature', 'Control', 'Low', 'Middle',
                         'High', 'verbose', 'dose_level', 'individual_id',
                         'duration_h', 'viability', '.SD', 'durations',
                         'tSetName'))

#' Define an S4 Generic for the methods::as function
#'
#' This will allow creation of new definitons for object conversions
#'
#' @param object An object to coerce to a different class
#' @param ... To allow new parameters for this generc
#' @param value The type to coerce the object to
#'
#' @export
setGeneric('as', function(object, ..., value) methods::as(object, ..., value))

#' Coerce a SummarizedExperiment object to an ExpressionSet
#'
#' @section Warning:
#' This method assumes that all slots not present in a SummarizedExperiment
#'   were moved into the metadata list from the original ExpressionSet
#'
#' @param object An ExpressionSet to coerce to a SummarizedExperiment
#' @param value A \code{character} vector specifying the type to coerce to
#'
#' @importFrom SummarizedExperiment colData rowData assays assay
#' @importFrom S4Vectors metadata
#' @export
setMethod('as',
          signature('SummarizedExperiment'),
          function(object, value) {
            if (value != 'ExpressionSet')
              as(object, value)
            else
              Biobase::ExpressionSet(
                as.environment(as.list(assays(object))),
                phenoData=as(colData(object), 'AnnotatedDataFrame'),
                featureData=as(rowData(object), 'AnnotatedDataFrame'),
                experimentData=metadata(object)$experimentData,
                annotation=metadata(object)$annotation,
                protocolData=metadata(object)$protocolData
              )
          })

#' Coerce a DFrame object to an AnnotatedDataFrame
#'
#' @param object A DFrame object to coerce to AnnotatedDataFrame
#' @param value The type to convert to; if not "AnnotatedDataFrame" this
#'   method falls back to as
#'
#' @importFrom Biobase AnnotatedDataFrame
#' @export
setMethod('as',
          signature('DFrame'),
          function(object, value) {
            if (value != "AnnotatedDataFrame")
              as(object, value)
            else
              AnnotatedDataFrame(as.data.frame(object))
          })