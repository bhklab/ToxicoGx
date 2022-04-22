#' @include ToxicoSet-class.R
NULL

# Navigating this file:
# - Slot section names start with ----
# - Method section names start with ==
#
# As a result, you can use Ctrl + f to find the slot or method you are looking
# for quickly, assuming you know its name.
#
# For example Ctrl + f '== molecularProfiles' would take you the molecularProfiles
# method, while Ctrl + f '---- molecularProfiles' would take you to the slot
# section.


#### CoreGx dynamic documentation
####
#### Warning: for dynamic docs to work, you must set
#### Roxygen: list(markdown = TRUE, r6=FALSE)
#### in the DESCRPTION file!


# =======================================
# Accessor Method Documentation Object
# ---------------------------------------

#' @name ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_accessors(class_='CoreSet')
#' @eval CoreGx:::.parseToRoxygen("@examples data({data_})", data_=.local_data)
NULL



# ======================================
# Accessor Methods
# --------------------------------------


## ==============
## ---- drug slot
## --------------


##
## == drugInfo

#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_treatmentInfo(class_=.local_class,
#' data_=.local_data)
#' @importMethodsFrom CoreGx treatmentInfo
#' @aliases drugInfo
#' @export
drugInfo <- function(...) treatmentInfo(...)

#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_treatmentInfo(class_=.local_class,
#' data_=.local_data)
#' @importMethodsFrom CoreGx treatmentInfo<-
#' @aliases drugInfo<-
#' @export
`drugInfo<-` <- function(..., value) `treatmentInfo<-`(..., value=value)



##
## == drugNames


#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_treatmentNames(class_=.local_class,
#' data_=.local_data)
#' @importMethodsFrom CoreGx treatmentNames
#' @aliases drugNames
#' @export
drugNames <- function(...) treatmentNames(...)


#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_treatmentNames(class_=.local_class,
#' data_=.local_data)
#' @importMethodsFrom CoreGx treatmentNames<-
#' @aliases drugNames<-
#' @export
`drugNames<-` <- function(..., value) `treatmentNames<-`(..., value=value)


## ====================
## ---- annotation slot
## --------------------


##
## == annotation


#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_annotation(class_=.local_class, data_=.local_data)
#' @importMethodsFrom CoreGx annotation
setMethod('annotation', signature("ToxicoSet"), function(object) {
    callNextMethod(object=object)
})

#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_annotation(class_=.local_class, data_=.local_data)
#' @importMethodsFrom CoreGx annotation<-
setReplaceMethod("annotation", signature("ToxicoSet", "list"),
        function(object, value) {
    callNextMethod(object=object, value=value)
})


##
## == dateCreated


#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_dateCreated(class_=.local_class, data_=.local_data)
#' @importMethodsFrom CoreGx dateCreated
setMethod('dateCreated', signature("ToxicoSet"), function(object) {
    callNextMethod(object=object)
})

#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_dateCreated(class_=.local_class, data_=.local_data)
#' @importMethodsFrom CoreGx dateCreated<-
setReplaceMethod('dateCreated', signature(object="ToxicoSet", value="character"),
    function(object, value)
{
    callNextMethod(object=object, value=value)
})


##
## === name


#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_name(class_=.local_class, data_=.local_data)
#' @importMethodsFrom CoreGx name
setMethod('name', signature("ToxicoSet"), function(object) {
    callNextMethod(object)
})

#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_name(class_=.local_class, data_=.local_data)
#' @importMethodsFrom CoreGx name<-
setReplaceMethod('name', signature("ToxicoSet"), function(object, value) {
    object <- callNextMethod(object, value=value)
    return(invisible(object))
})

## ==============
## ---- sample slot
## --------------


##
## == sampleInfo

#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_sampleInfo(class_=.local_class,
#' sample_=.local_sample)
#' @importMethodsFrom CoreGx sampleInfo
#' @importFrom CoreGx cellInfo
#' @export
setMethod("sampleInfo", "ToxicoSet", function(object) {
    callNextMethod(object)
})


#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_sampleInfo(class_=.local_class,
#' data_=.local_data, sample_="cell")
#' @importMethodsFrom CoreGx sampleInfo<-
#' @importFrom CoreGx cellInfo<-
#' @export
setReplaceMethod("sampleInfo", signature(object="ToxicoSet",
        value="data.frame"), function(object, value) {
    callNextMethod(object, value=value)
})


##
## == sampleNames


#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_sampleNames(class_=.local_class,
#' data_=.local_data, sample_=.local_sample)
#' @importMethodsFrom CoreGx sampleNames
setMethod("sampleNames", signature("ToxicoSet"), function(object) {
    callNextMethod(object)
})


#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_sampleNames(class_=.local_class,
#' data_=.local_data, sample_=.local_sample)
#' @importMethodsFrom CoreGx sampleNames<-
setReplaceMethod("sampleNames", signature(object="ToxicoSet", value="character"),
        function(object, value) {
    callNextMethod(object=object, value=value)
})



## ------------------
## ---- curation slot


##
## == curation


#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_curation(class_=.local_class,
#' data_=.local_data, details_="Contains three `data.frame`s, 'cell' with
#' cell-line ids and 'tissue' with tissue ids and 'drug' with drug ids.")
#' @importMethodsFrom CoreGx curation
setMethod('curation', signature(object="ToxicoSet"), function(object) {
    callNextMethod(object=object)
})

#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_curation(class_=.local_class,
#' data_=.local_data, details_="For a `ToxicoSet` object the slot should
#' contain tissue, cell-line and drug id `data.frame`s.")
#' @importMethodsFrom CoreGx curation<-
setReplaceMethod("curation", signature(object="ToxicoSet", value="list"),
    function(object, value)
{
    callNextMethod(object=object, value=value)
})


## ----------------------
## ---- datasetType slot


#
# == datasetType


#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_datasetType(class_=.local_class,
#' data_=.local_data)
#' @importMethodsFrom CoreGx datasetType
setMethod("datasetType", signature("ToxicoSet"), function(object) {
    callNextMethod(object)
})

#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_datasetType(class_=.local_class,
#' data_=.local_data)
#' @importMethodsFrom CoreGx datasetType<-
setReplaceMethod("datasetType", signature(object="ToxicoSet",
    value='character'), function(object, value)
{
    callNextMethod(object=object, value=value)
})


## ---------------------------
## ---- molecularProfiles slot


##
## == molecularProfiles


#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_molecularProfiles(class_=.local_class,
#' data_=.local_data)
#' @importMethodsFrom CoreGx molecularProfiles
setMethod(molecularProfiles, "ToxicoSet", function(object, mDataType, assay)
{
    callNextMethod(object=object, mDataType=mDataType, assay=assay)
})

#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_molecularProfiles(class_=.local_class,
#' data_=.local_data)
#' @importMethodsFrom CoreGx molecularProfiles<-
setReplaceMethod("molecularProfiles", signature(object="ToxicoSet",
    mDataType ="character", assay="character", value="matrix"),
    function(object, mDataType, assay, value)
{
    callNextMethod(object=object, mDataType=mDataType, assay=assay, value=value)
})
setReplaceMethod("molecularProfiles",
    signature(object="ToxicoSet", mDataType ="character", assay="missing",
        value="matrix"), function(object, mDataType, assay, value)
{
    callNextMethod(object=object, mDataType=mDataType, assay=assay, value=value)
})


##
## == featureInfo


#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_featureInfo(class_=.local_class,
#' data_=.local_data)
#' @importMethodsFrom CoreGx featureInfo
setMethod(featureInfo, "ToxicoSet", function(object, mDataType) {
    callNextMethod(object=object, mDataType=mDataType)
})

#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_featureInfo(class_=.local_class,
#' data_=.local_data, mDataType_='rna')
#' @importMethodsFrom CoreGx featureInfo<-
setReplaceMethod("featureInfo", signature(object="ToxicoSet",
    mDataType ="character",value="data.frame"),
    function(object, mDataType, value)
{
    callNextMethod(object=object, mDataType=mDataType, value=value)
})
setReplaceMethod("featureInfo", signature(object="ToxicoSet",
    mDataType ="character",value="DataFrame"),
    function(object, mDataType, value)
{
    callNextMethod(object=object, mDataType=mDataType, value=value)
})



##
## == phenoInfo


#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_phenoInfo(class_=.local_class,
#' data_=.local_data, mDataType_='rna')
#' @importMethodsFrom CoreGx phenoInfo
setMethod('phenoInfo', signature(object='ToxicoSet', mDataType='character'),
    function(object, mDataType)
{
    callNextMethod(object=object, mDataType=mDataType)
})

#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_phenoInfo(class_=.local_class,
#' data_=.local_data, mDataType_='rna')
#' @importMethodsFrom CoreGx phenoInfo<-
setReplaceMethod("phenoInfo", signature(object="ToxicoSet",
    mDataType ="character", value="data.frame"),
    function(object, mDataType, value)
{
    callNextMethod(object=object, mDataType=mDataType, value=value)
})
setReplaceMethod("phenoInfo", signature(object="ToxicoSet",
    mDataType ="character", value="DataFrame"),
    function(object, mDataType, value)
{
    callNextMethod(object=object, mDataType=mDataType, value=value)
})


##
## == fNames


#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_fNames(class_=.local_class,
#' data_=.local_data, mDataType_='rna')
#' @importMethodsFrom CoreGx fNames
setMethod('fNames', signature(object='ToxicoSet', mDataType='character'),
    function(object, mDataType)
{
    callNextMethod(object=object, mDataType=mDataType)
})

#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_fNames(class_=.local_class,
#' data_=.local_data, mDataType_='rna')
#' @importMethodsFrom CoreGx fNames<-
setReplaceMethod('fNames', signature(object='ToxicoSet', mDataType='character',
    value='character'), function(object, mDataType, value)
{
    callNextMethod(object=object, mDataType=mDataType, value=value)
})


##
## == mDataNames


#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_mDataNames(class_=.local_class,
#' data_=.local_data)
#' @importMethodsFrom CoreGx mDataNames
setMethod("mDataNames", "ToxicoSet", function(object){
    callNextMethod(object=object)
})

#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_mDataNames(class_=.local_class,
#' data_=.local_data)
#' @importMethodsFrom CoreGx mDataNames<-
setReplaceMethod("mDataNames", "ToxicoSet", function(object, value){
    callNextMethod(object=object, value=value)
})



##
## == molecularProfilesSlot


#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_molecularProfilesSlot(class_=.local_class,
#' data_=.local_data)
#' @importMethodsFrom CoreGx molecularProfilesSlot
setMethod("molecularProfilesSlot", signature("ToxicoSet"), function(object) {
    callNextMethod(object=object)
})

#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_molecularProfilesSlot(class_=.local_class,
#' data_=.local_data)
#' @importMethodsFrom CoreGx molecularProfilesSlot<-
setReplaceMethod("molecularProfilesSlot", signature("ToxicoSet", "list_OR_MAE"),
    function(object, value)
{
    callNextMethod(object=object, value=value)
})


# ---------------------
## ---- sensitivity slot


##
## == sensitivityInfo

#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_sensitivityInfo(class_=.local_class,
#' data_=.local_data)
#' @importMethodsFrom CoreGx sensitivityInfo
setMethod('sensitivityInfo', signature("ToxicoSet"),
    function(object, dimension, ...)
{
    callNextMethod(object=object, dimension=dimension, ...)
})

#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_sensitivityInfo(class_=.local_class,
#' data_=.local_data)
#' @importMethodsFrom CoreGx sensitivityInfo<-
setReplaceMethod("sensitivityInfo", signature(object="ToxicoSet",
    value="data.frame"), function(object, dimension, ..., value)
{
    callNextMethod(object=object, dimension=dimension, ..., value=value)
})


##
## == sensitvityMeasures


#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_sensitivityMeasures(class_=.local_class,
#' data_=.local_data)
#' @importMethodsFrom CoreGx sensitivityMeasures
setMethod('sensitivityMeasures', signature(object="ToxicoSet"),
    function(object)
{
    callNextMethod(object=object)
})

#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_sensitityMeasures(class_=.local_class,
#' data_=.local_data)
setReplaceMethod('sensitivityMeasures',
    signature(object='ToxicoSet', value='character'), function(object, value)
{
    callNextMethod(object=object, value=value)
})


##
## == sensitivityProfiles


#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_sensitivityProfiles(class_=.local_class,
#' data_=.local_data)
#' @importMethodsFrom CoreGx sensitivityProfiles
setMethod('sensitivityProfiles', signature(object="ToxicoSet"), function(object)
{
    callNextMethod(object=object)
})

#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_sensitivityProfiles(class_=.local_class,
#' data_=.local_data)
#' @importMethodsFrom CoreGx sensitivityProfiles<-
setReplaceMethod("sensitivityProfiles",
    signature(object="ToxicoSet", value="data.frame"),
    function(object, value)
{
    callNextMethod(object=object, value=value)
})


#
# == sensitivityRaw


#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_sensitivityRaw(class_=.local_class,
#' data_=.local_data)
#' @importMethodsFrom CoreGx sensitivityRaw
setMethod("sensitivityRaw", signature("ToxicoSet"), function(object) {
    callNextMethod(object=object)
})

#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_sensitivityRaw(class_=.local_class,
#' data_=.local_data)
#' @importMethodsFrom CoreGx sensitivityRaw<-
setReplaceMethod('sensitivityRaw', signature("ToxicoSet", "array"),
    function(object, value)
{
    callNextMethod(object=object, value=value)
})


#
# == sensitivitySlot


#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_sensitivitySlot(class_=.local_class,
#'   data_=.local_data)
#' @importMethodsFrom CoreGx sensitivitySlot
setMethod("sensitivitySlot", signature("ToxicoSet"), function(object) {
    callNextMethod(object=object)
})



#' @rdname ToxicoSet-accessors
#' @importMethodsFrom CoreGx sensitivitySlot<-
#' @eval CoreGx:::.docs_CoreSet_set_sensitivitySlot(class_=.local_class,
#' data_=.local_data)
setReplaceMethod('sensitivitySlot', signature(object='ToxicoSet',
    value='list_OR_LongTable'), function(object, value)
{
    callNextMethod(object=object, value=value)
})


##
## == sensNumber


#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_sensNumber(class_=.local_class,
#' data_=.local_data)
#' @importMethodsFrom CoreGx sensNumber
setMethod('sensNumber', "ToxicoSet", function(object){
    callNextMethod(object=object)
})

#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_sensNumber(class_=.local_class,
#' data_=.local_data)
#' @importMethodsFrom CoreGx sensNumber<-
setReplaceMethod('sensNumber', signature(object="ToxicoSet", value="matrix"),
        function(object, value) {
    callNextMethod(object=object, value=value)
})


## ======================
## ---- perturbation slot


##
## == pertNumber


#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_pertNumber(class_=.local_class,
#' data_=.local_data)
#' @importMethodsFrom CoreGx pertNumber
setMethod('pertNumber', signature(object='ToxicoSet'), function(object) {
    callNextMethod(object=object)
})

#' @rdname ToxicoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_pertNumber(class_=.local_class,
#' data_=.local_data)
#' @importMethodsFrom CoreGx pertNumber<-
setReplaceMethod('pertNumber', signature(object='ToxicoSet', value="array"),
        function(object, value) {
    callNextMethod(object=object, value=value)
})