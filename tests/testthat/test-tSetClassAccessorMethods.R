library(testthat)
library(ToxicoGx)

#### Checking example tSet structure ####
context("Testing TGGATESsmall object validity...")

test_that("TGGATESsmall tSet has correct structure", {
    data(TGGATESsmall)
    expect_error(checkTSetStructure(TGGATESsmall), NA)
})

#### tSet Acessor methods ####
context("Testing tSet Class Accessor Methods...")

# @treatment Slot
test_that("@treatment slot accessors produce expected results", {
    data("TGGATESsmall")
    context("External validation...")
    expect_equal_to_reference(treatmentInfo(TGGATESsmall),
        "drugInfo.TGGATESsmall.rds")
    expect_equal_to_reference(treatmentNames(TGGATESsmall),
        "drugNames.TGGATESsmall.rds")

    context("Internal validation...")
    expect_equal(treatmentInfo(TGGATESsmall), TGGATESsmall@treatment)
    expect_equal(treatmentNames(TGGATESsmall), TGGATESsmall@treatment$treatmentid)
})

# @annotation Slot
test_that("@annotation slot accessors produce expected results", {
    data("TGGATESsmall")

    context("External validation...")
    expect_equal_to_reference(TGGATESsmall@annotation,
        "annotation.TGGATESsmall.rds")
    expect_equal_to_reference(name(TGGATESsmall), "name.TGGATESsmall.rds")

    context("Internal validation...")
    expect_equal(name(TGGATESsmall), TGGATESsmall@annotation$name)
})

# @molecularProfile Slot
test_that("@molecularProfiles slot accessors produce expected results", {
    data("TGGATESsmall")

    context("External validation...")
    expect_equal_to_reference(mDataNames(TGGATESsmall),
        "mDataNames.TGGATESsmall.rds")
    context("Internal validation...")
    expect_equal(mDataNames(TGGATESsmall),
        names(TGGATESsmall@molecularProfiles))

    ## TODO:: Test this with incorrect tSet structure to determine if error
    ##>messages print in the correct order
    for (name in names(TGGATESsmall@molecularProfiles)) {
        context("External validation...")
        expect_equal_to_reference(
            molecularProfiles(TGGATESsmall, name)[, 1:100],
            paste0(name, ".molecularProfiles.TGGATESsmall.rds")
        )
        expect_equal_to_reference(featureInfo(TGGATESsmall, name),
            paste0(name, ".featureInfo.TGGATESsmall.rds"))
        expect_equal_to_reference(fNames(TGGATESsmall, name),
            paste0(name, ".fNames.TGGATESsmall.rds"))
        expect_equal_to_reference(phenoInfo(TGGATESsmall, name)[1:100, ],
            paste0(name, ".phenoData.TGGATESsmall.rds"))
        context("Internal validation...")
        # expect_equal(molecularProfiles(TGGATESsmall, name),
        #     assay(TGGATESsmall@molecularProfiles[[name]], 1))
        expect_equal(featureInfo(TGGATESsmall, name),
            rowData(TGGATESsmall@molecularProfiles[[name]]))
        expect_equal(fNames(TGGATESsmall, name),
            rownames(rowData(TGGATESsmall@molecularProfiles[[name]])))
        expect_equal(phenoInfo(TGGATESsmall, name),
            colData(TGGATESsmall@molecularProfiles[[name]]))
    }
})

# @sample Slot
test_that("@sample slot accessors produce expected results", {
    data("TGGATESsmall")

    context("External validation...")
    expect_equal_to_reference(sampleInfo(TGGATESsmall),
        "cellInfo.TGGATESsmall.rds")
    expect_equal_to_reference(sampleNames(TGGATESsmall),
        "cellNames.TGGATESsmall.rds")

    context("Internal validation...")
    expect_equal(sampleInfo(TGGATESsmall), TGGATESsmall@sample)
    expect_equal(sampleNames(TGGATESsmall), TGGATESsmall@sample$sampleid)
})

# @sensitivty Slot
test_that("@sensitivity slot accessors produce expected results", {
    data("TGGATESsmall")

    context("External validation...")
    expect_equal_to_reference(sensitivityInfo(TGGATESsmall),
        "sensitivityInfo.TGGATESsmall.rds")
    expect_equal_to_reference(sensitivityProfiles(TGGATESsmall),
        "sensitivitProfiles.TGGATESsmall.rds")
    expect_equal_to_reference(sensitivityMeasures(TGGATESsmall),
        "sensitivityMeasures.TGGATESsmall.rds")
    expect_equal_to_reference(sensNumber(TGGATESsmall),
        "sensNumber.TGGATESsmall.rds")
    context("Internal validation...")
    expect_equal(sensitivityInfo(TGGATESsmall), TGGATESsmall@sensitivity$info)
    expect_equal(sensitivityProfiles(TGGATESsmall),
        TGGATESsmall@sensitivity$profiles)
    expect_equal(sensitivityMeasures(TGGATESsmall),
        colnames(TGGATESsmall@sensitivity$profiles))
    expect_equal(sensNumber(TGGATESsmall), TGGATESsmall@sensitivity$n)
})

# @perturbation Slot
test_that("@perturbation slot accessors produce expected results", {
    data("TGGATESsmall")
    context("External validation...")
    expect_equal_to_reference(TGGATESsmall@perturbation,
        "perturbation.TGGATESsmall.rds")
    expect_equal_to_reference(pertNumber(TGGATESsmall),
        "pertNumber.TGGATESsmall.rds")
    context("Internal validation...")
    expect_equal(pertNumber(TGGATESsmall), TGGATESsmall@perturbation$n)
})

# @curation Slot
test_that("@curation slot accessors produce expected results", {
    data("TGGATESsmall")
    context("External validation...")
    expect_equal_to_reference(TGGATESsmall@curation,
        "curation.TGGATESsmall.rds")
})

# subsetTo Method
test_that("subsetTo() class method produces expected results", {
    data("TGGATESsmall")

    ## TODO:: Add unit tests for `[` subset operator
    ## TODO:: Change context() messages to be more informative when
    ##>running devtools::test()
    context("External validation...")
    expect_equal_to_reference(
        subsetTo(TGGATESsmall, drugs = treatmentNames(TGGATESsmall)[1],
            cell_lines=sampleNames(TGGATESsmall)[1]),
        "subsetTo.TGGATESsmall.rds")
    context("Internal validation...")
    ## Tests that subsetting molecularProfiles on duration works
    expect_equal(all(
        sensitivityInfo(subsetTo(TGGATESsmall, duration = "2"))$duration_h
            %in% "2"),
        TRUE)
    # Tests that relationship between sensitivity experiments and
    #>molecularProfiles is preserved
    #>(4 molecular Profiles / 1 sensitivity experiment)
    for (name in names(molecularProfilesSlot(TGGATESsmall))) {
        testthat::context(paste0("Testing subsetTo on molecularProfile for ",
            name))
        ## TODO:: Generalize duration arguement so that it uses the first
        ##>unique duration value in tSet (replace "8" with this)
        testthat::expect_equal(all(
            SummarizedExperiment::colData(
                ToxicoGx::subsetTo(TGGATESsmall, duration = "8"
                    )@molecularProfiles[[name]])$duration %in% "8"),
            TRUE)
    }
})
