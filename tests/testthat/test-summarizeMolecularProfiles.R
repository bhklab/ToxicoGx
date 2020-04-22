library(testthat)
library(ToxicoGx)

## TODO:: Can probably rewrite this using an apply function?

context("Testing if summarizeMolecularProfiles error handling works correctly...")

# tSet
context("...Checking for correct tSet param errors...")
test_that("Errors if given more than one tSet as parameter.", { expect_error(
  summarizeMolecularProfiles(
    c(TGGATESsmall, TGGATESsmall), mDataType = "rna", cell_lines=cellNames(TGGATESsmall),
    drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall,"rna"), duration = "8",
    dose = c("Control", "High"), summary.stat = "median",
    fill.missing = TRUE, verbose=TRUE
    )
  )
})

# mDataTypes
context("...Checking for correct mDataType param errors...")
test_that("Errors if given more than one mDataType.", { expect_error(
  summarizeMolecularProfiles(
    TGGATESsmall, mDataType = c("rna", "cnv"), cell_lines=cellNames(TGGATESsmall),
    drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall,"rna"), duration = "8",
    dose = c("Control", "High"), summary.stat = "median",
    fill.missing = FALSE, verbose=TRUE
    )
  )
})
test_that("Errors if given mDataType as type other than character.", { expect_error(
  summarizeMolecularProfiles(
    TGGATESsmall, mDataType = 1, cell_lines=cellNames(TGGATESsmall),
    drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall,"rna"), duration = "8",
    dose = c("Low", "Medium"), summary.stat = "mean",
    fill.missing = FALSE, verbose=FALSE
    )
  )
})
test_that("Errors if specified mDataType is not in the tSet.", { expect_error(
  summarizeMolecularProfiles(
    TGGATESsmall, mDataType = "cnv", cell_lines=cellNames(TGGATESsmall),
    drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall,"rna"), duration = "8",
    dose = c("Control", "High"), summary.stat = "first",
    fill.missing = TRUE, verbose=TRUE
    )
  )
})

# cell_lines
context("...Checking for correct cell_lines param errors...")
test_that("Errors if given cell_lines as type other than character.", { expect_error(
  summarizeMolecularProfiles(
    TGGATESsmall, mDataType = "rna", cell_lines=5,
    drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall,"rna"), duration = "8",
    dose = c("Control", "High"), summary.stat = "last",
    fill.missing = TRUE, verbose=TRUE
    )
  )
})
test_that("Errors if specified cell_lines are not in the tSet", { expect_error(
  summarizeMolecularProfiles(
    TGGATESsmall, mDataType = "rna", cell_lines='NOTINtSET',
    drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall,"rna"), duration = "8",
    dose = c("Control", "High"), summary.stat = "median",
    fill.missing = TRUE, verbose=TRUE
    )
  )
})

# drugs
context("...Checking for correct drugs param errors...")
test_that("Errors if given drugs are type other than character.", { expect_error(
  summarizeMolecularProfiles(
    TGGATESsmall, mDataType = "rna", cell_lines=cellNames(TGGATESsmall),
    drugs = 5,
    features = fNames(TGGATESsmall,"rna"), duration = "8",
    dose = c("Control", "High"), summary.stat = "mean",
    fill.missing = TRUE, verbose=TRUE
    )
  )
})
test_that("Errors if specified drugs are not in the tSet.", { expect_error(
  summarizeMolecularProfiles(
    TGGATESsmall, mDataType = "rna", cell_lines=cellNames(TGGATESsmall),
    drugs = "NOTINtSET",
    features = fNames(TGGATESsmall,"rna"), duration = "8",
    dose = c("Control", "High"), summary.stat = "first",
    fill.missing = TRUE, verbose=TRUE
    )
  )
})

# features
context("...Checking for correct features param errors....")
test_that("Errors if given features as type other than character", { expect_error(
  summarizeMolecularProfiles(
    TGGATESsmall, mDataType = "rna", cell_lines=cellNames(TGGATESsmall),
    drugs = head(drugNames(TGGATESsmall)),
    features = c(5), duration = "8",
    dose = c("Control", "High"), summary.stat = "last",
    fill.missing = TRUE, verbose=TRUE
    )
  )
})
test_that("Errors if given features as type other than character", { expect_error(
  summarizeMolecularProfiles(
    TGGATESsmall, mDataType = "rna", cell_lines=cellNames(TGGATESsmall),
    drugs = head(drugNames(TGGATESsmall)),
    features = c(5), duration = "8",
    dose = c("Control", "High"), summary.stat = "median",
    fill.missing = TRUE, verbose=TRUE
    )
  )
})

# duration
context("...Checking for correct duration param errors")
test_that("Errors if given duration as type other than character", { expect_error(
  summarizeMolecularProfiles(
    TGGATESsmall, mDataType = "rna", cell_lines=cellNames(TGGATESsmall),
    drugs = head(drugNames(TGGATESsmall)),
    features = c(5), duration = 8,
    dose = c("Control", "High"), summary.stat = "mean",
    fill.missing = TRUE, verbose=TRUE
    )
  )
})
test_that("Errors if given features as type other than character", { expect_error(
  summarizeMolecularProfiles(
    TGGATESsmall, mDataType = "rna", cell_lines=cellNames(TGGATESsmall),
    drugs = head(drugNames(TGGATESsmall)),
    features = c(5), duration = "NOTINtSET",
    dose = c("Control", "High"), summary.stat = "first",
    fill.missing = TRUE, verbose=TRUE
    )
  )
})

# dose
context("...Checking for correct dose param errors")
test_that("Errors if given features as type other than character", { expect_error(
  summarizeMolecularProfiles(
    TGGATESsmall, mDataType = "rna", cell_lines=cellNames(TGGATESsmall),
    drugs = head(drugNames(TGGATESsmall)),
    features = c(5), duration = "8",
    dose = c(1, 2), summary.stat = "last",
    fill.missing = TRUE, verbose=TRUE
    )
  )
})
test_that("Errors if specified doses are not in the tSet", { expect_error(
  summarizeMolecularProfiles(
    TGGATESsmall, mDataType = "rna", cell_lines=cellNames(TGGATESsmall),
    drugs = head(drugNames(TGGATESsmall)),
    features = c(5), duration = "8",
    dose = "NOTINTtSET", summary.stat = "median",
    fill.missing = TRUE, verbose=TRUE
    )
  )
})


