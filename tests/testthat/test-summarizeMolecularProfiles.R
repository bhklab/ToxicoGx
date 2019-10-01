library(testthat)
library(ToxicoGx)

context("Testing if summarizeMolecularProfiles error handling works correctly...")

context("Checking for correct tSet param errors...")
test_that("Errors if given more than one tSet as parameter.", { expect_error(
  summarizeMolecularProfiles(
    c(TGGATESsmall, TGGATESsmall), mDataType = "rna", cell.lines=cellNames(TGGATESsmall),
    drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall,"rna"), duration = "8",
    dose = c("Control", "High"), summary.stat = "median",
    fill.missing = TRUE, verbose=TRUE
    )
  )
})

context("Checking for correct mDataType param errors...")
test_that("Errors if given more than one mDataType.", { expect_error(
  summarizeMolecularProfiles(
    TGGATESsmall, mDataType = c("rna", "cnv"), cell.lines=cellNames(TGGATESsmall),
    drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall,"rna"), duration = "8",
    dose = c("Control", "High"), summary.stat = "median",
    fill.missing = FALSE, verbose=TRUE
    )
  )
})

context("Checking correct tSet param errors...")
test_that("", { expect_error(
  summarizeMolecularProfiles(
    TGGATESsmall, mDataType = "rna", cell.lines=cellNames(TGGATESsmall),
    drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall,"rna"), duration = "8",
    dose = c("Low", "Medium"), summary.stat = "median",
    fill.missing = FALSE, verbose=FALSE
  ))
})

context("Checking correct tSet param errors...")
test_that("Errors if given more than one tSet as parameter", { expect_error(
  summarizeMolecularProfiles(
    c(TGGATESsmall, TGGATESsmall), mDataType = "rna", cell.lines=cellNames(TGGATESsmall),
    drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall,"rna"), duration = "8",
    dose = c("Control", "High"), summary.stat = "median",
    fill.missing = TRUE, verbose=TRUE
  ))
})

context("Checking correct tSet param errors...")
test_that("Errors if given more than one tSet as parameter", { expect_error(
  summarizeMolecularProfiles(
    c(TGGATESsmall, TGGATESsmall), mDataType = "rna", cell.lines=cellNames(TGGATESsmall),
    drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall,"rna"), duration = "8",
    dose = c("Control", "High"), summary.stat = "median",
    fill.missing = TRUE, verbose=TRUE
  ))
})

context("Checking correct tSet param errors...")
test_that("Errors if given more than one tSet as parameter", { expect_error(
  summarizeMolecularProfiles(
    c(TGGATESsmall, TGGATESsmall), mDataType = "rna", cell.lines=cellNames(TGGATESsmall),
    drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall,"rna"), duration = "8",
    dose = c("Control", "High"), summary.stat = "median",
    fill.missing = TRUE, verbose=TRUE
  ))
})

test_that("Errors if given more than one tSet as parameter", { expect_error(
  summarizeMolecularProfiles(
    c(TGGATESsmall, TGGATESsmall), mDataType = "rna", cell.lines=cellNames(TGGATESsmall),
    drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall,"rna"), duration = "8",
    dose = c("Control", "High"), summary.stat = "median",
    fill.missing = TRUE, verbose=TRUE
  ))
})

test_that("Errors if given more than one tSet as parameter", { expect_error(
  summarizeMolecularProfiles(
    c(TGGATESsmall, TGGATESsmall), mDataType = "rna", cell.lines=cellNames(TGGATESsmall),
    drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall,"rna"), duration = "8",
    dose = c("Control", "High"), summary.stat = "median",
    fill.missing = TRUE, verbose=TRUE
  ))
})

## TODO:: Add function return validation
#context("Testing if summarizeMolecularProfiles returns the correct values...")



