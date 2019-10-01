library(ToxicoGx)

# Testing Function Returns
context("Testing drugPertubrationSig...")

test_that("Summary statistics don't change based on number of drugs", {
  context("Testing summary stats are independent of number of drugs...")
  expect_equal(
    as.data.frame(drugPerturbationSig(TGGATESsmall, mDataType = "rna", drugs =  drugNames(TGGATESsmall)[1:3], cells = "Hepatocyte", features = fNames(TGGATESsmall, "rna")[1:200], duration = c("2"), verbose = T, nthread = 3)@.Data)$estimate[,1:3],
    as.data.frame(drugPerturbationSig(TGGATESsmall, mDataType = "rna", drugs =  drugNames(TGGATESsmall)[1:5], cells = "Hepatocyte", features = fNames(TGGATESsmall, "rna")[1:200], duration = c("2"), verbose = T, nthread = 3)@.Data)$estimate[,1:3]
    )
})

# Testing Error Handling Works Correctly
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
