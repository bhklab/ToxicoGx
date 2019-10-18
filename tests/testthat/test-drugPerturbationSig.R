library(ToxicoGx)

## TODO:: Specify what the error messages should be to improve robustness of tests

# Testing Function Returns
context("Testing drugPertubrationSig...")

test_that("Summary statistics don't change based on number of drugs", {
  context("...Testing summary stats are independent of number of drugs...")
  expect_equal(
    suppressWarnings(drugPerturbationSig(TGGATESsmall, mDataType = "rna", drugs =  drugNames(TGGATESsmall)[1:3], cell.lines = "Hepatocyte", features = fNames(TGGATESsmall, "rna")[1:200], duration = c("2"), verbose = T, nthread = 2)@.Data[1,1:2,1]),
    suppressWarnings(drugPerturbationSig(TGGATESsmall, mDataType = "rna", drugs =  drugNames(TGGATESsmall)[1:5], cell.lines = "Hepatocyte", features = fNames(TGGATESsmall, "rna")[1:200], duration = c("2"), verbose = T, nthread = 2)@.Data[1,1:2,1])
    )
})

# Testing Error Handling Works Correctly

context("...Testing if drugPerturbationSig error handling works correctly...")

# tSet
context("...Checking for correct tSet param errors...")
test_that("Errors if given more than one tSet as parameter.", { expect_error(
  drugPerturbationSig(
    c(TGGATESsmall, TGGATESsmall), mDataType="rna", nthread=1,
    duration = "2", drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall, "rna")[1:2]
    )
  )
})

# mDataTypes
context("...Checking for correct mDataType param errors...")
test_that("Errors if given more than one mDataType.", { expect_error(
  drugPerturbationSig(
    TGGATESsmall, mDataType=c("rna", "cnv"), nthread=1,
    duration = "2", drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall, "rna")[1:2]
  )
)
})
test_that("Errors if given mDataType as type other than character.", { expect_error(
  drugPerturbationSig(
    TGGATESsmall, mDataType=c(TRUE, FALSE), nthread=1,
    duration = "2", drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall, "rna")[1:2]
  )
)
})
test_that("Errors if specified mDataType is not in the tSet.", { expect_error(
  drugPerturbationSig(
    TGGATESsmall, mDataType="cnv", nthread=1,
    duration = "2", drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall, "rna")[1:2]
  )
)
})

# cell.lines
context("...Checking for correct cells param errors...")
test_that("Errors if given cell.lines as type other than character.", { expect_error(
  drugPerturbationSig(
    TGGATESsmall, mDataType="rna", nthread=1, cell.lines = c(TRUE, FALSE),
    duration = "2", drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall, "rna")[1:2]
  )
)
})
test_that("Errors if specified cells are not in the tSet", { expect_error(
  drugPerturbationSig(
    TGGATESsmall, mDataType="rna", nthread=1, cell.lines = "NOTINtSET",
    duration = "2", drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall, "rna")[1:2]
  )
)
})

# drugs
context("...Checking for correct drugs param errors...")
test_that("Errors if given drugs are type other than character.", { expect_error(
  drugPerturbationSig(
    TGGATESsmall, mDataType="rna", nthread=1,
    duration = "2", drugs = c(1, 2),
    features = fNames(TGGATESsmall, "rna")[1:2]
  )
)
})
test_that("Errors if specified drugs are not in the tSet.", { expect_error(
  drugPerturbationSig(
    TGGATESsmall, mDataType="rna", nthread=1,
    duration = "2", drugs = "NOTINtSET",
    features = fNames(TGGATESsmall, "rna")[1:2]
    )
  )
})

# features
context("...Checking for correct features param errors....")
test_that("Errors if given less than two features", { expect_error(
  drugPerturbationSig(
    TGGATESsmall, mDataType="rna", nthread=1,
    duration = "2", drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall, "rna")[1]
    )
  )
})
test_that("Errors if given features as type other than character", { expect_error(
  drugPerturbationSig(
    TGGATESsmall, mDataType="rna", nthread=1,
    duration = "2", drugs = head(drugNames(TGGATESsmall)),
    features = c(1, 2)
    )
  )
})
test_that("Errors given features are not in the tSet", { expect_error(
  drugPerturbationSig(
    TGGATESsmall, mDataType="rna", nthread=1,
    duration = "2", drugs = head(drugNames(TGGATESsmall)),
    features = "NOTINtSET"
    )
  )
})

# duration
context("...Checking for correct duration param errors")
test_that("Errors if given duration as type other than character", { expect_error(
  drugPerturbationSig(
    TGGATESsmall, mDataType="rna", nthread=1,
    duration = 2, drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall, "rna")[1:2]
    )
  )
})
test_that("Errors if specified duration is not in the tSet", { expect_error(
  drugPerturbationSig(
    TGGATESsmall, mDataType="rna", nthread=1,
    duration = "NOTINTtSET", drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall, "rna")[1:2]
    )
  )
})

# dose
context("...Checking for correct dose param errors")
test_that("Errors if given dose as type other than character", { expect_error(
  drugPerturbationSig(
    TGGATESsmall, mDataType = "rna", nthread = 1,
    duration = "2", drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall, "rna")[1:2], dose = c(1, 2)
    )
  )
})
test_that("Errors if specified doses are not in the tSet", {expect_error(
  drugPerturbationSig(
    TGGATESsmall, mDataType="rna", nthread = 1,
    duration = "2", drugs = head(drugNames(TGGATESsmall)),
    features = fNames(TGGATESsmall, "rna")[1:2], dose = c("Too Much", "Lots")
    )
  )
})
