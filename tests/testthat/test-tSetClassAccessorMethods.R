context("Testing tSet Class Accessor Methods...")

library(ToxicoGx)
library(parallel)
data("TGGATES_small")

test_that("@molecularProfiles slot accessors produce expected results", {
  # External validation
  #expect_equal_to_reference(featureInfo(TGGATES_small, "rna"), "featureData.TGGATES_small.rds")
  #expect_equal_to_reference(phenoInfo(TGGATES_small, "rna"), "phenoData.TGGATES_small.rds")
  # Internal validation
  parallel::mclapply(names(TTGATES_small@molecularProfiles),
                     function(x) {
                       expect_equal(molecularProfiles(TGGATES_small, x),
                                    TGGATES_small@molecularProfiles[x])
                       })
})

test_that("@cell slot accessors produce expected results", {
  # External validation
  expect_equal_to_reference(cellInfo(TGGATES_small), "cell.TGGATES_small.rds")
  # Internal validation
  expect_equal(cellInfo(TGGATES_small), TGGATES_small@cell)
})

test_that("@drug slot accessors produce expected results", {
  # External validation
  expect_equal_to_reference(drugInfo(TGGATES_small), "drug.TGGATES_small.rds")
  # Internal validation
})


#test_that("", {
#
#})

#test_that("", {
#
#})

#test_that("", {
#
#})

#test_that("", {
#
#})
