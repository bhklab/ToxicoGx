#context("Testing tSet Class Accessor Methods...")
#
#library(ToxicoGx)
#library(parallel)
#
#test_that("@molecularProfiles slot accessors produce expected results", {
#  data("TGGATESsmall")
#  # External validation
#  expect_equal_to_reference(featureInfo(TGGATESsmall, "rna"), "featureData.TGGATESsmall.rds")
#  expect_equal_to_reference(phenoInfo(TGGATESsmall, "rna"), "phenoData.TGGATESsmall.rds")
#  # Internal validation
#  parallel::mclapply(names(TTGATES_small@molecularProfiles),
#                     function(x) {
#                       expect_equal(molecularProfiles(TGGATESsmall, x),
#                                    TGGATESsmall@molecularProfiles[x])
#                       })
#})
#
#test_that("@cell slot accessors produce expected results", {
#  data("TGGATESsmall")
#  # External validation
#  expect_equal_to_reference(cellInfo(TGGATESsmall), "cell.TGGATESsmall.rds")
#  # Internal validation
#  expect_equal(cellInfo(TGGATESsmall), TGGATESsmall@cell)
#})
#
#test_that("@drug slot accessors produce expected results", {
#  data("TGGATESsmall")
#  # External validation
#  expect_equal_to_reference(drugInfo(TGGATESsmall), "drug.TGGATESsmall.rds")
#  # Internal validation
#})
#
#
#test_that("@annotation slot accessors produce expected results", {
#  data("TGGATESsmall")
#  # External validation
#  expect_equal_to_reference(drugInfo(TGGATESsmall), "drug.TGGATESsmall.rds")
#  # Internal validation
#
#})

#test_that("", {
#  data("TGGATESsmall")
#
#})
#
#test_that("", {
#
#})
#
#test_that("", {
#
#})
