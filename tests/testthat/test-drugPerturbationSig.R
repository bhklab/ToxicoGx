library(ToxicoGx)

context("Testing drugPertubrationSig...")

test_that("Summary statistics don't change based on number of drugs", {
  context("Testing summary stats are independent of number of drugs...")
  expect_equal(
    drugPerturbationSig(TGGATESsmall,mDataType = "rna", drugs =  drugNames(TGGATESsmall)[1:3], cells = "Hepatocyte", features = nrf2_module, duration = c("2"), verbose = T, nthread = 3)@.Data[1:3],
    drugPerturbationSig(TGGATESsmall,mDataType = "rna", drugs =  drugNames(TGGATESsmall)[1:5], cells = "Hepatocyte", features = nrf2_module, duration = c("2"), verbose = T, nthread = 3)@.Data[1:3]
    )
})
