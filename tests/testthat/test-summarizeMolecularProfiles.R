context("Tesing summarizeMolecularProfiles function...")

test_that("Error handling works correctly...", {
  data
  expect_error(summarizeMolecularProfiles())
})
