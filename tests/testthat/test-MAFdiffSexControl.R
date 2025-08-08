test_that("MAFdiffSexControl produces the correct number of files", {
  skip_on_bioc()
  DataDir <- GXwasR:::GXwasR_data()
  ResultDir <- tempdir()
  finput <- "GXwasR_example"
  foutput <- "Test_output"
  x <- MAFdiffSexControl(DataDir, ResultDir, finput, filterSNP = TRUE, foutput = foutput)

  expect_equal(list.files(ResultDir, pattern = '^GXwasR_data_') %>% length(), 1)

  unlink(ResultDir, recursive = TRUE)
})

