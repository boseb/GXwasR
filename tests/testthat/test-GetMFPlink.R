test_that("GetMFPlink creates the correct number of output files", {
  skip_on_ci()
  skip_on_bioc()
  DataDir <- GXwasR:::GXwasR_data()
  ResultDir <- tempdir()
  finput <- "GXwasR_example"
  foutput <- "Test_output"
  sex <- "females"
  x <- GetMFPlink(
      DataDir = DataDir, ResultDir = ResultDir,
      finput = finput, foutput = foutput, sex = sex,
      xplink = FALSE, autoplink = FALSE
  )
  expect_equal(list.files(ResultDir) %>% length(), 5)
  unlink(ResultDir, recursive = TRUE)
})

