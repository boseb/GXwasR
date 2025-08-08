test_that("FilterPlinkSample generates correct output files", {
  skip_on_bioc()
  DataDir <- GXwasR:::GXwasR_data()
  ResultDir <- tempdir()
  finput <- "GXwasR_example"
  foutput <- "casesPlink"
  filter_sample <- "cases"
  keep_remove_sample_file <- "samples_example"
  keep <- FALSE

  FilterPlinkSample(
      DataDir = DataDir, ResultDir = ResultDir,
      finput = finput, foutput = foutput, keep_remove_sample_file = keep_remove_sample_file,
      keep = keep
  )

  expect_equal(length(list.files(ResultDir)), 2)
  unlink(ResultDir, recursive = TRUE)
})

