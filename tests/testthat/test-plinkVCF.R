test_that("plinkVCF creates the correct number of output files", {
  skip_on_bioc()
  finput <- "GXwasR_example" # Plink file
  foutput <- "GXwasR_example1"
  DataDir <- GXwasR:::GXwasR_data()
  ResultDir <- tempdir()
  PtoV <- TRUE
  VtoP <- FALSE
  Famfile <- NULL
  PVbyCHR <- FALSE
  plinkVCF(DataDir, ResultDir, finput, foutput, VtoP, PtoV, Famfile, PVbyCHR)
  expect_equal(list.files(ResultDir) %>% length(), 4)
  unlink(ResultDir, recursive = TRUE)
})
