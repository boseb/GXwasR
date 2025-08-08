test_that("LDPrune returns expected output", {
  skip_on_bioc()
  DataDir <- GXwasR:::GXwasR_data()
  ResultDir <- tempdir()
  finput <- "GXwasR_example"
  prunedSNPs <- LDPrune(DataDir, finput, ResultDir, 50, 5, 0.2)
  expect_equal(head(prunedSNPs), c("rs143773730", "rs75530702", "rs147281566", "rs35854196", "rs115490086", "rs12041521"))
  unlink(ResultDir, recursive = TRUE)
})

