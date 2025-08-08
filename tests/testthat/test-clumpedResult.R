test_that("ClumpLD returns expected results", {
  skip_on_bioc()
  data("Summary_Stat_Ex1", package = "GXwasR")
  data("Summary_Stat_Ex2", package = "GXwasR")
  DataDir <- GXwasR:::GXwasR_data()
  ResultDir <- tempdir()
  finput <- "GXwasR_example"
  SNPdata <- list(Summary_Stat_Ex1, Summary_Stat_Ex2)
  clump_p1 <- 0.0001
  clump_p2 <- 0.001
  clump_r2 <- 0.5
  clump_kb <- 250
  byCHR <- TRUE
  clumpedResult <- ClumpLD(
      DataDir, finput, SNPdata, ResultDir, clump_p1,
      clump_p2, clump_r2, clump_kb, byCHR
  )

  expect_type(clumpedResult, "list")
  expect_equal(length(clumpedResult), 2)

  expected <- data.table::data.table(
    INDEX = c("rs6529954", "rs12858640", "rs5962098"),
    PSNP = c("rs6529954", "rs12858640", "rs5962098"),
    RSQ = rep("*", 3),
    KB = c("0", "0", "0"),
    P = c("4.41e-09", "3.89e-09", "1.71e-06"),
    ALLELES = c("AA/GG", "CC/TT", "AA/GG"),
    F = c("2", "2", "2"),
    CHR = c("23", "23", "23")
  )
  expect_equal(clumpedResult$BestClump, expected)
  unlink(ResultDir, recursive = TRUE)
})

