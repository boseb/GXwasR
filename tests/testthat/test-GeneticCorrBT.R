test_that("GeneticCorrBT returns expected output", {
  skip_on_bioc()
  data("Example_phenofile", package = "GXwasR")
  DataDir <- GXwasR:::GXwasR_data()
  ResultDir <- tempdir()
  finput <- "GXwasR_example"
  byCHR <- TRUE
  REMLalgo <- 0
  nitr <- 3
  ncores <- 3
  phenofile <- Example_phenofile # Cannot be NULL
  cat_covarfile <- NULL
  quant_covarfile <- NULL
  partGRM <- FALSE # Partition the GRM into m parts (by row),
  autosome <- TRUE
  Xsome <- TRUE
  cripticut <- 0.025
  minMAF <- 0.01 # if MAF filter apply
  maxMAF <- 0.04
  excludeResidual <- TRUE
  #'
  genetic_correlation <- GeneticCorrBT(
      DataDir = DataDir, ResultDir = ResultDir, finput = finput, byCHR = byCHR,
      REMLalgo = 0, nitr = 10, phenofile = phenofile, cat_covarfile = NULL, quant_covarfile = NULL,
      partGRM = FALSE, autosome = TRUE, Xsome = TRUE, nGRM = 3,
      cripticut = 0.025, minMAF = NULL, maxMAF = NULL, excludeResidual = TRUE, ncores = ncores
  )

  expected_genetic_correlation <- data.table(
    chromosome = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2),
    Source = c("V.G._tr1", "V.G._tr2", "C.G._tr12", "V.e._tr1", "V.e._tr2", "X", "V.G._tr1", "V.G._tr2", "C.G._tr12", "V.e._tr1"),
    Variance = c(0, 0, 0, 0.23783, 0.02517, NA, 0, 0, 0, 0.23880),
    SE = c(NA_real_)
  )
  
  expect_equal(head(genetic_correlation, 10), expected_genetic_correlation, tolerance = 1e-6)
  unlink(ResultDir, recursive = TRUE)
})
