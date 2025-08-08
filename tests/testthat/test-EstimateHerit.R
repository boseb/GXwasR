test_that("EstimateHerit returns expected output", {
  skip_on_ci()
  skip_on_bioc()
  data("Summary_Stat_Ex1", package = "GXwasR")
  data("highLD_hg19", package = "GXwasR")
  DataDir <- GXwasR:::GXwasR_data()
  ResultDir <- tempdir()
  precomputedLD <- NULL
  finput <- "GXwasR_example"
  test.sumstats <- na.omit(Summary_Stat_Ex1[Summary_Stat_Ex1$TEST == "ADD", c(seq_len(4), 6:8)])
  colnames(test.sumstats) <- c("chr", "rsid", "pos", "a1", "n_eff", "beta", "beta_se")
  summarystat <- test.sumstats
  ncores <- 3
  model <- "GREML"
  byCHR <- FALSE
  r2_LD <- 0
  LDSC_blocks <- 20
  REMLalgo <- 0
  nitr <- 3
  cat_covarfile <- NULL
  quant_covarfile <- NULL
  prevalence <- 0.01
  partGRM <- FALSE
  autosome <- TRUE
  Xsome <- TRUE
  nGRM <- 3
  cripticut <- 0.025
  minMAF <- NULL
  maxMAF <- NULL
  hg <- "hg19"
  PlotIndepSNP <- TRUE
  IndepSNP_window_size <- 50
  IndepSNP_step_size <- 5
  IndepSNP_r2_threshold <- 0.02
  highLD_regions <- highLD_hg19
  H2 <- EstimateHerit(
      DataDir = DataDir, ResultDir = ResultDir, finput = finput,
      summarystat = NULL, ncores, model = "GREML", byCHR = TRUE, r2_LD = 0,
      LDSC_blocks = 20, REMLalgo = 0, nitr = 100, cat_covarfile = NULL, quant_covarfile = NULL,
      prevalence = 0.01, partGRM = FALSE, autosome = TRUE, Xsome = TRUE, nGRM = 3,
      cripticut = 0.025, minMAF = NULL, maxMAF = NULL, hg = "hg19", PlotIndepSNP = TRUE,
      IndepSNP_window_size = 50, IndepSNP_step_size = 5, IndepSNP_r2_threshold = 0.02,
      highLD_regions = highLD_hg19
  )
  expect_s3_class(H2, "data.frame")
  unlink(ResultDir, recursive = TRUE)
})
