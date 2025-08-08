test_that("MetaGWAS produces the correct output", {
  skip_on_ci()
  skip_on_bioc()
  data("Summary_Stat_Ex1", package = "GXwasR")
  data("Summary_Stat_Ex2", package = "GXwasR")
  DataDir <- GXwasR:::GXwasR_data()
  ResultDir <- tempdir()
  SummData <- list(Summary_Stat_Ex1, Summary_Stat_Ex2)
  SNPfile <- "UniqueLoci"
  useSNPposition <- FALSE
  UseA1 <- TRUE
  GCse <- TRUE
  byCHR <- FALSE
  pval_filter <- "R"
  top_snp_pval <- 1e-08
  max_top_snps <- 10
  chosen_snps_file <- NULL
  pval_threshold_manplot <- 1e-05
  plotname <- "Meta_Analysis.plot"
  x <- MetaGWAS(
      DataDir = DataDir, SummData = SummData, ResultDir = ResultDir,
      SNPfile = NULL, useSNPposition = TRUE, UseA1 = UseA1, GCse = GCse,
      plotname = "Meta_Analysis.plot", pval_filter, top_snp_pval, max_top_snps,
      chosen_snps_file = NULL, byCHR, pval_threshold_manplot
  )
  expect_type(x, "list")
  expect_equal(length(x), 5)
  expect_equal(x$Resultfixed %>% head(), data.frame(
    CHR = c(1, 1, 1, 1, 1, 1),
    BP = c(73841, 775125, 863863, 928969, 1109154, 1127860),
    SNP = c("rs143773730", "rs147281566", "rs35854196", "rs115490086", "rs12041521", "rs148527527"),
    A1 = c("T", "T", "A", "T", "A", "G"),
    A2 = c("?", "?", "?", "?", "?", "?"),
    Q = c(0.2588, 0.7118, 0.2820, NA, 0.6258, 0.5778),
    I = c(21.57, 0.00, 13.60, NA, 0.00, 0.00),
    P = c(0.5535000, 0.9414000, 0.5093000, 0.6896000, 0.0159200, 0.0008237),
    ES = c(0.1108, -0.0534, 0.3621, 0.5649, -0.5550, 1.3978),
    SE = c(0.1869968, 0.7264282, 0.5486962, 1.4143788, 0.2302193, 0.4179142),
    CI_L = c(-0.2557137, -1.4771993, -0.7133445, -2.2072824, -1.0062299, 0.5786881),
    CI_U = c(0.4773137, 1.3703993, 1.4375445, 3.3370824, -0.1037701, 2.2169119),
    stringsAsFactors = FALSE
  ), tolerance = 1e-4)
  unlink(ResultDir, recursive = TRUE)
})

