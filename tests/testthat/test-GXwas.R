test_that("GXwas function returns expected results for binary trait with FMcombx02 model", {
  skip_on_bioc()
  DataDir <- GXwasR:::GXwasR_data()
  ResultDir <- tempdir()
  finput <- "GXwasR_example"
  standard_beta <- TRUE
  xsex <- FALSE
  sex <- TRUE
  Inphenocov <- NULL
  covartest <- NULL
  interaction <- FALSE
  MF.na.rm <- FALSE
  B <- 10000
  MF.zero.sub <- 0.00001
  trait <- "binary"
  xmodel <- "FMcombx02"
  combtest <- "fisher.method"
  snp_pval <- 1e-08
  covarfile <- NULL
  ncores <- 0
  MF.mc.cores <- 1
  ResultGXwas <- GXwas(
      DataDir = DataDir, ResultDir = ResultDir,
      finput = finput, xmodel = xmodel, trait = trait, covarfile = covarfile,
      sex = sex, xsex = xsex, combtest = combtest, MF.p.corr = "none",
      snp_pval = snp_pval, plot.jpeg = TRUE, suggestiveline = 5, genomewideline = 7.3,
      MF.mc.cores = 1, ncores = ncores
  )

  expected <- data.table(
    CHR = c(1L, 1L, 1L, 1L, 1L, 1L),
    SNP = c("rs143773730", "rs143773730", "rs75530702", "rs75530702", "rs147281566", "rs147281566"),
    BP = c(73841L, 73841L, 720797L, 720797L, 775125L, 775125L),
    A1 = c("T", "T", "G", "G", "T", "T"),
    TEST = c("ADD", "SEX", "ADD", "SEX", "ADD", "SEX"),
    NMISS = c(276L, 276L, 276L, 276L, 276L, 276L),
    BETA = c(0.1185, 0.2412, 0.4350, 0.2501, -0.05108, 0.25020),
    SE = c(0.1923, 0.2483, 1.4220, 0.2477, 0.7425, 0.2478),
    L95 = c(-0.2584, -0.2456, -2.3530, -0.2355, -1.5060, -0.2355),
    U95 = c(0.4953, 0.7279, 3.2230, 0.7357, 1.4040, 0.7358),
    STAT = c(0.6161, 0.9711, 0.3058, 1.0100, -0.0688, 1.0100),
    P = c(0.5378, 0.3315, 0.7598, 0.3127, 0.9452, 0.3127)
  )

  expect_equal(head(ResultGXwas), expected, tolerance = 1e-4)
  unlink(ResultDir, recursive = TRUE)
})

