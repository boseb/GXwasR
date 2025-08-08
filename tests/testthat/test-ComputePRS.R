test_that("ComputePRS returns expected results", {
    skip_on_bioc()
    data("Summary_Stat_Ex1", package = "GXwasR")
    data("Example_phenofile", package = "GXwasR")
    data("Example_covarfile", package = "GXwasR")
    data("Example_pthresoldfile", package = "GXwasR")
    data("highLD_hg19", package = "GXwasR")
    DataDir <- GXwasR:::GXwasR_data()
    ResultDir <- tempdir()
    finput <- "GXwasR_example"
    summarystat <- Summary_Stat_Ex1[, c(2, 4, 7, 1, 3, 12)]
    phenofile <- Example_phenofile # Cannot be NULL
    covarfile <- Example_covarfile
    clump_p1 <- 0.0001
    clump_p2 <- 0.0001
    clump_kb <- 500
    clump_r2 <- 0.5
    byCHR <- TRUE
    pthreshold <- Example_pthresoldfile$Threshold
    ld_prunning <- TRUE
    highLD_regions <- highLD_hg19
    window_size <- 50
    step_size <- 5
    r2_threshold <- 0.02
    nPC <- 6
    pheno_type <- "binary"

    PRSresult <- ComputePRS(DataDir, ResultDir, finput, summarystat, phenofile, covarfile,
        effectsize = "BETA", LDreference = "GXwasR_example", ldclump = FALSE, clump_p1, clump_p2,
        clump_r2, clump_kb, byCHR = TRUE, pthreshold = pthreshold, highLD_regions = highLD_regions,
        ld_prunning = TRUE, window_size = 50, step_size = 5, r2_threshold = 0.02, nPC = 6,
        pheno_type = "binary"
    )

    expect_type(PRSresult, "list")
    expect_equal(length(PRSresult), 4)

    PRS <- PRSresult$PRS
    expect_s3_class(PRS, "data.frame")

    BestPvalue <- PRSresult$BestP$Threshold
    expect_equal(BestPvalue, 0.05)
    unlink(ResultDir, recursive = TRUE)
})
