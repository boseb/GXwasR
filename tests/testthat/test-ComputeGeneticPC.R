test_that("ComputeGeneticPC returns expected output", {
    skip_on_bioc()
    data("highLD_hg19", package = "GXwasR")
    DataDir <- GXwasR:::GXwasR_data()
    ResultDir <- tempdir()
    finput <- "GXwasR_example"
    highLD_regions <- highLD_hg19
    ld_prunning <- "TRUE"
    window_size <- 50
    step_size <- 5
    r2_threshold <- 0.02
    countPC <- 20
    ## Genetic PC
    GP <- ComputeGeneticPC(
        DataDir = DataDir, ResultDir = ResultDir,
        finput = finput, highLD_regions = highLD_hg19, countPC = 20
    )
    expect_type(GP, "list")
    expect_equal(length(GP), 2)
    expect_s3_class(GP$PCs1, "data.frame")
    unlink(ResultDir, recursive = TRUE)
})
