test_that("QCsample returns expected output", {
    skip_on_bioc()
    DataDir <- GXwasR:::GXwasR_data()
    ResultDir <- tempdir()
    finput <- "GXwasR_example"
    foutput <- "Test_output"
    imiss <- 0.01
    het <- 2
    small_sample_mod <- FALSE
    IBD <- 0.2
    IBDmatrix <- FALSE
    ambi_out <- TRUE
    #'
    x <- QCsample(
        DataDir = DataDir, ResultDir = ResultDir, finput = finput,
        foutput = foutput, imiss = imiss, het = het, IBD = IBD,
        ambi_out = ambi_out
    )

    expect_equal(length(x), 8)
    expect_equal(x$IBD_results, data.table::data.table(
        IID1 = "HG00119",
        IID2 = "HG00124",
        PI_HAT = 0.3245
    ))
    unlink(ResultDir, recursive = TRUE)
})
