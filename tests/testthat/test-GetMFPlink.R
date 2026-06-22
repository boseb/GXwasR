test_that("GetMFPlink creates the correct number of output files", {
    skip_on_ci()
    skip_on_bioc()
    DataDir <- system.file("extdata", package = "GXwasR")
    ResultDir <- tempdir()
    finput <- "GXwasR_example"
    foutput <- "Test_output"
    sex <- "females"
    x <- GetMFPlink(
        DataDir = DataDir, ResultDir = ResultDir,
        finput = finput, foutput = foutput, sex = sex,
        xplink = FALSE, autoplink = FALSE
    )
    expect_equal(list.files(ResultDir, pattern = '^Test_output') %>% length(), 4)
    unlink(ResultDir, recursive = TRUE)
})
