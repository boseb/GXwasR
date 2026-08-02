test_that("Xhwe returns a character vector of length 3", {
    DataDir <- system.file("extdata", package = "GXwasR")
    ResultDir <- tempdir()
    finput <- "GXwasR_example"
    foutput <- "Test_output"
    x <- Xhwe(
        DataDir = DataDir, ResultDir = ResultDir,
        finput = finput, foutput = foutput, filterSNP = TRUE
    )
    expect_type(x, "character")
    expect_equal(length(x), 3)
    expect_equal(x, c("rs56053951", "rs12353847", "rs5940058"))
    unlink(ResultDir, recursive = TRUE)
})
