test_that("Xhwe returns a character vector of length 3", {
    skip_on_ci()
    skip_on_bioc()
    DataDir <- GXwasR:::GXwasR_data()
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
    expect_equal(list.files(ResultDir) %>% length(), 5)
    unlink(ResultDir, recursive = TRUE)
})
