test_that("MergeRegion creates the correct number of output files", {
    skip_on_bioc()
    DataDir <- GXwasR:::GXwasR_data()
    ResultDir <- tempdir()
    finput1 <- "GXwasR_example"
    finput2 <- "GXwasR_example_imputed"
    foutput <- "Test_output"
    y <- MergeRegion(DataDir, ResultDir, finput1, finput2, foutput, use_common_snps = TRUE)
    expect_equal(list.files(ResultDir, pattern = "^Test_output") %>% length(), 4)
    unlink(ResultDir, recursive = TRUE)
})
