test_that("MAFdiffSexControl produces the messaging for example", {
    skip_on_bioc()
    DataDir <- system.file("extdata", package = "GXwasR")
    ResultDir <- tempdir()
    finput <- "GXwasR_example"
    foutput <- "Test_output"
    
    expect_message(MAFdiffSexControl(DataDir, ResultDir, finput, filterSNP = TRUE, foutput = foutput), 'No SNP to be flagged or excluded.')
    unlink(ResultDir, recursive = TRUE)
})
