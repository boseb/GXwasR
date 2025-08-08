test_that("FilterAllele produces the correct number of files", {
    DataDir <- GXwasR:::GXwasR_data()
    ResultDir <- tempdir()
    finput <- "GXwasR_example"
    foutput <- "Filter_Test"
    x <- FilterAllele(DataDir, ResultDir, finput, foutput)

    expect_equal(list.files(ResultDir, pattern = "^GXwasR_data") %>% length(), 1)
    unlink(ResultDir, recursive = TRUE)
})
