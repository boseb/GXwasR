test_that("PlinkSummary produces the correct number of files", {
    skip_on_bioc()
    DataDir <- GXwasR:::GXwasR_data()
    ResultDir <- tempdir()
    finput <- "GXwasR_example"
    #'
    x <- PlinkSummary(DataDir, ResultDir, finput)

    expect_equal(list.files(ResultDir, pattern = "^GXwasR_data") %>% length(), 1)
    unlink(ResultDir, recursive = TRUE)
})
