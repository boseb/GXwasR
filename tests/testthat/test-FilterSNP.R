test_that("FilterSNP generates correct output files", {
    skip_on_bioc()
    DataDir <- system.file("extdata", package = "GXwasR")
    ResultDir <- tempdir()
    SNPvec <- c("rs6529954", "rs12858640", "rs5962098")
    finput <- "GXwasR_example"
    foutput <- "Filter_Test"
    FilterSNP(DataDir, ResultDir, finput, foutput, SNPvec = SNPvec, extract = TRUE)

    expect_equal(list.files(ResultDir, pattern = "^Filter") %>% length(), 4)

    unlink(ResultDir, recursive = TRUE)
})
