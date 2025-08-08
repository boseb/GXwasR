test_that("FilterSNP generates correct output files", {
    skip_on_bioc()
    DataDir <- GXwasR:::GXwasR_data()
    ResultDir <- tempdir()
    SNPvec <- c("rs6529954", "rs12858640", "rs5962098")
    finput <- "GXwasR_example"
    foutput <- "Filter_Test"
    FilterSNP(DataDir, ResultDir, finput, foutput, SNPvec = SNPvec, extract = TRUE)

    expect_equal(list.files(ResultDir) %>% length(), 6)

    unlink(ResultDir, recursive = TRUE)
})
