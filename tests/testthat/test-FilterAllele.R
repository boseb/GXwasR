test_that("FilterAllele produces emits informative message when no multi-allelic snps present", {
    DataDir <- system.file("extdata", package = "GXwasR")
    ResultDir <- tempdir()
    finput <- "GXwasR_example"
    foutput <- "Filter_Test"

    expect_message(
        FilterAllele(DataDir, ResultDir, finput, foutput),
        "There are no multi-allelic SNPs present in the input dataset."
    )
    unlink(ResultDir, recursive = TRUE)
})
