test_that("PlinkSummary produces the correct output", {
    DataDir <- system.file("extdata", package = "GXwasR")
    ResultDir <- tempdir()
    finput <- "GXwasR_example"
    #'
    x <- PlinkSummary(DataDir, ResultDir, finput)

    expected_result <- dplyr::tibble(
        n_chr = 12L,
        unique_chr = list(chr_unique = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 23L, 24L)),
        n_snps = 26527L,
        n_samples = 276L
    )
    expect_equal(x, expected_result)
    unlink(ResultDir, recursive = TRUE)
})
