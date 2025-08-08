test_that("executePlinkMAF returns expected output", {
    skip_on_bioc()
    DataDir <- GXwasR:::GXwasR_data()
    ResultDir <- tempdir()
    finput <- "GXwasR_example"
    maf_data <- executePlinkMAF(DataDir, ResultDir, finput)
    expected_maf_data <- data.frame(
        CHR = rep(1, 6),
        SNP = c("rs143773730", "rs75530702", "rs147281566", "rs35854196", "rs115490086", "rs12041521"),
        A1 = c("T", "G", "T", "A", "T", "A"),
        A2 = c("C", "A", "C", "G", "C", "G"),
        MAF = c(0.233700, 0.003623, 0.014490, 0.027170, 0.007246, 0.202900),
        NCHROBS = rep(552, 6)
    )

    expect_equal(head(maf_data), expected_maf_data, tolerance = 1e-4)
    unlink(ResultDir, recursive = TRUE)
})
