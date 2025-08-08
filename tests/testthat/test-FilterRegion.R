test_that("FilterRegion returns expected output", {
    skip_on_bioc()
    DataDir <- GXwasR:::GXwasR_data()
    ResultDir <- tempdir()
    finput <- "GXwasR_example"
    foutput <- "PostimputeEX_QC1"
    x <- FilterRegion(
        DataDir = DataDir, ResultDir = ResultDir,
        finput = finput, foutput = foutput, CHRX = TRUE, CHRY = FALSE,
        filterPAR = TRUE, filterXTR = TRUE, filterAmpliconic = TRUE,
        regionfile = FALSE, filterCHR = NULL, Hg = "38", exclude = TRUE
    )
    expect_type(x, "list")
    expect_equal(length(x), 3)
    expected_output <- data.frame(
        CHR = rep(23, 9),
        SNP = c("rs62602496", "rs6612314", "rs151231489", "rs73498395", "rs6527", "rs4907822", "rs142219143", "rs5987512", "rs139353379"),
        START = rep(0, 9),
        END = c(48349540, 55442087, 55445831, 71740842, 73066891, 102264099, 102354248, 102483549, 103970756),
        A1 = c("G", "A", "C", "T", "A", "A", "G", "G", "G"),
        A2 = c("A", "C", "A", "C", "C", "G", "A", "C", "A")
    )
    expect_equal(x$Ampliconic, expected_output)
    unlink(ResultDir, recursive = TRUE)
})
