test_that("SexDiff returns the correct output", {
    data("Mfile", package = "GXwasR")
    data("Ffile", package = "GXwasR")
    Difftest <- SexDiff(Mfile, Ffile)
    expect_s3_class(Difftest, "data.frame")
})
