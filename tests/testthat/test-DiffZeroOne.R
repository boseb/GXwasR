test_that("DiffZeroOne returns expected output", {
    data("Example_rgdata", package = "GXwasR")
    inputdata <- Example_rgdata
    colnames(inputdata) <- c("Trait", "Stat", "SE")
    x <- DiffZeroOne(inputdata, FALSE, TRUE)

    expected <- data.frame(
        Trait = c("RTB", "EA", "BD", "CUE", "NEU", "SCZ", "AFB", "SMKP", "INS", "ANX", "ALCC", "NEB", "MDD", "OCD", "SMKC", "ADHD"),
        Stat = c(0.811, 0.919, 0.862, 0.769, 0.928, 0.923, 0.886, 0.915, 0.786, 0.678, 0.911, 0.903, 0.996, 1.034, 1.022, 1.209),
        SE = c(0.041, 0.022, 0.057, 0.091, 0.029, 0.033, 0.057, 0.044, 0.125, 0.275, 0.090, 0.097, 0.305, 0.564, 0.050, 0.127),
        p1 = c(2.015708e-06, 1.157883e-04, 7.737818e-03, 5.567052e-03, 6.518470e-03, 9.815329e-03, 2.275013e-02, 2.669098e-02, 4.344833e-02, 1.208177e-01, 1.613588e-01, 1.586553e-01, 4.947681e-01, 5.240352e-01, 6.700314e-01, 9.500841e-01)
    )

    expect_equal(x, expected, tolerance = 1e-4)
})
