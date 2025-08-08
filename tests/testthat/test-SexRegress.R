test_that("SexRegress returns the correct values", {
    data("Regression_Ex", package = "GXwasR")
    fdata <- Regression_Ex
    fdata$SEX <- as.factor(as.character(fdata$SEX))
    response_index <- 1
    regressor_index <- 2
    #'
    x <- SexRegress(fdata, regressor_index, response_index)
    expected <- c(
        Estimate = 0.001305843,
        `Std. Error` = 0.004337332,
        `t value` = 0.301070570,
        `Pr(>|t|)` = 0.763594133
    )
    expect_equal(x, expected, tolerance = 1e-7)
})
