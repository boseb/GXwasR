test_that("DummyCovar() generates correct output", {
    DataDir <- system.file("extdata", package = "GXwasR")
    ResultDir <- tempdir()
    bfile <- "GXwasR_example"
    incovar <- "covarfile_w_pc_age.txt"
    outcovar <- "dummycovarfile"
    dummy_covars <- DummyCovar(
        DataDir = DataDir, ResultDir = ResultDir,
        bfile = bfile, incovar = incovar,
        outcovar = outcovar
    )
    expected_result <- data.frame(
        FID = c("EUR_FIN", "EUR_FIN", "EUR_FIN", "EUR_FIN", "EUR_FIN", "EUR_FIN"),
        IID = c("HG00171", "HG00173", "HG00174", "HG00176", "HG00177", "HG00178"),
        Categorical = c(-9, -9, -9, -9, -9, -9)
    )
    expect_equal(dummy_covars %>% head(), expected_result)
    unlink(ResultDir, recursive = TRUE)
})
