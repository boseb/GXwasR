test_that("plinkVCF creates the correct number of output files", {
    finput <- "GXwasR_example" # Plink file
    foutput <- "GXwasR_example1"
    DataDir <- system.file("extdata", package = "GXwasR")
    ResultDir <- tempdir()
    PtoV <- TRUE
    VtoP <- FALSE
    Famfile <- NULL
    PVbyCHR <- FALSE
    plinkVCF(DataDir, ResultDir, finput, foutput, VtoP, PtoV, Famfile, PVbyCHR)
    expect_equal(list.files(ResultDir, pattern = "^GXwasR_example") %>% length(), 3)
    unlink(ResultDir, recursive = TRUE)
})
