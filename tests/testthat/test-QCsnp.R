test_that("QCsnp returns a list of length 2", {
    skip_on_bioc()
    DataDir <- GXwasR:::GXwasR_data()
    ResultDir <- tempdir()
    finput <- "GXwasR_example"
    foutput <- "Test_output"
    geno <- NULL
    maf <- 0.05
    casecontrol <- FALSE
    hweCase <- NULL
    hweControl <- NULL
    hweCase <- NULL
    monomorphicSNPs <- FALSE
    caldiffmiss <- FALSE
    ld_prunning <- FALSE
    x <- QCsnp(
        DataDir = DataDir, ResultDir = ResultDir, finput = finput, foutput = foutput,
        geno = geno, maf = maf, hweCase = hweCase, hweControl = hweControl,
        ld_prunning = ld_prunning, casecontrol = casecontrol, monomorphicSNPs = monomorphicSNPs,
        caldiffmiss = caldiffmiss
    )
    expect_equal(length(x), 2)
    unlink(ResultDir, recursive = TRUE)
})
