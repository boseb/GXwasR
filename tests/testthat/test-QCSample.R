test_that("QCsample returns expected output", {
    skip_on_bioc()
    DataDir <- system.file("extdata", package = "GXwasR")
    ResultDir <- tempdir()
    finput <- "GXwasR_example"
    foutput <- "Test_output"
    imiss <- 0.01
    het <- 2
    small_sample_mod <- FALSE
    IBD <- 0.2
    IBDmatrix <- FALSE
    ambi_out <- TRUE
    #'
    x <- QCsample(
        DataDir = DataDir, ResultDir = ResultDir, finput = finput,
        foutput = foutput, imiss = imiss, het = het, IBD = IBD,
        ambi_out = ambi_out
    )

    expect_equal(length(x), 8)
    expect_equal(x$IBD_results, data.table::data.table(
        IID1 = "HG00119",
        IID2 = "HG00124",
        PI_HAT = 0.3245
    ))
    unlink(ResultDir, recursive = TRUE)
})

test_that("QCsample returns foutput when het = NULL", {
    skip_on_bioc()
    ## Use example from preimputationQC vignette to ensure
    ## all necessary intermediate files are present.
    DataDir <- system.file("extdata", package = "GXwasR")
    ResultDir <- tempdir()
    finput <- "GXwasR_example"
    foutput <- "PreimputeEX_QC1"
    x <- FilterAllele(DataDir, ResultDir, finput, foutput)
    foutput <- "PreimputeEX_QC1"
    geno <- 0.2
    maf <- NULL
    casecontrol <- FALSE
    caldiffmiss <- FALSE
    diffmissFilter <- FALSE
    dmissX <- FALSE
    dmissAutoY <- FALSE
    monomorphicSNPs <- TRUE
    ld_prunning <- FALSE
    casecontrol <- FALSE
    hweCase <- NULL
    hweControl <- NULL
    monomorphicSNPs <- FALSE
    ld_prunning <- FALSE
    x <- QCsnp(
        DataDir = DataDir, ResultDir = ResultDir, finput = finput,
        foutput = foutput, geno = geno, maf = maf, hweCase = hweCase,
        hweControl = hweControl, ld_prunning = ld_prunning,
        casecontrol = casecontrol, monomorphicSNPs = monomorphicSNPs,
        caldiffmiss = caldiffmiss, dmissX = dmissX,
        dmissAutoY = dmissAutoY, diffmissFilter = diffmissFilter
    )
    ftemp <- list.files(ResultDir, pattern = "PreimputeEX_QC1", full.names = TRUE)
    file.copy(ftemp, DataDir)
    finput <- "PreimputeEX_QC1"
    foutput <- "PreimputeEX_QC2"
    imiss <- 0.2
    het <- NULL
    IBD <- NULL
    filterSample <- TRUE
    ambi_out <- TRUE
    x <- QCsample(
        DataDir = DataDir, ResultDir = ResultDir, finput = finput, foutput = foutput,
        imiss = imiss, het = het, IBD = NULL, filterSample = filterSample,
        ambi_out = ambi_out
    )
    output_files <- list.files(ResultDir, pattern = foutput)
    expect_equal(length(output_files), 7)
    unlink(ResultDir, recursive = TRUE)
    cleanupDataDir <- list.files(DataDir, pattern = 'PreimputeEX_QC', full.names = TRUE)
    unlink(cleanupDataDir)
})
