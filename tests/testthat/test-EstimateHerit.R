test_that("EstimateHerit returns expected output", {
    withr::local_envvar(
        GENEINFO_HG19 = system.file(
            "extdata",
            "GXwasR_geneinfo_hg19_test.txt",
            package = "GXwasR"
        )
    )
    data("Summary_Stat_Ex1", package = "GXwasR")
    data("highLD_hg19", package = "GXwasR")
    data("PrecomputedLD_Ex1", package = "GXwasR")

    test.sumstats <- na.omit(
        Summary_Stat_Ex1[
            Summary_Stat_Ex1$TEST == "ADD",
            c(seq_len(4), 6:8)
        ]
    )

    colnames(test.sumstats) <- c(
        "chr", "rsid", "pos", "a1",
        "n_eff", "beta", "beta_se"
    )

    H2ldsc <- EstimateHerit(
        summarystat = test.sumstats,
        precomputedLD = PrecomputedLD_Ex1,
        ncores = 1,
        model = "LDSC",
        byCHR = TRUE,
        hg = "hg19",
        plotjpeg = FALSE
    )

    expect_s3_class(H2ldsc, "data.frame")
})
