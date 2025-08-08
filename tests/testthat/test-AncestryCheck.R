test_that("AncestryCheck returns expected output", {
    skip_on_ci()
    skip_on_bioc()
    data("highLD_hg19", package = "GXwasR")
    data("example_data_study_sample_ancestry", package = "GXwasR")
    DataDir <- GXwasR:::GXwasR_data()
    ResultDir <- tempdir()
    finput <- "GXwasR_example"
    reference <- "HapMapIII_NCBI36"
    highLD_regions <- highLD_hg19
    study_pop <- example_data_study_sample_ancestry # PreimputeEX
    studyLD_window_size <- 50
    studyLD_step_size <- 5
    studyLD_r2_threshold <- 0.02
    filterSNP <- TRUE
    studyLD <- FALSE
    referLD <- FALSE
    referLD_window_size <- 50
    referLD_step_size <- 5
    referLD_r2_threshold <- 0.02
    outlier <- TRUE
    outlier_threshold <- 3
    result <- AncestryCheck(
        DataDir = DataDir, ResultDir = ResultDir, finput = finput,
        reference = reference, highLD_regions = highLD_regions,
        study_pop = study_pop, studyLD = studyLD, referLD = referLD,
        outlierOf = "EUR", outlier = outlier, outlier_threshold = outlier_threshold
    )
    expect_type(result, "list")
    expect_equal(length(result), 4)
    expect_s3_class(result$Outlier_samples, "data.frame")
    expect_s3_class(result$Samples_with_predicted_ancestry, "data.frame")
    expect_s3_class(result$NonOutlier_samples, "data.frame")
    unlink(ResultDir, recursive = TRUE)
})
