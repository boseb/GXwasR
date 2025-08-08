test_that("SexCheck returns a data.frame", {
  skip_on_bioc()
  DataDir <- GXwasR:::GXwasR_data()
  ResultDir <- tempdir()
  finput <- "GXwasR_example"
  LD <- TRUE
  LD_window_size <- 50
  LD_step_size <- 5
  LD_r2_threshold <- 0.02
  fmax_F <- 0.2
  mmin_F <- 0.8
  impute_sex <- FALSE
  compute_freq <- FALSE
  x <- SexCheck(
      DataDir = DataDir, ResultDir = ResultDir, finput = finput, impute_sex = impute_sex,
      compute_freq = compute_freq, LD_window_size = LD_window_size, LD_step_size = LD_step_size,
      LD_r2_threshold = 0.02, fmax_F = 0.2, mmin_F = 0.8
  )
  
  expect_s3_class(x, "data.frame")
  problematic_sex <- x[x$STATUS != "OK", ]
  rownames(problematic_sex) <- NULL
  expect_equal(problematic_sex, data.frame(
    FID = c("EUR_FIN", "EUR_TSI", "EUR_TSI", "EUR_TSI"),
    IID = c("HG00361", "NA20506", "NA20530", "NA20533"),
    PEDSEX = c(2L, 2L, 2L, 2L),
    SNPSEX = c(0L, 0L, 0L, 0L),
    STATUS = c("PROBLEM", "PROBLEM", "PROBLEM", "PROBLEM"),
    F = c(0.4890, 0.7969, 0.7478, 0.3285)
  ))
  unlink(ResultDir, recursive = TRUE)
})

