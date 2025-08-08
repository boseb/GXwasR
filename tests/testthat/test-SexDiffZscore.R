test_that("SexDiffZscore returns expected output", {
  data("Example_h2data", package = "GXwasR")
  inputdata <- Example_h2data
  x <- SexDiffZscore(inputdata)
  expected <- data.frame(
    ID = c("ADHD", "AFB", "ALCC", "ALCD", "ANX", "ASD", "BD", "CUE"),
    Mstat = c(0.249, 0.113, 0.086, 0.113, 0.009, 0.227, 0.364, 0.076),
    Mse = c(0.021, 0.010, 0.010, 0.031, 0.003, 0.019, 0.032, 0.010),
    Fstat = c(0.135, 0.052, 0.082, 0.014, 0.005, -0.075, 0.339, 0.067),
    Fse = c(0.028, 0.004, 0.011, 0.052, 0.002, 0.034, 0.027, 0.008),
    Zscore = c(3.2571429, 5.6637078, 0.2690691, 1.6353029, 1.1094004, 7.7537921, 0.5971027, 0.7027819),
    p = c(1.125398e-03, 1.481366e-08, 7.878765e-01, 1.019856e-01, 2.672575e-01, 8.918830e-15, 5.504388e-01, 4.821917e-01),
    adjP = c(9.003182e-03, 1.185093e-07, 1.000000e+00, 8.158845e-01, 1.000000e+00, 7.135064e-14, 1.000000e+00, 1.000000e+00)
  )
  
  expect_equal(x, expected, tolerance = 1e-4)
})

