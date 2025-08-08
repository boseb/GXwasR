test_that("ComputeLD function returns expected output", {
  skip_on_bioc()
  ResultDir <- tempdir()
  snpld <- ComputeLD(
    DataDir = GXwasR:::GXwasR_data(), ResultDir = ResultDir,
    finput = "GXwasR_example", ByCHR = TRUE, CHRnum = 1, r2_LD = 0.2
  )
  expected <- data.frame(
    CHR_A = c(1L, 1L, 1L, 1L, 1L, 1L),
    BP_A = c(1617443L, 1821625L, 2793754L, 2838513L, 4520285L, 4558972L),
    SNP_A = c("rs112033089", "rs6603805", "rs113572417", "rs78083445", "rs6426412", "rs13374152"),
    CHR_B = c(1L, 1L, 1L, 1L, 1L, 1L),
    BP_B = c(1711414L, 1875530L, 2806390L, 2845752L, 4549482L, 4580135L),
    SNP_B = c("rs867207", "rs2803333", "rs112163476", "rs10797367", "rs667500", "rs13375764"),
    R2 = c(0.247273, 0.220796, 1, 0.476771, 0.228998, 0.249483)
  )
  expect_equal(head(snpld), expected)
  unlink(ResultDir, recursive = TRUE)
})

