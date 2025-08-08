test_that("Miami plot is generated", {
  skip_on_ci()
  skip_on_bioc()
  data("Ffile", package = "GXwasR")
  data("Mfile", package = "GXwasR")
  ResultDir <- tempdir()
  FemaleWAS <- na.omit(Ffile[, c("SNP", "CHR", "BP", "P")])
  colnames(FemaleWAS) <- c("SNP", "CHR", "POS", "pvalue")
  MaleWAS <- na.omit(Mfile[, c("SNP", "CHR", "BP", "P")])
  colnames(MaleWAS) <- c("SNP", "CHR", "POS", "pvalue")
  #'
  GXWASmiami(ResultDir = ResultDir, FemaleWAS = FemaleWAS, MaleWAS = MaleWAS, snp_pval = 0.05)
  expect_equal(length(list.files(ResultDir, pattern = '^Stratified')), 1)
  unlink(ResultDir, recursive = TRUE)
})

