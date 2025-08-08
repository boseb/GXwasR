test_that("PvalComb returns expected results", {
  data("Mfile", package = "GXwasR")
  data("Ffile", package = "GXwasR")
  PlotDir <- tempdir()
  SumstatMale <- Mfile
  colnames(SumstatMale)[3] <- "POS"
  SumstatFemale <- Ffile
  colnames(SumstatFemale)[3] <- "POS"
  PvalComb_Result <- PvalComb(
      SumstatMale = SumstatMale, SumstatFemale = SumstatFemale,
      combtest = "fisher.method", MF.mc.cores = 1, snp_pval = 0.001, plot.jpeg = FALSE,
      suggestiveline = 3, genomewideline = 5.69897, ncores = 1
  )
  expected_result <- data.table::data.table(
    SNP = c("rs1015493", "rs1028717", "rs1030166", "rs10416111", "rs10494251", "rs10518607"),
    CHR = c(5, 14, 5, 19, 1, 4),
    POS = c(11359627, 85392037, 140165657, 8970692, 145490518, 133249707),
    P = c(0.13051019, 0.92974281, 0.94368070, 0.83297152, 0.05014729, 0.10323853)
  )
  data.table::setkey(expected_result, SNP)

  expect_equal(head(PvalComb_Result), expected_result)
  expect_equal(length(list.files(PlotDir, pattern = '^Stratified')), 2)
  unlink(PlotDir)
})

