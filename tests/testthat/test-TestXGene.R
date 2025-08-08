test_that("TestXGene returns expected output", {
  skip_on_bioc()
  skip_on_ci()
  data("XWAS_Summary_Example", package = "GXwasR")
  DataDir <- GXwasR:::GXwasR_data()
  ResultDir <- tempdir()
  finput <- "GXwasR_example"
  sumstat <- XWAS_Summary_Example
  ref_data <- NULL
  gene_file <- "Xlinkedgenes_hg19.txt"
  gene_range <- 500000
  max_gene <- 10
  gene_approximation <- TRUE
  beta_par <- c(1, 25)
  weights_function <- NULL
  geno_variance_weights <- "se.beta"
  method <- "kuonen"
  acc_devies <- 1e-8
  lim_devies <- 1e+6
  rho <- TRUE
  skato_p_threshold <- 0.8
  mac_threshold <- 3
  sample_size <- 4000
  reference_matrix_used <- FALSE
  regularize_fun <- "LH"
  pca_var_fraction <- 0.85
  flm_basis_function <- "fourier"
  flm_num_basis <- 25
  flm_poly_order <- 4
  flip_genotypes <- FALSE
  omit_linear_variant <- FALSE
  kernel_p_method <- "kuonen"
  anno_type <- ""
  GenetestResult <- TestXGene(DataDir, ResultDir, finput, sumstat, gene_file,
      gene_range, score_file, ref_data, max_gene, sample_size,
      genebasedTest = "SKAT",
      gene_approximation, beta_par, weights_function, geno_variance_weights,
      kernel_p_method, acc_devies, lim_devies, rho, skato_p_threshold, anno_type,
      mac_threshold, reference_matrix_used, regularize_fun, pca_var_fraction,
      flm_basis_function, flm_num_basis, flm_poly_order, flip_genotypes,
      omit_linear_variant
  )
  rownames(GenetestResult) <- NULL
  expected_output <- data.frame(
    gene = c("uc004cqy.3", "uc004crr.3", "uc004csf.3", "uc004csr.3", "uc004cst.2", "uc004csu.1", "uc004csx.4", "uc004cti.4", "uc004cup.1", "uc004cuw.3"),
    chrom = rep(23, 10),
    start = c(2822010, 5808082, 8496914, 9433200, 9693452, 9754495, 9983794, 10413349, 11155662, 11776277),
    end = c(2847416, 6146706, 8700227, 9687780, 9734005, 9917481, 10112518, 10851809, 11683821, 11793872),
    markers = c(0, 1, 0, 0, 0, 0, 0, 1, 1, 0),
    filtered.markers = c(0, 1, 0, 0, 0, 0, 0, 1, 1, 0),
    pvalue = c(NA, 0.9365, NA, NA, NA, NA, NA, 0.1846, 0.9751, NA)
  )
  expect_s3_class(GenetestResult, "data.frame")
  expect_equal(GenetestResult, expected_output)
  unlink(ResultDir, recursive = TRUE)
})

