#' TestXGene: Performing gene-based association test using GWAS/XWAS summary statistics.
#'
#' @author Banabithi Bose
#'
#' @description
#' This function performs gene-based association tests using GWAS/XWAS summary statistics and SNP-SNP correlation
#' matrices. For  SNP-SNP correlation matrices, users have the flexibility to use either the base genotype data or 1000 Genomes
#' Phase 3 reference genotype data. Users also have options to define the regional positions of genes to include the SNPs according
#' to their investigation.
#'
#' This function computes gene-wise SNP-SNP correlation matrices and can perform nine different gene-based tests, such as, “BT" (burden test),
#' "SKAT" (sequence kernel association test), "SKATO" (combination of BT and SKAT), "sumchi" (sum of \eqn{\chi^2}-statistics), "ACAT" (aggregated
#' Cauchy association test for combining P values), "PCA"(principal component approach), "FLM"( functional multiple linear regression model),
#' "simpleM" (Bonferroni correction test), "minp" (minimum P-value) leveraging PLINK1.9 \insertCite{Purcell2007}{GXwasR} and sumFREGAT
#' \insertCite{Svishcheva2019,Belonogova2022}{GXwasR} tools.
#'
#' Though this function implicitly performs X-linked gene-based test, it is flexible to perform this analysis genome-wide.
#' For the details about the different tests, please follow the associated paper.
#'
#'
#' @param DataDir
#' A character string for the file path of the all the input files.
#'
#' @param ResultDir
#' A character string for the file path where all output files will be stored. The default is `tempdir()`.
#'
#' @param finput
#' Character string, specifying the prefix of the input PLINK binary files for the genotype data. This file is used to compute the
#' correlation between the SNPs. This file needs to be in DataDir. If the base genotype data is unavailable, then users can use the
#' 1000 Genomes Project samples. Users should use the population that most closely represents the base sample. For ACAT model, this
#' parameter is not mandatory and could be set `NULL`.
#'
#' @param sumstat
#' A dataframe object with GWAS summary statistics. When the base-genotype data is used to compute genetic correlations, the mandatory
#' columns are:
#'
#' * Column 1: `CHROM` (i.e., chromosome code),
#' * Column 2: `POS` (i.e., base-pair position),
#' * Column 3: `ID` (i.e. SNP IDs),
#' * Column 4: `P` (i.e., p-values),
#' * Column 5: `BETA` (i.e., effect-size),
#' * Column 6: `A1` (i.e., effect allele),
#' * Column 7: `A2` (i.e., alternative allele) and
#' * Column 8: `EAF` (i.e., the effect allele frequency)
#'
#' These are mandatory when base-genotype data is used to compute genetic correlations. Otherwise, if the users are using reference data,
#' then columns 5 to 8 are optional. Also, in that case, columns, such as `REF` (i.e., reference allele), and `ALT` (i.e., alternative allele)
#' could be present to compare alleles with those in the reference file and exclude genetic variants if alleles do not match. There could be an
#' additional column, `ANNO` with functional annotations (like "intron_variant", "synonymous", "missense" etc.)
#'
#' @param gene_file
#' Character string, specifying the prefix of the name of a .txt file listing genes in refFlat format. This file needs to be in `DataDir`. The
#' X-linked gene files, "Xlinkedgenes_hg19.txt" and "Xlinkedgenes_hg38.txt" and autosomal gene files, “Autosomes_hg19.txt” and “Autosomes_hg38.txt”
#' can be specified. The default is "Xlinkedgenes_hg19.txt". The genome built should be in agreement with the analysis.
#'
#' @param gene_range
#' Integer value, specifying the up stream and down stream range (in kilo base) of a gene for SNPs to be considered. The default is 500000.
#'
#' @param score_file
#' Character string, specifying the prefix of a file which will be used to produce score files with Z scores from P values and beta input
#' from GWAS summary statistics.
#'
#' @param ref_data
#' Character string, specifying the path to a reference dataframe with additional data needed to recode user data according to correlation matrices
#' that will be used. It contains `ID` column with names of  SNPs,  `REF` and `ALT` columns with alleles that were coded as 0 and 1, respectively.
#' Effect sizes from data will be inverted for variants with effect alleles different from `ALT` alleles in reference data. If presented, `REF` and
#' `ALT` columns from the input data will be used to sort out variants with alleles different from those in reference data. This dataframe  can also
#' be a source of map data and allele frequencies if they are not present in data. `AF` column in the reference file represents the allele frequency
#' of `ALT` allele. The default is "ref1KG.MAC5.EUR_AF.RData".
#'
#' @param max_gene
#' Positive integer value, specifying the number of genes for which the gene-based test will be performed. The default is NULL to consider
#' all the genes.
#'
#' @param sample_size
#' Positive integer value, specifying the sample size of the GWAS. Only needed for FLM and PCA models.
#'
#' @param genebasedTest
#' Character string, specifying the name of the gene-based test. Nine different tests can be specified, "SKAT","SKATO","sumchi","ACAT","BT","PCA",
#' "FLM","simpleM","minp". The default is "SKAT".
#'
#' @param beta_par
#' Boolean value, `TRUE` or `FALSE`, specifying whether approximation for large genes (>= 500 SNPs) should be used. Applicable for SKAT, SKATO,
#' sumchi, PCA, FLM (default = `TRUE` for these methods).
#'
#' @param weights_function
#' A function of MAF to assign weights for each genetic variant. By default is `NULL`. In this case the weights will be calculated using
#' the beta distribution.
#'
#' @param geno_variance_weights
#' Character string, indicating whether scores should be weighted by the variance of genotypes: "none" (i.e., no weights applied, resulting
#' in a sum chi-square test); "se.beta" (i.e., scores weighted by variance of genotypes estimated from P values and effect sizes); "af"
#' (i.e., scores weighted by variance of genotypes calculated as \eqn{AF * (1 - AF)}, where AF is allele frequency.
#'
#' @param kernel_p_method
#' Character string, specifying the method for computing P value in kernel-based tests, such as SKAT, SKATO and sumchi. Available methods
#' are "kuonen" \insertCite{Belonogova2022}{GXwasR} "davies" \insertCite{Belonogova2022}{GXwasR} and "hybrid" \insertCite{Belonogova2022}{GXwasR}.
#' The default is "kuonen".
#'
#' @param acc_devies
#' Positive numeric value, specifying the accuracy parameter for "davies" method. The default is 1e-8.
#'
#' @param lim_devies
#' Positive numeric value, specifying the limit parameter for "davies" method. The default is 1e+6.
#'
#' @param rho
#' Logical value, 'TRUE' or 'FALSE' or can be a vector of grid values from 0 to 1. If TRUE, the optimal test (SKAT-O) is performed (12).
#' The default grid is c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1).
#'
#' @param skato_p_threshold
#' Positive numeric value, specifying the largest P value that will be considered as important when performing computational optimization in
#' SKAT-O. All P values larger than skato_p_threshold will be processed via burden test. The default is 0.8
#'
#' @param anno_type
#' A character (or character vector) indicating annotation types to be used. The default is "" (i.e, nothing).
#'
#' @param mac_threshold
#' Integer value, specifying the threshold of MACs (Minor allele content) calculated from MAFs. In ACAT, scores with MAC <= 10 will be
#' combined using Burden test.
#'
#' @param regularize_fun
#' Character string, specifying the one of two regularization algorithms if ‘reference_matrix’ is TRUE:  'LH' (default) or 'derivLH'.
#' Currently, both give similar results.
#'
#' @param pca_var_fraction P
#' ositive numeric value, specifying the minimal proportion of genetic variance within the region that should be explained by principal
#' components used in PCA method. This is also valid in 'simpleM'. The default is 0.85.
#'
#' @param flm_basis_function
#' Character string, specifying the name of a basis function type for beta-smooth in FLM method. Can be set to "bspline" (B-spline basis)
#' or "fourier" (Fourier basis, default).
#'
#' @param flm_num_basis
#' Positive integer value, specifying the number of basis functions to be used for beta-smooth in FLM method. The default is 25.
#'
#' @param flm_poly_order
#' Positive integer value, specifying the polynomial order to be used in "bspline" for FLM model. The default = 4, which corresponds to
#' the cubic B-splines. This has no effect if only Fourier bases are used
#'
#' @param flip_genotypes L
#' ogical value, `TRUE` or `FALSE`, indicating whether the genotypes of some genetic variants should be flipped (relabelled) for their
#' better functional representation (13). The default is `FALSE`.
#'
#' @param omit_linear_variant
#' Logical value, `TRUE` or `FALSE`, indicating whether to omit linearly dependent genetic variants. It was done in the FLM test (4).
#' The default is `FALSE`.
#'
#' @param gene_approximation
#' Boolean value, `TRUE` or `FALSE`, specifying whether approximation for large genes (>= 500 SNPs) should be used. Applicable for SKAT,
#' SKATO, sumchi, PCA, FLM. The default is `TRUE` for these methods).
#'
#' @param reference_matrix_used
#' Boolean value, `TRUE` or `FALSE` logical indicating whether the correlation matrices were generated using the reference matrix. The
#' default is `FALSE`. If  `TRUE`, regularization algorithms will be applied to ensure the invertibility and numerical stability of
#' the matrices.
#'
#' @returns
#' A data frame with columns:
#'
#' * gene
#' * chrom
#' * start
#' * end
#' * markers (i.e., numbers of SNPs),
#' * filtered.markers (i.e. filtered SNPs)
#' * pvalue (i.e., p-value).
#'
#' Additionally, for  “BT”, there will be “beta” (i.e., gene-level estimates of betas) and “beta.se” (i.e., standard errors of betas).
#' For “FLM”, there will be the “model” column with the names of the functional models used for each region. Names shortly describe the
#' functional basis and the number of basis functions used. E.g., "F25" means 25 Fourier basis functions, "B15" means 15 B-spline basis
#' functions. For “PCA”, there will be the “ncomponents” (the number of the principal components used for each region) and
#' “explained.variance.fraction” (i.e., the proportion of genetic variance they make up) columns.
#'
#' @references
#' \insertAllCited{}
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom magrittr %>%
#' @importFrom regioneR toGRanges
#' @importFrom plyranges join_overlap_intersect
#' @importFrom sumFREGAT SKAT SKATO sumchi ACAT BT PCA FLM simpleM minp
#' @importFrom rlang inform format_error_bullets
#'
#' @export
#'
#' @examples
#' if (!(Sys.getenv("CI") == "true" && Sys.info()[["sysname"]] == "Darwin")) {
#'     data("XWAS_Summary_Example", package = "GXwasR")
#'     DataDir <- GXwasR:::GXwasR_data()
#'     ResultDir <- tempdir()
#'     finput <- "GXwasR_example"
#'     sumstat <- XWAS_Summary_Example
#'     ref_data <- NULL
#'     gene_file <- "Xlinkedgenes_hg19.txt"
#'     gene_range <- 500000
#'     max_gene <- 10
#'     gene_approximation <- TRUE
#'     beta_par <- c(1, 25)
#'     weights_function <- NULL
#'     geno_variance_weights <- "se.beta"
#'     method <- "kuonen"
#'     acc_devies <- 1e-8
#'     lim_devies <- 1e+6
#'     rho <- TRUE
#'     skato_p_threshold <- 0.8
#'     mac_threshold <- 3
#'     sample_size <- 4000
#'     reference_matrix_used <- FALSE
#'     regularize_fun <- "LH"
#'     pca_var_fraction <- 0.85
#'     flm_basis_function <- "fourier"
#'     flm_num_basis <- 25
#'     flm_poly_order <- 4
#'     flip_genotypes <- FALSE
#'     omit_linear_variant <- FALSE
#'     kernel_p_method <- "kuonen"
#'     anno_type <- ""
#'     GenetestResult <- TestXGene(DataDir, ResultDir, finput, sumstat, gene_file,
#'         gene_range, score_file, ref_data, max_gene, sample_size,
#'         genebasedTest = "SKAT",
#'         gene_approximation, beta_par, weights_function, geno_variance_weights,
#'         kernel_p_method, acc_devies, lim_devies, rho, skato_p_threshold, anno_type,
#'         mac_threshold, reference_matrix_used, regularize_fun, pca_var_fraction,
#'         flm_basis_function, flm_num_basis, flm_poly_order, flip_genotypes,
#'         omit_linear_variant
#'     )
#' }
TestXGene <- function(
      DataDir,
      ResultDir = tempdir(),
      finput,
      sumstat,
      gene_file,
      gene_range = 500000,
      score_file,
      ref_data = NULL,
      max_gene = NULL,
      sample_size = NULL,
      genebasedTest = c(
          "SKAT",
          "SKATO",
          "sumchi",
          "ACAT",
          "BT",
          "PCA",
          "FLM",
          "simpleM",
          "minp"
      ),
      gene_approximation = TRUE,
      beta_par,
      weights_function,
      geno_variance_weights,
      kernel_p_method = "kuonen",
      acc_devies = 1e-8,
      lim_devies = 1e+6,
      rho = TRUE,
      skato_p_threshold = 0.8,
      anno_type = "",
      mac_threshold,
      reference_matrix_used,
      regularize_fun,
      pca_var_fraction = 0.85,
      flm_basis_function = "fourier",
      flm_num_basis = 25,
      flm_poly_order = 4,
      flip_genotypes = FALSE,
      omit_linear_variant = FALSE
) {
    tryCatch(
        withCallingHandlers(
            {
                if (!checkFiles(DataDir, finput)) {
                    stop("Missing required Plink files in the specified DataDir.")
                }

                input.dat <- sumstat[, c("CHROM", "POS", "ID", "A1", "P", "BETA", "EAF")]
                colnames(input.dat) <- c("CHROM", "POS", "ID", "EA", "P", "BETA", "EAF")
                ref.data <- sumstat[, c("CHROM", "POS", "ID", "A2", "A1", "EAF")]
                colnames(ref.data) <-
                    c("CHROM", "POS", "ID", "REF", "ALT", "AF") ## Following convention for reference data

                if (is.null(ref_data)) {
                    ref_data <- ref.data
                } else {
                    ref_data <- ref_data
                }

                geneTestScoreFile(
                    ResultDir = ResultDir, data = input.dat,
                    reference = ref_data,
                    output.file.prefix = "gene.test.score.file"
                ) ## use suppressWarnings()

                genes <- read.table(normalizePath(file.path(DataDir, gene_file), mustWork = FALSE))
                colnames(genes) <- c(c("gene_name", "X", "chr", "Y", "start", "end"))
                genes$up_Mb <- genes$start - gene_range
                genes$down_Mb <- genes$end + gene_range
                genes.gr <- GenomicRanges::makeGRangesFromDataFrame(genes, keep.extra.columns = TRUE)

                suppressWarnings(SNPfile <- read.table(
                    file = paste0(file.path(DataDir, finput), ".bim"),
                    header = FALSE,
                    # na = "NA",
                    na.strings = "NA"
                ))

                SNPfile$chr <- SNPfile$V1
                SNPfile$start <- SNPfile$V4
                SNPfile$end <- SNPfile$V4
                SNPfile$SNP <- SNPfile$V2
                snp_data <- SNPfile %>% select("chr", "start", "end", "SNP")
                snp.gr <- regioneR::toGRanges(snp_data)
                gene_snp_intersect <-
                    as.data.frame(plyranges::join_overlap_intersect(genes.gr, snp.gr))
                rlang::inform(
                    rlang::format_error_bullets(
                        c("i" = paste0(
                            length(
                                unique(
                                    gene_snp_intersect$gene_name
                                )
                            ), " genes are having ", length(
                                unique(
                                    gene_snp_intersect$SNP
                                )
                            ), " SNPs"
                        ))
                    )
                )
                gene_snp <- unique(gene_snp_intersect[, c(6, 11)])
                snpcount <- as.data.frame(table(gene_snp$gene_name))

                g <- as.character(snpcount[, 1])
                dir.create(path = file.path(ResultDir, "cormatrix"))
                rlang::inform(rlang::format_error_bullets("SNP-SNP correlation matrices are being created..."))

                snpcorrFun <- function(g) {
                    snps <- gene_snp[gene_snp$gene_name == g, 2, drop = FALSE]
                    write.table(
                        snps,
                        file = normalizePath(file.path(ResultDir, "cor_snps.txt"), mustWork = FALSE),
                        quote = FALSE,
                        row.names = FALSE,
                        col.names = FALSE,
                        eol = "\r\n"
                    )

                    invisible(sys::exec_wait(
                        plink(),
                        args = c(
                            "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE),
                            "--r2", "square",
                            "--extract", normalizePath(file.path(ResultDir, "cor_snps.txt"), mustWork = FALSE),
                            "--out", normalizePath(file.path(ResultDir, "snpcorr"), mustWork = FALSE),
                            "--silent"
                        ),
                        std_out = FALSE,
                        std_err = FALSE
                    ))

                    snpcorr <- as.matrix(read.table(file = normalizePath(file.path(ResultDir, "snpcorr.ld"), mustWork = FALSE)))
                    colnames(snpcorr) <- snps$SNP
                    rownames(snpcorr) <- snps$SNP
                    snpcorr <- as.data.frame(snpcorr)
                    save(snpcorr, file = normalizePath(file.path(ResultDir, "cormatrix", paste0(g, ".RData")), mustWork = FALSE))
                    return()
                }

                invisible(lapply(g, snpcorrFun))
                rlang::inform(rlang::format_error_bullets(c("v" = "SNP-SNP correlation matrices are done.")))

                score_file <- normalizePath(file.path(ResultDir, "gene.test.score.file.vcf.gz"), mustWork = FALSE)
                gene.file <- gene_file

                if (is.null(max_gene)) {
                    genes1 <- as.vector(g)
                } else {
                    maxgene <- max_gene
                    genes1 <- as.vector(g[seq_len(max_gene)])
                }

                if (genebasedTest == "SKAT") {
                    return(runSKAT(
                        score.file = score_file,
                        gene.file = normalizePath(file.path(DataDir, gene_file), mustWork = FALSE),
                        genes = genes1,
                        cor.path = normalizePath(file.path(ResultDir, "cormatrix"), mustWork = FALSE),
                        gene_approximation = gene_approximation,
                        anno.type = anno_type,
                        beta.par = beta_par,
                        weights.function = weights_function,
                        geno_variance_weights = geno_variance_weights,
                        kernel_p_method = kernel_p_method,
                        acc_devies = acc_devies,
                        lim_devies = lim_devies,
                        rho = rho,
                        skato_p_threshold = skato_p_threshold
                    ))
                } else if (genebasedTest == "SKATO") {
                    return(runSKATO(
                        score_file,
                        normalizePath(file.path(DataDir, gene_file), mustWork = FALSE),
                        genes1,
                        normalizePath(file.path(ResultDir, "cormatrix"), mustWork = FALSE),
                        anno_type, gene_approximation, beta_par, weights_function, kernel_p_method, acc_devies, lim_devies, rho, skato_p_threshold
                    ))
                } else if (genebasedTest == "sumchi") {
                    return(runSumChi(
                        score_file,
                        normalizePath(file.path(DataDir, gene_file), mustWork = FALSE),
                        genes1,
                        normalizePath(file.path(ResultDir, "cormatrix"), mustWork = FALSE),
                        gene_approximation, anno_type, kernel_p_method, acc_devies, lim_devies
                    ))
                } else if (genebasedTest == "ACAT") {
                    return(runACAT(
                        score_file,
                        normalizePath(file.path(DataDir, gene_file), mustWork = FALSE),
                        genes1, anno_type, beta_par, weights_function, geno_variance_weights, mac_threshold, sample_size
                    ))
                } else if (genebasedTest == "BT") {
                    return(runBT(
                        score_file,
                        normalizePath(file.path(DataDir, gene_file), mustWork = FALSE),
                        genes1, normalizePath(file.path(ResultDir, "cormatrix"), mustWork = FALSE),
                        anno_type, beta_par, weights_function
                    ))
                } else if (genebasedTest == "PCA") {
                    return(runPCA(
                        score_file,
                        normalizePath(file.path(DataDir, gene_file), mustWork = FALSE),
                        genes1,
                        normalizePath(file.path(ResultDir, "cormatrix"), mustWork = FALSE),
                        gene_approximation, anno_type, sample_size, beta_par, weights_function, reference_matrix_used, regularize_fun, pca_var_fraction
                    ))
                } else if (genebasedTest == "FLM") {
                    return(runFLM(
                        score_file,
                        normalizePath(file.path(DataDir, gene_file), mustWork = FALSE),
                        genes1, normalizePath(file.path(ResultDir, "cormatrix"), mustWork = FALSE),
                        gene_approximation, anno_type, sample_size, beta_par, weights_function, flm_basis_function, flm_num_basis, flm_poly_order,
                        flip_genotypes, omit_linear_variant, reference_matrix_used, regularize_fun
                    ))
                } else if (genebasedTest == "simpleM") {
                    return(runSimpleM(
                        score_file, gene_file, genes1, anno_type, pca_var_fraction
                    ))
                } else if (genebasedTest == "minp") {
                    return(runMinP(
                        score_file, gene_file, genes1, anno_type
                    ))
                }
            },
            error = function(e) {
                message("An error occurred: ", e$message)
                return(NULL)
            },
            warning = function(w) {
                rlang::inform(
                    rlang::format_error_bullets(
                        c("!" = paste("Warning:", conditionMessage(w)))
                    )
                )
                invokeRestart("muffleWarning")
            }
        )
    )
}

#' MetaGWAS: Combining summary-level results from two or more GWA studies into a single estimate.
#'
#' @description
#' This function combine K sets of GWAS association statistics on same (or at least similar) phenotype. This function employs
#' PLINK's \insertCite{Purcell2007}{GXwasR} inverse variance-based analysis to run a number of models, including a)
#' Fixed-effect model and b) Random-effect model, assuming there may be variation between the genuine underlying effects,
#' i.e., effect size beta. 'This function also calculates weighted Z-score-based p-values after METAL \insertCite{Willer2010}{GXwasR}.
#' For more information about the algorithms, please see the associated paper.
#'
#' @param DataDir
#' A character string for the file path of the input files needed for `SummData` and `SNPfile` arguments.
#'
#' @param SummData
#' Vector value containing the name(s) of the .Rda file(s) with GWAS summary statistics, with ‘SNP’
#' (i.e., SNP identifier), ‘BETA’ (i.e., effect-size or logarithm of odds ratio), ‘SE’ (i.e., standard error of BETA),
#' ‘P’ (i.e., p-values), 'NMISS' (i.e., effective sample size), 'L95' (i.e., lower limit of 95% confidence interval) and
#' 'U95' (i.e., upper limit of 95% confidence interval) are in mandatory column headers. These files needed to be in DataDir.
#' If the numbers of cases and controls are unequal, effective sample size should be \eqn{4 / (1/<qty of cases> + 1/<qty of controls>)}.
#' A smaller "effective" sample size may be used for samples that include related individuals, however simulations indicate
#' that small changes in the effective sample size have relatively little effect on the final p-value
#' \insertCite{Willer2010}{GXwasR}. Columns, such as, `CHR` (Chromosome code), `BP` (Basepair position), `A1` (First allele code),
#' `A2` (Second allele code) columns are optional. If these are present, setting `useSNPposition` to `FALSE`, causes `CHR`, `BP`
#' and `A1` to be ignored and setting `UseA1` to be `FALSE` causes `A1` to be ignored. If, both these arguments are `TRUE`, this
#' function takes care of A1/A2 allele flips properly. Otherwise, A1 mismatches are thrown out. Values of CHR/BP are allowed
#' to vary.
#'
#' @param ResultDir
#' A character string for the file path where all output files will be stored. The default is `tempdir()`.
#'
#' @param SNPfile
#' Character string specifying the name of the plain-text file with a column of SNP names. These could be LD clumped SNPs or
#' any other list of chosen SNPs for Meta analysis. This file needs to be in `DataDir`.
#'
#' @param useSNPposition
#' Boolean value, `TRUE` or `FALSE` for using `CHR`, `BP`, and `A1` or not. The default is `FALSE.` Note: if
#' this is `FALSE` then there will be no Manhattan and QQ plot will be generated.
#'
#' @param UseA1
#' Boolean value, `TRUE` or `FALSE` for `A1` to be used or not. The default is `FALSE`.
#'
#' @param GCse
#' Boolean value, `TRUE` or `FALSE` for applying study specific genomic control to adjust each study for potential population
#' structure for all the SNPs. The default is `TRUE`. If users would want to apply genomic control separately for directly
#' genotyped and imputed SNPs prior using the function, set this parameter as `FALSE`.
#'
#' @param plotname
#' Character string, specifying the plot name of the file containing forest plots for the SNPs. The default is
#' “Meta_Analysis.plot”.
#'
#' @param pval_filter
#' Character value as "R","F" or "W", specifying whether p-value threshold should be chosen based on “Random”, “Fixed” or
#' “Weighted” effect model for the SNPs to be included in the forest plots.
#'
#' @param top_snp_pval
#' Numeric value, specifying the threshold to be used to filter the SNPs for the forest plots. The default is 1e-08.
#'
#' @param max_top_snps
#' Integer value, specifying the maximum number of top SNPs (SNPs with the lowest p-values) to be ploted in the forest
#' plot file. The default is 6.
#'
#' @param chosen_snps_file
#' Character string specifying the name of the plain-text file with a column of SNP names for the forest plots.
#' The default is NULL.
#'
#' @param byCHR
#' Boolean value, `TRUE` or `FALSE`, specifying whether the meta analysis will be performed chromosome wise or not.
#' The default is `FALSE`.
#'
#' @param pval_threshold_manplot
#' Numeric value, specifying the p-value threshold for plotting Manhattan plots.
#'
#' @returns
#' A list object containing five dataframes and a list of forest plots. The first three dataframes, such as Mfixed, Mrandom and Mweighted contain results
#' for fixed effect, random effect and weighted model. Each of these dataframes can have maximum 12 columns, such as:
#' * `CHR` (Chromosome code)
#' * `BP` (Basepair position)
#' * `SNP` (SNP identifier)
#' * `A1` (First allele code)
#' * `A2` (Second allele code)
#' * `Q` (p-value for Cochrane's Q statistic)
#' * `I` (I^2 heterogeneity index (0-100))
#' * `P` (P-value from mata analysis)
#' * `ES` (Effect-size estimate from mata analysis)
#' * `SE` (Standard Error from mata analysis)
#' * `CI_L` (Lower limit of confidence interval)
#' * `CI_U` (Uper limit of confidence interval)
#'
#' The fourth dataframe contains the same columns `CHR`, `BP`, `SNP`, `A1`, `A2`, `Q`, `I`", with column `N`' ( Number of
#' valid studies for this SNP), P (Fixed-effects meta-analysis p-value), and other columns as `Fx...` (Study x (0-based input file
#' indices) effect estimate, Examples: F0, F1 etc.).
#'
#' The fifth dataframe, ProblemSNP has three columns, such as
#' * `File` (file name of input data),
#' * `SNP` (Problematic SNPs that are thrown)
#' * `Problem` (Problem code)
#'
#' Problem codes are:
#' * BAD_CHR (Invalid chromosome code)
#' * BAD_BP  (Invalid base-position code), BAD_ES (Invalid effect-size)
#' * BAD_SE (Invalid standard error), MISSING_A1 (Missing allele 1 label)
#' * MISSING_A2 (Missing allele 2 label)
#' * ALLELE_MISMATCH (Mismatching allele codes across files)
#'
#' A .pdf file comprising the forest plots of the SNPs is produced in the ResultDir with Plotname as prefix.
#' If `useSNPposition` is set `TRUE`, a .jpeg file with Manhattan Plot and Q-Q plot will be in the `ResultDir` with Plotname
#' as prefix.
#'
#' @references
#' \insertAllCited{}
#'
#' @importFrom qqman manhattan qq
#' @importFrom graphics par
#' @importFrom grDevices jpeg pdf
#'
#' @export
#'
#' @examples
#' data("Summary_Stat_Ex1", package = "GXwasR")
#' data("Summary_Stat_Ex2", package = "GXwasR")
#' DataDir <- GXwasR:::GXwasR_data()
#' ResultDir <- tempdir()
#' SummData <- list(Summary_Stat_Ex1, Summary_Stat_Ex2)
#' SNPfile <- "UniqueLoci"
#' useSNPposition <- FALSE
#' UseA1 <- TRUE
#' GCse <- TRUE
#' byCHR <- FALSE
#' pval_filter <- "R"
#' top_snp_pval <- 1e-08
#' max_top_snps <- 10
#' chosen_snps_file <- NULL
#' pval_threshold_manplot <- 1e-05
#' plotname <- "Meta_Analysis.plot"
#' x <- MetaGWAS(
#'     DataDir = DataDir, SummData = SummData, ResultDir = ResultDir,
#'     SNPfile = NULL, useSNPposition = TRUE, UseA1 = UseA1, GCse = GCse,
#'     plotname = "Meta_Analysis.plot", pval_filter, top_snp_pval, max_top_snps,
#'     chosen_snps_file = NULL, byCHR, pval_threshold_manplot
#' )
MetaGWAS <- function(DataDir, SummData = c(""), ResultDir = tempdir(), SNPfile = NULL,
    useSNPposition = TRUE,
    UseA1 = FALSE, GCse = TRUE,
    plotname = "Meta_Analysis.plot", pval_filter = "R",
    top_snp_pval = 1e-08, max_top_snps = 6, chosen_snps_file = NULL,
    byCHR = FALSE, pval_threshold_manplot = 1e-05) {
    # Validate input parameters
    validateInputForMetaGWAS(DataDir, ResultDir, SummData, SNPfile, useSNPposition, UseA1, GCse, plotname, pval_filter, top_snp_pval, max_top_snps, chosen_snps_file, byCHR, pval_threshold_manplot)

    tryCatch(
        {
            # Setup nominal parameters
            if (useSNPposition == TRUE) {
                nomap <- NULL
            } else {
                nomap <- "no-map"
            }

            if (UseA1 == TRUE) {
                UseA1v <- NULL
            } else {
                UseA1v <- "no-allele"
            }

            if (is.null(SNPfile) | byCHR == TRUE) {
                extract <- NULL
                SNPfilev <- NULL
            } else {
                extract <- "--extract"
                SNPfilev <- SNPfile
            }
            # Write summary data files to ResultDir
            for (i in seq_along(SummData)) {
                rlang::inform(rlang::format_error_bullets(c("i" = paste0("Processing file number ", i))))
                write.table(SummData[[i]], normalizePath(file.path(ResultDir, paste0("SNPdata_", i)), mustWork = FALSE), row.names = FALSE, col.names = TRUE, quote = FALSE)
            }

            SummData <- list.files(normalizePath(file.path(ResultDir), mustWork = FALSE), pattern = "SNPdata_")

            # Apply genomic control
            if (GCse == TRUE) {
                invisible(lapply(SummData, getGCse, ResultDir = ResultDir))
            } else {
                rlang::inform(rlang::format_error_bullets(c("i" = "No study-specific genomic control was applied.")))
            }

            if (is.null(SNPfile)) {
                SNPfilev <- NULL
            } else {
                SNPfilev <- normalizePath(file.path(DataDir, SNPfile), mustWork = FALSE)
            }

            if (byCHR == FALSE) {
                MR <- metaFun(
                    DataDir = DataDir, ResultDir = ResultDir, SummData = SummData,
                    CHR = NULL, chromosome = NULL, nomap = nomap,
                    UseA1v = UseA1v, extract = extract, SNPfilev = SNPfilev
                )
            } else {
                chromosome <- seq_len(23)
                MR <- data.table::rbindlist(lapply(chromosome, metaFun,
                    DataDir = DataDir, ResultDir = ResultDir, SummData = SummData,
                    CHR = "--chr", nomap = NULL,
                    UseA1v = NULL, extract = extract, SNPfilev = SNPfilev
                ))
            }


            # Calculate standard error and 95% confidence intervals

            # Standard error
            MR$SEfixed <- abs(MR$BETA / qnorm(MR$P / 2))
            MR$SErandom <- abs(MR$BETA.R. / qnorm(MR$P.R. / 2))
            MR$SEweighted <- abs(MR$WEIGHTED_Z / qnorm(MR$P.WZ. / 2))

            # 95% Confidence interval
            MR$CIfixedLL <- MR$BETA - 1.96 * MR$SEfixed
            MR$CIfixedUL <- MR$BETA + 1.96 * MR$SEfixed

            MR$CIrandomLL <- MR$BETA.R. - 1.96 * MR$SErandom
            MR$CIrandomUL <- MR$BETA.R. + 1.96 * MR$SErandom

            MR$CIweightedLL <- MR$WEIGHTED_Z - 1.96 * MR$SEweighted
            MR$CIweightedUL <- MR$WEIGHTED_Z + 1.96 * MR$SEweighted

            # Get effect size and confidence interval for studies - only for SNPs in MR
            MRsnps <- unique(MR[, "SNP", drop = FALSE])
            Sbeta <- data.table::rbindlist(lapply(SummData, getStudyCI, MRsnps = MRsnps, ResultDir = ResultDir))


            # Filter and prepare SNPs for forest plots
            top_snp_pval <- adjustPvalThreshold(top_snp_pval, MR, pval_filter)
            MRfiltered <- filterSNPsForForestPlot(MR, top_snp_pval, pval_filter)


            # Update MRfiltered if a specific SNP file is provided
            if (!is.null(chosen_snps_file)) {
                chosenS <- read.table(normalizePath(file.path(DataDir, chosen_snps_file), mustWork = FALSE), header = TRUE)
                colnames(chosenS) <- "SNP"
                MRfiltered <- merge(chosenS, MR, by = "SNP")
            }

            # Visualization
            forest_plots <- generatePlots(MRfiltered, Sbeta, ResultDir, plotname, useSNPposition, pval_threshold_manplot, chosen_snps_file)


            ## Produce all forest plots in .pdf
            grDevices::pdf(normalizePath(file.path(ResultDir, paste0(plotname, ".pdf")), mustWork = FALSE), width = 10, height = 5)
            i <- seq_len(length(MRfiltered$SNP))
            invisible(suppressWarnings(lapply(i, allForestplot, MR2 = MRfiltered, Sbeta = Sbeta)))
            dev.off()

            rlang::inform(rlang::format_error_bullets(c(
                "v" = paste0("Forest plots for SNPS have compiled and saved as ", plotname, ".pdf"),
                "i" = paste("You can find them in the directory:", ResultDir)
            )))

            # Check for problem SNPs
            problemFile <- normalizePath(file.path(ResultDir, "MetaResult.prob"), mustWork = FALSE)
            if (file.exists(problemFile)) {
                MP <- read.table(problemFile)
                colnames(MP) <- c("File", "SNP", "Problem")
            } else {
                MP <- data.frame(File = "None", SNP = "None", Problem = "None")
            }

            # Cleanup temporary files
            removeTempFiles(ResultDir, "SNPdata")
            removeTempFiles(ResultDir, "MetaResult")

            if (useSNPposition == TRUE) {
                Mfixed <- MR[, c("CHR", "BP", "SNP", "A1", "A2", "Q", "I", "P", "BETA", "SEfixed", "CIfixedLL", "CIfixedUL")]
                colnames(Mfixed) <- c("CHR", "BP", "SNP", "A1", "A2", "Q", "I", "P", "ES", "SE", "CI_L", "CI_U")
                Mrandom <- MR[, c("CHR", "BP", "SNP", "A1", "A2", "Q", "I", "P.R.", "BETA.R.", "SErandom", "CIrandomLL", "CIrandomUL")]
                colnames(Mrandom) <- c("CHR", "BP", "SNP", "A1", "A2", "Q", "I", "P", "ES", "SE", "CI_L", "CI_U")
                Mweighted <- MR[, c("CHR", "BP", "SNP", "A1", "A2", "Q", "I", "P.WZ.", "WEIGHTED_Z", "SEweighted", "CIweightedLL", "CIweightedUL")]
                colnames(Mweighted) <- c("CHR", "BP", "SNP", "A1", "A2", "Q", "I", "P", "ES", "SE", "CI_L", "CI_U")
                Msummdata <- MR[, !names(MR) %in%
                    c("P", "BETA", "SEfixed", "CIfixedLL", "CIfixedUL", "P.R.", "BETA.R.", "SErandom", "CIrandomLL", "CIrandomUL", "P.WZ.", "WEIGHTED_Z", "SEweighted", "CIweightedLL", "CIweightedUL")]

                ## Plot Mahattan and QQ plots
                options(bitmapType = "cairo")
                grDevices::jpeg(normalizePath(file.path(ResultDir, paste0(plotname, ".jpeg")), mustWork = FALSE),
                    width = 20,
                    height = 10,
                    units = "in",
                    res = 300
                )
                graphics::par(mfrow = c(3, 2))
                mR <- na.omit(Mfixed[, c("SNP", "CHR", "BP", "P")])
                # From p-values, calculate chi-squared statistic
                chisq <- qchisq(1 - na.omit(mR$P), 1)
                lamdaGC <- median(chisq) / qchisq(0.5, 1)

                invisible(suppressWarnings(qqman::manhattan(mR, ylim = c(0, 10), annotatePval = pval_threshold_manplot, annotateTop = FALSE, main = "Manhattan plot of fixed effect meta GWAS")))
                invisible(suppressWarnings(qqman::qq(mR$P, main = paste0(("Q-Q plot of fixed effect meta GWAS p-values with GIF = "), lamdaGC))))

                mR <- na.omit(Mrandom[, c("SNP", "CHR", "BP", "P")])
                # From p-values, calculate chi-squared statistic
                chisq <- qchisq(1 - na.omit(mR$P), 1)
                lamdaGC <- median(chisq) / qchisq(0.5, 1)
                invisible(suppressWarnings(qqman::manhattan(mR, ylim = c(0, 10), annotatePval = pval_threshold_manplot, annotateTop = FALSE, main = "Manhattan plot of random effect meta GWAS")))
                invisible(suppressWarnings(qqman::qq(mR$P, main = paste0(("Q-Q plot of random effect meta GWAS p-values with GIF = "), lamdaGC))))

                mR <- na.omit(Mweighted[, c("SNP", "CHR", "BP", "P")])
                # From p-values, calculate chi-squared statistic
                chisq <- qchisq(1 - na.omit(mR$P), 1)
                lamdaGC <- median(chisq) / qchisq(0.5, 1)
                invisible(suppressWarnings(qqman::manhattan(mR, ylim = c(0, 10), annotatePval = pval_threshold_manplot, annotateTop = TRUE, main = "Manhattan plot of weighted Z-score meta GWAS")))
                invisible(suppressWarnings(qqman::qq(mR$P, main = paste0(("Q-Q plot of weighted Z-score meta GWAS p-values with GIF = "), lamdaGC))))
                dev.off()
                rlang::inform(rlang::format_error_bullets(c(
                    "v" = paste0("Manhattan and QQ plots for SNPS have been compiled and saved as ", plotname, ".jpeg"),
                    "i" = paste("You can find them in the directory:", ResultDir)
                )))

                #####
            } else {
                Mfixed <- MR[, c("SNP", "Q", "I", "P", "BETA", "SEfixed", "CIfixedLL", "CIfixedUL")]
                colnames(Mfixed) <- c("SNP", "Q", "I", "P", "ES", "SE", "CI_L", "CI_U")
                Mrandom <- MR[, c("SNP", "Q", "I", "P.R.", "BETA.R.", "SErandom", "CIrandomLL", "CIrandomUL")]
                colnames(Mrandom) <- c("SNP", "Q", "I", "P", "ES", "SE", "CI_L", "CI_U")
                Mweighted <- MR[, c("SNP", "Q", "I", "P.WZ.", "WEIGHTED_Z", "SEweighted", "CIweightedLL", "CIweightedUL")]
                colnames(Mweighted) <- c("SNP", "Q", "I", "P", "ES", "SE", "CI_L", "CI_U")
                Msummdata <- MR[, !names(MR) %in%
                    c("P", "BETA", "SEfixed", "CIfixedLL", "CIfixedUL", "P.R.", "BETA.R.", "SErandom", "CIrandomLL", "CIrandomUL", "P.WZ.", "WEIGHTED_Z", "SEweighted", "CIweightedLL", "CIweightedUL")]

                rlang::inform(rlang::format_error_bullets(c("i" = "Since useSNPposition = FALSE, there will be no Manhattan and QQ plot will be generated.")))
            }
            return(list(Resultfixed = Mfixed, Resultrandom = Mrandom, Resultweighted = Mweighted, Metadata = Msummdata, ProblemSNP = MP, forestPlots = forest_plots))
        },
        error = function(e) {
            message("An error occurred: ", e$message)
            return(NULL)
        },
        warning = function(w) {
            message("Warning: ", w$message)
        }
    )
}

#' ComputePGS: Computing polygenic score (PGS)
#'
#' @author Banabithi Bose
#'
#' @description This function calculates the polygenic risk score, which summarizes the estimated effect of many genetic variants on an
#' individual’s phenotype. It is calculated  as the sum of the allele counts (genotypes), each weighted by their estimated phenotypic
#' effect sizes from genome-wide association studies. It uses C+T filtering techniques. Users can perform the clumping procedure
#' chromosome-wise or genome-wide. Also, the function offers the choice of including genetic principal components and other covariates.
#' Using this function, users have freedom to experiment with various clumping and thresholding arrangements to test a wide range of
#' various parameter values.
#'
#' @param DataDir
#' A character string specifying file path of the all the input files.
#'
#' @param ResultDir
#' A character string specifying file path where all output files will be stored. The default is tempdir().
#'
#' @param finput
#' Character string, specifying the prefix of the input PLINK binary files containing the genotype data i.e., the target data based on
#' which the clumping procedure will be performed. This file needs to be in DataDir. If the sample size of the target dataset is small
#' (e.g., N < 500 individuals) then users can utilize the 1000 Genomes Project samples, ensuring the use of the population that most
#' closely represents the target sample.
#'
#' @param summarystat
#' A dataframe object containing GWAS summary statistics.
#'
#' The mandatory column headers in this dataframe are:
#' * `CHR`(Chromosome code)
#' * `BP`(Basepair position)
#' * `A1` (effect allele)
#' * `SNP` (i.e., SNP identifier)
#' * `BETA` or `OR` (i.e., effect-size or logarithm of odds ratio)
#' * `P` (i.e., p-values).
#'
#' Special Notes: The first three columns needed to be `SNP`, `A1` and `BETA` or `OR`.
#'
#' @param phenofile
#' A character string, specifying the name of the mandatory phenotype file. This is a plain text file with no header line; columns are:
#' family ID, individual ID, and phenotype.  For a binary trait, the phenotypic value should be coded as 0 or 1, then it will be recognized
#' as a case-control study (0 for controls and 1 for cases). The missing value should be represented by "-9" or "NA". The desired phenotype
#' column should be labeled as "Pheno1". This file needs to be in 'DataDir'.
#'
#' @param covarfile
#' A character string, specifying the name of the covariate file which is a plain .text file with no header line; columns are family ID,
#' individual ID, and the covariates. The default is 'NULL'. This file needs to be in 'DataDir'.
#'
#' @param pheno_type
#' Boolean value, ‘binary’ or ‘quantitative’, specifying the type of the trait. The default is ‘binary’.
#'
#' @param effectsize
#' Boolean value, 'BETA' or 'OR', specifying the type of the GWAS effect size. The default is 'BETA'.
#'
#' @param ldclump
#' Boolean value, `TRUE` or `FALSE`, specifying whether to perform clumping or not.
#'
#' @param LDreference
#' A character string, specifying the prefix of the PLINK files of the genetic similarity reference panel, (ideally the same that was used
#' to impute the target dataset). These files should be in 'DataDir'.
#'
#' @param clump_p1
#' Numeric value, specifying the significance threshold for index SNPs if `ldclump` was set to be `TRUE`. The default is 0.0001.
#'
#' @param clump_p2
#' Numeric value, specifying the secondary significance threshold for clumped SNPs if `ldclump` was set to be `TRUE`. The default is 0.01
#'
#' @param clump_r2
#' Numeric value, specifying the linkage disequilibrium (LD) threshold for clumping if `ldclump` was set to be `TRUE`. The default is 0.50.
#'
#' @param clump_kb
#' Integer value, specifying the physical distance threshold in base-pair for clumping if `ldclump` was set to be `TRUE`. The default is 250.
#'
#' @param byCHR
#' Boolean value, 'TRUE' or 'FALSE', specifying chromosome-wise clumping if 'ldclump' was set to be 'TRUE'. The default is 'TRUE'.
#'
#' @param pthreshold
#' Numeric vector, containing several p value thresholds to maximize predictive ability of the derived polygenic scores.
#'
#' @param ld_prunning
#' Boolean value, `TRUE` or `FALSE` for LD-based filtering for computing genetic PC as covariates.
#'
#' @param nPC
#' Positive integer value, specifying the number of genetic PCs to be included as predictor in the PGS model fit. The default is 6.
#'
#' @param window_size
#' Integer value, specifying a window size in variant count or kilobase for LD-based filtering in computing genetic PC. The default is 50.
#'
#' @param step_size
#' Integer value, specifying a variant count to shift the window at the end of each step for LD filtering in computing genetic PCs. The default is 5.
#'
#' @param r2_threshold
#' Numeric value between 0 to 1 of pairwise \eqn{r^2} threshold for LD-based filtering in computing genetic PCs. The default is 0.02.
#'
#' @param highLD_regions
#' Character string, specifying the .txt file name with known genomic regions with high LD. The default is `NULL`.
#'
#' @return
#' A list object containing a dataframe, a numeric value, a GeneticPC plot (if requested), and a PGS plot. The dataframe, PGS, contains four
#' mandatory columns: IID (i.e., Individual ID), FID (i.e., Family ID), Pheno1 (i.e., the trait for PGS) and Score (i.e., the best PGS). Other
#' columns of covariates could be there. The numeric value, BestP contains the threshold of the best p-value for the best PGS model fit. Also,
#' the function produces several plots, including p-value thresholds vs PGS model fit and PGS distribution among male and females. For
#' case-control data, it also plots the PGS distribution among cases and controls, and plots ROC curves.
#'
#' Also, the function produces several plots such as p-value thresholds vs PGS model fit and PGS distribution among male and females.
#' For case-control data, it shows PGS distribution among cases and controls and ROC curves as well.
#'
#' @importFrom dplyr distinct
#' @importFrom stats lm predict logLik
#' @importFrom ggplot2 theme_classic theme element_text ggtitle geom_density xlab aes scale_y_continuous geom_bar scale_fill_gradient2
#' @importFrom grid grid.newpage grid.draw
#' @importFrom data.table as.data.table
#' @importFrom ggpubr ggarrange
#' @export
#'
#' @examples
#' data("Summary_Stat_Ex1", package = "GXwasR")
#' data("Example_phenofile", package = "GXwasR")
#' data("Example_covarfile", package = "GXwasR")
#' data("Example_pthresoldfile", package = "GXwasR")
#' data("highLD_hg19", package = "GXwasR")
#' DataDir <- GXwasR:::GXwasR_data()
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' summarystat <- Summary_Stat_Ex1[, c(2, 4, 7, 1, 3, 12)]
#' phenofile <- Example_phenofile # Cannot be NULL
#' # The interested phenotype column should be labeled as "Pheno1".
#' covarfile <- Example_covarfile
#' clump_p1 <- 0.0001
#' clump_p2 <- 0.0001
#' clump_kb <- 500
#' clump_r2 <- 0.5
#' byCHR <- TRUE
#' pthreshold <- Example_pthresoldfile$Threshold
#' ld_prunning <- TRUE
#' highLD_regions <- highLD_hg19
#' window_size <- 50
#' step_size <- 5
#' r2_threshold <- 0.02
#' nPC <- 6 # We can incorporate PCs into our PGS analysis to account for population stratification.
#' pheno_type <- "binary"
#'
#' PGSresult <- ComputePGS(DataDir, ResultDir, finput, summarystat, phenofile, covarfile,
#'     effectsize = "BETA", LDreference = "GXwasR_example", ldclump = FALSE, clump_p1, clump_p2,
#'     clump_r2, clump_kb, byCHR = TRUE, pthreshold = pthreshold, highLD_regions = highLD_regions,
#'     ld_prunning = TRUE, window_size = 50, step_size = 5, r2_threshold = 0.02, nPC = 6,
#'     pheno_type = "binary"
#' )
#'
#' ## This table shows 10 samples with phenotype, covariates and a PGS column.
#' PGS <- PGSresult$PGS
#' PGS[seq_len(10), ]
#' ## The best threshold
#' BestPvalue <- PGSresult$BestP$Threshold
#' BestPvalue
ComputePGS <- function(
      DataDir, ResultDir = tempdir(), finput, summarystat, phenofile, covarfile = NULL,
      effectsize = c("BETA", "OR"), ldclump = FALSE, LDreference, clump_p1, clump_p2, clump_r2, clump_kb, byCHR = TRUE,
      pthreshold = c(0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5), highLD_regions, ld_prunning = FALSE,
      window_size = 50, step_size = 5, r2_threshold = 0.02, nPC = 6, pheno_type = "binary"
) {
    # Validate inputs
    if (!validateInputForComputePGS(DataDir, ResultDir, finput, summarystat, phenofile, covarfile, effectsize, ldclump, LDreference, clump_p1, clump_p2, clump_r2, clump_kb, byCHR, pthreshold, highLD_regions, ld_prunning, window_size, step_size, r2_threshold, nPC, pheno_type)) {
        stop("Please validate all inputs")
    }

    if (!checkFiles(DataDir, finput)) {
        stop("Missing required Plink files in the specified DataDir.")
    }

    tryCatch(
        {
            if (effectsize == "OR") {
                summarystat$OR <- log(summarystat$OR)
            }

            summarystat <- as.data.frame(dplyr::distinct(summarystat, summarystat$SNP, .keep_all = TRUE))

            write.table(summarystat, file = normalizePath(file.path(ResultDir, "pgssummarystat"), mustWork = FALSE), quote = FALSE, row.names = FALSE)
            SNP.pvalue <- unique(summarystat[, c("SNP", "P")])
            write.table(SNP.pvalue, file = normalizePath(file.path(ResultDir, "SNP.pvalue"), mustWork = FALSE), quote = FALSE, row.names = FALSE)


            # LD Clumping
            clumpResults <- performLDClumping(ldclump, DataDir, ResultDir, LDreference, summarystat, clump_p1, clump_p2, clump_r2, clump_kb, byCHR)
            clumpExtract <- clumpResults$clumpExtract
            clumpSNP <- clumpResults$clumpSNP

            # Prepare Phenotype Data
            GP <- preparePhenotypeData(phenofile, nPC, DataDir, ResultDir, finput, highLD_regions, ld_prunning, window_size, step_size, r2_threshold)

            # Read in the phenotype file
            phenotype <- cbind(phenofile[, seq_len(2)], phenofile[, "Pheno1"])
            colnames(phenotype) <- c("FID", "IID", "Pheno1")

            # Read the covariates (here, it is sex)
            if (is.null(covarfile)) {
                pheno <- merge(phenotype, GP$PCs1, by = c("FID", "IID"))
            } else {
                covariate <- covarfile
                # Now merge the files
                pheno <- merge(merge(phenotype, covariate, by = c("FID", "IID")), GP$PCs1, by = c("FID", "IID"))
            }

            # We can then calculate the null model (model with PGS) using a linear regression
            # (as height is quantitative) ## Check for binary
            # Compute Null Model
            null_model_result <- computeNullModel(pheno, pheno_type)

            # Accessing the null model and R-squared value
            null_model <- null_model_result$model
            # null_r2 <- null_model_result$r_squared

            ## PGS using thresholding
            ## Best fit PGS
            pgsResult <- data.table::rbindlist(lapply(pthreshold, function(pt) {
                pgsFun(pt, ResultDir, DataDir, finput, clumpExtract, clumpSNP, pheno, pheno_type, null_model)
            }))


            # Generate a pretty format for p-value output
            pgsResult$WriteP <- round(pgsResult$P, digits = 3)
            pgsResult$WriteP[!is.na(pgsResult$WriteP) & pgsResult$WriteP == 0] <- format(pgsResult$P[!is.na(pgsResult$WriteP) & pgsResult$WriteP == 0], digits = 2)
            pgsResult$WriteP <- sub("e", "*x*10^", pgsResult$WriteP)

            p1 <- createPGSPlot(pgsResult)

            # Best result is:
            bestP <- pgsResult[which.max(pgsResult$R2), "Threshold"]
            # Getting PGS score with best p-value threshold
            pt <- cbind(bestP, 0, bestP)
            write.table(pt, file = normalizePath(file.path(ResultDir, "range_list"), mustWork = FALSE), quote = FALSE, row.names = FALSE)
            # By default, if a genotype in the score is missing for a particular individual, then the expected value is imputed, i.e. based on the sample allele frequency. To change this behavior, add the flag --score-no-mean-imputation
            invisible(sys::exec_wait(
                plink(),
                args = c(
                    "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE),
                    "--score", normalizePath(file.path(ResultDir, "pgssummarystat"), mustWork = FALSE), 1, 2, 3, "header",
                    "--q-score-range", normalizePath(file.path(ResultDir, "range_list"), mustWork = FALSE), normalizePath(file.path(ResultDir, "SNP.pvalue"), mustWork = FALSE),
                    clumpExtract, clumpSNP,
                    "--out", normalizePath(file.path(ResultDir, "PGS"), mustWork = FALSE),
                    "--silent"
                ),
                std_out = FALSE,
                std_err = FALSE
            ))

            pgs <- read.table(normalizePath(file.path(ResultDir, paste0("PGS.", bestP, ".profile")), mustWork = FALSE), header = TRUE)
            pheno.pgs <- merge(pheno, pgs[, c("FID", "IID", "SCORE")], by = c("FID", "IID"))

            ## PGS with sex
            d1 <- pheno.pgs[, c("FID", "IID", "Pheno1"), drop = FALSE]
            famfile <- read.table(normalizePath(file.path(DataDir, paste0(finput, ".fam")), mustWork = FALSE), header = FALSE)
            sex <- famfile[!famfile$V5 == 0, c(1, 2, 5)]
            colnames(sex) <- c("FID", "IID", "SEX")
            dat <- merge(d1, sex, by = c("FID", "IID"))

            # Rename the sex
            dat$SEX[dat$SEX == 1] <- "Male"
            dat$SEX[dat$SEX == 2] <- "Female"
            dat$SEX <- as.factor(as.character(dat$SEX))

            # Merge the files
            dat <- merge(dat, pgs, by = c("FID", "IID"))

            # Basic density plot with custom color
            p2 <- createSexDistributionPlot(dat)

            if (pheno_type == "binary") {
                plot_out <- createBinaryPhenotypePlots(dat, p1, p2)
            } else {
                plot_out <- ggpubr::ggarrange(p1, p2)
            }

            # Define patterns for files to be removed
            filePatternsToRemove <- c("pruned_", "PGS", "Clump", "pcfile")

            # Remove files for each pattern
            for (pattern in filePatternsToRemove) {
                ftemp <- list.files(ResultDir, pattern = pattern)
                removeFiles(ftemp, ResultDir)
            }

            # Additional specific files to remove
            additionalFilesToRemove <- c("range_list", "Valid.SNP", "SNPdata_1", "SNP.pvalue", "pgssummarystat")
            removeFiles(additionalFilesToRemove, ResultDir)
        },
        error = function(e) {
            message("An error occurred: ", e$message)
            return(NULL)
        },
        warning = function(w) {
            message("Warning: ", w$message)
        }
    )

    return(list(PGS = pheno.pgs, BestP = bestP, GeneticPC_plot = GP$plot, PGS_plot = plot_out))
}

## Function 144
## Added in 3.0
validateInputForGeneticCorrBT <- function(DataDir, ResultDir, finput, byCHR, REMLalgo, nitr, phenofile, cat_covarfile, quant_covarfile, partGRM, autosome, Xsome, nGRM, cripticut, minMAF, maxMAF, excludeResidual, ncores) {
    # Validate DataDir and ResultDir
    if (!dir.exists(DataDir)) {
        stop("Error in DataDir: Directory does not exist.")
    }
    if (!dir.exists(ResultDir)) {
        stop("Error in ResultDir: Directory does not exist.")
    }

    # Validate finput
    if (!is.character(finput)) {
        stop("Error in finput: Must be a character string.")
    }

    # Validate phenofile - must be a dataframe with exactly four columns
    if (!is.null(phenofile)) {
        if (!is.data.frame(phenofile)) {
            stop("Error in phenofile: Must be a dataframe.")
        }
        if (ncol(phenofile) != 4) {
            stop("Error in phenofile: Dataframe must contain exactly four columns.")
        }
    }


    # Validate categorical and quantitative covariate files
    if (!is.null(cat_covarfile) && !file.exists(normalizePath(file.path(DataDir, cat_covarfile), mustWork = FALSE))) {
        stop("Error in cat_covarfile: Specified file does not exist in DataDir.")
    }
    if (!is.null(quant_covarfile) && !file.exists(normalizePath(file.path(DataDir, quant_covarfile), mustWork = FALSE))) {
        stop("Error in quant_covarfile: Specified file does not exist in DataDir.")
    }

    # Validate boolean parameters
    if (!is.logical(byCHR) || !is.logical(partGRM) || !is.logical(autosome) || !is.logical(Xsome) || !is.logical(excludeResidual)) {
        stop("Error in boolean parameters: All must be TRUE or FALSE.")
    }

    # Validate REMLalgo
    if (!is.numeric(REMLalgo) || !all(REMLalgo %in% c(0, 1, 2))) {
        stop("Error in REMLalgo: Must be 0, 1, or 2.")
    }

    # Validate nitr
    if (!is.numeric(nitr) || nitr <= 0 || nitr != as.integer(nitr)) {
        stop("Error in nitr: Must be a positive whole number.")
    }

    # Validate nGRM
    if (!is.numeric(nGRM) || nGRM <= 0 || nGRM != as.integer(nGRM)) {
        stop("Error in nGRM: Must be a positive whole number.")
    }

    # Validate cripticut, minMAF, maxMAF
    if (!is.null(cripticut) && (!is.numeric(cripticut) || cripticut < 0 || cripticut > 1)) {
        stop("Error in cripticut: Must be within the range [0, 1].")
    }
    if (!is.null(minMAF) && (!is.numeric(minMAF) || minMAF < 0 || minMAF > 1)) {
        stop("Error in minMAF: Must be within the range [0, 1].")
    }
    if (!is.null(maxMAF) && (!is.numeric(maxMAF) || maxMAF < 0 || maxMAF > 1)) {
        stop("Error in maxMAF: Must be within the range [0, 1].")
    }

    # Validate ncores
    if (!is.numeric(ncores) || ncores <= 0 || ncores != as.integer(ncores)) {
        stop("Error in ncores: Must be a positive whole number.")
    }

    return(TRUE)
}

#' GeneticCorrBT: Computing genetic correlation between two traits.
#'
#' @description
#' This function computes genetic correlation, a quantitative genetic measure that describes the genetic link between two
#' traits and has been predicted to indicate pleiotropic gene activity or correlation between causative loci in two traits.
#' For example, it does a bivariate GREML analysis to determine the genetic association between two quantitative traits,
#' two binary disease traits from case-control studies, and between a quantitative trait and a binary disease trait
#' following \insertCite{Yang2011,Lee2012}{GXwasR}. If users want, this function gives the flexibility to compute the genetic
#' correlation chromosome-wise.
#'
#' @param DataDir
#' A character string for the file path of the all the input files.
#'
#' @param ResultDir
#' A character string for the file path where all output files will be stored. The default is `tempdir()`.
#'
#' @param finput
#' Character string, specifying the prefix of the input PLINK binary files for the genotype data. This file needs to be in `DataDir`.
#'
#' @param byCHR
#' Boolean value, `TRUE` or `FALSE`, specifying whether the analysis will be performed chromosome wise or not. The default is `FALSE`.
#'
#' @param REMLalgo
#' Integer value of 0, 1 or 2, specifying the algorithm to run REML iterations, 0 for average information (AI), 1 for Fisher-scoring and
#' 2 for EM. The default option is 0, i.e. AI-REML (1).
#'
#' @param nitr
#' Integer value, specifying the number of iterations for performing the REML. The default is 100.
#'
#' @param phenofile
#' A dataframe for Bivar RELM has four columns `family ID`, `individual ID` and two trait columns. For binary trait, the phenotypic value
#' should be coded as 0 or 1, then it will be recognized as a case-control study (0 for controls and 1 for cases). Missing value should
#' be represented by "-9" or "NA".
#'
#' @param cat_covarfile
#' A character string, specifying the name of the categorical covariate file which is a plain text file with no header line; columns are
#' `family ID`, `individual ID` and discrete covariates. The default is `NULL`. This file needs to be in `DataDir`.
#'
#' @param quant_covarfile
#' A character string, specifying the name of the quantitative covariate file which is a plain text file with no header line; columns
#' are `family ID`, `individual ID` and continuous covariates. The default is `NULL`. This file needs to be in `DataDir`.
#'
#' @param computeGRM
#' Boolean value, `TRUE` or `FALSE`, specifying whether to compute GRM matrices or not. The default is `TRUE`.
#'
#' @param grmfile_name
#' A string of characters specifying the prefix of autosomal .grm.bin file. Users need to provide separate GRM files for autosomes
#' and X chromosome in `ResultDir`. The X chromosomal GRM file should have "x" added in the autosomal prefix as file name.
#'
#' For instance, if autosomal file is "ABC.grm.bin", then X chromosomal file should be "xABC.grm.bim". If you are providing
#' chromosome-wise GRMs, then the prefix should add "ChrNumber_" at the starting of the prefix like, "Chr1_ABC.grm.bin".
#' The default is `NULL`.
#'
#' @param partGRM
#' Boolean value, `TRUE` or `FALSE`, specifying whether the GRM will be partitioned into n parts (by row) in GREML model. The default is `FALSE`.
#'
#' @param autosome
#' Boolean value, `TRUE` or `FALSE`, specifying whether estimate of heritability will be done for autosomes or not. The default is `TRUE`.
#'
#' @param Xsome
#' Boolean value, `TRUE` or `FALSE`, specifying whether estimate of heritability will be done for X chromosome or not. The default is `TRUE`.
#'
#' @param nGRM
#' Integer value, specifying the number of the partition of the GRM in GREML model. The default is 3.
#'
#' @param cripticut
#' Numeric value, specifying the threshold to create a new GRM of "unrelated" individuals in GREML model. The default is arbitrary
#' chosen as 0.025 following \insertCite{Yang2011}{GXwasR}.
#'
#' @param minMAF
#' Positive numeric value (< maxMAF), specifying the minimum threshold for the MAF filter of the SNPs in the Bivariate GREML model.
#'
#' @param maxMAF
#' Positive numeric value (minMAF,1), specifying the maximum threshold for the MAF filter of the SNPs in the Bivariate GREML model.
#'
#' @param excludeResidual
#' Boolean value, `TRUE` or `FALSE`, specifying whether to drop the residual covariance from the model. Recommended to set this `TRUE`
#' if the traits were measured on different individuals. The default is `FALSE`.
#'
#' @param ncores
#' Integer value, specifying the number of cores to be used.
#'
#' @importFrom dplyr everything
#' @importFrom stringr str_remove
#' @importFrom tidyr pivot_longer
#'
#' @return
#' A dataframe with minimum three columns:
#'
#' * Source" (i.e., source of heritability)
#' * Variance" (i.e. estimated heritability)
#' * SE" (i.e., standard error of the estimated heritability)
#'
#' Source column will have rows, such as V(G)_tr1 (genetic variance for trait 1), V(G)_tr2 (genetic variance for trait 2),
#' C(G)_tr12 (genetic covariance between traits 1 and 2),V(e)_tr1 (residual variance for trait 1), V(e)_tr2 (residual variance for trait 2),
#' C(e)_tr12 (residual covariance between traits 1 and 2), Vp_tr1 (proportion of variance explained by all SNPs for trait 1),
#' Vp_tr2 (proportion of variance explained by all SNPs for trait 2), V(G)/Vp_tr1 (phenotypic variance for trait 1),
#' V(G)/Vp_tr2 (phenotypic variance for trait 2), rG (genetic correlation) and n (sample size). In case of chromosome-wise
#' analysis, there will be 'chromosome' column for chromosome code.
#'
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' data("Example_phenofile", package = "GXwasR")
#' DataDir <- GXwasR:::GXwasR_data()
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' byCHR <- TRUE
#' REMLalgo <- 0
#' nitr <- 3
#' ncores <- 3
#' phenofile <- Example_phenofile # Cannot be NULL
#' cat_covarfile <- NULL
#' quant_covarfile <- NULL
#' partGRM <- FALSE # Partition the GRM into m parts (by row),
#' autosome <- TRUE
#' Xsome <- TRUE
#' cripticut <- 0.025
#' minMAF <- 0.01 # if MAF filter apply
#' maxMAF <- 0.04
#' excludeResidual <- TRUE
#'
#' genetic_correlation <- GeneticCorrBT(
#'     DataDir = DataDir, ResultDir = ResultDir, finput = finput, byCHR = byCHR,
#'     REMLalgo = 0, nitr = 10, phenofile = phenofile, cat_covarfile = NULL, quant_covarfile = NULL,
#'     partGRM = FALSE, autosome = TRUE, Xsome = TRUE, nGRM = 3,
#'     cripticut = 0.025, minMAF = NULL, maxMAF = NULL, excludeResidual = TRUE, ncores = ncores
#' )
GeneticCorrBT <- function(
      DataDir, ResultDir, finput, byCHR = FALSE,
      REMLalgo = c(0, 1, 2), nitr = 100, phenofile, cat_covarfile = NULL, quant_covarfile = NULL,
      computeGRM = TRUE, grmfile_name = NULL,
      partGRM = FALSE, autosome = TRUE, Xsome = TRUE, nGRM = 3,
      cripticut = 0.025, minMAF = NULL, maxMAF = NULL,
      excludeResidual = FALSE, ncores = 2
) {
    # Validate input parameters
    validateInputForGeneticCorrBT(
        DataDir, ResultDir, finput, byCHR, REMLalgo, nitr, phenofile,
        cat_covarfile, quant_covarfile, partGRM, autosome, Xsome, nGRM,
        cripticut, minMAF, maxMAF, excludeResidual, ncores
    )

    if (!checkFiles(DataDir, finput)) {
        stop("Missing required Plink files in the specified DataDir.")
    }

    tryCatch(
        withCallingHandlers(
            {
                ## ComputeBivarREMLone phenofile
                write.table(phenofile, file = normalizePath(file.path(ResultDir, "GCphenofile"), mustWork = FALSE), row.names = FALSE, quote = FALSE)

                if (byCHR == FALSE) {
                    if (autosome == TRUE && Xsome == FALSE) {
                        ## Compute GRM
                        if (computeGRM == TRUE) {
                            ComputeGRMauto(
                                DataDir = DataDir, ResultDir = ResultDir, finput = finput,
                                partGRM = partGRM, nGRM = nGRM, cripticut = cripticut, minMAF = minMAF, maxMAF = maxMAF, ncores = ncores
                            )

                            grmfile_name <- "GXwasR"
                        } else {
                            grmfile_name <- grmfile_name
                        }
                        ## Compute REML
                        herit_result <- ComputeBivarREMLone(
                            DataDir = DataDir, ResultDir = ResultDir, REMLalgo = REMLalgo, nitr = nitr, phenofile = "GCphenofile", cat_covarfile = cat_covarfile,
                            quant_covarfile = quant_covarfile, excludeResidual = excludeResidual, chr = "chromosome", grmfile = grmfile_name, ncores = ncores
                        )

                        return(herit_result)
                    } else if (autosome == TRUE && Xsome == TRUE) {
                        if (computeGRM == TRUE) {
                            ## Compute GRM Autosome
                            ComputeGRMauto(
                                DataDir = DataDir, ResultDir = ResultDir, finput = finput,
                                partGRM = partGRM, nGRM = nGRM, cripticut = cripticut, minMAF = minMAF, maxMAF = maxMAF, ncores = ncores
                            )
                            ## Compute GRM X
                            ComputeGRMX(
                                DataDir = DataDir, ResultDir = ResultDir, finput = finput,
                                partGRM = partGRM, nGRM = nGRM, minMAF = minMAF, maxMAF = maxMAF, ncores = ncores
                            )
                        } else {
                            rlang::inform(rlang::format_error_bullets(c("i" = "Skipping GRM computation.")))
                        }

                        ## Compute REML
                        herit_result <- ComputeBivarREMLmulti(
                            DataDir = DataDir, ResultDir = ResultDir, REMLalgo = REMLalgo, nitr = nitr, phenofile = "GCphenofile", cat_covarfile = cat_covarfile,
                            quant_covarfile = quant_covarfile, excludeResidual = excludeResidual, grmfile = "multi_GRMs.txt", computeGRM = computeGRM, grmfile_name = grmfile_name, ncores = ncores
                        )

                        return(herit_result)
                    } else if (autosome == FALSE && Xsome == TRUE) {
                        if (computeGRM == TRUE) {
                            ## Compute GRM X
                            ComputeGRMX(
                                DataDir = DataDir, ResultDir = ResultDir, finput = finput,
                                partGRM = partGRM, nGRM = nGRM, minMAF = minMAF, maxMAF = maxMAF, ncores = ncores
                            )


                            grmfile_name <- "xGXwasR"
                        } else {
                            grmfile_name <- paste0("x", grmfile_name)
                        }
                        ## Compute REML X
                        herit_result <- ComputeBivarREMLone(
                            DataDir = DataDir, ResultDir = ResultDir, REMLalgo = REMLalgo, nitr = nitr, phenofile = "GCphenofile", cat_covarfile = cat_covarfile,
                            quant_covarfile = quant_covarfile, excludeResidual = excludeResidual, chr = "chromosome", grmfile = grmfile_name, ncores = ncores
                        )

                        return(herit_result)
                    } else {
                        rlang::inform(rlang::format_error_bullets(c("x" = "autosome and Xsome cannot both be set as FALSE.")))
                    }
                } else {
                    bimfile <- read.table(normalizePath(file.path(DataDir, paste0(finput, ".bim")), mustWork = FALSE))

                    chrnum <- seq_len(length(unique(bimfile$V1)))

                    chrwiseRELM <- function(chrnum, ncores) {
                        chromosome <- as.integer(unique(bimfile$V1)[chrnum])

                        rlang::inform(paste0("Processing chromosome ", chromosome))


                        if (chromosome == 23) {
                            if (computeGRM == TRUE) {
                                ## Compute GRM X

                                ComputeGRMX(
                                    DataDir = DataDir, ResultDir = ResultDir, finput = finput,
                                    partGRM = partGRM, nGRM = nGRM, minMAF = minMAF, maxMAF = maxMAF, ncores = ncores
                                )

                                grmfile_name <- "xGXwasR"
                            } else {
                                grmfile_name <- paste0("x", grmfile_name)
                            }

                            ## Compute REML X
                            herit_result <- ComputeBivarREMLone(
                                DataDir = DataDir, ResultDir = ResultDir, REMLalgo = REMLalgo, nitr = nitr, phenofile = "GCphenofile", cat_covarfile = cat_covarfile,
                                quant_covarfile = quant_covarfile, excludeResidual = excludeResidual, chr = 23, grmfile = grmfile_name, ncores = ncores
                            )

                            herit_result <- data.table::as.data.table(cbind(chromosome, herit_result))
                        } else {
                            if (computeGRM == TRUE) {
                                ## Compute GRM
                                ComputeGRMauto(
                                    DataDir = DataDir, ResultDir = ResultDir, finput = finput,
                                    partGRM = partGRM, nGRM = nGRM, cripticut = cripticut,
                                    minMAF = minMAF, maxMAF = maxMAF, ByCHR = byCHR, CHRnum = chromosome, ncores = ncores
                                )

                                grmfile_name <- "GXwasR"
                            } else {
                                grmfile_name <- paste0("Chr", chromosome, "_", grmfile_name)
                            }

                            ## Compute REML
                            herit_result <- ComputeBivarREMLone(
                                DataDir = DataDir, ResultDir = ResultDir, REMLalgo = REMLalgo, nitr = nitr, phenofile = "GCphenofile", cat_covarfile = cat_covarfile,
                                quant_covarfile = quant_covarfile, excludeResidual = excludeResidual, chr = chromosome, grmfile = grmfile_name, ncores = ncores
                            )


                            herit_result <- data.table::as.data.table(cbind(chromosome, herit_result))
                        }

                        return(herit_result)
                    }
                    result <- data.table::rbindlist(lapply(chrnum, function(x) chrwiseRELM(x, ncores)), fill = TRUE)
                    # result <- na.omit(result)

                    # Gather files matching the patterns
                    patterns <- c("bireml", "grm", "test", "gcta", "GCphenofile", "multi_GRMs.txt")
                    patterns_regex <- paste0(patterns, collapse = "|")
                    files_to_remove <- list.files(ResultDir, pattern = patterns_regex, full.names = TRUE)
                    # Remove Files
                    invisible(file.remove(files_to_remove))

                    return(result)
                }
            },
            error = function(e) {
                message("An error occurred: ", e$message)
                return(NULL)
            },
            warning = function(w) {
                message("Warning: ", w$message)
                invokeRestart("muffleWarning")
            }
        )
    )
}


## Function 129
## Added in 3.0
validateInputForEstimateHerit <- function(DataDir = NULL, ResultDir = tempdir(), finput = NULL,
    summarystat = NULL, ncores = parallel::detectCores(),
    model = c("LDSC", "GREML"), byCHR = FALSE, r2_LD = 0,
    LDSC_blocks = 20, REMLalgo = c(0, 1, 2), nitr = 100,
    cat_covarfile = NULL, quant_covarfile = NULL,
    prevalence = 0.01, partGRM = FALSE, autosome = TRUE,
    Xsome = TRUE, nGRM = 3, cripticut = 0.025,
    minMAF = NULL, maxMAF = NULL, hg = c("hg19", "hg38"),
    PlotIndepSNP = c(TRUE, FALSE), IndepSNP_window_size = 50,
    IndepSNP_step_size = 5, IndepSNP_r2_threshold = 0.02,
    highLD_regions = NULL) {
    # Validate directories
    if (!is.null(DataDir) && !dir.exists(DataDir)) {
        stop("Error in DataDir: Directory does not exist.")
    }
    if (!is.null(ResultDir) && !dir.exists(ResultDir)) {
        stop("Error in ResultDir: Directory does not exist.")
    }


    # Validate summarystat if not NULL
    if (!is.null(summarystat) && !is.data.frame(summarystat)) {
        stop("Error in summarystat: Must be a dataframe.")
    }

    # Validate ncores
    if (!is.numeric(ncores) || ncores < 0) {
        stop("Error in ncores: Must be a non-negative integer.")
    }

    # Validate model
    if (!is.character(model) || !all(model %in% c("LDSC", "GREML"))) {
        stop("Error in model: Must be either 'LDSC' or 'GREML', or both.")
    }

    # Validate file prefix
    if (model == "GREML") {
        if (is.null(finput) || !is.character(finput)) {
            stop("Error in finput: Must be a non-null character string.")
        }
    }

    # Validate boolean parameters
    if (!is.logical(byCHR) || !is.logical(partGRM) || !is.logical(autosome) ||
        !is.logical(Xsome) || !is.logical(PlotIndepSNP)) {
        stop("Error in Boolean parameters: Must be TRUE or FALSE.")
    }

    # Validate numerical parameters
    if (!is.numeric(r2_LD) || r2_LD < 0) {
        stop("Error in r2_LD: Must be a non-negative number.")
    }
    if (!is.null(prevalence) && (!is.numeric(prevalence) || prevalence < 0)) {
        stop("Error in prevalence: Must be a non-negative number or NULL.")
    }
    if (!is.numeric(cripticut) || cripticut < 0) {
        stop("Error in cripticut: Must be a non-negative number.")
    }

    # Validate integer-like parameters
    if (!is.null(LDSC_blocks) && (!is.numeric(LDSC_blocks) || LDSC_blocks <= 0 ||
        LDSC_blocks != as.integer(LDSC_blocks))) {
        stop("Error in LDSC_blocks: Must be a positive whole number.")
    }
    if (!is.null(nitr) && (!is.numeric(nitr) || nitr <= 0 || nitr != as.integer(nitr))) {
        stop("Error in nitr: Must be a positive whole number.")
    }
    if (!is.null(nGRM) && (!is.numeric(nGRM) || nGRM <= 0 || nGRM != as.integer(nGRM))) {
        stop("Error in nGRM: Must be a positive whole number.")
    }
    if (!is.numeric(IndepSNP_window_size) || IndepSNP_window_size <= 0 ||
        IndepSNP_window_size != as.integer(IndepSNP_window_size)) {
        stop("Error in IndepSNP_window_size: Must be a positive whole number.")
    }
    if (!is.numeric(IndepSNP_step_size) || IndepSNP_step_size <= 0 ||
        IndepSNP_step_size != as.integer(IndepSNP_step_size)) {
        stop("Error in IndepSNP_step_size: Must be a positive whole number.")
    }

    # Validate REMLalgo
    if (!is.numeric(REMLalgo) || !all(REMLalgo %in% c(0, 1, 2))) {
        stop("Error in REMLalgo: Must be 0, 1, or 2.")
    }

    # Validate MAF parameters
    if ((!is.null(minMAF) && (!is.numeric(minMAF) || minMAF < 0 || minMAF > 1)) ||
        (!is.null(maxMAF) && (!is.numeric(maxMAF) || maxMAF < 0 || maxMAF > 1))) {
        stop("Error in MAF parameters: Must be within the range [0, 1].")
    }

    # Validate hg parameter
    if (!is.character(hg) || !all(hg %in% c("hg19", "hg38"))) {
        stop("Error in hg: Must be 'hg19', 'hg38', or both.")
    }

    # Validate highLD_regions if not NULL
    if (!is.null(highLD_regions) && !is.data.frame(highLD_regions)) {
        stop("Error in highLD_regions: Must be a dataframe.")
    }

    return(TRUE)
}

#' EstimateHerit: Computing SNP heritability i.e., the proportion of phenotypic variance explained by SNPs.
#'
#' @author Banabithi Bose
#'
#' @description This functions performs two types of heritability estimation, (i)GREML:Genomic relatedness matrix (GRM) restricted
#' maximum likelihood-based method following GCTA \insertCite{Yang2011}{GXwasR} and (ii)LDSC: LD score regression-based method
#' following \insertCite{Bulik-Sullivan2014,Prive2020}{GXwasR}. For the details, please follow the associated paper.
#'
#' Prior to using this function, it is recommended to apply QCsnp and QCsample to ensure data quality control.
#'
#' @param DataDir
#' A character string for the file path of the all the input files. The default is `NULL`.
#'
#' @param ResultDir
#' A character string for the file path where all output files will be stored. The default is `tempdir()`.
#'
#' @param finput
#' Character string, specifying the prefix of the input PLINK binary files for the genotype data. This file needs to be in `DataDir`.
#' For LDSC model, if the original genotype data is not available, Hapmap 3 or 1000Genome data can be used. If use NULL, then you need
#' to provide `precomputedLD` argument. See below.
#'
#' @param precomputedLD
#' A dataframe object as LD matrix with columns: `CHR`, `SNP`, `BP`, `ld_size`, `MAF`, `ld_score`. . The default is `NULL`.
#'
#' @param chi2_thr1
#' Numeric value for threshold on chi2 in step 1 of LDSC regression. Default is 30.
#'
#' @param chi2_thr2
#' Numeric value for threshold on chi2 in step 2. Default is `Inf` (none).
#'
#' @param intercept
#' Numeric value to constrain the intercept to some value (e.g. 1) in LDSC regression. Default is `NULL`.
#'
#' @param summarystat
#' A dataframe object with GWAS summary statistics. The mandatory column headers in this dataframe are
#' * `chr` (Chromosome code),
#' * `pos` (Basepair position)
#' * `a1` (First allele code)
#' * `rsid` (i.e., SNP identifier)
#' * `beta` (i.e., effect-size or logarithm of odds ratio)
#' * `beta_se` (i.e., standard error of beta)
#' * `P` (i.e., p-values)
#' * `n_eff` (i.e., effective sample size)
#'
#' For case-control study, effective sample size should be \eqn{4 / (1/<qty of cases> + 1/<qty of controls>)}. The default is `NULL`.
#'
#' @param ncores
#' Integer value, specifying the number of cores to be used for running LDSC model. The default is 2.
#'
#' @param model
#' Character string, specifying the heritability estimation model. There are two options, “GREML” or “LDSC”. The default is “GREML”.
#'
#' Note: argument For LDSC, DataDir and finput can be `NULL`.
#'
#' @param byCHR
#' Boolean value, `TRUE` or `FALSE`, specifying whether the analysis will be performed chromosome wise or not. The default is `FALSE`.
#'
#' @param r2_LD
#' Numeric value, specifying the LD threshold for clumping in LDSC model. The default is 0.
#'
#' @param LDSC_blocks
#' Integer value, specifying the block size for performing jackknife variance estimator in LDSC model following \insertCite{Prive2020}{GXwasR}.
#' The default is 200.
#'
#' @param REMLalgo
#' Integer value of 0, 1 or 2, specifying the algorithm to run REML iterations, 0 for average information (AI), 1 for Fisher-scoring and
#' 2 for EM. The default option is 0, i.e. AI-REML \insertCite{Yang2011}{GXwasR}.
#'
#' @param nitr
#' Integer value, specifying the number of iterations for performing the REML. The default is 100.
#'
#' @param cat_covarfile
#' A character string, specifying the name of the categorical covariate file which is a plain text file with no header line; columns
#' are family ID, individual ID and discrete covariates. The default is `NULL`. This file needs to be in `DataDir`.
#'
#' @param quant_covarfile
#' A character string, specifying the name of the quantitative covariate file which is a plain text file with no header line;
#' columns are family ID, individual ID and continuous covariates. The default is `NULL`. This file needs to be in `DataDir`.
#'
#' @param prevalence
#' Numeric value, specifying the disease prevalence. The default is `NULL`.
#'
#' Note: for the continuous trait value, users should use the default.
#'
#' @param computeGRM
#' Boolean value, `TRUE` or `FALSE`, specifying whether to compute GRM matrices or not. The default is `TRUE`.
#'
#' @param grmfile_name
#' A string of characters specifying the prefix of autosomal .grm.bin file. Users need to provide separate GRM files
#' for autosomes and X chromosome in `ResultDir`.
#'
#' The X chromosomal GRM file should have "x" added in the autosomal prefix as file name. For instance, if autosomal file
#' is "ABC.grm.bin", then X chromosomal file should be "xABC.grm.bim".
#'
#' If you are providing chromosome-wise GRMs, then the prefix should add "ChrNumber_" at the start of the prefix like,
#' "Chr1_ABC.grm.bin". The default is NULL.
#'
#' @param partGRM
#' Boolean value, `TRUE` or `FALSE`, specifying whether the GRM will be partitioned into n parts (by row) in GREML model. The default is `FALSE`.
#'
#' @param autosome
#' Boolean value, `TRUE` or `FALSE`, specifying whether estimate of heritability will be done for autosomes or not. The default is `TRUE`.
#'
#' @param Xsome
#' Boolean value, `TRUE` or `FALSE`, specifying whether estimate of heritability will be done for X chromosome or not. The default is `TRUE`.
#'
#' @param nGRM
#' Integer value, specifying the number of the partition of the GRM in GREML model. The default is 3.
#'
#' @param cripticut
#' Numeric value, specifying the threshold to create a new GRM of "unrelated" individuals in GREML model. The default is arbitrary chosen
#' as 0.025 following \insertCite{Yang2011}{GXwasR}.
#'
#' @param minMAF
#' Positive numeric value (0,1), specifying the minimum threshold for the MAF filter of the SNPs in the GREML model. This value cannot be
#' greater than `maxMAF` parameter. The default is `NULL`. For `NULL`, maximum MAF value of the genotype data will be computed and printed
#' on the plot.
#'
#' @param maxMAF
#' Positive numeric value (0,1), specifying the maximum threshold for the MAF filter of the SNPs in the GREML model. This value cannot be less
#' than `minMAF` parameter. The default is `NULL`. For `NULL`, minimum MAF value of the genotype data will be computed and printed on the plot.
#'
#' @param hg
#' Boolean value, specifying the genome built, “hg19” or “hg38” to use chromosome length from UCSC genome browser and getting genes and proteins
#' according to this built. The default is “hg19”.
#'
#' @param PlotIndepSNP
#' Boolean value, `TRUE` or `FALSE`, specifying whether to use independent SNPs i.e., chromosome-wise LD pruned SNPs in the plots or not.
#' The default is `TRUE`.
#'
#' @param indepSNPs
#' A dataframe with independent SNP ids with column name "rsid". The default is `NULL`.
#'
#' @param IndepSNP_window_size
#' Integer value, specifying a window size in variant count or kilobase for LD-based filtering. The default is 50.
#'
#' @param IndepSNP_step_size
#' Integer value, specifying a variant count to shift the window at the end of each step for LD filtering for pruned SNPs in the plots.
#' The default is 5.
#'
#' @param IndepSNP_r2_threshold
#' Numeric value between 0 to 1 of pairwise \eqn{r^2} threshold for LD-based filtering for pruned SNPs in the plots. The default is 0.02.
#'
#' @param highLD_regions
#' Character string, specifying the .txt file name with genomic regions with high LD for using in finding pruned SNPs in the plots.
#' This file needs to be in `DataDir`.
#'
#' @param plotjpeg
#' Boolean value, `TRUE` or `FALSE`, specifying whether to save the plots in jpeg file in `ResultDir`. The default is `TRUE`.
#'
#' @param plotname
#' String of character value specifying the name of the jpeg file with the plots. The default is "Heritability_Plots".
#'
#' @returns
#' A dataframe with maximum eight columns for GREML (here, three columns if running genome-wide) and ten columns for LDSC model if byCHR is `TRUE`.
#' The columns, such as, "chromosome"(i.e., chromosome code),"snp_proportion" (i.e.,chromosome-wise SNP proportion)", "no.of.genes" (i.e., number of genes per chromosome),
#' "no.of.proteins" (i.e., number of genes per chromosome),"size_mb" (i.e., chromosome length), "Source" (i.e., source of heritability),
#' "Variance" (i.e., estimated heritability), and "SE" (i.e., standard error of the estimated heritability) are common for both GREML and LDSC model.
#' The column, "Intercept" (i.e., LDSC regression intercept) and "Int_SE" (i.e., standard error of the intercept) will be two extra columns for LDSC models.
#' Source column will have rows, such as `V(1)` (i.e., name of genetic variance), `V(e)` (i.e., residual variance), `V(p)` (i.e., phenotypic variance), `V(1)/Vp` (i.e.,
#' ratio of genetic variance to phenotypic variance), and `V(1)/Vp_L` (i.e., ratio of genetic variance to phenotypic variance in liability scale for binary phenotypes).
#' If `byCHR` is `FALSE`, then the first five columns will not be reported in the dataframe.
#'
#' @references
#' \insertAllCited{}
#'
#' @importFrom bigsnpr snp_readBed snp_attach snp_match coef_to_liab
#' @importFrom data.table as.data.table rbindlist
#'
#' @export
#'
#' @examples
#' data("Summary_Stat_Ex1", package = "GXwasR")
#' data("highLD_hg19", package = "GXwasR")
#' DataDir <- GXwasR:::GXwasR_data()
#' ResultDir <- tempdir()
#' precomputedLD <- NULL
#' finput <- "GXwasR_example"
#' test.sumstats <- na.omit(Summary_Stat_Ex1[Summary_Stat_Ex1$TEST == "ADD", c(seq_len(4), 6:8)])
#' colnames(test.sumstats) <- c("chr", "rsid", "pos", "a1", "n_eff", "beta", "beta_se")
#' summarystat <- test.sumstats
#' ncores <- 3
#' model <- "GREML"
#' byCHR <- FALSE
#' r2_LD <- 0
#' LDSC_blocks <- 20
#' REMLalgo <- 0
#' nitr <- 3
#' cat_covarfile <- NULL
#' quant_covarfile <- NULL
#' prevalence <- 0.01
#' partGRM <- FALSE
#' autosome <- TRUE
#' Xsome <- TRUE
#' nGRM <- 3
#' cripticut <- 0.025
#' minMAF <- NULL
#' maxMAF <- NULL
#' hg <- "hg19"
#' PlotIndepSNP <- TRUE
#' IndepSNP_window_size <- 50
#' IndepSNP_step_size <- 5
#' IndepSNP_r2_threshold <- 0.02
#' highLD_regions <- highLD_hg19
#' H2 <- EstimateHerit(
#'     DataDir = DataDir, ResultDir = ResultDir, finput = finput,
#'     summarystat = NULL, ncores, model = "GREML", byCHR = TRUE, r2_LD = 0,
#'     LDSC_blocks = 20, REMLalgo = 0, nitr = 100, cat_covarfile = NULL, quant_covarfile = NULL,
#'     prevalence = 0.01, partGRM = FALSE, autosome = TRUE, Xsome = TRUE, nGRM = 3,
#'     cripticut = 0.025, minMAF = NULL, maxMAF = NULL, hg = "hg19", PlotIndepSNP = TRUE,
#'     IndepSNP_window_size = 50, IndepSNP_step_size = 5, IndepSNP_r2_threshold = 0.02,
#'     highLD_regions = highLD_hg19
#' )
EstimateHerit <- function(DataDir = NULL, ResultDir = tempdir(), finput = NULL, precomputedLD = NULL,
    indepSNPs = NULL, summarystat = NULL, ncores = 2, model = c("LDSC", "GREML"),
    computeGRM = TRUE, grmfile_name = NULL, byCHR = FALSE,
    r2_LD = 0, LDSC_blocks = 20, intercept = NULL, chi2_thr1 = 30,
    chi2_thr2 = Inf, REMLalgo = c(0, 1, 2), nitr = 100, cat_covarfile = NULL,
    quant_covarfile = NULL, prevalence = NULL, partGRM = FALSE, autosome = TRUE,
    Xsome = TRUE, nGRM = 3, cripticut = 0.025, minMAF = NULL, maxMAF = NULL,
    hg = c("hg19", "hg38"), PlotIndepSNP = TRUE, IndepSNP_window_size = 50,
    IndepSNP_step_size = 5, IndepSNP_r2_threshold = 0.02, highLD_regions = NULL,
    plotjpeg = TRUE, plotname = "Heritability_Plots") {
    # Validate inputs
    if (!validateInputForEstimateHerit(DataDir, ResultDir, finput, summarystat, ncores, model, byCHR, r2_LD, LDSC_blocks, REMLalgo, nitr, cat_covarfile, quant_covarfile, prevalence, partGRM, autosome, Xsome, nGRM, cripticut, minMAF, maxMAF, hg, PlotIndepSNP, IndepSNP_window_size, IndepSNP_step_size, IndepSNP_r2_threshold, highLD_regions)) {
        return(NULL)
    }

    if (is.null(precomputedLD)) {
        if (checkFiles(DataDir, finput) == TRUE) {
            rlang::inform(rlang::format_error_bullets(c("v" = "Input genotype files are present in specified directory.")))
        } else {
            stop("Missing required Plink files in the specified DataDir.")
        }
    }


    tryCatch(
        {
            if (is.null(precomputedLD)) {
                maf_range <- computeMAFRange(DataDir, ResultDir, finput, minMAF, maxMAF)
                miMAF <- maf_range$miMAF
                maMAF <- maf_range$maMAF
            } else {
                miMAF <- minMAF
                maMAF <- maxMAF
            }

            if (model == "LDSC") {
                heritability_results <- processLDSCModel(DataDir, ResultDir, finput, precomputedLD, IndepSNPs = indepSNPs, summarystat, byCHR, r2_LD, LDSC_blocks, chi2_thr1, chi2_thr2, intercept, ncores, prevalence, PlotIndepSNP, highLD_regions, IndepSNP_window_size, IndepSNP_step_size, IndepSNP_r2_threshold, hg, miMAF, maMAF, plotjpeg = plotjpeg, plotname = plotname)

                return(heritability_results)
            } else if (model == "GREML") {
                greml_results <- processGREMLModel(DataDir, ResultDir, finput, byCHR, autosome, Xsome, partGRM, nGRM, computeGRM, grmfile_name, cripticut, minMAF, maxMAF, REMLalgo, nitr, cat_covarfile, quant_covarfile, prevalence, PlotIndepSNP, highLD_regions, IndepSNP_window_size, IndepSNP_step_size, IndepSNP_r2_threshold, hg, miMAF, maMAF, ncores, plotjpeg, plotname)

                # Gather files matching the patterns

                patterns_to_remove <- c("PLINK", "test_reml", "LDfiltered", "HumanGenome", "LDsnp", "test", "gcta", "MAF", "GRM", "phenofile.phen", "sink")

                # Use removeFiles helper function to delete the files
                for (pattern in patterns_to_remove) {
                    removeTempFiles(ResultDir, pattern)
                }

                rlang::inform(rlang::format_error_bullets(c("v" = paste0("All GRM related files are in ", ResultDir))))

                return(greml_results)
            }
        },
        error = function(e) {
            message("An error occurred: ", e$message)
            return(NULL)
        },
        warning = function(w) {
            message("Warning: ", w$message)
        }
    )
}

