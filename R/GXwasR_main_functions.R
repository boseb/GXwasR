#' AncestryCheck: Evaluation of the samples' ancestry label.
#'
#' @author Banabithi Bose
#'
#' @description
#' This function displays the result of the ancestry analysis in a color-coded scatter plot of the first two principal
#' components for samples of the reference populations and the study population. Specifically, it compares the study samples'
#' ancestry labels to a panel representing a reference population, and it also flags the outlier samples with respect to a
#' chosen reference population.
#'
#' The function first filters the reference and study data for non-A-T or G-C SNPs. It next conducts
#' LD pruning, fixes the chromosome mismatch between the reference and study datasets, checks for allele flips, updates the
#' locations, and flips the alleles. The two datasets are then joined, and the resulting genotype dataset is subjected to Principal
#' Component Analysis (PCA).
#'
#' The detection of population structure down to the level of the reference dataset can then be accomplished
#' using PCA on this combined genotyping panel. For instance, the center of the European reference samples is determined using the
#' data from principal components 1 and 2 (median(PC1 europeanRef), median(PC2 europeanRef)). It determines the European reference
#' samples' maximum Euclidean distance (maxDist) from this center.
#'
#' All study samples that are non-European, or outliers, are those whose Euclidean distances from the center are more than or
#' equal to the radius r= outlier threshold* maxDist. This function utilizes the HapMap phase 3 data in NCBI 36 and 1000GenomeIII
#' in CGRCh37. Both study and reference datasets should be of the same genome build. If not, users need to lift over one of the
#' datasets to the same build.
#'
#' @param DataDir
#' A character string for the file path of the input PLINK binary files.
#'
#' @param ResultDir
#' A character string for the file path where all output files will be stored. The default is `tempdir()`.
#'
#' @param finput
#' Character string, specifying the prefix of the input PLINK binary files for the study samples.
#'
#' @param reference
#' Boolean value,'HapMapIII_NCBI36' and 'ThousandGenome', specifying Hapmap Phase3 \insertCite{HapMap2010}{GXwasR} and 1000 Genomes
#' phase III \insertCite{1000Genomes2015}{GXwasR} reference population, respectively. The default is 'HapMapIII_NCBI36'.
#'
#' @param filterSNP
#' Boolean value, `TRUE` or `FALSE` for filtering out the SNPs. The default is `TRUE`. We recommend setting it `FALSE`
#' only when the users are sure that they could join the study and reference samples directly.
#'
#' @param studyLD
#' Boolean value, `TRUE` or `FALSE` for applying linkage disequilibrium (LD)-based filtering on study genotype data.
#'
#' @param studyLD_window_size
#' Integer value, specifying a window size in variant count or kilobase for LD-based filtering of the variants for the study data.
#'
#' @param studyLD_step_size
#' Integer value, specifying a variant count to shift the window at the end of each step for LD filtering for the study data.
#'
#' @param studyLD_r2_threshold
#' Numeric value between 0 to 1 of pairwise \eqn{r^2} threshold for LD-based filtering for the study data.
#'
#' @param referLD
#' Boolean value, 'TRUE' or 'FALSE' for applying linkage disequilibrium (LD)-based filtering on reference genotype data.
#'
#' @param referLD_window_size
#' Integer value, specifying a window size in variant count or kilobase for LD-based filtering of the variants for the reference data.
#'
#' @param referLD_step_size
#' Integer value, specifying a variant count to shift the window at the end of each step for LD filtering for the reference data.
#'
#' @param referLD_r2_threshold
#' Numeric value between 0 to 1 of pairwise \eqn{r^2} threshold for LD-based filtering for the reference data.
#'
#' @param highLD_regions
#' A dataframe with known high LD regions \insertCite{Anderson2010}{GXwasR} is provided with the package.
#'
#' @param study_pop
#' A dataframe containing two columns for study in first column, sample ID (i.e., IID) and in second column, the ancestry label.
#'
#' @param outlier
#' Boolean value, `TRUE` or `FALSE`, specifying outlier detection will be performed or not.
#'
#' @param outlierOf
#' Chracter string, specifying the reference ancestry name for detecting outlier samples. The default is "outlierOf = "EUR".
#'
#' @param outlier_threshold
#' Numeric value, specifying the threshold to be be used to detect outlier samples. This threshold will be multiplied with the
#' Eucledean distance from the center of the PC 1 and PC2 to the maximum Euclidean distance of the reference samples. Study samples
#' outside this distance will be considered as outlier. The default is 3.
#'
#' @importFrom data.table as.data.table
#' @importFrom ggplot2 ggplot aes geom_hline geom_vline guides geom_point guide_legend scale_shape_manual
#'
#' @return A dataframe with the IDs of non-European samples as outliers.
#'
#' @references
#' \insertAllCited{}
#'
#' @export
#'
#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' reference <- "HapMapIII_NCBI36"
#' data("GXwasRData")
#' highLD_regions <- GXwasR::highLD_hg19 # See: ?GXwasR::highLD_hg19
#' study_pop <- example_data_study_sample_ancestry # PreimputeEX
#' studyLD_window_size <- 50
#' studyLD_step_size <- 5
#' studyLD_r2_threshold <- 0.02
#' filterSNP <- TRUE
#' studyLD <- TRUE
#' referLD <- TRUE
#' referLD_window_size <- 50
#' referLD_step_size <- 5
#' referLD_r2_threshold <- 0.02
#' outlier <- TRUE
#' outlier_threshold <- 3
#' x <- AncestryCheck(
#'   DataDir = DataDir, ResultDir = ResultDir, finput = finput,
#'   reference = "HapMapIII_NCBI36", highLD_regions = highLD_regions,
#'   study_pop = study_pop, studyLD = studyLD, referLD = referLD,
#'   outlierOf = "EUR", outlier = outlier, outlier_threshold = outlier_threshold
#' )
AncestryCheck <-
  function(DataDir,
           ResultDir = tempdir(),
           finput,
           reference = c("HapMapIII_NCBI36", "ThousandGenome"),
           filterSNP = TRUE,
           studyLD = TRUE,
           studyLD_window_size = 50,
           studyLD_step_size = 5,
           studyLD_r2_threshold = 0.02,
           referLD = FALSE,
           referLD_window_size = 50,
           referLD_step_size = 5,
           referLD_r2_threshold = 0.02,
           highLD_regions,
           study_pop,
           outlier = FALSE,
           outlierOf = "EUR",
           outlier_threshold = 3) {
    tryCatch(
      {
        # Validate inputs
        validateAncestryCheckInputs(DataDir, ResultDir, finput, reference, filterSNP, studyLD, studyLD_window_size, studyLD_step_size, studyLD_r2_threshold, referLD, referLD_window_size, referLD_step_size, referLD_r2_threshold, highLD_regions, study_pop, outlier, outlierOf, outlier_threshold)

        if (!checkFiles(DataDir, finput)) {
          stop("Missing required Plink files in the specified DataDir.")
        }

        # Copying the input data and changing the SNP ids in study data.
        # Copy the file in the same directory with a new name
        file.copy(from = paste0(DataDir, "/", finput, ".fam"), to = paste0(DataDir, "/", finput, "_new.fam"))
        file.copy(from = paste0(DataDir, "/", finput, ".bim"), to = paste0(DataDir, "/", finput, "_new.bim"))
        file.copy(from = paste0(DataDir, "/", finput, ".bed"), to = paste0(DataDir, "/", finput, "_new.bed"))

        sbim <- read.table(paste0(DataDir, "/", finput, "_new.bim"))
        # Assuming sbim is your data frame
        sbim$V2 <- paste(sbim$V1, sbim$V4, sep = ":")
        # Replace the new input .bim file with new snp ids.
        write.table(sbim, file = paste0(DataDir, "/", finput, "_new.bim"), quote = FALSE, row.names = FALSE, col.names = FALSE)

        # Changing the input file
        finput <- paste0(finput, "_new")

        # setupPlink(ResultDir)

        Download_reference(refdata = reference, wdir = ResultDir)

        # Changing the snp ids in reference .bim file
        rbim <- read.table(paste0(ResultDir, "/", reference, ".bim"))
        # Assuming sbim is your data frame
        rbim$V2 <- paste(rbim$V1, rbim$V4, sep = ":")
        # Replace the new input .bim file with new snp ids.
        write.table(rbim, file = paste0(ResultDir, "/", reference, ".bim"), quote = FALSE, row.names = FALSE, col.names = FALSE)


        if (!is.null(highLD_regions)) {
          write.table(highLD_regions, file = paste0(ResultDir, "/high-LD-regions-temp.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

          highLD_regions <- paste0(ResultDir, "/high-LD-regions-temp.txt")
        } else {
          highLD_regions <- NULL
        }

        if (filterSNP == TRUE) {
          # Filter AT-GC
          filterATGCSNPs(DataDir, ResultDir, finput, reference)

          # LD Pruning for Study and Reference Data
          studyLDMessage <- processLDstudyData(ResultDir, highLD_regions, studyLD, studyLD_window_size, studyLD_step_size, studyLD_r2_threshold)
          referenceLDMessage <- processLDreferenceData(ResultDir, highLD_regions, referLD, referLD_window_size, referLD_step_size, referLD_r2_threshold)

          # Printing the messages returned by the functions
          rlang::inform(studyLDMessage)
          rlang::inform(referenceLDMessage)

          # Define patterns for files to remove

          # Remove files with specified patterns
          file.remove(list.files(ResultDir, pattern = "filtered_study_temp1", full.names = TRUE))
          file.remove(list.files(ResultDir, pattern = "filtered_ref_temp1", full.names = TRUE))
        } else if (filterSNP == FALSE) {
          executePlinkForUnfilteredData(DataDir, ResultDir, finput, reference)
        }

        # Find common SNPs between study and reference data
        commonSNPResults <- findCommonSNPs(ResultDir)
        common_snps <- commonSNPResults$common_snps
        pruned_study <- commonSNPResults$pruned_study
        pruned_ref <- commonSNPResults$pruned_ref


        # Process common SNPs
        processCommonSNPs(ResultDir, common_snps)

        # Process for correcting chromosome mismatches
        correctedData <- correctChromosomeMismatches(ResultDir, common_snps, pruned_study, pruned_ref)
        S1 <- correctedData$S1
        S2 <- correctedData$S2
        snpSameNameDiffPos <- correctedData$snpSameNameDiffPos

        # Finding mis-matching allele positions
        updated_ref <-
          read.table(
            file = paste0(ResultDir, "/", "filtered_ref_temp4", ".bim"),
            stringsAsFactors = FALSE
          )
        S3 <-
          updated_ref[match(common_snps, updated_ref[, "V2"]), , drop = FALSE]
        colnames(S1) <- c("V1", "V2", "V3", "V4", "Sa", "Sb")
        colnames(S3) <- c("V1", "V2", "V3", "V4", "Ra", "Rb")

        S1 <- data.table::as.data.table(S1)
        S3 <- data.table::as.data.table(S3)
        # Using SNP name and chr no. for merging, not using base-pair position
        # S4 <- merge(S1, S3, by = c("V1", "V2"))
        # Updating it in final
        S4 <- merge(S1, S3, by = c("V1", "V4")) # using chr and position
        snps_flips <- S4[which(S4[, 5] != S4[, 9] & S4[, 6] != S4[, 10]), ]
        snp_allele_flips <- unique(snps_flips$V2)

        ## Allele Flips
        handleSnpAlleleFlips(ResultDir, snp_allele_flips)

        # Checking allele flips again after correcting
        flipped_ref <-
          read.table(
            file = paste0(ResultDir, "/", "filtered_ref_temp5", ".bim"),
            stringsAsFactors = FALSE
          )
        S5 <-
          flipped_ref[match(common_snps, flipped_ref[, "V2"]), , drop = FALSE]
        colnames(S5) <- c("V1", "V2", "V3", "V4", "Ra", "Rb")
        S5 <- data.table::as.data.table(S5)
        # S6 <- merge(S1, S5, by = c("V1", "V2"))
        # Updating it in final
        S6 <- merge(S1, S5, by = c("V1", "V4"))
        snps_flips_wrong <- S6[which(S6[, 5] != S6[, 9] &
          S6[, 6] != S6[, 10]), ]
        allele_flips_wrong <- unique(snps_flips_wrong$V2)

        # Handel Allele flips wrong
        handleAlleleFlipsWrong(ResultDir, allele_flips_wrong)

        # Checking number of SNPs in reference after clean-up.
        cleaned_ref <-
          read.table(
            file = paste0(ResultDir, "/", "filtered_ref_temp6", ".bim"),
            stringsAsFactors = FALSE
          )

        snps_final <- unique(cleaned_ref$V2)

        # Merge study and reference to perform PCA

        mergeDatasetsAndPerformPCA(ResultDir)

        # Process Reference
        ref_ancestry_EUR_AFR_ASIAN <- loadAndProcessReferenceAncestry(reference)

        # Plot PCA
        combined_pop <- prepareAncestryData(study_pop, ref_ancestry_EUR_AFR_ASIAN)
        tab <- loadPCAData(ResultDir, combined_pop)
        pop_type <- createPopulationTypeData(tab)

        plotPCA(tab, pop_type)

        reportAlleleFlips(snp_allele_flips, ResultDir)

        # Example of using the function
        Outlier_samples1 <- detectOutliers(tab, ResultDir, outlier, outlierOf, outlier_threshold)

        removeTempFiles(ResultDir, "study_ref_merge")
        removeTempFiles(DataDir, "_new")

        # Define patterns for other files to remove
        patterns_to_remove <- c(
          c("study_SNP", "ref_SNP", "common_snps", "snp_allele_flips", "allele_flips_wrong"),
          "Outlier_ancestry",
          "snpSameNameDiffPos",
          reference
        )
        for (pattern in patterns_to_remove) {
          removeTempFiles(ResultDir, pattern)
        }

        return(Outlier_samples1)
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
#' "SKAT" (sequence kernel association test), "SKATO" (combination of BT and SKAT), "sumchi" (sum of χ2-statistics), "ACAT" (aggregated
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
#' @importFrom regioneR toGRanges
#' @importFrom plyranges join_overlap_intersect
#' @importFrom sumFREGAT SKAT SKATO sumchi ACAT BT PCA FLM simpleM minp
#' @importFrom rlang inform format_error_bullets
#'
#' @export
#'
#' @examples
#' \donttest{
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' data("GXwasRData")
#' sumstat <- XWAS_Summary_Example
#' ref_data <- NULL
#' gene_file <- "Xlinkedgenes_hg19.txt"
#' gene_range <- 500000
#' max_gene <- 10
#' gene_approximation <- TRUE
#' beta_par <- c(1, 25)
#' weights_function <- NULL
#' geno_variance_weights <- "se.beta"
#' method <- "kuonen"
#' acc_devies <- 1e-8
#' lim_devies <- 1e+6
#' rho <- TRUE
#' skato_p_threshold <- 0.8
#' mac_threshold <- 3
#' sample_size <- 4000
#' reference_matrix_used <- FALSE
#' regularize_fun <- "LH"
#' pca_var_fraction <- 0.85
#' flm_basis_function <- "fourier"
#' flm_num_basis <- 25
#' flm_poly_order <- 4
#' flip_genotypes <- FALSE
#' omit_linear_variant <- FALSE
#' kernel_p_method <- "kuonen"
#' anno_type <- ""
#' GenetestResult <- TestXGene(DataDir, ResultDir, finput, sumstat, gene_file,
#'   gene_range, score_file, ref_data, max_gene, sample_size,
#'   genebasedTest = "SKAT",
#'   gene_approximation, beta_par, weights_function, geno_variance_weights,
#'   kernel_p_method, acc_devies, lim_devies, rho, skato_p_threshold, anno_type,
#'   mac_threshold, reference_matrix_used, regularize_fun, pca_var_fraction,
#'   flm_basis_function, flm_num_basis, flm_poly_order, flip_genotypes,
#'   omit_linear_variant
#' )
#' }

TestXGene <- function(DataDir,
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
                      omit_linear_variant = FALSE) {
  tryCatch(
    withCallingHandlers(
      {
        if (!checkFiles(DataDir, finput)) {
          stop("Missing required Plink files in the specified DataDir.")
        }

        # setupPlink(ResultDir)

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

      genes <- read.table(paste0(DataDir, "/", gene_file))
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
        snp_data <- SNPfile %>% select(.data$chr, .data$start, .data$end, .data$SNP)
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
              ), " SNPs")
            )
          )
        )
        gene_snp <- unique(gene_snp_intersect[, c(6, 11)])
        snpcount <- as.data.frame(table(gene_snp$gene_name))

        g <- as.character(snpcount[, 1])
        dir.create(path = file.path(ResultDir, "cormatrix"))
        rlang::inform(rlang::format_error_bullets("SNP-SNP correlation matrices are being created..."))

        snpcorrFun <- function(g) {
          # g <- gene_snp$gene_name
          snps <- gene_snp[gene_snp$gene_name == g, 2, drop = FALSE]
          # print(g)
          write.table(
            snps,
            file = paste0(ResultDir, "/cor_snps.txt"),
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            eol = "\r\n"
          )

          invisible(sys::exec_wait(
            plink(),
            args = c(
              "--bfile",
              paste0(DataDir, "/", finput),
              "--r2",
              "square",
              "--extract",
              paste0(ResultDir, "/cor_snps.txt"),
              "--out",
              paste0(ResultDir, "/snpcorr"),
              "--silent"
            ),
            std_out = FALSE,
            std_err = FALSE
          ))

          snpcorr <- as.matrix(read.table(file = paste0(ResultDir, "/snpcorr.ld")))
          colnames(snpcorr) <- snps$SNP
          rownames(snpcorr) <- snps$SNP
          snpcorr <- as.data.frame(snpcorr)
          save(snpcorr, file = paste0(ResultDir, "/cormatrix/", g, ".RData"))
          return()
        }

        invisible(lapply(g, snpcorrFun))
        rlang::inform(rlang::format_error_bullets(c('v' ="SNP-SNP correlation matrices are done.")))

        score_file <- paste0(ResultDir, "/gene.test.score.file.vcf.gz")
        gene.file <- gene_file

        if (is.null(max_gene)) {
          genes1 <- as.vector(g)
        } else {
          maxgene <- max_gene
          genes1 <- as.vector(g[1:max_gene])
        }

        if (genebasedTest == "SKAT") {
          return(runSKAT(
            score.file = score_file, 
            gene.file = paste0(DataDir, "/", gene_file), 
            genes = genes1, 
            cor.path = paste0(ResultDir, "/cormatrix"), 
            gene_approximation = gene_approximation, 
            anno.type = anno_type, 
            beta.par = beta_par, 
            weights.function = weights_function, 
            geno_variance_weights = geno_variance_weights, 
            kernel_p_method = kernel_p_method, 
            acc_devies = acc_devies, 
            lim_devies = lim_devies, 
            rho = rho, 
            skato_p_threshold = skato_p_threshold))
        } else if (genebasedTest == "SKATO") {
          return(runSKATO(score_file, paste0(DataDir, "/", gene_file), genes1, paste0(ResultDir, "/cormatrix"), anno_type, gene_approximation, beta_par, weights_function, kernel_p_method, acc_devies, lim_devies, rho, skato_p_threshold))
        } else if (genebasedTest == "sumchi") {
          return(runSumChi(score_file, paste0(DataDir, "/", gene_file), genes1, paste0(ResultDir, "/cormatrix"), gene_approximation, anno_type, kernel_p_method, acc_devies, lim_devies))
        } else if (genebasedTest == "ACAT") {
          return(runACAT(score_file, paste0(DataDir, "/", gene_file), genes1, anno_type, beta_par, weights_function, geno_variance_weights, mac_threshold, sample_size))
        } else if (genebasedTest == "BT") {
          return(runBT(score_file, paste0(DataDir, "/", gene_file), genes1, paste0(ResultDir, "/cormatrix"), anno_type, beta_par, weights_function))
        } else if (genebasedTest == "PCA") {
          return(runPCA(score_file, paste0(DataDir, "/", gene_file), genes1, paste0(ResultDir, "/cormatrix"), gene_approximation, anno_type, sample_size, beta_par, weights_function, reference_matrix_used, regularize_fun, pca_var_fraction))
        } else if (genebasedTest == "FLM") {
          return(runFLM(score_file, paste0(DataDir, "/", gene_file), genes1, paste0(ResultDir, "/cormatrix"), gene_approximation, anno_type, sample_size, beta_par, weights_function, flm_basis_function, flm_num_basis, flm_poly_order, flip_genotypes, omit_linear_variant, reference_matrix_used, regularize_fun))
        } else if (genebasedTest == "simpleM") {
          return(runSimpleM(score_file, gene_file, genes1, anno_type, pca_var_fraction))
        } else if (genebasedTest == "minp") {
          return(runMinP(score_file, gene_file, genes1, anno_type))
        }
      },
      error = function(e) {
        message("An error occurred: ", e$message)
        return(NULL)
      },
      warning = function(w) {
        rlang::inform(
          rlang::format_error_bullets(
            c('!' = paste("Warning:", conditionMessage(w)))
          )
        )
        invokeRestart("muffleWarning") 
      }
    )
  )
}


#' SexDiff: Sex difference in effect size for each SNP using t-test.
#'
#' @author Banabithi Bose
#'
#' @description
#' This function uses the GWAS summary statistics from sex-stratified tests like FM01comb or FM02comb, to evaluate
#' the difference in effect size between males and females at each SNP using a t-test.
#'
#' The input dataframes should only include X-chromosome in order to obtain results for sex differences based solely
#' on X-linked loci.
#'
#' @param Mfile
#' R dataframe of summary statistics of GWAS or XWAS of male samples with six mandatory columns, SNP(Variant),CHR(Chromosome number),
#' BP(Base pair position),A1(Minor allele),BETA_M(Effect size) and SE_M(Standard error). This can be generated by running FM01comb or
#' FM02comb model with GXWAS function.
#'
#' @param Ffile
#' R dataframe of summary statistics of GWAS or XWAS of male samples with six mandatory columns, SNP(Variant),CHR(Chromosome number),
#' BP(Base pair position),A1(Minor allele),BETA_F(Effect size) and SE_F(Standard error). This can be generated by running FM01comb or
#' FM02comb model with GXWAS function.
#'
#' @return
#' R dataframe with seven columns:
#'
#' * `SNP` (Variant)
#' * `CHR` (Chromosome number)
#' * `BP` (Base pair position)
#' * `A1` (Minor allele)
#' * `tstat` (t-statistics for effect-size test)
#' * `P` (p-value) and
#' * `adjP` (Bonferroni corrected p-value)
#'
#' @importFrom stats cor pt
#'
#' @export
#'
#' @examples
#' data("GXwasRData")
#' Difftest <- SexDiff(Mfile, Ffile)
#' significant_snps <- Difftest[Difftest$adjP < 0.05, ] # 9
SexDiff <- function(Mfile, Ffile) {
  # Validate input data
  if (!validateInputDataSexDiff(Mfile, Ffile)) {
    return(NULL)
  }

  tryCatch(
    {
      MFWAS <- merge(na.omit(Mfile), na.omit(Ffile), by = c("SNP", "CHR", "BP", "A1"))
      gc(reset = TRUE)

      r <-
        stats::cor(MFWAS$BETA_M, MFWAS$BETA_F, method = "spearman", use = "pairwise.complete.obs")

      MFWAS$tstat <-
        (MFWAS$BETA_M - MFWAS$BETA_F) / sqrt(((MFWAS$SE_M)^2) + (MFWAS$SE_F^2) - 2 * r * (MFWAS$SE_M) * MFWAS$SE_F)

      gc(reset = TRUE)
      MFWAS$P <- stats::pt(q = abs(MFWAS$tstat), df = 2 * nrow(MFWAS) - 2, lower.tail = FALSE)
      MFWAS <- data.table::as.data.table(MFWAS)
      MFWAS$adjP <- p.adjust(MFWAS$P, method = "bonferroni", n = nrow(MFWAS))
      gc(reset = TRUE)
      x <- MFWAS[, c("SNP", "CHR", "BP", "A1", "tstat", "P", "adjP")]
      y <- x[order(x$adjP), ]
      # qq plot
      chisq <- qchisq(1 - y$P, 1)
      lamdaGC <- median(chisq) / qchisq(0.5, 1)
      gc(reset = TRUE)
      par(mar = c(1, 1, 1, 1))
      qqman::qq(y$P, main = paste0(("QQ-plots for test of sex-differentiated/n effect size with GIF = "), round(lamdaGC, 3)))
      gc(reset = TRUE)
      return(y)
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


#' QCsnp: Quality control (QC) for SNPs.
#'
#' @author Banabithi Bose
#'
#' @description
#' This function performs quality control on genetic variant data from PLINK binary files. It filters variants based on parameters
#' like minor allele frequency, Hardy-Weinberg equilibrium, call rate, and differential missingness between cases and controls.
#' It also handles monomorphic SNPs and applies linkage disequilibrium-based filtering if specified.
#'
#' @param DataDir
#' A character string for the file path of the input PLINK binary files.
#'
#' @param ResultDir
#' A character string for the file path where all output files will be stored. The default is `tempdir()`.
#'
#' @param finput
#' Character string, specifying the prefix of the input PLINK binary files with both male and female samples. This file needs
#' to be in `DataDir`.
#'
#' @param foutput
#' Character string, specifying the prefix of the output PLINK binary files if the filtering option for the SNPs is chosen.
#' The default is "FALSE".
#'
#' @param casecontrol
#' Boolean value, `TRUE` or `FALSE` indicating if the input plink files has cases-control status or not.
#' The default is `FALSE`.
#'
#' @param hweCase
#' Numeric value between 0 to 1 or `NULL` for removing SNPs which fail Hardy-Weinberg equilibrium for cases.
#' The default is `NULL`.
#'
#' @param hweControl
#' Numeric value between 0 to 1 or `NULL` for removing SNPs which fail Hardy-Weinberg equilibrium for controls.
#' The default is `NULL`.
#'
#' @param hwe
#' Numeric value between 0 to 1 or `NULL` for removing SNPs which fail Hardy-Weinberg equilibrium for entire dataset.
#' The default is `NULL`.
#'
#' @param maf
#' Numeric value between 0 to 1 for removing SNPs with minor allele frequency less than the specified threshold.
#' The default is 0.05.
#'
#' @param geno
#' Numeric value between 0 to 1 for removing SNPs that have less than the specified call rate. The default is 0.05.
#'
#' Users can set this as `NULL` to not apply this filter.
#'
#' @param monomorphicSNPs
#' Boolean value, `TRUE` or `FALSE` for filtering out monomorphic SNP. The default is `TRUE`.
#'
#' @param caldiffmiss
#' Boolean value, `TRUE` or `FALSE`, specifying whether to compute differential missingness between cases and controls
#' for each SNP (threshold is \eqn{0.05/length(unique(No. of. SNPs in the test))}). The default is `TRUE.`
#'
#' @param diffmissFilter
#' Boolean value, `TRUE` or `FALSE`, specifying whether to filter out the SNPs or only flagged them for differential
#' missingness in cases vs contols. The deafailt is `TRUE`.
#'
#' @param dmissX
#' Boolean value, `TRUE` or `FALSE` for computing differential missingness between cases and controls for X chromosome
#' SNPs only. The default is `FALSE`. The diffmissFilter will work for all these SNPs.
#'
#' @param dmissAutoY
#' Boolean value, `TRUE` or `FALSE` for computing differential missingness between cases and controls for SNPs on autosomes
#' and Y chromosome only. The default is `FALSE`.
#'
#' If `dmissX` and `dmissAutoY` are both `FALSE`, then this will be computed genome-wide. The `diffmissFilter` will work
#' for all these SNPs.
#'
#' @param ld_prunning
#' Boolean value, `TRUE` or `FALSE` for applying linkage disequilibrium (LD)-based filtering.
#'
#' @param highLD_regions
#' A dataframe with known high LD regions \insertCite{Anderson2010}{GXwasR} is provided with the package.
#'
#' @param window_size
#' Integer value, specifying a window size in the variant counts for LD-based filtering. The default is 50.
#'
#' @param step_size
#' Integer value, specifying a variant count to shift the window at the end of each step for LD filtering. The default is 5.
#'
#' @param r2_threshold
#' Numeric value between 0 to 1 of pairwise \eqn{r^2} threshold for LD-based filtering. The default is 0.02.
#'
#' @references
#' \insertAllCited{}
#'
#' @return
#' A list of two objects, namely, `MonomorSNPs` and `DiffMissSNPs` containing monomorphic SNPs and SNPs with differential missingness
#' in cases vs controls, respectively. Output plink binary files in the working directory.
#' @export
#'
#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' foutput <- "Test_output"
#' geno <- NULL
#' maf <- 0.05
#' casecontrol <- FALSE
#' hweCase <- NULL
#' hweControl <- NULL
#' hweCase <- NULL
#' monomorphicSNPs <- FALSE
#' caldiffmiss <- FALSE
#' ld_prunning <- FALSE
#' x <- QCsnp(
#'   DataDir = DataDir, ResultDir = ResultDir, finput = finput, foutput = foutput,
#'   geno = geno, maf = maf, hweCase = hweCase, hweControl = hweControl,
#'   ld_prunning = ld_prunning, casecontrol = casecontrol, monomorphicSNPs = monomorphicSNPs,
#'   caldiffmiss = caldiffmiss
#' )
QCsnp <-
  function(DataDir,
           ResultDir = tempdir(),
           finput,
           foutput = "FALSE",
           casecontrol = TRUE,
           hweCase = NULL,
           hweControl = NULL,
           hwe = NULL,
           maf = 0.05,
           geno = 0.1,
           monomorphicSNPs = FALSE,
           caldiffmiss = FALSE,
           diffmissFilter = FALSE,
           dmissX = FALSE,
           dmissAutoY = FALSE,
           highLD_regions = NULL,
           ld_prunning = FALSE,
           window_size = 50,
           step_size = 5,
           r2_threshold = 0.02) {
    if (!validateInputForQCsnp(DataDir, ResultDir, finput, foutput, casecontrol, hweCase, hweControl, hwe, maf, geno, monomorphicSNPs, caldiffmiss, diffmissFilter, dmissX, dmissAutoY, highLD_regions, ld_prunning, window_size, step_size, r2_threshold)) {
      return(NULL)
    }

    if (!checkFiles(DataDir, finput)) {
      stop("There are no Plink files in DataDir. Please specify correct directory path with input Plink files.")
    }

    tryCatch(
      {
        # setupPlink(ResultDir)

        fam <-
          as.data.frame(utils::read.table(file = paste0(DataDir, "/", finput, ".fam")))

        fam$V6 <- as.numeric(as.character(fam$V6))
        fam <- stats::na.omit(fam)
        fam1 <- fam[fam$V5 != 0, ]
        fam2 <- fam1[fam1$V6 != 0, ]
        fam4 <- fam2[fam2$V6 != -9, ]

        plinkFlags <- setPlinkFlags(maf, geno, hwe, hweCase, hweControl)
        MAF <- plinkFlags$MAF
        GENO <- plinkFlags$GENO
        HWE <- plinkFlags$HWE
        HWECase <- plinkFlags$HWECase
        HWECon <- plinkFlags$HWECon

        # Remove ambiguous SNPs
        removedSNPCount <- removeAmbiguousSNPs(DataDir, ResultDir, finput)
        rlang::inform(rlang::format_error_bullets(c("i" = paste0(removedSNPCount, " Ambiguous SNPs (A-T/G-C), indels etc. were removed."))))


        ## This will be done for the entire file irrespective of case-control status. This will create "filtered_temp1".
        applyFiltersWithPlink(ResultDir, DataDir, finput, MAF, maf, GENO, geno, HWE, hwe)

        ## Apply HWE filters and monomprpic snps. This will create "filtered_temp4"
        casecontrol <- applyCaseControlFilters(ResultDir, fam4, casecontrol, HWECase, hweCase, HWECon, hweControl)

        freq <-
          read.table(file = paste0(ResultDir, "/", "filtered_temp4", ".frq"), stringsAsFactors = FALSE, header = TRUE)

        mmSNPs <- freq[which(freq[, 5] == 0), 2, drop = FALSE]

        write.table(
          mmSNPs,
          file = paste0(ResultDir, "/", "monomorphicSNPs"),
          quote = FALSE,
          row.names = FALSE,
          col.names = FALSE,
          eol = "\r\n",
          sep = " "
        )

        ## Handle Monomorphic SNPs
        monomorphicSNPSettings <- handleMonomorphicSNPs(monomorphicSNPs, mmSNPs, ResultDir)
        mmSNP1 <- monomorphicSNPSettings$mmSNP1
        exclude <- monomorphicSNPSettings$exclude
        excludemono <- monomorphicSNPSettings$excludemono


        ## Handle LD Pruning
        ldPruningSettings <- handleLDPruning(ld_prunning, highLD_regions, ResultDir, window_size, step_size, r2_threshold)
        excluderange <- ldPruningSettings$excluderange
        highLD_regions <- ldPruningSettings$highLD_regions
        indep <- ldPruningSettings$indep
        window_size <- ldPruningSettings$window_size
        step_size <- ldPruningSettings$step_size
        r2_threshold <- ldPruningSettings$r2_threshold

        executePlinkWithParams(ResultDir, "filtered_temp4", exclude, excludemono, excluderange, highLD_regions, indep, window_size, step_size, r2_threshold)

        SNPmissCC <- NULL

        # Filter for case-control differential missingness
        SNPmissCC <- handleCaseControlFiltering(ResultDir, casecontrol, dmissX, dmissAutoY, caldiffmiss, SNPmissCC, diffmissFilter, foutput)

        rlang::inform(rlang::format_error_bullets(c("v" = paste0("Output plink files prefixed as ,", foutput, ", with passed SNPs are saved in ResultDir."))))

        # Remove filtered_temp* files
        removeTempFiles(ResultDir, "filtered_temp")

        # Remove specific files if they exist
        removeTempFiles(ResultDir, "SNPdifCallrate")
        removeTempFiles(ResultDir, "monomorphicSNPs")

        # Remove NoAmbiguousSNP* files
        removeTempFiles(ResultDir, "NoAmbiguousSNP")

        # Remove plink file
        removeTempFiles(ResultDir, "plink")

        # Remove other files
        removeTempFiles(ResultDir, "study_SNP")
        removeTempFiles(ResultDir, "sink_file.txt")

        resultbim <- read.table(paste0(ResultDir, "/", foutput, ".bim"))
        inputbim <- read.table(paste0(DataDir, "/", finput, ".bim"))

        rlang::inform(rlang::format_error_bullets(c("i" = paste0("Input file has ", length(unique(inputbim$V2)), " SNPs."))))
        rlang::inform(rlang::format_error_bullets(c("i" = paste0("Output file has ", length(unique(resultbim$V2)), " SNPs after filtering."))))


        return(list(MonomorSNPs = mmSNP1, DiffMissSNPs = SNPmissCC))
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
#' * `rsid` (i.e., SNP idenitifier)
#' * `beta` (i.e., effect-size or logarithm of odds ratio)
#' * `beta_se` (i.e., standard error of beta)
#' * `P` (i.e., p-values)
#' * `n_eff` (i.e., effective sample size)
#'
#' For case-control study, effective sample size should be \eqn{4 / (1/<# of cases> + 1/<# of controls>)}. The default is `NULL`.
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
#' Integer value, specifying the number of the partision of the GRM in GREML model. The default is 3.
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
#' String of character value specifying the name of the jpeg file with the plots. The deafult is "Heritability_Plots".
#'
#' @returns
#' A dataframe with maximum eight columns for GREML (here, three columns if running genome-wide) and ten columns for LDSC model if byCHR is `TRUE`.
#' The columns, such as, "chromosome"(i.e., chromosome code),"snp_proportion" (i.e.,chromosome-wise SNP propotion)", "no.of.genes" (i.e., number of genes per chromosome),
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
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' precomputedLD <- NULL
#' finput <- "GXwasR_example"
#' data("GXwasRData")
#' test.sumstats <- na.omit(Summary_Stat_Ex1[Summary_Stat_Ex1$TEST == "ADD", c(1:4, 6:8)])
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
#' H2 <- EstimateHerit(
#'   DataDir = DataDir, ResultDir = ResultDir, finput = finput,
#'   summarystat = NULL, ncores, model = "GREML", byCHR = TRUE, r2_LD = 0,
#'   LDSC_blocks = 20, REMLalgo = 0, nitr = 100, cat_covarfile = NULL, quant_covarfile = NULL,
#'   prevalence = 0.01, partGRM = FALSE, autosome = TRUE, Xsome = TRUE, nGRM = 3,
#'   cripticut = 0.025, minMAF = NULL, maxMAF = NULL, hg = "hg19", PlotIndepSNP = TRUE,
#'   IndepSNP_window_size = 50, IndepSNP_step_size = 5, IndepSNP_r2_threshold = 0.02,
#'   highLD_regions = highLD_hg19
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
        # setupPlink(ResultDir)

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

        patterns <- c("plink", "test_reml", "LDfiltered", "HumanGenome", "LDsnp", "test", "gcta", "MAF", "GRM", "phenofile.phen", "sink")

        files_to_remove <- unlist(lapply(patterns, function(pattern) {
          list.files(ResultDir, pattern = patterns, full.names = TRUE)
        }))

        # Use removeFiles helper function to delete the files
        suppressWarnings(removeFiles1(files_to_remove))

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



#' ComputeGeneticPC: Computing principal components from genetic relationship matrix
#'
#' @author Banabithi Bose
#' @description
#' This function performs principal components analysis (PCA) based on the variance-standardized relationship
#' matrix \insertCite{Purcell2007}{GXwasR}.
#'
#' Top principal components are generally used as covariates in association analysis regressions to help correct for
#' population stratification
#'
#'
#' @param DataDir
#' A character string for the file path of the input PLINK binary files.
#'
#' @param finput
#' Character string, specifying the prefix of the input PLINK binary files. This file needs to be in `DataDir.`
#'
#' @param ResultDir
#' A character string for the file path where all output files will be stored. The default is `tempdir()`.
#'
#' @param countPC
#' Integer value, specifying the number of principal components. The default is 10.
#'
#' @param plotPC
#' Boolean value, `TRUE` or `FALSE`, specifying whether to plot the first two PCs.
#'
#' @param highLD_regions
#' A R dataframe with genomic regions with high LD for using in finding pruned SNPs in the plots. The default is `NULL`.
#'
#' @param ld_prunning
#' Numeric value between 0 to 1 of pairwise \eqn{r^2} threshold for LD-based filtering for pruned SNPs in the plots.
#' The default is 0.02.
#'
#' @param window_size
#' Integer value, specifying a window size in variant count or kilobase for LD-based filtering. The default is 50.
#'
#' @param step_size
#' Integer value, specifying a variant count to shift the window at the end of each step for LD filtering for pruned SNPs
#' in the plots. The default is 5.
#'
#' @param r2_threshold
#' Numeric value between 0 to 1 of pairwise \eqn{r^2} threshold for LD-based filtering for pruned SNPs in the plots.
#' The default is 0.02.
#'
#' @return A dataframe with genetic principal components. The first two columns are IID (i.e., Individual Id) and
#' FID (i.e., Family ID). The other columns are PCs.
#'
#' @references
#' \insertAllCited{}
#'
#' @importFrom ggplot2 ggplot geom_bar ylab xlab theme_classic geom_point theme_light coord_equal aes_string
#' @importFrom ggpubr ggarrange
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' data("GXwasRData")
#' highLD_regions <- highLD_hg19
#' ld_prunning <- "TRUE"
#' window_size <- 50
#' step_size <- 5
#' r2_threshold <- 0.02
#' countPC <- 20
#' ## Genetic PC
#' GP <- ComputeGeneticPC(
#'   DataDir = DataDir, ResultDir = ResultDir,
#'   finput = finput, highLD_regions = highLD_hg19, countPC = 20
#' )
ComputeGeneticPC <- function(DataDir, ResultDir = tempdir(), finput, countPC = 10, plotPC = TRUE,
                             highLD_regions = NULL, ld_prunning = TRUE,
                             window_size = 50, step_size = 5, r2_threshold = 0.02) {
  # Validate inputs
  if (!validateInputForComputeGeneticPC(DataDir, ResultDir, finput, countPC, plotPC, highLD_regions, ld_prunning, window_size, step_size, r2_threshold)) {
    return(NULL)
  }

  if (!checkFiles(DataDir, finput)) {
    stop("Missing required Plink files in the specified DataDir.")
  }

  tryCatch(
    {
      # setupPlink(ResultDir)
      

      ## Handle high LD regions exclusion
      processed_file <- paste0(DataDir, "/", finput)
      if (!is.null(highLD_regions)) {
        # Write high LD regions to a temporary file
        options(scipen = 100)
        write.table(highLD_regions,
          file = paste0(ResultDir, "/", "highLD_regions_temp"),
          quote = FALSE, row.names = FALSE, col.names = FALSE
        )
        highLD_regions_file <- paste0(ResultDir, "/", "highLD_regions_temp")
        options(scipen = 0)

        # Exclude high LD regions
        invisible(sys::exec_wait(
          plink(),
          args = c(
            "--bfile",
            processed_file,
            "--exclude", "range",
            highLD_regions_file,
            "--make-bed",
            "--out",
            paste0(ResultDir, "/", "no_highLD_", finput),
            "--silent"
          ),
          std_out = FALSE,
          std_err = FALSE
        ))

        # Update the processed file to the one without high LD regions
        processed_file <- paste0(ResultDir, "/", "no_highLD_", finput)
      }

      ## LD pruning if enabled
      if (ld_prunning == TRUE) {
        # Perform LD pruning
        invisible(sys::exec_wait(
          plink(),
          args = c(
            "--bfile",
            processed_file,
            "--indep-pairwise",
            window_size,
            step_size,
            r2_threshold,
            "--allow-no-sex",
            "--out",
            paste0(ResultDir, "/", "pruned_", finput),
            "--silent"
          ),
          std_out = FALSE,
          std_err = FALSE
        ))

        # Extract pruned SNPs
        invisible(sys::exec_wait(
          plink(),
          args = c(
            "--bfile",
            processed_file,
            "--extract",
            paste0(ResultDir, "/", "pruned_", finput, ".prune.in"),
            "--make-bed",
            "--out",
            paste0(ResultDir, "/", "final_pruned_", finput),
            "--silent"
          ),
          std_out = FALSE,
          std_err = FALSE
        ))

        # Update the processed file to the pruned dataset
        processed_file <- paste0(ResultDir, "/", "final_pruned_", finput)
      }

      ## PCA calculation
      invisible(sys::exec_wait(
        plink(),
        args = c(
          "--bfile",
          processed_file,
          "--pca", countPC,
          "--out",
          paste0(ResultDir, "/", "pcfile"),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))

      # Process PCA results
      PCs1 <- read.table(paste0(ResultDir, "/pcfile.eigenvec"))
      PCs <- PCs1[, -c(1:2)]
      names(PCs) <- paste0("PC", 1:ncol(PCs))
      EV <- scan(paste0(ResultDir, "/pcfile.eigenval"))
      Percent.var <- data.frame(PC = 1:ncol(PCs), Percent.var = EV / sum(EV) * 100)

      # Plot PCA results
      if (plotPC) {
        p1 <- ggplot2::ggplot(data = Percent.var, ggplot2::aes(x = .data$PC, y = .data$Percent.var)) +
          ggplot2::geom_bar(stat = "identity") +
          ggplot2::ylab("Percent Variance Explained") +
          ggplot2::theme_classic()

        p2 <- ggplot2::ggplot(data = PCs, ggplot2::aes(x = .data$PC1, y = .data$PC2)) +
          ggplot2::geom_point() +
          ggplot2::xlab(paste0("PC1 (", signif(Percent.var$Percent.var[1]), "%)")) +
          ggplot2::ylab(paste0("PC2 (", signif(Percent.var$Percent.var[2]), "%)")) +
          ggplot2::theme_light() +
          ggplot2::coord_equal()

        print(ggpubr::ggarrange(p1, p2, labels = c("A", "B"), ncol = 2, nrow = 1))
      }

      # Clean up temporary files
      if (file.exists(paste0(ResultDir, "/pcfile.log"))) {
        file.remove(list.files(ResultDir, pattern = "pcfile", full.names = TRUE))
      }
      file.remove(list.files(ResultDir, pattern = "temp", full.names = TRUE))

      return(PCs1)
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

#' ComputePRS: Computing polygenic risk score (PRS)
#'
#' @author Banabithi Bose
#'
#' @description This function calculates the polygenic risk score, which is the total of allele counts (genotypes) weighted by estimated
#' effect sizes from genome-wide association studies. It uses C+T filtering techniques. The users could perform clumping procedure
#' choromosome-wise and genome-wide. Also, the function offers the choice of including several genetic principal components along with
#' other covariates. Using this function, users have the freedom to experiment with various clumping and thresholding arrangements to
#' test a wide range of various parameter values.
#'
#'
#' @param DataDir
#' A character string for the file path of the all the input files.
#'
#' @param ResultDir
#' A character string for the file path where all output files will be stored. The default is tempdir().
#'
#' @param finput
#' Character string, specifying the prefix of the input PLINK binary files for the genotype data i.e., the target data based on which
#' clumping procedure will be performed. This file needs to be in DataDir. If your target data are small (e.g. N < 500) then you can use
#' the 1000 Genomes Project samples. Make sure to use the population that most closely reflects represents the base sample.
#'
#' @param summarystat
#' A dataframe object with GWAS summary statistics.
#'
#' The mandatory column headers in this dataframe are:
#' * `CHR`(Chromosome code)
#' * `BP`(Basepair position)
#' * `A1` (effect allele)
#' * `SNP` (i.e., SNP idenitifier)
#' * `BETA` or `OR` (i.e., effect-size or logarithm of odds ratio)
#' * `P` (i.e., p-values).
#'
#' Special Notes: The first three columns needed to be `SNP`, `A1` and `BETA` or `OR`.
#'
#' @param phenofile
#' A character string, specifying the name of the mandatory phenotype file. This is a plain text file with no header line; columns
#' family ID, individual ID and phenotype columns. For binary trait, the phenotypic value should be coded as 0 or 1, then it will be
#' recognized as a case-control study (0 for controls and 1 for cases). Missing value should be represented by "-9" or "NA". The
#' interested phenotype column should be labeled as "Pheno1". This file needs to be in `DataDir`.
#'
#' @param covarfile
#' A character string, specifying the name of the covariate file which is a plain .text file with no header line; columns are family ID,
#' individual ID and the covariates. The default is `NULL`. This file needs to be in `DataDir`.
#'
#' @param pheno_type
#' Boolean value, ‘binary’ or ‘quantitative’, specifying the type of the trait. The default is ‘binary’.
#'
#' @param effectsize
#' Boolean value, `BETA` or `OR`, specifying the type of the GWAS effectsize. The default is `BETA`.
#'
#' @param ldclump
#' Boolean value, `TRUE` or `FALSE`, specifying whether to perform clumping or not.
#'
#' @param LDreference
#' A character string, specifying the  prefix of the PLINK files of the population reference panel of the same ancestry, and ideally
#' the one that was used for imputing your target dataset. These files should be in `DataDir`.
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
#' Boolean value, 'TRUE' or 'FALSE', specifying chromosome-wise clumping procedure if `ldclump` was set to be `TRUE`. The default is `TRUE`
#'
#' @param pthreshold
#' Numeric vector, containing several p value thresholds to maximize predictive ability of the derived polygenic scores.
#'
#' @param ld_prunning
#' Boolean value, `TRUE` or `FALSE` for LD-based filtering for computing genetic PC as covariates.
#'
#' @param nPC
#' Positive integer value, specifying the number of genetic PCs to be included as predictor in the PRS model fit. The default is 6.
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
#' A list object containing a dataframe and a numeric value. The dataframe,PRS, contains four mandatory columns, such as,
#' IID (i.e., Individual ID), FID (i.e., Family ID), Pheno1 (i.e., the trait for PRS) and Score (i.e., the best PRS).
#' Other columns of covariates could be there. The numeric value, BestP contains the threshold of
#' of the best p-value for the best pRS model fit.
#'
#' Also, the function produces several plots such as p-value thresholds vs PRS model fit and PRS distribution among male and females.
#' For case-control data, it shows PRS distribution among cases and controls and ROC curves as well.
#'
#' @importFrom dplyr distinct
#' @importFrom stats lm predict logLik
#' @importFrom ggplot2 theme_classic theme element_text ggtitle geom_density xlab aes scale_y_continuous geom_bar scale_fill_gradient2
#' @importFrom data.table as.data.table
#' @importFrom ggpubr ggarrange
#' @export
#'
#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' data("GXwasRData")
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
#' nPC <- 6 # We can incorporate PCs into our PRS analysis to account for population stratification.
#' pheno_type <- "binary"
#'
#' PRSresult <- ComputePRS(DataDir, ResultDir, finput, summarystat, phenofile, covarfile,
#'   effectsize = "BETA", LDreference = "GXwasR_example", ldclump = FALSE, clump_p1, clump_p2,
#'   clump_r2, clump_kb, byCHR = TRUE, pthreshold = pthreshold, highLD_regions = highLD_regions,
#'   ld_prunning = TRUE, window_size = 50, step_size = 5, r2_threshold = 0.02, nPC = 6,
#'   pheno_type = "binary"
#' )
#'
#' ## This table shows 10 samples with phenotype, covariates and a PRS column.
#' PRS <- PRSresult$PRS
#' PRS[1:10, ]
#' ## The best threshold
#' BestPvalue <- PRSresult$BestP$Threshold
#' BestPvalue
ComputePRS <- function(DataDir, ResultDir = tempdir(), finput, summarystat, phenofile, covarfile = NULL,
                       effectsize = c("BETA", "OR"), ldclump = FALSE, LDreference, clump_p1, clump_p2, clump_r2, clump_kb, byCHR = TRUE,
                       pthreshold = c(0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5), highLD_regions, ld_prunning = FALSE,
                       window_size = 50, step_size = 5, r2_threshold = 0.02, nPC = 6, pheno_type = "binary") {
  # Validate inputs
  if (!validateInputForComputePRS(DataDir, ResultDir, finput, summarystat, phenofile, covarfile, effectsize, ldclump, LDreference, clump_p1, clump_p2, clump_r2, clump_kb, byCHR, pthreshold, highLD_regions, ld_prunning, window_size, step_size, r2_threshold, nPC, pheno_type)) {
    return(NULL)
  }

  if (!checkFiles(DataDir, finput)) {
    stop("Missing required Plink files in the specified DataDir.")
  }

  tryCatch(
    {
      # setupPlink(ResultDir)
      

      if (effectsize == "OR") {
        summarystat$OR <- log(summarystat$OR)
      }

      summarystat <- as.data.frame(dplyr::distinct(summarystat, summarystat$SNP, .keep_all = TRUE))

      write.table(summarystat, file = paste0(ResultDir, "/", "prssummarystat"), quote = FALSE, row.names = FALSE)
      SNP.pvalue <- unique(summarystat[, c("SNP", "P")])
      write.table(SNP.pvalue, file = paste0(ResultDir, "/", "SNP.pvalue"), quote = FALSE, row.names = FALSE)


      # LD Clumping
      clumpResults <- performLDClumping(ldclump, DataDir, ResultDir, LDreference, summarystat, clump_p1, clump_p2, clump_r2, clump_kb, byCHR)
      clumpExtract <- clumpResults$clumpExtract
      clumpSNP <- clumpResults$clumpSNP

      # Prepare Phenotype Data
      GP <- preparePhenotypeData(phenofile, nPC, DataDir, ResultDir, finput, highLD_regions, ld_prunning, window_size, step_size, r2_threshold)

      # Read in the phenotype file
      phenotype <- cbind(phenofile[, 1:2], phenofile[, "Pheno1"])
      colnames(phenotype) <- c("FID", "IID", "Pheno1")

      # Read the covariates (here, it is sex)
      if (is.null(covarfile)) {
        pheno <- merge(phenotype, GP, by = c("FID", "IID"))
      } else {
        covariate <- covarfile
        # Now merge the files
        pheno <- merge(merge(phenotype, covariate, by = c("FID", "IID")), GP, by = c("FID", "IID"))
      }

      # We can then calculate the null model (model with PRS) using a linear regression
      # (as height is quantitative) ## Check for binary
      # Compute Null Model
      null_model_result <- computeNullModel(pheno, pheno_type)

      # Accessing the null model and R-squared value
      null_model <- null_model_result$model
      # null_r2 <- null_model_result$r_squared

      ## PRS using thresholding
      ## Best fit PRS
      prsResult <- data.table::rbindlist(lapply(pthreshold, function(pt) {
        prsFun(pt, ResultDir, DataDir, finput, clumpExtract, clumpSNP, pheno, pheno_type, null_model)
      }))


      # Generate a pretty format for p-value output
      prsResult$WriteP <- round(prsResult$P, digits = 3)
      prsResult$WriteP[!is.na(prsResult$WriteP) & prsResult$WriteP == 0] <- format(prsResult$P[!is.na(prsResult$WriteP) & prsResult$WriteP == 0], digits = 2)
      prsResult$WriteP <- sub("e", "*x*10^", prsResult$WriteP)

      p1 <- createPRSPlot(prsResult)

      # Best result is:
      bestP <- prsResult[which.max(prsResult$R2), "Threshold"]
      # Getting PRS score with best p-value threshold
      pt <- cbind(bestP, 0, bestP)
      write.table(pt, file = paste0(ResultDir, "/range_list"), quote = FALSE, row.names = FALSE)
      # By default, if a genotype in the score is missing for a particular individual, then the expected value is imputed, i.e. based on the sample allele frequency. To change this behavior, add the flag --score-no-mean-imputation
      invisible(sys::exec_wait(
        plink(),
        args = c(
          "--bfile", paste0(DataDir, "/", finput),
          "--score", paste0(ResultDir, "/", "prssummarystat"), 1, 2, 3, "header",
          "--q-score-range", paste0(ResultDir, "/range_list"), paste0(ResultDir, "/", "SNP.pvalue"),
          clumpExtract, clumpSNP,
          "--out",
          paste0(ResultDir, "/", "PRS"),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))

      prs <- read.table(paste0(ResultDir, "/", "PRS.", bestP, ".profile"), header = TRUE)
      pheno.prs <- merge(pheno, prs[, c("FID", "IID", "SCORE")], by = c("FID", "IID"))

      ## PRS with sex
      d1 <- pheno.prs[, c("FID", "IID", "Pheno1"), drop = FALSE]
      famfile <- read.table(paste0(DataDir, "/", finput, ".fam"), header = FALSE)
      sex <- famfile[!famfile$V5 == 0, c(1, 2, 5)]
      colnames(sex) <- c("FID", "IID", "SEX")
      dat <- merge(d1, sex, by = c("FID", "IID"))

      # Rename the sex
      dat$SEX[dat$SEX == 1] <- "Male"
      dat$SEX[dat$SEX == 2] <- "Female"
      dat$SEX <- as.factor(as.character(dat$SEX))

      # Merge the files
      dat <- merge(dat, prs, by = c("FID", "IID"))

      # Basic density plot with custom color
      p2 <- createSexDistributionPlot(dat)

      if (pheno_type == "binary") {
        rlang::inform(rlang::format_error_bullets(c("i" = "Plots are initiated.")))
        print(createBinaryPhenotypePlots(dat, p1, p2))
        rlang::inform(rlang::format_error_bullets(c("v" = "Plots are printed.")))
      } else {
        rlang::inform(rlang::format_error_bullets(c("i" = "Plots are initiated.")))
        print(ggpubr::ggarrange(p1, p2))
        rlang::inform(rlang::format_error_bullets(c("v" = "Plots are printed.")))
      }

      # Define patterns for files to be removed
      filePatternsToRemove <- c("pruned_", "PRS", "Clump", "pcfile")

      # Remove files for each pattern
      for (pattern in filePatternsToRemove) {
        ftemp <- list.files(ResultDir, pattern = pattern)
        removeFiles(ftemp, ResultDir)
      }

      # Additional specific files to remove
      additionalFilesToRemove <- c("range_list", "Valid.SNP", "SNPdata_1", "SNP.pvalue", "prssummarystat")
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

  return(list(PRS = pheno.prs, BestP = bestP))
}


#' MergeRegion: Merging two sets of plink binary files.
#'
#' @author Banabithi Bose
#'
#' @description This function combines the two genotype datasets based on either common SNPs or all the SNPs between them.
#'
#' @param DataDir A character string for the file path of the input plink binary files.
#' @param ResultDir A character string for the file path where all output files will be stored. The default is tempdir().
#' @param finput1 Character string, specifying the prefix of the first input PLINK binary files.
#' @param finput2 Character string, specifying the prefix of the first input PLINK binary files.
#' @param foutput Character string, specifying the prefix of the output PLINK binary files if filtering option for the SNPs is chosen. The default is "FALSE".
#' @param use_common_snps Boolean value, TRUE or FALSE, specifying to use common SNPs for merging or to use all the SNPs.
#'
#' @return NULL
#'
#' The output plink files will be saved in ResultDir.
#'
#' @export
#'
#' @examples
#'
#' # Not Run
#' # DataDir <- system.file("extdata", package = "GXwasR")
#' # ResultDir <- tempdir()
#' # finput1 <- "GXwasR_example"
#' # finput2 <- "GXwasR_example_imputed"
#' # foutput <- "Test_output"
#' # Not Run
#' # y <- MergeRegion(DataDir, ResultDir, finput1, finput2, foutput,  use_common_snps = TRUE)
MergeRegion <- function(DataDir, ResultDir, finput1, finput2, foutput, use_common_snps = TRUE) {
  # Validate inputs
  if (!validateInputForMergeRegion(DataDir, ResultDir, finput1, finput2, foutput, use_common_snps)) {
    return(NULL)
  }

  if (!(checkFiles(DataDir, finput1) && checkFiles(DataDir, finput2))) {
    stop("Missing required Plink files in the specified DataDir.")
  }

  tryCatch(
    {
      # setupPlink(ResultDir)
      


      if (use_common_snps) {
        bim1 <- read.table(file.path(DataDir, paste0(finput1, ".bim")))
        bim2 <- read.table(file.path(DataDir, paste0(finput2, ".bim")))
        common_snps <- intersect(bim1$V2, bim2$V2)
        write.table(common_snps,
          file = file.path(ResultDir, paste0("common_snps_", foutput)),
          quote = FALSE, row.names = FALSE, col.names = FALSE
        )

        args1 <- c(
          "--bfile", file.path(DataDir, finput1), "--extract", file.path(ResultDir, paste0("common_snps_", foutput)), "--allow-no-sex",
          "--make-bed", "--out", file.path(ResultDir, paste0("new", finput1)), "--silent"
        )
        executePlink(args1)

        args2 <- c(
          "--bfile", file.path(DataDir, finput2), "--extract", file.path(ResultDir, paste0("common_snps_", foutput)), "--allow-no-sex",
          "--make-bed", "--out", file.path(ResultDir, paste0("new", finput2)), "--silent"
        )
        executePlink(args2)

        merge_args <- c(
          "--bfile", file.path(ResultDir, paste0("new", finput1)), "--bmerge", file.path(ResultDir, paste0("new", finput2)),
          "--allow-no-sex", "--make-bed", "--out", file.path(ResultDir, foutput), "--silent"
        )
        executePlink(merge_args)

        rlang::inform(rlang::format_error_bullets(c("v" = "Merging is done using the common SNPs between the input genotype files.")))

        # Clean-up
        ftemp <- c(list.files(ResultDir, pattern = "new"), list.files(ResultDir, pattern = "common_snps"))
        invisible(file.remove(file.path(ResultDir, ftemp)))
      } else {
        merge_args <- c(
          "--bfile", file.path(DataDir, finput1), "--bmerge", file.path(DataDir, finput2),
          "--allow-no-sex", "--make-bed", "--out", file.path(ResultDir, foutput), "--silent"
        )
        executePlink(merge_args)
        rlang::inform(rlang::format_error_bullets(c("v" = "Merging is done with all the SNPs i.e., union of the SNPs.")))
      }

      rlang::inform(rlang::format_error_bullets(c("v" = paste0("Plink files with merged regions are in ", ResultDir, " prefixed as ", foutput))))
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

#' plinkVCF: Converting VCF files to plink binary files and vice-versa.
#'
#' @author Banabithi Bose
#'
#' @description
#' This function performs the conversion between VCF files to plink binary formats.
#'
#' For VCF to plink files conversion, if you do not specify any FAM file when you are converting from VCF to plink
#' format, then plink will just create a 'dummy' FAM file with the same name as your dataset with missing phenotypes
#' and missing sex.
#'
#' @param DataDir
#' A character string for the file path of the input plink binary files and all other input files.
#'
#' @param ResultDir
#' A character string for the file path where all output files will be stored. The default is `tempdir()`.
#'
#' @param finput
#' Character string, specifying the prefix of the input plink binary files or vcf files. This file needs to be in `DataDir`.
#'
#' @param foutput
#' Character string, specifying the prefix of the output plink binary files if filtering option for the SNPs is chosen.
#' The default is "FALSE".
#'
#' @param VtoP
#' Boolean value, `TRUE` or `FALSE`, specifying the conversion of VCF files to plink binary files or not. The default is `TRUE`.
#'
#' @param PtoV
#' Boolean value, `TRUE` or `FALSE`, specifying the conversion of plink binary files to VCF  files or not. The default is `TRUE`.
#'
#' @param Famfile
#' Character string, specifying the name of the original .fam file if VtoP was set to be `TRUE`. This file needs to be in `DataDir`.
#' The default is `NULL`.
#'
#' @param PVbyCHR
#' Boolean value, `TRUE` or `FALSE` specifying to do the plink to vcf conversion chromosome-wise or not. The default is `TRUE`.
#'
#' @importFrom Rsamtools bgzip indexTabix 
#' 
#' @return
#' `NULL`
#'
#' The output files will be saved in `ResultDir`.
#'
#' @export
#'
#' @examples
#' finput <- "GXwasR_example" # Plink file
#' foutput <- "GXwasR_example1"
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' PtoV <- TRUE
#' VtoP <- FALSE
#' Famfile <- NULL
#' PVbyCHR <- FALSE
#' plinkVCF(DataDir, ResultDir, finput, foutput, VtoP, PtoV, Famfile, PVbyCHR)

plinkVCF <- function(DataDir, ResultDir = tempdir(), finput, foutput, 
                     VtoP = FALSE, PtoV = TRUE, Famfile = NULL, PVbyCHR = TRUE) {
  # Validate Inputs
  if (!validateInputForPlinkVCF(DataDir, ResultDir, finput, foutput, VtoP, PtoV, Famfile, PVbyCHR)) {
    return(NULL)
  }
  
  tryCatch({
    if (PtoV) {
      if (!checkFiles(DataDir, finput)) {
        stop("There are no Plink files in DataDir. Please specify correct DataDir path with input Plink files.")
      }
      
      convertPlinkToVCF <- function(prefix, chr = NULL) {
        args <- c("--bfile", file.path(DataDir, finput),
        "--recode", "vcf", "--allow-extra-chr",
        "--out", file.path(ResultDir, prefix), "--silent")
        if (!is.null(chr)) {
          args <- c(args, "--chr", chr)
        }
        executePlinkAd(ResultDir, args)
        
        # Compress and index
        vcf_path <- file.path(ResultDir, paste0(prefix, ".vcf"))
        vcf_gz_path <- Rsamtools::bgzip(
          file = vcf_path,
          dest = paste0(vcf_path, ".gz"),
          overwrite = TRUE
        )
        Rsamtools::indexTabix(vcf_gz_path, format = "vcf")
      }
      
      if (PVbyCHR) {
        bimfile <- read.table(file.path(DataDir, paste0(finput, ".bim")))
        chrs <- unique(bimfile$V1)
        chrs <- gsub("^chr", "", chrs)  # Normalize chromosome names
        invisible(lapply(chrs, function(chr) {
          message("Processing chromosome: ", chr)
          convertPlinkToVCF(paste0(foutput, "_chr", chr), chr)
        }))
      } else {
        convertPlinkToVCF(paste0(foutput, "_vcf"))
      }
      removeTempFiles(ResultDir, "log")
    }
    
    if (VtoP) {
      vcf_file <- file.path(DataDir, paste0(finput, ".vcf"))
      if (!file.exists(vcf_file)) {
        stop("VCF file not found in DataDir. Please specify correct directory path with input VCF files.")
      }
      
      executePlinkAd(ResultDir, c(
        "--vcf", vcf_file,
        "--keep-allele-order", "--allow-extra-chr",
        "--make-bed", "--const-fid", "1",
        "--out", file.path(ResultDir, foutput), "--silent"
      ))
      
      if (!is.null(Famfile)) {
        fam <- read.table(file.path(DataDir, Famfile))
        write.table(fam, file = file.path(ResultDir, paste0(foutput, ".fam")),
        col.names = FALSE, row.names = FALSE, quote = FALSE)
      } else {
        message("Famfile is NULL. The generated .fam file will have missing phenotypes.")
      }
    }
    
    rlang::inform(
      format_error_bullets(c(
        "v" = paste("Output files created in ResultDir:", ResultDir)
      )
    ))
  },
  error = function(e) {
    message("An error occurred: ", e$message)
    return(NULL)
  },
  warning = function(w) {
    message("Warning: ", w$message)
  })
}


#' SexCheck: Compare sex assignments in the input plink files with those imputed from X chromosome inbreeding coefficients
#'
#' @author Banabithi Bose
#'
#' @description
#' This function compares sex assignments in the input dataset with those predicted from X chromosome inbreeding coefficients \insertCite{Purcell2007}{GXwasR},
#' and gives the option to convert the sex assignments to the predicted values. Implicitly, this function computes observed and expected autosomal homozygous
#' genotype counts for each sample and reports method-of-moments F coefficient estimates (i.e., observed hom. \eqn{count - expected count) / (total observations - expected count)}).
#' The expected counts will be based on loaded or imputed minor allele frequencies.  Since imputed MAFs are highly inaccurate when there are few samples,
#' the 'compute freq' parameter should be set to TRUE to compute MAF implicitly.
#'
#' Due to the use of allele frequencies, if a cohort is comprised of individuals of different ancestries, users may need to process any samples with rare
#' ancestry individually if the dataset has a very unbalanced ancestry distribution. It is advised to run this function with all the parameters set to zero,
#' then examine the distribution of the F estimates (there should be a clear gap between a very tight male clump on the right side of the distribution and the
#' females everywhere else). Then, rerun the function with the parameters that correspond to this gap.
#'
#'
#' @param DataDir
#' Character string for the file path of the input PLINK binary files.
#'
#' @param ResultDir
#' A character string for the file path where all output files will be stored. The default is `tempdir()`.
#'
#' @param finput
#' Character string, specifying the prefix of the input PLINK binary files. Note: Input dataset should contain X and Y regions.
#'
#' @param impute_sex
#' Boolean value, `TRUE` or `FALSE`, specifying sex to be imputed or not. If `TRUE` then sex-imputed PLINK files, prefixed, 'seximputed_plink', will
#' be produced in `DataDir`.
#'
#' @param compute_freq
#' Boolean value, `TRUE` or `FALSE`, specifying minor allele frequency (MAF). This function requires reasonable MAF estimates, so it is essential
#' to use `compute_freq` = `TRUE` for computing MAF from an input PLINK file if there are very few samples in the input dataset. The default is `FALSE`.
#'
#' @param LD
#' Boolean value, `TRUE` or `FALSE` for applying linkage disequilibrium (LD)-based filtering. The default is `TRUE`.
#'
#' @param LD_window_size
#' Integer value, specifying a window size in variant count for LD-based filtering. The default is 50.
#'
#' @param LD_step_size
#' Integer value, specifying a variant count to shift the window at the end of each step for LD filtering. The default is 5.
#'
#' @param LD_r2_threshold
#' Numeric value between 0 to 1 of pairwise \eqn{r^2} threshold for LD-based filtering. The default is 0.02.
#'
#' @param fmax_F
#' Numeric value between 0 to 1. Samples with F estimates smaller than this value will be labeled as females. The default is 0.2.
#'
#' @param mmin_F
#' Numeric value between 0 to 1. Samples with F estimates larger than this value will be labeled as males. The default is 0.8.
#'
#' @return
#' A dataframe with six columns:
#'
#' * `FID` (Family ID)
#' * `IID` (Individual ID)
#' * `PEDSEX `(Sex as determined in pedigree file (1=male, 2=female))
#' * `SNPSEX` (Sex as determined by X chromosome)
#' * `STATUS` (Displays "PROBLEM" or "OK" for each individual)
#' * `F` (The actual X chromosome inbreeding (homozygosity) estimate)
#'
#' A PROBLEM arises if the two sexes do not match, or if the SNP data or pedigree data are ambiguous with regard to sex.
#'
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' LD <- TRUE
#' LD_window_size <- 50
#' LD_step_size <- 5
#' LD_r2_threshold <- 0.02
#' fmax_F <- 0.2
#' mmin_F <- 0.8
#' impute_sex <- FALSE
#' compute_freq <- FALSE
#'
#' x <- SexCheck(
#'   DataDir = DataDir, ResultDir = ResultDir, finput = finput, impute_sex = impute_sex,
#'   compute_freq = compute_freq, LD_window_size = LD_window_size, LD_step_size = LD_step_size,
#'   LD_r2_threshold = 0.02, fmax_F = 0.2, mmin_F = 0.8
#' )
#'
#' # Checking if there is any wrong sex assignment
#' problematic_sex <- x[x$STATUS != "OK", ]
SexCheck <-
  function(DataDir,
           ResultDir = tempdir(),
           finput,
           impute_sex = FALSE,
           compute_freq = FALSE,
           LD = TRUE,
           LD_window_size = 50,
           LD_step_size = 5,
           LD_r2_threshold = 0.02,
           fmax_F = 0.2,
           mmin_F = 0.8) {
    # Validate inputs
    if (!validateInputForSexCheck(DataDir, ResultDir, finput, impute_sex, compute_freq, LD, LD_window_size, LD_step_size, LD_r2_threshold, fmax_F, mmin_F)) {
      return(NULL)
    }

    # Check if required files exist using checkFiles helper function
    if (!checkFiles(DataDir, finput)) {
      stop("Required PLINK files are missing in the specified DataDir.")
    }

    tryCatch(
      {
        # setupPlink(ResultDir)
        

        # Read BIM file to check for X and Y chromosomes
        bim_file <- read.table(file.path(DataDir, paste0(finput, ".bim")))
        xChr <- nrow(subset(bim_file, bim_file$V1 == 23 | bim_file$V1 == "X"))
        yChr <- nrow(subset(bim_file, bim_file$V1 == 24 | bim_file$V1 == "Y"))

        if (xChr == 0) {
          stop("There are no X chromosomes in the input PLINK files.")
        }

        if (yChr == 0) {
          rlang::inform(rlang::format_error_bullets(c("i" = "There are no Y chromosomes in the input PLINK files. Estimates will be based solely on the X chromosome.")))
        }

        if (impute_sex == FALSE) {
          if (compute_freq == TRUE) {
            freq_file_args <- c(
              "--bfile", paste0(DataDir, "/", finput),
              "--freq",
              "--out", paste0(ResultDir, "/freq_file"),
              "--silent"
            )
            executePlink(freq_file_args)

            if (LD == TRUE) {
              # LD pruning and sex check using executePlink
              ld_check_sex_args <- c(
                "--bfile", paste0(DataDir, "/", finput),
                "--indep-pairwise", LD_window_size, LD_step_size, LD_r2_threshold,
                "--read-freq", paste0(ResultDir, "/freq_file.frq"),
                "--check-sex", fmax_F, mmin_F,
                "--out", paste0(ResultDir, "/cs"),
                "--silent"
              )
              executePlink(ld_check_sex_args)
            } else {
              # Sex check using executePlink
              sex_check_args <- c(
                "--bfile", paste0(DataDir, "/", finput),
                "--read-freq", paste0(ResultDir, "/freq_file.frq"),
                "--check-sex", fmax_F, mmin_F,
                "--out", paste0(ResultDir, "/cs"),
                "--silent"
              )
              executePlink(sex_check_args)
            }
          } else if (compute_freq == FALSE) {
            if (LD == TRUE) {
              # LD pruning and sex check using executePlink
              ld_prune_sex_check_args <- c(
                "--bfile", paste0(DataDir, "/", finput),
                "--indep-pairwise", LD_window_size, LD_step_size, LD_r2_threshold,
                "--check-sex", fmax_F, mmin_F,
                "--out", paste0(ResultDir, "/cs"),
                "--silent"
              )
              executePlink(ld_prune_sex_check_args)
            } else {
              # Execute sex check using executePlink
              sex_check_args <- c(
                "--bfile", paste0(DataDir, "/", finput),
                "--check-sex", fmax_F, mmin_F,
                "--out", paste0(ResultDir, "/cs"),
                "--silent"
              )
              executePlink(sex_check_args)
            }

            check_sex <-
              read.table(
                file = paste0(ResultDir, "/", "cs.sexcheck"),
                stringsAsFactors = FALSE,
                header = TRUE
              )
          }
        } else if (impute_sex == TRUE) {
          if (LD == TRUE) {
            plink_args <- c(
              "--bfile", paste0(DataDir, "/", finput),
              "--indep-pairwise", LD_window_size, LD_step_size, LD_r2_threshold,
              "--make-bed",
              "--out", paste0(ResultDir, "/", "csLD"),
              "--silent"
            )

            executePlink(plink_args)

            plink_args_sex_imputation <- c(
              "--bfile", paste0(ResultDir, "/", "csLD"),
              "--impute-sex",
              "--make-bed",
              "--out", paste0(ResultDir, "/", "seximputed_plink"),
              "--silent"
            )

            executePlink(plink_args_sex_imputation)
          } else {
            plink_args_sex_imputation <- c(
              "--bfile", paste0(DataDir, "/", finput),
              "--impute-sex",
              "--make-bed",
              "--out", paste0(ResultDir, "/", "seximputed_plink"),
              "--silent"
            )

            executePlink(plink_args_sex_imputation)
          }

          check_sex <-
            read.table(
              file = paste0(ResultDir, "/", "seximputed_plink.sexcheck"),
              stringsAsFactors = FALSE,
              header = TRUE
            )
          rlang::inform(rlang::format_error_bullets(c("v" = "The output plink files with imputed sex, prefixed, seximputed_plink, are available in the ResultDir.")))
        }

        # Cleanup temporary files
        removeTempFiles(ResultDir, "cs")

        return(check_sex)
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



#' FilterPlinkSample: Making PLINK files with desired samples.
#'
#' @author Banabithi Bose
#'
#' @description
#' This function prepares PLINK binary files with the desired samples.
#'
#' @param DataDir
#' Character string for the file path of the all input files.
#'
#' @param ResultDir
#' character string for the file path where the output PLINK files will be stored.
#'
#' @param finput
#' Character string, specifying the prefix of the input PLINK binary files.
#'
#' @param foutput
#' Character string, specifying the prefix of the output PLINK binary files.
#'
#' @param filter_sample
#' Character string, specifying the sample type to be retained. The choices are, "cases", "controls", "males" and "females".
#' The default is "cases".
#'
#' @param keep_remove_sample_file
#' Character string, specifying the prefix of a space/tab-delimited text file with no header. For the samples that we want
#' to keep or remove, the family IDs should be in the first column and within-family IDs in the second column. This file
#' needs to be in the `DataDir`. The default is `NULL`.
#'
#' @param keep
#' Boolean value, `TRUE` or `FALSE` for specifying desired samples to keep or remove. The default is `TRUE`.
#'
#' @return
#' `NULL`
#'
#' The output plink files with passed samples will be saved in ResultDir.
#'
#' @export
#'
#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' foutput <- "casesPlink"
#' filter_sample <- "cases"
#' keep_remove_sample_file <- "samples_example"
#' keep <- FALSE
#'
#' FilterPlinkSample(
#'   DataDir = DataDir, ResultDir = ResultDir,
#'   finput = finput, foutput = foutput, keep_remove_sample_file = keep_remove_sample_file,
#'   keep = keep
#' )
FilterPlinkSample <- function(DataDir, ResultDir,
                              finput,
                              foutput = NULL,
                              filter_sample = "cases",
                              keep_remove_sample_file = NULL,
                              keep = TRUE) {
  # Validate inputs
  if (!validateInputForFilterPlinkSample(DataDir, ResultDir, finput, foutput, filter_sample, keep_remove_sample_file, keep)) {
    return(NULL)
  }

  if (!checkFiles(DataDir, finput)) {
    stop("Missing required Plink files in the specified DataDir.")
  }

  tryCatch(
    {
      # setupPlink(ResultDir)
      

      if (is.null(keep_remove_sample_file)) {
        invisible(sys::exec_wait(
          plink(),
          args = c(
            "--bed",
            paste0(DataDir, "/", finput, ".bed"),
            "--bim",
            paste0(DataDir, "/", finput, ".bim"),
            "--fam",
            paste0(DataDir, "/", finput, ".fam"),
            paste0("--filter-", filter_sample),
            "--allow-no-sex", # 4.0
            "--make-bed",
            "--out",
            paste0(ResultDir, "/", foutput),
            "--silent"
          ),
          std_out = FALSE,
          std_err = FALSE
        ))
      } else {
        if (keep == TRUE) {
          invisible(sys::exec_wait(
            plink(),
            args = c(
              "--bed",
              paste0(DataDir, "/", finput, ".bed"),
              "--bim",
              paste0(DataDir, "/", finput, ".bim"),
              "--fam",
              paste0(DataDir, "/", finput, ".fam"),
              "--keep", paste0(DataDir, "/", keep_remove_sample_file),
              "--allow-no-sex", # 4.0
              "--make-bed",
              "--out",
              paste0(ResultDir, "/", foutput),
              "--silent"
            ),
            std_out = FALSE,
            std_err = FALSE
          ))
        } else if (keep == FALSE) {
          invisible(sys::exec_wait(
            plink(),
            args = c(
              "--bed",
              paste0(DataDir, "/", finput, ".bed"),
              "--bim",
              paste0(DataDir, "/", finput, ".bim"),
              "--fam",
              paste0(DataDir, "/", finput, ".fam"),
              "--remove", paste0(DataDir, "/", keep_remove_sample_file),
              "--allow-no-sex", # 4.0
              "--make-bed",
              "--out",
              paste0(ResultDir, "/", foutput),
              "--silent"
            ),
            std_out = FALSE,
            std_err = FALSE
          ))
        }
      }
      rlang::inform(rlang::format_error_bullets(c("v" = paste0(foutput, " plink files with desired samples are in ", ResultDir))))
      return()
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


#' GetMFPlink: Getting male and female PLINK binary files.
#'
#' @author Banabithi Bose
#'
#' @description
#' This function prepares separate male and female PLINK binary files from combined PLINK files.
#'
#'
#' @param DataDir
#' Character string for the file path of the input PLINK binary files.
#'
#' @param ResultDir
#' Character string for the file path where all output files will be stored. The default is `tempdir()`.
#'
#' @param finput
#' Character string, specifying the prefix of the input PLINK binary files.
#'
#' @param foutput
#' Character string, specifying the prefix of the output PLINK binary files.
#'
#' @param sex
#' Boolean value, 'males' or 'females', specifying output plink binary files with male or female samples.
#'
#' @param xplink
#' Boolean value, `TRUE` or `FALSE`, specifying output plink binary files with only X chromosome or not. Default is `FALSE.`
#'
#' @param autoplink
#' Boolean value, `TRUE` or `FALSE`, specifying output plink binary files with only autosome or not. Default is `FALSE.`
#'
#' @return
#' None
#'
#' @export
#'
#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' foutput <- "Test_output"
#' sex <- "females"
#' x <- GetMFPlink(
#'   DataDir = DataDir, ResultDir = ResultDir,
#'   finput = finput, foutput = foutput, sex = sex,
#'   xplink = FALSE, autoplink = FALSE
#' )
GetMFPlink <- function(DataDir,
                       ResultDir = tempdir(),
                       finput,
                       foutput,
                       sex,
                       xplink = FALSE,
                       autoplink = FALSE) {
  if (!checkFiles(DataDir, finput)) {
    stop("Missing required Plink files in the specified DataDir.")
  }
  # setupPlink(ResultDir)
  
  # Validate inputs
  if (!validateInputForGetMFPlink(DataDir, ResultDir, finput, foutput, sex, xplink, autoplink)) {
    return(NULL)
  }

  tryCatch(
    {
      if (xplink == FALSE && autoplink == FALSE) {
        invisible(sys::exec_wait(
          plink(),
          args = c(
            "--bed",
            paste0(DataDir, "/", finput, ".bed"),
            "--bim",
            paste0(DataDir, "/", finput, ".bim"),
            "--fam",
            paste0(DataDir, "/", finput, ".fam"),
            paste0("--filter-", sex),
            "--make-bed",
            "--out",
            paste0(ResultDir, "/", foutput),
            "--silent"
          ),
          std_out = FALSE,
          std_err = FALSE
        ))
      } else if (xplink == TRUE && autoplink == FALSE) {
        invisible(sys::exec_wait(
          plink(),
          args = c(
            "--bed",
            paste0(DataDir, "/", finput, ".bed"),
            "--bim",
            paste0(DataDir, "/", finput, ".bim"),
            "--fam",
            paste0(DataDir, "/", finput, ".fam"),
            paste0("--filter-", sex),
            "--chr",
            23,
            "--make-bed",
            "--out",
            paste0(ResultDir, "/", foutput),
            "--silent"
          ),
          std_out = FALSE,
          std_err = FALSE
        ))
      } else if (xplink == FALSE && autoplink == TRUE) {
        invisible(sys::exec_wait(
          plink(),
          args = c(
            "--bed",
            paste0(DataDir, "/", finput, ".bed"),
            "--bim",
            paste0(DataDir, "/", finput, ".bim"),
            "--fam",
            paste0(DataDir, "/", finput, ".fam"),
            paste0("--filter-", sex),
            "--not-chr",
            23,
            "--make-bed",
            "--out",
            paste0(ResultDir, "/", foutput),
            "--silent"
          ),
          std_out = FALSE,
          std_err = FALSE
        ))
      }

      rlang::inform(rlang::format_error_bullets(c("v" = paste0("Output plink files, prefixed as ", foutput, ", are in ", ResultDir))))
      return()
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



#' Xhwe: Filter X-chromosome variants for HWE in females.
#'
#' @author Banabithi Bose
#'
#' @description
#' This function is a part of the post-imputation quality control process prior to GWAS. This tests for Hardy-Weinberg
#' Equilibrium (HWE) for X-chromosome variants in females. Males' hemizygous X chromosome prevents testing for HWE on
#' their haploid X calls, and testing for HWE across all samples would have a high failure rate. This function will check
#' for HWE across the X in females (cases and controls combined), following the recommendation in Khramtsova et al., 2023,
#' and can remove these regions from analysis in all samples. The p-value threshold for filtering out SNPs is 0.05/no.of.
#' X-chromosome variants.
#'
#' @param DataDir
#' Character string for the file path of the input PLINK binary files.
#'
#' @param ResultDir
#' Character string for the file path where all output files will be stored. The default is `tempdir()`.
#'
#' @param finput
#' Character string, specifying the prefix of the input PLINK binary files.
#'
#' @param foutput
#' Character string, specifying the prefix of the output PLINK binary files if filtering option for the SNPs is chosen.
#' The default is "FALSE".
#'
#' @param filterSNP
#' Boolean value, `TRUE` or `FALSE` for filtering out the X-chromosome variants i.e., SNPs from the input file or not
#' (i.e., only flagged). The default is `FALSE`.
#'
#' @return A list object containing SNPs. If `filterSNP` = `TRUE`, the output filtered PLINK binary files will be
#' produced inside `DataDir`.
#'
#' @export
#'
#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' foutput <- "Test_output"
#' x <- Xhwe(
#'   DataDir = DataDir, ResultDir = ResultDir,
#'   finput = finput, foutput = foutput, filterSNP = TRUE
#' )
#' x
Xhwe <- function(DataDir, ResultDir = tempdir(), finput, filterSNP = TRUE, foutput) {
  # Validate inputs
  if (!validateInputForXhwe(DataDir, ResultDir, finput, foutput, filterSNP)) {
    return(NULL)
  }

  if (!checkFiles(DataDir, finput)) {
    stop("Missing required Plink files in the specified DataDir.")
  }

  tryCatch(
    {
      # setupPlink(ResultDir)
      

      ## Getting female Plink file
      invisible(GetMFPlink(DataDir = DataDir, ResultDir, finput = finput, foutput = "female", sex = "females", xplink = FALSE, autoplink = FALSE))

      fam <-
        as.data.frame(utils::read.table(file = paste0(ResultDir, "/", "female.fam")))

      # Check for case control status in input plink file
      fam <- na.omit(fam)
      fam$V6 <- as.numeric(as.character(fam$V6))
      fam2 <- fam[fam$V6 != 0, ]
      fam3 <- fam2[fam2$V6 != -9, ]

      # Updated this warning part
      if (length(unique(fam3$V6)) != 2) {
        writeLines(
          "There is not both case-control status for females in input Plink files."
        )
      } else if (length(unique(fam3$V6)) == 2) {
        writeLines(
          "This test is running on a case-control dataset with female samples."
        )
      }
      ######
      invisible(sys::exec_wait(
        plink(),
        args = c(
          "--bfile",
          paste0(ResultDir, "/", "female"),
          "--chr", 23,
          "--hardy",
          "--out",
          paste0(ResultDir, "/", "Xhwe"),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))


      x <-
        as.data.frame(
          read.table(
            file = paste0(ResultDir, "/", "Xhwe.hwe"),
            header = TRUE,
            sep = ""
          )
        )

      ## Bonferroni-corrected pvalue threshold set to 0.05/(number of X chromosome variants)
      p <- 0.05 / length(unique(x$SNP))
      snp <- x[x$P < p, 2, drop = TRUE]
      X_excluded_SNPs <- unique(snp)

      if (length(X_excluded_SNPs) == 0) {
        rlang::inform(
          rlang::format_error_bullets(c(
            "i" = "No SNP to be excluded.",
            "i" = "Input plink files are unchanged. No output plink files are produced."
          )))
        return()
      } else {
        if (length(X_excluded_SNPs) != 0 & filterSNP == "TRUE") {
          utils::write.table(
            X_excluded_SNPs,
            file = paste0(ResultDir, "/", "XhweSNPs"),
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            eol = "\r\n",
            sep = " "
          )


          invisible(sys::exec_wait(
            plink(),
            args = c(
              "--bfile",
              paste0(DataDir, "/", finput),
              "--exclude",
              paste0(ResultDir, "/", "XhweSNPs"),
              "--allow-no-sex", ## 4.0
              "--make-bed",
              "--out",
              paste0(ResultDir, "/", foutput),
              "--silent"
            ),
            std_out = FALSE,
            std_err = FALSE
          ))


          rlang::inform(rlang::format_error_bullets(c("i" = paste0("Failed SNPs are excluded from the output plink files prefixed as ", foutput, " is in ", ResultDir))))

          ftemp <- list.files(paste0(ResultDir, "/"), pattern = "hwe")
          invisible(file.remove(paste0(ResultDir, "/", ftemp)))
          ftemp <- list.files(paste0(ResultDir, "/"), pattern = "female")
          invisible(file.remove(paste0(ResultDir, "/", ftemp)))
          return(X_excluded_SNPs)
        } else {
          rlang::inform(rlang::format_error_bullets(c("i" = "SNPs are flagged.")))
          ftemp <- list.files(paste0(ResultDir, "/"), pattern = "hwe")
          invisible(file.remove(paste0(ResultDir, "/", ftemp)))
          ftemp <- list.files(paste0(ResultDir, "/"), pattern = "female")
          invisible(file.remove(paste0(ResultDir, "/", ftemp)))
          return(X_excluded_SNPs)
        }
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



#' MAFdiffSexControl: Test for significantly different minor allele frequency (MAF) between sexes in control samples
#'
#' @author Banabithi Bose
#'
#' @description
#' With parameters to filter out SNPs and/or flag the SNPs, this function tests for significantly different MAF
#' (p-value < 0.05/no. of SNPs) between sexes in control samples solely for binary phenotypes. Since the disparities
#' may be caused by technical confounding or sample biases for the research cohorts, it is advised that any SNPs in
#' the controls with a sex difference in MAF be carefully evaluated and identified for further examination
#' (Khramtsova et. al., 2023). In autosomal allele frequencies, sex differences are not anticipated.
#'
#' @param DataDir
#' Character string for the file path of the input PLINK binary files.
#'
#' @param ResultDir
#' Character string for the file path where all output files will be stored. The default is `tempdir()`.
#'
#' @param finput
#' Character string, specifying the prefix of the input PLINK binary files with both male and female samples.
#' This file needs to be in `DataDir`.
#'
#' @param foutput
#' Character string, specifying the prefix of the output PLINK binary files if filtering option for the SNPs
#' is chosen. The default is NULL.
#'
#' @param filterSNP
#' Boolean value, `TRUE` or `FALSE` for filtering out the SNPs or not (i.e., only flagged). The default is `FALSE`.
#'
#' @return
#' A list object containing excluded or flagged SNPs. If `filterSNP` = `TRUE`, the output filtered PLINK binary
#' files will be produced inside `DataDir`.
#'
#' @importFrom stats na.omit
#' @importFrom utils download.file read.table unzip write.table
#' @importFrom sys exec_wait
#' @export
#'
#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' foutput <- "Test_output"
#' x <- MAFdiffSexControl(DataDir, ResultDir, finput, filterSNP = TRUE, foutput = foutput)
MAFdiffSexControl <- function(DataDir,
                              ResultDir = tempdir(),
                              finput,
                              filterSNP = FALSE,
                              foutput = NULL) {
  if (!validateInputForMAFdiffSexControl(DataDir, ResultDir, finput, filterSNP, foutput)) {
    return(NULL)
  }

  if (!checkFiles(DataDir, finput)) {
    stop("Missing required Plink files in the specified DataDir.")
  }

  tryCatch(
    {
      # setupPlink(ResultDir)
      

      fam <-
        as.data.frame(utils::read.table(file = paste0(DataDir, "/", finput, ".fam")))

      # Check for sex
      fam$V6 <- as.numeric(as.character(fam$V6))
      fam <- stats::na.omit(fam)
      fam1 <- fam[fam$V5 != 0, ]
      fam2 <- fam1[fam1$V6 != 0, ]
      fam4 <- fam2[fam2$V6 != -9, ]


      if (length(unique(fam1$V5)) == 2 &&
        length(unique(fam4$V6)) == 2) {
        # For having phenotype for control sample only
        fam$V7 <- 0
      } else {
        writeLines(
          "There is incorrect male-female sex status or incorrect case-control status in input Plink files.\nNeeds both male and female samples with both case and control status to run this function."
        )
      }

      fam$V7[fam$V5 == 1 & fam$V6 == 1] <- 1 # for male and control
      fam[fam$V5 == 2 &
        fam$V6 == 1, 7] <- 2 # for female and control
      fam[fam$V5 == 1 & fam$V6 == 2, 7] <- -9
      fam[fam$V5 == 2 & fam$V6 == 2, 7] <- -9

      phenofile <- unique(fam[, c(1, 2, 7)])
      utils::write.table(
        phenofile,
        file = paste0(ResultDir, "/", "phenofile"),
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE,
        eol = "\r\n",
        sep = " "
      )

      invisible(sys::exec_wait(
        plink(),
        args = c(
          "--bfile",
          paste0(DataDir, "/", finput),
          "--logistic",
          "--pheno",
          paste0(ResultDir, "/", "phenofile"),
          "--out",
          paste0(ResultDir, "/", "OUTPUT"),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))

      x <-
        utils::read.table(
          file = paste0(ResultDir, "/", "OUTPUT.assoc.logistic"),
          header = TRUE
        )

      # Filter for SEX chromosome (i.e. Filter for X chromosome)

      y <- x[, c(1, 2, 9)]
      y <- na.omit(y)
      bf <-
        0.05 / length(unique(y$SNP)) # taking the SNPs for which we have finite p-values.
      y$P <- as.numeric(as.character(y$P))

      flaggedSnps <- unique(y[y$P < bf, 2, drop = TRUE])

      if (length(flaggedSnps) == 0) {
        rlang::inform(rlang::format_error_bullets(c("i" = "No SNP to be flagged or excluded.")))
        flaggedSnps <- NULL
      } else if (length(flaggedSnps) != 0 & filterSNP == TRUE) {
        utils::write.table(
          flaggedSnps,
          file = paste0(ResultDir, "/", "flaggedSnpsSexMafDiff"),
          quote = FALSE,
          row.names = FALSE,
          col.names = FALSE,
          eol = "\r\n",
          sep = " "
        )


        invisible(sys::exec_wait(
          plink(),
          args = c(
            "--bfile",
            paste0(DataDir, "/", finput),
            "--exclude",
            paste0(ResultDir, "/", "flaggedSnpsSexMafDiff"),
            "--allow-no-sex", # 4.0
            "--make-bed",
            "--out",
            paste0(ResultDir, "/", foutput),
            "--silent"
          ),
          std_out = FALSE,
          std_err = FALSE
        ))


        writeLines(
          paste0("SNPs with significantly MAF difference are excluded.\nFiltered plink files are saved in ", ResultDir)
        )
        return(as.list(flaggedSnps))
      } else if (length(flaggedSnps) != 0 & filterSNP == FALSE) {
        rlang::inform(rlang::format_error_bullets(c("i" = "SNPs are flagged.")))
      }

      ftemp <- list.files(paste0(ResultDir, "/"), pattern = "OUTPUT")
      invisible(do.call(file.remove, list(paste0(ResultDir, "/", ftemp))))
      invisible(do.call(file.remove, list(paste0(ResultDir, "/", "phenofile"))))

      return(flaggedSnps)
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



#' QCsample: Quality control for samples in the plink binary files.
#'
#' @author Banabithi Bose
#'
#' @description
#' This function identifies outlier individuals for heterozygosity and/or missing genotype rates, which aids in the
#' detection of samples with subpar DNA quality and/or concentration that should be removed from the study. Individuals
#' missing more than 3-7% of their genotype calls are often excluded from the analysis.
#'
#' Having the correct designation of sex is important to obtain accurate genotype rate estimates, or avoid incorrectly
#' removing samples, etc. Details can be accessed from the paper.

#'
#' @param DataDir
#' Character string, specifying the file path of the input PLINK binary files. The default is `NULL`.
#'
#' @param ResultDir
#' A character string for the file path where all output files will be stored. The default is `tempdir()`.
#'
#' @param finput
#' Character string, specifying the prefix of the input PLINK binary files with both male and female samples.
#' This file needs to be in `DataDir`.
#'
#' @param foutput
#' Character string, specifying the prefix of the output PLINK binary files if filtering option for the samples is chosen.
#'
#' @param imiss
#' Numeric value between 0 to 1 for removing samples that have more than the specified missingness. The default is 0.03.
#'
#' @param het
#' Positive numeric value, specifying the standard deviation from the mean heterozygosity rate. The samples whose rates are more
#' than the specified sd from the mean heterozygosity rate are removed. The default is 3. With this default value, outlying
#' heterozygosity rates would remove individuals who are three sd away from the mean rate (1).
#'
#' @param small_sample_mod
#' Boolean value indicating whether to apply modifications for small sample sizes. Default is `FALSE`.
#'
#' @param IBD
#' Numeric value for setting the threshold for Identity by Descent (IBD) analysis. Default is `NULL`.
#'
#' @param IBDmatrix
#' Boolean value indicating whether to generate an entire IBD matrix. Default is `FALSE`. In this case filtered IBD
#' matrix will be stored.
#'
#' @param ambi_out
#' Boolean value indicating whether to process ambiguous samples.
#'
#' @param title_size
#' Integer, specifying the size of the title of the plot heterozygosity estimate vs missingness across samples.
#'
#' @param legend_text_size
#' Integer, specifying the size for legend text in the plot.
#'
#' @param legend_title_size
#' Integer, specifying the size for the legend title in the plot.
#'
#' @param axis_text_size
#' Integer, specifying the size for axis text in the plot.
#'
#' @param axis_title_size
#' Integer, specifying the size for the axis title in the plot.
#'
#' @param filterSample
#' Boolean value, `TRUE` or `FALSE` for filtering out the samples or not (i.e., only flagged). The default is `TRUE`.
#'
#' @importFrom stats sd
#' @importFrom ggplot2 ggplot
#'
#' @return
#' A plot of heterogysity estimate vs missingness accross sample and a list containing five R dataframe objects, namely,
#' `HM` (samples with outlying heterozygosity and/or missing genotype rates), `Failed_Missingness` (samples with missing genotype rates),
#' `Failed_heterozygosity` (samples with outlying heterozygosity), `Missingness_results` (missingness results) and `Heterozygosity_results`
#' (heterozygosity results) with output plink files in ResultDir if filtering out the samples option is chosen.
#'
#' `Missingness_results` contains missingness results for each individual, with six columns as `FID`, `IID`, `MISS_PHENO`, `N_MISS`, `N_GENO` and
#' `F_MISS` for Family ID, Within-family ID, Phenotype missing? (Y/N), Number of missing genotype call(s), not including obligatory missings
#' or heterozygous haploids, number of potentially valid call(s), and missing call rate, respectively.
#'
#' `Heterozygosity_results` contains heterozygosity results for each individual, with six columns as `FID`, `IID`, `O(HOM)`, `E(HOM)`, `N(NM)`,
#' and `F` for Family ID, Within-family ID, Observed number of homozygotes, Expected number of homozygotes, Number of (non-missing, non-monomorphic)
#' autosomal genotype observations and, Method-of-moments F coefficient estimate, respectively.
#' @export
#'
#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' foutput <- "Test_output"
#' imiss <- 0.01
#' het <- 2
#' small_sample_mod <- FALSE
#' IBD <- 0.2
#' IBDmatrix <- FALSE
#' ambi_out <- TRUE
#'
#' x <- QCsample(
#'   DataDir = DataDir, ResultDir = ResultDir, finput = finput,
#'   foutput = foutput, imiss = imiss, het = het, IBD = IBD,
#'   ambi_out = ambi_out
#' )
QCsample <- function(DataDir,
                     ResultDir,
                     finput,
                     foutput = NULL,
                     imiss,
                     het,
                     small_sample_mod = FALSE,
                     IBD,
                     IBDmatrix = FALSE,
                     ambi_out = TRUE,
                     legend_text_size = 8,
                     legend_title_size = 7,
                     axis_text_size = 5,
                     axis_title_size = 7,
                     title_size = 9,
                     filterSample = TRUE) {
  # Validate parameters
  validateInputForQCsample(DataDir, ResultDir, finput, foutput, imiss, het, small_sample_mod, IBD, IBDmatrix, ambi_out, legend_text_size, legend_title_size, axis_text_size, axis_title_size, title_size, filterSample = TRUE)

  if (!checkFiles(DataDir, finput)) {
    stop("There are no Plink files in DataDir. Please specify correct directory path with input Plink files.")
  }

  tryCatch(
    {
      # setupPlink(ResultDir)
      

      if (small_sample_mod == TRUE) {
        SSM <- "small-sample"
      } else {
        SSM <- NULL
      }


      ############ Adding this in version 3.0 #############
      ## Filter out samples with missing phenotype

      fam1 <- nrow(read.table(paste0(DataDir, "/", finput, ".fam"), header = FALSE))

      # finput <- if(ambi_out) processAmbiguousSamples(DataDir, ResultDir, finput, fam1) else finput ## Closing it

      # Prepare the arguments for the Plink command for missing data analysis
      missingDataArgs <- c(
        "--bfile", paste0(DataDir, "/", finput),
        "--missing",
        "--out", paste0(ResultDir, "/", "filtered_missing"),
        "--silent"
      )

      executePlink(missingDataArgs)


      # Prepare the arguments for the Plink command for heterozygosity analysis
      heterozygosityArgs <- c(
        "--bfile", paste0(DataDir, "/", finput),
        "--dog",
        "--het",
        SSM,
        "--out", paste0(ResultDir, "/", "filtered_hetero"),
        "--silent"
      )

      executePlink(heterozygosityArgs)


      miss <- readDataFile(paste0(ResultDir, "/", "filtered_missing.imiss"))
      heter <- readDataFile(paste0(ResultDir, "/", "filtered_hetero.het"))

      heter$F <- as.numeric(as.character(heter$F))
      miss$F_MISS <- as.numeric(as.character(miss$F_MISS))

      imissfail <- miss[miss$F_MISS > imiss, , drop = FALSE]


      hetfail <- heter[heter$F < (mean(heter$F) - het * stats::sd(heter$F)) |
        heter$F > (mean(heter$F) + het * stats::sd(heter$F)), , drop = FALSE]


      hetermiss <- merge(miss, heter, by = "IID")

      failed_het_imiss <-
        hetermiss[which(hetermiss$IID %in% union(hetfail$IID, imissfail$IID)), , drop = FALSE]

      write.table(
        failed_het_imiss[, c(2, 1)],
        file = paste0(ResultDir, "/", "failed_het_imiss"),
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE,
        eol = "\r\n",
        sep = " "
      )

      ## Updating it in 6.0
      if (!is.null(imiss) && !is.null(het)) {
        filterSamples(DataDir, ResultDir, finput, failed_het_imiss, filterSample)


        ## Plot
        print(createHeterozygosityPlot(hetermiss, hetfail, imissfail, het, imiss, legend_text_size, legend_title_size, axis_text_size, axis_title_size, title_size))

        printSampleFilterResults(imissfail, hetfail, failed_het_imiss)


        if (nrow(hetermiss) == 0) {
          hetermiss <- NULL
        } else {
          hetermiss <- hetermiss
        }

        if (nrow(hetermiss) == 0) {
          hetermiss1 <- NULL
        } else {
          hetermiss1 <- hetermiss[, 1:2]
        }

        if (nrow(imissfail) == 0) {
          imissfail1 <- NULL
        } else {
          imissfail1 <- imissfail[, 1:2]
        }

        if (nrow(hetfail) == 0) {
          hetfail1 <- NULL
        } else {
          hetfail1 <- hetfail[, 1:2]
        }

        ftemp <- c("failed_het_imiss", "filtered_hetero.log", "filtered_missing.lmiss", "filtered_missing.log")
        invisible(do.call(file.remove, list(paste0(ResultDir, "/", ftemp))))

        if (file.exists(paste0(ResultDir, "/filtered_missing.hh"))) {
          ftemp <- c("filtered_missing.hh")
          invisible(do.call(file.remove, list(paste0(ResultDir, "/", ftemp))))
        }
        if (file.exists(paste0(ResultDir, "/", "foutput", ".hh"))) {
          ftemp <- c(paste0("foutput", ".log"), paste0("foutput", ".hh"))
          invisible(do.call(file.remove, list(paste0(ResultDir, "/", ftemp))))
        }
        fmi <- read.table(paste0(ResultDir, "/filtered_missing.imiss"))
        fhh <- read.table(paste0(ResultDir, "/filtered_hetero.het"))

        ftemp <- c("filtered_missing.imiss", "filtered_hetero.het")
        invisible(do.call(file.remove, list(paste0(ResultDir, "/", ftemp))))

        hm <- failed_het_imiss[, 2:1]
        colnames(hm) <- c("FID", "IID")
      } else {
        excludeSamplesArgs <- c(
          "--bfile", paste0(DataDir, "/", finput),
          "--make-bed",
          "--out", paste0(ResultDir, "/", "foutput"),
          "--silent"
        )
        executePlink(excludeSamplesArgs)
        rlang::inform(rlang::format_error_bullets(c("i" = "Missingness and heterogygosity thresholds are NULL.")))
        hm <- NULL
        imissfail <- NULL
        hetfail <- NULL
        fmi <- NULL
        fhh <- NULL
      }

      ######## IBD########
      fd <- processIBDData(IBD, IBDmatrix, ResultDir, foutput, filterSample)
      failed_ibd <- fd$failed_ibd
      ibd <- fd$ibd

      if (is.null(IBD)) {
        failed_ibd <- NULL
      } else {
        failed_ibd <- failed_ibd
      }

      if (is.null(IBD)) {
        rlang::inform(rlang::format_error_bullets(c("i" = "There will be no sample to be filtered for IBD with 'IBD' threshold.")))
        ibd <- NULL
      } else if (!is.null(failed_ibd)) {
        if (filterSample == TRUE) {
          rlang::inform(rlang::format_error_bullets(c("i" = paste0("No. of samples marked to be filtered out for IDB after missingness and heterozygosity filter: ", nrow(failed_ibd)))))
        } else if (filterSample == FALSE) {
          rlang::inform(rlang::format_error_bullets(c("i" = paste0("No. of samples are flagged out for IDB after missingness and heterozygosity filter: ", nrow(failed_ibd)))))
        }
      } else if (is.null(failed_ibd)) {
        rlang::inform(rlang::format_error_bullets(c("i" = "No sample is filtered out for IDB after missingness and heterozygosity filter.")))
      }

      ###################

      Outputsample <- nrow(read.table(paste0(ResultDir, "/", foutput, ".fam"), header = FALSE))
      rlang::inform(rlang::format_error_bullets(c("i" = paste0("No. of samples in input plink files: ", fam1))))
      rlang::inform(rlang::format_error_bullets(c("i" = paste0("No. of samples in output plink files: ", Outputsample))))
      rlang::inform(rlang::format_error_bullets(c("v" = paste0("Output plink files, ", foutput, " with final samples are in ", ResultDir, "."))))
      if (filterSample == FALSE) {
        rlang::inform(rlang::format_error_bullets(c("i" = "Samples are flagged only.")))
      }
      return(list(
        HM = hm,
        Failed_Missingness = imissfail[, 1:2],
        Failed_heterozygosity = hetfail[, 1:2],
        Failed_IBD = failed_ibd[, 1:2],
        Missingness_results = fmi,
        Heterozygosity_results = fhh,
        IBD_results = ibd
      ))
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




#' Miami plot
#'
#' @description
#' This function generates Miami plots for GWAS and XWAS.
#'
#' @param ResultDir
#' Character string for the folder path where the outputs will be saved.
#'
#' @param FemaleWAS
#' R dataframe of summary statistics of GWAS or XWAS of female samples with four columns, SNP(Variant),
#' CHR(Chromosome number), POS(Base pair position) and pvalue(P-value of the test). This can be generated
#' by running FM01comb or FM02comb model with GXWAS function.
#'
#' @param MaleWAS
#' R dataframe of summary statistics of GWAS or XWAS of male samples with four columns, SNP(Variant),
#' CHR(Chromosome number), POS(Base pair position) and pvalue(P-value of the test). This can be generated
#' by running FM01comb or FM02comb model with GXWAS function.
#'
#' @param snp_pval
#' Numeric value as p-value threshold for annotation. SNPs below this p-value will be annotated on the plot.
#' The default is 1e-08.
#'
#' @param Xchr
#' Boolean value, `TRUE` or `FALSE`, specifying whether to generate Miami plot for stratified XWAS or not.
#' The default is `TRUE`.
#'
#' @return Invisibly returns `NULL`. Generates and saves Miami plots as a side effect.
#' @export
#'
#' @examples
#' data("GXwasRData")
#' FemaleWAS <- na.omit(Ffile[, c("SNP", "CHR", "BP", "P")])
#' colnames(FemaleWAS) <- c("SNP", "CHR", "POS", "pvalue")
#' MaleWAS <- na.omit(Mfile[, c("SNP", "CHR", "BP", "P")])
#' colnames(MaleWAS) <- c("SNP", "CHR", "POS", "pvalue")
#'
#' GXWASmiami(FemaleWAS = FemaleWAS, MaleWAS = MaleWAS, snp_pval = 0.05)
GXWASmiami <- function(ResultDir = tempdir(), FemaleWAS, MaleWAS, snp_pval = 1e-08, Xchr = FALSE) {
  # Validate input parameters
  validateInputForGXWASmiami(ResultDir, FemaleWAS, MaleWAS, snp_pval, Xchr)

  tryCatch(
    {
      rlang::inform(rlang::format_error_bullets("Generating Miami plots for stratified test."))
      suppressWarnings(invisible(gmirror(
        top = FemaleWAS, bottom = MaleWAS, tline = snp_pval, bline = snp_pval,
        toptitle = "GWAS of females", bottomtitle = "GWAS of males",
        highlight_p = c(snp_pval, snp_pval), highlighter = "green", chrblocks = TRUE, file = paste0(ResultDir, "/", "Stratified_GWAS")
      )))

      rlang::inform(rlang::format_error_bullets(c("v" = paste0("Miami plot of stratified GWAS is saved in ", ResultDir))))

      if (Xchr == TRUE) {
        FemaleWAS <- as.data.frame(FemaleWAS)
        FemaleWAS[FemaleWAS$CHR == "23", "CHR"] <- "X"
        MaleWAS <- as.data.frame(MaleWAS)
        MaleWAS[MaleWAS$CHR == "23", "CHR"] <- "X"

        # Stratified XWAS plot
        gwas.t2 <- FemaleWAS[FemaleWAS$CHR == "X", ]
        gwas.b2 <- MaleWAS[MaleWAS$CHR == "X", ]

        rm(FemaleWAS)
        rm(MaleWAS)
        suppressWarnings(invisible(gmirror(
          top = gwas.t2, bottom = gwas.b2, tline = snp_pval, bline = snp_pval,
          toptitle = "XWAS of females", bottomtitle = "XWAS of males",
          highlight_p = c(snp_pval, snp_pval), highlighter = "green", chrblocks = TRUE, file = paste0(ResultDir, "/", "Stratified_XWAS")
        )))
        gc(reset = TRUE)
        rlang::inform(rlang::format_error_bullets(c("v" = paste0("Miami plot of stratified XWAS is saved in ", ResultDir))))
      } else {
        return(invisible(NULL))
      }
      return(invisible(NULL))
    },
    error = function(e) {
      rlang::abort("An error occurred: ", e$message)
    },
    warning = function(w) {
      rlang::warn("Warning: ", w$message)
    }
  )
}

#' GXwas: Running genome-wide association study (GWAS) and X-chromosome-wide association study (XWAS) models.
#'
#' @author Banabithi Bose
#'
#' @description
#' This function runs GWAS models in autosomes with several alternative XWAS models.
#' Models such as "FMcombx01","FMcombx02",and "FMstatrified" can be applied to both binary and quantitative traits,
#' while "GWAcxci" can only be applied to a binary trait.
#'
#' For binary and quantitative features, this function uses logistic and linear regression,
#' allowing for multiple covariates and the interactions with those covariates in a multiple-regression approach.
#' These models are all run using the additive effects of SNPs, and each additional minor allele's influence
#' is represented by the direction of the regression coefficient.
#'
#' This function attempts to identify the multi-collinearity among predictors by displaying NA for the test statistic
#' and a p-value for all terms in the model. The more terms you add, the more likely you are to run into issues.
#'
#' For details about the different XWAS model, please follow the associated publication.
#'
#' @param DataDir
#' Character string for the file path of the input PLINK binary files.
#'
#' @param ResultDir
#' Character string for the folder path where the outputs will be saved.
#'
#' @param finput
#' Character string, specifying the prefix of the input PLINK binary files with both male and female samples.
#' This file needs to be in `DataDir`.
#'
#' Note: Case/control phenotypes are expected to be encoded as 1=unaffected (control), 2=affected (case); 0 is accepted as
#' an alternate missing value encoding. The missing case/control or quantitative phenotypes are expected to be encoded as 'NA'/'nan'
#' (any capitalization) or -9.
#' @param trait
#' Boolean value, 'binary' or 'quantitative' for the phenotype i.e. the trait.
#'
#' @param sex
#' Boolean value, `TRUE` or `FALSE` for using sex as covariate in association test. It is applicable genome-wide.
#'
#' The default is FALSE.
#' @param xsex
#' Boolean value, `TRUE` or `FALSE` for using sex as covariate in association test for X-chromosomal SNPs.
#' The default is FALSE. This will overwrite 'sex' argument for X-chromosome.
#'
#' @param standard_beta
#' Boolean value, `TRUE` or `FALSE` in case of quantitative trait for standardizing the trait or phenotype values
#' (mean 0, unit variance), so the resulting coefficients will be standardized. The default is `TRUE`.
#'
#' @param xmodel
#' Models "FMcombx01","FMcombx02",and "FMstatrified" can be chosen for both binary and quantitative traits
#' while "GWAcxci" can only apply to the binary trait. These models take care of the X-chromosomal marker.
#' Three female genotypes are coded by 0, 1, and 2 in FM01 and FM02. The two genotypes of males that follow the
#' X-chromosome inactivation (XCI) pattern as random (XCI-R) in the FM01 model are coded by 0 and 1, while the two
#' genotypes that follow the XCI is escaped (XCI-E) in the FM02 model are coded by 0 and 1. To reflect the dose
#' compensation connection between the sexes, FM02 treats men as homozygous females.
#'
#' In the FM01comb and FM01comb methods, associations are tested separately for males and females with the FM01 and FM02
#' models, respectively, and then the combined p values are computed the Fisher's method, Fisher's method with permutation,
#' or Stouffer's method(1,3-7]. An X-chromosome inactivation (XCI) pattern, or coding technique for X-chromosomal genotypes
#' between sexes, is not required for the XCGA. By simultaneously accounting for four distinct XCI patterns, namely XCI-R, XCI-E,
#' XCI-SN (XCI fully toward normal allele), and XCI-SR (XCI fully toward risk allele), this model may maintain a
#' respectably high power \insertCite{Su2022}{GXwasR}.
#'
#' Note: `sex` shouldn't be provided as a covariate in the XCGA model.
#'
#' @param covarfile
#' Character string for the full name of the covariate file in .txt format. This file should be placed in `DataDir`.
#'
#' Note about the covariate file: The first column of this file should be `FID`, the second column should be `IID` and
#' the other columns should be covariates. The primary header line should be there starting with “FID”, and “IID”
#' followed by covariate names. If an individual is not present in the covariate file, or if the individual has a
#' missing phenotype value (i.e. -9 by default) for the covariate, then that individual is set to missing (i.e. will
#' be excluded from association analysis). It is important to note that for statrified GWAS model, if PCs are included
#' as covar then it should be generated separately for each cohort and then included in the covarfile. Use the function
#' \code{\link{DummyCovar}} to generate a new covariate file with categorical variables down-coded as binary dummy variables for
#' the covariate file with categorical variables. For instance, if a variable has K categories, K-1 new dummy variables
#' are constructed, and the original covariate is now estimated with a coefficient for each category.
#'
#' @param covartest
#' Vector value with `NULL`,"ALL" or covarite name/names to be included in the test. The default is `NULL.` For instance,
#' the user can choose “AGE” and “SEX” as covartest = c(“AGE”, “SEX”) or all the covariates as covartest = c(“ALL”).
#'
#' @param interaction
#' Boolean value, `TRUE` or `FALSE` for including SNP x covariate interaction term/terms from the association analysis.
#' The default is `FALSE`. If a permutation procedure is chosen then the interaction will be automatically `FALSE`. For the
#' interaction with the two covariates COV1 and COV2, the model will look like: \eqn{Y = b0 + b1.ADD + b2.COV1 + b3.COV2 +
#' b4.ADD x COV1 + b5.ADD x COV2 + e}. When interaction factors are incorporated into the model, the main effects'
#' significance is not always determined simply; rather, it depends on the arbitrary coding of the variables. To put it
#' another way, you should probably just interpret the p-value for the interaction. Also, The p-values for the covariates
#' do not represent the test for the SNP-phenotype association after controlling for the covariate. That is the first row
#' (ADD). Rather, the covariate term is the test associated with the covariate-phenotype association. These p-values might
#' be extremely significant (e.g. if one covaries for smoking in an analysis of heart disease, etc) but this does not mean
#' that the SNP has a highly significant effect necessarily. Note that, this feature is not valid for XCGA model for XWAS part.
#'
#' @param Inphenocov Vector of integer values starting from 1 to extract the terms which user wants from the above model:
#' \eqn{Y = b0 + b1.ADD + b2.COV1 + b3.COV2 + b4.ADDxCOV1 + b5.ADDxCOV2 + e}. The terms will appear in order as
#' \insertCite{Purcell2007}{GXwasR} for ADD, \insertCite{Su2022}{GXwasR} for COV1, \insertCite{Rhodes2002}{GXwasR} for ADD x COV1,
#' and \insertCite{Moreau2003}{GXwasR} for ADD x COV2. If the user wants to extract the terms for COV1 and ADD x COV1, they need to specify it as c(2,4).
#' The default is `c(“ALL”)`.
#'
#' Note: This feature is not valid for the XCGA model for the XWAS part.
#'
#' @param combtest
#' Character vector specifying method for combining p-values for stratified GWAS with FM01comb and FM02comb XWAS models.
#' Choices are “stouffer.method”, "fisher.method" and "fisher.method.perm". For fisher.method the function for combining
#' p-values uses a statistic, \eqn{S = -2 ∑^k /log p}, which follows a \eqn{χ^2} distribution with 2k degrees of freedom \insertCite{Fisher1925}{GXwasR}.
#'
#' For fisher.method.perm, using p-values from stratified tests, the summary statistic for combining p-values is \eqn{S = -2 ∑ /log p}.
#' A p-value for this statistic can be derived by randomly generating summary statistics \insertCite{Rhodes2002}{GXwasR}. Therefore, a p-value is randomly
#' sampled from each contributing study, and a random statistic is calculated. The fraction of random statistics greater or
#' equal to S then gives the final p-value.
#'
#' For stouffer.method ,the function applies Stouffer’s method \insertCite{Stouffer1949}{GXwasR} to the p-values assuming that the p-values to be combined are
#' independent. Letting p1, p2, . . . , pk denote the individual (one- or two-sided) p-values of the k hypothesis tests to be
#' combined, the test statistic is then computed with \eqn{$z = ∑^k_{1}frac{z_{i}}{sqrt(k)}$} where \eqn{$z_{i}$ = Φ−1 (1 – $p_{i}$)} and
#' \eqn{Φ −1 (·)} denotes the inverse of the cumulative distribution function of a standard normal distribution. Under the joint null
#' hypothesis, the test statistic follows a standard normal distribution which is used to compute the combined p-value. This
#' functionality is taken from the R package poolr \insertCite{Cinar2022}{GXwasR}.
#'
#' Note that only p-values between 0 and 1 are allowed to be passed to these methods.
#'
#' @param MF.zero.sub
#' Small numeric value for substituting p-values of 0 in in stratified GWAS with FM01comb and FM02comb XWAS models.
#' The default is 0.00001. As log(0) results in Inf this replaces p-values of 0 by default with a small float.
#'
#' @param B
#' Integer value specifying the number of permutation in case of using fisher.method.perm method in stratified GWAS with
#' FM01comb and FM02comb XWAS models. The default is 10000.
#'
#' @param MF.na.rm
#' Boolean value, `TRUE` or `FALSE` for removing p-values of NA in stratified GWAS with FM01comb and FM02comb XWAS
#' in case of using Fisher’s and Stouffer’s methods. The default is FALSE.
#'
#' @param MF.p.corr
#' Character vector specifying method for correcting the summary p-values for FMfcomb and FMscomb models. Choices
#' are "bonferroni", "BH" and "none" for Bonferroni,  Benjamini-Hochberg and none, respectively. The default is "none".
#'
#' @param MF.mc.cores
#' Number of cores used for fisher.method.perm in stratified GWAS with FM01comb and FM02comb XWAS models.
#'
#' @param plot.jpeg
#' Boolean value, `TRUE` or `FALSE` for saving the plots in .jpeg file. The default is TRUE.
#'
#' @param plotname
#' A character string specifying the prefix of the file for plots. This file will be saved in DataDir. The default is
#' "GXwas.plot".
#'
#' @param snp_pval
#' Numeric value as p-value threshold for annotation. SNPs below this p-value will be annotated on the plot. The default
#' is 1e-08.
#'
#' @param annotateTopSnp
#' Boolean value, `TRUE` or 'FALSE. If TRUE, it only annotates the top hit on each chromosome that is below the
#' snp_pval threshold. The default is FALSE.
#'
#' @param suggestiveline
#' The default is 5 (for p-value 1e-05).
#'
#' @param genomewideline
#' The default is 7.3 (for p-value 5e-08).
#'
#' @param ncores
#' Integer value, specifying the number of cores for parallel processing. The default is 0 (no parallel computation).
#'
#' @importFrom progress progress_bar
#' @importFrom Rdpack reprompt
#'
#' @return A dataframe with GWAS (with XWAS for X-chromosomal variants) along with Manhattan and Q-Q plots.
#' In the case of the stratified test, the return is a list containing three dataframes, namely, FWAS, MWAS, and MFWAS with association
#' results in only female, only male, and both cohorts, respectively. This will be accompanied by Miami and Q-Q plots. The individual manhattan
#' and Q-Q-plots for stratified tests prefixed with xmodel type will be in the DataDir.
#'
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#'
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' standard_beta <- TRUE
#' xsex <- FALSE
#' sex <- TRUE
#' Inphenocov <- NULL
#' covartest <- NULL
#' interaction <- FALSE
#' MF.na.rm <- FALSE
#' B <- 10000
#' MF.zero.sub <- 0.00001
#' trait <- "binary"
#' xmodel <- "FMcombx02"
#' combtest <- "fisher.method"
#' snp_pval <- 1e-08
#' covarfile <- NULL
#' ncores <- 0
#' MF.mc.cores <- 1
#' ResultGXwas <- GXwas(
#'   DataDir = DataDir, ResultDir = ResultDir,
#'   finput = finput, xmodel = xmodel, trait = trait, covarfile = covarfile,
#'   sex = sex, xsex = xsex, combtest = combtest, MF.p.corr = "none",
#'   snp_pval = snp_pval, plot.jpeg = TRUE, suggestiveline = 5, genomewideline = 7.3,
#'   MF.mc.cores = 1, ncores = ncores
#' )
GXwas <- function(DataDir, ResultDir, finput, trait = c("binary", "quantitative"), standard_beta = TRUE,
                  xmodel = c("FMcombx01", "FMcombx02", "FMstatrified", "GWAScxci"), sex = FALSE, xsex = FALSE,
                  covarfile = NULL, interaction = FALSE, covartest = c("ALL"), Inphenocov = c("ALL"), combtest = c("fisher.method", "fisher.method.perm", "stouffer.method"),
                  MF.zero.sub = 0.00001, B = 10000, MF.mc.cores = 1, MF.na.rm = FALSE,
                  MF.p.corr = "none", plot.jpeg = FALSE, plotname = "GXwas.plot", snp_pval = 1e-08,
                  annotateTopSnp = FALSE, suggestiveline = 5, genomewideline = 7.3, ncores = 0) {
  # Initialize progress bar

  pb <- progress::progress_bar$new(
    format = "[:bar] :percent in :elapsed",
    total = 100, clear = FALSE, width = 60
  )

  # Validate inputs
  validationError <- validateGXwasInputs(DataDir, ResultDir, finput, trait, standard_beta, xmodel, sex, xsex, covarfile, interaction, covartest, Inphenocov, combtest, MF.zero.sub, B, MF.mc.cores, MF.na.rm, MF.p.corr, plot.jpeg, plotname, snp_pval, annotateTopSnp, suggestiveline, genomewideline, ncores)

  if (!is.null(validationError)) {
    stop(validationError)
  }

  tryCatch(
    {
      pb$tick(5)

      # setupPlink(ResultDir)
      

      pb$tick(6)

      if (xmodel[1] == "FMcombx01") {
        rlang::inform(rlang::format_error_bullets("Running FMcombx01 model"))

        x <- suppressWarnings(FMmain(
          DataDir = DataDir, ResultDir = ResultDir, finput = finput, trait = trait, standard_beta = standard_beta, xmodel = xmodel,
          sex = sex, xsex = xsex, covarfile = covarfile, interaction = interaction, covartest = covartest, Inphenocov = Inphenocov, plot.jpeg = plot.jpeg, plotname = plotname, snp_pval = snp_pval, annotateTopSnp = annotateTopSnp, suggestiveline = suggestiveline, genomewideline = genomewideline, ncores = ncores
        ))


        pb$tick(20)
      } else if (xmodel[1] == "FMcombx02") {
        rlang::inform(rlang::format_error_bullets("Running FMcombx02 model"))


        x <- suppressWarnings(FMmain(
          DataDir = DataDir, ResultDir = ResultDir, finput = finput, trait = trait, standard_beta = standard_beta, xmodel = xmodel,
          sex = sex, xsex = xsex, covarfile = covarfile, interaction = interaction, covartest = covartest, Inphenocov = Inphenocov, plot.jpeg = plot.jpeg, plotname = plotname, snp_pval = snp_pval, annotateTopSnp = annotateTopSnp, suggestiveline = suggestiveline, genomewideline = genomewideline, ncores = ncores
        ))


        pb$tick(20)
      } else if (xmodel[1] == "FMstatrified") {
        rlang::inform(rlang::format_error_bullets("Running FMstatrified model"))

        ## Making male and female files in ResultDir

        MFsplitPlink(DataDir = DataDir, ResultDir = ResultDir, finput = finput, foutput = "finput.female", sex = "females")
        pb$tick(10)
        gc(reset = TRUE)
        MFsplitPlink(DataDir = DataDir, ResultDir = ResultDir, finput = finput, foutput = "finput.male", sex = "males")
        gc(reset = TRUE)
        pb$tick(15)

        x <- suppressWarnings(FMcomb(
          DataDir = DataDir, ResultDir = ResultDir, trait = trait, standard_beta = standard_beta, xmodel = xmodel,
          covarfile = covarfile, interaction = interaction, covartest = covartest, Inphenocov = Inphenocov,
          plot.jpeg = plot.jpeg, plotname = plotname, snp_pval = snp_pval, annotateTopSnp = annotateTopSnp,
          combtest = combtest, B = B, MF.p.corr = MF.p.corr, MF.zero.sub = MF.zero.sub, MF.na.rm = MF.na.rm, MF.mc.cores = MF.mc.cores, suggestiveline = suggestiveline, genomewideline = genomewideline, ncores = ncores
        ))
        pb$tick(20)
        ftemp <- list.files(paste0(ResultDir, "/"), pattern = "finput")
        invisible(file.remove(paste0(ResultDir, "/", ftemp)))
      } else if (xmodel[1] == "GWAScxci") {
        if (trait[1] == "quantitative") {
          return(rlang::inform(rlang::format_error_bullets(c("x" = "For GWAScxci model, trait needs to be quantitative. Please correct the input file."))))
        } else {
          rlang::inform(rlang::format_error_bullets("Running GWAScxci model"))
        }
        x <- suppressWarnings(XCMAFun(
          DataDir = DataDir, ResultDir = ResultDir, finput = finput, standard_beta = standard_beta,
          sex = sex, covarfile = covarfile, interaction = interaction, covartest = covartest, Inphenocov = Inphenocov, plot.jpeg = plot.jpeg, plotname = plotname, snp_pval = snp_pval, annotateTopSnp = annotateTopSnp, suggestiveline = suggestiveline, genomewideline = genomewideline, ncores = ncores
        ))

        pb$tick(20)
      }

      patterns <- c("_ss", ".logistic", "_snps", "plink", "allsnpsresults.rda", "all_snps_results")
      removePatternFiles(ResultDir = ResultDir, patterns = patterns)
      pb$tick(100)
      return(x)
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
#' (i.e., SNP idenitifier), ‘BETA’ (i.e., effect-size or logarithm of odds ratio), ‘SE’ (i.e., standard error of BETA),
#' ‘P’ (i.e., p-values), 'NMISS' (i.e., effective sample size), 'L95' (i.e., lower limit of 95% confidence interval) and
#' 'U95' (i.e., upper limit of 95% confidence interval) are in mandatory column headers. These files needed to be in DataDir.
#' If the numbers of cases and controls are unequal, effective sample size should be \eqn{4 / (1/<# of cases> + 1/<# of controls>)}.
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
#' Character string, spycifying the plot name of the file containg forest plots for the SNPs. The default is
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
#' Character string specifing the name of the plain-text file with a column of SNP names for the forest plots.
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
#' A list object containing five dataframes. The first three dataframes, such as Mfixed, Mrandom and Mweighted contain results
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
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' data("GXwasRData")
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
#'   DataDir = DataDir, SummData = SummData, ResultDir = ResultDir,
#'   SNPfile = NULL, useSNPposition = TRUE, UseA1 = UseA1, GCse = GCse,
#'   plotname = "Meta_Analysis.plot", pval_filter, top_snp_pval, max_top_snps,
#'   chosen_snps_file = NULL, byCHR, pval_threshold_manplot
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
        write.table(SummData[[i]], paste0(ResultDir, "/", "SNPdata_", i), row.names = FALSE, col.names = TRUE, quote = FALSE)
      }

      SummData <- list.files(paste0(ResultDir, "/"), pattern = "SNPdata_")

      # Apply genomic control
      if (GCse == TRUE) {
        invisible(lapply(SummData, getGCse, ResultDir = ResultDir))
      } else {
        rlang::inform(rlang::format_error_bullets(c("i" = "No study-specific genomic control was applied.")))
        # file.copy(from=paste0(DataDir,"/",SummData),to=paste0(ResultDir,"/",SummData),overwrite = TRUE,copy.mode = TRUE)
      }

      if (is.null(SNPfile)) {
        SNPfilev <- NULL
      } else {
        SNPfilev <- paste0(DataDir, "/", SNPfile)
      }

      if (byCHR == FALSE) {
        MR <- metaFun(
          DataDir = DataDir, ResultDir = ResultDir, SummData = SummData,
          CHR = NULL, chromosome = NULL, nomap = nomap,
          UseA1v = UseA1v, extract = extract, SNPfilev = SNPfilev
        )
      } else {
        chromosome <- 1:23
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
        chosenS <- read.table(paste0(DataDir, "/", chosen_snps_file), header = TRUE)
        colnames(chosenS) <- "SNP"
        MRfiltered <- merge(chosenS, MR, by = "SNP")
      }

      # Visualization
      generatePlots(MRfiltered, Sbeta, ResultDir, plotname, useSNPposition, pval_threshold_manplot, chosen_snps_file)


      ## Produce all forest plots in .pdf
      grDevices::pdf(paste0(ResultDir, "/", plotname, ".pdf"), width = 10, height = 5)
      # graphics::par(mar = c(3, 2, 2, 3),oma=c(3,0,3,3))
      i <- 1:length(MRfiltered$SNP)
      invisible(suppressWarnings(lapply(i, allForestplot, MR2 = MRfiltered, Sbeta = Sbeta)))
      dev.off()

      rlang::inform(rlang::format_error_bullets(c("v" = paste0(plotname, " files containing the forest plots of the SNPs are produced in the directory ", ResultDir, "."))))


      # Check for problem SNPs
      problemFile <- paste0(ResultDir, "/MetaResult.prob")
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
        grDevices::jpeg(paste0(ResultDir, "/", plotname, ".jpeg"),
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
        rlang::inform(rlang::format_error_bullets(c("v" = paste0(plotname, " files containing the forest plots of the SNPs are produced in the directory ", ResultDir, "."))))

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
      return(list(Resultfixed = Mfixed, Resultrandom = Mrandom, Resultweighted = Mweighted, Metadata = Msummdata, ProblemSNP = MP))
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



#' ClumpLD: Clumping SNPs using linkage disequilibrium between SNPs
#'
#' @description
#' This function, which is based on empirical estimations of linkage disequilibrium between SNPs, groups the SNP-based results
#' across one or more datasets or analysis. This approach can be used in two basic scenarios: (i) To summarize the top X single
#' SNP findings from a genome-wide scan as fewer clusters of connected SNPs (i.e., to assess how many independent loci are
#' associated). (ii) To give researchers a simple approach to merge sets of data from multiple studies when those studies may
#' have used various marker sets for genotyping.
#'
#' The clumping process begins with the index SNPs that are significant at threshold p1 and have not yet been clumped. It then
#' creates clumps of all additional SNPs that are within a specified kb of the index SNP and that are in linkage disequilibrium
#' with the index SNP based on an r-squared threshold. Following that, these SNPs are filtered based on the outcome for that SNP.
#' As this method is greedy \insertCite{Purcell2007}{GXwasR}, each SNP will, at most, only appear in one clump. The P value and
#' ALLELES would always, at random, be chosen from the first input file if the same SNP appeared in several input files in SNPdata
#' argument. Instead of the best p-value, the function refer to the SNP that has the strongest LD to the index as the best proxy.
#' Based on the genotype data, the SNP with the highest LD will be the same for all input files.
#'
#' @author Banabithi Bose
#'
#' @param DataDir
#' A character string for the file path of the input PLINK binary files.
#'
#' @param finput
#' Character string, specifying the prefix of the input PLINK binary files which will be used to calculate linkage disequilibrium
#' between the SNPs. This actual genotype data may or may not be the same dataset that was used to generate the summary statistics.
#' This file needs to be in `DataDir`.
#'
#' @param SNPdata
#' A list of R dataframes containing a single or multiple summary statistics with SNP and P (i.e., p-values) in mandatory column
#' headers. Other columns could be present.
#'
#' @param ResultDir
#' A character string for the file path where all output files will be stored. The default is `tempdir()`.
#'
#' @param clump_p1
#' Numeric value, specifying the significance threshold for index SNPs. The default is 0.0001.
#'
#' @param clump_p2
#' Numeric value, specifying the secondary significance threshold for clumped SNPs. The default is 0.01
#'
#' @param clump_r2
#' Numeric value, specifying the LD threshold for clumping. The default is 0.50.
#'
#' @param clump_kb
#' Integer value, specifying the physical distance threshold in base-pair for clumping. The default is 250.
#'
#' @param clump_index_first
#' Boolean value, `TRUE` or `FALSE`, specifying whether to force the index SNP to appear first in each clump. This option should
#' typically be `TRUE` if clump_best is `TRUE.` Default is `TRUE`.
#'
#' @param clump_best
#' Boolean value, `TRUE` or `FALSE`, specifying whether to select and output the best SNP from each clump. Default is `TRUE`.
#'
#' @param byCHR
#' Boolean value, `TRUE` or `FALSE`, specifying whether to perform the clumping chromosome-wise.
#'
#' @return
#' A list with two dataframes.
#'
#' BestClump: a dataframe with eight columns showing the single best proxy SNP for each index SNP with
#' coulmns "INDEX"(Index SNP identifier), "PSNP"(Best proxy SNP), "RSQ LD"(r-squared) between index and proxy,
#' "KB"(Physical distance between index and proxy), P(p-value for proxy SNP), "ALLELES"(The associated haplotypes for the index and proxy SNP),
#' and "F"(Which file used for clumping from which this result came from).
#'
#' AllClump: a dataframe with eight columns providing a detailed summary
#' of each clump identified by PLINK. It includes "INDEX_SNP" (the identifier for the index SNP that represents the clump), "SNP"
#' (the SNP being reported, which for the index SNP is the same as INDEX_SNP), "DISTANCE" (the physical distance in base pairs between the index
#' SNP and the reported SNP, with 0.0 indicating the index itself), "RSQ" (the r-squared value showing the degree of linkage disequilibrium between
#' the index SNP and the SNP in the clump), "ALLELES" (the allele information, which in some cases may appear misaligned if the data isn’t formatted
#' as expected), "F" (a statistic or indicator related to the association test, which may be NA when not applicable), "P" (the p-value for the
#' association test of the SNP), and "CHR" (the chromosome on which the SNP is located).
#' @export
#'
#' @references
#' \insertAllCited{}

#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' data("GXwasRData")
#' SNPdata <- list(Summary_Stat_Ex1, Summary_Stat_Ex2)
#' clump_p1 <- 0.0001
#' clump_p2 <- 0.001
#' clump_r2 <- 0.5
#' clump_kb <- 250
#' byCHR <- TRUE
#' clumpedResult <- ClumpLD(
#'   DataDir, finput, SNPdata, ResultDir, clump_p1,
#'   clump_p2, clump_r2, clump_kb, byCHR
#' )
ClumpLD <- function(DataDir, finput, SNPdata, ResultDir = tempdir(),
                    clump_p1, clump_p2, clump_r2, clump_kb, byCHR = TRUE,
                    clump_best = TRUE, clump_index_first = TRUE) {
  # Check for an unsupported combination:
  if (clump_best == TRUE && clump_index_first == FALSE) {
    warning("The combination clump_best = TRUE and clump_index_first = FALSE is not recommended. Enforcing clump_index_first = TRUE.")
    clump_index_first <- TRUE
  }

  # Validate input parameters
  validateInputForClumpLD(
    DataDir, finput, SNPdata, ResultDir,
    clump_p1, clump_p2, clump_r2, clump_kb, byCHR
  )

  if (!checkFiles(DataDir, finput)) {
    stop("Missing required Plink files in the specified DataDir.")
  }

  tryCatch(
    {
      # setupPlink(ResultDir)
      

      for (i in seq_along(SNPdata)) {
        rlang::inform(rlang::format_error_bullets(paste0("Processing summary statistics ", i)))
        write.table(SNPdata[[i]],
          paste0(ResultDir, "/", "SNPdata_", i),
          row.names = FALSE, col.names = TRUE, quote = FALSE
        )
      }

      SumData <- list.files(paste0(ResultDir, "/"), pattern = "SNPdata_")

      if (byCHR == TRUE) {
        bimfile <- read.table(paste0(DataDir, "/", finput, ".bim"))
        chrnum <- 1:length(unique(bimfile$V1))

        chrwiseLD <- function(chrnum) {
          chromosome <- unique(bimfile$V1)[chrnum]
          rlang::inform(rlang::format_error_bullets(paste0("Running LD clumping for chromosome ", chromosome)))
          # Ensure 'plink' is executable and in the correct directory
          # plink_executable <- paste0(ResultDir, "/plink")

          # Construct the PLINK command arguments with toggles
          plink_args <- c(
            "--bfile", paste0(DataDir, "/", finput),
            "--chr", as.character(chromosome),
            "--clump", paste0(ResultDir, "/", SumData),
            "--clump-p1", as.character(clump_p1),
            "--clump-p2", as.character(clump_p2),
            "--clump-r2", as.character(clump_r2),
            "--clump-kb", as.character(clump_kb)
          )
          if (clump_best) {
            plink_args <- c(plink_args, "--clump-best")
          }
          if (clump_index_first) {
            plink_args <- c(plink_args, "--clump-index-first")
          }
          plink_args <- c(
            plink_args,
            "--clump-allow-overlap",
            "--clump-snp-field", "SNP",
            "--clump-field", "P",
            "--clump-verbose",
            "--out", paste0(ResultDir, "/ClumpLD")
          )

          # Execute the PLINK command
          invisible(sys::exec_wait(
            plink(),
            args = plink_args,
            std_out = FALSE, # Do not capture standard output
            std_err = FALSE # Do not capture standard error
          ))

          if (file.exists(paste0(ResultDir, "/", "ClumpLD.clumped"))) {
            # If both flags are ON, try reading the produced best file.
            # Otherwise, create a dummy BestClump dataframe.
            if (clump_best && clump_index_first &&
              file.exists(paste0(ResultDir, "/", "ClumpLD.clumped.best"))) {
              ldc2 <- as.data.frame(data.table::fread(
                paste0(ResultDir, "/", "ClumpLD.clumped.best"),
                header = TRUE
              ))
            } else {
              ldc2 <- data.frame("", "", "", "", "", "", "", "")
              colnames(ldc2) <- c("INDEX", "PSNP", "RSQ", "KB", "P", "ALLELES", "F", "CHR")
            }

            ldc2$CHR <- chromosome
            ldcall <- suppressWarnings(read_plink_clumped_clean(
              resultDir = ResultDir,
              filename = "ClumpLD.clumped"
            ))
            ldcall$CHR <- chromosome
            save(ldcall, file = paste0(ResultDir, "/ldcall_", chromosome, ".Rda"))
            invisible(do.call(file.remove, list(paste0(ResultDir, "/", "ClumpLD.clumped"))))
            return(ldc2)
          } else {
            rlang::inform(rlang::format_error_bullets(c("i" = paste0("No significant clump results for chromosome ", chromosome))))
            ldc2 <- data.frame("", "", "", "", "", "", "", "")
            colnames(ldc2) <- c("INDEX", "PSNP", "RSQ", "KB", "P", "ALLELES", "F", "CHR")
            return(ldc2)
          }
        }

        LDC <- data.table::rbindlist(lapply(chrnum, chrwiseLD))
        LDC <- LDC[!apply(LDC == "", 1, all), ]

        ldfiles <- list.files(paste0(ResultDir, "/"), pattern = "ldcall_")
        ldf <- function(ldfiles) {
          load(paste0(ResultDir, "/", ldfiles))
          return(ldcall)
        }
        All_ldc <- data.table::rbindlist(lapply(ldfiles, ldf), fill = TRUE)
      } else {
        # plink_executable <- paste0(ResultDir, "/plink")

        # Construct the PLINK command arguments with toggles
        plink_args <- c(
          "--bfile", paste0(DataDir, "/", finput),
          "--clump", paste0(ResultDir, "/", SumData),
          "--clump-p1", as.character(clump_p1),
          "--clump-p2", as.character(clump_p2),
          "--clump-r2", as.character(clump_r2),
          "--clump-kb", as.character(clump_kb)
        )
        if (clump_best) {
          plink_args <- c(plink_args, "--clump-best")
        }
        if (clump_index_first) {
          plink_args <- c(plink_args, "--clump-index-first")
        }
        plink_args <- c(
          plink_args,
          "--clump-allow-overlap",
          "--clump-snp-field", "SNP",
          "--clump-field", "P",
          "--clump-verbose",
          "--out", paste0(ResultDir, "/ClumpLD")
        )

        # Execute the PLINK command
        invisible(sys::exec_wait(
          plink(),
          args = plink_args,
          std_out = FALSE, # Do not capture standard output
          std_err = FALSE # Do not capture standard error
        ))

        if (file.exists(paste0(ResultDir, "/", "ClumpLD.clumped.best"))) {
          ldc2 <- read.table(paste0(ResultDir, "/", "ClumpLD.clumped.best"), header = TRUE)
          ldc2 <- ldc2[!apply(ldc2 == "", 1, all), ]
          ldcall <- suppressWarnings(read_plink_clumped_clean(
            resultDir = ResultDir,
            filename = "ClumpLD.clumped"
          ))
          save(ldcall, file = paste0(ResultDir, "/ldcall_genome.Rda"))
          invisible(do.call(file.remove, list(paste0(ResultDir, "/", "ClumpLD.clumped"))))

          LDC <- ldc2
          All_ldc <- ldcall
        } else {
          rlang::inform(rlang::format_error_bullets(c("i" = "No significant clump results")))
          # Create dummy dataframes if no best clump file exists.
          LDC <- data.frame("", "", "", "", "", "", "", "")
          colnames(LDC) <- c("INDEX", "PSNP", "RSQ", "KB", "P", "ALLELES", "F", "CHR")
          All_ldc <- data.frame()
        }
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

  ftemp <- list.files(paste0(ResultDir, "/"), pattern = "SNPdata_")
  invisible(do.call(file.remove, list(paste0(ResultDir, "/", ftemp))))
  ftemp <- list.files(paste0(ResultDir, "/"), pattern = "Clump")
  invisible(do.call(file.remove, list(paste0(ResultDir, "/", ftemp))))
  ftemp <- list.files(paste0(ResultDir, "/"), pattern = "ldcall")
  invisible(do.call(file.remove, list(paste0(ResultDir, "/", ftemp))))

  return(list(BestClump = LDC, AllClump = All_ldc))
}

#' DiffZeroOne: Assessing the Z-score for deviation from one and zero.
#'
#' @description
#' This function tests the null hypothesis that a measured statistics (example: genetic correlation,
#' rg for a trait) < 1 using a 1-tailed test compared with a normal distribution (z = (1 − measure statistics)/Standard error).
#' For multiple tests, users are encouraged to apply a Bonferroni multiple-testing correction.
#'
#' @param inputdata
#' A dataframe object, contaning three columns:
#' * `Trait` (i.e., the phenotype of interest)
#' * `Stat` (i.e., the measured statistics)
#' * `SE` (i.e., the standard error of the measured statistics)
#'
#' @param diffzero
#' Boolean value, `TRUE` or `FALSE`, specifying to perform diviation from 0 test.
#'
#' @param diffone
#' Boolean value, `TRUE` or `FALSE`, specifying to perform diviation from 1 test.
#'
#' @return
#' A dataframe with columns:
#' * `Trait`
#' * `Stat`
#' * `SE`
#' * `P0` (i.e, p-value for deviation from zero test)
#' * `P1` (i.e., p-value for deviation from 1 test)
#'
#' @export
#'
#' @examples
#' data("GXwasRData")
#' colnames(Example_rgdata) <- c("Trait", "Stat", "SE")
#' inputdata <- Example_rgdata
#' x <- DiffZeroOne(inputdata, FALSE, TRUE)
DiffZeroOne <- function(inputdata, diffzero = TRUE, diffone = TRUE) {
  validateInputForDiffZeroOne(inputdata, diffzero, diffone)

  tryCatch(
    {
      if (diffzero == FALSE & diffone == FALSE) {
        rlang::inform(rlang::format_error_bullets(c("x" = "Both diffzero and diffone cannot set to be FALSE.")))
        return()
      } else if (diffzero == TRUE && diffone == TRUE) {
        inputdata1 <- inputdata
        inputdata$x <- as.numeric(inputdata$Stat)
        inputdata$y <- as.numeric(inputdata$SE)
        inputdata$Zscore <- inputdata$x / inputdata$y
        inputdata1$p0 <- pnorm(inputdata$Zscore, mean = 0, sd = 1, lower.tail = FALSE) # one-tailed p-value

        inputdata$x <- 1 - as.numeric(inputdata$Stat)
        inputdata$y <- as.numeric(inputdata$SE)
        inputdata$Zscore <- inputdata$x / inputdata$y
        inputdata1$p1 <- pnorm(inputdata$Zscore, mean = 0, sd = 1, lower.tail = FALSE) # one-tailed p-value
      } else if (diffzero == FALSE && diffone == TRUE) {
        inputdata1 <- inputdata
        inputdata$x <- 1 - as.numeric(inputdata$Stat)
        inputdata$y <- as.numeric(inputdata$SE)
        inputdata$Zscore <- inputdata$x / inputdata$y
        inputdata1$p1 <- pnorm(inputdata$Zscore, mean = 0, sd = 1, lower.tail = FALSE) # one-tailed p-value
      } else if (diffzero == TRUE && diffone == FALSE) {
        inputdata1 <- inputdata
        inputdata$x <- as.numeric(inputdata$Stat)
        inputdata$y <- as.numeric(inputdata$SE)
        inputdata$Zscore <- inputdata$x / inputdata$y
        inputdata1$p0 <- pnorm(inputdata$Zscore, mean = 0, sd = 1, lower.tail = FALSE) # one-tailed p-value
      }

      return(inputdata1)
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


#' SexDiffZscore: Z-score-based sex difference test.
#'
#' @description
#' This function calculates the difference in any kind of measured entities,(example: including SNP heritability estimate,
#' genetic correlation, and GWAS β values) between sexes using a Z-score and its associated p-value statistic.
#' When STAT/SE is normally distributed and the test statistics are independent in sex, the test is well calibrated. If
#' the statistics are positively correlated, this test is conservative (1).
#'
#' We could define SNPs with SDEs as those variants at the extreme ends of the distribution with an absolute value of the
#' Z-score greater than 3(|Z-score| > 3), which is roughly equivalent to p <10−3, and represents 0.3% of all tested SNPs.
#' The input dataframes should only include X-chromosome in order to obtain results for sex differences based solely on
#' X-linked loci.
#'
#' @param inputdata
#' A dataframe with five columns:
#' * `ID` (i.e., SNP ID or the phenotype of interest, etc.)
#' * `Fstat` (i.e., the measured statistics in females)
#' * `Fse` (i.e., the standard error of the measured statistics in females)
#' * `Mstat` (i.e., the measured statistics in males)
#' * `Mse` (i.e., the standard error of the measured statistics in males)
#'
#'
#' @return
#' Original input dataframe with:
#' * `Zscore` (i.e., Z-score),
#' * `p` (i.e., p-value) and
#' * `adjP` (i.e., Bonferroni corrected p-value)
#' columns added.
#'
#' @export
#'
#' @examples
#' data("GXwasRData")
#' inputdata <- Example_h2data
#' x <- SexDiffZscore(inputdata)
#'
SexDiffZscore <- function(inputdata) {
  # Validate input parameters
  validateInputForSexDiffZscore(inputdata)

  tryCatch(
    {
      inputdata$x <- as.numeric(inputdata$Mstat) - as.numeric(inputdata$Fstat)
      inputdata$y <- sqrt((as.numeric(inputdata$Mse))^2 + (as.numeric(inputdata$Fse))^2)
      inputdata$Zscore <- inputdata$x / inputdata$y

      # Calculate p-values
      inputdata$p <- 2 * pnorm(-abs(inputdata$Zscore))
      inputdata$adjP <- p.adjust(inputdata$p, method = "bonferroni")

      inputdata <- inputdata[, c("ID", "Mstat", "Mse", "Fstat", "Fse", "Zscore", "p", "adjP")]
      return(inputdata)
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
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' byCHR <- TRUE
#' REMLalgo <- 0
#' nitr <- 3
#' ncores <- 3
#' data("GXwasRData")
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
#'   DataDir = DataDir, ResultDir = ResultDir, finput = finput, byCHR = byCHR,
#'   REMLalgo = 0, nitr = 10, phenofile = phenofile, cat_covarfile = NULL, quant_covarfile = NULL,
#'   partGRM = FALSE, autosome = TRUE, Xsome = TRUE, nGRM = 3,
#'   cripticut = 0.025, minMAF = NULL, maxMAF = NULL, excludeResidual = TRUE, ncores = ncores
#' )
GeneticCorrBT <- function(DataDir, ResultDir, finput, byCHR = FALSE,
                          REMLalgo = c(0, 1, 2), nitr = 100, phenofile, cat_covarfile = NULL, quant_covarfile = NULL,
                          computeGRM = TRUE, grmfile_name = NULL,
                          partGRM = FALSE, autosome = TRUE, Xsome = TRUE, nGRM = 3,
                          cripticut = 0.025, minMAF = NULL, maxMAF = NULL,
                          excludeResidual = FALSE, ncores = 2) {
  # Validate input parameters
  validateInputForGeneticCorrBT(DataDir, ResultDir, finput, byCHR, REMLalgo, nitr, phenofile, cat_covarfile, quant_covarfile, partGRM, autosome, Xsome, nGRM, cripticut, minMAF, maxMAF, excludeResidual, ncores)

  if (!checkFiles(DataDir, finput)) {
    stop("Missing required Plink files in the specified DataDir.")
  }

  tryCatch(
    withCallingHandlers(
      {
        # setupGCTA(ResultDir)
        ## ComputeBivarREMLone phenofile
        write.table(phenofile, file = paste0(ResultDir, "/", "GCphenofile"), row.names = FALSE, quote = FALSE)

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
          bimfile <- read.table(paste0(DataDir, "/", finput, ".bim"))

          chrnum <- 1:length(unique(bimfile$V1))
          # chrnum <- 1:3

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



#' SexRegress: Performing linear regression analysis with quantitative response variable.
#'
#' @description
#' This function could be used to check association of two variables. For instance, PRS with sex.
#'
#' @param fdata
#' R dataframe object. The column with header `response` should contain the response variable. All other column are the regressor.
#'
#' @param regressor_index
#' Integer value, specifying the column number of the main regressor variable.
#'
#' @param response_index
#' Integer value, specifying the column number of the response variable.
#'
#' @return
#' Numeric value containing the regression estimate ("Estimate"), standard error ("Std. Error"), statistics ("t value") and
#' p-value (\eqn{Pr(>|t|)})
#'
#' @importFrom stats lm
#'
#' @export
#'
#' @examples
#'
#' data("GXwasRData")
#' fdata <- Regression_Ex
#' fdata$SEX <- as.factor(as.character(fdata$SEX))
#' response_index <- 1
#' regressor_index <- 2
#'
#' x <- SexRegress(fdata, regressor_index, response_index)
SexRegress <- function(fdata, regressor_index, response_index) {
  # Validate input parameters
  validateInputForSexRegress(fdata, regressor_index, response_index)

  tryCatch(
    {
      names(fdata)[response_index] <- "response"
      nullfdata <- fdata[, -regressor_index]

      null.model <- stats::lm(nullfdata$response ~ ., data = nullfdata)
      model <- stats::lm(fdata$response ~ ., data = fdata)
      # model R2 is obtained as
      null.r2 <- summary(null.model)$r.squared
      model.r2 <- summary(model)$r.squared

      # R2 of response is simply calculated as the model R2 minus the null R2
      response.r2 <- model.r2 - null.r2
      model.result <- summary(model)$coeff[regressor_index, ]
      return(model.result)
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

#' FilterRegion: Filter chromosomal regions.
#'
#' @author Banabithi Bose
#'
#' @description
#' Filtering Pseudo-Autosomal Region (PAR), X-transposed region (XTR), Ampliconic, filter based on chromosome code or
#' user-defined regions from input plink files. Only one type of filtering can be done from three types, either by region
#' (using `regionfile` = `TRUE`), by chromosome (`filterCHR`) or by any combination of these three, `filterPAR`,
#' `filterXTR` and `filterAmpliconic.`

#'
#' @param DataDir
#' A character string for the file path of the input PLINK binary files.
#'
#' @param ResultDir
#' A character string for the file path where all output files will be stored. The default is `tempdir()`.
#'
#' @param finput
#' Character string, specifying the prefix of the input PLINK binary files.
#'
#' @param foutput
#' Character string, specifying the prefix of the output PLINK binary files if the filtering option for the SNPs is chosen.
#' The default is "FALSE".
#'
#' @param CHRX
#' Boolean value, `TRUE` or `FALSE` to filter/flag regions from chromosome X. The default is `TRUE`. Note: `CHRX` only in effect if
#' one of `filterPAR`, `filterXTR` or `filterAmpliconic` filter is in effect.
#'
#' @param CHRY
#' Boolean value, `TRUE` or `FALSE` to filter/flag regions from chromosome X. The default is `FALSE`. Note: CHRY only in effect
#' if one of `filterPAR`, `filterXTR` or `filterAmpliconic` filter is in effect.
#'
#' @param filterPAR
#' Boolean value, `TRUE` or `FALSE` to filter out PARs from input plink file. The default is `TRUE`.
#'
#' @param filterXTR
#' Boolean value, `TRUE` or `FALSE` to filter out XTRs from input plink file. The default is `TRUE`.
#'
#' @param filterAmpliconic
#' Boolean value, `TRUE` or `FALSE` to filter out Ampliconic regions from input plink file. The default is `TRUE`.
#'
#' @param regionfile
#' Character string, specifying the name of the .txt file containing the user-defined regions to be filtered out from input plink
#' file in bed format. The default is `FALSE`. If `regionfile` = `TRUE`, only this filtering will be in effect. Also, PAR, XTR and
#' Ampliconic SNPs from X-chomosome will be flagged and returned.
#'
#' @param filterCHR
#' Vector value with positive integer, specifying the chromosome code to filter/flag the SNPs. The default is 0, means no filtering
#' based on chromosome code. For non-zero values of this argument, the function will only consider the chromosome code to filter or
#' flag. All other filtering will not work. If filterCHR = TRUE, only this filtering will be in effect. Also, PAR, XTR and Ampliconic
#' SNPs from X-chomosome will be flagged and returned.
#'
#' @param Hg
#' Character value, '19', or '38', specifying which genome build to use for PAR, XTR and Ampliconic regions. The default is Hg = "19".
#'
#' @param exclude
#' Boolean value, `TRUE` or `FALSE` to filter and flag or only flag the SNPs. The default is `TRUE`.
#'
#' @return
#' A list of three dataframes: PAR containing SNPs from PAR regions; XTR containing SNPs from XTR region and Ampliconic containing
#' SNPs from Ampliconic region.
#'
#' For non-zero value of `filterCHR`, a dataframe containing the excluded/flagged SNPs will be returned.
#'
#' For `exclude` = `TRUE`, two sets of plink binary files will be produced in ResultDir. One set will have the remaining SNPs after
#' filtering and other one will have the discarded SNPs.
#'
#' @export
#'
#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' foutput <- "PostimputeEX_QC1"
#' x <- FilterRegion(
#'   DataDir = DataDir, ResultDir = ResultDir,
#'   finput = finput, foutput = foutput, CHRX = TRUE, CHRY = FALSE,
#'   filterPAR = TRUE, filterXTR = TRUE, filterAmpliconic = TRUE,
#'   regionfile = FALSE, filterCHR = NULL, Hg = "38", exclude = TRUE
#' )
FilterRegion <-
  function(DataDir,
           ResultDir,
           finput,
           foutput,
           CHRX = TRUE,
           CHRY = FALSE,
           filterPAR = TRUE,
           filterXTR = TRUE,
           filterAmpliconic = TRUE,
           regionfile = FALSE,
           filterCHR = NULL,
           Hg = "19",
           exclude = TRUE) {
    # Validate parameters
    validateFilterRegionParams(DataDir, ResultDir, finput, foutput, CHRX, CHRY, filterPAR, filterXTR, filterAmpliconic, regionfile, filterCHR, Hg, exclude)

    if (!checkFiles(DataDir, finput)) {
      stop("Missing required Plink files in the specified DataDir.")
    }

    tryCatch(
      {
        # setupPlink(ResultDir)
        

        DataDir1 <- system.file("extdata", package = "GXwasR")

        # Set filter parameters
        para <- setFilterParameters(CHRX, CHRY, filterCHR, regionfile, filterPAR, filterXTR, filterAmpliconic)
        CHRX <- para$CHRX
        fch <- para$fch
        rf <- para$rf

        if (is.null(filterCHR)) {
          if (Hg == "19") {
            HG <- "hg19"
          } else {
            HG <- "GRCh38"
          }

          if (CHRX == TRUE) {
            CHR <- "chrX"
          } else if (CHRY == TRUE) {
            CHR <- "chrY"
          }

          # print("line 3491")
          rlang::inform(rlang::format_error_bullets(c("i" = CHR)))
          if (exclude == TRUE) {
            if (regionfile == FALSE) {
              x <- readGenomicFeatures(DataDir1, CHRX, CHRY, CHR, HG)

              snps <- processRegionFilter(x, filterPAR, filterXTR, filterAmpliconic, ResultDir, DataDir, finput, foutput)
              par_snps <- snps[[1]]
              xtr_snps <- snps[[2]]
              ampliconic_snps <- snps[[3]]
              filterPAR <- snps[[4]]
              filterXTR <- snps[[5]]
              filterAmpliconic <- snps[[6]]

              y <- filterGenomicFeatures(x, filterPAR, filterXTR, filterAmpliconic)

              if (is.null(y)) {
                rangefile <- NULL
              } else {
                write.table(
                  y,
                  file = paste0(ResultDir, "/region.txt"),
                  quote = FALSE,
                  row.names = FALSE,
                  col.names = FALSE,
                  eol = "\r\n",
                  sep = " "
                )

                rangefile <- paste0(ResultDir, "/region.txt")
              }
            } else {
              rangefile <- paste0(DataDir, "/", regionfile)
              par_snps <- NULL
              xtr_snps <- NULL
              ampliconic_snps <- NULL
            }

            if (!is.null(rangefile)) {
              executePlinkExcludeExtract(ResultDir, DataDir, finput, rangefile, foutput)

              bim <- read.table(paste0(ResultDir, "/", foutput, ".bim"))
              bim1 <- read.table(paste0(DataDir, "/", finput, ".bim"))

              num_marker_excluded <- nrow(bim1) - nrow(bim)
              rlang::inform(
                rlang::format_error_bullets(c(
                  "i" = paste0(num_marker_excluded, " SNPs are discarded."),
                  "v" = paste0("Plink files with passed SNPs are in ", ResultDir, " prefixed as ", foutput),
                  "v" = paste0("Plink files with discarded SNPs are in ", ResultDir, " prefixed as ", foutput, "_snps_extracted")
                )))
            } else {
              rlang::inform(rlang::format_error_bullets(c("i" = "No SNPs to be discarded or flagged. No output plink files are created.")))
            }
            ## Modified in V7
          } else if (exclude == FALSE) {
            x <- readGenomicFeatures(DataDir1, CHRX, CHRY, CHR, HG)
            snps <- processRegionFilter(x, filterPAR, filterXTR, filterAmpliconic, ResultDir, DataDir, finput, foutput)
            par_snps <- snps[[1]]
            xtr_snps <- snps[[2]]
            ampliconic_snps <- snps[[3]]
            filterPAR <- snps[[4]]
            filterXTR <- snps[[5]]
            filterAmpliconic <- snps[[6]]

            rlang::inform(rlang::format_error_bullets(c("i" = "SNPs are only flagged for the desired region.")))
          }

          if (regionfile == FALSE) {
            ftemp <- list.files(paste0(ResultDir, "/"), pattern = "region")
            invisible(do.call(file.remove, list(paste0(ResultDir, "/", ftemp))))
          }

          return(list(PAR = par_snps, XTR = xtr_snps, Ampliconic = ampliconic_snps))
        } else {
          executePlinkChrFilter(ResultDir, DataDir, finput, filterCHR, foutput)

          return(NULL) ## Added in 5.0
        }
        # return(list(PAR = par_snps, XTR = xtr_snps, Ampliconic = ampliconic_snps ))

        if (exclude == TRUE) {
          bim <- read.table(paste0(ResultDir, "/", foutput, ".bim"))
          bim1 <- read.table(paste0(DataDir, "/", finput, ".bim"))

          num_marker_excluded <- nrow(bim1) - nrow(bim)
          rlang::inform(
            rlang::format_error_bullets(c(
              "i" = paste0(num_marker_excluded, " SNPs are discarded."),
              "v" = paste0("Plink files with passed SNPs are in ", ResultDir, " prefixed as ", foutput), 
              "v" = paste0("Plink files with discarded SNPs are in ", ResultDir, " prefixed as ", foutput, "_snps_extracted")
            )))
        } else if (exclude == FALSE) {
          bim <- read.table(paste0(ResultDir, "/", foutput, ".bim"))
          colnames(bim) <- c("CHR", "SNP", "START", "END", "A1", "A2")
          Flagged_SNPs <- bim
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



# Updated in 5.0
#' PlinkSummary: Summary of plink format genotype dataset
#'
#' @param DataDir A character string for the file path of the input PLINK binary files.
#' @param ResultDir A character string for the file path of the plink program to be set up.
#' @param finput Character string, specifying the prefix of the input PLINK binary files. This file needs to be in DataDir.
#'
#' @return This function is called for its side effect: printing summary statistics to the console. It returns `NULL` invisibly.
#' @export
#'
#' @importFrom rlang abort warn
#' 
#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#'
#' x <- PlinkSummary(DataDir, ResultDir, finput)
#'
PlinkSummary <- function(DataDir, ResultDir = tempdir(), finput) {
  # Validate DataDir
  if (!is.character(DataDir) || !dir.exists(DataDir)) {
    stop("DataDir must be a valid directory path.")
  }

  # Validate ResultDir
  if (!is.character(ResultDir) || (!dir.exists(ResultDir) && ResultDir != tempdir())) {
    stop("ResultDir must be a valid directory path or the default tempdir().")
  }

  # Validate finput and foutput
  if (!is.character(finput)) {
    stop("finput must be character strings.")
  }

  if (!checkFiles(DataDir, finput)) {
    stop("Missing required Plink files in the specified DataDir.")
  }

  tryCatch(
    {
      # setupPlink(ResultDir)
      

      fam <- as.data.frame(utils::read.table(file.path(DataDir, paste0(finput, ".fam"))))
      if (ncol(fam) == 5) {
        fam$V6 <- fam$V1
        fam <- fam[, c(6, 1:5)]
        colnames(fam) <- c("V1", "V2", "V3", "V4", "V5", "V6")
        rlang::inform(rlang::format_error_bullets(c("x" = ".fam file has five columns, please provide six columns in this to utilize GXwasR.")))
      }

      fam$V6 <- as.numeric(as.character(fam$V6))
      fam <- stats::na.omit(fam)
      fam4 <- fam[fam$V5 != 0 & fam$V6 != 0 & fam$V6 != -9, ]

      rlang::inform(rlang::format_error_bullets(c("i" = paste("Dataset:", finput))))

      # Analyze phenotype data
      analyzePhenotypeData(fam, fam4)

      # Process SNP data
      bim <- as.data.frame(utils::read.table(file.path(DataDir, paste0(finput, ".bim"))))
      No.of.chr <- length(unique(bim$V1))
      No.of.snps <- length(unique(bim$V2))
      No.of.samples <- length(unique(fam$V2))

      rlang::inform(
        rlang::format_error_bullets(c(
          "i" = paste("Number of chromosomes:", No.of.chr),
          " " = paste("  - Chr:", unique(bim$V1)),
          "i" = paste("Total number of SNPs:", No.of.snps),
          "i" = paste("Total number of samples:", No.of.samples)
        )))
      return(invisible(NULL))
    },
    error = function(e) {
      rlang::abort("An error occurred: ", e$message)
    },
    warning = function(w) {
      rlang::warn("Warning: ", w$message)
    }
  )
}


#' FilterAllele: Filtering out the multi-allelic variants
#'
#' @author Banabithi Bose
#'
#' @description
#' This function filters out the multi-allelic SNPs from the input dataset.
#'
#' @param DataDir
#' A character string for the file path of the input PLINK binary files.
#'
#' @param ResultDir
#' A character string for the file path where all output files will be stored. The default is `tempdir()`.
#'
#' @param finput
#' Character string, specifying the prefix of the input PLINK binary files.
#'
#' @param foutput
#' Character string, specifying the prefix of the output PLINK binary files. If multi-allelic variants are present,
#' this file will be produced after filtering out these variants.
#'
#' @return
#' `NULL`. After multi-allelic variant filtering, the filtered plink files with only biallelic SNPs will be saved in `ResultDir`.
#' @export
#'
#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' foutput <- "Filter_Test"
#' x <- FilterAllele(DataDir, ResultDir, finput, foutput)
FilterAllele <- function(DataDir, ResultDir, finput, foutput) {
  # Validate DataDir
  if (!is.character(DataDir) || !dir.exists(DataDir)) {
    stop("DataDir must be a valid directory path.")
  }

  # Validate ResultDir
  if (!is.character(ResultDir) || (!dir.exists(ResultDir) && ResultDir != tempdir())) {
    stop("ResultDir must be a valid directory path or the default tempdir().")
  }

  # Validate finput and foutput
  if (!is.character(finput) || !is.character(foutput)) {
    stop("finput and foutput must be character strings.")
  }

  if (!checkFiles(DataDir, finput)) {
    stop("Missing required Plink files in the specified DataDir.")
  }

  tryCatch(
    {
      # setupPlink(ResultDir)
      


      bimf <- read.table(paste0(DataDir, "/", finput, ".bim"))
      x1 <- bimf[nchar(bimf[, 5]) > 1 | nchar(bimf[, 6]) > 1, , drop = FALSE]

      if (nrow(x1) != 0) {
        write.table(x1$V2, file = paste0(ResultDir, "/snps_multiallelic"), quote = FALSE, col.names = FALSE, row.names = FALSE)
      } else {
        rlang::inform(rlang::format_error_bullets(c("i" = "There is no multi-allelic SNP present in the input dataset.")))
        return()
      }

      invisible(sys::exec_wait(
        plink(),
        args = c(
          "--bfile",
          paste0(DataDir, "/", finput),
          "--exclude",
          paste0(ResultDir, "/", "snps_multiallelic"),
          "--allow-no-sex", # 4.0
          "--make-bed",
          "--out",
          paste0(ResultDir, "/", foutput),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))

      bimf1 <- read.table(paste0(ResultDir, "/", foutput, ".bim"))

      rlang::inform(
        rlang::format_error_bullets(c(
          "i" = paste0("Input dataset has ", nrow(bimf), " SNPs."),
          "i" = paste0("Output dataset has ", nrow(bimf1), " SNPs."),
          "v" = paste0("Plink files with only biallelic SNPs are in ", ResultDir, " prefixed as ", foutput)
        )))
      return()
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


#' PvalComb
#'
#' @description
#' This function combines the p-values of two separate GWAS summary statistics (for instance male and female populations),
#' merges them, and then applies various statistical methods (like Stouffer's method, Fisher's method) to integrate the
#' p-values. It also includes functionality for generating plots like Manhattan plots and Q-Q plots.
#'
#'
#' @param SumstatMale
#' R dataframe object of summary statistics of male GWAS with five mandatory columns:
#' * `CHR` (numeric chromosome code)
#' * `SNP` (variant id)
#' * `A1` (allele)
#' * `POS` (base-pair position)
#' * `P` (p-value).
#'
#' Other coulmns may present.
#'
#' @param SumstatFemale
#' R dataframe object of summary statistics of female GWAS with five mandatory columns:
#' * `SNP`
#' * `A1`
#' * `TEST`
#' * `POS`
#' * `P`
#'
#' Other coulmns may present.
#'
#' @param combtest
#' Character vector specifying method for combining p-values for stratified GWAS models. Choices are “stouffer.method”,
#' "fisher.method" and "fisher.method.perm". For fisher.method the function for combining p-values uses a statistic,
#' \eqn{S = -2 ∑^k /log p}, which follows a \eqn{χ^2} distribution with 2k degrees of freedom \insertCite{Fisher1925}{GXwasR}.
#' For fisher.method.perm, using p-values from stratified tests, the summary statistic for combining p-values
#' is \eqn{S = -2 ∑ /log p}. A p-value for this statistic can be derived by randomly generating summary statistics \insertCite{Rhodes2002}{GXwasR}.
#' Therefore, a p-value is randomly sampled from each contributing study, and a random statistic is calculated. The
#' fraction of random statistics greater or equal to S then gives the final p-value.
#'
#' @param MF.p.corr
#' Character vector specifying method for correcting the summary p-values for FMfcomb and FMscomb models. Choices are
#' "bonferroni", "BH" and "none" for Bonferroni,  Benjamini-Hochberg and none, respectively. The default is "none".
#'
#' @param MF.zero.sub
#' Small numeric value for substituting p-values of 0 in GWAS summary statistics. The default is 0.00001. As \eqn{log(0)}
#' results in Inf this replaces p-values of 0 by default with a small float.
#'
#' @param MF.na.rm
#' Boolean value, `TRUE` or `FALSE` for removing p-values of NA in stratified GWAS summary satistics in case of using Fisher’s
#' and Stouffer’s methods. The default is `TRUE`.
#'
#' @param MF.mc.cores
#' Number of cores used for fisher.method.perm for combining p-values. The default is 1.
#'
#' @param B
#' Integer value specifying the number of permutation in case of using fisher.method.perm method. The default is 10000.
#'
#' @param plot.jpeg
#' Boolean value, `TRUE` or `FALSE` for saving the plots in .jpeg file. The default is `TRUE`.
#'
#' @param plotname
#' A character string specifying the prefix of the file for plots. This file will be saved in DataDir.
#' The default is "GXwas.plot".
#'
#' @param PlotDir
#' A character string specifying the path of the directory where the plots will be saved. The default is `tempdir()`.
#'
#' @param snp_pval
#' Numeric value as p-value threshold for annotation. SNPs below this p-value will be annotated on the plot. The default is 1e-08.
#'
#' @param annotateTopSnp
#' Boolean value, `TRUE` or `FALSE.` If `TRUE`, it only annotates the top hit on each chromosome that is below the snp_pval threshold. The default is `FALSE`.
#'
#' @param suggestiveline
#' Numeric value for suggestive cut-off line in GWAS manhattan plot. The default is 5 (for p-value 1e-05).
#'
#' @param genomewideline
#' Numeric value for genome-wide significant cut-off line in GWAS manhattan plot. The default is 7.3 (for p-value 5e-08).
#'
#' @param ncores
#' Integer value, specifying the number of cores for parallel processing. The default is 0 (no parallel computation).

#' @return
#' A dataframe with GWAS summary statistics (with XWAS for X-chromosomal variants) along with Manhattan and Q-Q plots.
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' data("GXwasRData")
#' SumstatMale <- Mfile
#' colnames(SumstatMale)[3] <- "POS"
#' SumstatFemale <- Ffile
#' colnames(SumstatFemale)[3] <- "POS"
#' # Use the below datasets for tutorial document
#' # SumstatMale <- Summary_Stat_Ex1
#' # colnames(SumstatMale)[3]<- "POS"
#' # SumstatFemale <- Summary_Stat_Ex2
#' # colnames(SumstatFemale)[3]<- "POS"
#' PvalComb_Result <- PvalComb(
#'   SumstatMale = SumstatMale, SumstatFemale = SumstatFemale,
#'   combtest = "fisher.method", MF.mc.cores = 1, snp_pval = 0.001, plot.jpeg = FALSE,
#'   suggestiveline = 3, genomewideline = 5.69897, ncores = 1
#' )
#'
PvalComb <- function(SumstatMale, SumstatFemale,
                     combtest,
                     MF.p.corr = "none",
                     MF.zero.sub = 0.00001,
                     MF.na.rm = TRUE,
                     MF.mc.cores = 1,
                     B = 1000,
                     plot.jpeg = TRUE,
                     plotname = "GXwas.plot",
                     PlotDir = tempdir(),
                     snp_pval,
                     annotateTopSnp = FALSE,
                     suggestiveline = 5,
                     genomewideline = 7.3,
                     ncores = 0) {
  # Validate inputs
  validation_result <- validatePvalCombInputs(SumstatMale, SumstatFemale, combtest, MF.p.corr, MF.zero.sub, MF.na.rm, MF.mc.cores, B, plot.jpeg, plotname, snp_pval, annotateTopSnp, suggestiveline, genomewideline, ncores)
  if (!is.null(validation_result)) {
    stop(validation_result)
  }

  tryCatch(
    {
      MaleWAS <- data.table::as.data.table(SumstatMale)
      FemaleWAS <- data.table::as.data.table(SumstatFemale)

      MFWAS <- merge(FemaleWAS, MaleWAS, by = c("SNP", "A1"))
      gc(reset = TRUE)
      pvals <- as.data.frame(MFWAS[, c("P.x", "P.y")])

      if (combtest[1] == "stouffer.method") {
        Pnew <- applyStoufferMethod(pvals, MF.p.corr, MF.zero.sub, MF.na.rm, MF.mc.cores, ncores)
        Result <- cbind(MFWAS, Pnew[, 3:4])
        gc(reset = TRUE)
        Result <- Result[, c("SNP", "CHR.x", "POS.x", "p2")]
        gc(reset = TRUE)
        colnames(Result) <- c("SNP", "CHR", "POS", "P")

        XWAS_ADD_X <- Result[Result$CHR == 23, ]
      } else if (combtest[1] == "fisher.method") {
        Pnew <-
          fisher.method(
            pvals = pvals,
            p.corr = MF.p.corr,
            zero.sub = MF.zero.sub,
            na.rm = MF.na.rm,
            mc.cores = MF.mc.cores
          )
        Result <- cbind(MFWAS, Pnew[, 3:4])
        gc(reset = TRUE)
        Result <- Result[, c("SNP", "CHR.x", "POS.x", "p.value")] # we could choose "p.adj" as well.
        colnames(Result) <- c("SNP", "CHR", "POS", "P")
        gc(reset = TRUE)
        XWAS_ADD_X <- Result[Result$CHR == 23, ]
      } else if (combtest[1] == "fisher.method.perm") {
        Pnew <-
          fisher.method.perm(
            pvals = pvals,
            p.corr = MF.p.corr,
            zero.sub = MF.zero.sub,
            B = B,
            mc.cores = MF.mc.cores,
            blinker = 1000
          )
        Result <- cbind(MFWAS, Pnew[, 3:4])
        gc(reset = TRUE)
        Result <- Result[, c("SNP", "CHR", "POS", "P")]
        gc(reset = TRUE)
        XWAS_ADD_X <- Result[Result$CHR == 23, ]
      }

      # From p-values, calculate chi-squared statistic
      chisq <- qchisq(1 - Result$P, 1)
      lamdaGC <- median(chisq) / qchisq(0.5, 1)
      chisq1 <- qchisq(1 - XWAS_ADD_X$P, 1)
      lamdaGC1 <- median(chisq1) / qchisq(0.5, 1)

      # Stratified GWAS plot
      # Manhattan and QQ-plots will be produced using P values from additive effect only. For all other tests, please use the final output.
      FemaleWAS <- na.omit(FemaleWAS[, c("SNP", "CHR", "POS", "P")])
      gc(reset = TRUE)
      MaleWAS <- na.omit(MaleWAS[, c("SNP", "CHR", "POS", "P")])
      gc(reset = TRUE)
      colnames(FemaleWAS) <- c("SNP", "CHR", "POS", "pvalue")
      colnames(MaleWAS) <- c("SNP", "CHR", "POS", "pvalue")
      FemaleWAS <- as.data.frame(FemaleWAS)
      FemaleWAS[FemaleWAS$CHR == "23", "CHR"] <- "X"
      MaleWAS <- as.data.frame(MaleWAS)
      MaleWAS[MaleWAS$CHR == "23", "CHR"] <- "X"

      # Stratified XWAS plot
      gwas.t2 <- FemaleWAS[FemaleWAS$CHR == "X", ]
      gwas.b2 <- MaleWAS[MaleWAS$CHR == "X", ]

      Result1 <- Result
      XWAS_ADD_X1 <- XWAS_ADD_X
      colnames(Result1) <- c("SNP", "CHR", "BP", "P")
      colnames(XWAS_ADD_X1) <- c("SNP", "CHR", "BP", "P")

      # Call to the helper function for generating plots
      generateGWASPlots(plot.jpeg, plotname, FemaleWAS, MaleWAS, gwas.t2, gwas.b2, Result1, XWAS_ADD_X1, snp_pval, annotateTopSnp, suggestiveline, genomewideline, lamdaGC, lamdaGC1, PlotDir)

      gc(reset = TRUE)
      return(na.omit(Result))
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



#' FilterSNP: Filter out SNPs.
#'
#' @param DataDir
#' A character string for the file path of the input PLINK binary files.
#'
#' @param ResultDir
#' A character string for the file path where all output files will be stored. The default is `tempdir()`.
#'
#' @param finput
#' Character string, specifying the prefix of the input PLINK binary files.
#'
#' @param foutput
#' Character string, specifying the prefix of the output PLINK binary files if the filtering option for the SNPs is chosen.
#' The default is "FALSE".
#'
#' @param SNPvec
#' R dataframe with SNP names to be excluded.
#'
#' @param extract
#' Boolean value, `TRUE` or `FALSE`, specifying whether to extract the snps or discard the snps. The default is `FALSE`.
#'
#' @return `NULL`. The filtered file will be saved in `ResultDir`.
#' @export
#'
#' @examples
#'
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' SNPvec <- c("rs6529954", "rs12858640", "rs5962098")
#' finput <- "GXwasR_example"
#' foutput <- "Filter_Test"
#' FilterSNP(DataDir, ResultDir, finput, foutput, SNPvec = SNPvec, extract = TRUE)
FilterSNP <- function(DataDir, ResultDir, finput, foutput, SNPvec, extract = FALSE) {
  # Validate inputs
  validation_result <- validateFilterSNPInputs(DataDir, finput, SNPvec, extract)
  if (!is.null(validation_result)) {
    stop(validation_result)
  }

  if (!checkFiles(DataDir, finput)) {
    stop("Missing required Plink files in the specified DataDir.")
  }

  tryCatch(
    {
      # setupPlink(ResultDir)

      if (extract == TRUE) {
        remov <- "--extract"
      } else if (extract == FALSE) {
        remov <- "--exclude"
      }
      write.table(SNPvec, file = paste0(ResultDir, "/snplist"), sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)

      # Exclude SNPs
      invisible(sys::exec_wait(
        plink(),
        args = c(
          "--bfile",
          paste0(DataDir, "/", finput),
          remov,
          paste0(ResultDir, "/snplist"),
          "--allow-no-sex", # 4.0
          "--make-bed",
          "--out",
          paste0(ResultDir, "/", foutput),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))

      bim <- read.table(paste0(ResultDir, "/", foutput, ".bim"))

      rlang::inform(
        rlang::format_error_bullets(c(
          "i" = paste0(nrow(bim), " SNPs are extracted"),
          "v" = paste0("Plink files with extracted SNPs are in ", ResultDir, " prefixed as ", foutput)
        )))
      return(NULL)
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


#' Download Hapmap phase 3 and 1000 Genome phase 3 Data
#'
#' @description
#' Downloads reference data sets from specified URLs based on the reference dataset name and working directory.
#' Currently supports 'HapMapIII_NCBI36' and 'ThousandGenome'.
#'
#' @param refdata
#' A character string specifying the reference dataset to download. Should be one of 'HapMapIII_NCBI36' or 'ThousandGenome'.
#'
#' @param wdir
#' A character string specifying the working directory where the reference data will be downloaded and extracted.
#'
#' @return
#' Invisible. The function prints a message upon successful download and extraction of the reference data.
#'
#' @export
#' @examples
#' \dontrun{
#' Download_reference("ThousandGenome", tempdir())
#' Download_reference("HapMapIII_NCBI36", tempdir())
#' }
Download_reference <- function(refdata, wdir = tempdir()) {
  # Input validation
  if (!is.character(refdata) || !is.character(wdir)) {
    stop("Both 'refdata' and 'wdir' must be character strings.")
  }


  if (!refdata %in% c("HapMapIII_NCBI36", "ThousandGenome")) {
    stop("Invalid 'refdata'. Choose either 'HapMapIII_NCBI36' or 'ThousandGenome'.")
  }

  if (!dir.exists(wdir)) {
    stop("The specified working directory does not exist.")
  }

  tryCatch(
    {
      os_type <- detect_os_type()
      options(timeout = 200)

      if (refdata == "HapMapIII_NCBI36") {
        url <- "https://figshare.com/ndownloader/files/40585145"
        file_ext <- ".zip"
      } else if (refdata == "ThousandGenome") {
        # url <- "https://figshare.com/ndownloader/files/40728539"
        url <- "https://figshare.com/ndownloader/files/46552177"
        file_ext <- ".tar.gz"
      }

      destfile <- paste0(wdir, "/", refdata, file_ext)

      # Downloading based on OS
      if (os_type == "unix") {
        utils::download.file(url, destfile, quiet = TRUE)
      } else if (os_type == "windows") {
        utils::download.file(url, destfile, quiet = TRUE, mode = "wb")
      } else {
        stop("Unsupported Operating System.")
      }

      # Extracting files
      if (file_ext == ".zip") {
        utils::unzip(destfile, exdir = wdir, junkpaths = TRUE)
      } else if (file_ext == ".tar.gz") {
        utils::untar(destfile, exdir = wdir)
      }

      invisible(file.remove(destfile))

      # List the extracted files in the working directory

      ##### Added in V7

      # Define new names

      if (refdata == "ThousandGenome") {
        new_names <- c("ThousandGenome.bed", "ThousandGenome.bim", "ThousandGenome.fam")

        # Original file paths
        extracted_files <- as.list(paste0(wdir, "/", c("Ref10Kgenome.bed", "Ref10Kgenome.bim", "Ref10Kgenome.fam")))

        # Rename files
        for (i in seq_along(extracted_files)) {
          old_name <- extracted_files[[i]]
          new_name <- file.path(wdir, new_names[i])
          file.rename(old_name, new_name)
        }
      }
      #########
      rlang::inform(rlang::format_error_bullets(c("i" = paste0("Reference data '", refdata, "' downloaded and extracted in ", wdir, "."))))
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


############ New Function added in 5.0
#' DummyCovar: Recode a categorical variable to a set of binary dummy variables.
#'
#' @description
#' When dealing with categorical variables in genetic analysis using, a common approach is to convert these into dummy variables
#' for proper analysis \insertCite{Purcell2007}{GXwasR}. This function creates K-1 new dummy variables for a variable with K categories.
#' One level is automatically excluded from the dummy variables which serves as the reference category for subsequent analyses. This setup
#' implicitly sets the excluded that category as the baseline against which other categories are compared.
#'
#' @param DataDir
#' A character string for the file path of the input PLINK binary files.
#'
#' @param bfile
#' Character string, specifying the prefix of the input PLINK binary files for which covariate file will be generated.
#'
#' @param incovar
#' Character string, specifying the prefix of the input covariate file. First two columns will be, FID (i.e., Family ID) and IID (i.e.,
#' Sample ID) and rest of the columns are covariates.
#'
#' @param outcovar
#' Character string, specifying the prefix of the Output covariate file
#'
#' @return
#' R dataframe object with covariates
#'
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' bfile <- "GXwasR_example"
#' incovar <- "covarfile_w_pc_age.txt"
#' outcovar <- "dummycovarfile"
#' dummy_covars <- DummyCovar(DataDir = DataDir, ResultDir = ResultDir,
#'                            bfile = bfile, incovar = incovar, 
#'                            outcovar = outcovar)
DummyCovar <- function(DataDir, ResultDir = DataDir, bfile, incovar, outcovar) {
  if (
    file.exists(paste0(DataDir, "/", bfile, ".bed")) &&
    file.exists(paste0(DataDir, "/", bfile, ".bim")) &&
    file.exists(paste0(DataDir, "/", bfile, ".fam"))
  ) {
    invisible(
      sys::exec_wait(
        plink(),
        args = c(
          "--bed",
          paste0(DataDir, "/", bfile, ".bed"),
          "--bim",
          paste0(DataDir, "/", bfile, ".bim"),
          "--fam",
          paste0(DataDir, "/", bfile, ".fam"),
          "--covar",
          paste0(DataDir, "/", incovar),
          "--write-covar",
          "--dummy-coding",
          "--out",
          paste0(ResultDir, "/", outcovar),
          "--silent"
        ),
    std_out = FALSE,
    std_err = FALSE
  ))} else {
    stop("There are no Plink files in DataDir.\nPlease specify correct directory path with input Plink files.")
  }

  if(file.exists(paste0(ResultDir, "/", outcovar, ".cov"))){
    x <- read.table(paste0(ResultDir, "/", outcovar, ".cov"), header = TRUE)
    rlang::inform(rlang::format_error_bullets(c("i" = paste0("Covariate file: ", outcovar, ".cov is in ", ResultDir))))
    return(x)
  }
  
}


## New function in V7
#' executePlinkMAF: Execute PLINK to Calculate Minor Allele Frequencies (MAF)
#'
#' @description
#' This function executes PLINK to calculate minor allele frequencies (MAF) for a given dataset. It sets up the necessary PLINK environment,
#' runs the PLINK command, and returns the MAF results as a DataFrame. Intermediate files generated by PLINK are cleaned up after execution.
#'
#' @param DataDir
#' Character. Directory containing the input PLINK files (.bed, .bim, .fam).
#'
#' @param ResultDir
#' Character. Directory to store the output files generated by PLINK.
#'
#' @param finput
#' Character. Base name of the PLINK input files (without extensions).
#'
#' @return
#' DataFrame containing the minor allele frequency (MAF) results.
#'
#' @export
#'
#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' maf_data <- executePlinkMAF(DataDir, ResultDir, finput)
executePlinkMAF <- function(DataDir, ResultDir, finput) {
  if (!checkFiles(DataDir, finput)) {
    stop("Missing required Plink files in the specified DataDir.")
  }

  # setupPlink(ResultDir)
  # Compute MAF using PLINK and return the results as a DataFrame, and clean up intermediate files afterward.

  # Generate a unique output file prefix based on the timestamp to prevent any overwrite
  foutput <- paste0("maf_output_", format(Sys.time(), "%Y%m%d_%H%M%S"))

  # Execute PLINK using --bfile for simplified input file specification
  tryCatch(
    {
      invisible(sys::exec_wait(
        plink(), # Path to the PLINK executable
        args = c(
          "--bfile", paste0(DataDir, "/", finput), # Base filename for .bed, .bim, and .fam
          "--freq", # Command to calculate allele frequencies
          "--out", paste0(ResultDir, "/", foutput),
          "--silent" # Suppress output to standard output
        ),
        std_out = TRUE,
        std_err = TRUE
      ))
    },
    error = function(e) {
      stop("An error occurred while executing Plink: ", e$message)
    }
  )

  # Define the output filename for the frequency result
  freqOutputFile <- paste0(ResultDir, "/", foutput, ".frq")

  # Read the .frq output file into R as a DataFrame
  if (file.exists(freqOutputFile)) {
    maf_data <- read.table(freqOutputFile, header = TRUE, stringsAsFactors = FALSE)
  } else {
    stop("The expected PLINK output file does not exist: ", freqOutputFile)
  }

  # Clean up intermediate files created by PLINK
  file.remove(list.files(ResultDir, pattern = paste0(foutput, "\\."), full.names = TRUE))


  # Return the MAF data as a DataFrame
  return(maf_data)
}


## Added in V7
#' LDPrune: Performs LD pruning on SNP data using PLINK
#' @description
#' This function utilizes PLINK to perform LD pruning on genetic data. It identifies and removes SNPs that are in high
#' linkage disequilibrium with each other within specified windows.
#'
#' @param DataDir
#' Character string representing the file path of the input PLINK binary files.
#'
#' @param finput
#' Character string specifying the prefix of the input PLINK binary files.
#'
#' @param ResultDir
#' Character string for the file path where all output files will be stored, defaulting to a temporary directory.
#'
#' @param window_size
#' Integer, specifying the number of SNPs to include in the sliding window.
#'
#' @param step_size
#' Integer, specifying the number of SNPs the window moves over in each step.
#'
#' @param r2_threshold
#' Numeric, specifying the R^2 threshold for LD pruning.
#'
#' @return Returns a character vector of SNP identifiers that remain after LD pruning or NULL if an error occurs.
#' @export
#'
#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' prunedSNPs <- LDPrune(DataDir, finput, ResultDir, 50, 5, 0.2)
LDPrune <- function(DataDir, finput, ResultDir = tempdir(), window_size = 50, step_size = 5, r2_threshold = 0.2) {
  # Validate input parameters
  validateInputForLDPrune(DataDir, finput, ResultDir, window_size, step_size, r2_threshold)

  # Check if required PLINK files are available
  if (!checkFiles(DataDir, finput)) {
    stop("Missing required PLINK files in the specified DataDir.")
  }


  # Constructing the PLINK command arguments
  args <- c(
    "--bfile", paste0(DataDir, "/", finput),
    "--indep-pairwise",
    window_size,
    step_size,
    r2_threshold,
    "--out", paste0(ResultDir, "/ld_prune"),
    "--silent"
  )

  # Execute the PLINK command using the helper function within tryCatch
  tryCatch(
    {
      # setupPlink(ResultDir)
      executePlinkAd(ResultDir, args)

      # Load the list of pruned SNPs if available
      prune_in_file <- paste0(ResultDir, "/ld_prune.prune.in")
      if (file.exists(prune_in_file)) {
        prunedSNPs <- read.table(prune_in_file, col.names = "SNP", stringsAsFactors = FALSE)
        # Return the list of pruned SNP identifiers
        return(prunedSNPs$SNP)
      } else {
        stop("LD pruning did not generate any output. Check PLINK logs for details.")
      }
    },
    error = function(e) {
      # Handle errors by returning a more user-friendly message
      rlang::inform(rlang::format_error_bullets(c("x" = paste0("Error during LD pruning: ", e$message))))
      return(NULL) # Return NULL to indicate failure
    }
  )
}


## New function in 7.0

#' SumstatGenCorr: Genetic Correlation Calculation from GWAS Summary Statistics
#'
#' @description
#' This function calculates the genetic correlation between two summary statistics
#' using a specified reference linkage disequilibrium (LD) matrix from the UK Biobank.
#'
#' @param ResultDir
#' Directory where results should be saved.
#'
#' @param referenceLD
#' Reference LD matrix identifier. These are the LD matrices and their eigen-decomposition from 335,265 genomic
#' British UK Biobank individuals. Two sets of reference panel are provided:
#' 1) 307,519 QCed UK Biobank Axiom Array SNPs. The size is about 7.5 GB after unzipping.
#' 2) 1,029,876 QCed UK Biobank imputed SNPs. The size is about 31 GB after unzipping. Although it takes more time,
#' using the imputed panel provides more accurate estimates of genetic correlations.
#' Therefore if the GWAS includes most of the HapMap3 SNPs, then it is recommend using the imputed reference panel.
#'
#' @param sumstat1
#' Data frame for the first set of summary statistics.
#' The input data frame should include following columns: SNP, SNP ID; A1, effect allele; A2, reference allele;
#' N, sample size; Z, z-score; If Z is not given, alternatively, you may provide: b, estimate of marginal effect in GWAS; se,
#' standard error of the estimates of marginal effects in GWAS.
#'
#' @param sumstat2
#' Data frame for the second set of summary statistics.
#' The input data frame should include following columns: SNP, SNP ID; A1, effect allele; A2, reference allele;
#' N, sample size; Z, z-score; If Z is not given, alternatively, you may provide: b, estimate of marginal effect in GWAS; se,
#' standard error of the estimates of marginal effects in GWAS.
#'
#' @param Nref
#' Sample size of the reference sample where LD is computed. If the default UK Biobank reference sample is used, Nref = 335265
#'
#' @param N0
#' Number of individuals included in both cohorts. The estimated genetic correlation is usually robust against misspecified N0.
#' If not given, the default value is set to the minimum sample size across all SNPs in cohort 1 and cohort 2.
#'
#' @param eigen.cut
#' Which eigenvalues and eigenvectors in each LD score matrix should be used for HDL.
#' Users are allowed to specify a numeric value between 0 and 1 for eigen.cut. For example, eigen.cut = 0.99 means using the
#' leading eigenvalues explaining 99% of the variance
#' and their correspondent eigenvectors. If the default 'automatic' is used, the eigen.cut gives the most stable heritability
#' estimates will be used.
#'
#' @param lim
#' Tolerance limitation, default lim = exp(-18).
#'
#' @param parallel
#' Boolean value, TRUE or FALSE for whether to perform parallel computation. The default is FALSE
#'
#' @param numCores
#' The number of cores to be used. The default is 2.
#' 
#' @details
#' This function requires access to the reference LD data via an environment variable.
#' You must set one of the following environment variables to the appropriate directory:
#' 
#' - `UKB_ARRAY_PATH` for the Axiom Array reference (`UKB_array_SVD_eigen90_extraction`)
#' - `UKB_IMPUTED_PATH` for the full imputed reference (`UKB_imputed_SVD_eigen99_extraction`)
#' - `UKB_IMPUTED_HAPMAP2_PATH` for the imputed HapMap2 subset (`UKB_imputed_hapmap2_SVD_eigen99_extraction`)
#'
#'
#' @return A list is returned with:
#' \itemize{
#' \item{`rg`}: The estimated genetic correlation.
#' \item{`rg.se`}: The standard error of the estimated genetic correlation.
#' \item{`P`}: P-value based on Wald test.
#' \item{`estimates.df`}: A detailed matrix includes the estimates and standard errors of heritabilities, genetic covariance
#' and genetic correlation.
#' \item{`eigen.use`}: The eigen.cut used in computation.
#' }
#'
#' @references
#' \insertRef{Ning2020}{GXwasR}
#' 
#' @export
#'
#' @examples
#' sumstat1 <- GXwasR:::simulateSumstats()
#' sumstat2 <- GXwasR:::simulateSumstats()
#' if (nzchar(Sys.getenv("UKB_IMPUTED_HAPMAP2_PATH"))) {
#'  res <- SumstatGenCorr(
#'    ResultDir = tempdir(), 
#'    referenceLD = "UKB_imputed_hapmap2_SVD_eigen99_extraction",
#'    sumstat1 = sumstat1,
#'    sumstat2 = sumstat2, 
#'    parallel = FALSE
#'  )
#' }

SumstatGenCorr <- function(
  ResultDir = tempdir(),
  referenceLD,
  sumstat1,
  sumstat2,
  Nref = 335265,
  N0 = min(sumstat1$N),
  eigen.cut = "automatic",
  lim = exp(-18),
  parallel = FALSE,
  numCores = 2
) {
  reference_paths <- list(
    UKB_imputed_hapmap2_SVD_eigen99_extraction = Sys.getenv("UKB_IMPUTED_HAPMAP2_PATH", unset = NA),
    UKB_imputed_SVD_eigen99_extraction = Sys.getenv("UKB_IMPUTED_PATH", unset = NA),
    UKB_array_SVD_eigen90_extraction = Sys.getenv("UKB_ARRAY_PATH", unset = NA)
  )

  reference_urls <- c(
    UKB_imputed_hapmap2_SVD_eigen99_extraction = "https://www.dropbox.com/s/kv5zhgu274wg9z5/UKB_imputed_hapmap2_SVD_eigen99_extraction.tar.gz?dl=1",
    UKB_imputed_SVD_eigen99_extraction = "https://www.dropbox.com/s/6js1dzy4tkc3gac/UKB_imputed_SVD_eigen99_extraction.tar.gz?dl=1",
    UKB_array_SVD_eigen90_extraction = "https://www.dropbox.com/s/fuvpwsf6r8tjd6c/UKB_array_SVD_eigen90_extraction.tar.gz?dl=1"
  )

  if (!referenceLD %in% names(reference_paths)) {
    stop("Invalid referenceLD. Choose from: ", paste(names(reference_paths), collapse = ", "))
  }

  LD_path <- reference_paths[[referenceLD]]

  if (is.na(LD_path) || !dir.exists(LD_path)) {
    stop(
      "The environment variable for '", referenceLD, "' is not set or the directory does not exist.\n",
      "Please set the appropriate environment variable before running this function.\n",
      "Reference data can be obtained from: ", reference_urls[[referenceLD]]
    )
  }

  # Compute genetic correlation
  res.HDL <- tryCatch(
    {
      if (!parallel) {
        HDL.rg(gwas1.df = sumstat1, gwas2.df = sumstat2,
               LD.path = LD_path, Nref = Nref, N0 = N0,
               eigen.cut = eigen.cut, lim = lim)
      } else {
        HDL.rg.parallel(gwas1.df = sumstat1, gwas2.df = sumstat2,
                        LD.path = LD_path, Nref = Nref, N0 = N0,
                        eigen.cut = eigen.cut, lim = lim, numCores = numCores)
      }
    },
    error = function(e) {
      message("Error in estimating Genetic Correlation: ", e$message)
      return(NULL)
    }
  )

  return(res.HDL)
}

# New Function Addedin 10.0
#' ComputeLD: Compute Linkage Disequilibrium (LD) for SNP Data
#'
#' @description
#' This function computes linkage disequilibrium (LD) statistics for SNP data using PLINK. It allows for computation across all
#' SNPs or within specific chromosomes.
#'
#' @param DataDir Character string representing the file path of the input PLINK binary files.
#' @param finput Character string specifying the prefix of the input PLINK binary files.
#' @param ResultDir Character string for the file path where all output files will be stored, defaulting to a temporary directory.
#' @param ByCHR Logical indicating whether to perform the computation by chromosome. The default is FALSE.
#' @param CHRnum If ByCHR is TRUE, specifies the chromosome number for which LD should be computed. The default is NULL.
#' @param r2_LD The threshold for r-squared LD values to report in the output.
#'
#' @return Returns a data frame containing the computed LD values among SNPs, read from the output file generated by PLINK.
#'
#' @export
#'
#' @examples
#' \donttest{
#' ComputeLD(
#'   DataDir = system.file("extdata", package = "GXwasR"), ResultDir = tempdir(),
#'   finput = "GXwasR_example", ByCHR = TRUE, CHRnum = 1, r2_LD = 0.2
#' )
#' }
ComputeLD <- function(DataDir, ResultDir, finput, ByCHR = FALSE, CHRnum = NULL, r2_LD) {
  if (ByCHR == FALSE) {
    chr <- NULL
    CHRnum <- NULL
  } else {
    chr <- "--chr"
    CHRnum <- CHRnum
  }
  invisible(sys::exec_wait(
    plink(),
    args = c(
      "--bfile",
      paste0(DataDir, "/", finput),
      chr, CHRnum,
      "--r2",
      "--ld-window-r2", r2_LD,
      "--out",
      paste0(ResultDir, "/", "snpcorr"),
      "--silent"
    ),
    std_out = FALSE,
    std_err = FALSE
  ))
  snpld <- read.table(paste0(ResultDir, "/snpcorr.ld"), header = TRUE)
  return(snpld)
}
