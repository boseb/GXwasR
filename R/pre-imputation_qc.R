## Function 128
## Added in 3.0
validateInputForQCsnp <- function(DataDir, ResultDir = tempdir(), finput, foutput, casecontrol = FALSE, hweCase = NULL, hweControl = NULL, hwe = NULL, maf = 0.05, geno = 0.1, monomorphicSNPs = FALSE, caldiffmiss = FALSE, diffmissFilter = FALSE, dmissX = FALSE, dmissAutoY = FALSE, highLD_regions, ld_prunning = FALSE, window_size = 50, step_size = 5, r2_threshold = 0.02) {
    # Validate directories
    if (!dir.exists(DataDir)) {
        stop("Error in DataDir: Directory does not exist.")
    }
    if (!dir.exists(ResultDir)) {
        stop("Error in ResultDir: Directory does not exist.")
    }

    # Validate file prefixes
    if (!is.character(finput)) {
        stop("Error in finput: Must be a character string.")
    }
    if (!is.character(foutput)) {
        stop("Error in foutput: Must be a character string.")
    }

    # Validate boolean parameters
    boolean_params <- list(casecontrol = casecontrol, monomorphicSNPs = monomorphicSNPs, caldiffmiss = caldiffmiss, diffmissFilter = diffmissFilter, dmissX = dmissX, dmissAutoY = dmissAutoY, ld_prunning = ld_prunning)
    for (param_name in names(boolean_params)) {
        param_value <- boolean_params[[param_name]]
        if (!is.logical(param_value)) {
            stop("Error in ", param_name, ": Must be a boolean value.")
        }
    }

    # Validate numeric parameters
    numeric_params <- list(hweCase = hweCase, hweControl = hweControl, hwe = hwe, maf = maf, geno = geno, r2_threshold = r2_threshold)
    for (param_name in names(numeric_params)) {
        param_value <- numeric_params[[param_name]]
        if (!is.null(param_value) && (!is.numeric(param_value) || param_value < 0 || param_value > 1)) {
            stop("Error in ", param_name, ": Must be a numeric value between 0 and 1.")
        }
    }

    # Validate highLD_regions if not NULL
    if (!is.null(highLD_regions) && !is.data.frame(highLD_regions)) {
        stop("Error in highLD_regions: Must be a dataframe.")
    }

    # Validate integer-like parameters
    if (!is.numeric(window_size) || window_size <= 0 || window_size != as.integer(window_size)) {
        stop("Error in window_size: Must be a positive whole number.")
    }
    if (!is.numeric(step_size) || step_size <= 0 || step_size != as.integer(step_size)) {
        stop("Error in step_size: Must be a positive whole number.")
    }

    return(TRUE)
}

## Function 138
## Added in 3.0
validateInputForQCsample <- function(DataDir, ResultDir, finput, foutput, imiss, het, small_sample_mod, IBD, IBDmatrix, ambi_out, legend_text_size, legend_title_size, axis_text_size, axis_title_size, title_size, filterSample) {
    # Validation for directory paths (should be strings)
    if (!is.character(DataDir) || !is.character(ResultDir)) {
        stop("DataDir and ResultDir must be strings representing directory paths.")
    }

    # Validation for file names (should be strings)
    if (!is.character(finput) || (!is.null(foutput) && !is.character(foutput))) {
        stop("finput and foutput must be strings representing file names.")
    }

    # Validation for numeric parameters
    numeric_params <- list(imiss = imiss, het = het)

    # Validate 'imiss' to be between 0 and 1 if not NULL
    param_value_imiss <- numeric_params[["imiss"]]
    if (!is.null(param_value_imiss)) {
        if (!is.numeric(param_value_imiss) || param_value_imiss < 0 || param_value_imiss > 1) {
            stop("imiss must be NULL or a numeric value between 0 and 1.")
        }
    }

    # Validate 'het' to be numeric
    param_value_het <- numeric_params[["het"]]
    if (!is.null(param_value_het)) {
        if (!is.null(param_value_het) && !is.numeric(param_value_het)) {
            stop("het must be a numeric value or NULL.")
        }
    }

    # Validation for size parameters (should be positive integers)
    size_params <- list(legend_text_size = legend_text_size, legend_title_size = legend_title_size, axis_text_size = axis_text_size, axis_title_size = axis_title_size, title_size = title_size)
    for (param_name in names(size_params)) {
        if (!is.numeric(size_params[[param_name]]) || size_params[[param_name]] <= 0) {
            stop(param_name, " must be a positive integer.")
        }
    }
}

#' QCsample: Quality control for samples in the PLINK binary files.
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
#' A plot of heterogysity estimate vs missingness across sample and a list containing five R dataframe objects, namely,
#' `HM` (samples with outlying heterozygosity and/or missing genotype rates), `Failed_Missingness` (samples with missing genotype rates),
#' `Failed_heterozygosity` (samples with outlying heterozygosity), `Missingness_results` (missingness results) and `Heterozygosity_results`
#' (heterozygosity results) with output PLINK files in ResultDir if filtering out the samples option is chosen.
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
#'     DataDir = DataDir, ResultDir = ResultDir, finput = finput,
#'     foutput = foutput, imiss = imiss, het = het, IBD = IBD,
#'     ambi_out = ambi_out
#' )
#' cleanupDataDir <- list.files(DataDir, pattern = 'PreimputeEX_QC', full.names = TRUE)
#' unlink(cleanupDataDir)
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
            if (small_sample_mod == TRUE) {
                SSM <- "small-sample"
            } else {
                SSM <- NULL
            }


            ############ Adding this in version 3.0 #############
            ## Filter out samples with missing phenotype

            fam1 <- nrow(read.table(normalizePath(file.path(DataDir, paste0(finput, ".fam")), mustWork = FALSE), header = FALSE))

            # finput <- if(ambi_out) processAmbiguousSamples(DataDir, ResultDir, finput, fam1) else finput ## Closing it

            # Prepare the arguments for the Plink command for missing data analysis
            missingDataArgs <- c(
                "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE),
                "--missing",
                "--out", normalizePath(file.path(ResultDir, "filtered_missing"), mustWork = FALSE),
                "--silent"
            )

            executePlink(missingDataArgs)


            # Prepare the arguments for the Plink command for heterozygosity analysis
            heterozygosityArgs <- c(
                "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE),
                "--dog",
                "--het",
                SSM,
                "--out", normalizePath(file.path(ResultDir, "filtered_hetero"), mustWork = FALSE),
                "--silent"
            )

            executePlink(heterozygosityArgs)

            miss <- readDataFile(normalizePath(file.path(ResultDir, "filtered_missing.imiss"), mustWork = FALSE))
            heter <- readDataFile(normalizePath(file.path(ResultDir, "filtered_hetero.het"), mustWork = FALSE))

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
                file = normalizePath(file.path(ResultDir, "failed_het_imiss"), mustWork = FALSE),
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
                het_plot <- createHeterozygosityPlot(
                    hetermiss, hetfail, imissfail, het, imiss,
                    legend_text_size, legend_title_size,
                    axis_text_size, axis_title_size, title_size
                )

                printSampleFilterResults(imissfail, hetfail, failed_het_imiss)
              
                if (nrow(hetermiss) == 0) {
                    hetermiss <- NULL
                } else {
                    hetermiss <- hetermiss
                }

                if (nrow(hetermiss) == 0) {
                    hetermiss1 <- NULL
                } else {
                    hetermiss1 <- hetermiss[, seq_len(2)]
                }

                if (nrow(imissfail) == 0) {
                    imissfail1 <- NULL
                } else {
                    imissfail1 <- imissfail[, seq_len(2)]
                }

                if (nrow(hetfail) == 0) {
                    hetfail1 <- NULL
                } else {
                    hetfail1 <- hetfail[, seq_len(2)]
                }

                ftemp <- c("failed_het_imiss", "filtered_hetero.log", "filtered_missing.lmiss", "filtered_missing.log")
                invisible(do.call(file.remove, list(normalizePath(file.path(ResultDir, ftemp), mustWork = FALSE))))

                if (file.exists(normalizePath(file.path(ResultDir, "filtered_missing.hh"), mustWork = FALSE))) {
                    ftemp <- c("filtered_missing.hh")
                    invisible(do.call(file.remove, list(normalizePath(file.path(ResultDir, ftemp), mustWork = FALSE))))
                }
                if (file.exists(normalizePath(file.path(ResultDir, "foutput.hh"), mustWork = FALSE))) {
                    ftemp <- c("foutput.log", "foutput.hh")
                    invisible(do.call(file.remove, list(normalizePath(file.path(ResultDir, ftemp), mustWork = FALSE))))
                }
                fmi <- read.table(normalizePath(file.path(ResultDir, "filtered_missing.imiss"), mustWork = FALSE))
                fhh <- read.table(normalizePath(file.path(ResultDir, "filtered_hetero.het"), mustWork = FALSE))

                ftemp <- c("filtered_missing.imiss", "filtered_hetero.het")
                invisible(do.call(file.remove, list(normalizePath(file.path(ResultDir, ftemp), mustWork = FALSE))))

                hm <- failed_het_imiss[, 2:1]
                colnames(hm) <- c("FID", "IID")
            } else {
                excludeSamplesArgs <- c(
                    "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE),
                    "--make-bed",
                    "--out", normalizePath(file.path(ResultDir, foutput), mustWork = FALSE),
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

            Outputsample <- nrow(read.table(normalizePath(file.path(ResultDir, paste0(foutput, ".fam")), mustWork = FALSE), header = FALSE))
            rlang::inform(rlang::format_error_bullets(c("i" = paste0("No. of samples in input PLINK files: ", fam1))))
            rlang::inform(rlang::format_error_bullets(c("i" = paste0("No. of samples in output PLINK files: ", Outputsample))))
            rlang::inform(rlang::format_error_bullets(c("v" = paste0("Output PLINK files, ", foutput, " with final samples are in ", ResultDir, "."))))
            if (filterSample == FALSE) {
                rlang::inform(rlang::format_error_bullets(c("i" = "Samples are flagged only.")))
            }
            return(list(
                HM = hm,
                Failed_Missingness = imissfail[, seq_len(2)],
                Failed_heterozygosity = hetfail[, seq_len(2)],
                Failed_IBD = failed_ibd[, seq_len(2)],
                Missingness_results = fmi,
                Heterozygosity_results = fhh,
                IBD_results = ibd,
                het_plot = if (!is.null(imiss) && !is.null(het)) {
                    het_plot
                } else {
                    NULL
                }
            ))
        },
        error = function(e) {
            rlang::abort(
                message = e$message, 
                class = "QCSample_error"
            )
        },
        warning = function(w) {
            rlang::warn(
                message = w$message, 
                class = "QCSample_warning", 
                .frequency = "regularly", 
                .frequency_id = "QCSample_warning"
            )
        }
    )
}

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
#' Users are required to provide SNPs ids or rsids in the input PLINK files.
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
#' Eucledean distance from the center of the PC1 and PC2 to the maximum Euclidean distance of the reference samples. Study samples
#' outside this distance will be considered as outlier. The default is 3.
#'
#' @importFrom data.table as.data.table .SD
#' @importFrom ggplot2 ggplot aes geom_hline geom_vline guides geom_point guide_legend scale_shape_manual
#' @importFrom vroom vroom
#'
#' @return A list containing three data frames: one with the IDs of outlier samples (Outlier_samples), another with samples
#' annotated with predicted ancestry (Samples_with_predicted_ancestry), and one with the IDs of non-outlier samples (Non_outlier_samples).
#' A PCA plot is also returned.
#'
#' @references
#' \insertAllCited{}
#'
#' @export
#'
#' @examples
#' data("highLD_hg19", package = "GXwasR")
#' data("example_data_study_sample_ancestry", package = "GXwasR")
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' reference <- "HapMapIII_NCBI36"
#' highLD_regions <- highLD_hg19
#' study_pop <- example_data_study_sample_ancestry # PreimputeEX
#' studyLD_window_size <- 50
#' studyLD_step_size <- 5
#' studyLD_r2_threshold <- 0.02
#' filterSNP <- TRUE
#' studyLD <- FALSE
#' referLD <- FALSE
#' referLD_window_size <- 50
#' referLD_step_size <- 5
#' referLD_r2_threshold <- 0.02
#' outlier <- TRUE
#' outlier_threshold <- 3
#' x <- AncestryCheck(
#'     DataDir = DataDir, ResultDir = ResultDir, finput = finput,
#'     reference = reference, highLD_regions = highLD_regions,
#'     study_pop = study_pop, studyLD = studyLD, referLD = referLD,
#'     outlierOf = "EUR", outlier = outlier, outlier_threshold = outlier_threshold
#' )
AncestryCheck <- function(
      DataDir,
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
      outlier_threshold = 3
) {
    tryCatch(
        {
            # Validate inputs
            validateAncestryCheckInputs(
                DataDir, ResultDir, finput, reference, filterSNP, studyLD, studyLD_window_size, studyLD_step_size,
                studyLD_r2_threshold, referLD, referLD_window_size, referLD_step_size, referLD_r2_threshold,
                highLD_regions, study_pop, outlier, outlierOf, outlier_threshold
            )

            if (!checkFiles(DataDir, finput)) {
                stop("Missing required Plink files in the specified DataDir.")
            }
            # Read study bim file
            sbim <- vroom::vroom(
                file = normalizePath(file.path(DataDir, paste0(finput, ".bim")), mustWork = FALSE),
                col_names = FALSE,
                delim = "\t",
                show_col_types = FALSE
            )
            names(sbim) <- paste0("V", seq_len(ncol(sbim)))

            # Verify SNP format of study data
            study_snp_format <- verify_snp_format(sbim)

            # Verify existence of required reference data
            ref_path <- validate_reference_data(reference)
            if (reference == "ThousandGenome") {
                reference <- "Ref10Kgenome"
            }

            # Read reference .bim file
            rbim <- vroom::vroom(
                file = normalizePath(file.path(ref_path, paste0(reference, ".bim")), mustWork = FALSE),
                col_names = FALSE,
                delim = "\t",
                show_col_types = FALSE
            )
            names(rbim) <- paste0("V", seq_len(ncol(rbim)))

            # Verify SNP format of reference data
            ref_snp_format <- verify_snp_format(rbim)

            if (study_snp_format == "chr:pos" & ref_snp_format == "rsID") {
                rbim <-
                    rbim %>%
                    mutate(V2 = paste0(.data$V1, ":", .data$V4))
            }

            if (!is.null(highLD_regions)) {
                write.table(
                    highLD_regions,
                    file = normalizePath(file.path(ResultDir, "high-LD-regions-temp.txt"), mustWork = FALSE),
                    quote = FALSE, row.names = FALSE, col.names = FALSE
                )

                highLD_regions <- normalizePath(file.path(ResultDir, "high-LD-regions-temp.txt"), mustWork = FALSE)
            } else {
                highLD_regions <- NULL
            }

            if (filterSNP == TRUE) {
                # Filter AT-GC
                filterATGCSNPs(
                    study_bim_path = normalizePath(file.path(DataDir, paste0(finput, ".bim")), mustWork = FALSE),
                    study_bed_path = normalizePath(file.path(DataDir, paste0(finput, ".bed")), mustWork = FALSE),
                    study_fam_path = normalizePath(file.path(DataDir, paste0(finput, ".fam")), mustWork = FALSE),
                    ref_bim_path = normalizePath(file.path(ref_path, paste0(reference, ".bim")), mustWork = FALSE), # make sure you've copied all 3!
                    ref_bed_path = normalizePath(file.path(ref_path, paste0(reference, ".bed")), mustWork = FALSE),
                    ref_fam_path = normalizePath(file.path(ref_path, paste0(reference, ".fam")), mustWork = FALSE),
                    ResultDir = ResultDir
                )

                # LD Pruning for Study and Reference Data
                studyLDMessage <- processLDstudyData(
                    ResultDir,
                    highLD_regions = normalizePath(file.path(ResultDir, "high-LD-regions-temp.txt"), mustWork = FALSE),
                    studyLD, studyLD_window_size, studyLD_step_size, studyLD_r2_threshold
                )
                referenceLDMessage <- processLDreferenceData(
                    ResultDir, highLD_regions, referLD, referLD_window_size, referLD_step_size, referLD_r2_threshold
                )

                # Printing the messages returned by the functions
                rlang::inform(studyLDMessage)
                rlang::inform(referenceLDMessage)
            } else if (filterSNP == FALSE) {
                executePlinkForUnfilteredData(DataDir, ResultDir, finput, reference)
            }

            # Find common SNPs between study and reference data
            commonSNPResults <- findCommonSNPs(ResultDir)
            common_snps <- commonSNPResults$common_snps
            pruned_study <- commonSNPResults$pruned_study
            pruned_ref <- commonSNPResults$pruned_ref

            # Process common SNPs
            processCommonSNPs(ResultDir)

            # Process for correcting chromosome mismatches
            correctedData <- correctChromosomeMismatches(ResultDir, common_snps, pruned_study, pruned_ref)
            S1 <- correctedData$S1
            S2 <- correctedData$S2
            snpSameNameDiffPos <- correctedData$snpSameNameDiffPos

            # Finding mis-matching allele positions
            updated_ref <-
                vroom::vroom(
                    file = normalizePath(file.path(ResultDir, "filtered_ref_temp4.bim"), mustWork = FALSE),
                    col_names = FALSE,
                    delim = "\t",
                    show_col_types = FALSE
                )
            names(updated_ref) <- paste0("V", seq_len(ncol(updated_ref)))
            S3 <- updated_ref[match(common_snps, updated_ref$V2), , drop = FALSE]
            colnames(S1) <- c("V1", "V2", "V3", "V4", "Sa", "Sb")
            colnames(S3) <- c("V1", "V2", "V3", "V4", "Ra", "Rb")

            S1 <- data.table::as.data.table(S1)
            S3 <- data.table::as.data.table(S3)
            # Using SNP name and chr no. for merging, not using base-pair position
            S4 <- merge(S1, S3, by = c("V1", "V4")) # using chr and position
            snps_flips <- S4[which(S4[, 5] != S4[, 9] & S4[, 6] != S4[, 10]), ]
            snp_allele_flips <- unique(snps_flips$V2)

            ## Allele Flips
            handleSnpAlleleFlips(ResultDir, snp_allele_flips)

            # Checking allele flips again after correcting
            flipped_ref <-
                vroom::vroom(
                    file = normalizePath(file.path(ResultDir, "filtered_ref_temp5.bim"), mustWork = FALSE),
                    col_names = FALSE,
                    delim = "\t",
                    show_col_types = FALSE
                )
            names(flipped_ref) <- paste0("V", seq_len(ncol(flipped_ref)))
            S5 <- flipped_ref[match(common_snps, flipped_ref$V2), , drop = FALSE]
            colnames(S5) <- c("V1", "V2", "V3", "V4", "Ra", "Rb")
            S5 <- data.table::as.data.table(S5)
            S6 <- merge(S1, S5, by = c("V1", "V4"))
            snps_flips_wrong <- S6[which(S6[, 5] != S6[, 9] &
                S6[, 6] != S6[, 10]), ]
            allele_flips_wrong <- unique(snps_flips_wrong$V2)

            # Handel Allele flips wrong
            handleAlleleFlipsWrong(ResultDir, allele_flips_wrong)

            # Checking number of SNPs in reference after clean-up.
            cleaned_ref <-
                vroom::vroom(
                    file = normalizePath(file.path(ResultDir, "filtered_ref_temp6.bim"), mustWork = FALSE),
                    col_names = FALSE,
                    delim = "\t",
                    show_col_types = FALSE
                )
            names(cleaned_ref) <- paste0("V", seq_len(ncol(cleaned_ref)))

            snps_final <- unique(cleaned_ref$V2)

            # Merge study and reference to perform PCA

            mergeDatasetsAndPerformPCA(ResultDir)

            # Process Reference
            ref_ancestry_EUR_AFR_ASIAN <- loadAndProcessReferenceAncestry(ResultDir, reference)
            # Plot PCA
            combined_pop <- prepareAncestryData(study_pop, ref_ancestry_EUR_AFR_ASIAN)
            tab <- loadPCAData(ResultDir, combined_pop)
            pop_type <- createPopulationTypeData(tab)

            pca_plot <- plotPCA(tab, pop_type)

            reportAlleleFlips(snp_allele_flips, ResultDir)

            # Example of using the function
            Outlier_samples1 <- detectOutliers(
                tab = tab, ResultDir = ResultDir, DataDir = DataDir,
                finput = finput, outlier = outlier, outlierOf = outlierOf,
                outlier_threshold = outlier_threshold
            )

            removeTempFiles(ResultDir, "study_ref_merge")

            # Define patterns for other files to remove
            patterns_to_remove <- c(
                "study_SNP", "ref_SNP", "common_snps",
                "snp_allele_flips", "allele_flips_wrong",
                "Outlier_ancestry", "snpSameNameDiffPos", "temp"
            )
            for (pattern in patterns_to_remove) {
                removeTempFiles(ResultDir, pattern)
            }

            len <- length(Outlier_samples1)
            Outlier_samples1[["pca_plot"]] <- pca_plot

            return(Outlier_samples1)
        },
        error = function(e) {
            rlang::abort(
                message = e$message, 
                class = "AncestryCheck_error"
            )
        },
        warning = function(w) {
            rlang::warn(
                message = w$message, 
                class = "AncestryCheck_warning", 
                .frequency = "regularly", 
                .frequency_id = "AncestryCheck_warning"
            )
        }
    )
}

## Function 134
## Added in 3.0
validateInputForSexCheck <- function(DataDir, ResultDir = tempdir(), finput, infer_sex = FALSE, compute_freq = FALSE, LD = TRUE, LD_window_size = 50, LD_step_size = 5, LD_r2_threshold = 0.02, fmax_F = 0.2, mmin_F = 0.8) {
    # Validate directories
    if (!dir.exists(DataDir)) {
        stop("Error in DataDir: Directory does not exist.")
    }
    if (!dir.exists(ResultDir)) {
        stop("Error in ResultDir: Directory does not exist.")
    }

    # Validate file prefix
    if (!is.character(finput)) {
        stop("Error in finput: Must be a character string.")
    }

    # Validate boolean parameters
    boolean_params <- list(infer_sex = infer_sex, compute_freq = compute_freq, LD = LD)
    for (param_name in names(boolean_params)) {
        param_value <- boolean_params[[param_name]]
        if (!is.logical(param_value)) {
            stop("Error in ", param_name, ": Must be a boolean value.")
        }
    }

    # Validate integer-like parameters for LD_window_size and LD_step_size
    if (!is.numeric(LD_window_size) || LD_window_size <= 0 || (LD_window_size != as.integer(LD_window_size) && LD_window_size != 50)) {
        stop("Error in LD_window_size: Must be a positive whole number or the default value of 50.")
    }
    if (!is.numeric(LD_step_size) || LD_step_size <= 0 || (LD_step_size != as.integer(LD_step_size) && LD_step_size != 5)) {
        stop("Error in LD_step_size: Must be a positive whole number or the default value of 5.")
    }

    # Validate numeric parameters
    numeric_params <- list(LD_r2_threshold = LD_r2_threshold, fmax_F = fmax_F, mmin_F = mmin_F)
    for (param_name in names(numeric_params)) {
        param_value <- numeric_params[[param_name]]
        if (!is.numeric(param_value) || param_value < 0 || param_value > 1) {
            stop("Error in ", param_name, ": Must be a numeric value between 0 and 1.")
        }
    }

    return(TRUE)
}

#' SexCheck: Compare sex assignments in the input PLINK files with those imputed from X chromosome inbreeding coefficients
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
#' @param infer_sex
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
#' infer_sex <- FALSE
#' compute_freq <- FALSE
#'
#' x <- SexCheck(
#'     DataDir = DataDir, ResultDir = ResultDir, finput = finput, infer_sex = infer_sex,
#'     compute_freq = compute_freq, LD_window_size = LD_window_size, LD_step_size = LD_step_size,
#'     LD_r2_threshold = 0.02, fmax_F = 0.2, mmin_F = 0.8
#' )
#'
#' # Checking if there is any wrong sex assignment
#' problematic_sex <- x[x$STATUS != "OK", ]
SexCheck <-
    function(DataDir,
    ResultDir = tempdir(),
    finput,
    infer_sex = FALSE,
    compute_freq = FALSE,
    LD = TRUE,
    LD_window_size = 50,
    LD_step_size = 5,
    LD_r2_threshold = 0.02,
    fmax_F = 0.2,
    mmin_F = 0.8) {
        # Validate inputs
        if (!validateInputForSexCheck(DataDir, ResultDir, finput, infer_sex, compute_freq, LD, LD_window_size, LD_step_size, LD_r2_threshold, fmax_F, mmin_F)) {
            return(NULL)
        }

        # Check if required files exist using checkFiles helper function
        if (!checkFiles(DataDir, finput)) {
            stop("Required PLINK files are missing in the specified DataDir.")
        }

        tryCatch(
            {
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

                if (infer_sex == FALSE) {
                    if (compute_freq == TRUE) {
                        freq_file_args <- c(
                            "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE),
                            "--freq",
                            "--out", normalizePath(file.path(ResultDir, "freq_file"), mustWork = FALSE),
                            "--silent"
                        )
                        executePlink(freq_file_args)

                        if (LD == TRUE) {
                            # LD pruning and sex check using executePlink
                            ld_check_sex_args <- c(
                                "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE),
                                "--indep-pairwise", LD_window_size, LD_step_size, LD_r2_threshold,
                                "--read-freq", normalizePath(file.path(ResultDir, "freq_file.frq"), mustWork = FALSE),
                                "--check-sex", fmax_F, mmin_F,
                                "--out", normalizePath(file.path(ResultDir, "cs"), mustWork = FALSE),
                                "--silent"
                            )
                            executePlink(ld_check_sex_args)
                        } else {
                            # Sex check using executePlink
                            sex_check_args <- c(
                                "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE),
                                "--read-freq", normalizePath(file.path(ResultDir, "freq_file.frq"), mustWork = FALSE),
                                "--check-sex", fmax_F, mmin_F,
                                "--out", normalizePath(file.path(ResultDir, "cs"), mustWork = FALSE),
                                "--silent"
                            )
                            executePlink(sex_check_args)
                        }
                    } else if (compute_freq == FALSE) {
                        if (LD == TRUE) {
                            # LD pruning and sex check using executePlink
                            ld_prune_sex_check_args <- c(
                                "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE),
                                "--indep-pairwise", LD_window_size, LD_step_size, LD_r2_threshold,
                                "--check-sex", fmax_F, mmin_F,
                                "--out", normalizePath(file.path(ResultDir, "cs"), mustWork = FALSE),
                                "--silent"
                            )
                            executePlink(ld_prune_sex_check_args)
                        } else {
                            # Execute sex check using executePlink
                            sex_check_args <- c(
                                "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE),
                                "--check-sex", fmax_F, mmin_F,
                                "--out", normalizePath(file.path(ResultDir, "cs"), mustWork = FALSE),
                                "--silent"
                            )
                            executePlink(sex_check_args)
                        }

                        check_sex <-
                            read.table(
                                file = normalizePath(file.path(ResultDir, "cs.sexcheck"), mustWork = FALSE),
                                stringsAsFactors = FALSE,
                                header = TRUE
                            )
                    }
                } else if (infer_sex == TRUE) {
                    if (LD == TRUE) {
                        plink_args <- c(
                            "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE),
                            "--indep-pairwise", LD_window_size, LD_step_size, LD_r2_threshold,
                            "--make-bed",
                            "--out", normalizePath(file.path(ResultDir, "csLD"), mustWork = FALSE),
                            "--silent"
                        )

                        executePlink(plink_args)

                        plink_args_sex_imputation <- c(
                            "--bfile", normalizePath(file.path(ResultDir, "csLD"), mustWork = FALSE),
                            "--impute-sex",
                            "--make-bed",
                            "--out", normalizePath(file.path(ResultDir, "seximputed_plink"), mustWork = FALSE),
                            "--silent"
                        )

                        executePlink(plink_args_sex_imputation)
                    } else {
                        plink_args_sex_imputation <- c(
                            "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE),
                            "--impute-sex",
                            "--make-bed",
                            "--out", normalizePath(file.path(ResultDir, "seximputed_plink"), mustWork = FALSE),
                            "--silent"
                        )

                        executePlink(plink_args_sex_imputation)
                    }

                    check_sex <-
                        read.table(
                            file = normalizePath(file.path(ResultDir, "seximputed_plink.sexcheck"), mustWork = FALSE),
                            stringsAsFactors = FALSE,
                            header = TRUE
                        )
                    rlang::inform(rlang::format_error_bullets(c("v" = "The output PLINK files with imputed sex, prefixed, seximputed_plink, are available in the ResultDir.")))
                }

                # Cleanup temporary files
                removeTempFiles(ResultDir, "cs")

                return(check_sex)
            },
            error = function(e) {
                rlang::abort(
                    message = e$message,
                    class = 'SexCheck_error'
                )
            },
            warning = function(w) {
                rlang::warn(
                    message = w$message,
                    class = 'PreImputationQC_warning',
                    .frequency = "regularly",
                    .frequency_id = "SexCheck_warning"
                )
            }
        )
    }