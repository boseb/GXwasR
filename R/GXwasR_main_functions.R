#' QCsnp: Quality control (QC) for SNPs.
#'
#' @author Banabithi Bose
#'
#' @description
#' This function performs QC of genotype data from PLINK binary files.
#' It can filter based on minor allele frequency, Hardy-Weinberg equilibrium, call rate, and
#' differential missingness between cases and controls. It can also perform linkage disequilibrium-based filtering.
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
#' Boolean value, `TRUE` or `FALSE` indicating if the input PLINK files has cases-control status or not.
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
#' missingness in cases vs controls. The default is `TRUE`.
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
#' in cases vs controls, respectively. Output PLINK binary files in the working directory.
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
#'     DataDir = DataDir, ResultDir = ResultDir, finput = finput, foutput = foutput,
#'     geno = geno, maf = maf, hweCase = hweCase, hweControl = hweControl,
#'     ld_prunning = ld_prunning, casecontrol = casecontrol, monomorphicSNPs = monomorphicSNPs,
#'     caldiffmiss = caldiffmiss
#' )
QCsnp <-
    function(
      DataDir,
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
      r2_threshold = 0.02
    ) {
        if (!validateInputForQCsnp(DataDir, ResultDir, finput, foutput, casecontrol, hweCase, hweControl, hwe, maf, geno, monomorphicSNPs, caldiffmiss, diffmissFilter, dmissX, dmissAutoY, highLD_regions, ld_prunning, window_size, step_size, r2_threshold)) {
            return(NULL)
        }

        if (!checkFiles(DataDir, finput)) {
            stop("There are no Plink files in DataDir. Please specify correct directory path with input Plink files.")
        }

        tryCatch(
            {
                fam <-
                    as.data.frame(utils::read.table(file = normalizePath(file.path(DataDir, paste0(finput, ".fam")), mustWork = FALSE)))

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
                    read.table(file = normalizePath(file.path(ResultDir, "filtered_temp4.frq"), mustWork = FALSE), stringsAsFactors = FALSE, header = TRUE)

                mmSNPs <- freq[which(freq[, 5] == 0), 2, drop = FALSE]

                write.table(
                    mmSNPs,
                    file = normalizePath(file.path(ResultDir, "monomorphicSNPs"), mustWork = FALSE),
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

                rlang::inform(rlang::format_error_bullets(c("v" = paste0("Output PLINK files prefixed as ,", foutput, ", with passed SNPs are saved in ResultDir."))))

                # Remove filtered_temp* files
                removeTempFiles(ResultDir, "filtered_temp")

                # Remove specific files if they exist
                removeTempFiles(ResultDir, "SNPdifCallrate")
                removeTempFiles(ResultDir, "monomorphicSNPs")

                # Remove NoAmbiguousSNP* files
                removeTempFiles(ResultDir, "NoAmbiguousSNP")

                # Remove PLINK file
                removeTempFiles(ResultDir, "PLINK")

                # Remove other files
                removeTempFiles(ResultDir, "study_SNP")
                removeTempFiles(ResultDir, "sink_file.txt")

                resultbim <- read.table(normalizePath(file.path(ResultDir, paste0(foutput, ".bim")), mustWork = FALSE))
                inputbim <- read.table(normalizePath(file.path(DataDir, paste0(finput, ".bim")), mustWork = FALSE))

                rlang::inform(rlang::format_error_bullets(c("i" = paste0("Input file has ", length(unique(inputbim$V2)), " SNPs."))))
                rlang::inform(rlang::format_error_bullets(c("i" = paste0("Output file has ", length(unique(resultbim$V2)), " SNPs after filtering."))))


                return(list(MonomorSNPs = mmSNP1, DiffMissSNPs = SNPmissCC))
            },
            error = function(e) {
                rlang::abort(
                    message = e$message,
                    class = "QCsnp_error"
                )
            },
            warning = function(w) {
                rlang::warn(
                    message = w$message,
                    .frequency = "regularly",
                    .frequency_id = "QCsnp_warning"
                )
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
#' data("highLD_hg19", package = "GXwasR")
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' highLD_regions <- highLD_hg19
#' ld_prunning <- "TRUE"
#' window_size <- 50
#' step_size <- 5
#' r2_threshold <- 0.02
#' countPC <- 20
#' ## Genetic PC
#' GP <- ComputeGeneticPC(
#'     DataDir = DataDir, ResultDir = ResultDir,
#'     finput = finput, highLD_regions = highLD_hg19, countPC = 20
#' )
ComputeGeneticPC <- function(DataDir, ResultDir = tempdir(), finput, countPC = 10, plotPC = TRUE,
    highLD_regions = NULL, ld_prunning = TRUE,
    window_size = 50, step_size = 5, r2_threshold = 0.02) {
    # Validate inputs
    if (!validateInputForComputeGeneticPC(DataDir, ResultDir, finput, countPC, plotPC, highLD_regions, ld_prunning, window_size, step_size, r2_threshold)) {
        stop("Please verify all inputs.")
    }

    if (!checkFiles(DataDir, finput)) {
        stop("Missing required Plink files in the specified DataDir.")
    }

    tryCatch(
        {
            ## Handle high LD regions exclusion
            processed_file <- normalizePath(file.path(DataDir, finput), mustWork = FALSE)
            if (!is.null(highLD_regions)) {
                # Write high LD regions to a temporary file
                options(scipen = 100)
                write.table(highLD_regions,
                    file = normalizePath(file.path(ResultDir, "highLD_regions_temp"), mustWork = FALSE),
                    quote = FALSE, row.names = FALSE, col.names = FALSE
                )
                highLD_regions_file <- normalizePath(file.path(ResultDir, "highLD_regions_temp"), mustWork = FALSE)
                options(scipen = 0)

                # Exclude high LD regions
                invisible(sys::exec_wait(
                    plink(),
                    args = c(
                        "--bfile", processed_file,
                        "--exclude", "range", highLD_regions_file,
                        "--make-bed",
                        "--out", normalizePath(file.path(ResultDir, paste0("no_highLD_", finput)), mustWork = FALSE),
                        "--silent"
                    ),
                    std_out = FALSE,
                    std_err = FALSE
                ))

                # Update the processed file to the one without high LD regions
                processed_file <- normalizePath(file.path(ResultDir, paste0("no_highLD_", finput)), mustWork = FALSE)
            }

            ## LD pruning if enabled
            if (ld_prunning == TRUE) {
                # Perform LD pruning
                invisible(sys::exec_wait(
                    plink(),
                    args = c(
                        "--bfile", processed_file,
                        "--indep-pairwise",
                        window_size,
                        step_size,
                        r2_threshold,
                        "--allow-no-sex",
                        "--out", normalizePath(file.path(ResultDir, paste0("pruned_", finput)), mustWork = FALSE),
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
                        "--extract", normalizePath(file.path(ResultDir, paste0("pruned_", finput, ".prune.in")), mustWork = FALSE),
                        "--make-bed",
                        "--out", normalizePath(file.path(ResultDir, paste0("final_pruned_", finput)), mustWork = FALSE),
                        "--silent"
                    ),
                    std_out = FALSE,
                    std_err = FALSE
                ))

                # Update the processed file to the pruned dataset
                processed_file <- normalizePath(file.path(ResultDir, paste0("final_pruned_", finput)), mustWork = FALSE)
            }

            ## PCA calculation
            invisible(sys::exec_wait(
                plink(),
                args = c(
                    "--bfile", processed_file,
                    "--pca", countPC,
                    "--out", normalizePath(file.path(ResultDir, "pcfile"), mustWork = FALSE),
                    "--silent"
                ),
                std_out = FALSE,
                std_err = FALSE
            ))

            # Process PCA results
            PCs1 <- read.table(normalizePath(file.path(ResultDir, "pcfile.eigenvec"), mustWork = FALSE))
            PCs <- PCs1[, -c(seq_len(2))]
            names(PCs) <- paste0("PC", seq_len(ncol(PCs)))
            EV <- scan(normalizePath(file.path(ResultDir, "pcfile.eigenval"), mustWork = FALSE))
            Percent.var <- data.frame(PC = seq_len(ncol(PCs)), Percent.var = EV / sum(EV) * 100)

            # Plot PCA results if requested
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

                combined_plot <- ggpubr::ggarrange(p1, p2, labels = c("A", "B"), ncol = 2, nrow = 1)
            }

            # Clean up temporary files
            if (file.exists(normalizePath(file.path(ResultDir, "pcfile.log"), mustWork = FALSE))) {
                file.remove(list.files(normalizePath(file.path(ResultDir), mustWork = FALSE), pattern = "pcfile", full.names = TRUE))
            }
            file.remove(list.files(normalizePath(file.path(ResultDir), mustWork = FALSE), pattern = "temp", full.names = TRUE))

            return(
                list(
                    PCs1 = PCs1,
                    plot = if (plotPC) {
                        combined_plot
                    } else {
                        NA
                    }
                )
            )
        },
        error = function(e) {
            rlang::abort(
                message = e$message,
                class = "ComputeGeneticPC_error"
            )
        },
        warning = function(w) {
            rlang::warn(
                message = w$message,
                .frequency = "regularly",
                .frequency_id = "ComputeGeneticPC_warning"
            )
            if (str_detect(conditionMessage(w), "cannot remove")) {
                return(
                    list(
                        PCs1 = PCs1,
                        plot = if (plotPC) {
                            combined_plot
                        } else {
                            NA
                        }
                    )
                )
            }
        }
    )
}


#' MergeRegion: Merging two sets of PLINK binary files.
#'
#' @author Banabithi Bose
#'
#' @description This function combines the two genotype datasets based on either common SNPs or all the SNPs between them.
#'
#' @param DataDir A character string for the file path of the input PLINK binary files.
#' @param ResultDir A character string for the file path where all output files will be stored. The default is tempdir().
#' @param finput1 Character string, specifying the prefix of the first input PLINK binary files.
#' @param finput2 Character string, specifying the prefix of the first input PLINK binary files.
#' @param foutput Character string, specifying the prefix of the output PLINK binary files if filtering option for the SNPs is chosen. The default is "FALSE".
#' @param use_common_snps Boolean value, TRUE or FALSE, specifying to use common SNPs for merging or to use all the SNPs.
#'
#' @return NULL
#'
#' The output PLINK files will be saved in ResultDir.
#'
#' @export
#'
#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput1 <- "GXwasR_example"
#' finput2 <- "GXwasR_example_imputed"
#' foutput <- "Test_output"
#' y <- MergeRegion(DataDir, ResultDir, finput1, finput2, foutput, use_common_snps = TRUE)
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
            if (use_common_snps) {
                bim1 <- read.table(normalizePath(file.path(DataDir, paste0(finput1, ".bim")), mustWork = FALSE))
                bim2 <- read.table(normalizePath(file.path(DataDir, paste0(finput2, ".bim")), mustWork = FALSE))
                common_snps <- intersect(bim1$V2, bim2$V2)
                write.table(common_snps,
                    file = normalizePath(file.path(ResultDir, paste0("common_snps_", foutput)), mustWork = FALSE),
                    quote = FALSE, row.names = FALSE, col.names = FALSE
                )

                args1 <- c(
                    "--bfile", normalizePath(file.path(DataDir, finput1), mustWork = FALSE),
                    "--extract", normalizePath(file.path(ResultDir, paste0("common_snps_", foutput)), mustWork = FALSE),
                    "--allow-no-sex", "--make-bed",
                    "--out", normalizePath(file.path(ResultDir, paste0("new", finput1)), mustWork = FALSE),
                    "--silent"
                )
                executePlink(args1)

                args2 <- c(
                    "--bfile", normalizePath(file.path(DataDir, finput2), mustWork = FALSE),
                    "--extract", normalizePath(file.path(ResultDir, paste0("common_snps_", foutput)), mustWork = FALSE),
                    "--allow-no-sex", "--make-bed",
                    "--out", normalizePath(file.path(ResultDir, paste0("new", finput2)), mustWork = FALSE),
                    "--silent"
                )
                executePlink(args2)

                merge_args <- c(
                    "--bfile", normalizePath(file.path(ResultDir, paste0("new", finput1)), mustWork = FALSE),
                    "--bmerge", normalizePath(file.path(ResultDir, paste0("new", finput2)), mustWork = FALSE),
                    "--allow-no-sex", "--make-bed",
                    "--out", normalizePath(file.path(ResultDir, foutput), mustWork = FALSE),
                    "--silent"
                )
                executePlink(merge_args)

                rlang::inform(rlang::format_error_bullets(c("v" = "Merging is done using the common SNPs between the input genotype files.")))

                # Clean-up
                ftemp <- c(list.files(normalizePath(ResultDir, mustWork = FALSE), pattern = "new"), list.files(normalizePath(ResultDir, mustWork = FALSE), pattern = "common_snps"))
                invisible(file.remove(normalizePath(file.path(ResultDir, ftemp), mustWork = FALSE)))
            } else {
                merge_args <- c(
                    "--bfile", normalizePath(file.path(DataDir, finput1), mustWork = FALSE),
                    "--bmerge", normalizePath(file.path(DataDir, finput2), mustWork = FALSE),
                    "--allow-no-sex", "--make-bed",
                    "--out", normalizePath(file.path(ResultDir, foutput), mustWork = FALSE),
                    "--silent"
                )
                executePlink(merge_args)
                rlang::inform(rlang::format_error_bullets(c("v" = "Merging is done with all the SNPs i.e., union of the SNPs.")))
            }

            rlang::inform(rlang::format_error_bullets(c("v" = paste0("Plink files with merged regions are in ", ResultDir, " prefixed as ", foutput))))
        },
        error = function(e) {
            rlang::abort(
                message = e$message,
                class = "MergeRegion_error"
            )
        },
        warning = function(w) {
            rlang::warn(
                message = w$message,
                .frequency = "regularly",
                .frequency_id = "MergeRegion_warning"
            )
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
#' data("Ffile", package = "GXwasR")
#' data("Mfile", package = "GXwasR")
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
                highlight_p = c(snp_pval, snp_pval), highlighter = "green", chrblocks = TRUE, file = normalizePath(file.path(ResultDir, "Stratified_GWAS"), mustWork = FALSE)
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
                    highlight_p = c(snp_pval, snp_pval), highlighter = "green", chrblocks = TRUE, file = normalizePath(file.path(ResultDir, "Stratified_XWAS"), mustWork = FALSE)
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
#' Models such as "FMcombx01","FMcombx02",and "FMstratified" can be applied to both binary and quantitative traits,
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
#' Models "FMcombx01","FMcombx02",and "FMstratified" can be chosen for both binary and quantitative traits
#' while "GWAcxci" can only apply to the binary trait. These models take care of the X-chromosomal marker.
#' Three female genotypes are coded by 0, 1, and 2 in FM01 and FM02. The two genotypes of males that follow the
#' X-chromosome inactivation (XCI) pattern as random (XCI-R) in the FM01 model are coded by 0 and 1, while the two
#' genotypes that follow the XCI is escaped (XCI-E) in the FM02 model are coded by 0 and 1. To reflect the dose
#' compensation connection between the sexes, FM02 treats men as homozygous females.
#'
#' In the "FMstratified" associations are tested separately for males and females, and then the combined p values are computed the Fisher's method, Fisher's method with permutation,
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
#' be excluded from association analysis). It is important to note that for stratified GWAS model, if PCs are included
#' as covar then it should be generated separately for each cohort and then included in the covarfile. Use the function
#' \code{\link{DummyCovar}} to generate a new covariate file with categorical variables down-coded as binary dummy variables for
#' the covariate file with categorical variables. For instance, if a variable has K categories, K-1 new dummy variables
#' are constructed, and the original covariate is now estimated with a coefficient for each category.
#'
#' @param covartest
#' Vector value with `NULL`,"ALL" or covariate name/names to be included in the test. The default is `NULL.` For instance,
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
#' Character vector specifying method for combining p-values after stratified GWAS/XWAS models.
#' Choices are “stouffer.method”, "fisher.method" and "fisher.method.perm". For fisher.method the function for combining
#' p-values uses a statistic, \eqn{S = -2 \sum_{i=1}^{k} \log(p_i)}, which follows a \eqn{\chi^2} distribution with 2k degrees of freedom \insertCite{Fisher1925}{GXwasR}.
#'
#' For fisher.method.perm, using p-values from stratified tests, the summary statistic for combining p-values is \eqn{S = -2 \sum_{i=1}^{k} \log(p_i)}.
#' A p-value for this statistic can be derived by randomly generating summary statistics \insertCite{Rhodes2002}{GXwasR}. Therefore, a p-value is randomly
#' sampled from each contributing study, and a random statistic is calculated. The fraction of random statistics greater or
#' equal to S then gives the final p-value.
#'
#' For stouffer.method ,the function applies Stouffer’s method \insertCite{Stouffer1949}{GXwasR} to the p-values assuming that the p-values to be combined are
#' independent. Letting p1, p2, . . . , pk denote the individual (one- or two-sided) p-values of the k hypothesis tests to be
#' combined, the test statistic is then computed as \eqn{z = \frac{\sum_{i=1}^{k} z_i}{\sqrt{k}}}, where \eqn{z_i = \Phi^{-1}(1 - p_i)} and
#' \eqn{\Phi^{-1}(\cdot)} denotes the inverse of the cumulative distribution function of a standard normal distribution. Under the joint null
#' hypothesis, the test statistic follows a standard normal distribution which is used to compute the combined p-value. This
#' functionality is taken from the R package poolr \insertCite{Cinar2022}{GXwasR}.
#'
#' Note that only p-values between 0 and 1 are allowed to be passed to these methods.
#'
#' Note: Though this parameter is enabled for both autosome GWAS and XWAS, the combining pvalue after
#' sex-stratified test is recommended to ChrX only.
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
#' ncores <- 2
#' MF.mc.cores <- 1
#' ResultGXwas <- GXwas(
#'     DataDir = DataDir, ResultDir = ResultDir,
#'     finput = finput, xmodel = xmodel, trait = trait, covarfile = covarfile,
#'     sex = sex, xsex = xsex, combtest = combtest, MF.p.corr = "none",
#'     snp_pval = snp_pval, plot.jpeg = TRUE, suggestiveline = 5, genomewideline = 7.3,
#'     MF.mc.cores = 1, ncores = ncores
#' )
GXwas <- function(
      DataDir, ResultDir, finput, trait = c("binary", "quantitative"), standard_beta = TRUE,
      xmodel = c("FMcombx01", "FMcombx02", "FMstratified", "GWAScxci"), sex = FALSE, xsex = FALSE,
      covarfile = NULL, interaction = FALSE, covartest = c("ALL"), Inphenocov = c("ALL"), combtest = c("fisher.method", "fisher.method.perm", "stouffer.method"),
      MF.zero.sub = 0.00001, B = 10000, MF.mc.cores = 1, MF.na.rm = FALSE,
      MF.p.corr = "none", plot.jpeg = FALSE, plotname = "GXwas.plot", snp_pval = 1e-08,
      annotateTopSnp = FALSE, suggestiveline = 5, genomewideline = 7.3, ncores = 0
) {
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
            } else if (xmodel[1] == "FMstratified") {
                rlang::inform(rlang::format_error_bullets("Running FMstratified model"))

                ## Making male and female files in ResultDir

                MFsplitPlink(DataDir = DataDir, ResultDir = ResultDir, finput = finput, foutput = "finput.female", sex = "females")
                pb$tick(10)
                gc(reset = TRUE)
                MFsplitPlink(DataDir = DataDir, ResultDir = ResultDir, finput = finput, foutput = "finput.male", sex = "males")
                gc(reset = TRUE)
                pb$tick(15)
                x <-
                    suppressWarnings(
                        FMcomb(
                            DataDir = DataDir, ResultDir = ResultDir, trait = trait, standard_beta = standard_beta, xmodel = xmodel,
                            covarfile = covarfile, interaction = interaction, covartest = covartest, Inphenocov = Inphenocov,
                            plot.jpeg = plot.jpeg, plotname = plotname, snp_pval = snp_pval, annotateTopSnp = annotateTopSnp,
                            combtest = combtest, B = B, MF.p.corr = MF.p.corr, MF.zero.sub = MF.zero.sub, MF.na.rm = MF.na.rm,
                            MF.mc.cores = MF.mc.cores, suggestiveline = suggestiveline, genomewideline = genomewideline, ncores = ncores
                        )
                    )
                pb$tick(20)
                ftemp <- list.files(normalizePath(file.path(ResultDir), mustWork = FALSE), pattern = "finput")
                invisible(file.remove(normalizePath(file.path(ResultDir, ftemp), mustWork = FALSE)))
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

            patterns <- c("_ss", ".logistic", "_snps", "PLINK", "allsnpsresults.rda", "all_snps_results")
            removePatternFiles(ResultDir = ResultDir, patterns = patterns)
            pb$tick(100)
            return(x)
        },
        error = function(e) {
            rlang::abort(
                message = e$message,
                class = "GXwasR_error"
            )
        },
        warning = function(w) {
            rlang::warn(
                message = w$message,
                .frequency = "regularly",
                .frequency_id = "GXwasR_warning"
            )
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
#' columns "INDEX"(Index SNP identifier), "PSNP"(Best proxy SNP), "RSQ LD"(r-squared) between index and proxy,
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
#' data("Summary_Stat_Ex1", package = "GXwasR")
#' data("Summary_Stat_Ex2", package = "GXwasR")
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' SNPdata <- list(Summary_Stat_Ex1, Summary_Stat_Ex2)
#' clump_p1 <- 0.0001
#' clump_p2 <- 0.001
#' clump_r2 <- 0.5
#' clump_kb <- 250
#' byCHR <- TRUE
#' clumpedResult <- ClumpLD(
#'     DataDir, finput, SNPdata, ResultDir, clump_p1,
#'     clump_p2, clump_r2, clump_kb, byCHR
#' )
ClumpLD <- function(
      DataDir, finput, SNPdata, ResultDir = tempdir(),
      clump_p1, clump_p2, clump_r2, clump_kb, byCHR = TRUE,
      clump_best = TRUE, clump_index_first = TRUE
) {
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
            for (i in seq_along(SNPdata)) {
                rlang::inform(rlang::format_error_bullets(paste0("Processing summary statistics ", i)))
                write.table(SNPdata[[i]],
                    normalizePath(file.path(ResultDir, paste0("SNPdata_", i)), mustWork = FALSE),
                    row.names = FALSE, col.names = TRUE, quote = FALSE
                )
            }

            SumData <- list.files(normalizePath(file.path(ResultDir), mustWork = FALSE), pattern = "SNPdata_")

            if (byCHR == TRUE) {
                bimfile <- read.table(normalizePath(file.path(DataDir, paste0(finput, ".bim")), mustWork = FALSE))
                chrnum <- seq_len(length(unique(bimfile$V1)))

                chrwiseLD <- function(chrnum) {
                    chromosome <- unique(bimfile$V1)[chrnum]
                    rlang::inform(rlang::format_error_bullets(paste0("Running LD clumping for chromosome ", chromosome)))

                    # Construct the PLINK command arguments with toggles
                    plink_args <- c(
                        "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE),
                        "--chr", as.character(chromosome),
                        "--clump", normalizePath(file.path(ResultDir, SumData), mustWork = FALSE),
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
                        "--out", normalizePath(file.path(ResultDir, "ClumpLD"), mustWork = FALSE)
                    )

                    # Execute the PLINK command
                    invisible(sys::exec_wait(
                        plink(),
                        args = plink_args,
                        std_out = FALSE, # Do not capture standard output
                        std_err = FALSE # Do not capture standard error
                    ))

                    if (file.exists(normalizePath(file.path(ResultDir, "ClumpLD.clumped"), mustWork = FALSE))) {
                        # If both flags are ON, try reading the produced best file.
                        # Otherwise, create a dummy BestClump dataframe.
                        if (clump_best && clump_index_first &&
                            file.exists(normalizePath(file.path(ResultDir, "ClumpLD.clumped.best"), mustWork = FALSE))) {
                            ldc2 <- as.data.frame(data.table::fread(
                                normalizePath(file.path(ResultDir, "ClumpLD.clumped.best"), mustWork = FALSE),
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
                        save(ldcall, file = normalizePath(file.path(ResultDir, paste0("ldcall_", chromosome, ".Rda")), mustWork = FALSE))
                        invisible(do.call(file.remove, list(normalizePath(file.path(ResultDir, "ClumpLD.clumped"), mustWork = FALSE))))
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

                ldfiles <- list.files(normalizePath(file.path(ResultDir), mustWork = FALSE), pattern = "ldcall_")
                ldf <- function(ldfiles) {
                    load(normalizePath(file.path(ResultDir, ldfiles), mustWork = FALSE))
                    return(ldcall)
                }
                All_ldc <- data.table::rbindlist(lapply(ldfiles, ldf), fill = TRUE)
            } else {
                # Construct the PLINK command arguments with toggles
                plink_args <- c(
                    "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE),
                    "--clump", normalizePath(file.path(ResultDir, SumData), mustWork = FALSE),
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
                    "--out", normalizePath(file.path(ResultDir, "ClumpLD"), mustWork = FALSE)
                )

                # Execute the PLINK command
                invisible(sys::exec_wait(
                    plink(),
                    args = plink_args,
                    std_out = FALSE, # Do not capture standard output
                    std_err = FALSE # Do not capture standard error
                ))

                if (file.exists(normalizePath(file.path(ResultDir, "ClumpLD.clumped.best"), mustWork = FALSE))) {
                    ldc2 <- read.table(normalizePath(file.path(ResultDir, "ClumpLD.clumped.best"), mustWork = FALSE), header = TRUE)
                    ldc2 <- ldc2[!apply(ldc2 == "", 1, all), ]
                    ldcall <- suppressWarnings(read_plink_clumped_clean(
                        resultDir = ResultDir,
                        filename = "ClumpLD.clumped"
                    ))
                    save(ldcall, file = normalizePath(file.path(ResultDir, "ldcall_genome.Rda"), mustWork = FALSE))
                    invisible(do.call(file.remove, list(normalizePath(file.path(ResultDir, "ClumpLD.clumped"), mustWork = FALSE))))

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
            rlang::abort(
                message = e$message,
                class = "ClumpLD_error"
            )
        },
        warning = function(w) {
            rlang::warn(
                message = w$message,
                .frequency = "regularly",
                .frequency_id = "ClumpLD_warning"
            )
        }
    )

    ftemp <- list.files(normalizePath(file.path(ResultDir), mustWork = FALSE), pattern = "SNPdata_")
    invisible(do.call(file.remove, list(normalizePath(file.path(ResultDir, ftemp), mustWork = FALSE))))
    ftemp <- list.files(normalizePath(file.path(ResultDir), mustWork = FALSE), pattern = "Clump")
    invisible(do.call(file.remove, list(normalizePath(file.path(ResultDir, ftemp), mustWork = FALSE))))
    ftemp <- list.files(normalizePath(file.path(ResultDir), mustWork = FALSE), pattern = "ldcall")
    invisible(do.call(file.remove, list(normalizePath(file.path(ResultDir, ftemp), mustWork = FALSE))))

    return(list(BestClump = LDC, AllClump = All_ldc))
}


#' SexRegress: Performing linear regression analysis with quantitative response variable.
#'
#' @description
#' This function could be used to check association of two variables. For instance, PGS with sex.
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
#' data("Regression_Ex", package = "GXwasR")
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
            model.result <- summary(model)$coefficients[regressor_index, ]
            return(model.result)
        },
        error = function(e) {
            rlang::abort(
                message = e$message,
                class = "SexRegress_error"
            )
        },
        warning = function(w) {
            rlang::warn(
                message = w$message,
                .frequency = "regularly",
                .frequency_id = "SexRegress_warning"
            )
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
#' `NULL`. After multi-allelic variant filtering, the filtered PLINK files with only biallelic SNPs will be saved in `ResultDir`.
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
            bimf <- read.table(normalizePath(file.path(DataDir, paste0(finput, ".bim")), mustWork = FALSE))
            x1 <- bimf[nchar(bimf[, 5]) > 1 | nchar(bimf[, 6]) > 1, , drop = FALSE]

            if (nrow(x1) != 0) {
                write.table(x1$V2, file = normalizePath(file.path(ResultDir, "snps_multiallelic"), mustWork = FALSE), quote = FALSE, col.names = FALSE, row.names = FALSE)
            } else {
                rlang::inform(
                    rlang::format_error_bullets(c("i" = "There are no multi-allelic SNPs present in the input dataset."))
                )
            }
            exclude_arg <- if (nrow(x1) > 0) normalizePath(file.path(ResultDir, "snps_multiallelic"), mustWork = FALSE) else NULL
            invisible(sys::exec_wait(
                plink(),
                args = c(
                    "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE),
                    "--exclude", exclude_arg,
                    "--allow-no-sex", # 4.0
                    "--make-bed",
                    "--out", normalizePath(file.path(ResultDir, foutput), mustWork = FALSE),
                    "--silent"
                ),
                std_out = FALSE,
                std_err = FALSE
            ))
            if (nrow(x1) > 0) {
                bimf1 <- read.table(normalizePath(file.path(ResultDir, paste0(foutput, ".bim")), mustWork = FALSE))

                rlang::inform(
                    rlang::format_error_bullets(c(
                        "i" = paste0("Input dataset has ", nrow(bimf), " SNPs."),
                        "i" = paste0("Output dataset has ", nrow(bimf1), " SNPs."),
                        "v" = paste0("Plink files with only biallelic SNPs are in ", ResultDir, " prefixed as ", foutput)
                    ))
                )
                return(invisible(bimf1))
            }
        },
        error = function(e) {
            rlang::abort(
                message = e$message,
                class = "FilterAllele_error"
            )
        },
        warning = function(w) {
            rlang::warn(
                message = w$message,
                .frequency = "regularly",
                .frequency_id = "FilterAllele_warning"
            )
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
#' Other columns may present.
#'
#' @param SumstatFemale
#' R dataframe object of summary statistics of female GWAS with five mandatory columns:
#' * `SNP`
#' * `A1`
#' * `TEST`
#' * `POS`
#' * `P`
#'
#' Other columns may present.
#'
#' @param combtest
#' Character vector specifying method for combining p-values for stratified GWAS models. Choices are “stouffer.method”,
#' "fisher.method" and "fisher.method.perm". For fisher.method, the function for combining p-values uses a statistic,
#' \eqn{S = -2 \sum_{i=1}^{k} \log(p_i)}, which follows a \eqn{\chi^2} distribution with 2k degrees of freedom \insertCite{Fisher1925}{GXwasR}.
#' For fisher.method.perm, using p-values from stratified tests, the summary statistic for combining p-values
#' is \eqn{S = -2 \sum_{i=1}^{k} \log(p_i)}. A p-value for this statistic can be derived by randomly generating summary statistics \insertCite{Rhodes2002}{GXwasR}.
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
#' data("Mfile", package = "GXwasR")
#' data("Ffile", package = "GXwasR")
#' SumstatMale <- Mfile
#' colnames(SumstatMale)[3] <- "POS"
#' SumstatFemale <- Ffile
#' colnames(SumstatFemale)[3] <- "POS"
#' PvalComb_Result <- PvalComb(
#'     SumstatMale = SumstatMale, SumstatFemale = SumstatFemale,
#'     combtest = "fisher.method", MF.mc.cores = 1, snp_pval = 0.001, plot.jpeg = FALSE,
#'     suggestiveline = 3, genomewideline = 5.69897, ncores = 1
#' )
#'
PvalComb <- function(
      SumstatMale, SumstatFemale,
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
      ncores = 0
) {
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
            rlang::abort(message = e$message, class = "PvalComb_error")
        },
        warning = function(w) {
            rlang::warn(
                message = w$message,
                .frequency = "regularly",
                .frequency_id = "PvalComb_warning"
            )
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
            if (extract == TRUE) {
                remov <- "--extract"
            } else if (extract == FALSE) {
                remov <- "--exclude"
            }
            write.table(SNPvec, file = normalizePath(file.path(ResultDir, "snplist"), mustWork = FALSE), sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)

            # Exclude SNPs
            invisible(sys::exec_wait(
                plink(),
                args = c(
                    "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE),
                    remov, normalizePath(file.path(ResultDir, "snplist"), mustWork = FALSE),
                    "--allow-no-sex", # 4.0
                    "--make-bed",
                    "--out", normalizePath(file.path(ResultDir, foutput), mustWork = FALSE),
                    "--silent"
                ),
                std_out = FALSE,
                std_err = FALSE
            ))

            bim <- read.table(normalizePath(file.path(ResultDir, paste0(foutput, ".bim")), mustWork = FALSE))

            rlang::inform(
                rlang::format_error_bullets(c(
                    "i" = paste0(nrow(bim), " SNPs are extracted"),
                    "v" = paste0("Plink files with extracted SNPs are in ", ResultDir, " prefixed as ", foutput)
                ))
            )
            return(invisible(NULL))
        },
        error = function(e) {
            rlang::abort(
                message = e$message,
                class = "FilterSNP_error"
            )
        },
        warning = function(w) {
            rlang::warn(
                message = w$message,
                .frequency = "regularly",
                .frequency_id = "FilterSNP_warning"
            )
        }
    )
}


#' Validate Path to Reference Data Set
#'
#' @description
#' Validates that reference data for either 'HapMapIII_NCBI36' or 'ThousandGenome' is present
#' in the path specified by the appropriate environment variable:
#' - 'HAPMAPIII_NCBI36_DIR' for HapMapIII_NCBI36
#' - 'THOUSANDGENOME_DIR' for ThousandGenome
#'
#' If files are missing or paths are not set, informative guidance is provided.
#'
#' @param refdata
#' A character string specifying the reference dataset to validate. Should be one of
#' 'HapMapIII_NCBI36' or 'ThousandGenome'.
#'
#' @return
#' A normalized path to the directory containing the validated reference data files.
#'
#' @export
#' @examples
#' if (nzchar(Sys.getenv("HAPMAPIII_NCBI36_DIR"))) {
#'     validate_reference_data("HapMapIII_NCBI36")
#' }
#'
#' if (nzchar(Sys.getenv("HAPMAPIII_NCBI36_DIR"))) {
#'     validate_reference_data("ThousandGenome")
#' }
validate_reference_data <- function(refdata) {
    valid_refs <- c("HapMapIII_NCBI36", "ThousandGenome", "Ref10Kgenome")
    if (!is.character(refdata) || length(refdata) != 1 || !(refdata %in% valid_refs)) {
        stop("Invalid 'refdata'. Must be one of: ", paste(valid_refs, collapse = ", "))
    }

    # Define environment variable and expected files based on dataset
    env_var <- switch(refdata,
        "HapMapIII_NCBI36" = "HAPMAPIII_NCBI36_DIR",
        "ThousandGenome" = "THOUSANDGENOME_DIR",
        "Ref10Kgenome" = "THOUSANDGENOME_DIR"
    )

    expected_files <- switch(refdata,
        "HapMapIII_NCBI36" = c("HapMapIII_NCBI36.bed", "HapMapIII_NCBI36.bim", "HapMapIII_NCBI36.fam"),
        "ThousandGenome" = c("Ref10Kgenome.bed", "Ref10Kgenome.bim", "Ref10Kgenome.fam"),
        "Ref10Kgenome" = c("Ref10Kgenome.bed", "Ref10Kgenome.bim", "Ref10Kgenome.fam")
    )

    dir_path <- Sys.getenv(env_var, unset = NA)
    if (is.na(dir_path) || !dir.exists(dir_path)) {
        stop(
            "Environment variable '", env_var, "' is not set or points to an invalid directory.\n",
            "Set it using Sys.setenv(", env_var, " = '/path/to/your/data')"
        )
    }

    # Check presence of all expected files
    full_paths <- file.path(dir_path, expected_files)
    missing_files <- expected_files[!file.exists(full_paths)]

    if (length(missing_files) > 0) {
        stop(
            "Missing files for '", refdata, "' in ", normalizePath(dir_path), ":\n",
            paste("-", missing_files, collapse = "\n"), "\n\n",
            "Please download and extract the reference data manually.\n",
            "Download URLs:\n",
            "- HapMapIII_NCBI36: https://figshare.com/ndownloader/files/40585145\n",
            "- ThousandGenome: https://figshare.com/ndownloader/files/46552177"
        )
    }

    rlang::inform(
        rlang::format_error_bullets(c(
            "i" = paste0("'", refdata, "' reference data found at ", normalizePath(dir_path), ".")
        )),
        .frequency = "regularly", .frequency_id = "validate_reference"
    )

    return(normalizePath(dir_path))
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
#' @param ResultDir
#' A character string for the file path where all output files will be stored. The default is `tempdir()`.
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
#' dummy_covars <- DummyCovar(
#'     DataDir = DataDir, ResultDir = ResultDir,
#'     bfile = bfile, incovar = incovar,
#'     outcovar = outcovar
#' )
DummyCovar <- function(DataDir, ResultDir = DataDir, bfile, incovar, outcovar) {
    if (
        file.exists(normalizePath(file.path(DataDir, paste0(bfile, ".bed")), mustWork = FALSE)) &&
            file.exists(normalizePath(file.path(DataDir, paste0(bfile, ".bim")), mustWork = FALSE)) &&
            file.exists(normalizePath(file.path(DataDir, paste0(bfile, ".fam")), mustWork = FALSE))
    ) {
        invisible(
            sys::exec_wait(
                plink(),
                args = c(
                    "--bed", normalizePath(file.path(DataDir, paste0(bfile, ".bed")), mustWork = FALSE),
                    "--bim", normalizePath(file.path(DataDir, paste0(bfile, ".bim")), mustWork = FALSE),
                    "--fam", normalizePath(file.path(DataDir, paste0(bfile, ".fam")), mustWork = FALSE),
                    "--covar",
                    normalizePath(file.path(DataDir, incovar), mustWork = FALSE),
                    "--write-covar",
                    "--dummy-coding",
                    "--out", normalizePath(file.path(ResultDir, outcovar), mustWork = FALSE),
                    "--silent"
                ),
                std_out = FALSE,
                std_err = FALSE
            )
        )
    } else {
        stop("There are no PLINK files in DataDir.\nPlease specify correct directory path with input PLINK files.")
    }

    if (file.exists(normalizePath(file.path(ResultDir, paste0(outcovar, ".cov")), mustWork = FALSE))) {
        x <- read.table(normalizePath(file.path(ResultDir, paste0(outcovar, ".cov")), mustWork = FALSE), header = TRUE)
        rlang::inform(rlang::format_error_bullets(c("i" = paste0("Covariate file: ", outcovar, ".cov is in ", ResultDir))))
        return(x)
    }
}


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
        "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE),
        "--indep-pairwise",
        window_size,
        step_size,
        r2_threshold,
        "--out", normalizePath(file.path(ResultDir, "ld_prune"), mustWork = FALSE),
        "--silent"
    )

    # Execute the PLINK command using the helper function within tryCatch
    tryCatch(
        {
            executePlinkAd(ResultDir, args)

            # Load the list of pruned SNPs if available
            prune_in_file <- normalizePath(file.path(ResultDir, "ld_prune.prune.in"), mustWork = FALSE)
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
#' This function requires access to the \href{https://zenodo.org/records/16923484}{reference LD data} via an
#' environment variable. You must set one of the following environment variables to the appropriate directory:
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
#'     res <- SumstatGenCorr(
#'         ResultDir = tempdir(),
#'         referenceLD = "UKB_imputed_hapmap2_SVD_eigen99_extraction",
#'         sumstat1 = sumstat1,
#'         sumstat2 = sumstat2,
#'         parallel = TRUE
#'     )
#' }
SumstatGenCorr <- function(ResultDir = tempdir(),
    referenceLD,
    sumstat1,
    sumstat2,
    Nref = 335265,
    N0 = min(sumstat1$N),
    eigen.cut = "automatic",
    lim = exp(-18),
    parallel = FALSE,
    numCores = 2) {
    reference_paths <- list(
        UKB_imputed_hapmap2_SVD_eigen99_extraction = Sys.getenv("UKB_IMPUTED_HAPMAP2_PATH", unset = NA),
        UKB_imputed_SVD_eigen99_extraction = Sys.getenv("UKB_IMPUTED_PATH", unset = NA),
        UKB_array_SVD_eigen90_extraction = Sys.getenv("UKB_ARRAY_PATH", unset = NA)
    )

    reference_urls <- c(
        UKB_imputed_hapmap2_SVD_eigen99_extraction = "https://zenodo.org/records/16923484/files/UKB_imputed_hapmap2_SVD_eigen99_extraction.tar.gz?download=1",
        UKB_imputed_SVD_eigen99_extraction = "https://zenodo.org/records/16923484/files/UKB_imputed_hapmap2_SVD_eigen99_extraction.tar.gz?download=1",
        UKB_array_SVD_eigen90_extraction = "https://zenodo.org/records/16923484/files/UKB_array_SVD_eigen90_extraction.tar.gz?download=1"
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
                HDL.rg(
                    gwas1.df = sumstat1, gwas2.df = sumstat2,
                    LD.path = LD_path, Nref = Nref, N0 = N0,
                    eigen.cut = eigen.cut, lim = lim
                )
            } else {
                HDL.rg.parallel(
                    gwas1.df = sumstat1, gwas2.df = sumstat2,
                    LD.path = LD_path, Nref = Nref, N0 = N0,
                    eigen.cut = eigen.cut, lim = lim, numCores = numCores
                )
            }
        },
        error = function(e) {
            rlang::abort(
                message = glue::glue("Error in estimating Genetic Correlation: {e$message}"),
                class = "SumstatGenCorr_error"
            )
        }
    )

    return(res.HDL)
}

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
#' snpld <- ComputeLD(
#'     DataDir = system.file("extdata", package = "GXwasR"), ResultDir = tempdir(),
#'     finput = "GXwasR_example", ByCHR = TRUE, CHRnum = 1, r2_LD = 0.2
#' )
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
            "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE),
            chr, CHRnum,
            "--r2",
            "--ld-window-r2", r2_LD,
            "--out", normalizePath(file.path(ResultDir, "snpcorr"), mustWork = FALSE),
            "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
    ))
    snpld <- read.table(normalizePath(file.path(ResultDir, "snpcorr.ld"), mustWork = FALSE), header = TRUE)
    return(snpld)
}
