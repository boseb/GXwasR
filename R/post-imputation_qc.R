## Function 136
## Added in 3.0
validateInputForXhwe <- function(DataDir, ResultDir = tempdir(), finput, foutput, filterSNP = TRUE) {
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

    # Validate filterSNP
    if (!is.logical(filterSNP)) {
        stop("Error in filterSNP: Must be a boolean value.")
    }

    return(TRUE)
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
#' DataDir <- GXwasR:::GXwasR_data()
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' foutput <- "Test_output"
#' x <- Xhwe(
#'     DataDir = DataDir, ResultDir = ResultDir,
#'     finput = finput, foutput = foutput, filterSNP = TRUE
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
            ## Getting female Plink file
            invisible(GetMFPlink(DataDir = DataDir, ResultDir, finput = finput, foutput = "female", sex = "females", xplink = FALSE, autoplink = FALSE))

            fam <-
                as.data.frame(utils::read.table(file = normalizePath(file.path(ResultDir, "female.fam"), mustWork = FALSE)))

            # Check for case control status in input PLINK file
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
                    "--bfile", normalizePath(file.path(ResultDir, "female"), mustWork = FALSE),
                    "--chr", 23,
                    "--hardy",
                    "--out", normalizePath(file.path(ResultDir, "Xhwe"), mustWork = FALSE),
                    "--silent"
                ),
                std_out = FALSE,
                std_err = FALSE
            ))


            x <-
                as.data.frame(
                    read.table(
                        file = normalizePath(file.path(ResultDir, "Xhwe.hwe"), mustWork = FALSE),
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
                        "i" = "Input PLINK files are unchanged. No output PLINK files are produced."
                    ))
                )
                return()
            } else {
                if (length(X_excluded_SNPs) != 0 & filterSNP == "TRUE") {
                    utils::write.table(
                        X_excluded_SNPs,
                        file = normalizePath(file.path(ResultDir, "XhweSNPs"), mustWork = FALSE),
                        quote = FALSE,
                        row.names = FALSE,
                        col.names = FALSE,
                        eol = "\r\n",
                        sep = " "
                    )


                    invisible(sys::exec_wait(
                        plink(),
                        args = c(
                            "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE),
                            "--exclude", normalizePath(file.path(ResultDir, "XhweSNPs"), mustWork = FALSE),
                            "--allow-no-sex", ## 4.0
                            "--make-bed",
                            "--out", normalizePath(file.path(ResultDir, foutput), mustWork = FALSE),
                            "--silent"
                        ),
                        std_out = FALSE,
                        std_err = FALSE
                    ))


                    rlang::inform(rlang::format_error_bullets(c("i" = paste0("Failed SNPs are excluded from the output PLINK files prefixed as ", foutput, " is in ", ResultDir))))

                    ftemp <- list.files(normalizePath(file.path(ResultDir)), pattern = "hwe")
                    invisible(file.remove(normalizePath(file.path(ResultDir, ftemp))))
                    ftemp <- list.files(normalizePath(file.path(ResultDir)), pattern = "female")
                    invisible(file.remove(normalizePath(file.path(ResultDir, ftemp))))
                    return(X_excluded_SNPs)
                } else {
                    rlang::inform(rlang::format_error_bullets(c("i" = "SNPs are flagged.")))
                    ftemp <- list.files(normalizePath(file.path(ResultDir)), pattern = "hwe")
                    invisible(file.remove(normalizePath(file.path(ResultDir, ftemp))))
                    ftemp <- list.files(normalizePath(file.path(ResultDir)), pattern = "female")
                    invisible(file.remove(normalizePath(file.path(ResultDir, ftemp))))
                    return(X_excluded_SNPs)
                }
            }
        },
        error = function(e) {
            rlang::abort(
                message = e$message, 
                class = "Xhwe_error"
            )
        },
        warning = function(w) {
            rlang::warn(
                message = w$message, 
                class = "Xhwe_warning", 
                .frequency = "regularly", 
                .frequency_id = "Xhwe_warning"
            )
        }
    )
}

## Function 137
## Added in 3.0
validateInputForMAFdiffSexControl <- function(DataDir, ResultDir = tempdir(), finput, filterSNP = FALSE, foutput = NULL) {
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

    # Validate foutput if not NULL
    if (!is.null(foutput) && !is.character(foutput)) {
        stop("Error in foutput: Must be a character string.")
    }

    # Validate filterSNP
    if (!is.logical(filterSNP)) {
        stop("Error in filterSNP: Must be a boolean value.")
    }

    return(TRUE)
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
#' @importFrom utils read.table write.table
#' @importFrom sys exec_wait
#' @export
#'
#' @examples
#' DataDir <- GXwasR:::GXwasR_data()
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' foutput <- "Test_output"
#' x <- MAFdiffSexControl(DataDir, ResultDir, finput, filterSNP = TRUE, foutput = foutput)
MAFdiffSexControl <- function(
      DataDir,
      ResultDir = tempdir(),
      finput,
      filterSNP = FALSE,
      foutput = NULL
) {
    if (!validateInputForMAFdiffSexControl(DataDir, ResultDir, finput, filterSNP, foutput)) {
        return(NULL)
    }

    if (!checkFiles(DataDir, finput)) {
        stop("Missing required Plink files in the specified DataDir.")
    }

    tryCatch(
        {
            fam <-
                as.data.frame(utils::read.table(file = normalizePath(file.path(DataDir, paste0(finput, ".fam")), mustWork = FALSE)))

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
                file = normalizePath(file.path(ResultDir, "phenofile"), mustWork = FALSE),
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE,
                eol = "\r\n",
                sep = " "
            )

            invisible(sys::exec_wait(
                plink(),
                args = c(
                    "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE),
                    "--logistic",
                    "--pheno", normalizePath(file.path(ResultDir, "phenofile"), mustWork = FALSE),
                    "--out", normalizePath(file.path(ResultDir, "OUTPUT"), mustWork = FALSE),
                    "--silent"
                ),
                std_out = FALSE,
                std_err = FALSE
            ))

            x <-
                utils::read.table(
                    file = normalizePath(file.path(ResultDir, "OUTPUT.assoc.logistic"), mustWork = FALSE),
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
                    file = normalizePath(file.path(ResultDir, "flaggedSnpsSexMafDiff"), mustWork = FALSE),
                    quote = FALSE,
                    row.names = FALSE,
                    col.names = FALSE,
                    eol = "\r\n",
                    sep = " "
                )


                invisible(sys::exec_wait(
                    plink(),
                    args = c(
                        "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE),
                        "--exclude", normalizePath(file.path(ResultDir, "flaggedSnpsSexMafDiff"), mustWork = FALSE),
                        "--allow-no-sex", # 4.0
                        "--make-bed",
                        "--out", normalizePath(file.path(ResultDir, foutput), mustWork = FALSE),
                        "--silent"
                    ),
                    std_out = FALSE,
                    std_err = FALSE
                ))


                writeLines(
                    paste0("SNPs with significantly MAF difference are excluded.\nFiltered PLINK files are saved in ", ResultDir)
                )
                return(as.list(flaggedSnps))
            } else if (length(flaggedSnps) != 0 & filterSNP == FALSE) {
                rlang::inform(rlang::format_error_bullets(c("i" = "SNPs are flagged.")))
            }

            ftemp <- list.files(normalizePath(file.path(ResultDir), mustWork = FALSE), pattern = "OUTPUT")
            invisible(do.call(file.remove, list(normalizePath(file.path(ResultDir, ftemp), mustWork = FALSE))))
            invisible(do.call(file.remove, list(normalizePath(file.path(ResultDir, "phenofile"), mustWork = FALSE))))

            return(flaggedSnps)
        },
        error = function(e) {
            rlang::abort(
                message = e$message, 
                class = "MAFdiffSexControl_error")
        },
        warning = function(w) {
            rlang::warn(
                message = w$message, 
                class = "MAFdiffSexControl_warning", 
                .frequency = "regularly", 
                .frequency_id = "MAFdiffSexControl_warning"
            )
        }
    )
}

#' FilterRegion: Filter chromosomal regions.
#'
#' @author Banabithi Bose
#'
#' @description
#' Filtering Pseudo-Autosomal Region (PAR), X-transposed region (XTR), Ampliconic, filter based on chromosome code or
#' user-defined regions from input PLINK files. Only one type of filtering can be done from three types, either by region
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
#' Boolean value, `TRUE` or `FALSE` to filter out PARs from input PLINK file. The default is `TRUE`.
#'
#' @param filterXTR
#' Boolean value, `TRUE` or `FALSE` to filter out XTRs from input PLINK file. The default is `TRUE`.
#'
#' @param filterAmpliconic
#' Boolean value, `TRUE` or `FALSE` to filter out Ampliconic regions from input PLINK file. The default is `TRUE`.
#'
#' @param regionfile
#' Character string, specifying the name of the .txt file containing the user-defined regions to be filtered out from input PLINK
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
#' For `exclude` = `TRUE`, two sets of PLINK binary files will be produced in ResultDir. One set will have the remaining SNPs after
#' filtering and other one will have the discarded SNPs.
#'
#' @export
#'
#' @examples
#' DataDir <- GXwasR:::GXwasR_data()
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' foutput <- "PostimputeEX_QC1"
#' x <- FilterRegion(
#'     DataDir = DataDir, ResultDir = ResultDir,
#'     finput = finput, foutput = foutput, CHRX = TRUE, CHRY = FALSE,
#'     filterPAR = TRUE, filterXTR = TRUE, filterAmpliconic = TRUE,
#'     regionfile = FALSE, filterCHR = NULL, Hg = "38", exclude = TRUE
#' )
FilterRegion <-
    function(
      DataDir,
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
      exclude = TRUE
    ) {
        # Validate parameters
        validateFilterRegionParams(DataDir, ResultDir, finput, foutput, CHRX, CHRY, filterPAR, filterXTR, filterAmpliconic, regionfile, filterCHR, Hg, exclude)

        if (!checkFiles(DataDir, finput)) {
            stop("Missing required Plink files in the specified DataDir.")
        }

        tryCatch(
            {
                DataDir1 <- GXwasR_data()

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
                                    file = normalizePath(file.path(ResultDir, "region.txt"), mustWork = FALSE),
                                    quote = FALSE,
                                    row.names = FALSE,
                                    col.names = FALSE,
                                    eol = "\r\n",
                                    sep = " "
                                )

                                rangefile <- normalizePath(file.path(ResultDir, "region.txt"))
                            }
                        } else {
                            rangefile <- normalizePath(file.path(DataDir, regionfile), mustWork = FALSE)
                            par_snps <- NULL
                            xtr_snps <- NULL
                            ampliconic_snps <- NULL
                        }

                        if (!is.null(rangefile)) {
                            executePlinkExcludeExtract(ResultDir, DataDir, finput, rangefile, foutput)

                            bim <- read.table(normalizePath(file.path(ResultDir, paste0(foutput, ".bim")), mustWork = FALSE))
                            bim1 <- read.table(normalizePath(file.path(DataDir, paste0(finput, ".bim")), mustWork = FALSE))

                            num_marker_excluded <- nrow(bim1) - nrow(bim)
                            rlang::inform(
                                rlang::format_error_bullets(c(
                                    "i" = paste0(num_marker_excluded, " SNPs are discarded."),
                                    "v" = paste0("PLINK files with passed SNPs are in ", ResultDir, " prefixed as ", foutput),
                                    "v" = paste0("PLINK files with discarded SNPs are in ", ResultDir, " prefixed as ", foutput, "_snps_extracted")
                                ))
                            )
                        } else {
                            rlang::inform(rlang::format_error_bullets(c("i" = "No SNPs to be discarded or flagged. No output PLINK files are created.")))
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
                        ftemp <- list.files(normalizePath(file.path(ResultDir), mustWork = FALSE), pattern = "region")
                        invisible(do.call(file.remove, list(normalizePath(file.path(ResultDir, ftemp), mustWork = FALSE))))
                    }

                    return(list(PAR = par_snps, XTR = xtr_snps, Ampliconic = ampliconic_snps))
                } else {
                    executePlinkChrFilter(ResultDir, DataDir, finput, filterCHR, foutput)

                    return(NULL) ## Added in 5.0
                }

                if (exclude == TRUE) {
                    bim <- read.table(normalizePath(file.path(ResultDir, paste0(foutput, ".bim")), mustWork = FALSE))
                    bim1 <- read.table(normalizePath(file.path(DataDir, paste0(finput, ".bim")), mustWork = FALSE))

                    num_marker_excluded <- nrow(bim1) - nrow(bim)
                    rlang::inform(
                        rlang::format_error_bullets(c(
                            "i" = paste0(num_marker_excluded, " SNPs are discarded."),
                            "v" = paste0("Plink files with passed SNPs are in ", ResultDir, " prefixed as ", foutput),
                            "v" = paste0("Plink files with discarded SNPs are in ", ResultDir, " prefixed as ", foutput, "_snps_extracted")
                        ))
                    )
                } else if (exclude == FALSE) {
                    bim <- read.table(normalizePath(file.path(ResultDir, paste0(foutput, ".bim")), mustWork = FALSE))
                    colnames(bim) <- c("CHR", "SNP", "START", "END", "A1", "A2")
                    Flagged_SNPs <- bim
                }
            },
            error = function(e) {
                rlang::abort(
                    message = e$message, 
                    class = "FilterRegion_error"
                )
            },
            warning = function(w) {
                rlang::warn(message = w$message, 
                    class = "FilterRegion_warning", 
                    .frequency = "regularly", 
                    .frequency_id = "FilterRegion_warning"
                )
            }
        )
    }
