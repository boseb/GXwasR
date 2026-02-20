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
#' The output PLINK files with passed samples will be saved in ResultDir.
#'
#' @export
#'
#' @examples
#' DataDir <- GXwasR:::GXwasR_data()
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' foutput <- "casesPlink"
#' filter_sample <- "cases"
#' keep_remove_sample_file <- "samples_example"
#' keep <- FALSE
#'
#' FilterPlinkSample(
#'     DataDir = DataDir, ResultDir = ResultDir,
#'     finput = finput, foutput = foutput, keep_remove_sample_file = keep_remove_sample_file,
#'     keep = keep
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
            if (is.null(keep_remove_sample_file)) {
                invisible(sys::exec_wait(
                    plink(),
                    args = c(
                        "--bed", normalizePath(file.path(DataDir, paste0(finput, ".bed")), mustWork = FALSE),
                        "--bim", normalizePath(file.path(DataDir, paste0(finput, ".bim")), mustWork = FALSE),
                        "--fam", normalizePath(file.path(DataDir, paste0(finput, ".fam")), mustWork = FALSE),
                        paste0("--filter-", filter_sample),
                        "--allow-no-sex", # 4.0
                        "--make-bed",
                        "--out", normalizePath(file.path(ResultDir, foutput), mustWork = FALSE),
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
                            "--bed", normalizePath(file.path(DataDir, paste0(finput, ".bed")), mustWork = FALSE),
                            "--bim", normalizePath(file.path(DataDir, paste0(finput, ".bim")), mustWork = FALSE),
                            "--fam", normalizePath(file.path(DataDir, paste0(finput, ".fam")), mustWork = FALSE),
                            "--keep", normalizePath(file.path(DataDir, keep_remove_sample_file), mustWork = FALSE),
                            "--allow-no-sex", # 4.0
                            "--make-bed",
                            "--out", normalizePath(file.path(ResultDir, foutput), mustWork = FALSE),
                            "--silent"
                        ),
                        std_out = FALSE,
                        std_err = FALSE
                    ))
                } else if (keep == FALSE) {
                    invisible(sys::exec_wait(
                        plink(),
                        args = c(
                            "--bed", normalizePath(file.path(DataDir, paste0(finput, ".bed")), mustWork = FALSE),
                            "--bim", normalizePath(file.path(DataDir, paste0(finput, ".bim")), mustWork = FALSE),
                            "--fam", normalizePath(file.path(DataDir, paste0(finput, ".fam")), mustWork = FALSE),
                            "--remove", normalizePath(file.path(DataDir, keep_remove_sample_file), mustWork = FALSE),
                            "--allow-no-sex", # 4.0
                            "--make-bed",
                            "--out", normalizePath(file.path(ResultDir, foutput), mustWork = FALSE),
                            "--silent"
                        ),
                        std_out = FALSE,
                        std_err = FALSE
                    ))
                }
            }
            rlang::inform(rlang::format_error_bullets(c("v" = paste0(foutput, " PLINK files with desired samples are in ", ResultDir))))
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
#' Boolean value, 'males' or 'females', specifying output PLINK binary files with male or female samples.
#'
#' @param xplink
#' Boolean value, `TRUE` or `FALSE`, specifying output PLINK binary files with only X chromosome or not. Default is `FALSE.`
#'
#' @param autoplink
#' Boolean value, `TRUE` or `FALSE`, specifying output PLINK binary files with only autosome or not. Default is `FALSE.`
#'
#' @return
#' None
#'
#' @export
#'
#' @examples
#' DataDir <- GXwasR:::GXwasR_data()
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' foutput <- "Test_output"
#' sex <- "females"
#' x <- GetMFPlink(
#'     DataDir = DataDir, ResultDir = ResultDir,
#'     finput = finput, foutput = foutput, sex = sex,
#'     xplink = FALSE, autoplink = FALSE
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
                        "--bed", normalizePath(file.path(DataDir, paste0(finput, ".bed")), mustWork = FALSE),
                        "--bim", normalizePath(file.path(DataDir, paste0(finput, ".bim")), mustWork = FALSE),
                        "--fam", normalizePath(file.path(DataDir, paste0(finput, ".fam")), mustWork = FALSE),
                        paste0("--filter-", sex),
                        "--make-bed",
                        "--out", normalizePath(file.path(ResultDir, foutput), mustWork = FALSE),
                        "--silent"
                    ),
                    std_out = FALSE,
                    std_err = FALSE
                ))
            } else if (xplink == TRUE && autoplink == FALSE) {
                invisible(sys::exec_wait(
                    plink(),
                    args = c(
                        "--bed", normalizePath(file.path(DataDir, paste0(finput, ".bed")), mustWork = FALSE),
                        "--bim", normalizePath(file.path(DataDir, paste0(finput, ".bim")), mustWork = FALSE),
                        "--fam", normalizePath(file.path(DataDir, paste0(finput, ".fam")), mustWork = FALSE),
                        paste0("--filter-", sex),
                        "--chr",
                        23,
                        "--make-bed",
                        "--out", normalizePath(file.path(ResultDir, foutput), mustWork = FALSE),
                        "--silent"
                    ),
                    std_out = FALSE,
                    std_err = FALSE
                ))
            } else if (xplink == FALSE && autoplink == TRUE) {
                invisible(sys::exec_wait(
                    plink(),
                    args = c(
                        "--bed", normalizePath(file.path(DataDir, paste0(finput, ".bed")), mustWork = FALSE),
                        "--bim", normalizePath(file.path(DataDir, paste0(finput, ".bim")), mustWork = FALSE),
                        "--fam", normalizePath(file.path(DataDir, paste0(finput, ".fam")), mustWork = FALSE),
                        paste0("--filter-", sex),
                        "--not-chr",
                        23,
                        "--make-bed",
                        "--out", normalizePath(file.path(ResultDir, foutput), mustWork = FALSE),
                        "--silent"
                    ),
                    std_out = FALSE,
                    std_err = FALSE
                ))
            }

            rlang::inform(rlang::format_error_bullets(c("v" = paste0("Output PLINK files, prefixed as ", foutput, ", are in ", ResultDir))))
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

# Updated in 5.0
#' PlinkSummary: Summary of PLINK format genotype dataset
#'
#' @param DataDir A character string for the file path of the input PLINK binary files.
#' @param ResultDir A character string for the file path of the PLINK program to be set up.
#' @param finput Character string, specifying the prefix of the input PLINK binary files. This file needs to be in DataDir.
#'
#' @return This function is called for its side effect: printing summary statistics to the console. It returns `NULL` invisibly.
#' @export
#'
#' @importFrom rlang abort warn
#'
#' @examples
#' DataDir <- GXwasR:::GXwasR_data()
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
            fam <- as.data.frame(utils::read.table(file.path(DataDir, paste0(finput, ".fam"))))
            if (ncol(fam) == 5) {
                fam$V6 <- fam$V1
                fam <- fam[, c(6, seq_len(5))]
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
                ))
            )
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
#' DataDir <- GXwasR:::GXwasR_data()
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' maf_data <- executePlinkMAF(DataDir, ResultDir, finput)
executePlinkMAF <- function(DataDir, ResultDir, finput) {
    if (!checkFiles(DataDir, finput)) {
        stop("Missing required Plink files in the specified DataDir.")
    }

    # Compute MAF using PLINK and return the results as a DataFrame, and clean up intermediate files afterward.

    # Generate a unique output file prefix based on the timestamp to prevent any overwrite
    foutput <- paste0("maf_output_", format(Sys.time(), "%Y%m%d_%H%M%S"))

    # Execute PLINK using --bfile for simplified input file specification
    tryCatch(
        {
            invisible(sys::exec_wait(
                plink(), # Path to the PLINK executable
                args = c(
                    "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE), # Base filename for .bed, .bim, and .fam
                    "--freq", # Command to calculate allele frequencies
                    "--out", normalizePath(file.path(ResultDir, foutput), mustWork = FALSE),
                    "--silent" # Suppress output to standard output
                ),
                std_out = FALSE,
                std_err = FALSE
            ))
        },
        error = function(e) {
            stop("An error occurred while executing PLINK: ", e$message)
        }
    )

    # Define the output filename for the frequency result
    freqOutputFile <- normalizePath(file.path(ResultDir, paste0(foutput, ".frq")), mustWork = FALSE)

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

#' plinkVCF: Converting VCF files to PLINK binary files and vice-versa.
#'
#' @author Banabithi Bose
#'
#' @description
#' This function performs the conversion between VCF files to PLINK binary formats.
#'
#' For VCF to PLINK files conversion, if you do not specify any FAM file when you are converting from VCF to PLINK
#' format, then PLINK will just create a 'dummy' FAM file with the same name as your dataset with missing phenotypes
#' and missing sex.
#'
#' @param DataDir
#' A character string for the file path of the input PLINK binary files and all other input files.
#'
#' @param ResultDir
#' A character string for the file path where all output files will be stored. The default is `tempdir()`.
#'
#' @param finput
#' Character string, specifying the prefix of the input PLINK binary files or vcf files. This file needs to be in `DataDir`.
#'
#' @param foutput
#' Character string, specifying the prefix of the output PLINK binary files if filtering option for the SNPs is chosen.
#' The default is "FALSE".
#'
#' @param VtoP
#' Boolean value, `TRUE` or `FALSE`, specifying the conversion of VCF files to PLINK binary files or not. The default is `TRUE`.
#'
#' @param PtoV
#' Boolean value, `TRUE` or `FALSE`, specifying the conversion of PLINK binary files to VCF  files or not. The default is `TRUE`.
#'
#' @param Famfile
#' Character string, specifying the name of the original .fam file if VtoP was set to be `TRUE`. This file needs to be in `DataDir`.
#' The default is `NULL`.
#'
#' @param PVbyCHR
#' Boolean value, `TRUE` or `FALSE` specifying to do the PLINK to vcf conversion chromosome-wise or not. The default is `TRUE`.
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
#' DataDir <- GXwasR:::GXwasR_data()
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

    tryCatch(
        {
            if (PtoV) {
                if (!checkFiles(DataDir, finput)) {
                    stop("There are no Plink files in DataDir. Please specify correct DataDir path with input Plink files.")
                }

                convertPlinkToVCF <- function(prefix, chr = NULL) {
                    args <- c(
                        "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE),
                        "--recode", "vcf", "--allow-extra-chr",
                        "--out", normalizePath(file.path(ResultDir, prefix), mustWork = FALSE),
                        "--silent"
                    )
                    if (!is.null(chr)) {
                        args <- c(args, "--chr", chr)
                    }
                    executePlinkAd(ResultDir, args)

                    # Compress and index
                    vcf_path <- normalizePath(file.path(ResultDir, paste0(prefix, ".vcf")), mustWork = FALSE)
                    vcf_gz_path <- Rsamtools::bgzip(
                        file = vcf_path,
                        dest = paste0(vcf_path, ".gz"),
                        overwrite = TRUE
                    )
                    Rsamtools::indexTabix(vcf_gz_path, format = "vcf")
                }

                if (PVbyCHR) {
                    bimfile <- read.table(normalizePath(file.path(DataDir, paste0(finput, ".bim")), mustWork = FALSE))
                    chrs <- unique(bimfile$V1)
                    chrs <- gsub("^chr", "", chrs) # Normalize chromosome names
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
                vcf_file <- normalizePath(file.path(DataDir, paste0(finput, ".vcf")), mustWork = FALSE)
                if (!file.exists(vcf_file)) {
                    stop("VCF file not found in DataDir. Please specify correct directory path with input VCF files.")
                }

                executePlinkAd(ResultDir, c(
                    "--vcf", vcf_file,
                    "--keep-allele-order", "--allow-extra-chr",
                    "--make-bed", "--const-fid", "1",
                    "--out", normalizePath(file.path(ResultDir, foutput), mustWork = FALSE),
                    "--silent"
                ))

                if (!is.null(Famfile)) {
                    fam <- read.table(file.path(DataDir, Famfile))
                    write.table(fam,
                        file = file.path(ResultDir, paste0(foutput, ".fam")),
                        col.names = FALSE, row.names = FALSE, quote = FALSE
                    )
                } else {
                    message("Famfile is NULL. The generated .fam file will have missing phenotypes.")
                }
            }

            rlang::inform(
                format_error_bullets(c(
                    "v" = paste("Output files created in ResultDir:", ResultDir)
                ))
            )
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
