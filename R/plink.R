# PLINK Helpers 
## PLINK Dependency Check
verifyPlink <- function() {
    os_type <- detect_os_type()

    # Helper to determine if the found plink binary is the bioinformatics version
    is_bioinformatics_plink <- function(plink_path) {
        out <- tryCatch(
            {
                system2(plink_path, "--version", stdout = TRUE, stderr = TRUE)
            },
            error = function(e) {
                character()
            }
        )
        if (length(out) == 0) {
            return(FALSE)
        }

        has_plink_signature <- any(grepl("^PLINK v[0-9.]+", out))
        is_not_putty <- !any(grepl("plink:|Release|PuTTY|Build platform", out, ignore.case = TRUE))

        has_plink_signature && is_not_putty
    }

    # 1. Check PLINK_PATH environment variable
    env_path <- Sys.getenv("PLINK_PATH", unset = NA)
    if (!is.na(env_path)) {
        resolved_env_path <- normalizePath(env_path, mustWork = FALSE)
        if (file.exists(resolved_env_path) && is_bioinformatics_plink(resolved_env_path)) {
            return(resolved_env_path)
        } else if (file.exists(resolved_env_path)) {
            rlang::warn(
                message = glue::glue("PLINK_PATH is set but points to a non-bioinformatics 'plink' binary: {resolved_env_path}"), 
                .frequency = "regularly", 
                .frequency_id = "plink_check",
                class = "verifyPlink_env_warning"
            )
        }
    }

    # 2. Check system path via Sys.which
    sys_path <- Sys.which("plink")
    # On Windows, check for `.exe` if missing
    if (os_type == "windows" && !grepl("\\.exe$", sys_path, ignore.case = TRUE)) {
        sys_path <- paste0(sys_path, ".exe")
    }
    # Normalize and check existence
    resolved_sys_path <- normalizePath(sys_path, mustWork = FALSE)
    if (nzchar(resolved_sys_path) && file.exists(resolved_sys_path)) {
        if (is_bioinformatics_plink(resolved_sys_path)) {
            return(resolved_sys_path)
        } else {
            rlang::warn(
                message = glue::glue("System 'plink' found at {resolved_sys_path} but it does not appear to be the bioinformatics version."), 
                .frequency = "regularly", 
                .frequency_id = "plink_check",
                class = "verifyPlink_sys_warning"
            )
        }
    }

    # 3. Failure message
    rlang::abort(
        message = rlang::format_error_bullets(c(
            "PLINK binary not found.",
            "x" = "Attempted to locate the 'PLINK_PATH' environment variable and 'plink' in system PATH.",
            "i" = "Ensure PLINK is available and executable (bioinformatics version).",
            "i" = "You can permanently set the path to 'plink' by running:",
            " " = "usethis::edit_r_environ()  # then add a line like: PLINK_PATH=/full/path/to/plink"
        )),
        class = "plink_not_found"
    )
        }


## PLINK Binary Location
plink <- function() {
    plink_executable <- verifyPlink()
    plink_info <- system2(plink_executable, args = "--version", stdout = TRUE)
    rlang::inform(
        message = paste("Using", plink_info[[1]]),
        .frequency = "regularly",
        .frequency_id = "plink_check"
    )
    plink_executable
}

## Function 4
########## Added in 3.0
executePlink <- function(args, ResultDir) {
    # globalVariables("ResultDir")
    tryCatch(
        {
            # Redirect stderr to null to suppress warning messages
            stderr_dest <- ifelse(.Platform$OS.type == "windows", "NUL", "/dev/null")
            invisible(
                sys::exec_wait(
                    plink(),
                    args = args,
                    std_err = stderr_dest
                )
            )
        },
        error = function(e) {
            stop("An error occurred while executing Plink: ", e$message)
        }
    )
}

## Function 6
########## Added in 3.0
plinkExcludeExtract <- function(DataDir, finput, ResultDir, foutput, region_file_path) {
    # Plink command for excluding SNPs
    plinkArgsExclude <- c(
        "--bfile", file.path(DataDir, finput),
        "--exclude", "range", region_file_path,
        "--allow-no-sex", ## Adding in 4.0
        "--make-bed",
        "--out", file.path(ResultDir, foutput),
        "--silent"
    )
    executePlink(plinkArgsExclude, ResultDir)

    # Plink command for extracting SNPs
    plinkArgsExtract <- c(
        "--bfile", file.path(DataDir, finput),
        "--extract", "range", region_file_path,
        "--allow-no-sex",
        "--make-bed",
        "--out", file.path(ResultDir, paste0(foutput, "_snps_extracted")),
        "--silent"
    )
    executePlink(plinkArgsExtract, ResultDir)
}

## Function 19
######### Added in 3.0
executePlinkForIBD <- function(ResultDir, IBD, outFileName) {
    ####### Added in final version #######
    executePlinkAd(ResultDir, c(
        "--bfile", normalizePath(file.path(ResultDir, "foutput"), mustWork = FALSE),
        "--indep-pairwise", 50, 5, 0.02, # We made these as hard thresholds.
        "--allow-no-sex", ## Adding in 4.0
        "--out", normalizePath(file.path(ResultDir, "foutput"), mustWork = FALSE),
        "--silent"
    ))

    # Extract pruned SNPs based on the .prune.in file
    executePlinkAd(ResultDir, c(
        "--bfile", normalizePath(file.path(ResultDir, "foutput"), mustWork = FALSE), # Original data
        "--extract", normalizePath(file.path(ResultDir, "foutput.prune.in"), mustWork = FALSE), # Use pruned SNP list
        "--allow-no-sex",
        "--make-bed",
        "--out", normalizePath(file.path(ResultDir, "foutput1"), mustWork = FALSE),
        "--silent" # Final file with pruned SNPs
    ))

    ######################################


    ibdArgs <- c(
        "--bfile", normalizePath(file.path(ResultDir, "foutput1"), mustWork = FALSE),
        "--genome",
        "--out", normalizePath(file.path(ResultDir, outFileName), mustWork = FALSE),
        "--silent"
    )
    if (!is.null(IBD)) {
        ibdArgs <- c(ibdArgs, "--min", IBD)
    }
    executePlink(ibdArgs, ResultDir)
}

## Function 22
######### Added in 3.0
updatePlinkFilesWithIBDFilter <- function(ResultDir, foutput, failed_ibd) {
    if (!is.null(failed_ibd)) {
        write.table(failed_ibd, file = normalizePath(file.path(ResultDir, "samples_failed_ibd"), mustWork = FALSE), quote = FALSE, row.names = FALSE, col.names = FALSE)
        removeSamplesArgs <- c(
            "--bed", normalizePath(file.path(ResultDir, paste0("foutput", ".bed")), mustWork = FALSE),
            "--bim", normalizePath(file.path(ResultDir, paste0("foutput", ".bim")), mustWork = FALSE),
            "--fam", normalizePath(file.path(ResultDir, paste0("foutput", ".fam")), mustWork = FALSE),
            "--remove", normalizePath(file.path(ResultDir, "/samples_failed_ibd"), mustWork = FALSE),
            "--allow-no-sex", "--make-bed",
            "--out", normalizePath(file.path(ResultDir, foutput), mustWork = FALSE),
            "--silent"
        )
        executePlink(removeSamplesArgs, ResultDir)
    } else {
        removeSamplesArgs <- c(
            "--bed", normalizePath(file.path(ResultDir, paste0("foutput", ".bed")), mustWork = FALSE),
            "--bim", normalizePath(file.path(ResultDir, paste0("foutput", ".bim")), mustWork = FALSE),
            "--fam", normalizePath(file.path(ResultDir, paste0("foutput", ".fam")), mustWork = FALSE),
            "--make-bed", "--out", normalizePath(file.path(ResultDir, foutput), mustWork = FALSE),
            "--silent"
        )
        executePlink(removeSamplesArgs, ResultDir)
    }
}

## Function 23
######### Added in 3.0
executeMakeBed <- function(ResultDir, foutput) {
    makeBedArgs <- c(
        "--bfile", normalizePath(file.path(ResultDir, foutput), mustWork = FALSE),
        "--make-bed",
        "--out", normalizePath(file.path(ResultDir, foutput), mustWork = FALSE),
        "--silent"
    )
    executePlink(makeBedArgs, ResultDir)
}

## Function 24
######### Added in 3.0
# Adjusted Helper Function
executePlinkAd <- function(ResultDir, args) {
    tryCatch(
        {
            stderr_dest <- ifelse(.Platform$OS.type == "windows", "NUL", "/dev/null")
            sys::exec_wait(plink(), args = args, std_err = stderr_dest)
        },
        error = function(e) {
            stop("An error occurred while executing Plink: ", e$message)
        }
    )
}

## Function 25
######### Added in 3.0
# Helper Function to Set PLINK Flags
setPlinkFlags <- function(maf, geno, hwe, hweCase, hweControl) {
    MAF <- if (!is.null(maf)) "--maf" else NULL
    GENO <- if (!is.null(geno)) "--geno" else NULL

    if (!is.null(hwe)) {
        if (!is.null(hweControl)) {
            rlang::inform(rlang::format_error_bullets(c("i" = "Since hwe is not NULL, hweControl should be NULL. Setting hweControl = NULL implicitly.")))
            hweControl <- NULL
        } else if (!is.null(hweCase)) {
            rlang::inform(rlang::format_error_bullets(c("i" = "Since hwe is not NULL, hweCase should be NULL. Setting hweCase = NULL implicitly.")))
            hweCase <- NULL
        }
        HWE <- "--hwe"
    } else {
        HWE <- NULL
    }

    if (is.null(hweCase) && !is.null(hweControl)) {
        stop("hweControl cannot be NULL if hweCase is not NULL.")
    } else if (!is.null(hweCase) && is.null(hweControl)) {
        stop("hweCase cannot be NULL if hweControl is not NULL.")
    }

    HWECase <- if (!is.null(hweCase)) "--hwe" else NULL
    HWECon <- if (!is.null(hweControl)) "--hwe" else NULL

    if (!is.null(hweCase) || !is.null(hweControl)) {
        HWE <- NULL
        hwe <- NULL
    }

    return(list(MAF = MAF, GENO = GENO, HWE = HWE, HWECase = HWECase, HWECon = HWECon))
}

## Function 26
######### Added in 3.0
# Helper Function to Remove Ambiguous SNPs
removeAmbiguousSNPs <- function(DataDir, ResultDir, finput) {
    bimFilePath <- normalizePath(file.path(DataDir, paste0(finput, ".bim")), mustWork = FALSE)
    study <- read.table(file = bimFilePath, stringsAsFactors = FALSE)

    # Identifying Ambiguous SNPs
    study_AT <- study[study[, 5] == "A" & study[, 6] == "T", 2, drop = FALSE]
    study_TA <- study[study[, 5] == "T" & study[, 6] == "A", 2, drop = FALSE]
    study_GC <- study[study[, 5] == "G" & study[, 6] == "C", 2, drop = FALSE]
    study_CG <- study[study[, 5] == "C" & study[, 6] == "G", 2, drop = FALSE]

    # Identifying Indels
    study_indel1 <- study[!which(study[, 5] != "A" & study[, 5] != "T" & study[, 5] != "G" & study[, 5] != "C"), 2, drop = FALSE]
    study_indel2 <- study[!which(study[, 6] != "A" & study[, 6] != "T" & study[, 6] != "G" & study[, 6] != "C"), 2, drop = FALSE]

    study_SNP <- rbind(study_AT, study_TA, study_GC, study_CG, study_indel1, study_indel2)

    # Write Ambiguous SNPs to a file
    outputFile <- normalizePath(file.path(ResultDir, "study_SNP"), mustWork = FALSE)
    write.table(study_SNP, file = outputFile, quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\r\n", sep = " ")

    # Remove Ambiguous SNPs using PLINK
    executePlinkAd(ResultDir, args = c(
        "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE),
        "--exclude", outputFile, "--allow-no-sex",
        "--make-bed",
        "--out", normalizePath(file.path(ResultDir, "NoAmbiguousSNP"), mustWork = FALSE),
        "--silent"
    ))

    # Return the count of removed SNPs
    return(nrow(study_SNP))
}

## Function 27
######### Added in 3.0
# Helper Function to Apply Filters with PLINK
applyFiltersWithPlink <- function(ResultDir, DataDir, finput, MAF, maf, GENO, geno, HWE, hwe) {
    # Execute PLINK command
    executePlinkAd(ResultDir, args = c(
        "--bfile", normalizePath(file.path(ResultDir, "NoAmbiguousSNP"), mustWork = FALSE),
        MAF, maf,
        GENO, geno,
        HWE, hwe,
        "--allow-no-sex",
        "--make-bed",
        "--out", normalizePath(file.path(ResultDir, "filtered_temp1"), mustWork = FALSE),
        "--silent"
    ))

    # Check if the file exists and print relevant messages
    if (file.exists(normalizePath(file.path(ResultDir, "filtered_temp1.bed"), mustWork = FALSE))) {
        rlang::inform(rlang::format_error_bullets(c("v" = "Thresholds for maf, geno and hwe worked.")))
        logContents <- readLines(normalizePath(file.path(ResultDir, "filtered_temp1.log"), mustWork = FALSE))
        rlang::inform(rlang::format_error_bullets(c("i" = grep("variants removed", logContents, value = TRUE))))
    } else {
        rlang::inform(rlang::format_error_bullets(c("x" = "Error applying thresholds or file not found.")))
        executePlinkAd(ResultDir, args = c(
            "--bfile", normalizePath(file.path(DataDir, finput), mustWork = FALSE),
            "--make-bed",
            "--out", normalizePath(file.path(ResultDir, "filtered_temp1"), mustWork = FALSE),
            "--silent"
        ))
    }
}

## Function 28
######### Added in 3.0
# Helper Function to Apply HWE Filters and monomorpic part with PLINK
applyCaseControlFilters <- function(ResultDir, fam4, casecontrol, HWECase, hweCase, HWECon, hweControl) {
    nextFile <- "filtered_temp1" # Default next file

    if (length(unique(fam4$V6)) >= 2 && casecontrol) {
        # Filter for cases
        executePlinkAd(ResultDir, args = c(
            "--bfile", normalizePath(file.path(ResultDir, "filtered_temp1"), mustWork = FALSE),
            "--filter-cases", HWECase, hweCase,
            "--allow-no-sex",
            "--make-bed",
            "--out", normalizePath(file.path(ResultDir, "filtered_temp_hwe_case_filtered"), mustWork = FALSE),
            "--silent"
        ))
        printHWEMessages(ResultDir, "/filtered_temp_hwe_case_filtered.log", "In cases")

        # Filter for controls
        executePlinkAd(ResultDir, args = c(
            "--bfile", normalizePath(file.path(ResultDir, "/filtered_temp1"), mustWork = FALSE),
            "--filter-controls", HWECon, hweControl,
            "--allow-no-sex",
            "--make-bed",
            "--out", normalizePath(file.path(ResultDir, "filtered_temp_hwe_control_filtered"), mustWork = FALSE),
            "--silent"
        ))
        printHWEMessages(ResultDir, "/filtered_temp_hwe_control_filtered.log", "In controls")

        # Merge case and control filtered files
        MergeRegion(DataDir = ResultDir, ResultDir, "filtered_temp_hwe_case_filtered", "filtered_temp_hwe_control_filtered", "filtered_temp2", use_common_snps = TRUE)

        nextFile <- "filtered_temp2"
    } else {
        if (length(unique(fam4$V6)) == 1) {
            rlang::inform(rlang::format_error_bullets(c("i" = "There is no case-control status in the plink files. Setting casecontrol = FALSE implicitly.")))
            casecontrol <- FALSE
        } else {
            casecontrol <- FALSE
        }
    }

    # Processing for monomorphic SNPs
    if (file.exists(normalizePath(file.path(ResultDir, paste0(nextFile, ".bed")), mustWork = FALSE))) {
        executePlinkAd(ResultDir, args = c(
            "--bfile", normalizePath(file.path(ResultDir, nextFile), mustWork = FALSE),
            "--freq", "--make-bed", "--allow-no-sex",
            "--out", normalizePath(file.path(ResultDir, "filtered_temp4"), mustWork = FALSE),
            "--silent"
        ))
    } else {
        rlang::inform(rlang::format_error_bullets(c("x" = "Something went wrong.")))
        rlang::inform(rlang::format_error_bullets(c("x" = grep("Error", readLines(normalizePath(file.path(ResultDir, paste0(nextFile, ".log")), mustWork = FALSE)), value = TRUE))))
    }

    return(casecontrol)
}

## Function 32
executePlinkWithParams <- function(ResultDir, filtered_temp, exclude, excludemono, excluderange, highLD_regions, indep, window_size, step_size, r2_threshold) {
    if (is.null(excluderange)) {
        range <- NULL
    } else {
        range <- "range"
    }

    executePlinkAd(ResultDir, args = c(
        "--bfile", normalizePath(file.path(ResultDir, filtered_temp), mustWork = FALSE),
        exclude, excludemono, excluderange, range, highLD_regions,
        "--make-bed",
        "--allow-no-sex",
        "--out", normalizePath(file.path(ResultDir, paste0(filtered_temp, "dummy_processed")), mustWork = FALSE),
        "--silent"
    ))

    if (is.null(indep)) {
        # Extract prunned snps
        executePlinkAd(ResultDir, args = c(
            "--bfile", normalizePath(file.path(ResultDir, paste0(filtered_temp, "dummy_processed")), mustWork = FALSE),
            "--make-bed",
            "--out", normalizePath(file.path(ResultDir, paste0(filtered_temp, "_processed")), mustWork = FALSE),
            "--silent"
        ))
    } else {
        executePlinkAd(ResultDir, args = c(
            "--bfile", normalizePath(file.path(ResultDir, paste0(filtered_temp, "dummy_processed")), mustWork = FALSE),
            indep, window_size, step_size, r2_threshold,
            "--allow-no-sex",
            "--out", normalizePath(file.path(ResultDir, paste0(filtered_temp, "pruned_processed")), mustWork = FALSE),
            "--silent"
        ))
        # Extract prunned snps
        executePlinkAd(ResultDir, args = c(
            "--bfile", normalizePath(file.path(ResultDir, paste0(filtered_temp, "dummy_processed")), mustWork = FALSE),
            "--extract", normalizePath(file.path(ResultDir, paste0(filtered_temp, "pruned_processed.prune.in")), mustWork = FALSE),
            "--make-bed",
            "--out", normalizePath(file.path(ResultDir, paste0(filtered_temp, "_processed")), mustWork = FALSE),
            "--silent"
        ))
    }
}

## Function 35
######### Added in 3.0
# Helper function to apply SNP missingness filter
applySNPmissCCFilter <- function(ResultDir, SNPmissCC, diffmissFilter, foutput) {
    if (!is.null(SNPmissCC) && diffmissFilter) {
        write.table(
            SNPmissCC,
            file = normalizePath(file.path(ResultDir, "SNPdifCallrate"), mustWork = FALSE),
            quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\r\n", sep = " "
        )
        executePlinkAd(ResultDir, args = c(
            "--bfile", normalizePath(file.path(ResultDir, "filtered_temp4_processed"), mustWork = FALSE), ## TEST
            "--exclude", normalizePath(file.path(ResultDir, "SNPdifCallrate"), mustWork = FALSE),
            "--allow-no-sex",
            "--make-bed",
            "--out", normalizePath(file.path(ResultDir, foutput), mustWork = FALSE),
            "--silent"
        ))
    } else {
        rlang::inform(rlang::format_error_bullets(c("i" = "No SNP with differential missingness between cases and controls.")))
        executePlinkAd(ResultDir, args = c(
            "--bfile", normalizePath(file.path(ResultDir, "filtered_temp4_processed"), mustWork = FALSE), ## TEST
            "--allow-no-sex",
            "--make-bed",
            "--out", normalizePath(file.path(ResultDir, foutput), mustWork = FALSE),
            "--silent"
        ))
    }
}


# PLINK Functions
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
#' DataDir <- system.file("extdata", package = "GXwasR")
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
            return(invisible(NULL))
        },
        error = function(e) {
            rlang::abort(
                message = e$message, 
                class = "FilterPlinkSample_error"
            )
        },
        warning = function(w) {
            rlang::warn(message = w$message, 
                class = "FilterPlinkSample_warning", 
                .frequency = "regularly", 
                .frequency_id = "FilterPlinkSample_warning"
            )
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
#' DataDir <- system.file("extdata", package = "GXwasR")
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
            return(invisible(NULL))
        },
        error = function(e) {
            rlang::abort(
                message = e$message, 
                class = "GetMFPlink_error"
            )
        },
        warning = function(w) {
            rlang::warn(
                message = w$message, 
                class = "GetMFPlink_warning", 
                .frequency = "regularly", 
                .frequency_id = "GetMFPlink_warning"
            )
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
#' @importFrom dplyr tibble
#' @importFrom rlang abort warn
#' 
#' @return Invisible. tibble containing summary stats
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

            summary <- tibble(n_chr = No.of.chr, unique_chr = list(chr_unique = unique(bim$V1)), n_snps = No.of.snps, n_samples = No.of.samples)
            
            rlang::inform(
                rlang::format_error_bullets(c(
                    "i" = paste("Number of chromosomes:", No.of.chr),
                    " " = paste("  - Chr:", unique(bim$V1)),
                    "i" = paste("Total number of SNPs:", No.of.snps),
                    "i" = paste("Total number of samples:", No.of.samples)
                ))
            )
            return(invisible(summary))
        },
        error = function(e) {
            rlang::abort(
                e$message, 
                class = "PlinkSummary_error"
            )
        },
        warning = function(w) {
            rlang::warn(
                w$message, 
                class = "PlinkSummary_warning"
            )
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
#' DataDir <- system.file("extdata", package = "GXwasR")
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
                        rlang::inform(message = glue::glue("Processing chromosome: {chr}"))
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
                    rlang::inform(message = "Famfile is NULL. The generated .fam file will have missing phenotypes.")
                }
            }

            rlang::inform(
                format_error_bullets(c(
                    "v" = paste("Output files created in ResultDir:", ResultDir)
                ))
            )
        },
        error = function(e) {
            rlang::abort(
                message =  e$message, 
                class = "PlinkVCF_error"
            )
        },
        warning = function(w) {
            rlang::warn(
                message = w$message, 
                class = "PlinkVCF_warning", 
                .frequency = "regularly", 
                .frequency_id = "PlinkVCF_warning"
            )
        }
    )
}


MFsplitPlink <- function(DataDir, ResultDir, finput, foutput, sex, xplink = FALSE, autoplink = FALSE) {
    if (!checkFiles(DataDir, finput)) {
        stop("Missing required Plink files in the specified DataDir.")
    }
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
                "--chr", 23,
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
                "--not-chr", 23,
                "--make-bed",
                "--out", normalizePath(file.path(ResultDir, foutput), mustWork = FALSE),
                "--silent"
            ),
            std_out = FALSE,
            std_err = FALSE
        ))
    }

    rlang::inform(
        message = rlang::format_error_bullets(paste("Stratified test is running for", sex))
    )
}
