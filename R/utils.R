#' Detect the operating system type
#'
#' This helper function detects the current operating system and returns
#' a standardized string for branching logic.
#'
#' @return A character string: one of `"windows"`, `"unix"`, or `"unknown"`.
#' @examples
#' detect_os_type()
#' @keywords internal
#' @export
detect_os_type <- function() {
    os <- tolower(Sys.info()[["sysname"]])
    if (grepl("windows", os)) {
        return("windows")
    } else if (grepl("linux|darwin|unix|mac", os)) {
        return("unix")
    } else {
        return("unknown")
    }
}

#' Simulate mock summary statistics for SNPs
#'
#' Generates a mock data frame containing 100 SNPs with random identifiers,
#' distinct allele pairs, sample sizes, and Z-scores. Intended for internal
#' testing or demonstration purposes.
#'
#' @importFrom stats rnorm
#'
#' @return A data.frame with columns: SNP, A1, A2, N, Z
#' @keywords internal

simulateSumstats <- function() {
    # set.seed(123)

    n_snps <- 100000
    snp_ids <- paste0("rs", sample(1e6:2e6, n_snps))
    alleles <- c("A", "C", "G", "T")

    # Generate distinct mock allele pairs
    get_alleles <- function(n) {
        A1 <- sample(alleles, n, replace = TRUE)
        A2 <- vapply(A1, function(a) sample(setdiff(alleles, a), 1), character(1))
        list(A1 = A1, A2 = A2)
    }

    allele_data <- get_alleles(n_snps)

    data.frame(
        SNP = snp_ids,
        A1 = allele_data$A1,
        A2 = allele_data$A2,
        N = sample(50000:100000, n_snps, replace = TRUE),
        Z = rnorm(n_snps)
    )
}

#' @title Internal Logging Helper
#' @description
#' Appends one or more lines to a log file if a valid path is provided.
#' Optionally prepends a timestamp to each message.
#'
#' @param ... Character strings to be logged. Each will be combined with line breaks.
#' @param output.file Character string specifying the file path to write to. If empty (`""`), no logging occurs.
#' @param sep Character string to separate input lines. Defaults to `"\n"`.
#' @param timestamp Logical, whether to prepend a timestamp to the log entry. Default is `TRUE`.
#'
#' @return Invisible `NULL`. Used for side-effects only.
#'
#' @keywords internal
#' @noRd
log_output <- function(..., output.file, sep = "\n", timestamp = TRUE) {
    if (nzchar(output.file)) {
        con <- file(output.file, open = "a")
        on.exit(close(con))
        msg <- paste(..., sep = sep)
        if (timestamp) {
            msg <- paste(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg)
        }
        writeLines(msg, con = con)
    }
    invisible(NULL)
}
