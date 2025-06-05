#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL

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
