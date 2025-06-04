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
