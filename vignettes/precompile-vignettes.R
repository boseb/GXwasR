#' Knit and prepare vignettes
#'
#' This script is used to precompile some long running vignettes 
#' for the package. It takes the vignette source files, knits them 
#' into R Markdown files, and copies the figures into the 
#' \code{vignettes} directory.
#'
precompile_vignettes <- function() {
  # Knit vignettes
  knitr::knit(
    input = 'vignettes/decoding_ancestry.Rmd.orig',
    output = 'vignettes/decoding_ancestry.Rmd'
  )

  # Copy figures
  vignette_figures <- list.files('figures', recursive = TRUE, full.names = TRUE)
  dir.create('vignettes/figures', recursive = TRUE, showWarnings = FALSE)
  file.copy(vignette_figures, file.path('vignettes', vignette_figures))

  # Remove temporary figures directory
  unlink('figures', recursive = TRUE)
}

