% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{simulateSumstats}
\alias{simulateSumstats}
\title{Simulate mock summary statistics for SNPs}
\usage{
simulateSumstats()
}
\value{
A data.frame with columns: SNP, A1, A2, N, Z
}
\description{
Generates a mock data frame containing 100 SNPs with random identifiers,
distinct allele pairs, sample sizes, and Z-scores. Intended for internal
testing or demonstration purposes.
}
\keyword{internal}
