% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GXwasR_main_functions.R
\name{ComputeLD}
\alias{ComputeLD}
\title{ComputeLD: Compute Linkage Disequilibrium (LD) for SNP Data}
\usage{
ComputeLD(DataDir, ResultDir, finput, ByCHR = FALSE, CHRnum = NULL, r2_LD)
}
\arguments{
\item{DataDir}{Character string representing the file path of the input PLINK binary files.}

\item{ResultDir}{Character string for the file path where all output files will be stored, defaulting to a temporary directory.}

\item{finput}{Character string specifying the prefix of the input PLINK binary files.}

\item{ByCHR}{Logical indicating whether to perform the computation by chromosome. The default is FALSE.}

\item{CHRnum}{If ByCHR is TRUE, specifies the chromosome number for which LD should be computed. The default is NULL.}

\item{r2_LD}{The threshold for r-squared LD values to report in the output.}
}
\value{
Returns a data frame containing the computed LD values among SNPs, read from the output file generated by PLINK.
}
\description{
This function computes linkage disequilibrium (LD) statistics for SNP data using PLINK. It allows for computation across all
SNPs or within specific chromosomes.
}
\examples{
snpld <- ComputeLD(
    DataDir = system.file("extdata", package = "GXwasR"), ResultDir = tempdir(),
    finput = "GXwasR_example", ByCHR = TRUE, CHRnum = 1, r2_LD = 0.2
)
}
