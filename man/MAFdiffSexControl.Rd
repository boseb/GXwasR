% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GXwasR_main_functions.R
\name{MAFdiffSexControl}
\alias{MAFdiffSexControl}
\title{MAFdiffSexControl: Test for significantly different minor allele frequency (MAF) between sexes in control samples}
\usage{
MAFdiffSexControl(
  DataDir,
  ResultDir = tempdir(),
  finput,
  filterSNP = FALSE,
  foutput = NULL
)
}
\arguments{
\item{DataDir}{Character string for the file path of the input PLINK binary files.}

\item{ResultDir}{Character string for the file path where all output files will be stored. The default is \code{tempdir()}.}

\item{finput}{Character string, specifying the prefix of the input PLINK binary files with both male and female samples.
This file needs to be in \code{DataDir}.}

\item{filterSNP}{Boolean value, \code{TRUE} or \code{FALSE} for filtering out the SNPs or not (i.e., only flagged). The default is \code{FALSE}.}

\item{foutput}{Character string, specifying the prefix of the output PLINK binary files if filtering option for the SNPs
is chosen. The default is NULL.}
}
\value{
A list object containing excluded or flagged SNPs. If \code{filterSNP} = \code{TRUE}, the output filtered PLINK binary
files will be produced inside \code{DataDir}.
}
\description{
With parameters to filter out SNPs and/or flag the SNPs, this function tests for significantly different MAF
(p-value < 0.05/no. of SNPs) between sexes in control samples solely for binary phenotypes. Since the disparities
may be caused by technical confounding or sample biases for the research cohorts, it is advised that any SNPs in
the controls with a sex difference in MAF be carefully evaluated and identified for further examination
(Khramtsova et. al., 2023). In autosomal allele frequencies, sex differences are not anticipated.
}
\examples{
DataDir <- system.file("extdata", package = "GXwasR")
ResultDir <- tempdir()
finput <- "GXwasR_example"
foutput <- "Test_output"
x <- MAFdiffSexControl(DataDir, ResultDir, finput, filterSNP = TRUE, foutput = foutput)
}
\author{
Banabithi Bose
}
