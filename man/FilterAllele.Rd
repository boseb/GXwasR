% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GXwasR_main_functions.R
\name{FilterAllele}
\alias{FilterAllele}
\title{FilterAllele: Filtering out the multi-allelic variants}
\usage{
FilterAllele(DataDir, ResultDir, finput, foutput)
}
\arguments{
\item{DataDir}{A character string for the file path of the input PLINK binary files.}

\item{ResultDir}{A character string for the file path where all output files will be stored. The default is \code{tempdir()}.}

\item{finput}{Character string, specifying the prefix of the input PLINK binary files.}

\item{foutput}{Character string, specifying the prefix of the output PLINK binary files. If multi-allelic variants are present,
this file will be produced after filtering out these variants.}
}
\value{
\code{NULL}. After multi-allelic variant filtering, the filtered plink files with only biallelic SNPs will be saved in \code{ResultDir}.
}
\description{
This function filters out the multi-allelic SNPs from the input dataset.
}
\examples{
DataDir <- system.file("extdata", package = "GXwasR")
ResultDir <- tempdir()
finput <- "GXwasR_example"
foutput <- "Filter_Test"
x <- FilterAllele(DataDir, ResultDir, finput, foutput)
}
\author{
Banabithi Bose
}
