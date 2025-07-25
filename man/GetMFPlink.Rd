% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GXwasR_main_functions.R
\name{GetMFPlink}
\alias{GetMFPlink}
\title{GetMFPlink: Getting male and female PLINK binary files.}
\usage{
GetMFPlink(
  DataDir,
  ResultDir = tempdir(),
  finput,
  foutput,
  sex,
  xplink = FALSE,
  autoplink = FALSE
)
}
\arguments{
\item{DataDir}{Character string for the file path of the input PLINK binary files.}

\item{ResultDir}{Character string for the file path where all output files will be stored. The default is \code{tempdir()}.}

\item{finput}{Character string, specifying the prefix of the input PLINK binary files.}

\item{foutput}{Character string, specifying the prefix of the output PLINK binary files.}

\item{sex}{Boolean value, 'males' or 'females', specifying output plink binary files with male or female samples.}

\item{xplink}{Boolean value, \code{TRUE} or \code{FALSE}, specifying output plink binary files with only X chromosome or not. Default is \code{FALSE.}}

\item{autoplink}{Boolean value, \code{TRUE} or \code{FALSE}, specifying output plink binary files with only autosome or not. Default is \code{FALSE.}}
}
\value{
None
}
\description{
This function prepares separate male and female PLINK binary files from combined PLINK files.
}
\examples{
DataDir <- system.file("extdata", package = "GXwasR")
ResultDir <- tempdir()
finput <- "GXwasR_example"
foutput <- "Test_output"
sex <- "females"
x <- GetMFPlink(
    DataDir = DataDir, ResultDir = ResultDir,
    finput = finput, foutput = foutput, sex = sex,
    xplink = FALSE, autoplink = FALSE
)
}
\author{
Banabithi Bose
}
