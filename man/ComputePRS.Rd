% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GXwasR_main_functions.R
\name{ComputePRS}
\alias{ComputePRS}
\title{ComputePRS: Computing polygenic risk score (PRS)}
\usage{
ComputePRS(
  DataDir,
  ResultDir = tempdir(),
  finput,
  summarystat,
  phenofile,
  covarfile = NULL,
  effectsize = c("BETA", "OR"),
  ldclump = FALSE,
  LDreference,
  clump_p1,
  clump_p2,
  clump_r2,
  clump_kb,
  byCHR = TRUE,
  pthreshold = c(0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5),
  highLD_regions,
  ld_prunning = FALSE,
  window_size = 50,
  step_size = 5,
  r2_threshold = 0.02,
  nPC = 6,
  pheno_type = "binary"
)
}
\arguments{
\item{DataDir}{A character string for the file path of the all the input files.}

\item{ResultDir}{A character string for the file path where all output files will be stored. The default is tempdir().}

\item{finput}{Character string, specifying the prefix of the input PLINK binary files for the genotype data i.e., the target data based on which
clumping procedure will be performed. This file needs to be in DataDir. If your target data are small (e.g. N < 500) then you can use
the 1000 Genomes Project samples. Make sure to use the population that most closely reflects represents the base sample.}

\item{summarystat}{A dataframe object with GWAS summary statistics.

The mandatory column headers in this dataframe are:
\itemize{
\item \code{CHR}(Chromosome code)
\item \code{BP}(Basepair position)
\item \code{A1} (effect allele)
\item \code{SNP} (i.e., SNP idenitifier)
\item \code{BETA} or \code{OR} (i.e., effect-size or logarithm of odds ratio)
\item \code{P} (i.e., p-values).
}

Special Notes: The first three columns needed to be \code{SNP}, \code{A1} and \code{BETA} or \code{OR}.}

\item{phenofile}{A character string, specifying the name of the mandatory phenotype file. This is a plain text file with no header line; columns
family ID, individual ID and phenotype columns. For binary trait, the phenotypic value should be coded as 0 or 1, then it will be
recognized as a case-control study (0 for controls and 1 for cases). Missing value should be represented by "-9" or "NA". The
interested phenotype column should be labeled as "Pheno1". This file needs to be in \code{DataDir}.}

\item{covarfile}{A character string, specifying the name of the covariate file which is a plain .text file with no header line; columns are family ID,
individual ID and the covariates. The default is \code{NULL}. This file needs to be in \code{DataDir}.}

\item{effectsize}{Boolean value, \code{BETA} or \code{OR}, specifying the type of the GWAS effectsize. The default is \code{BETA}.}

\item{ldclump}{Boolean value, \code{TRUE} or \code{FALSE}, specifying whether to perform clumping or not.}

\item{LDreference}{A character string, specifying the  prefix of the PLINK files of the population reference panel of the same ancestry, and ideally
the one that was used for imputing your target dataset. These files should be in \code{DataDir}.}

\item{clump_p1}{Numeric value, specifying the significance threshold for index SNPs if \code{ldclump} was set to be \code{TRUE}. The default is 0.0001.}

\item{clump_p2}{Numeric value, specifying the secondary significance threshold for clumped SNPs if \code{ldclump} was set to be \code{TRUE}. The default is 0.01}

\item{clump_r2}{Numeric value, specifying the linkage disequilibrium (LD) threshold for clumping if \code{ldclump} was set to be \code{TRUE}. The default is 0.50.}

\item{clump_kb}{Integer value, specifying the physical distance threshold in base-pair for clumping if \code{ldclump} was set to be \code{TRUE}. The default is 250.}

\item{byCHR}{Boolean value, 'TRUE' or 'FALSE', specifying chromosome-wise clumping procedure if \code{ldclump} was set to be \code{TRUE}. The default is \code{TRUE}}

\item{pthreshold}{Numeric vector, containing several p value thresholds to maximize predictive ability of the derived polygenic scores.}

\item{highLD_regions}{Character string, specifying the .txt file name with known genomic regions with high LD. The default is \code{NULL}.}

\item{ld_prunning}{Boolean value, \code{TRUE} or \code{FALSE} for LD-based filtering for computing genetic PC as covariates.}

\item{window_size}{Integer value, specifying a window size in variant count or kilobase for LD-based filtering in computing genetic PC. The default is 50.}

\item{step_size}{Integer value, specifying a variant count to shift the window at the end of each step for LD filtering in computing genetic PCs. The default is 5.}

\item{r2_threshold}{Numeric value between 0 to 1 of pairwise \eqn{r^2} threshold for LD-based filtering in computing genetic PCs. The default is 0.02.}

\item{nPC}{Positive integer value, specifying the number of genetic PCs to be included as predictor in the PRS model fit. The default is 6.}

\item{pheno_type}{Boolean value, ‘binary’ or ‘quantitative’, specifying the type of the trait. The default is ‘binary’.}
}
\value{
A list object containing a dataframe and a numeric value. The dataframe,PRS, contains four mandatory columns, such as,
IID (i.e., Individual ID), FID (i.e., Family ID), Pheno1 (i.e., the trait for PRS) and Score (i.e., the best PRS).
Other columns of covariates could be there. The numeric value, BestP contains the threshold of
of the best p-value for the best pRS model fit.

Also, the function produces several plots such as p-value thresholds vs PRS model fit and PRS distribution among male and females.
For case-control data, it shows PRS distribution among cases and controls and ROC curves as well.
}
\description{
This function calculates the polygenic risk score, which is the total of allele counts (genotypes) weighted by estimated
effect sizes from genome-wide association studies. It uses C+T filtering techniques. The users could perform clumping procedure
choromosome-wise and genome-wide. Also, the function offers the choice of including several genetic principal components along with
other covariates. Using this function, users have the freedom to experiment with various clumping and thresholding arrangements to
test a wide range of various parameter values.
}
\examples{
data("Summary_Stat_Ex1", package = "GXwasR")
data("Example_phenofile", package = "GXwasR")
data("Example_covarfile", package = "GXwasR")
data("Example_pthresoldfile", package = "GXwasR")
data("highLD_hg19", package = "GXwasR")
DataDir <- system.file("extdata", package = "GXwasR")
ResultDir <- tempdir()
finput <- "GXwasR_example"
summarystat <- Summary_Stat_Ex1[, c(2, 4, 7, 1, 3, 12)]
phenofile <- Example_phenofile # Cannot be NULL
# The interested phenotype column should be labeled as "Pheno1".
covarfile <- Example_covarfile
clump_p1 <- 0.0001
clump_p2 <- 0.0001
clump_kb <- 500
clump_r2 <- 0.5
byCHR <- TRUE
pthreshold <- Example_pthresoldfile$Threshold
ld_prunning <- TRUE
highLD_regions <- highLD_hg19
window_size <- 50
step_size <- 5
r2_threshold <- 0.02
nPC <- 6 # We can incorporate PCs into our PRS analysis to account for population stratification.
pheno_type <- "binary"

PRSresult <- ComputePRS(DataDir, ResultDir, finput, summarystat, phenofile, covarfile,
    effectsize = "BETA", LDreference = "GXwasR_example", ldclump = FALSE, clump_p1, clump_p2,
    clump_r2, clump_kb, byCHR = TRUE, pthreshold = pthreshold, highLD_regions = highLD_regions,
    ld_prunning = TRUE, window_size = 50, step_size = 5, r2_threshold = 0.02, nPC = 6,
    pheno_type = "binary"
)

## This table shows 10 samples with phenotype, covariates and a PRS column.
PRS <- PRSresult$PRS
PRS[seq_len(10), ]
## The best threshold
BestPvalue <- PRSresult$BestP$Threshold
BestPvalue
}
\author{
Banabithi Bose
}
