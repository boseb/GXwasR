#' Example of covariate dataframe
#'
#' @name Example_covarfile
#' @docType data
#' @format A data frame with 276 rows and 4 variables:
#' \describe{
#'   \item{FID}{Family ID; character string}
#'   \item{IID}{Individual ID; character string}
#'   \item{AGE}{Age of individual; integer}
#'   \item{testcovar}{Example binary covariate; integer (0/1)}
#' }
#' @author Banabithi Bose
#' @keywords data
NULL

#' Example of ancestry dataframe of samples
#'
#' @name example_data_study_sample_ancestry
#' @docType data
#' @format A data frame with 276 rows and 2 variables:
#' \describe{
#'   \item{V1}{Individual ID; character string}
#'   \item{V2}{Population ancestry label; character string (e.g., "EUR_GBR")}
#' }
#' @author Banabithi Bose
#' @keywords data
NULL

#' Example dataset for h2 estimation between male and female for eight different traits
#'
#' @name Example_h2data
#' @docType data
#' @format A data frame with 8 rows and 5 variables:
#' \describe{
#'   \item{ID}{Trait identifier; character string (e.g., "ADHD", "AFB")}
#'   \item{Fstat}{Estimated heritability statistic for females; numeric}
#'   \item{Fse}{Standard error for female heritability estimate; numeric}
#'   \item{Mstat}{Estimated heritability statistic for males; numeric}
#'   \item{Mse}{Standard error for male heritability estimate; numeric}
#' }
#' @author Banabithi Bose
#' @keywords data
NULL

#' Example dataset for phenotype file
#'
#' @name Example_phenofile
#' @docType data
#' @format A data frame with 276 rows and 4 variables:
#' \describe{
#'   \item{FID}{Family ID; character string}
#'   \item{IID}{Individual ID; character string}
#'   \item{Pheno1}{Example phenotype 1; integer (e.g., case/control status)}
#'   \item{Pheno2}{Example phenotype 2; numeric (e.g., quantitative trait)}
#' }
#' @author Banabithi Bose
#' @keywords data
NULL

#' Example dataset with p-value threshold
#'
#' @name Example_pthresoldfile
#' @docType data
#' @format A data frame with 7 rows and 3 variables:
#' \describe{
#'   \item{Threshold}{Nominal p-value threshold; numeric}
#'   \item{LowerBound}{Lower bound of the p-value interval; integer}
#'   \item{UpperBound}{Upper bound of the p-value interval; numeric}
#' }
#' @author Banabithi Bose
#' @keywords data
NULL

#' Example dataset for running DiffZerOne function
#'
#' @name Example_rgdata
#' @docType data
#' @format A data frame with 16 rows and 3 variables:
#' \describe{
#'   \item{Trait}{Trait name or abbreviation; character string}
#'   \item{rg}{Estimated genetic correlation with target trait; numeric}
#'   \item{se}{Standard error of the genetic correlation estimate; numeric}
#' }
#' @author Banabithi Bose
#' @keywords data
NULL

#' Example dataset for MaleWAS
#'
#' @name Mfile
#' @docType data
#' @format A data frame with 120 rows and 7 variables:
#' \describe{
#'   \item{SNP}{SNP identifier (rsID); character string}
#'   \item{CHR}{Chromosome number; integer}
#'   \item{BP}{Base-pair position (GRCh37 or GRCh38 assumed); integer}
#'   \item{A1}{Effect allele; character string}
#'   \item{BETA_M}{Effect size estimate (beta coefficient) for males; numeric}
#'   \item{SE_M}{Standard error of the male beta estimate; numeric}
#'   \item{P}{P-value for the association; numeric}
#' }
#' @author Banabithi Bose
#' @keywords data
NULL

#' Example dataset for FemaleWAS
#'
#' @name Ffile
#' @docType data
#' @format A data frame with 120 rows and 7 variables:
#' \describe{
#'   \item{SNP}{SNP identifier (rsID); character string}
#'   \item{CHR}{Chromosome number; integer}
#'   \item{BP}{Base-pair position (GRCh37 or GRCh38 assumed); integer}
#'   \item{A1}{Effect allele; character string}
#'   \item{BETA_F}{Effect size estimate (beta coefficient) for females; numeric}
#'   \item{SE_F}{Standard error of the female beta estimate; numeric}
#'   \item{P}{P-value for the association; numeric}
#' }
#' @author Banabithi Bose
#' @keywords data
NULL

#' Example dataset for regression
#'
#' @name Regression_Ex
#' @docType data
#' @format A data frame with 276 rows and 8 variables:
#' \describe{
#'   \item{SCORE}{Polygenic score; numeric}
#'   \item{SEX}{Sex coded as 1 (male) or 2 (female); integer}
#'   \item{PC1}{First ancestry principal component; numeric}
#'   \item{PC2}{Second ancestry principal component; numeric}
#'   \item{PC3}{Third ancestry principal component; numeric}
#'   \item{PC4}{Fourth ancestry principal component; numeric}
#'   \item{PC5}{Fifth ancestry principal component; numeric}
#'   \item{PC6}{Sixth ancestry principal component; numeric}
#' }
#' @author Banabithi Bose
#' @keywords data
NULL

#' Example dataset for GWAS summary statistics
#'
#' @name Summary_Stat_Ex1
#' @docType data
#' @format A data frame with 1,500 rows and 12 variables:
#' \describe{
#'   \item{CHR}{Chromosome number; integer}
#'   \item{SNP}{SNP identifier (rsID); character string}
#'   \item{BP}{Base-pair position; integer}
#'   \item{A1}{Effect allele; character string}
#'   \item{TEST}{Test type (e.g., "ADD" for additive model); character string}
#'   \item{NMISS}{Number of non-missing samples for the SNP; integer}
#'   \item{BETA}{Effect size estimate (beta coefficient); numeric}
#'   \item{SE}{Standard error of the effect size estimate; numeric}
#'   \item{L95}{Lower bound of the 95% confidence interval; numeric}
#'   \item{U95}{Upper bound of the 95% confidence interval; numeric}
#'   \item{STAT}{Test statistic (e.g., Wald z-score); numeric}
#'   \item{P}{P-value for the association; numeric}
#' }
#' @author Banabithi Bose
#' @keywords data
NULL

#' Example dataset for GWAS summary statistics
#'
#' @name Summary_Stat_Ex2
#' @docType data
#' @format A data frame with 1,500 rows and 12 variables:
#' \describe{
#'   \item{CHR}{Chromosome number; integer}
#'   \item{SNP}{SNP identifier (rsID); character string}
#'   \item{BP}{Base-pair position; integer}
#'   \item{A1}{Effect allele; character string}
#'   \item{TEST}{Test type (e.g., "ADD" for additive model); character string}
#'   \item{NMISS}{Number of non-missing samples for the SNP; integer}
#'   \item{BETA}{Effect size estimate (beta coefficient); numeric}
#'   \item{SE}{Standard error of the effect size estimate; numeric}
#'   \item{L95}{Lower bound of the 95% confidence interval; numeric}
#'   \item{U95}{Upper bound of the 95% confidence interval; numeric}
#'   \item{STAT}{Test statistic (e.g., Wald z-score); numeric}
#'   \item{P}{P-value for the association; numeric}
#' }
#' @author Banabithi Bose
#' @keywords data
NULL

#' Example dataset for XWAS summary statistics
#'
#' @name XWAS_Summary_Example
#' @docType data
#' @format A data frame with 80 rows and 8 variables:
#' \describe{
#'   \item{CHROM}{Chromosome number (23 for X chromosome); integer}
#'   \item{POS}{Genomic position (base-pair); integer}
#'   \item{ID}{SNP identifier (rsID); character string}
#'   \item{P}{P-value for association; numeric}
#'   \item{BETA}{Effect size estimate (beta coefficient); numeric}
#'   \item{A1}{Effect allele; character string}
#'   \item{A2}{Other allele; character string}
#'   \item{EAF}{Effect allele frequency; numeric}
#' }
#' @author Banabithi Bose
#' @keywords data
NULL

#' Example dataset for high LD regions in hg19 built
#'
#' @name highLD_hg19
#' @docType data
#' @format A data frame with 20 rows and 4 variables:
#' \describe{
#'   \item{chr}{Chromosome name (e.g., "chr1"); character string}
#'   \item{start}{Start position of the high-LD region (base-pair); integer}
#'   \item{end}{End position of the high-LD region (base-pair); integer}
#'   \item{index}{Region index or identifier; integer}
#' }
#' @author Banabithi Bose
#' @keywords data
NULL
