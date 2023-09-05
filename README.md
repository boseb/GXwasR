# GXwasR
GXwasR: A tool for conducting sex-aware quality control, association analysis, and testing various models of sex-dependent genetic effects in complex traits

This package implements various statistical genetics models for Genome-Wide Association (GWA) and X-Chromosome Wide Association (XWA) analyses in a sex-combined or sex-stratified way considering X-Chromosome Inactivation (XCI) pattern. In addition to association analysis, the package also enables testing for sex differences in genetic effects, including the implementation of specific models and applying best practices for additional quality control (QC) of genetic data required for these tests. The package includes twenty-five different functions in six different categories (A-F) which enable a comprehensive pipeline for sex-aware genetic association analysis of common variants with unrelated individuals.

**(A)Pre-imputation QC:**  QCsnp(); QCsample(); AncestryCheck(); SexCheck()

**(B)Post-imputation QC:** QCsnp(); QCsample2(); Xhwe(); MAFdiffSexControl(); FilterRegion()

**(C)Sex-combined and sex-stratified GWAS with XWAS:** GXwasR()

**(D)Sex-differential test:** SexDiff(); SexDiffZscore; DiffZeroOne()

**(E)High level analysis:** TestXGene(); MetaGWAS(); ComputePRS(); ComputeCorrBT(); EstimateHerit()

**(F)Utility Functions:** FilterPlinkSample(); ComputeGeneticPC(); ClumpLD(); GetMFPlink(); plinkVCF(); MergeRegion(); FilterAllele(); PlinkSummary()

# Installation from R

install_github("boseb/GXwasR")
