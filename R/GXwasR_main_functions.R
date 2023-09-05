#' AncestryCheck: Evaluation of the samples' ancestry label.
#'
#' @author Banabithi Bose
#'
#' @description This function displays the result of the ancestry analysis in a color-coded scatter plot of the first two principal components for samples of the reference populations and the study population. Specifically, it compares the study samples' ancestry labels to a panel representing a reference population, and it also flags the outlier samples with respect to a chosen reference population. The function first filters the reference and study data for non-A-T or G-C SNPs. It next conducts LD pruning, fixes the chromosome mismatch between the reference and study datasets, checks for allele flips, updates the locations, and flips the alleles. The two datasets are then joined, and the resulting genotype dataset is subjected to Principal Component Analysis (PCA). The detection of population structure down to the level of the reference dataset can then be accomplished using PCA on this combined genotyping panel.  For instance, the center of the European reference samples is determined using the data from principal components 1 and 2 (median(PC1 europeanRef), median(PC2 europeanRef)). It determines the European reference samples' maximum Euclidean distance (maxDist) from this center. All study samples that are non-European, or outliers, are those whose Euclidean distances from the center are more than or equal to the radius r= outlier threshold* maxDist. This function utilizes the HapMap phase 3 data in NCBI 36 and 1000GenomeIII in CGRCh37. Both study and reference datasets should be of the same genome build. If not, users need to lift over one of the datasets to the same build.
#'
#' @param DataDir A character string for the file path of the input PLINK binary files.
#' @param ResultDir A character string for the file path where all output files will be stored. The default is tempdir().
#' @param finput Character string, specifying the prefix of the input PLINK binary files for the study samples.
#' @param reference Boolean value,'HapMapIII_NCBI36' and 'ThousandGenome', specifying Hapmap Phase3 (1) and 1000 Genomes phase III (2) reference population, respectively. The default is 'HapMapIII_NCBI36'.

#' @param filterSNP Boolean value, 'TRUE' or 'FALSE' for filtering out the SNPs. The default is TRUE. We recommend setting it FALSE only when the users are sure that they could join the study and reference samples directly.
#' @param studyLD Boolean value, 'TRUE' or 'FALSE' for applying linkage disequilibrium (LD)-based filtering on study genotype data.
#' @param studyLD_window_size Integer value, specifying a window size in variant count or kilobase for LD-based filtering of the variants for the study data.
#' @param studyLD_step_size Integer value, specifying a variant count to shift the window at the end of each step for LD filtering for the study data.
#' @param studyLD_r2_threshold Numeric value between 0 to 1 of pairwise r^2 threshold for LD-based filtering for the study data.
#' @param referLD Boolean value, 'TRUE' or 'FALSE' for applying linkage disequilibrium (LD)-based filtering on reference genotype data.
#' @param referLD_window_size Integer value, specifying a window size in variant count or kilobase for LD-based filtering of the variants for the reference data.
#' @param referLD_step_size Integer value, specifying a variant count to shift the window at the end of each step for LD filtering for the reference data.
#' @param referLD_r2_threshold Numeric value between 0 to 1 of pairwise r^2 threshold for LD-based filtering for the reference data.
#' @param highLD_regions A dataframe with known high LD regions. A dataframe from Anderson et. al, 2010 (3) is provided with the package.
#' @param study_pop A dataframe containing two columns for study in first column, sample ID (i.e., IID) and in second column, the ancestry label.
#' @param outlier Boolean value, 'TRUE' or 'FALSE', specifying outlier detection will be performed or not.
#' @param outlierOf Chracter string, specifying the reference ancestry name for detecting outlier samples. The default is "outlierOf = "CEU".
#' @param outlier_threshold Numeric value, specifying the threshold to be be used to detect outlier samples. This threshold will be multiplied with the Eucledean distance from the center of the PC 1 and PC2 to the maximum Euclidean distance of the reference samples. Study samples outside this distance will be considered as outlier. The default is 3.
#' @importFrom data.table as.data.table
#' @importFrom ggplot2 ggplot aes geom_hline geom_vline guides geom_point guide_legend scale_shape_manual
#' @importFrom DescTools %like%
#'
#' @return A dataframe with the IDs of non-European samples as outliers.
#'
#' @references (1) The International HapMap 3 Consortium. Integrating common and rare genetic variation in diverse human populations. Nature 467, 52–58 (2010). https://doi.org/10.1038/nature09298
#' (2) The 1000 Genomes Project Consortium. A global reference for human genetic variation. Nature 526, 68–74 (2015). https://doi.org/10.1038/nature15393
#' (3) Anderson CA, Pettersson FH, Clarke GM, Cardon LR, Morris AP, Zondervan KT. Data quality control in genetic case-control association studies. Nat Protoc. 2010 Sep;5(9):1564-73. doi: 10.1038/nprot.2010.116. Epub 2010 Aug 26. PMID: 21085122; PMCID: PMC3025522.
#' @export
#'
#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' reference <- "HapMapIII_NCBI36"
#' data("GXwasRData")
#' #load("~/b1137/BBose/GXwasR/data/highLD_hg19.Rda")
#' highLD_regions <- highLD_hg19 #"high-LD-regions-hg19-GRCh37.txt" # four columns, chr, start, end, index
#' #load("~/b1137/BBose/GXwasR/data/example_data_study_sample_ancestry.Rda")
#' study_pop <- example_data_study_sample_ancestry #PreimputeEX
#' studyLD_window_size = 50
#' studyLD_step_size = 5
#' studyLD_r2_threshold = 0.02
#' filterSNP = TRUE
#' studyLD = TRUE
#' referLD = TRUE
#' referLD_window_size = 50
#' referLD_step_size = 5
#' referLD_r2_threshold = 0.02
#' outlier = TRUE
#' outlier_threshold = 3
#' x <- AncestryCheck(DataDir = DataDir, ResultDir = ResultDir, finput = finput, reference = "HapMapIII_NCBI36",highLD_regions = highLD_regions,
#' study_pop = study_pop, studyLD = studyLD, referLD = referLD, outlierOf = "CEU", outlier = outlier, outlier_threshold = outlier_threshold )
#'
#' #x <- AncestryCheck(DataDir = DataDir, ResultDir = ResultDir, finput = finput, reference = "ThousandGenome",highLD_regions = highLD_regions,
#' #study_pop = study_pop, studyLD = studyLD, referLD = referLD, outlierOf = "EUR", outlier = outlier, outlier_threshold = outlier_threshold )

AncestryCheck <-
  function(DataDir,
           ResultDir = tempdir(),
           finput,
           reference = c("HapMapIII_NCBI36","ThousandGenome"),
           filterSNP = TRUE,
           studyLD = TRUE,
           studyLD_window_size = 50,
           studyLD_step_size = 5,
           studyLD_r2_threshold = 0.02,
           referLD = FALSE,
           referLD_window_size = 50,
           referLD_step_size = 5,
           referLD_r2_threshold = 0.02,
           highLD_regions,
           study_pop,
           outlier = FALSE,
           outlierOf = "CEU",
           outlier_threshold = 3) {

    if (file.exists(paste0(DataDir, "/", finput, ".bed")) &&
        file.exists(paste0(DataDir, "/", finput, ".bim")) &&
        file.exists(paste0(DataDir, "/", finput, ".fam"))) {

      setupPlink(ResultDir)

    } else{
      writeLines(
        "There are no Plink files in DataDir.\nPlease specify correct DataDir path with input Plink files."
      )
    }


    #Download_reference(refdata = reference, wdir = DataDir)
    Download_reference(refdata = reference, wdir = ResultDir)

    write.table(highLD_regions, file = paste0(ResultDir,"/high-LD-regions-temp.txt"),quote = FALSE, row.names = FALSE, col.names = FALSE)

    highLD_regions <- "high-LD-regions-temp.txt"

    if (filterSNP == TRUE) {
      ## Reading .bim file of study and reference data and filtered out AT, GC SNPs
      study <-
        read.table(file = paste0(DataDir, "/",finput, ".bim"),
                   stringsAsFactors = FALSE)
      study_AT <-
        study[which(study[, 5] == "A" &
                      study[, 6] == "T"), 2, drop = FALSE]
      study_TA <-
        study[which(study[, 5] == "T" &
                      study[, 6] == "A"), 2, drop = FALSE]
      study_GC <-
        study[which(study[, 5] == "G" &
                      study[, 6] == "C"), 2, drop = FALSE]
      study_CG <-
        study[which(study[, 5] == "C" &
                      study[, 6] == "G"), 2, drop = FALSE]

      study_SNP <- rbind(study_AT, study_TA, study_GC, study_CG)

      write.table(
        study_SNP,
        file = paste0(ResultDir, "/","study_SNP"),
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE,
        eol = "\r\n",
        sep = " "
      )
      invisible(sys::exec_wait(
        paste0(ResultDir, "/","./plink"),
        args = c(
          "--bfile",
          paste0(DataDir, "/",finput),
          "--exclude",
          paste0(ResultDir, "/","study_SNP"),
          "--make-bed",
          "--out",
          paste0(ResultDir, "/","filtered_study_temp1"),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE

      ))

      ref <-
        read.table(file = paste0(ResultDir, "/",reference, ".bim"),
                   stringsAsFactors = FALSE)
      ref_AT <-
        ref[which(ref[, 5] == "A" & ref[, 6] == "T"), 2, drop = FALSE]
      ref_TA <-
        ref[which(ref[, 5] == "T" & ref[, 6] == "A"), 2, drop = FALSE]
      ref_GC <-
        ref[which(ref[, 5] == "G" & ref[, 6] == "C"), 2, drop = FALSE]
      ref_CG <-
        ref[which(ref[, 5] == "C" & ref[, 6] == "G"), 2, drop = FALSE]

      ref_SNP <- rbind(ref_AT, ref_TA, ref_GC, ref_CG)

      write.table(
        ref_SNP,
        file = paste0(ResultDir, "/","ref_SNP"),
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE,
        eol = "\r\n",
        sep = " "
      )
      invisible(sys::exec_wait(
        paste0(ResultDir, "/","./plink"),
        args = c(
          "--bfile",
          paste0(ResultDir, "/",reference),
          "--exclude",
          paste0(ResultDir, "/","ref_SNP"),
          "--make-bed",
          "--out",
          paste0(ResultDir, "/","filtered_ref_temp1"),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE

      ))

      ## LD pruning
      if (studyLD == TRUE) {
        invisible(sys::exec_wait(
          paste0(ResultDir, "/","./plink"),
          args = c(
            "--bfile",
            paste0(ResultDir, "/","filtered_study_temp1"),
            "--exclude",
            "range",
            paste0(ResultDir, "/",highLD_regions),
            "--indep-pairwise",
            studyLD_window_size,
            studyLD_step_size,
            studyLD_r2_threshold,
            "--make-bed",
            "--out",
            paste0(ResultDir, "/","filtered_study_temp2"),
            "--silent"
          ),
          std_out = FALSE,
          std_err = FALSE
        ))

      } else{
        invisible(sys::exec_wait(
          paste0(ResultDir, "/","./plink"),
          args = c(
            "--bfile",
            paste0(ResultDir, "/","filtered_study_temp1"),
            "--make-bed",
            "--out",
            paste0(ResultDir, "/","filtered_study_temp2"),
            "--silent"
          ),
          std_out = FALSE,
          std_err = FALSE
        ))
      }

      if (referLD == TRUE) {
        invisible(sys::exec_wait(
          paste0(ResultDir, "/","./plink"),
          args = c(
            "--bfile",
            paste0(ResultDir, "/","filtered_ref_temp1"),
            "--exclude",
            "range",
            paste0(ResultDir, "/",highLD_regions),
            "--indep-pairwise",
            referLD_window_size,
            referLD_step_size,
            referLD_r2_threshold,
            "--make-bed",
            "--out",
            paste0(ResultDir, "/","filtered_ref_temp2"),
            "--silent"
          ),
          std_out = FALSE,
          std_err = FALSE
        ))

      } else{
        invisible(sys::exec_wait(
          paste0(ResultDir, "/","./plink"),
          args = c(
            "--bfile",
            paste0(ResultDir, "/","filtered_ref_temp1"),
            "--make-bed",
            "--out",
            paste0(ResultDir, "/","filtered_ref_temp2"),
            "--silent"
          ),
          std_out = FALSE,
          std_err = FALSE
        ))
      }

      ftemp <- list.files(paste0(ResultDir,"/"),pattern = "filtered_study_temp1")
      invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))
      ftemp <- list.files(paste0(ResultDir,"/"),pattern = "filtered_ref_temp1")
      invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))

    } else if (filterSNP == FALSE) {
      invisible(sys::exec_wait(
        paste0(ResultDir, "/","./plink"),
        args = c(
          "--bfile",
          paste0(DataDir, "/",finput),
          "--make-bed",
          "--out",
          paste0(ResultDir, "/","filtered_study_temp2"),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))

      invisible(sys::exec_wait(
        paste0(ResultDir, "/","./plink"),
        args = c(
          "--bfile",
          paste0(ResultDir, "/",reference),
          "--make-bed",
          "--out",
          paste0(ResultDir, "/","filtered_ref_temp2"),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))

    }

    ## common SNPs between study and reference data. This step will reduce the file size much, hence recommended now.
    pruned_study <-
      read.table(file = paste0(ResultDir, "/","filtered_study_temp2", ".bim"),
                 stringsAsFactors = FALSE)
    pruned_ref <-
      read.table(file = paste0(ResultDir, "/","filtered_ref_temp2", ".bim"),
                 stringsAsFactors = FALSE)
    common_snps <- intersect(pruned_study$V2, pruned_ref$V2)

    if (length(common_snps) == 0){
      print("There is no common SNPs found between study and reference data. This analysis cannot be done.")
      return()
    }

    write.table(
      common_snps,
      file = paste0(ResultDir, "/","common_snps"),
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      eol = "\r\n",
      sep = " "
    )

    invisible(sys::exec_wait(
      paste0(ResultDir, "/","./plink"),
      args = c("--bfile",
               paste0(ResultDir, "/","filtered_study_temp2"),
               "--extract",
               paste0(ResultDir, "/","common_snps"),
               "--make-bed",
               "--out",
               paste0(ResultDir, "/","filtered_study_temp3")
      ),
      std_out = FALSE,
      std_err = FALSE
    ))
    invisible(sys::exec_wait(
      paste0(ResultDir, "/","./plink"),
      args = c("--bfile",
               paste0(ResultDir, "/","filtered_ref_temp2"),
               "--extract",
               paste0(ResultDir, "/","common_snps"),
               "--make-bed",
               "--out",
               paste0(ResultDir, "/","filtered_ref_temp3")
      ),
      std_out = FALSE,
      std_err = FALSE
    ))

    ftemp <- list.files(paste0(ResultDir,"/"),pattern = "filtered_study_temp2")
    suppressWarnings(invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp)))))
    ftemp <- list.files(paste0(ResultDir,"/"),pattern = "filtered_ref_temp2")
    suppressWarnings(invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp)))))

    # Correcting the chromosome mismatch between study and reference data
    S1 <-
      pruned_study[match(common_snps, pruned_study[, 2]), , drop = FALSE]
    S2 <-
      pruned_ref[match(common_snps, pruned_ref[, "V2"]), , drop = FALSE]
    # check chr and bp position consistency
    whChrNotSame <- which(S1[, 1] != S2[, "V1"])
    whPosNotSame <- which(S1[, 4] != S2[, "V4"])
    whChrPosNotSame <- union(whChrNotSame, whPosNotSame)

    snpSameNameDiffPos <- S1[whChrPosNotSame, c(1, 4, 2), drop = FALSE]


    if (nrow(snpSameNameDiffPos) != 0) {
      write.table(
        snpSameNameDiffPos,
        file = paste0(ResultDir, "/","snpSameNameDiffPos"),
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE,
        eol = "\r\n",
        sep = " "
      )

      invisible(sys::exec_wait(
        paste0(ResultDir, "/","./plink"),
        args = c("--bfile",
                 paste0(ResultDir, "/","filtered_ref_temp3"),
                 "--update-chr",
                 paste0(ResultDir, "/","snpSameNameDiffPos"),
                 1,
                 3,
                 "--update-map",
                 paste0(ResultDir, "/","snpSameNameDiffPos"),
                 2,
                 3,
                 "--make-bed",
                 "--out",
                 paste0(ResultDir, "/","filtered_ref_temp4")
        ),
        std_out = FALSE,
        std_err = FALSE
      ))

    } else{
      invisible(sys::exec_wait(
        paste0(ResultDir, "/","./plink"),
        args = c(
          "--bfile",
          paste0(ResultDir, "/","filtered_ref_temp3"),
          "--make-bed",
          "--out",
          paste0(ResultDir, "/","filtered_ref_temp4"),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))

    }


    ftemp <- list.files(paste0(ResultDir,"/"),pattern = "filtered_ref_temp3")
    suppressWarnings(invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp)))))

    # Finding mis-matching allele positions
    updated_ref <-
      read.table(file = paste0(ResultDir, "/","filtered_ref_temp4", ".bim"),
                 stringsAsFactors = FALSE)
    S3 <-
      updated_ref[match(common_snps, updated_ref[, "V2"]), , drop = FALSE]
    colnames(S1) <- c("V1", "V2", "V3", "V4", "Sa", "Sb")
    colnames(S3) <- c("V1", "V2", "V3", "V4", "Ra", "Rb")

    #library(data.table)
    S1 <- data.table::as.data.table(S1)
    S3 <- data.table::as.data.table(S3)
    #S4 <- merge(S1, S3, by = c("V1", "V2", "V4"))
    #Using SNP name and chr no. for merging, not using base-pair position
    S4 <- merge(S1, S3, by = c("V1", "V2"))
    snps_flips <- S4[which(S4[, 5] != S4[, 9] & S4[, 6] != S4[, 10]), ]
    snp_allele_flips <- unique(snps_flips$V2)

    if (length(unique(snp_allele_flips)) != 0) {

      write.table(
        snp_allele_flips,
        file = paste0(ResultDir, "/","snp_allele_flips"),
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE,
        eol = "\r\n",
        sep = " "
      )

      invisible(sys::exec_wait(
        paste0(ResultDir, "/","./plink"),
        args = c("--bfile",
                 paste0(ResultDir, "/","filtered_ref_temp4"),
                 "--flip",
                 paste0(ResultDir, "/","snp_allele_flips"),
                 "--make-bed",
                 "--out",
                 paste0(ResultDir, "/","filtered_ref_temp5")
        ),
        std_out = FALSE,
        std_err = FALSE
      ))


    } else{
      invisible(sys::exec_wait(
        paste0(ResultDir, "/","./plink"),
        args = c(
          "--bfile",
          paste0(ResultDir, "/","filtered_ref_temp4"),
          "--make-bed",
          "--out",
          paste0(ResultDir, "/","filtered_ref_temp5"),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))

    }


    ftemp <- list.files(paste0(ResultDir,"/"),pattern = "filtered_ref_temp4")
    suppressWarnings(invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp)))))

    # Checking allele flips again after correcting
    flipped_ref <-
      read.table(file = paste0(ResultDir, "/","filtered_ref_temp5", ".bim"),
                 stringsAsFactors = FALSE)
    S5 <-
      flipped_ref[match(common_snps, flipped_ref[, "V2"]), , drop = FALSE]
    colnames(S5) <- c("V1", "V2", "V3", "V4", "Ra", "Rb")
    S5 <- data.table::as.data.table(S5)
    # S6 <- merge(S1, S5, by = c("V1", "V2", "V4"))
    # snps_flips_wrong <- S6[which(S6[, 5] != S6[, 8] &
    #                                S6[, 6] != S6[, 9]), ]
    S6 <- merge(S1, S5, by = c("V1", "V2"))
    snps_flips_wrong <- S6[which(S6[, 5] != S6[, 9] &
                                   S6[, 6] != S6[, 10]), ]
    allele_flips_wrong <- unique(snps_flips_wrong$V2)

    if (length(unique(allele_flips_wrong)) != 0) {
      write.table(
        allele_flips_wrong,
        file = paste0(ResultDir, "/","allele_flips_wrong"),
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE,
        eol = "\r\n",
        sep = " "
      )

      invisible(sys::exec_wait(
        paste0(ResultDir, "/","./plink"),
        args = c("--bfile",
                 paste0(ResultDir, "/","filtered_ref_temp5"),
                 "--exclude",
                 paste0(ResultDir, "/","allele_flips_wrong"),
                 "--make-bed",
                 "--out",
                 paste0(ResultDir, "/","filtered_ref_temp6")
        ),
        std_out = FALSE,
        std_err = FALSE
      ))


    } else{
      invisible(sys::exec_wait(
        paste0(ResultDir, "/","./plink"),
        args = c(
          "--bfile",
          paste0(ResultDir, "/","filtered_ref_temp5"),
          "--make-bed",
          "--out",
          paste0(ResultDir, "/","filtered_ref_temp6"),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))

    }

    ftemp <- list.files(paste0(ResultDir,"/"),pattern = "filtered_ref_temp5")
    suppressWarnings(invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp)))))

    # Checking number of SNPs in reference after clean-up.
    cleaned_ref <-
      read.table(file = paste0(ResultDir, "/","filtered_ref_temp6", ".bim"),
                 stringsAsFactors = FALSE)
    snps_final <- unique(cleaned_ref$V2)

    # Merge study and reference datasets

    invisible(sys::exec_wait(
      paste0(ResultDir, "/","./plink"),
      args = c(
        "--bfile",
        paste0(ResultDir, "/","filtered_study_temp3"),
        "--bmerge",
        paste0(ResultDir, "/","filtered_ref_temp6"),
        "--make-bed",
        "--out",
        paste0(ResultDir, "/","study_ref_merge"),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))

    ## Checking for strand inconsistency
    if (file.exists(paste0(ResultDir, "/study_ref_merge-merge.missnp"))){

      # Merge study and reference datasets
      invisible(sys::exec_wait(
        paste0(ResultDir, "/","./plink"),
        args = c("--bfile",
                 paste0(ResultDir, "/","filtered_ref_temp6"),
                 "--exclude",
                 paste0(ResultDir, "/","study_ref_merge-merge.missnp"),
                 "--make-bed",
                 "--out",
                 paste0(ResultDir, "/","filtered_ref_temp7")
        ),
        std_out = FALSE,
        std_err = FALSE
      ))

      invisible(sys::exec_wait(
        paste0(ResultDir, "/","./plink"),
        args = c(
          "--bfile",
          paste0(ResultDir, "/","filtered_study_temp3"),
          "--bmerge",
          paste0(ResultDir, "/","filtered_ref_temp7"),
          "--make-bed",
          "--out",
          paste0(ResultDir, "/","study_ref_merge"),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))
    }

    ftemp <- list.files(paste0(ResultDir,"/"),pattern = "filtered_study_temp3")
    suppressWarnings(invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp)))))
    ####check
    # d1 <- read.table(paste0(ResultDir,"/filtered_study_temp3.bim"))
    # x1 <- d1[nchar(d1[,5]) > 1 | nchar(d1[,6]) > 1, ]
    #
    # d2 <- read.table(paste0(ResultDir,"/filtered_ref_temp6.bim"))
    # x1 <- d2[nchar(d2[,5]) > 1 | nchar(d2[,6]) > 1, ]
    ####

    # PCA for ancestry: family and individual ID in columns 1 and 2, followed by the first 20 principal components.
    invisible(sys::exec_wait(
      paste0(ResultDir, "/","./plink"),
      args = c(
        "--bfile",
        paste0(ResultDir, "/","study_ref_merge"),
        "--pca",
        "--out",
        paste0(ResultDir, "/","study_ref_merge"),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))


    ftemp <- list.files(paste0(ResultDir,"/"),pattern = "filtered_ref_temp6")
    suppressWarnings(invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp)))))
    ftemp <- list.files(paste0(ResultDir,"/"),pattern = "filtered_ref_temp7")
    suppressWarnings(invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp)))))

    print("PCA done.")

    if (reference == "HapMapIII_NCBI36") {


      ref_ancestry <-
        read.table(file = paste0(system.file("extdata", package = "GXwasR"), "/","hapmap_relationships_w_pops_121708.txt"),
                   stringsAsFactors = FALSE,
                   header = TRUE)[, c("IID", "population")]
      colnames(ref_ancestry) <- c("ID", "Ancestry")
      ref_ancestry[which(ref_ancestry$Ancestry == "CHD"), "Ancestry"] <-
        "ASIAN"
      ref_ancestry[which(ref_ancestry$Ancestry == "CHB"), "Ancestry"] <-
        "ASIAN"
      ref_ancestry[which(ref_ancestry$Ancestry == "JPT"), "Ancestry"] <-
        "ASIAN"
      ref_ancestry[which(ref_ancestry$Ancestry == "GIH"), "Ancestry"] <-
        "ASIAN"
      ref_ancestry_CEU_YRI_ASIAN <-
        ref_ancestry[ref_ancestry$Ancestry == "CEU" |
                       ref_ancestry$Ancestry == "YRI" | ref_ancestry$Ancestry == "ASIAN", ]
      ref_ancestry_CEU_YRI_ASIAN[ref_ancestry_CEU_YRI_ASIAN$Ancestry == "CEU",2]<- "Ref_CEU"
      ref_ancestry_CEU_YRI_ASIAN[ref_ancestry_CEU_YRI_ASIAN$Ancestry == "YRI",2]<- "Ref_YRI"
      ref_ancestry_CEU_YRI_ASIAN[ref_ancestry_CEU_YRI_ASIAN$Ancestry == "ASIAN",2]<- "Ref_ASIAN"

    } else if (reference == "ThousandGenome") {

      ref_ancestry <-
        read.table(file = paste0(system.file("extdata", package = "GXwasR"), "/","1000genomesampleinfo.txt"),
                   stringsAsFactors = FALSE,
                   header = TRUE)[, c("Sample", "superpop")]
      ref_ancestry_CEU_YRI_ASIAN <- ref_ancestry
      ref_ancestry_CEU_YRI_ASIAN$Ancestry <- paste0("Ref_",ref_ancestry_CEU_YRI_ASIAN$superpop)
      colnames(ref_ancestry_CEU_YRI_ASIAN)<- c("ID", "Ancestry")
    }
    # study_ancestry <-
    #   read.table(
    #     file = paste0(DataDir, "/",study_pop),
    #     stringsAsFactors = FALSE,
    #     header = FALSE
    #   )
    study_ancestry <- study_pop
    colnames(study_ancestry) <- c("ID", "Ancestry")
    study_ancestry$Ancestry <-
      paste0("Study_", study_ancestry$Ancestry)

    pop <- rbind(ref_ancestry_CEU_YRI_ASIAN, study_ancestry)
    population = pop$Ancestry

    pca <-
      read.table(file = paste0(ResultDir, "/","study_ref_merge.eigenvec"),
                 stringsAsFactors = FALSE,
                 header = FALSE)
    tab <- data.frame(
      sample = pca$V2,
      pop = factor(population)[match(pca$V2, pop$ID)],
      PC1 = pca[, 3],
      # the first eigenvector
      PC2 = pca[, 4],
      # the second eigenvector
      stringsAsFactors = FALSE
    )
    tab <- na.omit(tab)
    #head(tab)
    pop_type <- as.data.frame(unique(tab$pop))
    colnames(pop_type) <- "type"
    pop_type$value <- 20
    pop_type[which(pop_type$type %like% "Study_.*"), "value"] <- 3
    #Setting reference population on top since labeling is alphabetical
    pop_type1 <- pop_type[which(pop_type$type %like% "Ref_.*"),]
    pop_type2 <- pop_type[which(pop_type$type %like% "Study_.*"),]
    pop_type <- rbind(pop_type1,pop_type2)

    #library(ggplot2)
    p <-
      ggplot2::ggplot(data = tab, ggplot2::aes(
        x = PC1,
        y = PC2,
        color = pop,
        shape = pop
      )) +

      ggplot2::geom_hline(yintercept = 0, lty = 2) +

      ggplot2::geom_vline(xintercept = 0, lty = 2) +

      ggplot2::guides(color = ggplot2::guide_legend(title = "Ancestry"),
                      shape = guide_legend(title = "Ancestry")) + ggplot2::scale_shape_manual(values = pop_type$value) +

      ggplot2::geom_point(alpha = 1, size = 4)

    print(p)

    if (outlier == TRUE) {
      # PC1 <- tab[tab$pop == "Ref_CEU", 3]
      # PC2 <- tab[tab$pop == "Ref_CEU", 4]
      # Eu <- tab[tab$pop == "Ref_CEU", ]
      PC1 <- tab[tab$pop == paste0("Ref_",outlierOf), 3]
      PC2 <- tab[tab$pop == paste0("Ref_",outlierOf), 4]
      Eu <- tab[tab$pop == paste0("Ref_",outlierOf), ]

      Eu_pc1_median <- median(PC1)
      Eu_pc2_median <- median(PC2)

      Eu$dis <- sqrt((Eu$PC1 - Eu_pc1_median) ^ 2 +
                       (Eu$PC2 - Eu_pc2_median) ^ 2)

      Eu_max <- max(Eu$dis)

      study <- tab[tab$pop %like% "Study_.*", ]
      study$dis <- sqrt((study$PC1 - Eu_pc1_median) ^ 2 +
                          (study$PC2 - Eu_pc2_median) ^ 2)

      Non_Eu <- study[study$dis > Eu_max * outlier_threshold, ]
      Outlier_samples <- unique(Non_Eu$sample)
      print("Ancestry check is performed.")

      if (length(Outlier_samples) != 0) {
        write.table(
          Outlier_samples,
          file = paste0(ResultDir, "/","Outlier_ancestry"),
          quote = FALSE,
          row.names = FALSE,
          col.names = FALSE,
          eol = "\r\n",
          sep = " "
        )
      }

    } else{
      Outlier_samples <- 0
    }


    ###############

    if (filterSNP == TRUE) {
      if (nrow(unique(study_SNP)) == 0){
        print(paste0(
          nrow(unique(study_SNP)),
          " No SNP had 'A-T' and 'G-C' in study data."
        ))
      }else{
        print(paste0(
          nrow(unique(study_SNP)),
          " SNPs had 'A-T' and 'G-C' in study data. These SNPs were removed."
        ))}

      if ( nrow(unique(ref_SNP)) == 0){
        print(
          " No SNP had 'A-T' and 'G-C' in reference data."
        )
      }else{
        print(
          paste0(
            nrow(unique(ref_SNP)),
            " SNPs were 'A-T' and 'G-C' in reference data. These SNPs were removed."
          )
        )}

    } else {
      print(
        "A-T and C-G SNPs recommended to remove from both the reference and study data set by setting filterSNP == TRUE"
      )
    }

    if (studyLD == TRUE) {
      print("LD prunning was done for study dataset.")
    } else {
      print("LD prunning for study dataset is recommended. Set studyLD == FALSE.")
    }

    if (referLD == TRUE) {
      print("LD prunning was done for reference dataset.")
    } else {
      print("LD prunning for reference dataset is recommended. Set referLD == TRUE.")
    }

    # if (nrow(snpSameNameDiffPos) != 0) {
    #   print(paste0(
    #     nrow(snpSameNameDiffPos)," SNPs that had the same name but different genomic position were updated in the reference data."
    #   ))
    #}
    if (length(unique(snp_allele_flips)) == 0) {
      print("No allele flips between study and reference data.")
    } else{
      print(paste0(
        length(unique(snp_allele_flips)),
        " allele flips identified between study and reference data."
      ))
      mergebim <- read.table(paste0(ResultDir, "/","study_ref_merge.bim"))
      print(
        paste0(
          length(unique(mergebim$V2)),
          " SNPs were finally retained in study and reference data after correcting for position mismatch and allele flips."
        )
      )
    }
    if (outlier == TRUE & length(Outlier_samples) != 0) {
      print(
        paste0(
          length(Outlier_samples),
          " Non-European samples as outlier."
        )
      )

    } else if (outlier == TRUE & length(Outlier_samples) == 0) {
      print("There is no non-European sample as outlier.")

    } else if (outlier == FALSE) {
      print("Non-European detection as outlier was not performed.")
    }

    ftemp <- list.files(paste0(ResultDir,"/"),pattern = "study_ref_merge")
    suppressWarnings(invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp)))))

    ftemp <- c("study_SNP","ref_SNP","common_snps","snp_allele_flips","allele_flips_wrong")
    suppressWarnings(invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp)))))

    if (file.exists(paste0(ResultDir,"/Outlier_ancestry"))) {
      ftemp <- c("Outlier_ancestry")
      suppressWarnings(invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp)))))
    }

    if (file.exists(paste0(ResultDir,"/snpSameNameDiffPos"))) {
      ftemp <- c("snpSameNameDiffPos")
      suppressWarnings(invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp)))))
    }

    ftemp <- list.files(paste0(ResultDir,"/"),pattern = paste0(reference))
    suppressWarnings(invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp)))))

    Outlier_samples1 <- as.data.frame(Non_Eu[,1:2])
    Outlier_samples1$pop <- as.character(Outlier_samples1$pop)
    Outlier_samples1 <- Outlier_samples1[Outlier_samples1$pop %like% "Study_.*",]
    return(Outlier_samples1)
  }


#' TestXGene: Performing gene-based association test using GWAS/XWAS summary statistics.
#'
#' @description This function performs gene-based association tests using GWAS/XWAS summary statistics and SNP-SNP correlation matrices. For  SNP-SNP correlation matrices, users have the flexibility to use either the base genotype data or 1000 Genomes Phase 3 reference genotype data. Users also have options to define the regional positions of genes to include the SNPs according to their investigation. This function computes gene-wise SNP-SNP correlation matrices and can perform nine different gene-based tests, such as, “BT" (burden test), "SKAT" (sequence kernel association test), "SKATO" (combination of BT and SKAT), "sumchi" (sum of χ2-statistics), "ACAT" (aggregated Cauchy association test for combining P values), "PCA"(principal component approach), "FLM"( functional multiple linear regression model), "simpleM" (Bonferroni correction test), "minp" (minimum P-value) leveraging PLINK1.9 (1) and sumFREGAT (2,3) tools.
#'
#' Though this function implicitly performs X-linked gene-based test, it is flexible to perform this analysis genome-wide. For the details about the different tests, please follow the associated paper.
#'
#' @param DataDir A character string for the file path of the all the input files.
#' @param ResultDir A character string for the file path where all output files will be stored. The default is tempdir().
#' @param finput Character string, specifying the prefix of the input PLINK binary files for the genotype data. This file is used to compute the correlation between the SNPs. This file needs to be in DataDir. If the base genotype data is unavailable, then users can use the 1000 Genomes Project samples. Users should use the population that most closely represents the base sample. For ACAT model, this parameter is not mandatory and could be set NULL.
#' @param sumstat A dataframe object with GWAS summary statistics. When the base-genotype data is used to compute genetic correlations, the mandatory columns are Column 1: "CHROM" (i.e., chromosome code), Column 2: "POS" (i.e., base-pair position), Column 3: "ID" (i.e. SNP IDs), Column 4: “P” (i.e., p-values), Column 5: “BETA” (i.e., effect-size), Column 6: "A1" (i.e., effect allele), Column 7: “A2” (i.e., alternative allele) and Column 8: "EAF" (i.e., the effect allele frequency) are mandatory when base-genotype data is used to compute genetic correlations. Otherwise, if the users are using reference data, then columns 5 to 8 are optional. Also, in that case, columns, such as "REF" (i.e., reference allele), and "ALT" (i.e., alternative allele) could be present to compare alleles with those in the reference file and exclude genetic variants if alleles do not match. There could be an additional column, “ANNO” with functional annotations (like "intron_variant", "synonymous", "missense" etc.)
#' @param gene_file Character string, specifying the prefix of the name of a .txt file listing genes in refFlat format. This file needs to be in DataDir. The X-linked gene files, "Xlinkedgenes_hg19.txt" and "Xlinkedgenes_hg38.txt" and autosomal gene files, “Autosomes_hg19.txt” and “Autosomes_hg38.txt” can be specified. The default is "Xlinkedgenes_hg19.txt". The genome built should be in agreement with the analysis.
#' @param gene_range Integer value, specifying the up stream and down stream range (in kilo base) of a gene for SNPs to be considered. The default is 500000.
#' @param score_file Character string, specifying the prefix of a file which will be used to produce score files with Z scores from P values and beta input from GWAS summary statistics.
#' @param ref_data Character string, specifying the path to a reference dataframe with additional data needed to recode user data according to correlation matrices that will be used. It contains "ID" column with names of  SNPs,  "REF" and "ALT" columns with alleles that were coded as 0 and 1, respectively. Effect sizes from data will be inverted for variants with effect alleles different from "ALT" alleles in reference data. If presented, "REF" and "ALT" columns from the input data will be used to sort out variants with alleles different from those in reference data. This dataframe  can also be a source of map data and allele frequencies if they are not present in data. "AF" column in the reference file represents the allele frequency of "ALT" allele. The default is "ref1KG.MAC5.EUR_AF.RData".
#' @param max_gene Positive integer value, specifying the number of genes for which the gene-based test will be performed. The default is NULL to consider all the genes.
#' @param sample_size Positive integer value, specifying the sample size of the GWAS. Only needed for FLM and PCA models.
#' @param genebasedTest Character string, specifying the name of the gene-based test. Nine different tests can be specified, "SKAT","SKATO","sumchi","ACAT","BT","PCA","FLM","simpleM","minp". The default is "SKAT".
#' @param beta_par Boolean value, 'TRUE' or 'FALSE', specifying whether approximation for large genes (>= 500 SNPs) should be used. Applicable for SKAT, SKATO, sumchi, PCA, FLM (default = TRUE for these methods).
#' @param weights_function A function of MAF to assign weights for each genetic variant. By default is NULL. In this case the weights will be calculated using the beta distribution.
#' @param geno_variance_weights Character string, indicating whether scores should be weighted by the variance of genotypes: "none" (i.e., no weights applied, resulting in a sum chi-square test); "se.beta" (i.e., scores weighted by variance of genotypes estimated from P values and effect sizes); "af" (i.e., scores weighted by variance of genotypes calculated as AF * (1 - AF), where AF is allele frequency.
#' @param kernel_p_method Character string, specifying the method for computing P value in kernel-based tests, such as SKAT, SKATO and sumchi. Available methods are "kuonen" (9), "davies" (10) and "hybrid" (11). The default is "kuonen".
#' @param acc_devies Positive numeric value, specifying the accuracy parameter for "davies" method. The default is 1e-8.
#' @param lim_devies Positive numeric value, specifying the limit parameter for "davies" method. The default is 1e+6.
#' @param rho Logical value, 'TRUE' or 'FALSE' or can be a vector of grid values from 0 to 1. If TRUE, the optimal test (SKAT-O) is performed (12). The default grid is c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1).
#' @param skato_p_threshold Positive numeric value, specifying the largest P value that will be considered as important when performing computational optimization in SKAT-O. All P values larger than skato_p_threshold will be processed via burden test. The default is 0.8
#' @param anno_type A character (or character vector) indicating annotation types to be used. The default is "" (i.e, nothing).
#' @param mac_threshold Integer value, specifying the threshold of MACs (Minor allele content) calculated from MAFs. In ACAT, scores with MAC <= 10 will be combined using Burden test.
#' @param regularize_fun Character string, specifying the one of two regularization algorithms if ‘reference_matrix’ is TRUE:  'LH' (default) or 'derivLH'. Currently, both give similar results.
#' @param pca_var_fraction Positive numeric value, specifying the minimal proportion of genetic variance within the region that should be explained by principal components used in PCA method. This is also valid in 'simpleM'. The default is 0.85.
#' @param flm_basis_function Character string, specifying the name of a basis function type for beta-smooth in FLM method. Can be set to "bspline" (B-spline basis) or "fourier" (Fourier basis, default).
#' @param flm_num_basis Positive integer value, specifying the number of basis functions to be used for beta-smooth in FLM method. The default is 25.
#' @param flm_poly_order Positive integer value, specifying the polynomial order to be used in "bspline" for FLM model. The default = 4, which corresponds to the cubic B-splines. This has no effect if only Fourier bases are used
#' @param flip_genotypes Logical value, 'TRUE' or 'FALSE', indicating whether the genotypes of some genetic variants should be flipped (relabelled) for their better functional representation (13). The default is FALSE.
#' @param omit_linear_variant Logical value, 'TRUE' or 'FALSE', indicating whether to omit linearly dependent genetic variants. It was done in the FLM test (4). The default is FALSE.
#' @param gene_approximation Boolean value, 'TRUE' or 'FALSE', specifying whether approximation for large genes (>= 500 SNPs) should be used. Applicable for SKAT, SKATO, sumchi, PCA, FLM. The default is TRUE for these methods).
#' @param reference_matrix_used Boolean value, 'TRUE' or'FALSE' logical indicating whether the correlation matrices were generated using the reference matrix. The default is FALSE. If  TRUE, regularization algorithms will be applied to ensure the invertibility and numerical stability of the matrices.
#'
#' @returns A data frame with columns "gene",  "chrom", "start", "end", "markers" (i.e., numbers of SNPs), "filtered.markers" (i.e. filtered SNPs) and "pvalue" (i.e., p-value). Additionally, for  “BT”, there will be “beta” (i.e., gene-level estimates of betas) and “beta.se” (i.e., standard errors of betas).  For “FLM”, there will be the “model” column with the names of the functional models used for each region. Names shortly describe the functional basis and the number of basis functions used. E.g., "F25" means 25 Fourier basis functions, "B15" means 15 B-spline basis functions. For “PCA”, there will be the “ncomponents” (the number of the principal components used for each region) and “explained.variance.fraction”(i.e., the proportion of genetic variance they make up) columns.
#'
#' @references
#' (1) Purcell et. al.,(2007). PLINK: A Tool Set for Whole-Genome Association and Population-Based Linkage Analyses. The American Journal of Human Genetics, 81(3), 559-575. https://doi.org/10.1086/519795.
#' (2) Svishcheva GR, Belonogova NM, Zorkoltseva IV, Kirichenko AV, Axenovich TI. Gene-based association tests using GWAS summary statistics. Bioinformatics. 2019 Oct 1;35(19):3701-3708. doi: 10.1093/bioinformatics/btz172. PMID: 30860568.
#' (3) Belonogova NM, Svishcheva GR, Kirichenko AV, Zorkoltseva IV, Tsepilov YA, Axenovich TI. sumSTAAR: A flexible framework for gene-based association studies using GWAS summary statistics. PLoS Comput Biol. 2022 Jun 2;18(6):e1010172. doi: 10.1371/journal.pcbi.1010172. PMID: 35653402; PMCID: PMC9197066.
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom regioneR toGRanges
#' @importFrom plyranges join_overlap_intersect
#' @importFrom sumFREGAT SKAT SKATO sumchi ACAT BT PCA FLM simpleM minp
#'
#' @export
#'
#' @examples
#'#Not Run
# DataDir <- system.file("extdata", package = "GXwasR")
# ResultDir <- tempdir()
# #DataDir <- "inst/extdata"
# #ResultDir <- "Data"
# #finput <- "ExamplePLINK"
# finput <- "GXwasR_example"
# data("GXwasRData")
# #load("~/b1137/BBose/GXwasR/inst/extdata/XWAS_Summary.Rda")
# #load("~/b1137/BBose/GXwasR/inst/extdata/XWAS_Summary_Example.Rda")
# #load("~/b1137/BBose/GXwasR/data/XWAS_Summary_Example.Rda")
# #sumstat <- XWAS_Summary
# sumstat <- XWAS_Summary_Example
# ref_data <- NULL
# gene_file <- "Xlinkedgenes_hg19.txt"
# gene_range <- 500000
# max_gene <- 10
# gene_approximation <- TRUE
# beta_par <- c(1, 25)
# weights_function <- NULL
# geno_variance_weights <- "se.beta"
# method = "kuonen"
# acc_devies = 1e-8
# lim_devies = 1e+6
# rho = TRUE
# skato_p_threshold = 0.8
# mac_threshold <- 3
# sample_size <- 4000
# reference_matrix_used <- FALSE
# regularize_fun <- "LH"
# pca_var_fraction <- 0.85
# flm_basis_function <- "fourier"
# flm_num_basis <- 25
# flm_poly_order <- 4
# flip_genotypes <- FALSE
# omit_linear_variant <- FALSE
# kernel_p_method <- "kuonen"
# anno_type <- ""
#
# GenetestResult <- TestXGene(
#   DataDir,
#   ResultDir,
#   finput,
#   sumstat,
#   gene_file,
#   gene_range,
#   score_file,
#   ref_data,
#   max_gene,
#   sample_size,
#   genebasedTest = "SKAT",
#   gene_approximation,
#   beta_par,
#   weights_function,
#   geno_variance_weights,
#   kernel_p_method,
#   acc_devies,
#   lim_devies,
#   rho,
#   skato_p_threshold,
#   anno_type,
#   mac_threshold,
#   reference_matrix_used,
#   regularize_fun,
#   pca_var_fraction,
#   flm_basis_function,
#   flm_num_basis,
#   flm_poly_order,
#   flip_genotypes,
#   omit_linear_variant
# )

TestXGene <-
  function(DataDir,
           ResultDir = tempdir(),
           finput,
           sumstat,
           gene_file,
           gene_range = 500000,
           score_file,
           ref_data = NULL,
           max_gene = NULL,
           sample_size = NULL,
           genebasedTest = c("SKAT",
                             "SKATO",
                             "sumchi",
                             "ACAT",
                             "BT",
                             "PCA",
                             "FLM",
                             "simpleM",
                             "minp"),
           gene_approximation = TRUE,
           beta_par,
           weights_function,
           geno_variance_weights,
           kernel_p_method = "kuonen",
           acc_devies = 1e-8,
           lim_devies = 1e+6,
           rho = TRUE,
           skato_p_threshold = 0.8,
           anno_type ="",
           mac_threshold,
           reference_matrix_used,
           regularize_fun,
           pca_var_fraction = 0.85,
           flm_basis_function = "fourier",
           flm_num_basis = 25,
           flm_poly_order = 4,
           flip_genotypes = FALSE,
           omit_linear_variant = FALSE) {

    if (file.exists(paste0(DataDir, "/", finput, ".bed")) &&
        file.exists(paste0(DataDir, "/", finput, ".bim")) &&
        file.exists(paste0(DataDir, "/", finput, ".fam"))) {

      setupPlink(ResultDir)

    } else{
      writeLines(
        "There are no Plink files in DataDir.\nPlease specify correct DataDir path with input Plink files."
      )
    }

    input.dat <- sumstat[, c("CHROM", "POS", "ID", "A1", "P", "BETA", "EAF")]
    colnames(input.dat) <- c("CHROM", "POS", "ID", "EA", "P", "BETA", "EAF")
    ref.data <- sumstat[, c("CHROM", "POS", "ID", "A2", "A1", "EAF")]
    colnames(ref.data) <-
      c("CHROM", "POS", "ID", "REF", "ALT", "AF") ## Following convention for reference data

    if (is.null(ref_data)) {
      ref_data <- ref.data
    } else{
      ref_data <- ref_data
    }

    #library(sumFREGAT)
    #library(GenomicRanges)
    #library(plyranges)
    #library(regioneR)

    geneTestScoreFile(ResultDir = ResultDir,data = input.dat,
                      reference = ref_data,
                      output.file.prefix = "gene.test.score.file") ## use suppressWarnings()

    genes <- read.table(paste0(DataDir,"/",gene_file))
    colnames(genes) <- c(c("gene_name", "X", "chr", "Y", "start", "end"))
    genes$up_Mb <- genes$start - gene_range
    genes$down_Mb <- genes$end + gene_range
    genes.gr <- GenomicRanges::makeGRangesFromDataFrame(genes, keep.extra.columns = T)

    suppressWarnings(SNPfile <- read.table(
      file = paste0(DataDir,"/",finput, ".bim"),
      header = FALSE,
      #na = "NA",
      na.strings = "NA"
    ))

    SNPfile$chr <- SNPfile$V1
    SNPfile$start <- SNPfile$V4
    SNPfile$end <- SNPfile$V4
    SNPfile$SNP <- SNPfile$V2
    snp_data <- SNPfile[, c(7, 8, 9, 10)]
    snp.gr <- regioneR::toGRanges(snp_data)
    gene_snp_intersect <-
      as.data.frame(plyranges::join_overlap_intersect(genes.gr, snp.gr))
    print(paste0(length(unique(
      gene_snp_intersect$gene_name
    )), " genes are having ", length(unique(
      gene_snp_intersect$SNP
    )), " SNPs"))
    gene_snp <- unique(gene_snp_intersect[, c(6, 11)])
    snpcount <- as.data.frame(table(gene_snp$gene_name))

    g <- as.character(snpcount[, 1])
    dir.create(path = paste0(ResultDir,"/cormatrix"))
    print("SNP-SNP correlation matrices are being created.")

    snpcorrFun <- function(g) {
      #g <- gene_snp$gene_name
      snps <- gene_snp[gene_snp$gene_name == g, 2, drop = FALSE]
      #print(g)
      write.table(
        snps,
        file = paste0(ResultDir,"/cor_snps.txt"),
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE,
        eol = "\r\n"
      )

      invisible(sys::exec_wait(
        paste0(ResultDir,"/","./plink"),
        args = c(
          "--bfile",
          paste0(DataDir,"/",finput),
          "--r2",
          "square",
          "--extract",
          paste0(ResultDir,"/cor_snps.txt"),
          "--out",
          paste0(ResultDir,"/snpcorr"),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))

      snpcorr <- as.matrix(read.table(file = paste0(ResultDir,"/snpcorr.ld")))
      colnames(snpcorr) <- snps$SNP
      rownames(snpcorr) <- snps$SNP
      snpcorr <- as.data.frame(snpcorr)
      save(snpcorr, file = paste0(ResultDir,"/cormatrix/", g, ".RData"))
      return()
    }

    invisible(lapply(g, snpcorrFun))
    print("SNP-SNP correlation matrices are done.")

    ## TESTS
    #score.file <- paste0(score_file, ".vcf.gz")
    score.file <- paste0(ResultDir,"/gene.test.score.file.vcf.gz")
    gene.file <- gene_file

    print("line 126")

    if (is.null(max_gene)) {
      genes1 <- as.vector(g)

    } else{
      maxgene <- max_gene
      genes1 <- as.vector(g[1:max_gene])

    }

    if (genebasedTest == "SKAT") {
      x <-
        sumFREGAT::SKAT(
          score.file = score.file,
          gene.file = paste0(DataDir,"/",gene_file),
          genes = genes1,
          cor.path = paste0(ResultDir,"/cormatrix"),
          approximation = gene_approximation,
          anno.type = anno_type,
          beta.par = beta_par,
          weights.function = weights_function,
          user.weights = FALSE,
          gen.var.weights = geno_variance_weights,
          method = kernel_p_method,
          acc = acc_devies,
          lim = lim_devies,
          rho = rho,
          p.threshold = skato_p_threshold,
          write.file = FALSE,
          quiet = TRUE
        )
      return(x)

    } else if (genebasedTest == "SKATO") {
      x <-
        sumFREGAT::SKATO(
          score.file = score.file,
          gene.file = paste0(DataDir,"/",gene_file),
          genes = genes1,
          cor.path = paste0(ResultDir,"/cormatrix"),
          anno.type = anno_type,
          approximation = gene_approximation,
          beta.par = beta_par,
          weights.function = weights_function,
          user.weights = FALSE,
          method = kernel_p_method,
          acc = acc_devies,
          lim = lim_devies,
          rho = rho,
          p.threshold = skato_p_threshold,
          write.file = FALSE,
          quiet = TRUE
        )
      return(x)

    } else if (genebasedTest == "sumchi") {
      x <-
        sumFREGAT::sumchi(
          score.file = score.file,
          gene.file = paste0(DataDir,"/",gene_file),
          genes = genes1,
          cor.path = paste0(ResultDir,"/cormatrix"),
          approximation = gene_approximation,
          anno.type = anno_type,
          method = kernel_p_method,
          acc = acc_devies,
          lim = lim_devies,
          write.file = FALSE,
          quiet = TRUE
        )
      return(x)

    } else if (genebasedTest == "ACAT") {
      x <-
        sumFREGAT::ACAT(
          score.file = score.file,
          gene.file = paste0(DataDir,"/",gene_file),
          genes = genes1,
          anno.type = anno_type,
          beta.par = beta_par,
          weights.function = weights_function,
          user.weights = FALSE,
          gen.var.weights = geno_variance_weights,
          mac.threshold = mac_threshold,
          n = sample_size,
          write.file = FALSE,
          quiet = TRUE
        )
      return(x)

    } else if (genebasedTest == "BT") {
      x <-
        sumFREGAT::BT(
          score.file = score.file,
          gene.file = paste0(DataDir,"/",gene_file),
          cor.path = paste0(ResultDir,"/cormatrix"),
          genes = genes1,
          anno.type = anno_type,
          beta.par = beta_par,
          weights.function = weights_function,
          user.weights = FALSE,
          write.file = FALSE,
          quiet = TRUE
        )
      return(x)

    } else if (genebasedTest == "PCA") {
      x <-
        sumFREGAT::PCA(
          score.file = score.file,
          gene.file = paste0(DataDir,"/",gene_file),
          genes = genes1,
          cor.path = paste0(ResultDir,"/cormatrix"),
          approximation = gene_approximation,
          anno.type = anno_type,
          n = sample_size,
          beta.par = beta_par,
          weights.function = weights_function,
          user.weights = FALSE,
          reference.matrix = reference_matrix_used,
          fun = regularize_fun,
          var.fraction = pca_var_fraction,
          write.file = FALSE,
          quiet = TRUE
        )
      return(x)

    } else if (genebasedTest == "FLM") {
      x <-
        sumFREGAT::FLM(
          score.file = score.file,
          gene.file = paste0(DataDir,"/",gene_file),
          genes = genes1,
          cor.path = paste0(ResultDir,"/cormatrix"),
          approximation = gene_approximation,
          anno.type = anno_type,
          n = sample_size,
          beta.par = beta_par,
          weights.function = weights_function,
          user.weights = FALSE,
          basis.function = flm_basis_function,
          k = flm_num_basis,
          order = flm_poly_order,
          flip.genotypes = flip_genotypes,
          Fan = omit_linear_variant,
          reference.matrix = reference_matrix_used,
          fun = regularize_fun,
          write.file = FALSE,
          quiet = TRUE
        )
      return(x)

    } else if (genebasedTest == "simpleM") {
      x <-
        sumFREGAT::simpleM(
          score.file = score_file,
          gene.file = gene_file,
          genes = genes1,
          #cor.path = cor_path,
          anno.type = anno_type,
          var.fraction = pca_var_fraction,
          write.file = FALSE,
          quiet = TRUE
        )
      return(x)

    } else if (genebasedTest == "minp") {
      x <-
        sumFREGAT::minp(
          score.file = score_file,
          gene.file = gene_file,
          genes = genes1,
          #cor.path = cor_path,
          anno.type = anno_type,
          write.file = FALSE,
          quiet = TRUE
        )

      return(x)

    }
  }


#' SexDiff: Sex difference in effect size for each SNP.
#'
#' @author Banabithi Bose
#'
#' @description This function uses the GWAS summary statistics from sex-stratified tests like FM01comb or FM02comb, to evaluate the difference in effect size between males and females at each SNP using a t-test.
#'
#' The input dataframes should only include X-chromosome in order to obtain results for sex differences based solely on X-linked loci.

#' @param Mfile R dataframe of summary statistics of GWAS or XWAS of male samples with six mandatory columns, SNP(Variant),CHR(Chromosome number)
#' ,BP(Base pair position),A1(Minor allele),BETA_M(Effect size) and SE_M(Standard error). This can be generated by running FM01comb or
#'  FM02comb model with GXWAS function.
#' @param Ffile R dataframe of summary statistics of GWAS or XWAS of male samples with six mandatory columns, SNP(Variant),CHR(Chromosome number)
#' ,BP(Base pair position),A1(Minor allele),BETA_F(Effect size) and SE_F(Standard error). This can be generated by running FM01comb or
#'  FM02comb model with GXWAS function.
#' @return R dataframe with seven columns, SNP(Variant),CHR(Chromosome number)
#' ,BP(Base pair position),A1(Minor allele), tstat(t-statistics for effect-size test), P(p-value) and adjP (Bonferroni corrected p-value).
#'
#' @importFrom stats cor pt
#'
#' @export
#'
#' @examples
#' #load("/projects/b1137/BBose/GXwasR/data/Mfile.Rda")
#' #load("/projects/b1137/BBose/GXwasR/data/Ffile.Rda")
#' data("GXwasRData")
#' Difftest <- SexDiff(Mfile,Ffile)
#' significant_snps <- Difftest[Difftest$adjP <0.05,]# 9

SexDiff <- function(Mfile,Ffile){

  #library(qqman)

  MFWAS <- merge(na.omit(Mfile),na.omit(Ffile),by = c("SNP","CHR","BP","A1"))
  gc(reset = TRUE)

  r <-
    stats::cor(log(abs(MFWAS$SE_M)), log(abs(MFWAS$SE_F)), method = "spearman", use = "pairwise.complete.obs")

  MFWAS$tstat <-
    (log(exp(MFWAS$BETA_M)) - log(exp(MFWAS$BETA_F))) / (sqrt(abs((MFWAS$SE_M) ^
                                                                    2 - (abs(MFWAS$SE_F)) ^ 2 - 2 * r * (abs(MFWAS$SE_M)) * (abs(MFWAS$SE_F))
    )))
  gc(reset = TRUE)
  MFWAS$P <- stats::pt(q=MFWAS$tstat, df= 2*nrow(MFWAS)-2, lower.tail=FALSE)
  MFWAS <- data.table::as.data.table(MFWAS)
  MFWAS$adjP <- p.adjust(MFWAS$P, method = "bonferroni", n = nrow(MFWAS))
  gc(reset = TRUE)
  #x <- MFWAS[,c(1:4,11:13)]
  x <- MFWAS[,c("SNP","CHR","BP","A1","tstat","P","adjP")]
  y <- x[order(x$adjP),]
  # qq plot
  chisq <- qchisq(1-y$P,1)
  lamdaGC <- median(chisq)/qchisq(0.5,1)
  gc(reset = TRUE)
  par(mar = c(1, 1, 1, 1))
  qqman::qq(y$P, main = paste0(("QQ-plots for test of sex-differentiated/n effect size with GIF = "), round(lamdaGC,3)))
  gc(reset = TRUE)
  return(y)
}


#' QCsnp: Quality control (QC) for SNPs.
#'
#' @author Banabithi Bose
#'
#' @description This function performs the quality control for the variants of plink binary files.
#'
#' @param DataDir A character string for the file path of the input PLINK binary files.
#' @param ResultDir A character string for the file path where all output files will be stored. The default is tempdir().
#' @param finput Character string, specifying the prefix of the input PLINK binary files with both male and female samples. This file needs to be in DataDir.
#' @param foutput Character string, specifying the prefix of the output PLINK binary files if the filtering option for the SNPs is chosen. The default is "FALSE".
#' @param casecontrol Boolean value, 'TRUE' or 'FALSE' indicating if the input plink files has cases-control status or not. The default is FALSE.
#' @param hweCase Numeric value between 0 to 1 or NULL for removing SNPs which fail Hardy-Weinberg equilibrium for cases. The default is NULL.
#' @param hweControl Numeric value between 0 to 1 or NULL for removing SNPs which fail Hardy-Weinberg equilibrium for controls. The default is NULL.
#' @param hwe Numeric value between 0 to 1 or NULL for removing SNPs which fail Hardy-Weinberg equilibrium for entire dataset. The default is NULL.
#' @param maf Numeric value between 0 to 1 for removing SNPs with minor allele frequency less than the specified threshold. The default is 0.05.
#' @param geno Numeric value between 0 to 1 for removing SNPs that have less than the specified call rate. The default is 0.05. Users can set this as NULL for not applying this filter.
#' @param monomorphicSNPs Boolean value, 'TRUE' or 'FALSE' for filtering out monomorphic SNP. The default is "TRUE".
#' @param caldiffmiss Boolean value, 'TRUE' or 'FALSE', specifying whether to compute differential missingness between cases and controls for each SNP (threshold is 0.05/length(unique(No. of. SNPs in the test))). The default is TRUE.
#' @param diffmissFilter Boolean value, 'TRUE' or 'FALSE', specifying whether to filter out the SNPs or only flagged them for differential missingness in cases vs contols. The deafailt is "TRUE".
#' @param dmissX Boolean value, 'TRUE' or 'FALSE' for computing differential missingness between cases and controls for X chromosome SNPs only. The default is "FALSE". The diffmissFilter will work for all these SNPs.
#' @param dmissAutoY Boolean value, 'TRUE' or 'FALSE' for computing differential missingness between cases and controls for SNPs on autosomes and Y chromosome only. The default is "FALSE". If dmissX and dmissAutoY are both FALSE, then this will be computed genome-wide. The diffmissFilter will work for all these SNPs.
#' @param ld_prunning Boolean value, 'TRUE' or 'FALSE' for applying linkage disequilibrium (LD)-based filtering.
#' @param highLD_regions NULL or Character string, specifying the .txt file name with genomic regions with high LD. Based on a genome build Hg19 and Hg38, two files such as "high-LD-regions-hg19-GRCh37.txt" and "" are provided in "extdata" folder.
#' @param window_size Integer value, specifying a window size in the variant counts for LD-based filtering. The default is 50.
#' @param step_size Integer value, specifying a variant count to shift the window at the end of each step for LD filtering. The default is 5.
#' @param r2_threshold Numeric value between 0 to 1 of pairwise r^2 threshold for LD-based filtering. The default is 0.02.
#'
#' @return A list of two objects, namely, MonomorSNPs and DiffMissSNPs containing monomorphic SNPs and SNPs with differential missingness in cases vs controls, respectively. Output plink binary files in the working directory.
#' @export
#'
#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' foutput <- "Test_output"
#' geno <- NULL
#' maf <- 0.05
#' casecontrol <- FALSE
#' hweCase <- NULL
#' hweControl <- NULL
#' hweCase <- NULL
#' monomorphicSNPs <- FALSE
#' caldiffmiss <- FALSE
#' ld_prunning <- FALSE
#' x <- QCsnp(DataDir = DataDir, ResultDir = ResultDir, finput = finput,foutput = foutput, geno = geno, maf = maf,hweCase = hweCase, hweControl = hweControl, ld_prunning = ld_prunning, casecontrol = casecontrol, monomorphicSNPs = monomorphicSNPs, caldiffmiss = caldiffmiss)

QCsnp <-
  function(DataDir,
           ResultDir = tempdir(),
           finput,
           foutput,
           casecontrol = TRUE,
           hweCase = NULL,
           hweControl = NULL,
           hwe = NULL,
           maf = 0.05,
           geno = 0.1,
           monomorphicSNPs = FALSE,
           caldiffmiss = FALSE,
           diffmissFilter = FALSE,
           dmissX = FALSE,
           dmissAutoY = FALSE,
           highLD_regions,
           ld_prunning = FALSE,
           window_size = 50,
           step_size = 5,
           r2_threshold = 0.02) {

    if (file.exists(paste0(DataDir, "/", finput, ".bed")) &&
        file.exists(paste0(DataDir, "/", finput, ".bim")) &&
        file.exists(paste0(DataDir, "/", finput, ".fam"))) {

      setupPlink(ResultDir)

    } else{
      writeLines(
        "There are no Plink files in DataDir.\nPlease specify correct directory path with input Plink files."
      )
    }

    fam <-
      as.data.frame(utils::read.table(file = paste0(DataDir, "/", finput, ".fam")))

    fam$V6 <- as.numeric(as.character(fam$V6))
    fam <- stats::na.omit(fam)
    fam1 <- fam[fam$V5 != 0,]
    fam2 <- fam1[fam1$V6 != 0,]
    fam4 <- fam2[fam2$V6 != -9,]

    if (is.null(maf)){
      MAF = NULL
    }else{
      MAF = "--maf"
    }

    if (is.null(geno)){
      GENO = NULL
    }else{
      GENO = "--geno"
    }


    if (is.null(hwe)){
      HWE = NULL
    }else{

      if (!is.null(hweControl)){
        print("Since hwe is not NULL, hweControl should be NULL. Setting hweControl = NULL implicitly.")
        hweControl = NULL
      }else if (!is.null(hweCase)){
        print("Since hwe is not NULL, hweCase should be NULL. Setting hweCase = NULL implicitly.")
        hweCase = NULL
      }

      HWE = "--hwe"
    }

    if (is.null(hweCase) && !is.null(hweControl)){
      print("hweControl cannot be NULL if hweCase is not NULL.")
      return()
    }else if (!is.null(hweCase) && is.null(hweControl)){
      print("hweCase cannot be NULL if hweControl is not NULL.")
      return()
    }

    if (is.null(hweCase) && is.null(hweControl)){
      HWECase = NULL
      HWECon = NULL

    }else{
      HWECase = "--hwe"
      HWECon = "--hwe"
      HWE = NULL
      hwe = NULL
    }



    ## This will be done for the entire file irrespective of case-control status
    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--bfile",
        paste0(DataDir,"/",finput),
        MAF,
        maf,
        GENO,
        geno,
        HWE,
        hwe,
        "--make-bed",
        "--out",
        paste0(ResultDir,"/","filtered_temp1"),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))

    ## ADD AN ERROR CATCH PART FOR filtered_temp1
    if (file.exists(paste0(ResultDir,"/","filtered_temp1.bed"))){
      print("Thresholds for maf, geno and hwe worked.")
    }else{

      ## Add read error line from file

      ##
      invisible(sys::exec_wait(
        paste0(ResultDir,"/","./plink"),
        args = c(
          "--bfile",
          paste0(DataDir,"/",finput),
          "--make-bed",
          "--out",
          paste0(ResultDir,"/","filtered_temp1"),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))
    }
    #######
    if (length(unique(fam4$V6)) >= 2 && casecontrol == TRUE) {

      invisible(sys::exec_wait(
        paste0(ResultDir,"/","./plink"),
        args = c(
          "--bfile",
          paste0(DataDir,"/",finput),
          "--filter-cases",
          HWECase,
          hweCase,
          "--make-bed",
          "--out",
          paste0(ResultDir,"/","filtered_temp_hwe_case_filtered"),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))

      invisible(sys::exec_wait(
        paste0(ResultDir,"/","./plink"),
        args = c(
          "--bfile",
          paste0(DataDir,"/",finput),
          "--filter-controls",
          HWECon,
          hweControl,
          "--make-bed",
          "--out",
          paste0(ResultDir,"/","filtered_temp_hwe_control_filtered"),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))



      invisible(sys::exec_wait(
        paste0(ResultDir,"/","./plink"),
        args = c(
          "--bfile",
          paste0(ResultDir,"/","filtered_temp_hwe_case_filtered"),
          "--bmerge",
          paste0(ResultDir,"/","filtered_temp_hwe_control_filtered"),
          "--freq",
          "--make-bed",
          "--out",
          paste0(ResultDir,"/","filtered_temp2"),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))

      invisible(sys::exec_wait(
        paste0(ResultDir,"/","./plink"),
        args = c(
          "--bfile",
          paste0(ResultDir,"/","filtered_temp1"),
          "--bmerge",
          paste0(ResultDir,"/","filtered_temp2"),
          "--freq",
          "--make-bed",
          "--out",
          paste0(ResultDir,"/","filtered_temp3"),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))

      invisible(sys::exec_wait(
        paste0(ResultDir,"/","./plink"),
        args = c(
          "--bfile",
          paste0(ResultDir,"/","filtered_temp3"),
          "--freq",
          "--make-bed",
          "--out",
          paste0(ResultDir,"/","filtered_temp4"),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))



    }else{

      if (length(unique(fam4$V6)) == 1){
        print("There is no case-control status in the plink files. Setting casecontrol = FALSE implicitly.")
        casecontrol = FALSE
      }else{
        casecontrol = FALSE}

      invisible(sys::exec_wait(
        paste0(ResultDir,"/","./plink"),
        args = c(
          "--bfile",
          paste0(ResultDir,"/","filtered_temp1"),
          "--freq",
          "--make-bed",
          "--out",
          paste0(ResultDir,"/","filtered_temp4"),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))

    }

    freq <-
      read.table(file = paste0(ResultDir,"/","filtered_temp4",".frq"), stringsAsFactors = FALSE, header = TRUE)
    mmSNPs <- freq[which(freq[, 5] == 0), 2,drop = FALSE]

    write.table(
      mmSNPs,
      file = paste0(ResultDir,"/","monomorphicSNPs"),
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      eol = "\r\n",
      sep = " "
    )

    if (monomorphicSNPs == TRUE && nrow(mmSNPs) != 0) {
      mmSNP1 <- mmSNPs$SNP
      exclude <- "--exclude"
      excludemono <- paste0(ResultDir,"/","monomorphicSNPs")

    }else if (monomorphicSNPs == FALSE | nrow(mmSNPs) == 0) {
      mmSNP1 <- NULL
      exclude <- NULL
      excludemono <- NULL
      print("There is no monomorphic SNPs.")
    }

    ## LD pruning
    if (ld_prunning == TRUE) {
      excluderange <- "--exclude"

      if (!is.null(highLD_regions)){
        write.table(highLD_regions, file = paste0(ResultDir,"/","highLD_regions_temp"), quote = FALSE, row.names = FALSE, col.names = FALSE)
        #highLD_regions <- paste0(DataDir,"/",highLD_regions)
        highLD_regions <- paste0(ResultDir,"/","highLD_regions_temp")

      }else{highLD_regions = NULL}

      indep <- "--indep-pairwise"
      window_size <- window_size
      step_size <- step_size
      r2_threshold <- r2_threshold

    }else if (ld_prunning == FALSE){
      excluderange <- NULL
      highLD_regions <- NULL
      indep <- NULL
      window_size <- NULL
      step_size <- NULL
      r2_threshold <- NULL
    }

    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--bfile",
        paste0(DataDir,"/",finput),
        MAF,
        maf,
        GENO,
        geno,
        HWE,
        hwe,
        "--make-bed",
        "--out",
        paste0(ResultDir,"/","filtered_temp1"),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))

    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--bfile",
        paste0(ResultDir,"/","filtered_temp4"),
        exclude,
        excludemono,
        excluderange,
        highLD_regions,
        indep,
        window_size,
        step_size,
        r2_threshold,
        "--make-bed",
        "--out",
        paste0(ResultDir,"/","filtered_temp5"),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))

    SNPmissCC = NULL

    if (casecontrol == TRUE){
      ## For case-control data:
      ## Write case/control missingness test to [ *.missing ]
      ## compute differential call rates
      ## between cases and controls for each SNP

      if (dmissX == TRUE & dmissAutoY == TRUE) {

        chrfilter <- NULL
        chrv <- NULL
      }else if (dmissX == TRUE & dmissAutoY == FALSE){

        chrfilter <- "--chr"
        chrv <- 23

      }else if(dmissX == FALSE & dmissAutoY == TRUE){

        chrfilter <- "--not-chr"
        chrv <- 23
      }else if (dmissX == FALSE & dmissAutoY == FALSE){

        chrfilter <- NULL
        chrv <- NULL
        print("There is no setting for filtering based on differential missingness.")
      }

      if (caldiffmiss == FALSE){

        print("Filtering for differential missingness between cases and controls is turned off by caldiffmiss = FALSE")
      }else{

        invisible(sys::exec_wait(
          paste0(ResultDir,"/","./plink"),
          args = c(
            "--bfile",
            paste0(ResultDir,"/","filtered_temp5"),
            chrfilter,
            chrv,
            "--test-missing",
            "--adjust",
            "--make-bed",
            "--out",
            paste0(ResultDir,"/","filtered_temp_casecontrol"),
            "--silent"
          ),
          std_out = FALSE,
          std_err = FALSE
        ))

      }

      if (file.exists(paste0(ResultDir,"/","filtered_temp_casecontrol.missing.adjusted"))){
        ccmissing <-
          read.table(file = paste0(ResultDir,"/","filtered_temp_casecontrol.missing.adjusted"),
                     header = TRUE,
                     sep = "")
      }else {
        ccmissing = NULL
      }

      if (!is.null(ccmissing)){

        diffmiss <- 0.05/length(unique(ccmissing$SNP)) ##DISCUSS THIS TO CONFIRM.


        SNPmissCC <-
          ccmissing[ccmissing$BONF < diffmiss, "SNP", drop = TRUE]

        write.table(
          SNPmissCC,
          file = paste0(ResultDir,"/","SNPdifCallrate"),
          quote = FALSE,
          row.names = FALSE,
          col.names = FALSE,
          eol = "\r\n",
          sep = " "
        )
      }else{
        SNPmissCC = NULL
      }
      if (length(SNPmissCC) == 0) {

        SNPmissCC = NULL
        print("No SNP with differential missingness between cases and controls.")

        invisible(sys::exec_wait(
          paste0(ResultDir,"/","./plink"),
          args = c(
            "--bfile",
            paste0(ResultDir,"/","filtered_temp5"),
            "--make-bed",
            "--out",
            paste0(ResultDir,"/",foutput),
            "--silent"
          ),
          std_out = FALSE,
          std_err = FALSE
        ))
      } else if (diffmissFilter == FALSE){

        print("No filter was applied differential missingness between cases and controls. You have set diffmissFilter = FALSE")
        print("Output plink files with all other filtering for the SNPs are saved in ResultDir.")

        invisible(sys::exec_wait(
          paste0(ResultDir,"/","./plink"),
          args = c(
            "--bfile",
            paste0(ResultDir,"/","filtered_temp5"),
            "--make-bed",
            "--out",
            paste0(ResultDir,"/",foutput),
            "--silent"
          ),
          std_out = FALSE,
          std_err = FALSE
        ))

      }else if(length(SNPmissCC) != 0 && diffmissFilter == TRUE){
        ## exclude SNPs
        invisible(sys::exec_wait(
          paste0(ResultDir,"/","./plink"),
          args = c(
            "--bfile",
            paste0(ResultDir,"/","filtered_temp5"),
            "--exclude",
            paste0(ResultDir,"/","SNPdifCallrate"),
            "--make-bed",
            "--out",
            paste0(ResultDir,"/",foutput),
            "--silent"
          ),
          std_out = FALSE,
          std_err = FALSE
        ))


      }

    }else{
      print("Since casecontrol = FALSE, no filter based on differential missigness will be applied.")

      invisible(sys::exec_wait(
        paste0(ResultDir,"/","./plink"),
        args = c(
          "--bfile",
          paste0(ResultDir,"/","filtered_temp5"),
          "--make-bed",
          "--out",
          paste0(ResultDir,"/",foutput),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))
    }


    print(paste0("Output plink files prefixed as ,",foutput,", with passed SNPs are saved in ResultDir."))


    ftemp <- list.files(paste0(ResultDir,"/"),pattern = "filtered_temp")
    invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))

    if (file.exists(paste0(ResultDir,"/SNPdifCallrate"))){
      invisible(do.call(file.remove,list(paste0(ResultDir,"/","SNPdifCallrate"))))
    }
    if (file.exists(paste0(ResultDir,"/monomorphicSNPs"))){
      invisible(do.call(file.remove,list(paste0(ResultDir,"/","monomorphicSNPs"))))
    }

    invisible(do.call(file.remove,list(paste0(ResultDir,"/","plink"))))

    resultbim <- read.table(paste0(ResultDir,"/",foutput,".bim"))

    inputbim <- read.table(paste0(DataDir,"/",finput,".bim"))

    print(paste0("Input file has ",length(unique(inputbim$V2))," SNPs."))
    print(paste0("Output file has ", length(unique(resultbim$V2))," SNPs after filtering."))


    return(list(MonomorSNPs = mmSNP1,DiffMissSNPs = SNPmissCC))
  } ## END

#' EstimateHerit: Computing SNP heritability i.e., the proportion of phenotypic variance explained by SNPs.
#'
#' @description This functions performs two types of heritability estimation, (i)GREML:Genomic relatedness matrix (GRM) restricted maximum likelihood-based method following GCTA (1) and (ii)LDSC: LD score regression-based method following (2,3). For the details, please follow the associated paper.
#'
#' Prior using this function, users are recommended to apply QCsnp and QCsample to ensure data quality control.
#'
#' @param DataDir A character string for the file path of the all the input files.
#' @param ResultDir A character string for the file path where all output files will be stored. The default is tempdir().
#' @param finput Character string, specifying the prefix of the input PLINK binary files for the genotype data. This file needs to be in DataDir. For LDSC model, if the original genotype data is not available, Hapmap 3 or 1000Genome data can be used.
#' @param summarystat A dataframe object with GWAS summary statistics. The mandatory column headers in this dataframe are 'chr'(Chromosome code), 'pos'(Basepair position), 'a1' (First allele code), ‘rsid’ (i.e., SNP idenitifier), ‘beta’ (i.e., effect-size or logarithm of odds ratio), ‘beta_se’ (i.e., standard error of beta), ‘P’ (i.e., p-values) and 'n_eff' (i.e., effective sample size).  For case-control study, effective sample size should be 4 / (1/<# of cases> + 1/<# of controls>). The default is NULL.
#' @param ncores Integer value, specifying the number of cores to be used for running LDSC model. The default is nb_cores().
#' @param model Character string, specifying the heritability estimation model. There are two options, “GREML” or “LDSC”. The default is “GREML”.
#' @param byCHR Boolean value, 'TRUE' or 'FALSE', specifying whether the analysis will be performed chromosome wise or not. The default is FALSE.
#' @param r2_LD Numeric value, specifying the LD threshold for clumping in LDSC model. The default is 0.
#' @param LDSC_blocks Integer value, specifying the block size for performing jackknife variance estimator in LDSC model following (3). The default is 200.
#' @param REMLalgo Integer value of 0, 1 or 2, specifying the algorithm to run REML iterations, 0 for average information (AI), 1 for Fisher-scoring and 2 for EM. The default option is 0, i.e. AI-REML (1).
#' @param nitr Integer value, specifying the number of iterations for performing the REML. The default is 100.
#' @param cat_covarfile A character string, specifying the name of the categorical covariate file which is a plain text file with no header line; columns are family ID, individual ID and discrete covariates. The default is NULL. This file needs to be in DataDir.
#' @param quant_covarfile A character string, specifying the name of the quantitative covariate file which is a plain text file with no header line; columns are family ID, individual ID and continuous covariates. The default is NULL. This file needs to be in DataDir.
#' @param prevalance Numeric value, specifying the disease prevalence. The default is 0.01.
#' @param partGRM Boolean value, 'TRUE' or 'FALSE', specifying whether the GRM will be partitioned into n parts (by row) in GREML model. The default is FALSE.
#' @param autosome Boolean value, 'TRUE' or 'FALSE', specifying whether estimate of heritability will be done for autosomes or not. The default is ‘TRUE’.
#' @param Xsome Boolean value, 'TRUE' or 'FALSE', specifying whether estimate of heritability will be done for X chromosome or not. The default is ‘TRUE’.
#' @param nGRM Integer value, specifying the number of the partision of the GRM in GREML model. The default is 3.
#' @param cripticut Numeric value, specifying the threshold to create a new GRM of "unrelated" individuals in GREML model. The default is arbitrary chosen as 0.025 following (1).
#' @param minMAF Positive numeric value (< maxMAF), specifying the minimum threshold for the MAF filter of the SNPs in the GREML model.
#' @param maxMAF Positive numeric value (minMAF,1), specifying the maximum threshold for the MAF filter of the SNPs in the GREML model.
#' @param hg Boolean value, specifying the genome built, “hg19” or “hg38” to use chromosome length from UCSC genome browser and getting genes and proteins according to this built. The default is “hg19”.
#' @param PlotIndepSNP Boolean value, 'TRUE' or 'FALSE', specifying whether to use independent SNPs i.e., chromosome-wise LD pruned SNPs in the plots or not. The default is TRUE.
#' @param IndepSNP_window_size Integer value, specifying a window size in variant count or kilobase for LD-based filtering. The default is 50.
#' @param IndepSNP_step_size Integer value, specifying a variant count to shift the window at the end of each step for LD filtering for pruned SNPs in the plots. The default is 5.
#' @param IndepSNP_r2_threshold Numeric value between 0 to 1 of pairwise r^2 threshold for LD-based filtering for pruned SNPs in the plots. The default is 0.02.
#' @param highLD_regions Character string, specifying the .txt file name with genomic regions with high LD for using in finding pruned SNPs in the plots. This file needs to be in DataDir.
#'
#' @returns A dataframe with eight columns for GREML and ten columns for LDSC model if byCHR is TRUE. The columns, such as, "chromosome"(i.e., chromosome code),"snp_proportion" (i.e.,chromosome-wise SNP propotion)", "no.of.genes" (i.e., number of genes per chromosome), "no.of.proteins" (i.e., number of genes per chromosome),"size_mb" (i.e., chromosome length), "Source" (i.e., source of heritability), "Variance" (i.e., estimated heritability), and "SE" (i.e., standard error of the estimated heritability) are common for both GREML and LDSC model. The column, "Intercept" (i.e., LDSC regression intercept) and "Int_SE" (i.e., standard error of the intercept) will be two extra columns for LDSC models. Source column will have rows, such as V(1) (i.e., name of genetic variance), V(e) (i.e., residual variance), V(p) (i.e., phenotypic variance), V(1)/Vp (i.e., ratio of genetic variance to phenotypic variance), and V(1)/Vp_L (i.e., ratio of genetic variance to phenotypic variance in liability scale for binary phenotypes). If byCHR is FALSE, then the first five columns will not be reported in the dataframe.
#'
#' @references
#'(1)	Yang et al. (2011) GCTA: a tool for Genome-wide Complex Trait Analysis. Am J Hum Genet. 88(1): 76-82.
#'(2)	Bulik-Sullivan B, et al. LD score regression distinguishes confounding from polygenicity in genome-wide association studies. Nat Genet. 2014;47:291–295.
#'(3)	Florian Privé, Julyan Arbel, Bjarni J Vilhjálmsson, LDpred2: better, faster, stronger, Bioinformatics, Volume 36, Issue 22-23, 1 December 2020, Pages 5424–5431, https://doi.org/10.1093/bioinformatics/btaa1029
#'
#' @importFrom bigsnpr snp_readBed snp_attach snp_match coef_to_liab
#' @importFrom data.table as.data.table rbindlist
#'
#'
#' @export
#'
#' @examples
#' #Not Run
# DataDir = system.file("extdata", package = "GXwasR")
# ResultDir = tempdir()
# finput <- "GXwasR_example"
# #load(paste0("/projects/b1137/BBose/GXwasR/data/Summary_Stat_Ex1.Rda"))
# data("GXwasRData")
# test.sumstats <- na.omit(Summary_Stat_Ex1[Summary_Stat_Ex1$TEST=="ADD",c(1:4,6:8)])
# colnames(test.sumstats) <- c("chr","rsid","pos","a1","n_eff","beta","beta_se")
# summarystat = test.sumstats
# ncores = 3
# model = "GREML"
# byCHR = FALSE
# r2_LD = 0
# LDSC_blocks = 20
# REMLalgo = 0
# nitr = 3
# cat_covarfile = NULL
# quant_covarfile = NULL
# prevalance = 0.01
# partGRM = FALSE
# autosome = TRUE
# Xsome = TRUE
# nGRM = 3
# cripticut = 0.025
# minMAF = NULL
# maxMAF = NULL
# hg = "hg19"
# PlotIndepSNP = TRUE
# IndepSNP_window_size = 50
# IndepSNP_step_size = 5
# IndepSNP_r2_threshold = 0.02
# #highLD_regions = "high-LD-regions-hg19-GRCh37.txt"
# H2 <- EstimateHerit(DataDir = DataDir, ResultDir = ResultDir, finput = finput, summarystat = NULL, ncores, model = "GREML", byCHR = "TRUE", r2_LD = 0, LDSC_blocks = 20,REMLalgo = 0, nitr = 100, cat_covarfile = NULL, quant_covarfile = NULL,prevalance = 0.01, partGRM = FALSE, autosome = TRUE, Xsome = TRUE, nGRM = 3,cripticut = 0.025, minMAF = NULL, maxMAF = NULL,hg = "hg19",PlotIndepSNP = TRUE, IndepSNP_window_size = 50,IndepSNP_step_size = 5,IndepSNP_r2_threshold = 0.02,highLD_regions = highLD_hg19)


EstimateHerit <- function(DataDir, ResultDir, finput, summarystat= NULL, ncores, model = c("LDSC","GREML"), byCHR = FALSE,
                          r2_LD = 0, LDSC_blocks = 20,
                          REMLalgo = c(0,1,2), nitr = 100, cat_covarfile = NULL, quant_covarfile = NULL,
                          prevalance = 0.01, partGRM = FALSE, autosome = TRUE, Xsome = TRUE, nGRM = 3,
                          cripticut = 0.025, minMAF = NULL, maxMAF = NULL, hg = c("hg19","hg38"),
                          PlotIndepSNP = c(TRUE,FALSE), IndepSNP_window_size = 50, IndepSNP_step_size = 5,
                          IndepSNP_r2_threshold = 0.02, highLD_regions = NULL){

  #library(bigsnpr)
  #library(data.table)
  #library(stringr)
  if (file.exists(paste0(DataDir, "/", finput, ".bed")) &&
      file.exists(paste0(DataDir, "/", finput, ".bim")) &&
      file.exists(paste0(DataDir, "/", finput, ".fam"))) {

    setupPlink(ResultDir)
  }else{
    writeLines(
      "There are no Plink files in DataDir.\nPlease specify correct directory path with input Plink files."
    )
    #return()
  }


  if (is.null(minMAF)|is.null(maxMAF)){
    #Compute min maf and max maf
    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--bfile",
        paste0(DataDir,"/",finput),
        "--freq",
        "--out",
        paste0(ResultDir,"/","MAF"),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))
    MAF <- na.omit(read.table(paste0(ResultDir,"/MAF.frq"),header = TRUE))
    minmaf <- min(MAF$MAF)
    maxmaf <- max(MAF$MAF)
    miMAF <- paste0("MAF = [",minmaf)
    maMAF <- paste0(maxmaf,"] ")
  }else{
    miMAF <- paste0("MAF  < ",minMAF)
    maMAF <- paste0(maxMAF," >")
  }

  if (model == "LDSC"){

    # if (file.exists(paste0(DataDir, "/", finput, ".bed")) &&
    #     file.exists(paste0(DataDir, "/", finput, ".bim")) &&
    #     file.exists(paste0(DataDir, "/", finput, ".fam"))) {
    #
    #   setupPlink(ResultDir)
    #
    # } else{
    #   writeLines(
    #     "There are no Plink files in DataDir.\nPlease specify correct directory path with input Plink files."
    #   )
    #   #return()
    # }

    gxwas_bedfile <- paste0(DataDir,"/",finput,".bed")
    # Reading the bedfile and storing the data in temporary directory
    rds <- snp_readBed(gxwas_bedfile, backingfile = tempfile())

    # Loading the data from backing files
    test.bigSNP <- snp_attach(rds)
    test.G <- test.bigSNP$genotypes

    y <- test.bigSNP$fam$affection - 1 # num[1:No. of. samples]
    # n_eff = 4 / (1 / sum(y == 0) + 1 / sum(y == 1))##This is for binary??
    # ##the effective sample sizes per variant for a GWAS using logistic regression, and simply the
    # ## sample size for continuous traits
    # test.sumstats$n_eff <- n_eff
    ## test statistics needs to have a1 and a0 both.
    test.map <- setNames(test.bigSNP$map[-c(3)], c("chr", "rsid", "pos", "a1", "a0"))
    summarystat1 <- merge(summarystat,test.map,by=c("chr", "rsid", "pos", "a1"))
    summarystat2 <- summarystat1[,c(1:4,8,5,6,7)]
    ## If no or few variants are actually flipped, you might want to disable the strand flipping
    ##option (strand_flip = FALSE) and maybe remove the few that were flipped (errors?)
    test.df_beta <- snp_match(summarystat2, test.map, join_by_pos = FALSE)  # use rsid instead of pos
    test.df_beta <- na.omit(test.df_beta)

    if (byCHR == FALSE){

      snpld <- ComputeLD(DataDir = DataDir,ResultDir = ResultDir,finput = finput, ByCHR = FALSE, CHRnum = NULL, r2_LD = r2_LD )
      result <- ComputeLDSC(snpld = snpld, test.df_beta = test.df_beta, ncores = ncores, LDSC_blocks = LDSC_blocks)
      herit_result <- data.table::as.data.table(t(as.data.frame(result)))
      herit_result$Source <- "V(G)/Vp"
      #Scaling coefficient to convert e.g. heritability to the liability scale.
      h2_liab <- herit_result$h2 * coef_to_liab(prevalance)# keeping the default value of K_gwas = 0.5, since n_eff is taken care of.
      herit_result2 <- herit_result
      herit_result2$h2 <- h2_liab
      herit_result2$Source <- "V(G)/Vp_L"
      herit_result3 <- rbind(herit_result,herit_result2)
      herit_result3 <- herit_result3[,c(5,3,4,1,2)]
      colnames(herit_result3) <- c("Source","Variance","SE","Intercept","Int_SE")
      return(herit_result3)

    }else{

      bimfile <- read.table(paste0(DataDir, "/", finput, ".bim"))

      chrnum <- 1:length(unique(bimfile$V1))
      chrnum <- 1:3

      chrwiseLD <- function(chrnum){
        chromosome <- unique(bimfile$V1)[chrnum]

        if (PlotIndepSNP == TRUE) {
          ## Do chromosome-wise LD prunning
          ChrwiseLDprun(DataDir=DataDir,ResultDir=ResultDir,chromosome=chromosome,
                        highLD_regions=highLD_regions, IndepSNP_window_size=IndepSNP_window_size,
                        IndepSNP_step_size=IndepSNP_step_size,IndepSNP_r2_threshold=IndepSNP_r2_threshold)

          ##
          bimfile1 <- read.table(paste0(ResultDir, "/LDfiltered.bim"))


        }else{
          bimfile1 <- bimfile[bimfile$V1==chromosome,]
        }

        snp_proportion <- nrow(bimfile1)/nrow(bimfile)

        ## Getting number of genes and proteins
        GP <- GeneProtein(hg = hg,chromosome = chromosome)

        print(paste0("Processing chromosome ",chromosome))
        snpld <- ComputeLD(DataDir = DataDir,ResultDir = ResultDir,finput = finput, ByCHR = byCHR, CHRnum = chromosome, r2_LD = r2_LD )
        result1 <- ComputeLDSC(snpld = snpld, test.df_beta = test.df_beta, ncores = ncores, LDSC_blocks = LDSC_blocks)
        result1 <- as.data.frame(result1)
        herit_result1 <- data.table::as.data.table(t(result1))
        herit_result <- data.table::as.data.table(cbind(chromosome,snp_proportion,GP,herit_result1))
        herit_result$Source <- "V(G)/Vp"
        #Scaling coefficient to convert e.g. heritability to the liability scale.
        h2_liab <- herit_result$h2 * coef_to_liab(prevalance)# keeping the default value of K_gwas = 0.5, since n_eff is taken care of.
        herit_result2 <- herit_result
        herit_result2$h2 <- h2_liab
        herit_result2$Source <- "V(G)/Vp_L"
        herit_result3 <- rbind(herit_result,herit_result2)
        return(herit_result3)
      }

      result <- data.table::rbindlist(lapply(chrnum,chrwiseLD))
      x <- na.omit(result)
      if (nrow(x)<3){
        print("Not enough data points for plots.")
      }else{
        result1 <- result[result$Source == "V(G)/Vp",]##Check if needed
        #load necessary libraries
        #library(ggplot2)
        #library(ggpubr)

        chr_mb <- ChrLength(hg=hg)
        result2 <- merge(result1,chr_mb, by.y = "chrom", by.x = "chromosome")
        colnames(result2) <- c("chromosome","snp_proportion","no.of.genes","no.of.proteins","Intercept","Int_SE","Variance","SE","Source","size_mb")
        #create plot with regression line, regression equation, Pearson correlation and p-value.
        PlotHeritability(Hdata = result2,miMAF = miMAF, maMAF = maMAF)##Check yval

        result3 <- result2
        colnames(result3) <- c("Chromosome","SNP_proportion","No.of.genes","No.of.proteins","Intercept","Int_SE","Variance","SE","Source","Size_mb")
        result3 <- result3[,c(1:4,10,9,7,8,5,6)]
        return(result3)
      }
    }

  }else if (model == "GREML"){

    if (file.exists(paste0(DataDir, "/", finput, ".bed")) &&
        file.exists(paste0(DataDir, "/", finput, ".bim")) &&
        file.exists(paste0(DataDir, "/", finput, ".fam"))) {

      setupGCTA(ResultDir)

    } else{
      writeLines(
        "There are no Plink files in DataDir.\nPlease specify correct directory path with input Plink files."
      )
      #return()
    }

    ## Making phenofile in ResultDir
    famfile <- read.table(paste0(DataDir, "/", finput, ".fam"))[,c(1,2,6)]
    write.table(famfile, file = paste0(ResultDir,"/phenofile.phen"),sep = " ",col.names = FALSE, row.names = FALSE, quote = FALSE)

    if (byCHR == FALSE){

      if (autosome == TRUE && Xsome == FALSE){

        ## Compute GRM
        ComputeGRMauto(DataDir = DataDir, ResultDir = ResultDir, finput = finput,
                       partGRM = partGRM, nGRM = nGRM, cripticut = cripticut,  minMAF = minMAF, maxMAF = maxMAF)

        ## Compute REML
        herit_result <- ComputeREMLone(DataDir = DataDir, ResultDir = ResultDir, REMLalgo = REMLalgo, nitr = nitr, phenofile = "phenofile.phen", cat_covarfile = cat_covarfile,
                                       quant_covarfile = quant_covarfile, prevalance = prevalance, chr = chromosome, grmfile = "test", ncores = ncores)

        return(herit_result)

      }else if (autosome == TRUE && Xsome == TRUE){
        ## Compute GRM Autosome
        ComputeGRMauto(DataDir = DataDir, ResultDir = ResultDir, finput = finput,
                       partGRM = partGRM, nGRM = nGRM, cripticut = cripticut,  minMAF = minMAF, maxMAF = maxMAF)
        ## Compute GRM X
        ComputeGRMX(DataDir = DataDir, ResultDir = ResultDir, finput = finput,
                    partGRM = partGRM, nGRM = nGRM, minMAF = minMAF, maxMAF = maxMAF)

        ## Compute REML
        herit_result <- ComputeREMLmulti(DataDir = DataDir, ResultDir = ResultDir, REMLalgo = REMLalgo, nitr = nitr, phenofile = "phenofile.phen", cat_covarfile = cat_covarfile,
                                         quant_covarfile = quant_covarfile, prevalance = prevalance, grmfile = "multi_GRMs.txt",ncores = ncores)

        return(herit_result)

      }else if (autosome == FALSE && Xsome == TRUE){
        ## Compute GRM X
        ComputeGRMX(DataDir = DataDir, ResultDir = ResultDir, finput = finput,
                    partGRM = partGRM, nGRM = nGRM, minMAF = minMAF, maxMAF = maxMAF)
        ## Compute REML X
        herit_result <- ComputeREMLone(DataDir = DataDir, ResultDir = ResultDir, REMLalgo = REMLalgo, nitr = nitr, phenofile = "phenofile.phen", cat_covarfile = cat_covarfile,
                                       quant_covarfile = quant_covarfile, prevalance = prevalance, chr = chromosome, grmfile = "xtest", ncores = ncores)

        return(herit_result)

      }else{
        print("autosome and Xsome cannot be set as FALSE together.")
      }

    }else{


      bimfile <- read.table(paste0(DataDir, "/", finput, ".bim"))

      chrnum <- 1:length(unique(bimfile$V1))
      chrnum <- 1:3

      chrwiseRELM <- function(chrnum){

        chromosome <- unique(bimfile$V1)[chrnum]

        if (PlotIndepSNP == TRUE) {
          ## Do chromosome-wise LD prunning
          ChrwiseLDprun(DataDir=DataDir,ResultDir=ResultDir,finput = finput, chromosome=chromosome,
                        highLD_regions=highLD_regions, IndepSNP_window_size=IndepSNP_window_size,
                        IndepSNP_step_size=IndepSNP_step_size,IndepSNP_r2_threshold=IndepSNP_r2_threshold)

          ##
          bimfile1 <- read.table(paste0(ResultDir, "/LDfiltered.bim"))


        }else{
          bimfile1 <- bimfile[bimfile$V1==chromosome,]
        }

        snp_proportion <- nrow(bimfile1)/nrow(bimfile)

        ## Getting number of genes and proteins
        GP <- GeneProtein(hg = hg,chromosome = chromosome)

        print(paste0("Processing chromosome ",chromosome))

        if (chromosome == 23){
          ## Compute GRM X

          ComputeGRMX(DataDir = DataDir, ResultDir = ResultDir, finput = finput,
                      partGRM = partGRM, nGRM = nGRM, minMAF = minMAF, maxMAF = maxMAF)
          ## Compute REML X
          herit_result <- ComputeREMLone(DataDir = DataDir, ResultDir = ResultDir, REMLalgo = REMLalgo, nitr = nitr, phenofile = "phenofile.phen", cat_covarfile = cat_covarfile,
                                         quant_covarfile = quant_covarfile, prevalance = prevalance, chr = chromosome, grmfile = "xtest", ncores = ncores)
          herit_result <- data.table::as.data.table(cbind(chromosome,snp_proportion,GP,herit_result))
          return(herit_result)

        }else{
          ## Compute GRM
          ComputeGRMauto(DataDir = DataDir, ResultDir = ResultDir, finput = finput,
                         partGRM = partGRM, nGRM = nGRM, cripticut = cripticut,
                         minMAF = minMAF, maxMAF = maxMAF, ByCHR = byCHR, CHRnum = chromosome)
          #grmfile <- paste0(chromosome,"test")
          #print(grmfile)
          ## Compute REML
          herit_result <- ComputeREMLone(DataDir = DataDir, ResultDir = ResultDir, REMLalgo = REMLalgo, nitr = nitr, phenofile = "phenofile.phen", cat_covarfile = cat_covarfile,
                                         quant_covarfile = quant_covarfile, prevalance = prevalance, chr = chromosome, grmfile = "test", ncores = ncores)


          herit_result <- data.table::as.data.table(cbind(chromosome,snp_proportion,GP,herit_result))

          return(herit_result)
        }


      }

      result <- data.table::rbindlist(lapply(chrnum,chrwiseRELM))

      x <- na.omit(result)
      if (nrow(x)<3){
        print("Not enough data points for plots.")
        return(result)
      }else{
        result1 <- result[result$Source == "V(G)/Vp",]
        #print(result1)
        #load necessary libraries
        #library(ggplot2)
        #library(ggpubr)

        chr_mb <- ChrLength(hg=hg)
        print("line 384")
        result2 <- merge(result1,chr_mb, by.y = "chrom", by.x = "chromosome")
        #print(result2)
        #create plot with regression line, regression equation, Pearson correlation and p-value.
        PlotHeritability(Hdata = result2, miMAF = miMAF, maMAF = maMAF)

        result3 <- merge(result,chr_mb, by.y = "chrom", by.x = "chromosome")
        result3 <- result3[,c(1:4,8,5:7)]
        colnames(result3) <- c("Chromosome","SNP_proportion","No.of.genes","No.of.proteins","Size_mb","Source","Variance","SE")

        return(result3)
      }
    }

  }
  ftemp <- list.files(paste0(ResultDir,"/"),pattern = "test_reml")
  invisible(file.remove(paste0(ResultDir,"/",ftemp)))
  ftemp <- list.files(paste0(ResultDir,"/"),pattern = "LDfiltered")
  invisible(file.remove(paste0(ResultDir,"/",ftemp)))
  ftemp <- list.files(paste0(ResultDir,"/"),pattern = "LDsnp")
  invisible(file.remove(paste0(ResultDir,"/",ftemp)))
  ftemp <- list.files(paste0(ResultDir,"/"),pattern = "grm")
  invisible(file.remove(paste0(ResultDir,"/",ftemp)))
  ftemp <- list.files(paste0(ResultDir,"/"),pattern = "test")
  invisible(file.remove(paste0(ResultDir,"/",ftemp)))
  ftemp <- list.files(paste0(ResultDir,"/"),pattern = "gcta")
  invisible(file.remove(paste0(ResultDir,"/",ftemp)))
  ftemp <- list.files(paste0(ResultDir,"/"),pattern = "MAF")
  invisible(file.remove(paste0(ResultDir,"/",ftemp)))
  ftemp <- list.files(paste0(ResultDir,"/"),pattern = "GRM")
  invisible(file.remove(paste0(ResultDir,"/",ftemp)))
  invisible(file.remove(paste0(ResultDir,"/","phenofile.phen")))
}


#' ComputeGeneticPC: Computing principal components from genetic relationship matrix
#'
#'@description  This function performs principal components analysis (PCA) based on the variance-standardized relationship matrix (1).
#' Top principal components are generally used as covariates in association analysis regressions to help correct for population stratification
#'
#' @param DataDir A character string for the file path of the input PLINK binary files.
#' @param finput Character string, specifying the prefix of the input PLINK binary files. This file needs to be in DataDir.
#' @param ResultDir A character string for the file path where all output files will be stored. The default is tempdir().
#' @param countPC Integer value, specifying the number of principal components. The default is 10.
#' @param plotPC Boolean value, 'TRUE' or 'FALSE', specifying whether to plot the first two PCs.
#' @param highLD_regions A R dataframe with genomic regions with high LD for using in finding pruned SNPs in the plots. The default is NULL.
#' @param ld_prunning Numeric value between 0 to 1 of pairwise r^2 threshold for LD-based filtering for pruned SNPs in the plots. The default is 0.02.
#' @param window_size Integer value, specifying a window size in variant count or kilobase for LD-based filtering. The default is 50.
#' @param step_size Integer value, specifying a variant count to shift the window at the end of each step for LD filtering for pruned SNPs in the plots. The default is 5.
#' @param r2_threshold Numeric value between 0 to 1 of pairwise r^2 threshold for LD-based filtering for pruned SNPs in the plots. The default is 0.02.
#'
#' @return A dataframe with genetic principal components. The first two columns are IID (i.e., Individual Id) and FID (i.e., Family ID). The other columns are PCs.
#'
#' @references (1) Purcell et. al.,(2007). PLINK: A Tool Set for Whole-Genome Association and Population-Based Linkage Analyses. The American Journal of Human Genetics, 81(3), 559-575. https://doi.org/10.1086/519795.
#'
#' @importFrom ggplot2 ggplot geom_bar ylab xlab theme_classic geom_point theme_light coord_equal aes_string
#' @importFrom ggpubr ggarrange
#'
#' @export
#'
#' @examples
#' #' # Not Run
# DataDir <- system.file("extdata", package = "GXwasR")
# ResultDir <- tempdir()
# finput <- "GXwasR_example"
# data("GXwasRData")
# #highLD_regions <- "high-LD-regions-hg19-GRCh37.txt"
# highLD_regions <- highLD_hg19
# ld_prunning <- "TRUE"
# window_size <- 50
# step_size <- 5
# r2_threshold <- 0.02
# countPC <- 20
# ## Genetic PC
# GP <- ComputeGeneticPC(DataDir = DataDir,ResultDir=ResultDir,finput=finput,highLD_regions = highLD_hg19, countPC = 20)

ComputeGeneticPC <- function(DataDir,ResultDir = tempdir(),finput,countPC = 10,plotPC = TRUE,
                             highLD_regions = NULL, ld_prunning = TRUE,
                             window_size = 50, step_size = 5,r2_threshold = 0.02){

  if (file.exists(paste0(DataDir, "/", finput, ".bed")) &&
      file.exists(paste0(DataDir, "/", finput, ".bim")) &&
      file.exists(paste0(DataDir, "/", finput, ".fam"))) {

    setupPlink(ResultDir)

  } else{
    writeLines(
      "There are no Plink files in DataDir.\nPlease specify correct DataDir path with input Plink files."
    )
  }
  ## LD pruning
  if (ld_prunning == TRUE && !is.null(highLD_regions)) {
    excluderange <- "--exclude"
    write.table(highLD_regions, file = paste0(ResultDir,"/","highLD_regions_temp"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    #highLD_regions <- paste0(DataDir,"/",highLD_regions)
    highLD_regions <- paste0(ResultDir,"/","highLD_regions_temp")
    indep <- "--indep-pairwise"
    window_size <- window_size
    step_size <- step_size
    r2_threshold <- r2_threshold

    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c("--bfile",
               paste0(DataDir,"/",finput),
               excluderange,
               highLD_regions,
               indep,
               window_size,
               step_size,
               r2_threshold,
               "--make-bed",
               "--out",
               paste0(ResultDir,"/","pruned_",finput),
               "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))
    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--bfile",
        paste0(ResultDir,"/","pruned_",finput),
        "--pca",countPC,
        "--out",
        paste0(ResultDir, "/","pcfile"),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))
  }else{

    # PCA for ancestry: family and individual ID in columns 1 and 2, followed by the principal components.
    invisible(sys::exec_wait(
      paste0(ResultDir, "/","./plink"),
      args = c(
        "--bfile",
        paste0(DataDir, "/",finput),
        "--pca",countPC,
        "--out",
        paste0(ResultDir, "/","pcfile"),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))
  }
  PCs1 <- read.table(paste0(ResultDir,"/pcfile.eigenvec"))
  PCs <- PCs1[,-c(1:2)]
  names(PCs)[1:ncol(PCs)] <- paste0("PC",(1:ncol(PCs)))
  EV <- scan(paste0(ResultDir,"/pcfile.eigenval"))
  Percent.var <- data.frame(PC = 1:ncol(PCs),Percent.var=EV/sum(EV)*100)
  #library(ggplot2)
  if (plotPC == TRUE){
    ##Plots
    p1<- ggplot2::ggplot(Percent.var,ggplot2::aes_string(Percent.var$PC,Percent.var$Percent.var)) + ggplot2::geom_bar(stat = "identity")+
      ggplot2::ylab("Percent Variance Explained")+ ggplot2::theme_classic()

    p2 <- ggplot2::ggplot(data=PCs, ggplot2::aes_string(PCs$PC1,PCs$PC2)) + ggplot2::geom_point() + ggplot2::xlab(paste0("PC1 (",signif(Percent.var$Percent.var[1]),"% )"))+
      ggplot2::ylab(paste0("PC2 (",signif(Percent.var$Percent.var[2]),"% )")) + ggplot2::theme_light() + ggplot2::coord_equal()
    #library(ggpubr)
    print(ggpubr::ggarrange(p1, p2,
                            labels = c("A", "B"),
                            ncol = 2, nrow = 1))

  }

  ftemp <- list.files(paste0(ResultDir,"/"),pattern = "pcfile")
  invisible(file.remove(paste0(ResultDir,"/",ftemp)))
  ftemp <- list.files(paste0(ResultDir,"/"),pattern = "temp")
  invisible(file.remove(paste0(ResultDir,"/",ftemp)))
  return(PCs1)
}

#' ComputePRS: Computing polygenic risk score (PRS)
#'
#' @description This function calculates the polygenic risk score, which is the total of allele counts (genotypes) weighted by estimated effect sizes from genome-wide association studies. It uses C+T filtering techniques. The users could perform clumping procedure choromosome-wise and genome-wide. Also, the function offers the choice of including several genetic principal components along with other covariates. Using this function, users have the freedom to experiment with various clumping and thresholding arrangements to test a wide range of various parameter values.
#'
#' @param DataDir A character string for the file path of the all the input files.
#' @param ResultDir A character string for the file path where all output files will be stored. The default is tempdir().
#' @param finput Character string, specifying the prefix of the input PLINK binary files for the genotype data i.e., the target data based on which clumping procedure will be performed. This file needs to be in DataDir. If your target data are small (e.g. N < 500) then you can use the 1000 Genomes Project samples. Make sure to use the population that most closely reflects represents the base sample.
#' @param summarystat A dataframe object with GWAS summary statistics. The mandatory column headers in this dataframe are 'CHR'(Chromosome code), 'BP'(Basepair position), 'A1' (effect allele), ‘SNP’ (i.e., SNP idenitifier), ‘BETA’ (i.e., effect-size or logarithm of odds ratio), and ‘P’ (i.e., p-values).
#' @param phenofile A character string, specifying the name of the mandatory phenotype file. This is a plain text file with no header line; columns family ID, individual ID and phenotype columns. For binary trait, the phenotypic value should be coded as 0 or 1, then it will be recognized as a case-control study (0 for controls and 1 for cases). Missing value should be represented by "-9" or "NA". The interested phenotype column should be labeled as "Pheno1". This file needs to be in DataDir.
#' @param covarfile A character string, specifying the name of the covariate file which is a plain .text file with no header line; columns are family ID, individual ID and the covariates. The default is NULL. This file needs to be in DataDir.
#' @param pheno_type Boolean value, ‘binary’ or ‘quantitative’, specifying the type of the trait. The default is ‘binary’.
#' @param effectsize Boolean value, ‘BETA’ or ‘OR’, specifying the type of the GWAS effectsize. The default is ‘BETA’.
#' @param ldclump Boolean value, 'TRUE' or 'FALSE', specifying whether to perform clumping or not.
#' @param LDreference A character string, specifying the  prefix of the PLINK files of the population reference panel of the same ancestry, and ideally the one that was used for imputing your target dataset. These files should be in DataDir.
#' @param clump_p1 Numeric value, specifying the significance threshold for index SNPs if 'ldclump' was set to be TRUE. The default is 0.0001.
#' @param clump_p2 Numeric value, specifying the secondary significance threshold for clumped SNPs if 'ldclump' was set to be TRUE. The default is 0.01
#' @param clump_r2 Numeric value, specifying the linkage disequilibrium (LD) threshold for clumping if 'ldclump' was set to be TRUE. The default is 0.50.
#' @param clump_kb Integer value, specifying the physical distance threshold in base-pair for clumping if 'ldclump' was set to be TRUE. The default is 250.
#' @param byCHR Boolean value, 'TRUE' or 'FALSE', specifying chromosome-wise clumping procedure if 'ldclump' was set to be TRUE. The default is TRUE.
#' @param pthreshold Numeric vector, containing several p value thresholds to maximize predictive ability of the derived polygenic scores.
#' @param ld_prunning Boolean value, 'TRUE' or 'FALSE' for LD-based filtering for computing genetic PC as covariates.
#' @param nPC Positive integer value, specifying the number of genetic PCs to be included as predictor in the PRS model fit. The default is 6.
#' @param window_size Integer value, specifying a window size in variant count or kilobase for LD-based filtering in computing genetic PC. The default is 50.
#' @param step_size Integer value, specifying a variant count to shift the window at the end of each step for LD filtering in computing genetic PCs. The default is 5.
#' @param r2_threshold Numeric value between 0 to 1 of pairwise r^2 threshold for LD-based filtering in computing genetic PCs. The default is 0.02.
#' @param highLD_regions Character string, specifying the .txt file name with known genomic regions with high LD. The default is NULL.
#'
#' @return A list object containing a dataframe and a numeric value. The dataframe,PRS, contains four mandatory columns, such as,
#' IID (i.e., Individual ID), FID (i.e., Family ID), Pheno1 (i.e., the trait for PRS) and Score (i.e., the best PRS). Other columns of covariates could be there. The numeric value, BestP contains the threshold of
#' of the best p-value for the best pRS model fit.
#'
#' Also, the function produces several plots such as p-value thresholds vs PRS model fit and PRS distribution among male and females. For case-control data, it shows PRS distribution among cases and controls and ROC curves as well.
#'
#' @importFrom dplyr distinct
#' @importFrom stats lm predict logLik
#' @importFrom ggplot2 theme_classic theme element_text ggtitle geom_density xlab aes scale_y_continuous geom_bar scale_fill_gradient2
#' @importFrom data.table as.data.table
#' @importFrom ggpubr ggarrange
#' @export
#'
#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' data("GXwasRData")
#' #load("/projects/b1137/BBose/GXwasR/data/Summary_Stat_Ex1.Rda")
#' summarystat <- Summary_Stat_Ex1[,c(2,4,7,1,3,12)]
#' #load("/projects/b1137/BBose/GXwasR/data/Example_phenofile.Rda")
#' phenofile <- Example_phenofile #Cannot be NULL, the interested phenotype column should be labeled as
#' #load("/projects/b1137/BBose/GXwasR/data/Example_covarfile.Rda") ## Add "Sex" as covar
#' covarfile <- Example_covarfile
#' clump_p1 = 0.0001
#' clump_p2 = 0.0001
#' clump_kb = 500
#' clump_r2 = 0.5
#' byCHR = TRUE
#' #load("/projects/b1137/BBose/GXwasR/data/Example_pthresoldfile.Rda")
#' pthreshold <- Example_pthresoldfile$Threshold
#' ld_prunning <- TRUE
#' highLD_regions <- "high-LD-regions-hg19-GRCh37.txt"
#' window_size <- 50
#' step_size <- 5
#' r2_threshold <- 0.02
#' nPC = 6 #We can incorporate PCs into our PRS analysis to account for population stratification.
#' pheno_type = "binary"
#'
#' PRSresult <- ComputePRS(DataDir,ResultDir,finput,summarystat,phenofile,covarfile,
#'                                     effectsize="BETA",LDreference = "GXwasR_example", ldclump = FALSE,clump_p1,clump_p2,clump_r2,clump_kb, byCHR = TRUE,
#'                                     pthreshold = pthreshold,highLD_regions = "high-LD-regions-hg19-GRCh37.txt", ld_prunning = TRUE,
#'                                     window_size = 50, step_size = 5,r2_threshold = 0.02, nPC = 6,pheno_type = "binary")
#' ## PRS
#' PRS <- PRSresult$PRS
#' ## The best threshold
#' BestPvalue <- PRSresult$BestP$Threshold

ComputePRS <- function(DataDir,ResultDir = tempdir(),finput,summarystat,phenofile,covarfile = NULL,
                       effectsize = c("BETA","OR"),ldclump = FALSE,LDreference,clump_p1,clump_p2,clump_r2,clump_kb, byCHR = TRUE,
                       pthreshold = c(0.001,0.05,0.1,0.2,0.3,0.4,0.5),highLD_regions = "high-LD-regions-hg19-GRCh37.txt", ld_prunning = FALSE,
                       window_size = 50, step_size = 5,r2_threshold = 0.02, nPC = 6,pheno_type = "binary"){

  #library(data.table)
  if (file.exists(paste0(DataDir, "/", finput, ".bed")) &&
      file.exists(paste0(DataDir, "/", finput, ".bim")) &&
      file.exists(paste0(DataDir, "/", finput, ".fam"))) {

    setupPlink(ResultDir)
  }else{
    writeLines(
      "There are no Plink files in DataDir.\nPlease specify correct directory path with input Plink files."
    )
  }

  if (effectsize == "OR"){
    summarystat$OR <- log(summarystat$OR)
  }
  #library(dplyr)

  summarystat <- as.data.frame(dplyr::distinct(summarystat,SNP,.keep_all=TRUE))
  write.table(summarystat, file=paste0(ResultDir,"/","prssummarystat"), quote=FALSE, row.names=FALSE)
  SNP.pvalue <- unique(summarystat[,c("SNP","P")])
  write.table(SNP.pvalue, file=paste0(ResultDir,"/","SNP.pvalue"), quote=FALSE, row.names=FALSE)


  ## Performing clumping
  #clumpedResult <- ClumpLD(DataDir = DataDir,ResultDir = ResultDir, finput = finput,SNPdata = "prssummarystat",
  #                        clump_p1 = clump_p1,clump_p2 = clump_p2,clump_r2 = clump_r2,clump_kb = clump_kb, byCHR = byCHR)

  if (ldclump == TRUE){
    clumpedResult <- ClumpLD(DataDir = DataDir,ResultDir = ResultDir, finput = LDreference,SNPdata = list(summarystat),
                             clump_p1 = clump_p1,clump_p2 = clump_p2,clump_r2 = clump_r2,clump_kb = clump_kb, byCHR = byCHR)

    write.table(clumpedResult$SNP, file = paste0(ResultDir,"/","Valid.SNP"), quote =  FALSE, row.names = FALSE)

    clumpExtract <- "--extract"
    clumpSNP <- paste0(ResultDir,"/","Valid.SNP")

  }else{
    clumpExtract <- NULL
    clumpSNP <- NULL

  }
  # Read in the phenotype file
  #phenotype <- read.table(paste0(DataDir,"/",phenofile), header= TRUE)
  #phenotype <- phenotype[,1:3]
  phenotype <- cbind(phenofile[,1:2],phenofile[,"Pheno1"])
  colnames(phenotype) <- c("FID","IID","Pheno1")

  # Genetic PC
  if (nPC > 0){
    # Compute PC
    GP <- ComputeGeneticPC(DataDir = DataDir,ResultDir = ResultDir, finput = finput,countPC = nPC,highLD_regions = highLD_regions,
                           ld_prunning = ld_prunning,window_size = window_size, step_size = step_size,r2_threshold = r2_threshold, plotPC = FALSE)
    # PCs
    colnames(GP) <- c("FID", "IID", paste0("PC",1:6))

  }else{
    print("Parameter 'nPC' is either zero or negative. Genetic PC will not be computed.")
    GP <- phenotype[,1:2]
  }

  # Read the covariates (here, it is sex)
  if (is.null(covarfile)){
    pheno <- merge(phenotype, GP, by=c("FID", "IID"))
  }else{
    #covariate <- read.table(paste0(DataDir,"/",covarfile), header= TRUE)
    covariate <- covarfile
    # Now merge the files
    pheno <- merge(merge(phenotype, covariate, by=c("FID", "IID")), GP, by=c("FID","IID"))
  }

  # We can then calculate the null model (model with PRS) using a linear regression
  # (as height is quantitative) ## Check for binary

  if (pheno_type == "binary"){
    phenoBin <- pheno
    phenoBin$Pheno1[phenoBin$Pheno1 == 2] <- 0
    phenoBin$Pheno1 <- as.integer(as.character(phenoBin$Pheno1))
    null.model2 <- glm(phenoBin$Pheno1~., data=phenoBin[,!colnames(phenoBin)%in%c("FID","IID")],family=binomial)
    # And the R2 of the null model is
    null.r2 <- summary(null.model2)$r.squared

  }else{
    null.model <- stats::lm(pheno$Pheno1~., data=pheno[,!colnames(pheno)%in%c("FID","IID")])
    # And the R2 of the null model is
    null.r2 <- summary(null.model)$r.squared
  }
  ## PRS using thresholding
  ## Best fit PRS
  prsFun <- function(pthreshold){
    print(pthreshold)
    print(paste0("Computing PRS for threshold ",pthreshold))
    #library(data.table)
    pt <- data.table::as.data.table(cbind(pthreshold,0,pthreshold))
    colnames(pt) <- c("Threshold","Lowerbound","UpperBound")
    write.table(pt, file = paste0(ResultDir,"/range_list"),quote = FALSE, row.names = FALSE)
    #By default, if a genotype in the score is missing for a particular individual, then the expected value is imputed, i.e. based on the sample allele frequency. To change this behavior, add the flag
    #--score-no-mean-imputation
    invisible(sys::exec_wait(
      paste0(ResultDir, "/","./plink"),
      args = c(
        "--bfile",paste0(DataDir, "/",finput),
        "--score",paste0(ResultDir,"/","prssummarystat"),1,2,3,"header",
        "--q-score-range",paste0(ResultDir,"/range_list"),paste0(ResultDir,"/","SNP.pvalue"),
        clumpExtract,clumpSNP,
        "--out",
        paste0(ResultDir, "/","PRS"),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))

    prs <- read.table(paste0(ResultDir,"/","PRS.",pthreshold,".profile"), header=TRUE)
    # Merge the prs with the phenotype matrix
    # We only want the FID, IID and PRS from the PRS file, therefore we only select the
    # relevant columns
    pheno.prs <- merge(pheno, prs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))
    # Now perform a logistic regression on trait with PRS and the covariates
    if (pheno_type == "binary"){
      phenoBin.prs <- pheno.prs
      phenoBin.prs$Pheno1[phenoBin.prs$Pheno1 == 2] <- 0
      phenoBin.prs$Pheno1 <- as.integer(as.character(phenoBin.prs$Pheno1))
      model2<-glm(phenoBin.prs$Pheno1~., data=phenoBin.prs[,!colnames(phenoBin.prs)%in%c("FID","IID")],
                  family=binomial)
      S=stats::predict(model2,type="response")
      McFaddenR2 <- 1-stats::logLik(model2)/stats::logLik(null.model2)
      prs.r2 <- McFaddenR2[1]
      prs.coef <- summary(model2)$coeff["SCORE",]
      prs.beta <- as.numeric(prs.coef[1])
      prs.se <- as.numeric(prs.coef[2])
      prs.p <- as.numeric(prs.coef[4])

    }else{
      # ignoring the FID and IID from our model
      model <- stats::lm(pheno$Pheno1~., data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")])
      # model R2 is obtained as
      model.r2 <- summary(model)$r.squared
      # R2 of PRS is simply calculated as the model R2 minus the null R2
      prs.r2 <- model.r2-null.r2
      # We can also obtain the coeffcient and p-value of association of PRS as follow
      prs.coef <- summary(model)$coeff["SCORE",]
      prs.beta <- as.numeric(prs.coef[1])
      prs.se <- as.numeric(prs.coef[2])
      prs.p <- as.numeric(prs.coef[4])
    }
    # We can then store the results
    prs.result <- rbind(data.frame(pthreshold, R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))
    colnames(prs.result) <- c("Threshold","R2","P","BETA","SE")
    return(prs.result)
  }
  prsResult <- data.table::rbindlist(lapply(pthreshold,prsFun))

  # Plots
  # library(ggplot2)
  # generate a pretty format for p-value output
  prsResult$WriteP <- round(prsResult$P, digits = 3)
  prsResult$WriteP[!is.na(prsResult$WriteP) & prsResult$WriteP == 0] <- format(prsResult$P[!is.na(prsResult$WriteP) & prsResult$WriteP == 0], digits = 2)
  prsResult$WriteP <- sub("e", "*x*10^", prsResult$WriteP)
  # Initialize ggplot, requiring the threshold as the x-axis (use factor so that it is uniformly distributed)
  p1 <- ggplot2::ggplot(data = prsResult, ggplot2::aes(x = factor(prsResult$Threshold), y = prsResult$R2)) +
    # Specify that we want to print p-value on top of the bars
    ggplot2::geom_text(
      ggplot2::aes(label = paste(prsResult$WriteP)),
      vjust = -1.5,
      hjust = 0,
      angle = 45,
      cex = 2,
      parse = T
    )  +
    # Specify the range of the plot, *1.25 to provide enough space for the p-values
    ggplot2::scale_y_continuous(limits = c(0, max(prsResult$R2) * 1.25)) +
    # Specify the axis labels
    ggplot2::xlab(expression(italic(prsResult$P) - value ~ threshold ~ (italic(prsResult$P)[T]))) +
    ggplot2::ylab(expression(paste("PRS model fit:  ", R ^ 2))) +
    # Draw a bar plot
    ggplot2::geom_bar(ggplot2::aes(fill = -log10(prsResult$P)), stat = "identity") +
    # Specify the colors
    ggplot2::scale_fill_gradient2(
      low = "dodgerblue",
      high = "firebrick",
      mid = "dodgerblue",
      midpoint = 1e-4,
      name = bquote(atop(-log[10] ~ model, italic(prsResult$P) - value),)
    ) +
    # Some beautification of the plot
    ggplot2::theme_classic() + ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold", size = 7),
      axis.text = ggplot2::element_text(size = 7),
      legend.title = ggplot2::element_text(face = "bold", size =
                                             7),
      legend.text = ggplot2::element_text(size = 7),
      axis.text.x = ggplot2::element_text(angle = 45, hjust =
                                            1)
    )+ ggplot2::ggtitle("P value thresolds vs PRS model fit") +
    ggplot2::theme(plot.title = element_text(size = 10,face = "bold"))

  # Best result is:
  bestP <- prsResult[which.max(prsResult$R2),"Threshold"]

  # Getting PRS score with best p-value threshold
  pt <- cbind(bestP,0,bestP)
  write.table(pt, file = paste0(ResultDir,"/range_list"),quote = FALSE, row.names = FALSE)
  #By default, if a genotype in the score is missing for a particular individual, then the expected value is imputed, i.e. based on the sample allele frequency. To change this behavior, add the flag
  #--score-no-mean-imputation
  invisible(sys::exec_wait(
    paste0(ResultDir, "/","./plink"),
    args = c(
      "--bfile",paste0(DataDir, "/",finput),
      "--score",paste0(ResultDir,"/","prssummarystat"),1,2,3,"header",
      "--q-score-range",paste0(ResultDir,"/range_list"),paste0(ResultDir,"/","SNP.pvalue"),
      clumpExtract,clumpSNP,
      "--out",
      paste0(ResultDir, "/","PRS"),
      "--silent"
    ),
    std_out = FALSE,
    std_err = FALSE
  ))

  prs <- read.table(paste0(ResultDir,"/","PRS.",bestP,".profile"), header=TRUE)
  pheno.prs <- merge(pheno, prs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))

  ## PRS with sex
  #library(ggplot2)
  # Read in the files
  #prs <- read.table("EUR.0.3.profile", header=T)
  d1 <- pheno.prs[,c("FID","IID","Pheno1"),drop = FALSE]
  famfile <- read.table(paste0(DataDir,"/",finput,".fam"), header=FALSE)
  sex <- famfile[!famfile$V5==0,c(1,2,5)]
  colnames(sex) <- c("FID","IID", "SEX")
  dat <- merge(d1,sex, by = c("FID","IID"))
  # Rename the sex
  dat$SEX[dat$SEX == 1] <- "Male"
  dat$SEX[dat$SEX == 2] <- "Female"
  dat$SEX <- as.factor(as.character(dat$SEX))
  #dat$Pheno1 <- as.factor(as.character(dat$Pheno1))
  #levels(hs$SEX) <- c("Male", "Female")
  # Merge the files
  dat <- merge(dat,prs,by=c("FID","IID"))
  # Start plotting
  # p2 <- ggplot(dat, aes(x=SCORE, y=Pheno1, color=SEX))+
  #   geom_point()+
  #   theme_classic()+
  #   labs(x="Polygenic Score", y="Pheno1")

  # Basic density plot with custom color
  p2 <- ggplot2::ggplot(dat, ggplot2::aes(x=dat$SCORE, color=dat$SEX)) +
    # color property for changing color of plot
    # geom_density() function plots the density plot
    ggplot2::geom_density()+ ggplot2::ggtitle("Best PRS distribution\n(males vs females)") +
    ggplot2::xlab("PRS")+
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10,face = "bold"))+
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = 8))

  if (pheno_type == "binary"){
    dat$Pheno1[dat$Pheno1 == 1] <- "control"
    dat$Pheno1[dat$Pheno1 == 2] <- "cases"

    dat$Pheno1 <- as.factor(as.character(dat$Pheno1))
    p3 <- ggplot2::ggplot(dat, ggplot2::aes(x=SCORE, color=Pheno1)) +
      # color property for changing color of plot
      # geom_density() function plots the density plot
      ggplot2::geom_density()+ggplot2::ggtitle("Best PRS distribution\n(cases vs controls)") +
      ggplot2::xlab("PRS")+
      ggplot2::theme(plot.title = ggplot2::element_text(size = 10,face = "bold"))+
      ggplot2::theme(axis.title.x = ggplot2::element_text(size = 8))

    dat1 <- dat[,c("Pheno1","SCORE")]

    mdat <- dat[dat$SEX == "Male",]
    p4 <- ggplot2::ggplot(mdat, ggplot2::aes(x=mdat$SCORE, color=mdat$Pheno1)) +
      # color property for changing color of plot
      # geom_density() function plots the density plot
      ggplot2::geom_density()+ggplot2::ggtitle("Best PRS distribution in males\n(cases vs controls)") +
      ggplot2::xlab("PRS")+
      ggplot2::theme(plot.title = ggplot2::element_text(size = 10,face = "bold"))+
      ggplot2::theme(axis.title.x = ggplot2::element_text(size = 8))

    fdat <- dat[dat$SEX == "Female",]
    p5 <- ggplot2::ggplot(fdat, ggplot2::aes(x=fdat$SCORE, color=fdat$Pheno1)) +
      # color property for changing color of plot
      # geom_density() function plots the density plot
      ggplot2::geom_density()+ ggplot2::ggtitle("Best PRS distribution in females\n(cases vs controls)") +
      ggplot2::xlab("PRS")+
      ggplot2::theme(plot.title = ggplot2::element_text(size = 10,face = "bold"))+
      ggplot2::theme(axis.title.x = ggplot2::element_text(size = 8))

    #library(pROC)
    #adjust plot margins
    par(mar = c(2, 4, 4, 2), cex = 0.7)
    range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    #library(grid)
    #library(ggplotify)
    #print("line 337")
    #dat1
    #x <- pROC::roc(as.vector((dat1$Pheno1)), range01(dat1$SCORE), direction=">",levels = c("cases","control"))
    #p6 <- as.grob(expression(plot.roc(x,col="red", tck = -.005,lwd=1,mgp=c(1.5, 0.5, 0), mar=c(6, 6, 1, 6)+.1,asp =0, cex.main= .8,cex.axis= .7,cex.lab=.8, main="ROC curve with best PRS",print.auc = TRUE, print.auc.cex = 0.5)))

    #ggpubr::ggarrange(p1,p2,p3,p4,p5,p6, ncol = 3, nrow = 2, labels=LETTERS[1:6],widths = c(1.5,1.5,1.5,1.5))
    print(ggpubr::ggarrange(p1,p2,p3,p4,p5, ncol = 3, nrow = 2, labels=LETTERS[1:5],widths = c(1.5,1.5,1.5,1.5)))


  }else{
    ggpubr::ggarrange(p1,p2)
  }

  return(list(PRS = pheno.prs, BestP = bestP))

  ftemp <- list.files(paste0(ResultDir,"/"),pattern = "pruned_")
  invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))
  ftemp <- list.files(paste0(ResultDir,"/"),pattern = "PRS")
  invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))
  ftemp <- list.files(paste0(ResultDir,"/"),pattern = "Clump")
  invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))
  ftemp <- list.files(paste0(ResultDir,"/"),pattern = "pcfile")
  invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))
  ftemp <- c("range_list","Valid.SNP","SNPdata_1","SNP.pvalue","prssummarystat")
  suppressWarnings(invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp)))))

}


#' MergeRegion: Merging two sets of plink binary files.
#'
#' @description This function combines the two genotype datasets based on either common SNPs or all the SNPs between them.
#' @author Banabithi Bose
#'
#' @param DataDir A character string for the file path of the input plink binary files.
#' @param ResultDir A character string for the file path where all output files will be stored. The default is tempdir().
#' @param finput1 Character string, specifying the prefix of the first input PLINK binary files.
#' @param finput2 Character string, specifying the prefix of the first input PLINK binary files.
#' @param foutput Character string, specifying the prefix of the output PLINK binary files if filtering option for the SNPs is chosen. The default is "FALSE".
#' @param use_common_snps Boolean value, TRUE or FALSE, specifying to use common SNPs for merging or to use all the SNPs.
#'
#' @return NULL
#'
#' The output plink files will be saved in ResultDir.
#'
#' @export
#'
#' @examples
#' # Not Run
# DataDir <- system.file("extdata", package = "GXwasR")
# ResultDir <- tempdir()
# finput1 <- "GXwasR_example"
# finput2 <- "GXwasR_example_imputed"
# foutput <- "Test_output"
#
# # Not Run
# #y <- MergeRegion(DataDir, ResultDir, finput1, finput2, foutput,  use_common_snps = TRUE)

MergeRegion <- function(DataDir, ResultDir, finput1, finput2, foutput, use_common_snps = TRUE){

  if (file.exists(paste0(DataDir, "/", finput1, ".bed")) &&
      file.exists(paste0(DataDir, "/", finput1, ".bim")) &&
      file.exists(paste0(DataDir, "/", finput1, ".fam")) &&
      file.exists(paste0(DataDir, "/", finput2, ".bed")) &&
      file.exists(paste0(DataDir, "/", finput2, ".bim")) &&
      file.exists(paste0(DataDir, "/", finput2, ".fam"))) {

    setupPlink(ResultDir)

  } else{
    writeLines(
      "There is missing plink files in specified DataDir.\nPlease specify correct directory path with input Plink files."
    )
  }

  if (use_common_snps == TRUE){

    bim1 <- read.table(paste0(DataDir, "/", finput1, ".bim"))
    bim2 <- read.table(paste0(DataDir, "/", finput2, ".bim"))
    common_snps <- intersect(bim1$V2,bim2$V2)
    write.table(common_snps,file = paste0(ResultDir, "/common_snps_", foutput),quote = FALSE, row.names = FALSE, col.names = FALSE)

    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--bfile",
        paste0(DataDir,"/",finput1),
        "--extract",
        paste0(ResultDir, "/common_snps_", foutput),
        "--make-bed",
        "--out",
        paste0(ResultDir,"/new",finput1),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))

    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--bfile",
        paste0(DataDir,"/",finput2),
        "--extract",
        paste0(ResultDir, "/common_snps_", foutput),
        "--make-bed",
        "--out",
        paste0(ResultDir,"/new",finput2),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))

    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--bfile",
        paste0(ResultDir,"/new",finput1),
        "--bmerge",
        paste0(ResultDir,"/new",finput2),
        "--make-bed",
        "--out",
        paste0(ResultDir,"/",foutput),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))

    print("Merging is done using the common SNPs between the input genotype files.")
  }else{
    print("Merging is done with all the SNPs i.e., union of the SNPs.")
    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--bfile",
        paste0(DataDir,"/",finput1),
        "--bmerge",
        paste0(DataDir,"/",finput2),
        "--make-bed",
        "--out",
        paste0(ResultDir,"/",foutput),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))
  }

  if (use_common_snps == TRUE){
    ftemp <- list.files(paste0(ResultDir,"/"),pattern = "new")
    invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))

    ftemp <- list.files(paste0(ResultDir,"/"),pattern = "common_snps")
    invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))
  }
  print(paste0("Plink files with merged regions are in ",ResultDir," prefixed as ",foutput))

}

#' plinkVCF: Converting VCF files to plink binary files and vice-versa.
#'
#' @author Banabithi Bose
#'
#' @description This function performs the conversion between VCF files to plink binary formats.
#'
#' For VCF to plink files conversion, if you do not specify any FAM file when you are converting from VCF to plink format, then plink will just create a 'dummy' FAM file with the same name as your dataset with missing phenotypes and missing sex.

#' @param DataDir A character string for the file path of the input plink binary files and all other input files.
#' @param ResultDir A character string for the file path where all output files will be stored. The default is tempdir().
#' @param finput Character string, specifying the prefix of the input plink binary files with both male and female samples. This file needs to be in DataDir.
#' @param foutput Character string, specifying the prefix of the output plink binary files if filtering option for the SNPs is chosen. The default is "FALSE".
#' @param VtoP Boolean value, TRUE or FALSE, specifying the conversion of VCF files to plink binary files or not. The default is TRUE.
#' @param PtoV Boolean value, TRUE or FALSE, specifying the conversion of plink binary files to VCF  files or not. The default is TRUE.
#' @param Famfile Character string, specifying the name of the original .fam file if VtoP was set to be TRUE. This file needs to be in DataDir. The default is NULL.
#' @param PVbyCHR Boolean value, TRUE or FALSE specifying to do the plink to vcf conversion chromosome-wise or not. The default is TRUE.
#' @return NULL
#'
#' The output files will be saved in ResultDir.
#'
#' @export
#'
#' @examples
#' # Not Run
# finput <- "GXwasR_example" #Plink file
# foutput <- "GXwasR_example1"
# DataDir <- system.file("extdata", package = "GXwasR")
# ResultDir <- tempdir()
# PtoV = TRUE
# VtoP = FALSE
# Famfile = NULL
# PVbyCHR = FALSE
# x <- plinkVCF(DataDir, ResultDir, finput, foutput, VtoP, PtoV, Famfile, PVbyCHR)
#'
plinkVCF <- function(DataDir, ResultDir = tempdir(), finput, foutput, VtoP = FALSE, PtoV = TRUE, Famfile = NULL, PVbyCHR = TRUE){

  if (PtoV == TRUE){

    if (file.exists(paste0(DataDir, "/", finput, ".bed")) &&
        file.exists(paste0(DataDir, "/", finput, ".bim")) &&
        file.exists(paste0(DataDir, "/", finput, ".fam"))) {

      setupPlink(ResultDir)

    } else{
      writeLines(
        "There are no Plink files in DataDir.\nPlease specify correct DataDir path with input Plink files."
      )
    }

    ## plink to vcf
    if (PVbyCHR == FALSE){

      invisible(sys::exec_wait(
        paste0(ResultDir,"/","./plink"),
        args = c(
          "--bfile",
          paste0(DataDir,"/",finput),
          "--recode","vcf",
          #"--make-bed",
          "--out",
          paste0(ResultDir,"/",foutput,"_vcf"),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))
      ###################
      #Generate vcf.gz file and its index file vcf.gz.tbi
      utils::download.file(destfile = paste0(ResultDir,"/","bgzip_tabix.zip"),
                           "https://github.com/boseb/bose_binaries/raw/main/bgzip_tabix.zip", quiet = TRUE,
      )

      utils::unzip(paste0(ResultDir,"/","bgzip_tabix.zip"), exdir = ResultDir)

      Sys.chmod(paste0(ResultDir,"/bgzip"), mode = "0777", use_umask = TRUE)
      #Sys.chmod(paste0(ResultDir,"/tabix"), mode = "0777", use_umask = TRUE)
      fn <- paste0(foutput,"_vcf.vcf")
      system(paste0(ResultDir,"/",'./bgzip ', ResultDir,"/",fn))
      #system(paste0(ResultDir,"/",'./tabix -p vcf ', ResultDir,"/",fn,".gz"))
      ftemp <- list.files(paste0(ResultDir,"/"),pattern = "bgzip")
      invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))
      ftemp <- c("tabix","plink")
      invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))
      ftemp <- list.files(paste0(ResultDir,"/"),pattern = "log")
      invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))
      ftemp <- list.files(paste0(ResultDir,"/"),pattern = "hh")
      invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))

    }else{

      bimfile <- read.table(paste0(DataDir, "/", finput, ".bim"))
      #chrnum <- 1:length(unique(bimfile$V1))
      chrnum <- as.numeric(gsub("chr","",unique(bimfile$V1)))
      #chrnum <- 1:3
      funVCFchr <- function(chrnum){
        print(paste0("Running for chromosome = ",chrnum))

        invisible(sys::exec_wait(
          paste0(ResultDir,"/","./plink"),
          args = c(
            "--bfile",
            paste0(DataDir,"/",finput),
            "--chr", chrnum,
            "--recode","vcf",
            #"--make-bed",
            "--out",
            paste0(ResultDir,"/",foutput,"_chr",chrnum),
            "--silent"
          ),
          std_out = FALSE,
          std_err = FALSE
        ))
        ###################
        #Generate vcf.gz file and its index file vcf.gz.tbi
        utils::download.file(destfile = paste0(ResultDir,"/","bgzip_tabix.zip"),
                             "https://github.com/boseb/bose_binaries/raw/main/bgzip_tabix.zip", quiet = TRUE,
        )

        utils::unzip(paste0(ResultDir,"/","bgzip_tabix.zip"), exdir = ResultDir)

        Sys.chmod(paste0(ResultDir,"/bgzip"), mode = "0777", use_umask = TRUE)
        #Sys.chmod(paste0(ResultDir,"/tabix"), mode = "0777", use_umask = TRUE)
        fn <- paste0(foutput,"_chr",chrnum,".vcf")
        system(paste0(ResultDir,"/",'./bgzip ', ResultDir,"/",fn))

        #system(paste0(ResultDir,"/",'./tabix -p vcf ', ResultDir,"/",fn,".gz"))
      }
      invisible(lapply(chrnum,funVCFchr))
      ftemp <- list.files(paste0(ResultDir,"/"),pattern = "bgzip")
      invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))
      ftemp <- list.files(paste0(ResultDir,"/"),pattern = "log")
      invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))
      ftemp <- c("tabix","plink")
      invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))
      ftemp <- list.files(paste0(ResultDir,"/"),pattern = "hh")
      invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))
    }
    #################

  }else if (VtoP == TRUE){

    if (file.exists(paste0(DataDir, "/", finput, ".vcf"))){

      setupPlink(ResultDir)

    } else{
      writeLines(
        "There is missing plink files in specified DataDir.\nPlease specify correct directory path with input Plink files."
      )
    }

    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--vcf",
        paste0(DataDir,"/",finput),
        "--keep-allele-order",
        "--make-bed",
        "--const-fid", 1,
        "--out",
        paste0(ResultDir,"/",foutput),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))

    ## Replacing the vcf.fam with original fam.
    fam <- read.table(paste0(DataDir,"/",Famfile))
    write.table(fam, file = paste0(ResultDir,"/",foutput,".fam"), col.names = FALSE, row.names = FALSE, quote = FALSE)
    invisible(do.call(file.remove,list(paste0(DataDir,"/","plink"))))
  }

  print("Output files are produced in DataDir.")
  return()

}


#' SexCheck: Compare sex assignments in the input plink files with those imputed from X chromosome inbreeding coefficients
#'
#' @author Banabithi Bose
#'
#' @description This function compares sex assignments in the input dataset with those predicted from X chromosome inbreeding coefficients following Purcell et. al.,(2007), and gives the option to convert the sex assignments to the predicted values. Implicitly, this function computes observed and expected autosomal homozygous genotype counts for each sample and reports method-of-moments F coefficient estimates (i.e., observed hom. count - expected count) / (total observations - expected count)). The expected counts will be based on loaded or imputed minor allele frequencies.  Since imputed MAFs are highly inaccurate when there are few samples, the 'compute freq' parameter should be set to TRUE to compute MAF implicitly.
#' Due to the use of allele frequencies, if a cohort is comprised of individuals of different ancestries, users may need to process any samples with rare ancestry individually if the dataset has a very unbalanced ancestry distribution. It is advised to run this function with all the parameters set to zero, then examine the distribution of the F estimates (there should be a clear gap between a very tight male clump on the right side of the distribution and the females everywhere else). Then, rerun the function with the parameters that correspond to this gap.
#'
#' @param DataDir Character string for the file path of the input PLINK binary files.
#' @param ResultDir A character string for the file path where all output files will be stored. The default is tempdir().
#' @param finput Character string, specifying the prefix of the input PLINK binary files. Note: Input dataset should contain X and Y regions.
#' @param impute_sex Boolean value, 'TRUE' or 'FALSE', specifying sex to be imputed or not. If 'TRUE' then sex-imputed PLINK files, prefixed, 'seximputed_plink', will be produced in DataDir.
#' @param compute_freq Boolean value, 'TRUE' or 'FALSE', specifying minor allele frequency (MAF). This function requires reasonable MAF estimates, so it is essential to use compute_freq = TRUE for computing MAF from an input PLINK file if there are very few samples in the input dataset. The default is FALSE.
#' @param LD Boolean value, 'TRUE' or 'FALSE' for applying linkage disequilibrium (LD)-based filtering. The default is TRUE.
#' @param LD_window_size Integer value, specifying a window size in variant count for LD-based filtering. The default is 50.
#' @param LD_step_size Integer value, specifying a variant count to shift the window at the end of each step for LD filtering. The default is 5.
#' @param LD_r2_threshold Numeric value between 0 to 1 of pairwise r^2 threshold for LD-based filtering. The default is 0.02.
#' @param fmax_F Numeric value between 0 to 1. Samples with F estimates smaller than this value will be labeled as females. The default is 0.2.
#' @param mmin_F Numeric value between 0 to 1. Samples with F estimates larger than this value will be labeled as males. The default is 0.8.
#'
#' @return A dataframe with six columns:FID(Family ID), IID(Individual ID), PEDSEX(Sex as determined in pedigree file (1=male, 2=female)), SNPSEX(Sex as determined by X chromosome), STATUS(Displays "PROBLEM" or "OK" for each individual) and F(The actual X chromosome inbreeding (homozygosity) estimate).A PROBLEM arises if the two sexes do not match, or if the SNP data or pedigree data are ambiguous with regard to sex.
#' @export
#'
#' @references (1) Purcell et. al.,(2007). PLINK: A Tool Set for Whole-Genome Association and Population-Based Linkage Analyses. The American Journal of Human Genetics, 81(3), 559-575. https://doi.org/10.1086/519795.
#'
#' @examples
#' # Not Run
# DataDir <- system.file("extdata", package = "GXwasR")
# ResultDir = tempdir()
# finput <- "GXwasR_example"
# LD = TRUE
# LD_window_size <- 50
# LD_step_size <- 5
# LD_r2_threshold <- 0.02
# fmax_F <- 0.2
# mmin_F <- 0.8
# impute_sex <- FALSE
# compute_freq <- FALSE
#
# x <-SexCheck(DataDir=DataDir,ResultDir = ResultDir, finput=finput,impute_sex=impute_sex,compute_freq =compute_freq,LD_window_size=LD_window_size,LD_step_size=LD_step_size,LD_r2_threshold=0.02,fmax_F=0.2,mmin_F=0.8)
#
# # Checking if there is any wrong sex assignment
# # Not Run
# # problematic_sex <- x[x$STATUS != "OK",]

SexCheck <-
  function(DataDir,
           ResultDir = tempdir(),
           finput,
           impute_sex = FALSE,
           compute_freq = FALSE,
           LD = TRUE,
           LD_window_size = 50,
           LD_step_size = 5,
           LD_r2_threshold = 0.02,
           fmax_F = 0.2 ,
           mmin_F = 0.8) {
    if (file.exists(paste0(DataDir, "/", finput, ".bed")) &&
        file.exists(paste0(DataDir, "/", finput, ".bim")) &&
        file.exists(paste0(DataDir, "/", finput, ".fam"))) {

      setupPlink(ResultDir)

    } else{
      writeLines(
        "There are no Plink files in DataDir.\nPlease specify correct directory path with input Plink files."
      )
      #return()
    }

    ## Checking for XY regions
    bim <-
      as.data.frame(utils::read.table(file = paste0(DataDir, "/", finput, ".bim")))
    xChr <- nrow(bim[bim$V1 == 23 | bim$V1 == "X", ])
    yChr <- nrow(bim[bim$V1 == 24 | bim$V1 == "Y", ])

    if (xChr == 0 | yChr == 0) {
      print(
        "There are no X and Y chromosomes in the input PLINK files. Please provide input with pre-existing XY regions."
      )
    } else {
      # cdir <- getwd()
      # if (cdir != DataDir){
      #
      #   #setwd(DataDir)
      # }else{
      #   print("Working directrory is same as DataDir")
      # }

      if (impute_sex == FALSE) {
        if (compute_freq == TRUE) {
          if (LD == TRUE) {
            invisible(sys::exec_wait(
              paste0(ResultDir,"/","./plink"),
              args = c(
                "--bfile",
                paste0(DataDir,"/",finput),
                "--freq",
                "--out",
                paste0(ResultDir,"/","freq_file"),
                "--silent"
              ),
              std_out = FALSE,
              std_err = FALSE
            ))

            invisible(sys::exec_wait(
              paste0(ResultDir,"/","./plink"),
              args = c(
                "--bfile",
                paste0(DataDir,"/",finput),
                "--indep-pairwise",
                LD_window_size,
                LD_step_size,
                LD_r2_threshold,
                "--read-freq",
                paste0(ResultDir,"/","freq_file.frq"),
                "--check-sex",
                fmax_F,
                mmin_F,
                "--out",
                paste0(ResultDir,"/","cs"),
                "--silent"
              ),
              std_out = FALSE,
              std_err = FALSE
            ))

            check_sex <-
              read.table(
                file = paste0(ResultDir,"/", "cs.sexcheck"),
                stringsAsFactors = FALSE,
                header = TRUE
              )


          } else{
            invisible(sys::exec_wait(
              paste0(ResultDir,"/","./plink"),
              args = c(
                "--bfile",
                paste0(DataDir,"/",finput),
                "--freq",
                "--out",
                paste0(ResultDir,"/","freq_file"),
                "--silent"
              ),
              std_out = FALSE,
              std_err = FALSE
            ))

            invisible(sys::exec_wait(
              paste0(ResultDir,"/","./plink"),
              args = c(
                "--bfile",
                paste0(DataDir,"/",finput),
                "--read-freq",
                paste0(ResultDir,"/","freq_file.frq"),
                "--check-sex",
                fmax_F,
                mmin_F,
                "--out",
                paste0(ResultDir,"/","cs"),
                "--silent"
              ),
              std_out = FALSE,
              std_err = FALSE
            ))

            check_sex <-
              read.table(
                file = paste0(ResultDir,"/", "cs.sexcheck"),
                stringsAsFactors = FALSE,
                header = TRUE
              )

          }
        } else if (compute_freq == FALSE) {
          if (LD == TRUE) {
            invisible(sys::exec_wait(
              paste0(ResultDir,"/","./plink"),
              args = c(
                "--bfile",
                paste0(DataDir,"/",finput),
                "--indep-pairwise",
                LD_window_size,
                LD_step_size,
                LD_r2_threshold,
                "--check-sex",
                fmax_F,
                mmin_F,
                "--out",
                paste0(ResultDir,"/","cs"),
                "--silent"
              ),
              std_out = FALSE,
              std_err = FALSE
            ))

            check_sex <-
              read.table(
                file = paste0(ResultDir,"/", "cs.sexcheck"),
                stringsAsFactors = FALSE,
                header = TRUE
              )

          } else{
            invisible(sys::exec_wait(
              paste0(ResultDir,"/","./plink"),
              args = c(
                "--bfile",
                paste0(DataDir,"/",finput),
                "--check-sex",
                fmax_F,
                mmin_F,
                "--out",
                paste0(ResultDir,"/","cs"),
                "--silent"
              ),
              std_out = FALSE,
              std_err = FALSE
            ))

            check_sex <-
              read.table(
                file = paste0(ResultDir,"/", "cs.sexcheck"),
                stringsAsFactors = FALSE,
                header = TRUE
              )
          }


        }
      } else if (impute_sex == TRUE) {
        if (LD == TRUE) {
          invisible(sys::exec_wait(
            paste0(ResultDir,"/","./plink"),
            args = c(
              "--bfile",
              paste0(DataDir,"/",finput),
              "--indep-pairwise",
              LD_window_size,
              LD_step_size,
              LD_r2_threshold,
              "--make-bed",
              "--out",
              paste0(ResultDir,"/","csLD"),
              "--silent"
            ),
            std_out = FALSE,
            std_err = FALSE
          ))
          invisible(sys::exec_wait(
            paste0(ResultDir,"/","./plink"),
            args = c(
              "--bfile",
              paste0(ResultDir,"/","csLD"),
              "--impute-sex",
              "--make-bed",
              "--out",
              paste0(ResultDir,"/","seximputed_plink"),
              "--silent"
            ),
            std_out = FALSE,
            std_err = FALSE

          ))
          check_sex <-
            read.table(
              file = paste0(ResultDir,"/","seximputed_plink.sexcheck"),
              stringsAsFactors = FALSE,
              header = TRUE
            )
          print(
            "The output plink files with imputed sex, prefixed, seximputed_plink, is in the working directory."
          )

        } else{
          invisible(sys::exec_wait(
            paste0(ResultDir,"/","./plink"),
            args = c(
              "--bfile",
              paste0(DataDir,"/",finput),
              "--impute-sex",
              "--make-bed",
              "--out",
              paste0(ResultDir,"/","seximputed_plink"),
              "--silent"
            ),
            std_out = FALSE,
            std_err = FALSE

          ))
          check_sex <-
            read.table(
              file = paste0(ResultDir,"/", "seximputed_plink.sexcheck"),
              stringsAsFactors = FALSE,
              header = TRUE
            )

          print(
            "The output plink files with imputed sex, prefixed as seximputed_plink, is in ResultDir."
          )
        }


      }

      ftemp <- list.files(paste0(ResultDir,"/"),pattern = "cs")
      invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))

      return(check_sex)
    }
  }


#' FilterPlinkSample: Making PLINK files with desired samples.
#'
#' @author Banabithi Bose
#'
#' @description This function prepares PLINK binary files with the desired samples.
#'
#' @param DataDir Character string for the file path of the all input files.
#' @param ResultDir character string for the file path where the output PLINK files will be stored.
#' @param finput Character string, specifying the prefix of the input PLINK binary files.
#' @param foutput Character string, specifying the prefix of the output PLINK binary files.
#' @param filter_sample Character string, specifying the sample type to be retained. The choices are, "cases","controls","males" and "females". The default is "cases".
#' @param keep_remove_sample_file Character string, specifying the prefix of a space/tab-delimited text file with no header. For the samples that we want to keep or remove, the family IDs should be in the first column and within-family IDs in the second column. This file needs to be in the DataDir. The default is NULL.
#' @param keep Boolean value, "TRUE" or "FALSE" for specifying desired samples to keep or remove. The default is "TRUE".
#'
#' @return NULL
#' The output plink files with passed samples will be saved in ResultDir.
#'
#' @export
#'
#' @examples
#' # Not Run
# DataDir <- system.file("extdata", package = "GXwasR")
# ResultDir <- tempdir()
# finput <- "GXwasR_example"
# foutput <- "casesPlink"
# filter_sample <- "cases"
# keep_remove_sample_file <- "samples_example"
# keep <- FALSE
#
# FilterPlinkSample(DataDir = DataDir,ResultDir = ResultDir,finput = finput,foutput = foutput,keep_remove_sample_file = keep_remove_sample_file, keep = keep)
#
FilterPlinkSample <- function(DataDir,ResultDir,
                              finput,
                              foutput = NULL,
                              filter_sample = c("cases","controls","males","females"),
                              keep_remove_sample_file = NULL,
                              keep = TRUE){

  if (file.exists(paste0(DataDir, "/", finput, ".bed")) &&
      file.exists(paste0(DataDir, "/", finput, ".bim")) &&
      file.exists(paste0(DataDir, "/", finput, ".fam"))) {

    setupPlink(ResultDir)

  } else{
    writeLines(
      "There are no Plink files in DataDir.\nPlease specify correct directory path with input Plink files."
    )
    #return()
  }

  if (is.null(keep_remove_sample_file)){

    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--bed",
        paste0(DataDir,"/",finput, ".bed"),
        "--bim",
        paste0(DataDir,"/",finput, ".bim"),
        "--fam",
        paste0(DataDir,"/",finput, ".fam"),
        paste0("--filter-", filter_sample),
        "--make-bed",
        "--out",
        paste0(ResultDir,"/",foutput),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))

  }else{

    if (keep == TRUE){

      invisible(sys::exec_wait(
        paste0(ResultDir,"/","./plink"),
        args = c(
          "--bed",
          paste0(DataDir,"/",finput, ".bed"),
          "--bim",
          paste0(DataDir,"/",finput, ".bim"),
          "--fam",
          paste0(DataDir,"/",finput, ".fam"),
          "--keep", paste0(DataDir,"/",keep_remove_sample_file),
          "--make-bed",
          "--out",
          paste0(ResultDir,"/",foutput),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))
    }else if (keep == FALSE){

      invisible(sys::exec_wait(
        paste0(ResultDir,"/","./plink"),
        args = c(
          "--bed",
          paste0(DataDir,"/",finput, ".bed"),
          "--bim",
          paste0(DataDir,"/",finput, ".bim"),
          "--fam",
          paste0(DataDir,"/",finput, ".fam"),
          "--remove", paste0(DataDir,"/",keep_remove_sample_file),
          "--make-bed",
          "--out",
          paste0(ResultDir,"/",foutput),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))
    }
  }
  print(paste0(foutput," plink files with desired samples are in ", ResultDir))
  return()
}

#' GetMFPlink: Getting male and female PLINK binary files.
#'
#' @author Banabithi Bose
#'
#' @description This function prepares separate male and female PLINK binary files from combined PLINK files.

#'
#' @param DataDir Character string for the file path of the input PLINK binary files.
#' @param ResultDir A character string for the file path where all output files will be stored. The default is tempdir().
#' @param finput Character string, specifying the prefix of the input PLINK binary files.
#' @param foutput Character string, specifying the prefix of the output PLINK binary files.
#' @param sex Boolean value, 'males' or 'females', specifying output plink binary files with male or female samples.
#' @param xplink Boolean value, 'TRUE' or 'FALSE', specifying output plink binary files with only X chromosome or not. Default is FALSE.
#' @param autoplink Boolean value, 'TRUE' or 'FALSE', specifying output plink binary files with only autosome or not. Default is FALSE.
#'
#' @return None
#' @export
#'
#' @examples
#' # Not Run
# DataDir <- system.file("extdata", package = "GXwasR")
# ResultDir <- tempdir()
# finput <- "GXwasR_example"
# foutput <- "Test_output"
# sex <- "females"
# x <- GetMFPlink(DataDir = DataDir, ResultDir = ResultDir, finput = finput, foutput = foutput,sex = sex, xplink = FALSE, autoplink = FALSE)

GetMFPlink <-
  function(DataDir,
           ResultDir = tempdir(),
           finput,
           foutput,
           sex,
           xplink = FALSE,
           autoplink = FALSE) {
    if (file.exists(paste0(DataDir, "/", finput, ".bed")) &&
        file.exists(paste0(DataDir, "/", finput, ".bim")) &&
        file.exists(paste0(DataDir, "/", finput, ".fam"))) {

      setupPlink(DataDir)

    } else{
      writeLines(
        "There are no Plink files in DataDir.\nPlease specify correct directory path with input Plink files."
      )

    }

    if (xplink == FALSE && autoplink == FALSE) {
      invisible(sys::exec_wait(
        paste0(DataDir,"/","./plink"),
        args = c(
          "--bed",
          paste0(DataDir,"/",finput, ".bed"),
          "--bim",
          paste0(DataDir,"/",finput, ".bim"),
          "--fam",
          paste0(DataDir,"/",finput, ".fam"),
          paste0("--filter-", sex),
          "--make-bed",
          "--out",
          paste0(ResultDir,"/",foutput),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))
    } else if (xplink == TRUE && autoplink == FALSE) {

      invisible(sys::exec_wait(
        paste0(DataDir,"/","./plink"),
        args = c(
          "--bed",
          paste0(DataDir,"/",finput, ".bed"),
          "--bim",
          paste0(DataDir,"/",finput, ".bim"),
          "--fam",
          paste0(DataDir,"/",finput, ".fam"),
          paste0("--filter-", sex),
          "--chr",
          23,
          "--make-bed",
          "--out",
          paste0(ResultDir,"/",foutput),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))

    } else if (xplink == FALSE && autoplink == TRUE) {
      invisible(sys::exec_wait(
        paste0(DataDir,"/","./plink"),
        args = c(
          "--bed",
          paste0(DataDir,"/",finput, ".bed"),
          "--bim",
          paste0(DataDir,"/",finput, ".bim"),
          "--fam",
          paste0(DataDir,"/",finput, ".fam"),
          paste0("--filter-", sex),
          "--not-chr",
          23,
          "--make-bed",
          "--out",
          paste0(ResultDir,"/",foutput),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))
    }

    print(paste0("Output plink files, prefixed as ",foutput,", are in ",ResultDir))
    return()
  }


#' Xhwe: Filter X-chromosome variants for HWE in females.
#'
#' @author Banabithi Bose
#'
#' @description This function is a part of the post-imputation quality control process prior to GWAS. This tests for Hardy-Weinberg Equilibrium (HWE) for X-chromosome variants in females. Males' hemizygous X chromosome prevents testing for HWE on their haploid X calls, and testing for HWE across all samples would have a high failure rate. This function will check for HWE across the X in females (cases and controls combined), following the recommendation in Khramtsova et al., 2023, and can remove these regions from analysis in all samples. The p-value threshold for filtering out SNPs is 0.05/no.of. X-chromosome variants.
#'
#' @param DataDir A character string for the file path of the input PLINK binary files.
#' @param ResultDir A character string for the file path where all output files will be stored. The default is tempdir().
#' @param finput Character string, specifying the prefix of the input PLINK binary files.
#' @param foutput Character string, specifying the prefix of the output PLINK binary files if filtering option for the SNPs is chosen. The default is "FALSE".
#' @param filterSNP Boolean value, 'TRUE' or 'FALSE' for filtering out the X-chromosome variants i.e., SNPs from the input file or not (i.e., only flagged). The default is "FALSE".

#' @return A list object containing SNPs. If filterSNP = TRUE, the output filtered PLINK binary files will be produced inside DataDir.
#' @export
#'
#' @examples
#' # Not Run
# DataDir <- system.file("extdata", package = "GXwasR")
# ResultDir <- tempdir()
# finput <- "GXwasR_example"
# foutput <- "Test_output"
# x <- Xhwe(DataDir = DataDir, ResultDir = ResultDir, finput = finput, foutput = foutput, filterSNP = TRUE)
# # Failed SNPs for HWE in X chromosome
# x
Xhwe <- function(DataDir, ResultDir = tempdir(), finput, filterSNP = TRUE, foutput){

  if (file.exists(paste0(DataDir, "/", finput, ".bed")) &&
      file.exists(paste0(DataDir, "/", finput, ".bim")) &&
      file.exists(paste0(DataDir, "/", finput, ".fam"))) {

    setupPlink(ResultDir)

  } else{
    writeLines(
      "There are no Plink files in DataDir.\nPlease specify correct directory path with input Plink files."
    )
  }

  ## Getting female Plink file
  invisible(GetMFPlink(DataDir = DataDir, ResultDir,finput = finput, foutput = "female",sex = "females", xplink = FALSE, autoplink = FALSE))

  fam <-
    as.data.frame(utils::read.table(file = paste0(ResultDir, "/", "female.fam")))

  # Check for case control status in input plink file
  fam <- na.omit(fam)
  fam$V6 <- as.numeric(as.character(fam$V6))
  fam2 <- fam[fam$V6 != 0,]
  fam3 <- fam2[fam2$V6 != -9,]

  if (length(unique(fam3$V6)) != 2) {
    writeLines(
      "There is not both case-control status for females in input Plink files. This test cannot be performed."
    )
  } else{

    #cdir <- getwd()
    #setwd(DataDir)

    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--bfile",
        paste0(ResultDir,"/","female"),
        "--chr", 23,
        "--hardy",
        "--out",
        paste0(ResultDir,"/","Xhwe"),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))


    x <-
      as.data.frame(
        read.table(file = paste0(ResultDir,"/", "Xhwe.hwe"),
                   header = TRUE,
                   sep = ""
        )
      )

    ## Bonferroni-corrected pvalue threshold set to 0.05/(number of X chromosome variants)
    p <- 0.05/length(unique(x$SNP))
    #system(paste0("awk -v bf=",p," '$9<bf {print $2}' Xhwe.hwe > Xsnps_to_exclude"))
    snp <- x[x$P < p, 2,drop = TRUE]
    X_excluded_SNPs <- unique(snp)

    if (length(X_excluded_SNPs) == 0) {
      print("No SNP to be excluded.")
      return()
    } else{
      if (length(X_excluded_SNPs) != 0 & filterSNP == "TRUE") {

        utils::write.table(
          X_excluded_SNPs,
          file = paste0(ResultDir,"/", "XhweSNPs"),
          quote = FALSE,
          row.names = FALSE,
          col.names = FALSE,
          eol = "\r\n",
          sep = " "
        )


        invisible(sys::exec_wait(
          paste0(ResultDir,"/","./plink"),
          args = c(
            "--bfile",
            paste0(DataDir,"/",finput),
            "--exclude",
            paste0(ResultDir,"/","XhweSNPs"),
            "--make-bed",
            "--out",
            paste0(ResultDir,"/",foutput),
            "--silent"
          ),
          std_out = FALSE,
          std_err = FALSE
        ))

        #setwd(cdir)
        print(
          paste0("Failed SNPs are excluded from the output plink files prefixed as ",foutput," is in ",ResultDir)
        )
        #setwd(cdir)
        ftemp <- list.files(paste0(ResultDir,"/"),pattern = "hwe")
        invisible(file.remove(paste0(ResultDir,"/",ftemp)))
        ftemp <- list.files(paste0(ResultDir,"/"),pattern = "female")
        invisible(file.remove(paste0(ResultDir,"/",ftemp)))
        return(X_excluded_SNPs)

      } else{
        #setwd(cdir)
        print("SNPs are flagged.")
        ftemp <- list.files(paste0(ResultDir,"/"),pattern = "hwe")
        invisible(file.remove(paste0(ResultDir,"/",ftemp)))
        ftemp <- list.files(paste0(ResultDir,"/"),pattern = "female")
        invisible(file.remove(paste0(ResultDir,"/",ftemp)))
        return(X_excluded_SNPs)

      }

    }
  }


}

#' QCsample2: Post-imputation quality control for samples in the plink binary files.
#'
#' @author Banabithi Bose
#'
#' @description This function finds and/or filters the samples with outlying heterozygosity, missing genotype rates and Identity By Descent (IBD).
#'
#' @param DataDir A character string for the file path of the input PLINK binary files.
#' @param ResultDir A character string for the file path where all output files will be stored. The default is tempdir().
#' @param finput Character string, specifying the prefix of the input PLINK binary files with both male and female samples. This file needs to be in DataDir.
#' @param foutput Character string, specifying the prefix of the output PLINK binary files if filtering option for the samples is chosen.
#' @param miss_geno Numeric value between 0 to 1 for removing samples that have more than the specified missing genotype call rate. The default is 0.1.
#' @param het_Fstat Numeric value between 0 to 1 for removing samples that have more than the specified absolute value of the heterozygosity F statistic. The default is 0.2.
#' @param small_sample_mod Boolean value, TRUE or FALSE, specifying whether to use small sample modifier in heterozygocity computation or not. The default is FALSE. Allele frequencies are used in this calculation. Ideally, the sample size should not be small.
#' @param IBD Numeric value between 0 to 1 for removing samples that have more than the specified pi-hat value as a measure of IBD. The default is 0.2.
#' @param keep Boolean value, TRUE or FALSE, specifying whether to filter out the failed samples or just flag them. The default is FALSE.
#' @param ld_prunning Boolean value, 'TRUE' or 'FALSE' for applying linkage disequilibrium (LD)-based filtering for IBD statistics computation. The default is FALSE. It is usually a good idea to perform some form of LD-based pruning before computing IBD statistics.
#' @param highLD_regions NULL or Character string, specifying the .txt file name with genomic regions with high LD. Based on a genome build Hg19 and Hg38, two files such as "high-LD-regions-hg19-GRCh37.txt" and "" are provided in "extdata" folder.
#' @param window_size Integer value, specifying a window size in variant count or kilobase for LD-based filtering. The default is 50.
#' @param step_size Integer value, specifying a variant count to shift the window at the end of each step for LD filtering. The default is 5.
#' @param r2_threshold Numeric value between 0 to 1 of pairwise r^2 threshold for LD-based filtering. The default is 0.02.

#' @return A list of three dataframes, namely, Failedheter, Failedmissingness and Failedibd, containing samples that failed heterozygosity test, genotype call rate and IBD test, respectively. Each of these dataframes will have two columns, containing family id (FID) and individual id (IID).
#'
#' If keep is FALSE, then the filtered output plink files with passed sample will be saved in the ResultDir.
#'
#' @export
#'
#' @examples
#' # Not Run
# # DataDir <- system.file("extdata", package = "GXwasR")ResultDir <- tempdir()
# finput <- "GXwasR_example"
# foutput <- "Test_output"
# miss_geno = 0.1
# het_Fstat = 0.2
# small_sample_mod = TRUE
# IBD = 0.2
# ld_prunning = FALSE
# keep = FALSE
#
# # Running the function
# x <- QCsample2(DataDir = DataDir,ResultDir = ResultDir, finput = finput,foutput = foutput,miss_geno = miss_geno,
#  het_Fstat = het_Fstat, small_sample_mod = small_sample_mod, IBD = IBD, ld_prunning = ld_prunning, keep = keep)

QCsample2 <- function(DataDir,
                      ResultDir = tempdir(),
                      finput,
                      foutput,
                      miss_geno = 0.1,
                      het_Fstat = 0.2,
                      small_sample_mod = FALSE,
                      IBD = 0.2,
                      ld_prunning = FALSE,
                      highLD_regions,
                      window_size = 50,
                      step_size = 5,
                      r2_threshold = 0.02,
                      keep = FALSE){


  if (file.exists(paste0(DataDir, "/", finput, ".bed")) &&
      file.exists(paste0(DataDir, "/", finput, ".bim")) &&
      file.exists(paste0(DataDir, "/", finput, ".fam"))) {

    setupPlink(ResultDir)

  } else{
    writeLines(
      "There are no Plink files in DataDir.\nPlease specify correct directory path with input Plink files."
    )
  }

  if (small_sample_mod == TRUE){
    SSM <- "small-sample"
  }else{
    SSM = NULL
  }

  if (is.null(miss_geno)){
    MIND = NULL
  }else{
    MIND = "--mind"
  }

  if (!is.null(het_Fstat)){

    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--bfile",
        paste0(DataDir,"/",finput),
        "--dog",
        "--het",
        SSM,
        "--out",
        paste0(ResultDir,"/","filtered_hetero"),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))

    heter <- read.table(
      paste0(ResultDir,"/","filtered_hetero.het"),
      sep = "",
      header = TRUE,
      as.is = TRUE,
      stringsAsFactors = FALSE
    )

    if(nrow(heter)==0){

      het_failed = NULL

      invisible(sys::exec_wait(
        paste0(ResultDir,"/","./plink"),
        args = c(
          "--bed",
          paste0(DataDir,"/",finput, ".bed"),
          "--bim",
          paste0(DataDir,"/",finput, ".bim"),
          "--fam",
          paste0(DataDir,"/",finput, ".fam"),
          "--make-bed",
          "--out",
          paste0(ResultDir,"/filterSam1"),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))

    }else{

      heter$F <- as.numeric(as.character(heter$F))
      samples_failed_heter <- heter[heter$F>het_Fstat,1:2]

      write.table(samples_failed_heter, file = paste0(ResultDir,"/","samples_failed_heter"),
                  quote = FALSE, row.names = FALSE, col.names = FALSE)

      het_failed <- samples_failed_heter


      if (nrow(het_failed) != 0){

        het_failed = NULL

        if (keep == FALSE){

          invisible(sys::exec_wait(
            paste0(ResultDir,"/","./plink"),
            args = c(
              "--bed",
              paste0(DataDir,"/",finput, ".bed"),
              "--bim",
              paste0(DataDir,"/",finput, ".bim"),
              "--fam",
              paste0(DataDir,"/",finput, ".fam"),
              "--remove", paste0(ResultDir,"/samples_failed_heter"),
              "--make-bed",
              "--out",
              paste0(ResultDir,"/filterSam1"),
              "--silent"
            ),
            std_out = FALSE,
            std_err = FALSE
          ))

        }else{


          invisible(sys::exec_wait(
            paste0(ResultDir,"/","./plink"),
            args = c(
              "--bed",
              paste0(DataDir,"/",finput, ".bed"),
              "--bim",
              paste0(DataDir,"/",finput, ".bim"),
              "--fam",
              paste0(DataDir,"/",finput, ".fam"),
              "--make-bed",
              "--out",
              paste0(ResultDir,"/filterSam1"),
              "--silent"
            ),
            std_out = FALSE,
            std_err = FALSE
          ))
        }
      }else{
        invisible(sys::exec_wait(
          paste0(ResultDir,"/","./plink"),
          args = c(
            "--bed",
            paste0(DataDir,"/",finput, ".bed"),
            "--bim",
            paste0(DataDir,"/",finput, ".bim"),
            "--fam",
            paste0(DataDir,"/",finput, ".fam"),
            "--make-bed",
            "--out",
            paste0(ResultDir,"/filterSam1"),
            "--silent"
          ),
          std_out = FALSE,
          std_err = FALSE
        ))
      }
    }

  }else{

    het_failed = NULL

    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--bed",
        paste0(DataDir,"/",finput, ".bed"),
        "--bim",
        paste0(DataDir,"/",finput, ".bim"),
        "--fam",
        paste0(DataDir,"/",finput, ".fam"),
        "--make-bed",
        "--out",
        paste0(ResultDir,"/filterSam1"),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))

  }

  if (keep == FALSE){

    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--bfile",
        paste0(ResultDir,"/filterSam1"),
        MIND,
        miss_geno,
        "--make-bed",
        "--out",
        paste0(ResultDir,"/","filterSam2"),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))

  }else{
    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--bfile",
        paste0(ResultDir,"/filterSam1"),
        "--make-bed",
        "--out",
        paste0(ResultDir,"/","filterSam2"),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))
  }

  fam1 <- read.table(paste0(ResultDir,"/filterSam1.fam"))
  fam2 <- read.table(paste0(ResultDir,"/filterSam2.fam"))

  failed_missingness1 <- as.data.frame(outersect(fam1$V2,fam2$V2))

  if(nrow(failed_missingness1)==0){
    failed_missingness = NULL
  }else{
    #dt[dt$fct %in% vc,]
    #failed_missingness <-fam1[failed_missingness1$`outersect(fam1$V2, fam2$V2)`,1:2]
    failed_missingness <-fam1[fam1$V2 %in% failed_missingness1$`outersect(fam1$V2, fam2$V2)`,1:2]
  }
  ### LD
  if (ld_prunning == TRUE) {
    excluderange <- "--exclude"

    if (!is.null(highLD_regions)){
      highLD_regions <- paste0(DataDir,"/",highLD_regions)
    }else{highLD_regions = NULL}

    indep <- "--indep-pairwise"
    window_size <- window_size
    step_size <- step_size
    r2_threshold <- r2_threshold

  }else if (ld_prunning == FALSE){
    excluderange <- NULL
    highLD_regions <- NULL
    indep <- NULL
    window_size <- NULL
    step_size <- NULL
    r2_threshold <- NULL
  }


  ###
  if (!is.null(IBD)){

    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--bfile",
        paste0(ResultDir,"/","filterSam2"),
        excluderange,
        highLD_regions,
        indep,
        window_size,
        step_size,
        r2_threshold,
        "--genome",
        "--out",
        paste0(ResultDir,"/","filtered_ibd"),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))

    ibd <- read.table(
      paste0(ResultDir,"/","filtered_ibd.genome"),
      sep = "",
      header = TRUE,
      as.is = TRUE,
      stringsAsFactors = FALSE
    )

    ibd$PI_HAT <- as.numeric(as.character(ibd$PI_HAT))
    samples_failed_ibd1 <- na.omit(ibd[ibd$PI_HAT>IBD,1:4,drop = FALSE])
    samples_failed_ibd2 <- as.data.frame(samples_failed_ibd1$IID1)
    colnames(samples_failed_ibd2) <- "V1"
    samples_failed_ibd3 <- as.data.frame(samples_failed_ibd1$IID2)
    colnames(samples_failed_ibd3) <- "V1"
    ibd_samples <- unique(rbind(samples_failed_ibd2,samples_failed_ibd3))

    fam2 <- read.table(paste0(ResultDir,"/filterSam2.fam"))

    failed_ibd1 <- as.data.frame(ibd_samples)


    if(nrow(failed_ibd1)==0){
      failed_ibd = NULL

      invisible(sys::exec_wait(
        paste0(ResultDir,"/","./plink"),
        args = c(
          "--bfile",
          paste0(ResultDir,"/","filterSam2"),
          "--make-bed",
          "--out",
          paste0(ResultDir,"/",foutput),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))

    }else{

      #failed_ibd <-fam1[failed_ibd1$V1,1:2]
      failed_ibd <-fam2[fam2$V2 %in% failed_ibd1$V1,1:2]

      write.table(failed_ibd, file = paste0(ResultDir,"/","samples_failed_ibd"),
                  quote = FALSE, row.names = FALSE, col.names = FALSE)

      if (keep == FALSE){

        invisible(sys::exec_wait(
          paste0(ResultDir,"/","./plink"),
          args = c(
            "--bed",
            paste0(ResultDir,"/filterSam2.bed"),
            "--bim",
            paste0(ResultDir,"/filterSam2.bim"),
            "--fam",
            paste0(ResultDir,"/filterSam2.fam"),

            "--remove", paste0(ResultDir,"/samples_failed_ibd"),
            "--make-bed",
            "--out",
            paste0(ResultDir,"/",foutput),
            "--silent"
          ),
          std_out = FALSE,
          std_err = FALSE
        ))

      }else{


        invisible(sys::exec_wait(
          paste0(ResultDir,"/","./plink"),
          args = c(
            "--bed",
            paste0(ResultDir,"/filterSam2.bed"),
            "--bim",
            paste0(ResultDir,"/filterSam2.bim"),
            "--fam",
            paste0(ResultDir,"/filterSam2.fam"),
            "--make-bed",
            "--out",
            paste0(ResultDir,"/",foutput),
            "--silent"
          ),
          std_out = FALSE,
          std_err = FALSE
        ))


      }
    }
  }else{


    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--bfile",
        paste0(ResultDir,"/","filterSam2"),
        "--make-bed",
        "--out",
        paste0(ResultDir,"/",foutput),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))
  }

  fam1 <- read.table(paste0(ResultDir,"/filterSam2.fam"))
  fam2 <- read.table(paste0(ResultDir,"/",foutput,".fam"))

  failed_ibd1 <- as.data.frame(outersect(fam1$V2,fam2$V2))


  if(nrow(failed_ibd1)==0){
    failed_ibd = NULL
  }else{
    #failed_ibd <-fam1[failed_ibd1$`outersect(fam1$V2, fam2$V2)`,1:2]
    failed_ibd <-fam1[fam1$V2 %in% failed_ibd1$`outersect(fam1$V2, fam2$V2)`,1:2]
  }

  if (nrow(het_failed)!=0 | !is.null(failed_missingness) | !is.null(failed_ibd)){
    if(keep == FALSE){

      print(paste0("The filtered plink files with passed samples are in ",ResultDir," prefixed as ", foutput))

    }else{

      print("The failed samples are flagged.")
      ftemp <- list.files(paste0(ResultDir,"/"),pattern = paste0(foutput))
      invisible(file.remove(paste0(ResultDir,"/",ftemp)))

    }
  }else {
    print("There is no outlier sample from these tests.")
  }

  ftemp <- list.files(paste0(ResultDir,"/"),pattern = "filter")
  invisible(file.remove(paste0(ResultDir,"/",ftemp)))

  if (file.exists(paste0(ResultDir,"/","samples_failed_ibd"))){
    invisible(file.remove(paste0(ResultDir,"/","samples_failed_ibd")))
  }
  if (file.exists(paste0(ResultDir,"/","samples_failed_heter"))){
    invisible(file.remove(paste0(ResultDir,"/","samples_failed_heter")))
  }

  print(paste0("Number of samples failed for heterogygosity test:",nrow(het_failed)))
  print(paste0("Number of samples failed for missingness test:",nrow(as.data.frame(failed_missingness))))
  print(paste0("Number of samples failed for IBD test:",nrow(as.data.frame(failed_ibd))))
  return(list(Failedheter = het_failed, Failedmissingness = failed_missingness, Failedibd = failed_ibd))

}


#' MAFdiffSexControl: Test for significantly different minor allele frequency (MAF) between sexes in control samples
#'
#' @author Banabithi Bose
#'
#' @description With parameters to filter out SNPs and/or flag the SNPs, this function tests for significantly different MAF (p-value < 0.05/no. of SNPs) between sexes in control samples solely for binary phenotypes. Since the disparities may be caused by technical confounding or sample biases for the research cohorts, it is advised that any SNPs in the controls with a sex difference in MAF be carefully evaluated and identified for further examination (Khramtsova et. al., 2023). In autosomal allele frequencies, sex differences are not anticipated.
#'
#' @param DataDir A character string for the file path of the input PLINK binary files.
#' @param ResultDir A character string for the file path where all output files will be stored. The default is tempdir().
#' @param finput Character string, specifying the prefix of the input PLINK binary files with both male and female samples. This file needs to be in DataDir.
#' @param foutput Character string, specifying the prefix of the output PLINK binary files if filtering option for the SNPs is chosen. The default is NULL.
#' @param filterSNP Boolean value, 'TRUE' or 'FALSE' for filtering out the SNPs or not (i.e., only flagged). The default is "FALSE".
#'
#' @return A list object containing excluded or flagged SNPs. If filterSNP = TRUE, the output filtered PLINK binary files will be produced inside DataDir.
#' @importFrom stats na.omit
#' @importFrom utils download.file read.table unzip write.table
#' @importFrom sys exec_wait
#' @export
#'
#' @examples
#' # Not Run
# DataDir <- system.file("extdata", package = "GXwasR")
# ResultDir <- tempdir()
# finput <- "GXwasR_example"
# foutput <- "Test_output"
# x <- MAFdiffSexControl(DataDir, ResultDir, finput,filterSNP = TRUE,foutput = foutput)

MAFdiffSexControl <-
  function(DataDir,
           ResultDir = tempdir(),
           finput,
           filterSNP = FALSE,
           foutput = NULL) {

    if (file.exists(paste0(DataDir, "/", finput, ".bed")) &&
        file.exists(paste0(DataDir, "/", finput, ".bim")) &&
        file.exists(paste0(DataDir, "/", finput, ".fam"))) {

      setupPlink(ResultDir)

    } else{
      writeLines(
        "There are no Plink files in DataDir.\nPlease specify correct directory path with input Plink files."
      )
    }

    fam <-
      as.data.frame(utils::read.table(file = paste0(DataDir, "/", finput, ".fam")))

    # Check for sex
    fam$V6 <- as.numeric(as.character(fam$V6))
    fam <- stats::na.omit(fam)
    fam1 <- fam[fam$V5 != 0,]
    fam2 <- fam1[fam1$V6 != 0,]
    fam4 <- fam2[fam2$V6 != -9,]
    #fam4 <- fam3[stats::na.omit(fam3$V6),]

    if (length(unique(fam1$V5)) == 2 &&
        length(unique(fam4$V6)) == 2) {
      # For having phenotype for control sample only
      fam$V7 <- 0
    } else{
      writeLines(
        "There is incorrect male-female sex status or incorrect case-control status in input Plink files.\nNeeds both male and female samples with both case and control status to run this function."
      )
    }

    fam$V7[fam$V5 == 1 & fam$V6 == 1] <- 1 # for male and control
    fam[fam$V5 == 2 &
          fam$V6 == 1, 7] <- 2   # for female and control
    fam[fam$V5 == 1 & fam$V6 == 2, 7] <- -9
    fam[fam$V5 == 2 & fam$V6 == 2, 7] <- -9

    phenofile <- unique(fam[, c(1, 2, 7)])
    utils::write.table(
      phenofile,
      file = paste0(ResultDir, "/", "phenofile"),
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      eol = "\r\n",
      sep = " "
    )

    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--bfile",
        paste0(DataDir,"/",finput),
        "--logistic",
        "--pheno",
        paste0(ResultDir,"/","phenofile"),
        "--out",
        paste0(ResultDir,"/","OUTPUT"),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))

    x <-
      utils::read.table(file = paste0(ResultDir, "/", "OUTPUT.assoc.logistic"),
                        header = TRUE)

    # Filter for SEX chromosome (i.e. Filter for X chromosome)
    #y <- x[x$CHR == 23,c(1,2,9)]
    y <- x[, c(1, 2, 9)]
    y <- na.omit(y)
    bf <-
      0.05 / length(unique(y$SNP)) # taking the SNPs for which we have finite p-values.
    y$P <- as.numeric(as.character(y$P))

    flaggedSnps <- unique(y[y$P < bf, 2, drop = TRUE])

    if (length(flaggedSnps) == 0) {
      print("No SNP to be flagged or excluded.")
      flaggedSnps = NULL


    }else if (length(flaggedSnps) != 0 & filterSNP == TRUE) {

      utils::write.table(
        flaggedSnps,
        file = paste0(ResultDir,"/","flaggedSnpsSexMafDiff"),
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE,
        eol = "\r\n",
        sep = " "
      )


      invisible(sys::exec_wait(
        paste0(ResultDir,"/","./plink"),
        args = c(
          "--bfile",
          paste0(DataDir,"/",finput),
          "--exclude",
          paste0(ResultDir,"/","flaggedSnpsSexMafDiff"),
          "--make-bed",
          "--out",
          paste0(ResultDir,"/",foutput),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))


      writeLines(
        paste0("SNPs with significantly MAF difference are excluded.\nFiltered plink files are saved in ",ResultDir)
      )
      return(as.list(flaggedSnps))

    }else if (length(flaggedSnps) != 0 & filterSNP == FALSE){


      print("SNPs are flagged.")


    }

    ftemp <- list.files(paste0(ResultDir,"/"),pattern = "OUTPUT")
    invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))
    invisible(do.call(file.remove,list(paste0(ResultDir,"/","phenofile"))))

    return(flaggedSnps)

  }

#' QCsample: Quality control for samples in the plink binary files.
#'
#' @description This function identifies outlier individuals for heterozygosity and/or missing genotype rates, which aids in the detection of samples with subpar DNA quality and/or concentration that should be removed from the study. Individuals missing more than 3-7% of their genotype calls are often excluded (1)from the analysis.
#' Having the correct designation of sex is important to obtain accurate genotype rate estimates, or avoid incorrectly removing samples, etc. Details can be accessed from the paper.

#' @param DataDir Character string, specifying the file path of the input PLINK binary files. The default is NULL.
#' @param ResultDir A character string for the file path where all output files will be stored. The default is tempdir().
#' @param finput Character string, specifying the prefix of the input PLINK binary files with both male and female samples. This file needs to be in DataDir.
#' @param foutput Character string, specifying the prefix of the output PLINK binary files if filtering option for the samples is chosen.
#' @param imiss Numeric value between 0 to 1 for removing samples that have more than the specified missingness. The default is 0.03.
#' @param het Positive numeric value, specifying the standard deviation from the mean heterozygosity rate. The samples whose rates are more than the specified sd from the mean heterozygosity rate are removed. The default is 3. With this default value, outlying heterozygosity rates would remove individuals who are three sd away from the mean rate (1).
#' @param title_size Integer, specifying the size of the title of the plot heterozygosity estimate vs missingness across samples.
#' @param legend_text_size Integer, specifying the size for legend text in the plot.
#' @param legend_title_size Integer, specifying the size for the legend title in the plot.
#' @param axis_text_size Integer, specifying the size for axis text in the plot.
#' @param axis_title_size Integer, specifying the size for the axis title in the plot.
#' @param filterSample Boolean value, 'TRUE' or 'FALSE' for filtering out the samples or not (i.e., only flagged). The default is "TRUE".
#' @importFrom stats sd
#' @importFrom ggplot2 ggplot
#'
#' @return A plot of heterogysity estimate vs missingness accross sample and a list containing five R dataframe objects, namely, HM (samples with outlying heterozygosity and/or missing genotype rates), Failed_Missingness (samples with missing genotype rates), Failed_heterozygosity (samples with outlying heterozygosity), Missingness_results (missingness results) and Heterozygosity_results (heterozygosity results) with output plink files in ResultDir if filtering out the samples option is chosen.
#' Missingness_results contains missingness results for each individual, with six columns as FID, IID, MISS_PHENO, N_MISS, N_GENO and F_MISS for Family ID, Within-family ID, Phenotype missing? (Y/N), Number of missing genotype call(s), not including obligatory missings or heterozygous haploids, number of potentially valid call(s), and missing call rate, respectively.
#' Heterozygosity_results contains heterozygosity results for each individual, with six columns as FID, IID, O(HOM), E(HOM), N(NM), and F for Family ID, Within-family ID, Observed number of homozygotes, Expected number of homozygotes, Number of (non-missing, non-monomorphic) autosomal genotype observations and, Method-of-moments F coefficient estimate, respectively.
#' @export
#'
#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' foutput <- "Test_output"
#' imiss = 0.01
#' het = 2
#' x = QCsample(DataDir = DataDir,ResultDir = ResultDir, finput = finput,foutput = foutput, imiss = imiss,het = het)

QCsample <- function(DataDir,
                     ResultDir,
                     finput,
                     foutput = NULL,
                     imiss,
                     het,
                     legend_text_size = 8,
                     legend_title_size = 7,
                     axis_text_size = 5,
                     axis_title_size = 7,
                     title_size = 9,
                     filterSample = TRUE) {


  if (file.exists(paste0(DataDir, "/", finput, ".bed")) &&
      file.exists(paste0(DataDir, "/", finput, ".bim")) &&
      file.exists(paste0(DataDir, "/", finput, ".fam"))) {

    setupPlink(ResultDir)

  } else{
    writeLines(
      "There are no Plink files in DataDir.\nPlease specify correct directory path with input Plink files."
    )
  }

  invisible(sys::exec_wait(
    paste0(ResultDir,"/","./plink"),
    args = c(
      "--bfile",
      paste0(DataDir,"/",finput),
      "--missing",
      "--out",
      paste0(ResultDir,"/","filtered_missing"),
      "--silent"
    ),
    std_out = FALSE,
    std_err = FALSE
  ))

  invisible(sys::exec_wait(
    paste0(ResultDir,"/","./plink"),
    args = c(
      "--bfile",
      paste0(DataDir,"/",finput),
      "--dog",
      "--het",
      "--out",
      paste0(ResultDir,"/","filtered_hetero"),
      "--silent"
    ),
    std_out = FALSE,
    std_err = FALSE
  ))

  miss <- read.table(
    paste0(ResultDir,"/","filtered_missing.imiss"),
    sep = "",
    header = TRUE,
    as.is = TRUE,
    stringsAsFactors = FALSE
  )
  heter <- read.table(
    paste0(ResultDir,"/","filtered_hetero.het"),
    sep = "",
    header = TRUE,
    as.is = TRUE,
    stringsAsFactors = FALSE
  )

  heter$F <- as.numeric(as.character(heter$F))
  miss$F_MISS <- as.numeric(as.character(miss$F_MISS))

  imissfail <- miss[miss$F_MISS > imiss, , drop = FALSE]


  hetfail = heter[heter$F < (mean(heter$F)  - het * stats::sd(heter$F)) |
                    heter$F > (mean(heter$F) + het * stats::sd(heter$F)), , drop = FALSE]


  hetermiss <- merge(miss, heter, by = "IID")

  failed_het_imiss <-
    hetermiss[which(hetermiss$IID %in% union(hetfail$IID, imissfail$IID)), , drop = FALSE]

  write.table(
    failed_het_imiss,
    file = paste0(ResultDir,"/","failed_het_imiss"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    eol = "\r\n",
    sep = " "
  )

  if (filterSample == TRUE){
    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--bfile",
        paste0(DataDir,"/",finput),
        "--exclude",
        paste0(ResultDir,"/","failed_het_imiss"),
        "--make-bed",
        "--out",
        paste0(ResultDir,"/",foutput),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))
  }else{

    print("There is no filtering of the samples since filterSample is FALSE.")

  }

  hetermiss$logF_MISS <- log10(hetermiss$F_MISS)

  hetermiss$type <- "Passed"
  hetermiss$type[hetermiss$IID %in% hetfail$IID] <- "Failed H"
  hetermiss$type[hetermiss$IID %in% imissfail$IID] <- "Failed M"
  hetermiss$type[hetermiss$IID %in%
                   intersect(hetfail$IID, imissfail$IID)] <-
    "Failed H & M"

  ## Standard deviations from mean for labeling
  minus_sd <- mean(hetermiss$F) - 1:5 * (stats::sd(hetermiss$F))
  plus_sd <- mean(hetermiss$F) + 1:5 * (stats::sd(hetermiss$F))


  #############
  colors <- c("#666666", "#1b9e77", "#d95f02", "#7570b3")
  names(colors) <-
    c("Passed", "Failed H", "Failed M", "Failed H & M")

  hetermiss$shape <- "general"
  shape_guide <- FALSE


  hetermiss$type <- factor(hetermiss$type, levels = names(colors))
  hetermiss$shape <- as.factor(hetermiss$shape)

  plot_hetimiss <- ggplot2::ggplot()
  plot_hetimiss <- plot_hetimiss + ggplot2::geom_point(data = hetermiss,
                                                       ggplot2::aes_string(
                                                         x = 'logF_MISS',
                                                         y = 'F',
                                                         color = 'type',
                                                         shape = "shape"
                                                       )) +
    ggplot2::scale_shape_manual(values = c(16, 17), guide = "none") +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(
      x = "Sample's missingness in log10",
      y = "Heterozygosity rate \n(and standard deviation)",
      color = "Sample Status",
      title = "Heterozygosity by missingness across samples"
    ) +
    ggplot2::geom_hline(
      yintercept = c(minus_sd[1:3], plus_sd[1:3]),
      lty = 2,
      col = "azure4"
    ) +
    ggplot2::scale_y_continuous(
      labels = c("-3", "-4", "-5" , "+3", "+4", "+5"),
      breaks = c(minus_sd[3:5], plus_sd[3:5])
    ) +
    ggplot2::scale_x_continuous(
      labels = c(0.0001, 0.001, 0.01, 0.03, 0.05, 0.01, 1),
      breaks = c(-4, -3, -2, log10(0.03), log10(0.05), -1, 0)
    ) +
    ggplot2::geom_hline(
      yintercept = mean(hetermiss$F) - (het * stats::sd(hetermiss$F)),
      col = "#e7298a",
      lty = 2
    ) +
    ggplot2::geom_hline(
      yintercept = mean(hetermiss$F) + (het * stats::sd(hetermiss$F)),
      col = "#e7298a",
      lty = 2
    ) +
    ggplot2::geom_vline(xintercept = log10(imiss),
                        col = "#e7298a",
                        lty = 2)


  plot_hetimiss <- plot_hetimiss +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.text = ggplot2::element_text(size = legend_text_size),
      legend.title = ggplot2::element_text(size = legend_title_size),
      axis.text = ggplot2::element_text(size = axis_text_size),
      title = ggplot2::element_text(size = title_size),
      axis.title = ggplot2::element_text(size = axis_title_size)
    )
  print(plot_hetimiss)
  if (nrow(imissfail) == 0) {
    print("There will be no sample to be filtered for missingness with 'miss' threshold.")
  } else if (nrow(imissfail) != 0) {
    ID_imissfail <- length(unique(imissfail$IID))
    print(paste0(ID_imissfail, " samples are filtered out for missingness."))
  }
  if (nrow(hetfail) == 0) {
    print("There will be no sample to be filtered for heterozygosity with 'het' threshold.")
  } else if (nrow(hetfail) != 0) {
    ID_hetfail <- length(unique(hetfail$IID))
    print(paste0(ID_hetfail, " samples are filtered out for heterozygosity."))
  }


  if (nrow(failed_het_imiss) == 0) {
    print("There will be no sample to filter for missingness and heterozygosity.")
  } else if (nrow(failed_het_imiss) != 0) {
    print(paste0(
      length(unique(failed_het_imiss$IID)),
      " samples are filtered out for missingness and heterozygosity."
    ))
  }

  if (nrow(hetermiss)==0){
    hetermiss <- NULL
  }else {
    hetermiss <- hetermiss
  }

  if (nrow(hetermiss)==0){
    hetermiss1 <- NULL
  }else {
    hetermiss1 <- hetermiss[,1:2]
  }

  if (nrow(imissfail)==0){
    imissfail1 <- NULL
  }else {
    imissfail1 <- imissfail[,1:2]
  }

  if (nrow(hetfail)==0){
    hetfail1 <- NULL
  }else {
    hetfail1 <- hetfail[,1:2]
  }


  print(paste0("Output plink files, ", foutput, " with final samples are in ", ResultDir,"."))

  ftemp <- c("failed_het_imiss","filtered_hetero.log","filtered_missing.lmiss","filtered_missing.log")
  invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))

  if (file.exists(paste0(ResultDir,"/filtered_missing.hh"))){
    ftemp <- c("filtered_missing.hh")
    invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))
  }
  if (file.exists(paste0(ResultDir,"/",foutput,".hh"))){
    ftemp <- c(paste0(foutput,".log"),paste0(foutput,".hh"))
    invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))
  }
  fmi <- read.table(paste0(ResultDir,"/filtered_missing.imiss"))
  fhh <- read.table(paste0(ResultDir,"/filtered_hetero.het"))

  ftemp <- c("filtered_missing.imiss","filtered_hetero.het")
  invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))

  hm <- failed_het_imiss[,2:1]
  colnames(hm)<- c("FID","IID")
  return(list(
    HM = hm,
    Failed_Missingness = imissfail[,1:2],
    Failed_heterozygosity = hetfail[,1:2],
    Missingness_results = fmi,
    Heterozygosity_results = fhh
  ))

}


#' Miami plot
#'
#' @description This function generates Miami plots for GWAS and XWAS
#' @param ResultDir Character string for the folder path where the outputs will be saved.
#' @param FemaleWAS R dataframe of summary statistics of GWAS or XWAS of female samples with four columns, SNP(Variant),CHR(Chromosome number)
#' ,POS(Base pair position) and pvalue(P-value of the test). This can be generated by running FM01comb or
#'  FM02comb model with GXWAS function.
#' @param MaleWAS R dataframe of summary statistics of GWAS or XWAS of male samples with four columns, SNP(Variant),CHR(Chromosome number)
#' ,POS(Base pair position) and pvalue(P-value of the test). This can be generated by running FM01comb or
#'  FM02comb model with GXWAS function.
#' @param snp_pval Numeric value as p-value threshold for annotation. SNPs below this p-value will be annotated on the plot. The default is 1e-08.
#' @param Xchr Boolean value, TRUE or FALSE, specifying whether to generate Miami plot for stratified XWAS or not. The default is TRUE.
#'
#' @return NULL
#' @export
#'
#' @examples
#' # Not Run
# data("GXwasRData")
# #load("~/b1137/BBose/Panscan_ResultDir/fmcomb2_fishertest_pc/FemaleWAS.Rda")
# #load("~/b1137/BBose/Panscan_ResultDir/fmcomb2_fishertest_pc/MaleWAS.Rda")
# #FemaleWAS <- na.omit(FemaleWAS[FemaleWAS$TEST=="ADD",c("SNP","CHR","BP","P")])
# #gc(reset=TRUE)
# #MaleWAS <- na.omit(MaleWAS[MaleWAS$TEST=="ADD",c("SNP","CHR","BP","P")])
# #gc(reset=TRUE)
# FemaleWAS <- na.omit(Ffile[,c("SNP","CHR","BP","P")])
# colnames(Ffile) <- c("SNP","CHR","POS","pvalue")
# MaleWAS <- na.omit(Mfile[,c("SNP","CHR","BP","P")])
# colnames(Mfile) <- c("SNP","CHR","POS","pvalue")
# ResultDir = tempdir()
# snp_pval = 0.05
# GXWASmiami(ResultDir,FemaleWAS,MaleWAS,snp_pval,Xchr = TRUE)

GXWASmiami <- function(ResultDir = tempdir(),FemaleWAS,MaleWAS,snp_pval = 1e-08,Xchr = FALSE){

  print("Generating Miami plots for stratified test.")
  suppressWarnings(invisible(gmirror(top=FemaleWAS, bottom=MaleWAS, tline= snp_pval, bline=snp_pval,
                                             toptitle="GWAS of females", bottomtitle = "GWAS of males",
                                             highlight_p = c(snp_pval,snp_pval), highlighter="green", chrblocks = TRUE,file = paste0(ResultDir,"/","Stratified_GWAS"))))

  print(paste0("Miami plot of stratified GWAS is saved in ", ResultDir))

  if (Xchr == TRUE){
    FemaleWAS <- as.data.frame(FemaleWAS)
    FemaleWAS[FemaleWAS$CHR == "23","CHR"]<-"X"
    MaleWAS <- as.data.frame(MaleWAS)
    MaleWAS[MaleWAS$CHR == "23","CHR"]<-"X"

    # Stratified XWAS plot
    gwas.t2 <- FemaleWAS[FemaleWAS$CHR=="X",]
    gwas.b2 <- MaleWAS[MaleWAS$CHR=="X",]

    rm(FemaleWAS)
    rm(MaleWAS)
    suppressWarnings(invisible(gmirror(top=gwas.t2, bottom=gwas.b2, tline=snp_pval, bline=snp_pval,
                                               toptitle="XWAS of females", bottomtitle = "XWAS of males",
                                               highlight_p = c(snp_pval,snp_pval), highlighter="green", chrblocks = TRUE, file = paste0(ResultDir,"/","Stratified_XWAS"))))
    gc(reset=TRUE)
    print(paste0("Miami plot of stratified XWAS is saved in ", ResultDir))
  }else{
    return()
  }
  return()
}

#' GXwas: Running genome-wide association study (GWAS) and X-chromosome-wide association study (XWAS) models.
#'
#' @author Banabithi Bose
#'
#' @description This function runs GWAS models in autosomes with several alternative XWAS models, such as "FM01", "FM02", "FM01comb" and "FM02comb" can be applied to both binary and quantitative traits while "XCGA" can only be applied to a binary trait. For binary and quantitative features, this function uses logistic and linear regression, allowing for multiple covariates and the interactions with those covariates in a multiple-regression approach. These models are all run using the additive effects of SNPs, and each additional minor allele's influence is represented by the direction of the regression coefficient. This function attempts to identify the multi-collinearity among predictors by displaying NA for the test statistic and a p-value for all terms in the model. The more terms you add, the more likely you are to run into issues. For details about the different XWAS model, please follow the associated publication.
#' @param DataDir Character string for the file path of the input PLINK binary files.
#' @param ResultDir Character string for the folder path where the outputs will be saved.
#' @param finput Character string, specifying the prefix of the input PLINK binary files with both male and female samples. This file needs to be in DataDir. Note: Case/control phenotypes are expected to be encoded as 1=unaffected (control), 2=affected (case); 0 is accepted as an alternate missing value encoding. The missing case/control or quantitative phenotypes are expected to be encoded as 'NA'/'nan' (any capitalization) or -9.
#' @param trait Boolean value, 'binary' or 'quantitative' for the phenotype i.e. the trait.
#' @param perm Boolean value, 'TRUE' or 'FALSE' for using Permutation procedures to generate significance levels empirically.According to the PLINK functionality, the first interval value indicates that this pruning will be carried out every 5 replicates, and the second pruning parameter (0.001) indicates that the rate of pruning will slow down as the number of replicates increases (i.e., pruning will be carried out every 5 + 0.001R replicates, where R is the current number of replicates). For each empirical p-value at each stage of pruning, a 100*(1 - beta / 2T)% confidence interval is calculated, where beta is, in this case, 0.01, and T is the total number of SNPs. We prune any SNP for which the lower confidence limit exceeds alpha or the upper confidence bound falls below alpha using the typical approximation to the binomial.
#' @param mperm Integer value for the threshold of max(T) permutation. The default is 0 i.e., no max(T) permutation. If 1000 permutations are specified, then all 1000 will be performed, for all SNPs. The benefit of doing this is that two sets of empirical significance values can then be calculated -- pointwise estimates of an individual SNPs significance, but also a value that controls for that fact that thousands of other SNPs were tested. This is achieved by comparing each observed test statistic against the maximum of all permuted statistics (i.e. over all SNPs) for each single replicate. In otherwords, the p-value now controls the familywise error rate, as the p-value reflects the chance of seeing a test statistic this large, given you've performed as many tests as you have. Because the permutation schemes preserve the correlational structure between SNPs, this provides a less stringent correction for multiple testing in comparison to the Bonferroni, which assumes all tests are independent. Because it is now the corrected p-value that is of interest, it is typically sufficient to perform a much smaller number of tests -- i.e. it is probably not necessary to demonstrate that something is genome-wide significant beyond 0.05 or 0.01. Note: For non zero mperm value, permutation procedure will take place irrespective of perm value.
#' @param genedrop Boolean value, 'TRUE' or 'FALSE' for whether to perform gene-dropping permutation. The default is FALSE. In its most basic form, this essentially entails flipping the allele that is 50/50 likely to be passed from parent to offspring (1). In this situation, a person must have both parents genotyped for genedropping. In all possible combinations of genedropping, the genotypes of founders and people without two genotyped parents remain constant. This strategy can be applied to generic pedigrees as well, passing founder-derived genes down the generations One can then perform a basic test of association for quantitative traits or samples that contain both affected and unaffected non-founders by treating the pedigree data as if all of the individuals were unrelated. However, permuted datasets made by genedropping will both control for stratification and the non-independence of related individuals (i.e., as these will also be properties of every permuted dataset). Applying the same series of 50:50 flip/no-flip decisions to all SNPs in the same permuted replicate for a certain transmission allows for the maintenance of LD between SNPs. Additionally, linking can be controlled for by giving each sibling in the same nuclear family an identical series of flip or no-flip selections. Following PLINK, both of these functionalities are automatically applied.
#' @param standard_beta Boolean value, 'TRUE' or 'FALSE' in case of quantitative trait for standardizing the trait or phenotype values (mean 0, unit variance), so the resulting coefficients will be standardized. The default is TRUE.
#' @param xmodel Models "FM01", "FM02", "FM01comb" and "FM02comb" can be chosen for both binary and quantitative traits while "XCGA" can only apply to the binary trait. These models take care of the X-chromosomal marker Three female genotypes are coded by 0, 1, and 2 in FM01 and FM02. The two genotypes of males that follow the X-chromosome inactivation (XCI) pattern as random (XCI-R) in the FM01 model are coded by 0 and 1, while the two genotypes that follow the XCI is escaped (XCI-E) in the FM02 model are coded by 0 and 1. To reflect the dose compensation connection between the sexes, FM02 treats men as homozygous females. In the FM01comb and FM01comb methods, associations are tested separately for males and females with the FM01 and FM02 models, respectively, and then the combined p values are computed the Fisher's method, Fisher's method with permutation, or Stouffer's method(1,3-7]. An X-chromosome inactivation (XCI) pattern, or coding technique for X-chromosomal genotypes between sexes, is not required for the XCGA. By simultaneously accounting for four distinct XCI patterns, namely XCI-R, XCI-E, XCI-SN (XCI fully toward normal allele), and XCI-SR (XCI fully toward risk allele), this model may maintain a respectably high power (2). Note that, sex shouldn't be provided as a covariate in the XCGA model.
#' @param noxsex Boolean value, 'TRUE' or 'FALSE' for eliminating sex or adding sex as an additional covariate in the XWAS. Sex will be added as defined in the .fam input PLINK file. This feature is not applicable to the XCGA model.
#' @param covarfile Character string for the full name of the covariate file in .txt format. This file should be placed in DataDir. Note about the covariate file: The first column of this file should be FID, the second column should be IID and the other columns should be covariates. The primary header line should be there starting with “FID”, and “IID” followed by covariate names. If an individual is not present in the covariate file, or if the individual has a missing phenotype value (i.e. -9 by default) for the covariate, then that individual is set to missing (i.e. will be excluded from association analysis). Use the function "DummyCovar()" to generate a new covariate file with categorical variables down-coded as binary dummy variables for the covariate file with categorical variables. For instance, if a variable has K categories, K-1 new dummy variables are constructed, and the original covariate is now estimated with a coefficient for each category.
#' @param covartest Vector value with "NULL","ALL" or covarite name/names to be included in the test. The default is NULL. For instance, the user can choose “AGE” and “SEX” as covartest = c(“AGE”, “SEX”) or all the covariates as covartest = c(“ALL”).
#' @param interaction Boolean value, 'TRUE' or 'FALSE' for including SNP x covariate interaction term/terms from the association analysis. The default is FALSE. If a permutation procedure is chosen then the interaction will be automatically FALSE. For the interaction with the two covariates COV1 and COV2, the model will look like: Y = b0 + b1.ADD + b2.COV1 + b3.COV2 + b4.ADD x COV1 + b5.ADD x COV2 + e. When interaction factors are incorporated into the model, the main effects' significance is not always determined simply; rather, it depends on the arbitrary coding of the variables. To put it another way, you should probably just interpret the p-value for the interaction. Also, The p-values for the covariates do not represent the test for the SNP-phenotype association after controlling for the covariate. That is the first row (ADD). Rather, the covariate term is the test associated with the covariate-phenotype association. These p-values might be extremely significant (e.g. if one covaries for smoking in an analysis of heart disease, etc) but this does not mean that the SNP has a highly significant effect necessarily. Note that, this feature is not valid for XCGA model for XWAS part
#' @param Inphenocov Vector of integer values starting from 1 to extract the terms which user wants from the above model: Y = b0 + b1.ADD + b2.COV1 + b3.COV2 + b4.ADDxCOV1 + b5.ADDxCOV2 + e. The terms will appear in order as (1) for ADD, (2) for COV1, (4) for ADD x COV1, and (5) for ADD x COV2. If the user wants to extract the terms for COV1 and ADD x COV1, she needs to specify it as c(2,4). The default is c(“ALL”). Note that, this feature is not valid for the XCGA model for the XWAS part.
#' @param combtest Character vector specifying method for combining p-values for stratified GWAS with FM01comb and FM02comb XWAS models. Choices are “stouffer.method”, "fisher.method" and "fisher.method.perm". For fisher.method the function for combining p-values uses a statistic, S = -2 ∑^k /log p, which follows a χ^2 distribution with 2k degrees of freedom (3).
#' For fisher.method.perm, using p-values from stratified tests, the summary statistic for combining p-values is S = -2 ∑ /log p. A p-value for this statistic can be derived by randomly generating summary statistics (4). Therefore, a p-value is randomly sampled from each contributing study, and a random statistic is calculated. The fraction of random statistics greater or equal to S then gives the final p-value.
#'
#' For stouffer.method ,the function applies Stouffer’s method (6]) to the p-values assuming that the p-values to be combined are independent. Letting p1, p2, . . . , pk denote the individual (one- or two-sided) p-values of the k hypothesis tests to be combined, the test statistic is then computed with $z = ∑^k_{1}frac{z_{i}}{sqrt(k)}$ where $z_{i}$ = Φ−1 (1 – $p_{i}$) and Φ −1 (·) denotes the inverse of the cumulative distribution function of a standard normal distribution. Under the joint null hypothesis, the test statistic follows a standard normal distribution which is used to compute the combined p-value. This functionality is taken from the R package poolr (7).
#' Note that only p-values between 0 and 1 are allowed to be passed to these methods.
#' @param MF.zero.sub Small numeric value for substituting p-values of 0 in in stratified GWAS with FM01comb and FM02comb XWAS models. The default is 0.00001. As log(0) results in Inf this replaces p-values of 0 by default with a small float.
#' @param B Integer value specifying the number of permutation in case of using fisher.method.perm method in stratified GWAS with FM01comb and FM02comb XWAS models. The default is 10000.
#' @param MF.na.rm Boolean value, 'TRUE' or 'FALSE' for removing p-values of NA in stratified GWAS with FM01comb and FM02comb XWAS in case of using Fisher’s and Stouffer’s methods. The default is FALSE.
#' @param MF.p.corr Character vector specifying method for correcting the summary p-values for FMfcomb and FMscomb models. Choices are "bonferroni", "BH" and "none" for Bonferroni,  Benjamini-Hochberg and none, respectively. The default is "none".
#' @param MF.mc.cores Number of cores used for fisher.method.perm in stratified GWAS with FM01comb and FM02comb XWAS models.
#' @param plot.jpeg Boolean value, 'TRUE' or 'FALSE' for saving the plots in .jpeg file. The default is TRUE.
#' @param plotname A character string specifying the prefix of the file for plots. This file will be saved in DataDir. The default is "GXwas.plot".
#' @param snp_pval Numeric value as p-value threshold for annotation. SNPs below this p-value will be annotated on the plot. The default is 1e-08.
#' @param annotateTopSnp Boolean value, 'TRUE' or 'FALSE. If TRUE, it only annotates the top hit on each chromosome that is below the snp_pval threshold. The default is FALSE.
#' @param suggestiveline The deafault is 5 (for p-value 1e-05).
#' @param genomewideline The default is 7.3 (for p-value 5e-08).
#' @param ncores Integer value, specifying the number of cores for parallel processing. The default is 0 (no parallel computation).
#'
#' @return A dataframe with GWAS (with XWAS for X-chromosomal variants) along with Manhattan and Q-Q plots.
#' In the case of the stratified test, the return is a list containing three dataframes, namely, FWAS, MWAS, and MFWAS with association results in only female, only male, and both cohorts, respectively. This will be accompanied by Miami and Q-Q plots. The individual manhattan and Q-Q-plots for stratified tests prefixed with xmodel type will be in the DataDir.
#' @export
#'
#' @references (1) Purcell et. al.,(2007). PLINK: A Tool Set for Whole-Genome Association and Population-Based Linkage Analyses. The American Journal of Human Genetics, 81(3), 559-575. https://doi.org/10.1086/519795.
#' (2) Su Y, Hu J, Yin P, Jiang H, Chen S, Dai M, Chen Z, Wang P. XCMAX4: A Robust X Chromosomal Genetic Association Test Accounting for Covariates. Genes (Basel). 2022 May 9;13(5):847. doi: 10.3390/genes13050847. PMID: 35627231; PMCID: PMC9141238.
#' (3) Fisher, R.A. (1925). Statistical Methods for Research Workers. Oliver and Boyd (Edinburgh);
#' (4) Rhodes, D. R., (2002). Meta-analysis of microarrays: interstudy validation of gene expression profiles reveals pathway dysregulation in prostate cancer. Cancer research, 62(15), 4427-33;
#' (5) Moreau, Y.et al. (2003). Comparison and meta-analysis of microarray data: from the bench to the computer desk. Trends in Genetics, 19(10), 570-577.
#' (6) Stouffer, S. A., Suchman, E. A., DeVinney, L. C., Star, S. A., & Williams, R. M., Jr. (1949). The American Soldier: Adjustment During Army Life (Vol. 1). Princeton, NJ: Princeton University Press.
#' (7) Cinar, O. & Viechtbauer, W. (2022). The poolr package for combining independent and dependent p values. Journal of Statistical Software, 101(1), 1–42. https://doi.org/10.18637/jss.v101.i01
#'
#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' #ResultGXwas <- GXwas(DataDir = DataDir,finput = finput, xmodel = "FM02comb", covarfile = NULL, covartest = NULL, combtest = "fisher.method", snp_pval = 0.001, plot.jpeg = FALSE, suggestiveline = 3, genomewideline = 5.69897)
#' ResultGXwas <- GXwas(DataDir = DataDir,ResultDir = ResultDir, finput = finput, xmodel = "FM02comb", covarfile = NULL, covartest = NULL, combtest = "fisher.method", snp_pval = 0.001, plot.jpeg = FALSE, suggestiveline = 3, genomewideline = 5.69897, ncores = 0)
#' #ResultGXwas$CombinedWAS[1:5,]
#' #ResultGXwas$MaleWAS[1:5,]
#' #ResultGXwas$FemaleWAS[1:5,]

# GXwas <- function(DataDir,ResultDir, finput, trait = c("binary","quantitative"), standard_beta = "TRUE", perm = "FALSE",
#                   mperm = 0, genedrop = "FALSE", xmodel = c("FM01","FM02","FM01comb","FM02comb","XCGA"), noxsex = "FALSE",
#                   covarfile = NULL, interaction = "FALSE", covartest = c("ALL"), Inphenocov = c("ALL"), combtest = c("fisher.method", "fisher.method.perm", "stouffer.method"),
#                   MF.zero.sub = 0.00001, B = 10000,  MF.mc.cores = NULL, MF.na.rm = FALSE,
#                   MF.p.corr = "none", plot.jpeg = FALSE, plotname = "GXwas.plot", snp_pval = 1e-08,
#                   annotateTopSnp = FALSE, suggestiveline = 5, genomewideline = 7.3, ncores = 0){
#
#   if (file.exists(paste0(DataDir, "/", finput, ".bed")) &&
#       file.exists(paste0(DataDir, "/", finput, ".bim")) &&
#       file.exists(paste0(DataDir, "/", finput, ".fam"))) {
#
#     setupPlink(ResultDir)
#
#   } else{
#     writeLines(
#       "There are no Plink files in DataDir.\nPlease specify correct directory path with input Plink files."
#     )
#   }
#
#
#   if (perm[1] == TRUE && interaction[1] == TRUE){
#
#     interaction <- FALSE
#
#   }else{
#
#     interaction <- interaction
#   }
#
#
#
#   if (mperm != 0){
#
#     perm <- TRUE
#
#   }else if (mperm == 0){
#
#     perm <- perm
#
#   }
#
#
#   if (xmodel[1] == "FM01"){
#
#     print("Running FM01 model")
#
#
#     x<- suppressWarnings(FMmain(DataDir = DataDir, ResultDir = ResultDir, finput = finput, trait = trait, standard_beta = standard_beta, perm = perm, mperm = mperm, genedrop = genedrop, xmodel = xmodel,
#                                 noxsex = noxsex, covarfile = covarfile, interaction = interaction, covartest = covartest, Inphenocov = Inphenocov, plot.jpeg = plot.jpeg, plotname = plotname, snp_pval = snp_pval, annotateTopSnp = annotateTopSnp, suggestiveline = suggestiveline, genomewideline = genomewideline, ncores = ncores))
#
#   }else if (xmodel[1] == "FM02"){
#
#     print("Running FM02 model")
#
#
#     x<- suppressWarnings(FMmain(DataDir = DataDir, ResultDir = ResultDir, finput = finput, trait = trait, standard_beta = standard_beta, perm = perm, mperm = mperm, genedrop = genedrop, xmodel = xmodel,
#                                 noxsex = noxsex, covarfile = covarfile, interaction = interaction, covartest = covartest, Inphenocov = Inphenocov, plot.jpeg = plot.jpeg, plotname = plotname, snp_pval = snp_pval, annotateTopSnp = annotateTopSnp, suggestiveline = suggestiveline, genomewideline = genomewideline))
#
#
#   }else if (xmodel[1] == "FM01comb" | xmodel[1] == "FM02comb"){
#
#     print("Running FMfcomb model")
#
#     ## Making male and female files in ResultDir
#     MFsplitPlink(DataDir = DataDir, ResultDir = ResultDir, finput = finput, foutput = "finput.female", sex = "females")
#     gc(reset = TRUE)
#     MFsplitPlink(DataDir = DataDir, ResultDir = ResultDir, finput = finput, foutput = "finput.male", sex = "males")
#     gc(reset = TRUE)
#
#     x<- suppressWarnings(FMcomb(DataDir = DataDir, ResultDir = ResultDir, trait = trait, standard_beta = standard_beta, perm = perm, mperm = mperm, genedrop = genedrop, xmodel = xmodel,
#                                 noxsex = noxsex, covarfile = covarfile, interaction = interaction, covartest = covartest, Inphenocov = Inphenocov
#                                 , plot.jpeg = plot.jpeg, plotname = plotname, snp_pval = snp_pval, annotateTopSnp = annotateTopSnp,
#                                 combtest = combtest, B = B, MF.p.corr = MF.p.corr, MF.zero.sub = MF.zero.sub, MF.na.rm = MF.na.rm, MF.mc.cores = MF.mc.cores, suggestiveline = suggestiveline, genomewideline = genomewideline, ncores = ncores))
#
#     ftemp <- list.files(paste0(ResultDir,"/"),pattern = "finput")
#     invisible(file.remove(paste0(ResultDir,"/",ftemp)))
#
#   }else if (xmodel[1] == "XCGA"){
#
#     if (trait[1] == "quantitative"){
#       return(print("For XCGA model, trait needs to be quantitative. Please correct the input file."))
#     }else{
#       print("Running XCGA model")}
#
#     x <- suppressWarnings(XCMAFun(DataDir = DataDir, ResultDir = ResultDir, finput = finput, trait = trait, standard_beta = standard_beta, perm = perm, mperm = mperm, genedrop = genedrop,
#                                   noxsex = noxsex, covarfile = covarfile, interaction = interaction, covartest = covartest, Inphenocov = Inphenocov, plot.jpeg = plot.jpeg, plotname = plotname, snp_pval = snp_pval, annotateTopSnp = annotateTopSnp, suggestiveline = suggestiveline, genomewideline = genomewideline, ncores = ncores))
#
#   }
#
#   return(x)
# }

GXwas <- function(DataDir,ResultDir, finput, trait = c("binary","quantitative"), standard_beta = "TRUE", perm = "FALSE",
                  mperm = 0, genedrop = "FALSE", xmodel = c("FM01","FM02","FM01comb","FM02comb","XCGA"), noxsex = "FALSE",
                  covarfile = NULL, interaction = "FALSE", covartest = c("ALL"), Inphenocov = c("ALL"), combtest = c("fisher.method", "fisher.method.perm", "stouffer.method"),
                  MF.zero.sub = 0.00001, B = 10000,  MF.mc.cores = NULL, MF.na.rm = FALSE,
                  MF.p.corr = "none", plot.jpeg = FALSE, plotname = "GXwas.plot", snp_pval = 1e-08,
                  annotateTopSnp = FALSE, suggestiveline = 5, genomewideline = 7.3, ncores = 0){

  if (file.exists(paste0(DataDir, "/", finput, ".bed")) &&
      file.exists(paste0(DataDir, "/", finput, ".bim")) &&
      file.exists(paste0(DataDir, "/", finput, ".fam"))) {

    setupPlink(ResultDir)

  } else{
    writeLines(
      "There are no Plink files in DataDir.\nPlease specify correct directory path with input Plink files."
    )
  }


  if (perm[1] == TRUE && interaction[1] == TRUE){

    interaction <- FALSE

  }else{

    interaction <- interaction
  }



  if (mperm != 0){

    perm <- TRUE

  }else if (mperm == 0){

    perm <- perm

  }


  if (xmodel[1] == "FM01"){

    print("Running FM01 model")


    x<- suppressWarnings(FMmain(DataDir = DataDir, ResultDir = ResultDir, finput = finput, trait = trait, standard_beta = standard_beta, perm = perm, mperm = mperm, genedrop = genedrop, xmodel = xmodel,
                                noxsex = noxsex, covarfile = covarfile, interaction = interaction, covartest = covartest, Inphenocov = Inphenocov, plot.jpeg = plot.jpeg, plotname = plotname, snp_pval = snp_pval, annotateTopSnp = annotateTopSnp, suggestiveline = suggestiveline, genomewideline = genomewideline, ncores = ncores))

  }else if (xmodel[1] == "FM02"){

    print("Running FM02 model")


    x<- suppressWarnings(FMmain(DataDir = DataDir, ResultDir = ResultDir, finput = finput, trait = trait, standard_beta = standard_beta, perm = perm, mperm = mperm, genedrop = genedrop, xmodel = xmodel,
                                noxsex = noxsex, covarfile = covarfile, interaction = interaction, covartest = covartest, Inphenocov = Inphenocov, plot.jpeg = plot.jpeg, plotname = plotname, snp_pval = snp_pval, annotateTopSnp = annotateTopSnp, suggestiveline = suggestiveline, genomewideline = genomewideline))


  }else if (xmodel[1] == "FM01comb" | xmodel[1] == "FM02comb"){

    print("Running FMfcomb model")

    ## Making male and female files in ResultDir
    MFsplitPlink(DataDir = DataDir, ResultDir = ResultDir, finput = finput, foutput = "finput.female", sex = "females")
    gc(reset = TRUE)
    MFsplitPlink(DataDir = DataDir, ResultDir = ResultDir, finput = finput, foutput = "finput.male", sex = "males")
    gc(reset = TRUE)

    x<- suppressWarnings(FMcomb(DataDir = DataDir, ResultDir = ResultDir, trait = trait, standard_beta = standard_beta, perm = perm, mperm = mperm, genedrop = genedrop, xmodel = xmodel,
                                noxsex = noxsex, covarfile = covarfile, interaction = interaction, covartest = covartest, Inphenocov = Inphenocov
                                , plot.jpeg = plot.jpeg, plotname = plotname, snp_pval = snp_pval, annotateTopSnp = annotateTopSnp,
                                combtest = combtest, B = B, MF.p.corr = MF.p.corr, MF.zero.sub = MF.zero.sub, MF.na.rm = MF.na.rm, MF.mc.cores = MF.mc.cores, suggestiveline = suggestiveline, genomewideline = genomewideline, ncores = ncores))

    ftemp <- list.files(paste0(ResultDir,"/"),pattern = "finput")
    invisible(file.remove(paste0(ResultDir,"/",ftemp)))

  }else if (xmodel[1] == "XCGA"){

    if (trait[1] == "quantitative"){
      return(print("For XCGA model, trait needs to be quantitative. Please correct the input file."))
    }else{
      print("Running XCGA model")}

    x <- suppressWarnings(XCMAFun(DataDir = DataDir, ResultDir = ResultDir, finput = finput, trait = trait, standard_beta = standard_beta, perm = perm, mperm = mperm, genedrop = genedrop,
                                  noxsex = noxsex, covarfile = covarfile, interaction = interaction, covartest = covartest, Inphenocov = Inphenocov, plot.jpeg = plot.jpeg, plotname = plotname, snp_pval = snp_pval, annotateTopSnp = annotateTopSnp, suggestiveline = suggestiveline, genomewideline = genomewideline, ncores = ncores))

  }

  return(x)
}


#' MetaGXwas: Combining summary-level results from two or more GWA studies into a single estimate.
#'
#' @description  This function combine K sets of GWAS association statistics on same (or at least similar) phenotype. This function employs PLINK's (1) inverse variance-based analysis to run a number of models, including a) Fixed-effect model and b) Random-effect model, assuming there may be variation between the genuine underlying effects, i.e., effect size beta. 'This function also calculates weighted Z-score-based p-values after METAL (3). For more information about the algorithms, please see the associated paper.
#'
#' @param DataDir A character string for the file path of the input files needed for ‘SummData’ and ‘SNPfile’ arguments.
#' @param SummData Vector value containing the name/names of the .Rda file/files with GWAS summary statistics, with ‘SNP’ (i.e., SNP idenitifier), ‘BETA’ (i.e., effect-size or logarithm of odds ratio), ‘SE’ (i.e., standard error of BETA), ‘P’ (i.e., p-values) and 'NMISS' (i.e., effective sample size) in mandatory column headers. These files needed to be in DataDir. If the numbers of cases and controls are unequal, effective sample size should be 4 / (1/<# of cases> + 1/<# of controls>). A smaller "effective" sample size may be used for samples that include related individuals, however simulations indicate that small changes in the effective sample size have relatively little effect on the final p-value (3). Columns, such as, 'CHR'(Chromosome code), 'BP'(Basepair position), 'A1' (First allele code), 'A2' (Second allele code) columns are optional. If these are present, setting ‘useSNPposition’ to FALSE, causes 'CHR', ‘BP’ and ‘A1’ to be ignored and setting ‘UseA1’ to be ‘FALSE’ causes A1 to be ignored. If, both these arguments are true, this function takes care of A1/A2 allele flips properly. Otherwise, A1 mismatches are thrown out. Values of CHR/BP are allowed to vary.
#' @param ResultDir A character string for the file path where all output files will be stored. The default is tempdir().
#' @param SNPfile Character string specifying the name of the plain-text file with a column of SNP names. These could be LD clumped SNPs or any other list of chosen SNPs for Meta analysis. This file needs to be in DataDir.
#' @param useSNPposition Boolean value, 'TRUE' or 'FALSE' for using 'CHR', 'BP', and ‘A1’ or not. The default is FALSE.
#' @param UseA1 Boolean value, 'TRUE' or 'FALSE' for ‘A1’ to be used or not. The default is FALSE.
#' @param GCse  Boolean value, 'TRUE' or 'FALSE' for applying study specific genomic control to adjust each study for potential population structure for all the SNPs. The default is TRUE. If users would want to apply genomic control separately for directly genotyped and imputed SNPs prior using the function, set this parameter as FALSE.
#' @param plotname Character string, spycifying the plot name of the file containg forest plots for the SNPs. The default is “Meta_Analysis.plot”.
#' @param pval_filter Character value as "R","F" or "W", specifying whether p-value threshold should be chosen based on “Random”, “Fixed” or “Weighted” effect model for the SNPs to be included in the forest plots.
#' @param top_snp_pval Numeric value, specifying the threshold to be used to filter the SNPs for the forest plots. The default is 1e-08.
#' @param max_top_snps Integer value, specifying the maximum number of top SNPs (SNPs with the lowest p-values) to be ploted in the forest plot file. The default is 6.
#' @param chosen_snps_file Character string specifing the name of the plain-text file with a column of SNP names for the forest plots. The default is NULL.
#' @param byCHR Boolean value, 'TRUE' or 'FALSE', specifying whether the meta analysis will be performed chromosome wise or not. The default is FALSE.
#' @param pval_threshold_manplot Numeric value, specifying the p-value threshold for plotting Manhattan plots.
#'
#' @returns A list object containing five dataframes. The first three dataframes, such as Mfixed, Mrandom and Mweighted contain results
#' for fixed effect, random effect and weighted model. Each of these dataframes can have maximum 12 columns, such as 'CHR'(Chromosome code),'BP'(Basepair position),'SNP'(SNP identifier), 'A1' (First allele code),
#' 'A2' (Second allele code), Q (p-value for Cochrane's Q statistic), I (I^2 heterogeneity index (0-100)),
#' P (P-value from mata analysis), ES (Effect-size estimate from mata analysis), SE (Standard Error from mata analysis),
#' CI_L (Lower limit of confidence interval) and CI_U (Uper limit of confidence interval).
#'
#' The fourth dataframe contains the same columns "CHR","BP","SNP","A1","A2","Q","I", with column N' ( Number of valid studies for this SNP), P (Fixed-effects meta-analysis p-value),
#' and other columns as Fx... (Study x (0-based input file indices) effect estimate, Examples: F0, F1 etc.).
#'
#' The fifth dataframe, ProblemSNP has three columns, such as, 'File' (file name of input data), 'SNP' (Problematic SNPs that are thrown)
#' and 'Problem' (Problem code). Problem codes are, BAD_CHR (Invalid chromosome code), BAD_BP (Invalid base-position code),
#' BAD_ES (Invalid effect-size), BAD_SE (Invalid standard error), MISSING_A1 (Missing allele 1 label),
#' MISSING_A2 (Missing allele 2 label), ALLELE_MISMATCH (Mismatching allele codes across files).
#'
#' A .pdf file comprising the forest plots of the SNPs is produced in the ResultDir with Plotname as prefix.
#' If useSNPposition is set TRUE, a .jpeg file with Manhattan Plot and Q-Q plot will be in the ResultDir with Plotname as prefix.
#'
#' @references (1) Purcell et. al.,(2007). PLINK: A Tool Set for Whole-Genome Association and Population-Based Linkage Analyses. The American Journal of Human Genetics, 81(3), 559-575. https://doi.org/10.1086/519795.
#' (2) Mägi, R., Morris, A.P. GWAMA: software for genome-wide association meta-analysis. BMC Bioinformatics 11, 288 (2010). https://doi.org/10.1186/1471-2105-11-288
#' (3) Willer CJ, Li Y, Abecasis GR. METAL: fast and efficient meta-analysis of genomewide association scans. Bioinformatics. 2010 Sep 1;26(17):2190-1. doi: 10.1093/bioinformatics/btq340. Epub 2010 Jul 8. PMID: 20616382; PMCID: PMC2922887.

#' @importFrom qqman manhattan qq
#' @importFrom graphics par
#' @importFrom grDevices jpeg pdf
#'
#' @export
#'
#' @examples
#' # Not Run
# DataDir = system.file("extdata", package = "GXwasR")
# ResultDir = tempdir()
# data("GXwasRData")
# #load(paste0("/projects/b1137/BBose/GXwasR/data/Summary_Stat_Ex1.Rda"))
# #load(paste0("/projects/b1137/BBose/GXwasR/data/Summary_Stat_Ex2.Rda"))
# #SummData = c("SNPtest1","SNPtest2")
# SummData <- list(Summary_Stat_Ex1,Summary_Stat_Ex2)
# #SNPfile = "UniqueLoci"
# useSNPposition = FALSE
# UseA1 = TRUE
# GCse = TRUE
# byCHR = FALSE
# pval_filter = "R"
# top_snp_pval = 1e-08
# max_top_snps = 10
# chosen_snps_file = NULL
# #chosen_snps_file = "MainSNP"
# pval_threshold_manplot = 1e-05
# x <- MetaGWAS(DataDir = DataDir, SummData = SummData,ResultDir=ResultDir, SNPfile = NULL, useSNPposition = TRUE, UseA1 = UseA1,GCse = GCse, plotname = "Meta_Analysis.plot", pval_filter, top_snp_pval, max_top_snps, chosen_snps_file = NULL, byCHR, pval_threshold_manplot)
# # Dataframe with the fixed effect result
# RandomResult <- x$Resultrandom
# # Dataframe with the random effect result
# RandomResult <- x$Resultrandom
# # Dataframe with the problem SNPs
# MetaProblem <- x$ProblemSNP
# # MetaData
# MetaData <- x$Metadata

MetaGWAS <- function(DataDir, SummData = c(""),ResultDir = tempdir(), SNPfile = NULL,
                     useSNPposition = TRUE,
                     UseA1 = FALSE,GCse = TRUE,
                     plotname = "Meta_Analysis.plot", pval_filter = "R",
                     top_snp_pval = 1e-08, max_top_snps = 6, chosen_snps_file = NULL,
                     byCHR = FALSE, pval_threshold_manplot = 1e-05){


  if (useSNPposition == TRUE){
    nomap = NULL
  }else{
    nomap = "no-map"
  }

  if (UseA1 == TRUE){
    UseA1v = NULL
  }else{
    UseA1v = "no-allele"
  }


  if (is.null(SNPfile)|byCHR == TRUE){
    extract = NULL
    SNPfilev = NULL
  }else{
    extract = "--extract"
    SNPfilev = SNPfile
  }

  for (i in 1:length(SummData)) {
    print(i)
    write.table(SummData[[i]], paste0(ResultDir,"/","SNPdata_",i), row.names=FALSE, col.names = TRUE, quote = FALSE)
  }

  SummData <- list.files(paste0(ResultDir,"/"),pattern = "SNPdata_")

  # Apply genomic control
  if (GCse == TRUE){
    getGCse <- function(SummData){
      s1 <- read.table(paste0(ResultDir,"/",SummData),header = TRUE)

      # From p-values, calculate chi-squared statistic
      chisq <- qchisq(1-s1$P,1)
      lamdaGC <- median(chisq)/qchisq(0.5,1)
      s1$SE <- s1$SE * sqrt(lamdaGC)
      write.table(s1,paste0(ResultDir,"/",SummData),col.names = TRUE, row.names = FALSE, quote = FALSE)
      return()
    }

    #lapply(as.list(SummData),getGCse)
    lapply(SummData,getGCse)

  }else{
    print("No study-specific genomic control was applied.")
    file.copy(from=paste0(DataDir,"/",SummData),to=paste0(ResultDir,"/",SummData),overwrite = TRUE,copy.mode = TRUE)
  }

  if (is.null(SNPfile)){
    SNPfilev = NULL
  }else{
    SNPfilev = paste0(DataDir,"/",SNPfile)
  }

  metaFun <- function(DataDir,ResultDir,SummData,CHR,chromosome,nomap,UseA1v,extract,SNPfilev){
    setupPlink(ResultDir)
    chromosomev = chromosome

    print(paste0("Processing chromosome ",chromosomev))

    invisible(sys::exec_wait(
      paste0(ResultDir, "/","./plink"),
      args = c(
        "--meta-analysis",paste0(ResultDir,"/",SummData),
        "+",
        "qt",
        "report-all",
        nomap,
        UseA1v,
        "report-all",
        "weighted-z",
        "study",
        CHR,chromosomev,
        extract,SNPfilev,
        "--out",
        paste0(ResultDir, "/","MetaResult"),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))

    metaresult <- read.table( paste0(ResultDir, "/","MetaResult.meta"),header = TRUE)
    return(metaresult)
  }

  if (byCHR == FALSE){
    MR <- metaFun(DataDir = DataDir,ResultDir = ResultDir,SummData = SummData,
                  CHR = NULL,chromosome = NULL,nomap = nomap,
                  UseA1v = UseA1v,extract = extract,SNPfilev = SNPfilev)
  }else{
    chromosome = 1:23
    MR <- data.table::rbindlist(lapply(chromosome,metaFun,DataDir = DataDir,ResultDir = ResultDir,SummData = SummData,
                                       CHR = "--chr",nomap = NULL,
                                       UseA1v = NULL,extract = extract, SNPfilev = SNPfilev))
  }


  # Standard error
  MR$SEfixed <- abs(MR$BETA/qnorm(MR$P/2))
  MR$SErandom <- abs(MR$BETA.R./qnorm(MR$P.R./2))
  MR$SEweighted <- abs(MR$WEIGHTED_Z/qnorm(MR$P.WZ./2))

  # 95% Confidence interval
  MR$CIfixedLL <- MR$BETA-1.96*MR$SEfixed
  MR$CIfixedUL <- MR$BETA+1.96*MR$SEfixed

  MR$CIrandomLL <- MR$BETA.R.-1.96*MR$SErandom
  MR$CIrandomUL <- MR$BETA.R.+1.96*MR$SErandom

  MR$CIweightedLL <- MR$WEIGHTED_Z-1.96*MR$SEweighted
  MR$CIweightedUL <- MR$WEIGHTED_Z+1.96*MR$SEweighted

  # Getting effect size and confidence interval for studies
  # Only for the SNPs in MR
  MRsnps <- MR[,"SNP",drop = FALSE]
  getStudyCI <- function(SummData,MRsnps){
    s1 <- read.table(paste0(ResultDir,"/",SummData),header = TRUE)
    s11 <- merge(MRsnps,s1, by = "SNP")
    S1 <- s11[,c("SNP","BETA","L95","U95")]
    return(S1)
  }

  Sbeta <- data.table::rbindlist(lapply(as.list(SummData),getStudyCI,MRsnps = MRsnps))


  # Selecting top SNPs for forest plot ()
  if (pval_filter == "R"){
    MR1 <- MR[MR$P<= top_snp_pval,,drop = FALSE]
    MR2 <- MR1[order(MR1$P),,drop = FALSE]
  }else if (pval_filter == "F"){
    MR1 <- MR[MR$P.F.<= top_snp_pval, ,drop = FALSE]
    MR2 <- MR1[order(MR1$P.F.),,drop = FALSE]
  }else if (pval_filter == "W"){
    MR1 <- MR[MR$P.WZ.<= top_snp_pval, ,drop = FALSE]
    MR2 <- MR1[order(MR1$P.WZ.),,drop = FALSE]
  }
  # If chosen_snps_file is provided, then MR2
  ##ADDING NEW 5/16/2023
  if (!is.null(chosen_snps_file)){
    chosenS <- read.table(paste0(DataDir,"/",chosen_snps_file),header = TRUE)
    colnames(chosenS) <- "SNP"
    MR2 <- merge(chosenS,MR, by = "SNP")
  }else{
    MR2 <- MR2
  }
  # Forest plots interactive
  #graphics::par(mfrow = c(2, 2),mar = c(1, 2, 5, 1),oma=c(3,3,3,3))
  # graphics::par(mfrow = c(2, 2))
  #adjust plot margins
  #par(mar = c(1, 1, 1, 1))
  if (length(unique(MR2$SNP))< 10){

    i = 1:length(unique(MR2$SNP))

  }else{
    print("Maximun 10 Forest plots of 10 chosen SNPs will be drawn in the plot window. For all other Forest plots, please check ResultDir.")
    i = 1:10
  }

  invisible(suppressWarnings(lapply(i,topForestplot,MR2=MR2,Sbeta=Sbeta)))

  options(warn=-1)
  invisible(suppressWarnings(graphics::mtext(text="Studies (black) and tests (blue)",side=4,line= 0, outer=TRUE)))
  invisible(suppressWarnings(graphics::mtext(text="Effect size (beta, 95% CI)",side=1,line= 1,outer=TRUE)))

  if (is.null(chosen_snps_file)){
    invisible(suppressWarnings(graphics::mtext(text="Forest plots of a few chosen SNPs",side=3,line=.5,outer=TRUE)))
  }else{
    invisible(suppressWarnings(graphics::mtext(text="Forest plots of the top SNPs",side=3,line= .5,outer=TRUE)))
  }
  #dev.off()
  if (useSNPposition == TRUE){
    Mfixed <- MR[,c("CHR","BP","SNP","A1","A2","Q","I","P","BETA","SEfixed","CIfixedLL","CIfixedUL")]
    colnames(Mfixed) <- c("CHR","BP","SNP","A1","A2","Q","I","P","ES","SE","CI_L","CI_U")
    Mrandom <- MR[,c("CHR","BP","SNP","A1","A2","Q","I","P.R.","BETA.R.","SErandom","CIrandomLL","CIrandomUL")]
    colnames(Mrandom) <- c("CHR","BP","SNP","A1","A2","Q","I","P","ES","SE","CI_L","CI_U")
    Mweighted <- MR[,c("CHR","BP","SNP","A1","A2","Q","I","P.WZ.","WEIGHTED_Z","SEweighted","CIweightedLL","CIweightedUL")]
    colnames(Mweighted) <- c("CHR","BP","SNP","A1","A2","Q","I","P","ES","SE","CI_L","CI_U")
    Msummdata <- MR[ , !names(MR) %in%
                       c("P","BETA","SEfixed","CIfixedLL","CIfixedUL","P.R.","BETA.R.","SErandom","CIrandomLL","CIrandomUL","P.WZ.","WEIGHTED_Z","SEweighted","CIweightedLL","CIweightedUL")]
  }else{
    Mfixed <- MR[,c("SNP","Q","I","P","BETA","SEfixed","CIfixedLL","CIfixedUL")]
    colnames(Mfixed) <- c("SNP","Q","I","P","ES","SE","CI_L","CI_U")
    Mrandom <- MR[,c("SNP","Q","I","P.R.","BETA.R.","SErandom","CIrandomLL","CIrandomUL")]
    colnames(Mrandom) <- c("SNP","Q","I","P","ES","SE","CI_L","CI_U")
    Mweighted <- MR[,c("SNP","Q","I","P.WZ.","WEIGHTED_Z","SEweighted","CIweightedLL","CIweightedUL")]
    colnames(Mweighted) <- c("SNP","Q","I","P","ES","SE","CI_L","CI_U")
    Msummdata <- MR[ , !names(MR) %in%
                       c("P","BETA","SEfixed","CIfixedLL","CIfixedUL","P.R.","BETA.R.","SErandom","CIrandomLL","CIrandomUL","P.WZ.","WEIGHTED_Z","SEweighted","CIweightedLL","CIweightedUL")]

  }

  if (useSNPposition == TRUE){

    options(bitmapType='cairo')
    grDevices::jpeg(paste0(ResultDir,"/",plotname,".jpeg"),  width = 20,
                    height = 10,
                    units = 'in',
                    res = 300)
    graphics::par(mfrow = c(3, 2))
    mR <- na.omit(Mfixed[,c("SNP","CHR","BP","P")])
    # From p-values, calculate chi-squared statistic
    chisq <- qchisq(1-mR$P,1)
    lamdaGC <- median(chisq)/qchisq(0.5,1)

    invisible(suppressWarnings(qqman::manhattan(mR, ylim = c(0,10),annotatePval = pval_threshold_manplot, annotateTop = FALSE, main = "Manhattan plot of fixed effect meta GWAS")))
    invisible(suppressWarnings(qqman::qq(mR$P, main = paste0(("Q-Q plot of fixed effect meta GWAS p-values with GIF = "), lamdaGC))))

    mR <- na.omit(Mrandom[,c("SNP","CHR","BP","P")])
    # From p-values, calculate chi-squared statistic
    chisq <- qchisq(1-mR$P,1)
    lamdaGC <- median(chisq)/qchisq(0.5,1)
    invisible(suppressWarnings(qqman::manhattan(mR, ylim = c(0,10),annotatePval = pval_threshold_manplot, annotateTop = FALSE, main = "Manhattan plot of random effect meta GWAS")))
    invisible(suppressWarnings(qqman::qq(mR$P, main = paste0(("Q-Q plot of random effect meta GWAS p-values with GIF = "), lamdaGC))))

    mR <- na.omit(Mweighted[,c("SNP","CHR","BP","P")])
    # From p-values, calculate chi-squared statistic
    chisq <- qchisq(1-mR$P,1)
    lamdaGC <- median(chisq)/qchisq(0.5,1)
    invisible(suppressWarnings(qqman::manhattan(mR, ylim = c(0,10),annotatePval = pval_threshold_manplot, annotateTop = TRUE, main = "Manhattan plot of weighted Z-score meta GWAS")))
    invisible(suppressWarnings(qqman::qq(mR$P, main = paste0(("Q-Q plot of weighted Z-score meta GWAS p-values with GIF = "), lamdaGC))))
    dev.off()
  }

  ## All forest plots in .pdf
  grDevices::pdf(paste0(ResultDir,"/",plotname,".pdf"),width=10,height = 5)
  #graphics::par(mar = c(3, 2, 2, 3),oma=c(3,0,3,3))
  i = 1:length(MR2$SNP)
  invisible(suppressWarnings(lapply(i,allForestplot,MR2=MR2,Sbeta = Sbeta)))
  dev.off()

  print(paste0(plotname," files containing Manhattan Plot, Q-Q plot and the forest plots of the SNPs are produced in the directory ", ResultDir,"."))


  if (file.exists(paste0(ResultDir, "/MetaResult.prob"))){
    MP <- read.table(paste0(ResultDir,"/","MetaResult.prob"))
    colnames(MP) <- c("File","SNP","Problem")
  }else{
    MP <- data.frame("None","None","None")
    colnames(MP) <- c("File","SNP","Problem")
  }
  ftemp <- list.files(paste0(ResultDir,"/"),pattern = "SNPdata")
  invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))
  ftemp <- list.files(paste0(ResultDir,"/"),pattern = "MetaResult")
  invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))

  return(list(Resultfixed = Mfixed, Resultrandom = Mrandom,Resultweighted = Mweighted, Metadata = Msummdata, ProblemSNP = MP))

}


#' ClumpLD: Clumping SNPs using linkage disequilibrium between SNPs
#' @description  This function, which is based on empirical estimations of linkage disequilibrium between SNPs, groups the SNP-based results across one or more datasets or analysis. This approach can be used in two basic scenarios: (i) To summarize the top X single SNP findings from a genome-wide scan as fewer clusters of connected SNPs (i.e., to assess how many independent loci are associated). (ii) To give researchers a simple approach to merge sets of data from multiple studies when those studies may have used various marker sets for genotyping.
#' The clumping process begins with the index SNPs that are significant at threshold p1 and have not yet been clumped. It then creates clumps of all additional SNPs that are within a specified kb of the index SNP and that are in linkage disequilibrium with the index SNP based on an r-squared threshold. Following that, these SNPs are filtered based on the outcome for that SNP. As this method is greedy (1), each SNP will, at most, only appear in one clump. The P value and ALLELES would always, at random, be chosen from the first input file if the same SNP appeared in several input files in SNPdata argument. Instead of the best p-value, the function refer to the SNP that has the strongest LD to the index as the best proxy. Based on the genotype data, the SNP with the highest LD will be the same for all input files.
#' @author Banabithi Bose
#'
#' @param DataDir A character string for the file path of the input PLINK binary files.
#' @param finput Character string, specifying the prefix of the input PLINK binary files which will be used to calculate linkage disequilibrium between the SNPs. This actual genotype data may or may not be the same dataset that was used to generate the summary statistics. This file needs to be in DataDir.
#' @param SNPdata a list of R dataframes containing a single or multiple summary statistics with SNP and P (i.e., p-values) in mandatory column headers. Other columns could be present.
#' @param ResultDir A character string for the file path where all output files will be stored. The default is tempdir().
#' @param clump_p1 Numeric value, specifying the significance threshold for index SNPs. The default is 0.0001.
#' @param clump_p2 Numeric value, specifying the secondary significance threshold for clumped SNPs. The default is 0.01
#' @param clump_r2 Numeric value, specifying the LD threshold for clumping. The default is 0.50.
#' @param clump_kb Integer value, specifying the physical distance threshold in base-pair for clumping. The default is 250.
#' @param byCHR Boolean value, 'TRUE' or 'FALSE', specifying whether to perform the clumping chromosome-wise.
#'
#' @return A dataframe with twelve columns:'CHR'(Chromosome code), 'F'(SNPdata fileset code as 1,2,...), 'SNP'(SNP identifier), 'BP'(Physical position of SNP inbase-pairs),'TOTAL'(Total number of other SNPs in clump i.e. passing clump_kb and clump_r2 thresholds), 'NSIG'(Number of clumped SNPs that are not significant p > 0.05), 'S05'(Number of clumped SNPs 0.01 < p < 0.05), 'S01'(Number of clumped SNPs 0.001 < p < 0.01), 'S001'(Number of clumped SNPs 0.0001 < p < 0.001), 'S0001'(Number of clumped SNPs p < 0.0001), 'SP2'(List of SNPs names and SNPdata code) clumped and significant at clump_p2 threshold.
#' @export
#'
#' @references (1) Purcell et. al.,(2007). PLINK: A Tool Set for Whole-Genome Association and Population-Based Linkage Analyses. The American Journal of Human Genetics, 81(3), 559-575. https://doi.org/10.1086/519795.

#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir = tempdir()
#' finput <- "GXwasR_example"
#' data("GXwasRData")
#' #load(paste0("/projects/b1137/BBose/GXwasR/data/Summary_Stat_Ex1.Rda"))
#' #load(paste0("/projects/b1137/BBose/GXwasR/data/Summary_Stat_Ex2.Rda"))
#' SNPdata <- list(Summary_Stat_Ex1,Summary_Stat_Ex2)
#' clump_p1 = 0.0001
#' clump_p2 = 0.001
#' clump_r2 = 0.5
#' clump_kb = 250
#' byCHR = TRUE
#' clumpedResult <- ClumpLD(DataDir,finput,SNPdata,ResultDir,clump_p1,clump_p2,clump_r2,clump_kb, byCHR)

ClumpLD <- function(DataDir,finput,SNPdata,ResultDir = tempdir(),clump_p1,clump_p2,clump_r2,clump_kb, byCHR = TRUE){

  if (file.exists(paste0(DataDir, "/", finput, ".bed")) &&
      file.exists(paste0(DataDir, "/", finput, ".bim")) &&
      file.exists(paste0(DataDir, "/", finput, ".fam"))) {

    setupPlink(ResultDir)
  }else{
    writeLines(
      "There are no Plink files in DataDir.\nPlease specify correct directory path with input Plink files."
    )
    #return()
  }

  for (i in 1:length(SNPdata)) {
    print(i)
    write.table(SNPdata[[i]], paste0(ResultDir,"/","SNPdata_",i), row.names=FALSE, col.names = TRUE, quote = FALSE)
  }

  SumData <- list.files(paste0(ResultDir,"/"),pattern = "SNPdata_")

  if (byCHR == TRUE){
    bimfile <- read.table(paste0(DataDir, "/", finput, ".bim"))

    chrnum <- 1:length(unique(bimfile$V1))

    chrwiseLD <- function(chrnum){

      chromosome <- unique(bimfile$V1)[chrnum]
      print(paste0("Running LD clumping for chromosome ",chromosome))

      invisible(sys::exec_wait(
        paste0(ResultDir, "/","./plink"),
        args = c(
          "--bfile",paste0(DataDir, "/",finput),
          "--chr",chromosome,
          "--clump", paste0(ResultDir,"/",SumData),
          "--clump-p1",clump_p1,
          "--clump-p2", clump_p2,
          "--clump-r2", clump_r2,
          "--clump-kb", clump_kb,
          "--clump-best",
          "--clump-index-first",
          "--clump-allow-overlap",
          "--clump-snp-field", "SNP",
          "--clump-field", "P",
          #"--clump-verbose",
          "--out",
          paste0(ResultDir, "/","ClumpLD"),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))

      if (file.exists(paste0(ResultDir,"/","ClumpLD.clumped"))){
        #ldc <- read.table(paste0(paste0(ResultDir, "/","ClumpLD.clumped.best")),header = T)
        #ldc1 <- as.data.frame(cbind(chromosome,ldc))
        #colnames(ldc1) <- c("CHROM","INDEX","PSNP","RSQ","KB","P","ALLELES","F" )
        ldc2 <- as.data.frame(read.table(paste0(paste0(ResultDir, "/","ClumpLD.clumped")),header = T))
        invisible(do.call(file.remove,list(paste0(ResultDir,"/","ClumpLD.clumped"))))
        return(ldc2)

      } else {
        print(paste0("No significant clump results for chromosome ",chromosome))

        ldc2 <- data.frame("","","","","","","","","","","","")
        colnames(ldc2) <- c("CHR","F","SNP","BP","P","TOTAL","NSIG","S05","S01","S001","S0001","SP2")
        return(ldc2)

      }

    }

    LDC <- data.table::rbindlist(lapply(chrnum,chrwiseLD))
    LDC <- LDC[!apply(LDC == "", 1, all),]
    return(LDC)

  }else{

    invisible(sys::exec_wait(
      paste0(ResultDir, "/","./plink"),
      args = c(
        "--bfile",
        paste0(DataDir, "/",finput),
        "--clump", paste0(ResultDir,"/",SumData),
        "--clump-p1",clump_p1,
        "--clump-p2", clump_p2,
        "--clump-r2", clump_r2,
        "--clump-kb", clump_kb,
        "--clump-best",
        "--clump-index-first",
        "--clump-allow-overlap",
        "--clump-snp-field", "SNP",
        "--clump-field", "P",
        "--out",
        paste0(ResultDir, "/","ClumpLD"),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))

    if(file.exists(paste0(ResultDir,"/","ClumpLD.clumped"))){
      #ldc <- read.table(paste0(paste0(ResultDir, "/","ClumpLD.clumped.best")),header = T)
      #ldc1 <- as.data.frame(cbind(chromosome,ldc))
      #colnames(ldc1) <- c("CHROM","INDEX","PSNP","RSQ","KB","P","ALLELES","F" )
      ldc2 <- read.table(paste0(paste0(ResultDir, "/","ClumpLD.clumped")),header = T)
      ldc2 <- ldc2[!apply(ldc2 == "", 1, all),]
      return(ldc2)
      #return(list(Clumped.Best = ldc1,Clumped.Detail = ldc2))
    } else {
      print(paste0("No significant clump results"))
    }
  }
  ftemp <- list.files(paste0(ResultDir,"/"),pattern = "SNPdata_")
  invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))
  ftemp <- list.files(paste0(ResultDir,"/"),pattern = "Clump")
  invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))

}


#' DiffZeroOne: Assessing the Z-score for deviation from one and zero.
#'
#' @description This function tests the null hypothesis that a measured statistics (example: genetic correlation, rg for a trait) < 1 using a 1-tailed test compared with a normal distribution (z = (1 − measure statistics)/Standard error).
#' For multiple tests, users are encouraged to apply a Bonferroni multiple-testing correction.
#'
#' @param inputdata A dataframe object, contaning three columns, namely, 'Trait' (i.e., the phenotype of interest), 'Stat' (i.e., the measured statistics) and 'SE' (i.e., the standard error of the measured statistics).
#' @param diffzero Boolean value, 'TRUE' or 'FALSE', specifying to perform diviation from 0 test.
#' @param diffone Boolean value, 'TRUE' or 'FALSE', specifying to perform diviation from 1 test.
#'
#' @return A dataframe with columns, 'Trait','Stat','SE', 'P0' (i.e, p-value for deviation from zero test) and 'P1' (i.e., p-value for deviation from 1 test)
#' @export
#'
#' @examples
#' #load("/projects/b1137/BBose/GXwasR/data/Example_rgdata.Rda")
#' data("GXwasRData")
#' colnames(Example_rgdata) <- c("Trait","Stat","SE")
#' inputdata = Example_rgdata
#' x <- DiffZeroOne(inputdata,FALSE,TRUE)
DiffZeroOne <- function(inputdata, diffzero = c('TRUE','FALSE'), diffone = c('TRUE','FALSE')){

  if (diffzero[1] == FALSE & diffone[1] == FALSE){
    print("Both diffzero and diffone cannot set to be FALSE.")
    return()
  } else if (diffzero[1] == TRUE && diffone[1] == TRUE){
    inputdata1 <- inputdata
    inputdata$x <- as.numeric(inputdata$Stat)
    inputdata$y <- as.numeric(inputdata$SE)
    inputdata$Zscore <- inputdata$x/inputdata$y
    inputdata1$p0 <- pnorm(inputdata$Zscore, mean = 0, sd = 1, lower.tail = FALSE) #one-tailed p-value

    inputdata$x <- 1-as.numeric(inputdata$Stat)
    inputdata$y <- as.numeric(inputdata$SE)
    inputdata$Zscore <- inputdata$x/inputdata$y
    inputdata1$p1 <- pnorm(inputdata$Zscore, mean = 0, sd = 1, lower.tail = FALSE) #one-tailed p-value

  } else if (diffzero[1] == FALSE && diffone[1] == TRUE){
    inputdata1 <- inputdata
    inputdata$x <- 1-as.numeric(inputdata$Stat)
    inputdata$y <- as.numeric(inputdata$SE)
    inputdata$Zscore <- inputdata$x/inputdata$y
    inputdata1$p1 <- pnorm(inputdata$Zscore, mean = 0, sd = 1, lower.tail = FALSE) #one-tailed p-value

  }else if (diffzero[1] == TRUE && diffone[1] == FALSE){
    inputdata1 <- inputdata
    inputdata$x <- as.numeric(inputdata$Stat)
    inputdata$y <- as.numeric(inputdata$SE)
    inputdata$Zscore <- inputdata$x/inputdata$y
    inputdata1$p0 <- pnorm(inputdata$Zscore, mean = 0, sd = 1, lower.tail = FALSE) #one-tailed p-value

  }

  return(inputdata1)

}


#' SexDiffZscore: Z-score-based sex difference test.
#'
#' @description This function calculates the difference in any kind of measured entities,(example: including SNP heritability estimate, genetic correlation, and GWAS β values) between sexes using a Z-score and its associated p-value statistic.
#' When STAT/SE is normally distributed and the test statistics are independent in sex, the test is well calibrated. If the statistics are positively correlated, this test is conservative (1).
#' We could define SNPs with SDEs as those variants at the extreme ends of the distribution with an absolute value of the Z-score greater than 3(|Z-score| > 3), which is roughly equivalent to p <10−3, and represents 0.3% of all tested SNPs.

#' @param inputdata A dataframe with five columns, such as, 'ID' (i.e., SNP ID or the phenotype of interest, etc.),
#' 'Fstat' (i.e., the measured statistics in females), 'Fse' (i.e., the standard error of the measured statistics in females),
#' 'Mstat' (i.e., the measured statistics in males), 'Mse' (i.e., the standard error of the measured statistics in males).
#'
#'
#' @return A dataframe with Zscore (i.e., Z-score), p (i.e., p-value) and adjP (i.e., Bonferroni corrected p-value) columns including other columns of the input dataframe.
#'
#' @export
#'
#' @examples
#' #load("/projects/b1137/BBose/GXwasR/data/Example_h2data.Rda")
#' data("GXwasRData")
#' inputdata = Example_h2data
#' x <- SexDiffZcrore(inputdata)
#'
SexDiffZcrore <- function(inputdata){

  inputdata$x <- as.numeric(inputdata$Fstat)-(as.numeric(inputdata$Mstat))

  inputdata$y <- sqrt((as.numeric(inputdata$Fse))^2+(as.numeric(inputdata$Mse))^2)
  inputdata$Zscore <- inputdata$x/inputdata$y
  #inputdata$p <- 2*pnorm(inputdata$Zscore, mean = 0, sd = 1, lower.tail = TRUE) #Two-tailed p-value

  inputdata <- inputdata[,c("ID","Fstat","Fse","Mstat","Mse","Zscore" )]
  #inputdata$adjP <- p.adjust(inputdata$p, method = "bonferroni", n = nrow(inputdata))
  return(inputdata)

}

#' GeneticCorrBT: Computing genetic correlation between two traits.
#'
#' @description This function computes genetic correlation, a quantitative genetic measure that describes the genetic link between two traits and has been predicted to indicate pleiotropic gene activity or correlation between causative loci in two traits. For example, it does a bivariate GREML analysis to determine the genetic association between two quantitative traits, two binary disease traits from case-control studies, and between a quantitative trait and a binary disease trait following (1,2). If users want, this function gives the flexibility to compute the genetic correlation chromosome-wise.
#'
#' @param DataDir A character string for the file path of the all the input files.
#' @param ResultDir A character string for the file path where all output files will be stored. The default is tempdir().
#' @param finput Character string, specifying the prefix of the input PLINK binary files for the genotype data. This file needs to be in DataDir.
#' @param byCHR Boolean value, 'TRUE' or 'FALSE', specifying whether the analysis will be performed chromosome wise or not. The default is FALSE.
#' @param REMLalgo Integer value of 0, 1 or 2, specifying the algorithm to run REML iterations, 0 for average information (AI), 1 for Fisher-scoring and 2 for EM. The default option is 0, i.e. AI-REML (1).
#' @param nitr Integer value, specifying the number of iterations for performing the REML. The default is 100.
#' @param phenofile A character string, specifying the name of the mandatory phenotype file for Bivar RELM. is a plain text file with no header line; columns family ID, individual ID and two trait columns. For binary trait, the phenotypic value should be coded as 0 or 1, then it will be recognized as a case-control study (0 for controls and 1 for cases). Missing value should be represented by "-9" or "NA".
#' @param cat_covarfile A character string, specifying the name of the categorical covariate file which is a plain text file with no header line; columns are family ID, individual ID and discrete covariates. The default is NULL. This file needs to be in DataDir.
#' @param quant_covarfile A character string, specifying the name of the quantitative covariate file which is a plain text file with no header line; columns are family ID, individual ID and continuous covariates. The default is NULL. This file needs to be in DataDir.
#' @param partGRM Boolean value, 'TRUE' or 'FALSE', specifying whether the GRM will be partitioned into n parts (by row) in GREML model. The default is FALSE.
#' @param autosome Boolean value, 'TRUE' or 'FALSE', specifying whether estimate of heritability will be done for autosomes or not. The default is ‘TRUE’.
#' @param Xsome Boolean value, 'TRUE' or 'FALSE', specifying whether estimate of heritability will be done for X chromosome or not. The default is ‘TRUE’.
#' @param nGRM Integer value, specifying the number of the partition of the GRM in GREML model. The default is 3.
#' @param cripticut Numeric value, specifying the threshold to create a new GRM of "unrelated" individuals in GREML model. The default is arbitrary chosen as 0.025 following (1).
#' @param minMAF Positive numeric value (< maxMAF), specifying the minimum threshold for the MAF filter of the SNPs in the Bivariate GREML model.
#' @param maxMAF Positive numeric value (minMAF,1), specifying the maximum threshold for the MAF filter of the SNPs in the Bivariate GREML model.
#' @param excludeResidual Boolean value, 'TRUE' or 'FALSE', specifying whether to drop the residual covariance from the model. Recommended to set this TRUE if the traits were measured on different individuals. The default is FALSE.
#' @param ncores Integer value, specifying the number of cores to be used.
#'
#' @return A dataframe with minimum three columns, such as "Source" (i.e., source of heritability),
#' "Variance" (i.e., estimated heritability), and "SE" (i.e., standard error of the estimated heritability).
#' Source column will have rows, such as V(G)_tr1 (genetic variance for trait 1),V(G)_tr2 (genetic variance for trait 2),
#' C(G)_tr12 (genetic covariance between traits 1 and 2),V(e)_tr1 (residual variance for trait 1), V(e)_tr2 (residual variance for trait 2),
#' C(e)_tr12 (residual covariance between traits 1 and 2), Vp_tr1 (proportion of variance explained by all SNPs for trait 1),
#' Vp_tr2 (proportion of variance explained by all SNPs for trait 2), V(G)/Vp_tr1 (phenotypic variance for trait 1),
#' V(G)/Vp_tr2 (phenotypic variance for trait 2), rG (genetic correlation) and n (sample size). In case of chromosome-wise
#' analysis, there will be 'chromosome' column for chromosome code.
#'
#'
#' @export
#'
#' @references
#'(1)	Yang et al. (2011) GCTA: a tool for Genome-wide Complex Trait Analysis. Am J Hum Genet. 88(1): 76-82.
#'(2) Lee et al. (2012) Estimation of pleiotropy between complex diseases using SNP-derived genomic relationships and restricted maximum likelihood. Bioinformatics, 28: 2540-2542. PubMed ID: 22843982.
#'
#' @examples
#' #Not Run
# DataDir <- system.file("extdata", package = "GXwasR")
# ResultDir <- tempdir()
# finput <- "GXwasR_example"
# byCHR = FALSE
# REMLalgo = 0
# nitr <- 3
# ncores <- 3
# #phenofile <- "GCphenofile.phen"
# #load("/projects/b1137/BBose/GXwasR/data/Example_phenofile.Rda")
# data("GXwasRData")
# phenofile <- Example_phenofile #Cannot be NULL, the interested phenotype column should be labeled as
#
# cat_covarfile <- NULL
# quant_covarfile <- NULL
# partGRM <- FALSE #Partition the GRM into m parts (by row),
# autosome = TRUE
# Xsome <- TRUE
# cripticut = 0.025
# minMAF <- 0.01 # if MAF filter apply
# maxMAF <- 0.04
# excludeResidual = "TRUE"
#
# GC <- GeneticCorrBT(DataDir = DataDir, ResultDir = ResultDir, finput = finput, byCHR = TRUE,
#                     REMLalgo = 0, nitr = 10, phenofile = phenofile, cat_covarfile = NULL, quant_covarfile = NULL,
#                     partGRM = FALSE, autosome = TRUE, Xsome = TRUE, nGRM = 3,
#                     cripticut = 0.025, minMAF = NULL, maxMAF = NULL,excludeResidual = TRUE, ncores = ncores)

GeneticCorrBT <- function(DataDir, ResultDir, finput,byCHR = FALSE,
                          REMLalgo = c(0,1,2), nitr = 100, phenofile, cat_covarfile = NULL, quant_covarfile = NULL,
                          partGRM = FALSE, autosome = TRUE, Xsome = TRUE, nGRM = 3,
                          cripticut = 0.025, minMAF = NULL, maxMAF = NULL,
                          excludeResidual = FALSE, ncores = 1){

  #library(data.table)

  if (file.exists(paste0(DataDir, "/", finput, ".bed")) &&
      file.exists(paste0(DataDir, "/", finput, ".bim")) &&
      file.exists(paste0(DataDir, "/", finput, ".fam"))) {

    setupGCTA(ResultDir)

  } else{
    writeLines(
      "There are no Plink files in DataDir.\nPlease specify correct directory path with input Plink files."
    )
    #return()
  }
  ##ComputeBivarREMLone phenofile
  write.table(phenofile,file = paste0(ResultDir,"/","GCphenofile"), row.names = FALSE, quote = FALSE)

  if (byCHR == FALSE){

    if (autosome == TRUE && Xsome == FALSE){

      ## Compute GRM
      ComputeGRMauto(DataDir = DataDir, ResultDir = ResultDir, finput = finput,
                     partGRM = partGRM, nGRM = nGRM, cripticut = cripticut,  minMAF = minMAF, maxMAF = maxMAF)

      ## Compute REML
      herit_result <- ComputeBivarREMLone(DataDir = DataDir, ResultDir = ResultDir, REMLalgo = REMLalgo, nitr = nitr, phenofile = "GCphenofile", cat_covarfile = cat_covarfile,
                                          quant_covarfile = quant_covarfile, excludeResidual = excludeResidual, grmfile = "test",ncores = ncores)

      return(herit_result)

    }else if (autosome == TRUE && Xsome == TRUE){
      ## Compute GRM Autosome
      ComputeGRMauto(DataDir = DataDir, ResultDir = ResultDir, finput = finput,
                     partGRM = partGRM, nGRM = nGRM, cripticut = cripticut,  minMAF = minMAF, maxMAF = maxMAF)
      ## Compute GRM X
      ComputeGRMX(DataDir = DataDir, ResultDir = ResultDir, finput = finput,
                  partGRM = partGRM, nGRM = nGRM, minMAF = minMAF, maxMAF = maxMAF)

      ## Compute REML
      herit_result <- ComputeBivarREMLmulti(DataDir = DataDir, ResultDir = ResultDir, REMLalgo = REMLalgo, nitr = nitr, phenofile = phenofile, cat_covarfile = cat_covarfile,
                                            quant_covarfile = quant_covarfile, excludeResidual = excludeResidual, grmfile = "multi_GRMs.txt",ncores = ncores)

      return(herit_result)

    }else if (autosome == FALSE && Xsome == TRUE){
      ## Compute GRM X
      ComputeGRMX(DataDir = DataDir, ResultDir = ResultDir, finput = finput,
                  partGRM = partGRM, nGRM = nGRM, minMAF = minMAF, maxMAF = maxMAF)
      ## Compute REML X
      herit_result <- ComputeBivarREMLone(DataDir = DataDir, ResultDir = ResultDir, REMLalgo = REMLalgo, nitr = nitr, phenofile = "GCphenofile", cat_covarfile = cat_covarfile,
                                          quant_covarfile = quant_covarfile, excludeResidual = excludeResidual, grmfile = "xtest", ncores = ncores)

      return(herit_result)

    }else{
      print("autosome and Xsome cannot be set as FALSE together.")
    }
  }else{


    bimfile <- read.table(paste0(DataDir, "/", finput, ".bim"))

    chrnum <- 1:length(unique(bimfile$V1))
    #chrnum <- 1:3

    chrwiseRELM <- function(chrnum){

      chromosome <- as.integer(unique(bimfile$V1)[chrnum])

      print(paste0("Processing chromosome ",chromosome))


      if (chromosome == 23){
        ## Compute GRM X

        ComputeGRMX(DataDir = DataDir, ResultDir = ResultDir, finput = finput,
                    partGRM = partGRM, nGRM = nGRM, minMAF = minMAF, maxMAF = maxMAF)

        ## Compute REML X
        herit_result <- ComputeBivarREMLone(DataDir = DataDir, ResultDir = ResultDir, REMLalgo = REMLalgo, nitr = nitr, phenofile = "GCphenofile", cat_covarfile = cat_covarfile,
                                            quant_covarfile = quant_covarfile, excludeResidual = excludeResidual, chr = 23, grmfile = "xtest", ncores = ncores)

        herit_result <- data.table::as.data.table(cbind(chromosome,herit_result))


      }else{

        ## Compute GRM
        ComputeGRMauto(DataDir = DataDir, ResultDir = ResultDir, finput = finput,
                       partGRM = partGRM, nGRM = nGRM, cripticut = cripticut,
                       minMAF = minMAF, maxMAF = maxMAF, ByCHR = byCHR, CHRnum = chromosome)

        ## Compute REML
        herit_result <- ComputeBivarREMLone(DataDir = DataDir, ResultDir = ResultDir, REMLalgo = REMLalgo, nitr = nitr, phenofile = "GCphenofile", cat_covarfile = cat_covarfile,
                                            quant_covarfile = quant_covarfile, excludeResidual = excludeResidual, chr = chromosome,grmfile = "test", ncores = ncores)


        herit_result <- data.table::as.data.table(cbind(chromosome,herit_result))


      }

      return(herit_result)
    }

    result <- rbindlist(lapply(chrnum,chrwiseRELM))
    result <- na.omit(result)
    ftemp <- list.files(paste0(ResultDir,"/"),pattern = "bireml")
    invisible(file.remove(paste0(ResultDir,"/",ftemp)))
    ftemp <- list.files(paste0(ResultDir,"/"),pattern = "grm")
    invisible(file.remove(paste0(ResultDir,"/",ftemp)))
    ftemp <- list.files(paste0(ResultDir,"/"),pattern = "test")
    invisible(file.remove(paste0(ResultDir,"/",ftemp)))
    ftemp <- list.files(paste0(ResultDir,"/"),pattern = "gcta")
    invisible(file.remove(paste0(ResultDir,"/",ftemp)))
    invisible(file.remove(paste0(ResultDir,"/","GCphenofile")))
    invisible(file.remove(paste0(ResultDir,"/","multi_GRMs.txt")))
    return(result)
  }

}


#' SexRegress: Performing linear regression analysis with quantitative response variable.
#'
#' @description This function could be used to check association of two variables. For instance, PRS with sex.
#'
#' @param fdata R dataframe object. The column with header "response" should contain the response variable. All other column are the regressor.
#' @param regressor_index Integer value, specifying the column number of the main regressor variable.
#' @param response_index Integer value, specifying the column number of the response variable.
#'
#' @return Numeric value containing the regression estimate ("Estimate"), standard error ("Std. Error"), statistics ("t value") and p-value ("Pr(>|t|)")
#'
#' @importFrom stats lm
#'
#' @export
#'
#' @examples
#' #load("/projects/b1137/BBose/GXwasR/data/Regression_Ex.Rda")
#' data("GXwasRData")
#' fdata = Regression_Ex
#' fdata$SEX <- as.factor(as.character(fdata$SEX))
#' response_index <- 1
#' regressor_index <- 2
#'
#' x <- SexRegress(fdata, regressor_index, response_index)
#'
SexRegress <- function(fdata, regressor_index, response_index){

  names(fdata)[response_index]<- "response"
  nullfdata <- fdata[,-regressor_index]

  null.model <- stats::lm(nullfdata$response~., data = nullfdata)
  model <- stats::lm(fdata$response~., data= fdata)
  # model R2 is obtained as
  null.r2 <- summary(null.model)$r.squared
  model.r2 <- summary(model)$r.squared

  # R2 of response is simply calculated as the model R2 minus the null R2
  response.r2 <- model.r2-null.r2
  model.result <- summary(model)$coeff[regressor_index,]
  return(model.result)
}

#' FilterRegion: Filter chromosomal regions.
#'
#' @author Banabithi Bose
#'
#' @description Filtering Pseudo-Autosomal Region (PAR), X-transposed region (XTR), Ampliconic, filter based on chromosome code or user-defined regions from input plink files.
#'
#' @param DataDir A character string for the file path of the input PLINK binary files.
#' @param ResultDir A character string for the file path where all output files will be stored. The default is tempdir().
#' @param finput Character string, specifying the prefix of the input PLINK binary files.
#' @param foutput Character string, specifying the prefix of the output PLINK binary files if the filtering option for the SNPs is chosen. The default is "FALSE".
#' @param CHRX Boolean value, 'TRUE' or 'FALSE' to filter/flag regions from chromosome X. The default is TRUE.
#' @param CHRY Boolean value, 'TRUE' or 'FALSE' to filter/flag regions from chromosome X. The default is FALSE.
#' @param filterPAR Boolean value, 'TRUE' or 'FALSE' to filter out PARs from input plink file. The default is TRUE.
#' @param filterXTR Boolean value, 'TRUE' or 'FALSE' to filter out XTRs from input plink file. The default is TRUE.
#' @param filterAmpliconic Boolean value, 'TRUE' or 'FALSE' to filter out Ampliconic regions from input plink file. The default is TRUE.
#' @param regionfile Character string, specifying the name of the .txt file containing the user-defined regions to be filtered out from input plink file in bed format. The default is FALSE.
#' @param filterCHR Vector value with positive integer, specifying the chromosome code to filter/flag the SNPs. The default is 0, means no filtering based on chromosome code. For non-zero values of this argument, the function will only consider the chromosome code to filter or flag. All other filtering will not work.
#' @param Hg Character value, '19', or '38', specifying which genome build to use for PAR, XTR and Ampliconic regions. The default is Hg = "19".
#' @param exclude Boolean value, 'TRUE' or 'FALSE' to filter and flag or only flag the SNPs. The default is TRUE.
#'
#' @return A list of three dataframes: PAR containing SNPs from PAR regions; XTR containing SNPs from XTR region and Ampliconic containing SNPs from Ampliconic region.
#'
#' For non-zero value of filterCHR, a dataframe containing the excluded/flagged SNPs will be returned.
#'
#' For exclude = TRUE, two sets of plink binary files will be produced in ResultDir. One set will have the remaining SNPs after filtering and other one will have the discarded SNPs.
#' @export
#'
#' @examples
#' # Not Run
# DataDir <- system.file("extdata", package = "GXwasR")
# ResultDir <- tempdir()
# finput <- "GXwasR_example_imputed"
# foutput <- "PostimputeEX_QC1"
# x <- FilterRegion(DataDir = DataDir, ResultDir = ResultDir, finput = finput, foutput = foutput, CHRX = TRUE, CHRY = FALSE, filterPAR = TRUE, filterXTR = TRUE, filterAmpliconic = TRUE, regionfile = FALSE, filterCHR = 0, Hg = "38",exclude = TRUE)

FilterRegion <-
  function(DataDir,
           ResultDir,
           finput,
           foutput,
           CHRX = TRUE,
           CHRY = FALSE,
           filterPAR = TRUE,
           filterXTR = TRUE,
           filterAmpliconic = TRUE,
           regionfile = FALSE,
           filterCHR = 0,
           Hg = "19",
           exclude = TRUE) {

    DataDir1 <- system.file("extdata", package = "GXwasR")

    if (file.exists(paste0(DataDir, "/", finput, ".bed")) &&
        file.exists(paste0(DataDir, "/", finput, ".bim")) &&
        file.exists(paste0(DataDir, "/", finput, ".fam"))) {

      setupPlink(ResultDir)

    } else{
      writeLines(
        "There are no Plink files in DataDir.\nPlease specify correct directory path with input Plink files."
      )
    }

    if (filterCHR[1] == 0){

      if (Hg == "19"){
        HG = "hg19"
      }else{
        HG = "GRCh38"
      }

      if (CHRX == TRUE){
        CHR <-"chrX"
      }else if (CHRY == TRUE){
        CHR <-"chrY"
      }

      if (exclude == TRUE){


        if (regionfile == FALSE) {

          if (CHRX == TRUE && CHRY == TRUE){

            x1 <-
              as.data.frame(
                read.table(
                  file = paste0(DataDir1,"/","X","_genomic_features_",HG,".bed.txt"),
                  header = FALSE,
                  sep = ""
                )
              )

            x2 <-
              as.data.frame(
                read.table(
                  file = paste0(DataDir1,"/","Y","_genomic_features_",HG,".bed.txt"),
                  header = FALSE,
                  sep = ""
                )
              )

            x <- rbind(x1, x2)

          }else{
            ## Need to use extdata for DataDir
            x <-
              as.data.frame(
                read.table(
                  file = paste0(DataDir1,"/",CHR,"_genomic_features_",HG,".bed.txt"),
                  header = FALSE,
                  sep = ""
                )
              )
          }

          #############
          #############
          if (filterPAR == TRUE){
            y <- x[x$V4 == "PAR1" | x$V4 == "PAR2", ]

            write.table(
              y,
              file = paste0(ResultDir,"/par_region.txt"),
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE,
              eol = "\r\n",
              sep = " "
            )

            invisible(sys::exec_wait(
              paste0(ResultDir,"/","./plink"),
              args = c(
                "--bfile",
                paste0(DataDir,"/",finput),
                "--extract",
                "range",
                paste0(ResultDir,"/","par_region.txt"),
                "--make-bed",
                "--out",
                paste0(ResultDir,"/",foutput,"_par_region"),
                "--silent"
              ),
              std_out = FALSE,
              std_err = FALSE
            ))

            if (file.exists(paste0(ResultDir,"/",foutput,"_par_region.bim"))){
              par_snps <-  read.table(paste0(ResultDir,"/",foutput,"_par_region.bim"))
              colnames(par_snps) <- c("CHR","SNP","START","END","A1","A2")
              print(paste0("PAR SNPs:",length(unique(par_snps$SNP))))

            }else{
              print("There is no PAR region in the input data. Argument filterPAR cannot set to be TRUE.")
              print("Changing it as filterPAR = FALSE")
              filterPAR = FALSE
              par_snps <- NULL
            }
          }else{
            par_snps <- NULL
          }

          if (filterXTR == TRUE){

            y <- x[x$V4 == "XTR", ]

            write.table(
              y,
              file = paste0(ResultDir,"/xtr_region.txt"),
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE,
              eol = "\r\n",
              sep = " "
            )

            invisible(sys::exec_wait(
              paste0(ResultDir,"/","./plink"),
              args = c(
                "--bfile",
                paste0(DataDir,"/",finput),
                "--extract",
                "range",
                paste0(ResultDir,"/","xtr_region.txt"),
                "--make-bed",
                "--out",
                paste0(ResultDir,"/",foutput,"_xtr_region"),
                "--silent"
              ),
              std_out = FALSE,
              std_err = FALSE
            ))

            if (file.exists(paste0(ResultDir,"/",foutput,"_xtr_region.bim"))){
              xtr_snps <-  read.table(paste0(ResultDir,"/",foutput,"_xtr_region.bim"))
              colnames(xtr_snps) <- c("CHR","SNP","START","END","A1","A2")
              print(paste0("XTR SNPs:",length(unique(xtr_snps$SNP))))
            }else{
              print("There is no XTR region in the input data. Argument filterXTR cannot set to be TRUE.")
              print("Changing it as filterXTR = FALSE")
              filterXTR = FALSE
              xtr_snps <- NULL
            }
          }else {
            xtr_snps <- NULL
          }

          if (filterAmpliconic == TRUE){

            y <- x[grep("^Ampliconic", x$V4),]

            write.table(
              y,
              file = paste0(ResultDir,"/ampliconic_region.txt"),
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE,
              eol = "\r\n",
              sep = " "
            )


            invisible(sys::exec_wait(
              paste0(ResultDir,"/","./plink"),
              args = c(
                "--bfile",
                paste0(DataDir,"/",finput),
                "--extract",
                "range",
                paste0(ResultDir,"/","ampliconic_region.txt"),
                "--make-bed",
                "--out",
                paste0(ResultDir,"/",foutput,"_ampliconic_region"),
                "--silent"
              ),
              std_out = FALSE,
              std_err = FALSE
            ))

            if (file.exists(paste0(ResultDir,"/",foutput,"_ampliconic_region.bim"))){
              ampliconic_snps <-  read.table(paste0(ResultDir,"/",foutput,"_ampliconic_region.bim"))
              colnames(ampliconic_snps) <- c("CHR","SNP","START","END","A1","A2")
              print(paste0("Ampliconic SNPs:",length(unique(ampliconic_snps$SNP))))

            }else{
              print("There is no ampliconic region in the input data. Argument filterAmpliconic cannot set to be TRUE.")
              print("Changing it as filterAmpliconic = FALSE")
              filterAmpliconic = FALSE
              ampliconic_snps <- NULL
            }
          }else{
            ampliconic_snps <- NULL
          }

          #############
          #############

          if (filterPAR == TRUE &&
              filterXTR == TRUE && filterAmpliconic == TRUE){
            y = x

          }else if (filterPAR == TRUE &&
                    filterXTR == TRUE && filterAmpliconic == FALSE){

            y <- x[x$V4 == "PAR1" | x$V4 == "PAR2" |x$V4 == "XTR", ]

          }else if (filterPAR == TRUE &&
                    filterXTR == FALSE && filterAmpliconic == TRUE){

            x1 <- x[grep("^Ampliconic", x$V4),]
            x2 <- x[x$V4 == "PAR1" | x$V4 == "PAR2", ]
            y <- rbind(x1,x2)

          }else if (filterPAR == FALSE &&
                    filterXTR == TRUE && filterAmpliconic == TRUE){
            x1 <- x[grep("^Ampliconic", x$V4),]
            x2 <- x[x$V4 == "XTR", ]
            y <- rbind(x1,x2)
          }else if (filterPAR == TRUE &&
                    filterXTR == FALSE && filterAmpliconic == FALSE){
            y <- x[x$V4 == "PAR1" | x$V4 == "PAR2", ]
          }else if (filterPAR == FALSE &&
                    filterXTR == TRUE && filterAmpliconic == FALSE){

            y <- x[x$V4 == "XTR", ]

          }else if (filterPAR == FALSE &&
                    filterXTR == FALSE && filterAmpliconic == TRUE){

            y <- x[grep("^Ampliconic", x$V4),]


          }

          write.table(
            y,
            file = paste0(ResultDir,"/region.txt"),
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            eol = "\r\n",
            sep = " "
          )

          rangefile = paste0(ResultDir,"/region.txt")

        }else{

          rangefile = paste0(DataDir,"/",regionfile)
          par_snps <- NULL
          xtr_snps <- NULL
          ampliconic_snps <- NULL
        }
        #}

        ## Exclude region
        invisible(sys::exec_wait(
          paste0(ResultDir,"/","./plink"),
          args = c(
            "--bfile",
            paste0(DataDir,"/",finput),
            "--exclude",
            "range",
            rangefile,
            "--make-bed",
            "--out",
            paste0(ResultDir,"/",foutput),
            "--silent"
          ),
          std_out = FALSE,
          std_err = FALSE
        ))

        ##Extract region
        #########

        invisible(sys::exec_wait(
          paste0(ResultDir,"/","./plink"),
          args = c(
            "--bfile",
            paste0(DataDir,"/",finput),
            "--extract",
            "range",
            rangefile,
            "--make-bed",
            "--out",
            paste0(ResultDir,"/",foutput,"_snps_extracted"),
            "--silent"
          ),
          std_out = FALSE,
          std_err = FALSE
        ))

        #ftemp <- list.files(paste0(ResultDir,"/"),pattern = "FF")
        #invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))

        bim <- read.table(paste0(ResultDir,"/",foutput,".bim"))
        bim1 <- read.table(paste0(DataDir,"/",finput,".bim"))

        num_marker_excluded <- nrow(bim1)-nrow(bim)
        print(paste0(num_marker_excluded," SNPs are discarded."))

        # marker_excluded <- outersect(bim1$V2,bim$V2)
        # bim2 <- bim1[bim1$V2 %in% marker_excluded,]

        print(paste0("Plink files with passed SNPs are in ",ResultDir," prefixed as ",foutput))
        print(paste0("Plink files with discarded SNPs are in ",ResultDir," prefixed as ",foutput,"_snps_extracted"))


      }else{

        print("SNPs are only flagged for the desired region.")

      }

      if (regionfile == FALSE){
        ftemp <- list.files(paste0(ResultDir,"/"),pattern = "region")
        invisible(do.call(file.remove,list(paste0(ResultDir,"/",ftemp))))
      }

      return(list(PAR = par_snps, XTR = xtr_snps, Ampliconic = ampliconic_snps ))

    }else{

      invisible(sys::exec_wait(
        paste0(ResultDir,"/","./plink"),
        args = c(
          "--bfile",
          paste0(DataDir,"/",finput),
          "--not-chr",
          filterCHR,
          "--make-bed",
          "--out",
          paste0(ResultDir,"/",foutput),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))

      invisible(sys::exec_wait(
        paste0(ResultDir,"/","./plink"),
        args = c(
          "--bfile",
          paste0(DataDir,"/",finput),
          "--chr",
          filterCHR,
          "--make-bed",
          "--out",
          paste0(ResultDir,"/",foutput,"_snps_extracted"),
          "--silent"
        ),
        std_out = FALSE,
        std_err = FALSE
      ))
    }


    if (exclude == TRUE){
      bim <- read.table(paste0(ResultDir,"/",foutput,".bim"))
      bim1 <- read.table(paste0(DataDir,"/",finput,".bim"))

      num_marker_excluded <- nrow(bim1)-nrow(bim)
      print(paste0(num_marker_excluded," SNPs are discarded."))

      print(paste0("Plink files with passed SNPs are in ",ResultDir," prefixed as ",foutput))
      print(paste0("Plink files with discarded SNPs are in ",ResultDir," prefixed as ",foutput,"_snps_extracted"))

      #return()

    }else if (exclude == FALSE){
      bim <- read.table(paste0(ResultDir,"/",foutput,".bim"))
      colnames(bim) <- c("CHR","SNP","START","END","A1","A2")
      Flagged_SNPs = bim
      #return(Flagged_SNPs)
    }
    #}

  }


#' PlinkSummary: Summary of plink format genotype dataset
#'
#' @param DataDir A character string for the file path of the input PLINK binary files.
#' @param ResultDir A character string for the file path of the plink program to be set up.
#' @param finput Character string, specifying the prefix of the input PLINK binary files. This file needs to be in DataDir.
#'
#' @return NULL
#' @export
#'
#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#'
#' x <- PlinkSummary(DataDir,ResultDir, finput)
#'
PlinkSummary <- function(DataDir, ResultDir = tempdir(), finput){

  if (file.exists(paste0(DataDir, "/", finput, ".bed")) &&
      file.exists(paste0(DataDir, "/", finput, ".bim")) &&
      file.exists(paste0(DataDir, "/", finput, ".fam"))) {

    setupPlink(ResultDir)

  } else{
    writeLines(
      "There are no Plink files in DataDir.\nPlease specify correct directory path with input Plink files."
    )
  }

  fam <-
    as.data.frame(utils::read.table(file = paste0(DataDir, "/", finput, ".fam")))

  fam$V6 <- as.numeric(as.character(fam$V6))
  fam <- stats::na.omit(fam)
  fam1 <- fam[fam$V5 != 0,]
  fam2 <- fam1[fam1$V6 != 0,]
  fam4 <- fam2[fam2$V6 != -9,]
  print(paste0("Dataset:",finput))
  if (length(unique(fam4$V6))==2){
    print("This is a case-control data.")
    No.of.cases <- nrow(fam4[fam4$V6==2,])#108
    print(paste0("Number of cases:",No.of.cases))
    No.of.controls <- nrow(fam4[fam4$V6==1,])#108
    print(paste0("Number of controls:",No.of.controls))
    No.of.missing.pheno <- nrow(fam[fam$V6==-9|fam$V6==0,])#108
    print(paste0("Number of missing phenotypes:",No.of.missing.pheno))
  }else if (length(unique(fam4$V6))==1){
    print("This dataset contains quantitative trait value.")
    No.of.missing.pheno <- nrow(fam[fam$V6==-9|fam$V6==0,])#108
    print(paste0("Number of missing phenotypes:",No.of.missing.pheno))
    print(paste0("Summary of the phenotype after removing missing values if present:",summary(fam4$V6)))
  }

  No.of.males <- nrow(fam[fam$V5==1,])
  print(paste0("Number of males:",No.of.males))
  No.of.females <- nrow(fam[fam$V5==2,])
  print(paste0("Number of females:",No.of.females))

  bim <-
    as.data.frame(utils::read.table(file = paste0(DataDir, "/", finput, ".bim")))


  No.of.chr <- length(unique(bim$V1))
  print(paste0("Number of chromosomes:",No.of.chr))
  print(paste0("Chr:",unique(bim$V1)))
  No.of.snps <- length(unique(bim$V2))
  print(paste0("Total number of SNPs:",No.of.snps))
  No.of.samples <- length(unique(fam$V2))
  print(paste0("Total number of samples:",No.of.samples))

  return(NULL)
}


#' FilterAllele: Filtering out the multi-allelic variants
#'
#' @author Banabithi Bose
#'
#' @description This function filters out the multi-allelic SNPs from the input dataset.
#' @param DataDir A character string for the file path of the input PLINK binary files.
#' @param ResultDir A character string for the file path where all output files will be stored. The default is tempdir().
#' @param finput Character string, specifying the prefix of the input PLINK binary files.
#' @param foutput Character string, specifying the prefix of the output PLINK binary files. If multi-allelic variants are present, this file will be produced after filtering out these variants.
#'
#' @return NULL
#' After multi-allelic variant filtering, the filtered plink files with only biallelic SNPs will be saved in ResultDir.
#' @export
#'
#' @examples
#' DataDir <- system.file("extdata", package = "GXwasR")
#' ResultDir <- tempdir()
#' finput <- "GXwasR_example"
#' foutput <- "Filter_Test"
#' x <- FilterAllele(DataDir, ResultDir, finput, foutput)

FilterAllele <- function(DataDir, ResultDir, finput, foutput){

  if (file.exists(paste0(DataDir, "/", finput, ".bed")) &&
      file.exists(paste0(DataDir, "/", finput, ".bim")) &&
      file.exists(paste0(DataDir, "/", finput, ".fam"))) {

    setupPlink(ResultDir)

  } else{
    writeLines(
      "There are no Plink files in DataDir.\nPlease specify correct directory path with input Plink files."
    )
  }

  bimf <- read.table(paste0(DataDir,"/", finput,".bim"))
  x1 <- bimf[nchar(bimf[,5]) > 1 | nchar(bimf[,6]) > 1, ,drop = FALSE]

  if (nrow(x1) != 0){

    write.table(x1$V2,file = paste0(ResultDir,"/snps_multiallelic"),quote = F,col.names = F,row.names = F)

  }else{
    print("There is no multi-allelic SNP present in the input dataset.")
    return()
  }

  invisible(sys::exec_wait(
    paste0(ResultDir,"/","./plink"),
    args = c(
      "--bfile",
      paste0(DataDir,"/",finput),
      "--exclude",
      paste0(ResultDir,"/","snps_multiallelic"),
      "--make-bed",
      "--out",
      paste0(ResultDir,"/",foutput),
      "--silent"
    ),
    std_out = FALSE,
    std_err = FALSE
  ))

  bimf1 <- read.table(paste0(ResultDir,"/", foutput,".bim"))

  print(paste0("Input dataset has ",nrow(bimf)," SNPs."))
  print(paste0("Plink files with only biallelic SNPs are in ",ResultDir," prefixed as ",foutput))
  print(paste0("Output dataset has ",nrow(bimf1)," SNPs."))

  return()

}









