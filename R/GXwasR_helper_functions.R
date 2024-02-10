## Function update in 3.0 again
# GXwasR helper functions
# @author Banabithi Bose
#' @importFrom stats na.omit median qchisq binomial dnorm glm na.exclude p.adjust pchisq pnorm lm pt qnorm setNames
#' @importFrom utils download.file read.table unzip write.table untar
#' @importFrom qqman manhattan qq
#' @importFrom grDevices dev.off jpeg
#' @importFrom poolr stouffer bonferroni
#' @importFrom SNPRelate snpgdsBED2GDS snpgdsGetGeno snpgdsOpen
#' @importFrom gdsfmt showfile.gds read.gdsn index.gdsn
#' @importFrom stringr str_sub
#' @importFrom GenomeInfoDb getChromInfoFromUCSC
#' @importFrom bigparallelr nb_cores
#' @importFrom data.table fread as.data.table data.table
#' @importFrom ggplot2 ylim xlim aes geom_smooth geom_point ggsave geom_text geom_rect labs geom_errorbarh scale_y_continuous geom_vline theme element_text ggtitle element_blank scale_fill_manual scale_y_reverse
#' @importFrom bigsnpr snp_ldsc2
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom regioneR toGRanges
#' @importFrom plyranges join_overlap_intersect
#' @importFrom ggpubr ggarrange annotate_figure text_grob stat_cor stat_regline_equation
#' @importFrom graphics abline arrows axis points par
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom Matrix sparseMatrix
#' @importFrom grid unit
#' @importFrom ggrepel geom_text_repel
#' @importFrom BiocStyle html_document

## Function 1
# Installing Plink
setupPlink <- function(wdir) {
  # Helper function to download and unzip files
  downloadAndUnzip <- function(file_url, dest_file, wdir) {
    utils::download.file(destfile = dest_file, file_url, quiet = TRUE)
    utils::unzip(dest_file, exdir = wdir)
  }

  # Helper function to remove files
  removeFiles <- function(files, wdir) {
    full_paths <- file.path(wdir, files)
    invisible(do.call(file.remove, list(full_paths)))
  }

  # Get operating system information
  OS <- Sys.info()['sysname']
  base_url <- "https://s3.amazonaws.com/plink1-assets/"
  os_specific_files <- list(
    Linux = list(file = "plink_linux_x86_64_20220402.zip", exec = "plink", remove = c("LICENSE", "prettify", "toy.map", "toy.ped")),
    Windows = list(file = "plink_win64_20230116.zip", exec = "plink.exe", remove = c("LICENSE", "prettify.exe", "toy.map", "toy.ped")),
    macOS = list(file = "plink_mac_20230116.zip", exec = "plink", remove = c("LICENSE", "prettify", "toy.map", "toy.ped"))
  )

  if (!OS %in% names(os_specific_files)) {
    stop("Unsupported operating system.")
  }

  # Perform the setup
  os_data <- os_specific_files[[OS]]
  file_url <- paste0(base_url, os_data$file)
  dest_file <- file.path(wdir, os_data$file)

  tryCatch({
    downloadAndUnzip(file_url, dest_file, wdir)
    Sys.chmod(file.path(wdir, os_data$exec), mode = "0777", use_umask = TRUE)
    removeFiles(c(os_data$file, os_data$remove), wdir)
    print("Program is set up.")
  }, error = function(e) {
    cat("An error occurred:", e$message, "\n")
  })
}

## Function 2
MFsplitPlink <- function(DataDir, ResultDir, finput, foutput, sex, xplink = FALSE, autoplink = FALSE){

  if (!checkFiles(DataDir,finput)) {
    stop("Missing required Plink files in the specified DataDir.")
  }

  setupPlink(ResultDir)

  if (xplink == FALSE && autoplink == FALSE){

    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--bed",
        paste0(DataDir,"/",finput,".bed"),
        "--bim",
        paste0(DataDir,"/",finput,".bim"),
        "--fam",
        paste0(DataDir,"/",finput,".fam"),
        paste0("--filter-",sex),
        "--make-bed",
        "--out",
        paste0(ResultDir,"/",foutput),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))
  }else if (xplink == TRUE && autoplink == FALSE){

    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--bed",
        paste0(DataDir,"/",finput,".bed"),
        "--bim",
        paste0(DataDir,"/",finput,".bim"),
        "--fam",
        paste0(DataDir,"/",finput,".fam"),
        paste0("--filter-",sex),
        "--chr",23,
        "--make-bed",
        "--out",
        paste0(ResultDir,"/",foutput),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))

  }else if (xplink == FALSE && autoplink == TRUE){
    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--bed",
        paste0(DataDir,"/",finput,".bed"),
        "--bim",
        paste0(DataDir,"/",finput,".bim"),
        "--fam",
        paste0(DataDir,"/",finput,".fam"),
        paste0("--filter-",sex),
        "--not-chr",23,
        "--make-bed",
        "--out",
        paste0(ResultDir,"/",foutput),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))
  }

  print("Stratified test is running")

}

## Function 3
########## Added in 3.0
checkFiles <- function(DataDir,finput) {
  all(file.exists(file.path(DataDir, paste0(finput, c(".bed", ".bim", ".fam")))))
}

## Function 4
########## Added in 3.0
executePlink <- function(args) {
  #globalVariables("ResultDir")
  tryCatch({
    # Redirect stderr to null to suppress warning messages
    stderr_dest <- ifelse(.Platform$OS.type == "windows", "NUL", "/dev/null")
    invisible(sys::exec_wait(file.path(ResultDir, "./plink"), args = args,
                             std_err = stderr_dest))
  }, error = function(e) {
    stop("An error occurred while executing Plink: ", e$message)
  })
}

## Function 5
########## Added in 3.0
analyzePhenotypeData <- function(fam, fam4) {
  unique_pheno <- unique(fam4$V6)
  if (length(unique_pheno) == 2) {
    No.of.cases <- nrow(fam4[fam4$V6 == 2, ])
    No.of.controls <- nrow(fam4[fam4$V6 == 1, ])
    message <- "This is a case-control data."
  } else if (length(unique_pheno) == 1 && unique_pheno != -9) {
    message <- paste0("This dataset contains a single trait value: ", unique_pheno)
  } else if (unique(fam$V6) == -9) {
    message <- "This dataset contains missing trait values for all samples."
  } else {
    message <- "This dataset contains quantitative trait value."
  }

  No.of.missing.pheno <- nrow(fam[fam$V6 == -9 | fam$V6 == 0, ])
  No.of.males <- nrow(fam[fam$V5 == 1, ])
  No.of.females <- nrow(fam[fam$V5 == 2, ])

  print(message)
  print(paste0("Number of males:", No.of.males))
  print(paste0("Number of females:", No.of.females))
  print(paste0("Number of missing phenotypes:", No.of.missing.pheno))
}

## Function 6
########## Added in 3.0
plinkExcludeExtract <- function(DataDir, finput, ResultDir, foutput, region_file_path) {
  # Plink command for excluding SNPs
  plinkArgsExclude <- c(
    "--bfile", file.path(DataDir, finput),
    "--exclude", "range", region_file_path,
    "--allow-no-sex",                       ## Adding in 4.0
    "--make-bed",
    "--out", file.path(ResultDir, foutput),
    "--silent"
  )
  executePlink(plinkArgsExclude)

  # Plink command for extracting SNPs
  plinkArgsExtract <- c(
    "--bfile", file.path(DataDir, finput),
    "--extract", "range", region_file_path,
    "--allow-no-sex",
    "--make-bed",
    "--out", file.path(ResultDir, paste0(foutput, "_snps_extracted")),
    "--silent"
  )
  executePlink(plinkArgsExtract)
}

## Function 7
########## Added in 3.0
# processRegionFile <- function(DataDir, finput, ResultDir, foutput, genomic_feature, region_file_path, HG, exclude) {
#   DataDir1 <- system.file("extdata", package = "GXwasR")
#   region_data <- read.table(file.path(DataDir1, paste0(genomic_feature, "_genomic_features_", HG, ".bed.txt")), header = FALSE, sep = "")
#   write.table(region_data, file = region_file_path, quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\r\n", sep = " ")
#
#   if (exclude) {
#     plinkExcludeExtract(DataDir, finput, ResultDir, foutput, region_file_path)
#   }
# }

## Function 8
########## Added in 3.0
# processCHRFilter <- function(DataDir, finput, ResultDir, foutput, filterCHR, exclude) {
#   if (exclude) {
#     plinkArgsExclude <- c("--bfile", file.path(DataDir, finput), "--not-chr", filterCHR, "--make-bed", "--out", file.path(ResultDir, foutput), "--silent")
#     plinkArgsExtract <- c("--bfile", file.path(DataDir, finput), "--chr", filterCHR, "--make-bed", "--out", file.path(ResultDir, paste0(foutput, "_snps_extracted")), "--silent")
#     executePlink(plinkArgsExclude)
#     executePlink(plinkArgsExtract)
#   }
# }

## Function 9
########## Added in 3.0
# finalProcessing <- function(DataDir, ResultDir, finput, foutput, exclude, regionfile) {
#   bim <- read.table(file.path(ResultDir, paste0(foutput, ".bim")))
#   bim1 <- read.table(file.path(DataDir, paste0(finput, ".bim")))
#   num_marker_excluded <- nrow(bim1) - nrow(bim)
#
#   print(paste0(num_marker_excluded, " SNPs are discarded."))
#   print(paste0("Plink files with passed SNPs are in ", ResultDir, " prefixed as ", foutput))
#   print(paste0("Plink files with discarded SNPs are in ", ResultDir, " prefixed as ", foutput, "_snps_extracted"))
#
#   if (!regionfile) {
#     ftemp <- list.files(ResultDir, pattern = "region")
#     invisible(file.remove(file.path(ResultDir, ftemp)))
#   }
# }

## Function 10
######### Added in 3.0
outersect <- function(x, y, ...) {
  big.vec <- c(x, y, ...)
  duplicates <- big.vec[duplicated(big.vec)]
  setdiff(big.vec, unique(duplicates))
}

## Function 11
######### Added in 3.0
removeTempFiles <- function(directory, pattern) {

  file_path <- paste0(directory,"/sink_file.txt")

  # Create an empty file (if it doesn't exist)
  if (!file.exists(file_path)) {
    file.create(file_path)
  }

  # Redirect output to the specified file
  sink(file_path)

  tempFiles <- list.files(directory, pattern = pattern, full.names = TRUE)
  if (length(tempFiles) > 0) {
    file.remove(tempFiles)
  } else {
    print(paste("No files found with pattern:", pattern, "in directory:", directory))
  }

  # Reset the sink to stop redirecting the output to the file
  if (sink.number() > 0) {
    sink()  # This line resets the output redirection
  }

}

## Function 12
######### Added in 3.0
copyTempFiles <- function(sourceDir, targetDir, pattern) {
  tempFiles <- list.files(sourceDir, pattern = pattern, full.names = TRUE)
  if (length(tempFiles) > 0) {
    sapply(tempFiles, function(file) file.copy(file, targetDir))
  }
}

## Function 13
######### Added in 3.0
readDataFile <- function(filePath) {
  data <- read.table(filePath, sep = "", header = TRUE, as.is = TRUE, stringsAsFactors = FALSE)
  return(data)
}

## Function 14
######### Added in 3.0
createHeterozygosityPlot <- function(hetermiss, hetfail, imissfail, het, imiss, legend_text_size, legend_title_size, axis_text_size, axis_title_size, title_size) {
  # Data preparation
  hetermiss$logF_MISS <- log10(hetermiss$F_MISS)
  hetermiss$type <- "Passed"
  hetermiss$type[hetermiss$IID %in% hetfail$IID] <- "Failed H"
  hetermiss$type[hetermiss$IID %in% imissfail$IID] <- "Failed M"
  hetermiss$type[hetermiss$IID %in% intersect(hetfail$IID, imissfail$IID)] <- "Failed H & M"

  minus_sd <- mean(hetermiss$F) - 1:5 * stats::sd(hetermiss$F)
  plus_sd <- mean(hetermiss$F) + 1:5 * stats::sd(hetermiss$F)

  colors <- c("#666666", "#1b9e77", "#d95f02", "#7570b3")
  names(colors) <- c("Passed", "Failed H", "Failed M", "Failed H & M")

  hetermiss$shape <- "general"
  hetermiss$type <- factor(hetermiss$type, levels = names(colors))
  hetermiss$shape <- as.factor(hetermiss$shape)

  # Plotting
  print("Plots are initiated.")
  plot_hetimiss <- ggplot2::ggplot(data = hetermiss, ggplot2::aes(x = with(hetermiss,logF_MISS), y = F, color = with(hetermiss,type), shape = with(hetermiss,shape))) +
    ggplot2::geom_point() +
    ggplot2::scale_shape_manual(values = c(16, 17), guide = "none") +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(x = "Sample's missingness in log10", y = "Heterozygosity rate \n(and standard deviation)", color = "Sample Status", title = "Heterozygosity by missingness across samples") +
    ggplot2::geom_hline(yintercept = c(minus_sd[1:3], plus_sd[1:3]), lty = 2, col = "azure4") +
    ggplot2::scale_y_continuous(labels = c("-3", "-4", "-5", "+3", "+4", "+5"), breaks = c(minus_sd[3:5], plus_sd[3:5])) +
    ggplot2::scale_x_continuous(labels = c(0.0001, 0.001, 0.01, 0.03, 0.05, 0.01, 1), breaks = c(-4, -3, -2, log10(0.03), log10(0.05), -1, 0)) +
    ggplot2::geom_hline(yintercept = mean(hetermiss$F) - (het * stats::sd(hetermiss$F)), col = "#e7298a", lty = 2) +
    ggplot2::geom_hline(yintercept = mean(hetermiss$F) + (het * stats::sd(hetermiss$F)), col = "#e7298a", lty = 2) +
    ggplot2::geom_vline(xintercept = log10(imiss), col = "#e7298a", lty = 2) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.text = ggplot2::element_text(size = legend_text_size),
      legend.title = ggplot2::element_text(size = legend_title_size),
      axis.text = ggplot2::element_text(size = axis_text_size),
      axis.title = ggplot2::element_text(size = axis_title_size),
      title = ggplot2::element_text(size = title_size)
    )

  return(plot_hetimiss)
}

## Function 15
######### Added in 3.0
processAmbiguousSamples <- function(DataDir, ResultDir, finput, fam1) {
  pruneArgs <- c(
    "--bfile", paste0(DataDir, "/", finput),
    "--prune",
    "--must-have-sex",
    "--make-bed",
    "--out", paste0(ResultDir, "/", "FINPUT"),
    "--silent"
  )
  executePlink(pruneArgs)

  fam2 <- nrow(read.table(paste0(ResultDir,"/","FINPUT",".fam"),header = FALSE))
  print(paste0("No. of ambiguous samples filtered out: ", fam1 - fam2))

  copyTempFiles(ResultDir, DataDir, "FINPUT")
  "FINPUT"
}

## Function 16
######### Added in 3.0
filterSamples <- function(DataDir, ResultDir, finput, failed_het_imiss) {

  excludeSamplesArgs <- c(
    "--bfile", paste0(DataDir, "/", finput),
    "--remove", paste0(ResultDir, "/", "failed_het_imiss"),
    "--allow-no-sex",                                      ## Adding in 4.0
    "--make-bed",
    "--out", paste0(ResultDir, "/", "foutput"),
    "--silent"
  )
  executePlink(excludeSamplesArgs)
}

## Function 17
######### Added in 3.0
printSampleFilterResults <- function(imissfail, hetfail, failed_het_imiss) {
  if (nrow(imissfail) == 0) {
    print("No samples filtered for missingness.")
  } else {
    print(paste0(length(unique(imissfail$IID)), " samples filtered for missingness."))
  }
  if (nrow(hetfail) == 0) {
    print("No samples filtered for heterozygosity.")
  } else {
    print(paste0(length(unique(hetfail$IID)), " samples filtered for heterozygosity."))
  }
  if (nrow(failed_het_imiss) == 0) {
    print("No samples filtered for missingness and heterozygosity.")
  } else {
    print(paste0(length(unique(failed_het_imiss$IID)), " samples filtered for missingness and heterozygosity."))
  }
}

## Function 18
######### Added in 3.0
processIBDData <- function(IBD, IBDmatrix, ResultDir,foutput) {

  if (!is.null(IBD)) {
    # Compute and save filtered IBD data
    executePlinkForIBD(ResultDir,IBD, "filtered_ibd")

    # Optionally save the entire IBD matrix
    if (IBDmatrix) {
      executePlinkForIBD(ResultDir,NULL, "Entire_ibd")
      print("Entire IBD matrix 'Entire_ibd.genome' saved in ResultDir.")
    }

    # Read filtered IBD data
    ibd <- readIBDData(ResultDir, "filtered_ibd.genome")
    failed_ibd <- identifyFailedSamplesFromIBD(ibd, ResultDir)

    # Update PLINK files based on IBD results
    updatePlinkFilesWithIBDFilter(ResultDir = ResultDir, foutput = foutput, failed_ibd = failed_ibd)
  } else {
    # Generate BED files without IBD filtering
    executeMakeBed(ResultDir = ResultDir, foutput = foutput)
    NULL
  }
  return(list(failed_ibd = failed_ibd,ibd = ibd))
}

## Function 19
######### Added in 3.0
executePlinkForIBD <- function(ResultDir, IBD, outFileName) {
  ibdArgs <- c(
    "--bfile", paste0(ResultDir, "/", "foutput"),
    "--genome",
    "--out", paste0(ResultDir, "/", outFileName),
    "--silent"
  )
  if (!is.null(IBD)) {
    ibdArgs <- c(ibdArgs, "--min", IBD)
  }
  executePlink(ibdArgs)
}

## Function 20
######### Added in 3.0
readIBDData <- function(ResultDir, fileName) {
  data.table::fread(paste0(ResultDir, "/", fileName), select = c("IID1", "IID2", "PI_HAT"))
}

## Function 21
######### Added in 3.0
identifyFailedSamplesFromIBD <- function(ibd, ResultDir) {
  failedSamples <- unique(c(ibd$IID1, ibd$IID2))
  if (length(failedSamples) > 0) {
    famData <- read.table(paste0(ResultDir, "/", "foutput", ".fam"))
    famData[famData$V2 %in% failedSamples, 1:2]
  } else {
    NULL
  }
}

## Function 22
######### Added in 3.0
updatePlinkFilesWithIBDFilter <- function(ResultDir, foutput, failed_ibd) {
  if (!is.null(failed_ibd)) {
    write.table(failed_ibd, file = paste0(ResultDir, "/samples_failed_ibd"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    removeSamplesArgs <- c("--bed", paste0(ResultDir, "/", "foutput", ".bed"), "--bim", paste0(ResultDir, "/", "foutput", ".bim"), "--fam", paste0(ResultDir, "/", "foutput", ".fam"), "--remove", paste0(ResultDir, "/samples_failed_ibd"), "--allow-no-sex", "--make-bed", "--out", paste0(ResultDir, "/", foutput), "--silent")
    executePlink(removeSamplesArgs)
  }else{
    removeSamplesArgs <- c("--bed", paste0(ResultDir, "/", "foutput", ".bed"), "--bim", paste0(ResultDir, "/", "foutput", ".bim"), "--fam", paste0(ResultDir, "/", "foutput", ".fam"), "--make-bed", "--out", paste0(ResultDir, "/", foutput), "--silent")
    executePlink(removeSamplesArgs)

  }

}

## Function 23
######### Added in 3.0
executeMakeBed <- function(ResultDir, foutput) {
  makeBedArgs <- c("--bfile", paste0(ResultDir, "/", "foutput"), "--make-bed", "--out", paste0(ResultDir, "/", foutput), "--silent")
  executePlink(makeBedArgs)
}

## Function 24
######### Added in 3.0
# Adjusted Helper Function
executePlinkAd <- function(ResultDir, args) {
  tryCatch({
    stderr_dest <- ifelse(.Platform$OS.type == "windows", "NUL", "/dev/null")
    invisible(sys::exec_wait(file.path(ResultDir, "./plink"), args = args, std_err = stderr_dest))
  }, error = function(e) {
    stop("An error occurred while executing Plink: ", e$message)
  })
}

## Function 25
######### Added in 3.0
# Helper Function to Set PLINK Flags
setPlinkFlags <- function(maf, geno, hwe, hweCase, hweControl) {
  MAF <- if (!is.null(maf)) "--maf" else NULL
  GENO <- if (!is.null(geno)) "--geno" else NULL

  if (!is.null(hwe)) {
    if (!is.null(hweControl)) {
      print("Since hwe is not NULL, hweControl should be NULL. Setting hweControl = NULL implicitly.")
      hweControl <- NULL
    } else if (!is.null(hweCase)) {
      print("Since hwe is not NULL, hweCase should be NULL. Setting hweCase = NULL implicitly.")
      hweCase <- NULL
    }
    HWE <- "--hwe"
  } else {
    HWE <- NULL
  }

  if (is.null(hweCase) && !is.null(hweControl)) {
    stop("hweControl cannot be NULL if hweCase is not NULL.")
  } else if (!is.null(hweCase) && is.null(hweControl)) {
    stop("hweCase cannot be NULL if hweControl is not NULL.")
  }

  HWECase <- if (!is.null(hweCase)) "--hwe" else NULL
  HWECon <- if (!is.null(hweControl)) "--hwe" else NULL

  if (!is.null(hweCase) || !is.null(hweControl)) {
    HWE <- NULL
    hwe <- NULL
  }

  return(list(MAF = MAF, GENO = GENO, HWE = HWE, HWECase = HWECase, HWECon = HWECon))
}


## Function 26
######### Added in 3.0
# Helper Function to Remove Ambiguous SNPs
removeAmbiguousSNPs <- function(DataDir, ResultDir, finput) {
  bimFilePath <- paste0(DataDir, "/", finput, ".bim")
  study <- read.table(file = bimFilePath, stringsAsFactors = FALSE)

  # Identifying Ambiguous SNPs
  study_AT <- study[study[, 5] == "A" & study[, 6] == "T", 2, drop = FALSE]
  study_TA <- study[study[, 5] == "T" & study[, 6] == "A", 2, drop = FALSE]
  study_GC <- study[study[, 5] == "G" & study[, 6] == "C", 2, drop = FALSE]
  study_CG <- study[study[, 5] == "C" & study[, 6] == "G", 2, drop = FALSE]

  # Identifying Indels
  study_indel1 <- study[!which(study[, 5] != "A" & study[, 5] != "T" & study[, 5] != "G" & study[, 5] != "C"), 2, drop = FALSE]
  study_indel2 <- study[!which(study[, 6] != "A" & study[, 6] != "T" & study[, 6] != "G" & study[, 6] != "C"), 2, drop = FALSE]

  study_SNP <- rbind(study_AT, study_TA, study_GC, study_CG, study_indel1, study_indel2)

  # Write Ambiguous SNPs to a file
  outputFile <- paste0(ResultDir, "/study_SNP")
  write.table(study_SNP, file = outputFile, quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\r\n", sep = " ")

  # Remove Ambiguous SNPs using PLINK
  executePlinkAd(ResultDir, args = c("--bfile", paste0(DataDir, "/", finput), "--exclude", outputFile, "--allow-no-sex", "--make-bed", "--out", paste0(ResultDir, "/NoAmbiguousSNP"), "--silent"))

  # Return the count of removed SNPs
  return(nrow(study_SNP))
}

## Function 27
######### Added in 3.0
# Helper Function to Apply Filters with PLINK
applyFiltersWithPlink <- function(ResultDir, DataDir, finput, MAF, maf, GENO, geno, HWE, hwe) {
  # Execute PLINK command
  executePlinkAd(ResultDir, args = c(
    "--bfile", paste0(ResultDir, "/NoAmbiguousSNP"),
    MAF, maf,
    GENO, geno,
    HWE, hwe,
    "--allow-no-sex",
    "--make-bed",
    "--out", paste0(ResultDir, "/filtered_temp1"),
    "--silent"
  ))

  # Check if the file exists and print relevant messages
  if (file.exists(paste0(ResultDir, "/filtered_temp1.bed"))) {
    print("Thresholds for maf, geno and hwe worked.")
    logContents <- readLines(paste0(ResultDir, "/filtered_temp1.log"))
    print(grep("variants removed", logContents, value = TRUE))
  } else {
    print("Error applying thresholds or file not found.")
    executePlinkAd(ResultDir, args = c(
      "--bfile", paste0(DataDir, "/", finput),
      "--make-bed",
      "--out", paste0(ResultDir, "/filtered_temp1"),
      "--silent"
    ))
  }
}

## Function 28
######### Added in 3.0
# Helper Function to Apply HWE Filters and monomorpic part with PLINK
applyCaseControlFilters <- function(ResultDir, fam4, casecontrol, HWECase, hweCase, HWECon, hweControl) {
  nextFile <- "filtered_temp1" # Default next file

  if (length(unique(fam4$V6)) >= 2 && casecontrol) {
    # Filter for cases
    executePlinkAd(ResultDir, args = c(
      "--bfile", paste0(ResultDir, "/filtered_temp1"),
      "--filter-cases", HWECase, hweCase,
      "--allow-no-sex",
      "--make-bed", "--out", paste0(ResultDir, "/filtered_temp_hwe_case_filtered"),
      "--silent"
    ))
    printHWEMessages(ResultDir, "/filtered_temp_hwe_case_filtered.log", "In cases")

    # Filter for controls
    executePlinkAd(ResultDir, args = c(
      "--bfile", paste0(ResultDir, "/filtered_temp1"),
      "--filter-controls", HWECon, hweControl,
      "--allow-no-sex",
      "--make-bed", "--out", paste0(ResultDir, "/filtered_temp_hwe_control_filtered"),
      "--silent"
    ))
    printHWEMessages(ResultDir, "/filtered_temp_hwe_control_filtered.log", "In controls")

    # Merge case and control filtered files
    MergeRegion(DataDir = ResultDir, ResultDir, "filtered_temp_hwe_case_filtered", "filtered_temp_hwe_control_filtered", "filtered_temp2", use_common_snps = TRUE)

    nextFile <- "filtered_temp2"
  } else {
    if (length(unique(fam4$V6)) == 1) {
      print("There is no case-control status in the plink files. Setting casecontrol = FALSE implicitly.")
      casecontrol <- FALSE
    } else {
      casecontrol <- FALSE
    }
  }

  # Processing for monomorphic SNPs
  if (file.exists(paste0(ResultDir, "/", nextFile, ".bed"))) {
    executePlinkAd(ResultDir, args = c(
      "--bfile", paste0(ResultDir, "/", nextFile),
      "--freq", "--make-bed", "--allow-no-sex",
      "--out", paste0(ResultDir, "/filtered_temp4"),
      "--silent"
    ))
  } else {
    print("Something went wrong.")
    print(grep("Error", readLines(paste0(ResultDir, "/", nextFile, ".log")), value = TRUE))
  }

  return(casecontrol)
}

## Function 29
######### Added in 3.0
# Helper function to print HWE messages from log files
printHWEMessages <- function(ResultDir, logFilePath, context) {
  if (file.exists(paste0(ResultDir, logFilePath))) {
    logContents <- readLines(paste0(ResultDir, logFilePath))
    message <- paste0(context, ", ", stringr::str_sub(grep("Hardy-Weinberg", logContents, value = TRUE), 7))
    print(message)
  } else {
    print(paste0(context, " log file not found."))
  }
}

## Function 30
######### Added in 3.0
# Helper Function to Handle Monomorphic SNPs
handleMonomorphicSNPs <- function(monomorphicSNPs, mmSNPs, ResultDir) {
  if (monomorphicSNPs == TRUE && nrow(mmSNPs) != 0) {
    mmSNP1 <- mmSNPs$SNP
    exclude <- "--exclude"
    excludemono <- paste0(ResultDir, "/monomorphicSNPs")
  } else {
    mmSNP1 <- NULL
    exclude <- NULL
    excludemono <- NULL
    if (monomorphicSNPs == TRUE || nrow(mmSNPs) != 0) {
      print("There are no monomorphic SNPs.")
    }
  }

  return(list(mmSNP1 = mmSNP1, exclude = exclude, excludemono = excludemono))
}

## Function 31
######### Added in 3.0
# Helper Function to Handle LD Pruning
handleLDPruning <- function(ld_prunning, highLD_regions, ResultDir, window_size, step_size, r2_threshold) {
  if (ld_prunning) {
    excluderange <- "--exclude"

    if (!is.null(highLD_regions)) {
      # Write high LD regions to a temporary file
      tempFile <- paste0(ResultDir, "/highLD_regions_temp")
      write.table(highLD_regions, file = tempFile, quote = FALSE, row.names = FALSE, col.names = FALSE)
      highLD_regions <- tempFile
    } else {
      highLD_regions <- NULL
    }

    indep <- "--indep-pairwise"
  } else {
    excluderange <- NULL
    highLD_regions <- NULL
    indep <- NULL
    window_size <- NULL
    step_size <- NULL
    r2_threshold <- NULL
  }

  return(list(excluderange = excluderange, highLD_regions = highLD_regions, indep = indep, window_size = window_size, step_size = step_size, r2_threshold = r2_threshold))
}

## Function 32
######### Added in 3.0
# Helper Function to Execute PLINK with Various Parameters
executePlinkWithParams <- function(ResultDir, filtered_temp, exclude, excludemono, excluderange, highLD_regions, indep, window_size, step_size, r2_threshold) {
  executePlinkAd(ResultDir, args = c(
    "--bfile", paste0(ResultDir, "/", filtered_temp),
    exclude, excludemono, excluderange, highLD_regions,
    indep, window_size, step_size, r2_threshold,
    "--make-bed",
    "--allow-no-sex",
    "--out", paste0(ResultDir, "/", filtered_temp, "_processed"),
    "--silent"
  ))
}

## Function 33
######### Added in 3.0
# Helper Function to Handle Differential Missingness Filtering for Case-Control Data
handleCaseControlFiltering <- function(ResultDir, casecontrol, dmissX, dmissAutoY, caldiffmiss, SNPmissCC, diffmissFilter, foutput) {

  SNPmissCC <- NULL  # Initialize

  if (casecontrol) {
    chrfilter <- NULL
    chrv <- NULL
    if (dmissX & dmissAutoY) {
      # No additional filters needed
    } else if (dmissX & !dmissAutoY) {
      chrfilter <- "--chr"
      chrv <- 23
    } else if (!dmissX & dmissAutoY) {
      chrfilter <- "--not-chr"
      chrv <- 23
    } else {
      print("Filtering for differential missingness between cases and controls is turned off.")
    }

    if (caldiffmiss) {
      executePlinkAd(ResultDir, args = c(
        "--bfile", paste0(ResultDir, "/filtered_temp4"),
        chrfilter, chrv,
        "--test-missing", "--adjust",
        "--make-bed", "--allow-no-sex",
        "--out", paste0(ResultDir, "/filtered_temp_casecontrol"),
        "--silent"
      ))

      # Process the differential missingness results
      SNPmissCC <- processDifferentialMissingnessResults(ResultDir)
    }

    applySNPmissCCFilter(ResultDir, SNPmissCC, diffmissFilter, foutput)

  } else {
    print("No filter based on differential missingness will be applied.")

    executePlinkAd(ResultDir, args = c(
      "--bfile", paste0(ResultDir, "/filtered_temp4"),
      "--make-bed", "--allow-no-sex",
      "--out", paste0(ResultDir, "/", foutput),
      "--silent"
    ))
  }

  return( SNPmissCC)
}

## Function 34
######### Added in 3.0
# Helper function to process differential missingness results
processDifferentialMissingnessResults <- function(ResultDir) {
  filePath <- paste0(ResultDir, "/filtered_temp_casecontrol.missing.adjusted")
  if (file.exists(filePath)) {
    ccmissing <- read.table(filePath, header = TRUE, sep = "")
    diffmiss <- 0.05 / length(unique(ccmissing$SNP))
    x <- ccmissing[ccmissing$BONF < diffmiss, "SNP", drop = TRUE]
    return(x)
  } else {
    return(NULL)
  }
}

## Function 35
######### Added in 3.0
# Helper function to apply SNP missingness filter
applySNPmissCCFilter <- function(ResultDir, SNPmissCC, diffmissFilter, foutput) {

  if (!is.null(SNPmissCC) && diffmissFilter) {
    write.table(
      SNPmissCC,
      file = paste0(ResultDir, "/SNPdifCallrate"),
      quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\r\n", sep = " "
    )
    executePlinkAd(ResultDir, args = c(
      "--bfile", paste0(ResultDir, "/filtered_temp4"),
      "--exclude", paste0(ResultDir, "/SNPdifCallrate"),
      "--allow-no-sex",
      "--make-bed",
      "--out", paste0(ResultDir, "/", foutput),
      "--silent"
    ))
  } else {
    print("No SNP with differential missingness between cases and controls.")
    executePlinkAd(ResultDir, args = c(
      "--bfile", paste0(ResultDir, "/filtered_temp4"),
      "--allow-no-sex",
      "--make-bed",
      "--out", paste0(ResultDir, "/", foutput),
      "--silent"
    ))
  }
}

## Function 36
######### Added in 3.0
#### gmirror function from hudson R package
# gmirror <- function(top, bottom, tline, bline, chroms = c(1:22, "X", "Y"),log10=TRUE,
#                     yaxis, opacity=1, annotate_snp, annotate_p, toptitle=NULL,
#                     bottomtitle=NULL, highlight_snp, highlight_p, highlighter="red",
#                     chrcolor1="#AAAAAA", chrcolor2="#4D4D4D", freey=FALSE,
#                     background="variegated", chrblocks=FALSE, file="gmirror",
#                     type="png", hgt=7, hgtratio=0.5, wi=12, res=300 ){
#
#   #Sort data
#   topn <- names(top)
#   bottomn <- names(bottom)
#   top$Location <- "Top"
#   bottom$Location <- "Bottom"
#
#   # Check file formats
#   if(!identical(topn, bottomn)){stop("Please ensure both inputs have the same metadata columns.")}
#
#   d <- as.data.frame(rbind(top, bottom))
#
#   d$POS <- as.numeric(as.character(d$POS))
#   d$CHR <- droplevels(factor(d$CHR, levels = as.character(chroms)))
#   d <- d[d$CHR %in% chroms, ]
#   d_order <- d[order(d$CHR, d$POS), ]
#   d_order$pos_index <- seq.int(nrow(d_order))
#   d_order_sub <- d_order[, c("SNP", "CHR", "POS", "pvalue", "pos_index")]
#
#   #Set up dataframe with color and position info
#   maxRows <- by(d_order_sub, d_order_sub$CHR, function(x) x[which.max(x$pos_index),])
#   minRows <- by(d_order_sub, d_order_sub$CHR, function(x) x[which.min(x$pos_index),])
#   milimits <- do.call(rbind, minRows)
#   malimits <- do.call(rbind, maxRows)
#   lims <- merge(milimits, malimits, by="CHR")
#   names(lims) <- c("Color", "snpx", "px", "posx", "posmin", "snpy", "py", "posy", "posmax")
#   lims$av <- (lims$posmin + lims$posmax)/2
#   lims <- lims[order(lims$Color),]
#   lims$shademap <- rep(c("shade_ffffff", "shade_ebebeb"), length.out=nrow(lims), each=1)
#
#   #Set up colors
#   nchrcolors <- nlevels(factor(lims$Color))
#
#   #Color by CHR
#   colnames(d_order)[2] <- "Color"
#   newcols <-c(rep(x=c(chrcolor1, chrcolor2), length.out=nchrcolors, each=1), "#FFFFFF", "#EBEBEB")
#   names(newcols) <-c(levels(factor(lims$Color)), "shade_ffffff", "shade_ebebeb")
#
#   #Info for y-axis
#   if(log10==TRUE){
#     d_order$pval <- -log10(d_order$pvalue)
#     yaxislab1 <- expression(paste("-log"[10], "(p-value)", sep=""))
#     yaxislab2 <- expression(paste("-log"[10], "(p-value)", sep=""))
#     if(!missing(tline)) {tredline <- -log10(tline)}
#     if(!missing(bline)) {bredline <- -log10(bline)}
#   } else {
#     d_order$pval <- d_order$pvalue
#     yaxislab1 <- yaxis[1]
#     yaxislab2 <- yaxis[2]
#     if(!missing(tline)) {tredline <- tline}
#     if(!missing(bline)) {bredline <- bline}
#   }
#   yaxismax1 <- ifelse(freey==FALSE, max(d_order$pval[which(d_order$pval< Inf)]), max(d_order$pval[which(d_order$pval< Inf) & d_order$Location=="Top"]))
#   yaxismax2 <- ifelse(freey==FALSE, max(d_order$pval[which(d_order$pval< Inf)]), max(d_order$pval[which(d_order$pval< Inf) & d_order$Location=="Bottom"]))
#   yaxismin1 <- ifelse(freey==FALSE, 0, min(d_order$pval[d_order$Location=="Top"]))
#   yaxismin2 <- ifelse(freey==FALSE, 0, min(d_order$pval[d_order$Location=="Bottom"]))
#
#   #Theme options
#   backpanel1 <- ifelse(background=="white", "NULL", "ggplot2::geom_rect(data = lims, ggplot2::aes(xmin = posmin-.5, xmax = posmax+.5, ymin = yaxismin1, ymax = Inf, fill=factor(shademap)), alpha = 0.5)" )
#   backpanel2 <- ifelse(background=="white", "NULL", "ggplot2::geom_rect(data = lims, ggplot2::aes(xmin = posmin-.5, xmax = posmax+.5, ymin = yaxismin2, ymax = Inf, fill=factor(shademap)), alpha = 0.5)" )
#
#   #Start plotting
#   #TOP PLOT
#   p1 <- ggplot2::ggplot() + eval(parse(text=backpanel1))
#   #Add shape info if available
#   if("Shape" %in% topn){
#     d1 <- d_order[d_order$Location=="Top",]
#     p1 <- p1 + ggplot2::geom_point(data=d1, ggplot2::aes(x=d1$pos_index, y=d1$pval, color=factor(d1$Color), shape=factor(d1$Shape)), alpha=opacity)
#   } else {
#     d1 <- d_order[d_order$Location=="Top",]
#     p1 <- p1 + ggplot2::geom_point(data=d1, ggplot2::aes(x=d1$pos_index, y=d1$pval, color=factor(d1$Color)), alpha=opacity)
#   }
#   p1 <- p1 + ggplot2::scale_x_continuous(breaks=lims$av, labels=lims$Color, expand=c(0,0))
#   if(chrblocks==TRUE){
#     p1 <- p1 + ggplot2::geom_rect(data = lims, ggplot2::aes(xmin = lims$posmin-.5, xmax = lims$posmax+.5, ymin = -Inf, ymax = min(d_order$pval), fill=as.factor(Color)), alpha = 1)
#   }
#   p1 <- p1 + ggplot2::scale_colour_manual(name = "Color", values = newcols) + ggplot2::scale_fill_manual(name = "Color", values = newcols)
#   p1 <- p1 + ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(), axis.title.x=ggplot2::element_blank(), legend.position="top", legend.title = ggplot2::element_blank())
#
#   #BOTTOM PLOT
#   p2 <- ggplot2::ggplot() + eval(parse(text=backpanel2))
#   #Add shape info if available
#   if("Shape" %in% bottomn){
#     d1 <- d_order[d_order$Location=="Bottom",]
#     p2 <- p2 + ggplot2::geom_point(data=d1, ggplot2::aes(x=d1$pos_index, y=d1$pval, color=factor(d1$Color), shape=factor(d1$Shape)), alpha=opacity)
#   } else {
#     d1 <- d_order[d_order$Location=="Bottom",]
#     p2 <- p2 + ggplot2::geom_point(data=d1, ggplot2::aes(x=d1, color=factor(d1$Color)), alpha=opacity)
#   }
#   p2 <- p2 + ggplot2::scale_x_continuous(breaks=lims$av, labels=lims$Color, expand=c(0,0))
#   if(chrblocks==TRUE){
#     p2 <- p2 + ggplot2::geom_rect(data = lims, ggplot2::aes(xmin = lims$posmin-.5, xmax = lims$posmax+.5, ymin = -Inf, ymax = min(d_order$pval), fill=as.factor(lims$Color)), alpha = 1)
#   }
#   p2 <- p2 + ggplot2::scale_colour_manual(name = "Color", values = newcols) + ggplot2::scale_fill_manual(name = "Color", values = newcols)
#   p2 <- p2 + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90), panel.grid.minor.x = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank(), legend.position="bottom", legend.title = ggplot2::element_blank())
#
#   #Highlight if given
#   if(!missing(highlight_snp)){
#     if("Shape" %in% topn){
#       d1 <- d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Top"]
#       p1 <- p1 + ggplot2::geom_point(data= d1, ggplot2::aes(x=d1$pos_index, y=d1$pval, shape=d1$Shape), colour=highlighter)
#       p1 <- p1 + ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(colour = "black")))
#     } else {
#       d1 <- d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Top", ]
#       p1 <- p1 + ggplot2::geom_point(data= d1, ggplot2::aes(x=d1$pos_index, y=d1$pval), colour=highlighter)
#     }
#     if("Shape" %in% bottomn){
#       d1 <- d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Bottom", ]
#       p2 <- p2 + ggplot2::geom_point(data=d1, ggplot2::aes(x=d1$pos_index, y=d1$pval, shape=d1$Shape), colour=highlighter)
#       p2 <- p2 + ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(colour = "black")))
#     } else {
#       d1 <- d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Bottom", ]
#       p2 <- p2 + ggplot2::geom_point(data=d1, ggplot2::aes(x=d1$pos_index, y=d1$pval), colour=highlighter)
#     }
#   }
#   if(!missing(highlight_p)){
#     if("Shape" %in% topn){
#       d1 <- d_order[d_order$pvalue < highlight_p[1] & d_order$Location=="Top", ]
#       p1 <- p1 + ggplot2::geom_point(data=d1, ggplot2::aes(x=d1$pos_index, y=d1$pval, shape=d1$Shape), colour=highlighter)
#       p1 <- p1 + ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(colour = "black")))
#     } else {
#       d1 <- d_order[d_order$pvalue < highlight_p[1] & d_order$Location=="Top", ]
#       p1 <- p1 + ggplot2::geom_point(data=d1, ggplot2::aes(x=d1$pos_index, y=d1$pval), colour=highlighter)
#     }
#     if("Shape" %in% bottomn){
#       d1 <- d_order[d_order$pvalue < highlight_p[2] & d_order$Location=="Bottom", ]
#       p2 <- p2 + ggplot2::geom_point(data=d1, ggplot2::aes(x=d1$pos_index, y=d1$pval, shape=d1$Shape), colour=highlighter)
#       p2 <- p2 + ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(colour = "black")))
#     } else {
#       d1 <- d_order[d_order$pvalue < highlight_p[2] & d_order$Location=="Bottom", ]
#       p2 <- p2 + ggplot2::geom_point(data=d1, ggplot2::aes(x=d1$pos_index, y=d1$pval), colour=highlighter)
#     }
#   }
#   #Add pvalue threshold line
#   if(!missing(tline)){
#     for(i in 1:length(tline)){
#       p1 <- p1 + ggplot2::geom_hline(yintercept = tredline[i], colour="red")
#     }
#   }
#   if(!missing(bline)){
#     for(i in 1:length(bline)){
#       p2 <- p2 + ggplot2::geom_hline(yintercept = bredline[i], colour="red")
#     }
#   }
#   #Annotate
#   if(!missing(annotate_p)){
#     if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE) {
#       print("Consider installing 'ggrepel' for improved text annotation")
#       d1 <- d_order[d_order$pvalue < annotate_p[1] & d_order$Location=="Top",]
#       p1 <- p1 + ggplot2::geom_text(data=d1, ggplot2::aes(d1$pos_index,d1$pval,label=d1$SNP))
#       p2 <- p2 + ggplot2::geom_text(data=d_order[d_order$pvalue < annotate_p[2] & d_order$Location=="Bottom",], ggplot2::aes(pos_index,pval,label=SNP))
#     } else {
#       d1 <- d_order[d_order$pvalue < annotate_p[1] & d_order$Location=="Top",]
#       p1 <- p1 + ggrepel::geom_text_repel(data=d1, ggplot2::aes(d1$pos_index,d1$pval,label=d1$SNP))
#       d2 <- d_order[d_order$pvalue < annotate_p[2] & d_order$Location=="Bottom",]
#       p2 <- p2 + ggrepel::geom_text_repel(data=d2, ggplot2::aes(d2$pos_index,d2$pval,label=d2$SNP))
#     }
#   }
#   if(!missing(annotate_snp)){
#     if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE){
#       print("Consider installing 'ggrepel' for improved text annotation")
#       d1 <- d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Top",]
#       p1 <- p1 + ggplot2::geom_text(data=d1, ggplot2::aes(d1$pos_index,d1$pval,label=d1$SNP))
#       d2 <- d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Bottom",]
#       p2 <- p2 + ggplot2::geom_text(data=d2, ggplot2::aes(d2$pos_index,d2$pval,label=d2$SNP))
#     } else {
#       d1 <- d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Top",]
#       p1 <- p1 + ggrepel::geom_text_repel(data=d1, ggplot2::aes(d1$pos_index,d1$pval,label=d1$SNP))
#       d2 <- d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Bottom",]
#       p2 <- p2 + ggrepel::geom_text_repel(data=d2, ggplot2::aes(d2$pos_index,d2$pval,label=d2$SNP))
#     }
#   }
#   #Add title and y axis title
#   p1 <- p1 + ggplot2::ylab(yaxislab1)
#   p2 <- p2 + ggplot2::ylab(yaxislab2)
#
#   #Format
#   if(chrblocks==TRUE){
#     if(freey==TRUE){
#       print("Sorry, drawing chrblocks with freey=TRUE is currently unsupported and will be ignored.")
#     } else {
#       p1 <- p1+ggplot2::theme(axis.text.x = ggplot2::element_text(vjust=1),axis.ticks.x = ggplot2::element_blank()) + ggplot2::ylim(c(yaxismin1,yaxismax1))
#       p2 <- p2+ggplot2::scale_y_reverse(limits=c(yaxismax2, yaxismin2)) + ggplot2::theme(axis.text.x = ggplot2::element_blank(),axis.ticks.x = ggplot2::element_blank())
#     }
#   } else {
#     p1 <- p1+ggplot2::theme(axis.text.x = ggplot2::element_text(vjust=1),axis.ticks.x = ggplot2::element_blank())+ ggplot2::scale_y_continuous(limits=c(yaxismin1, yaxismax1),expand=expansion(mult=c(0,0.1)))
#     p2 <- p2+ggplot2::scale_y_reverse(limits=c(yaxismax2,yaxismin2), expand=expansion(mult=c(0.1,0))) + ggplot2::theme(axis.text.x = ggplot2::element_blank(),axis.ticks.x = ggplot2::element_blank())
#   }
#
#   if(background=="white"){
#     p1 <- p1 + ggplot2::theme(panel.background = element_rect(fill="white"))
#     p2 <- p2 + ggplot2::theme(panel.background = element_rect(fill="white"))
#   }
#   p1 <- p1 + ggplot2::guides(fill="none", color="none")
#   p2 <- p2 + ggplot2::guides(fill="none", color="none")
#   #Save
#   print(paste0("Saving plot to ", file, ".", type))
#   p <- gridExtra::grid.arrange(gridExtra::arrangeGrob(p1, top=toptitle), gridExtra::arrangeGrob(p2, bottom=bottomtitle), padding=0, heights=c(hgtratio,1-hgtratio))
#   ggplot2::ggsave(p, filename=paste0(file, ".", type), dpi=res, units="in", height=hgt, width=wi)
#   return(p)
# }
gmirror <- function(top, bottom, tline, bline, chroms = c(1:22, "X", "Y"),log10=TRUE,
                    yaxis, opacity=1, annotate_snp, annotate_p, toptitle=NULL,
                    bottomtitle=NULL, highlight_snp, highlight_p, highlighter="red",
                    chrcolor1="#AAAAAA", chrcolor2="#4D4D4D", freey=FALSE,
                    background="variegated", chrblocks=FALSE, file="gmirror",
                    type="png", hgt=7, hgtratio=0.5, wi=12, res=300 ){

  #Sort data
  topn <- names(top)
  bottomn <- names(bottom)
  top$Location <- "Top"
  bottom$Location <- "Bottom"

  # Check file formats
  if(!identical(topn, bottomn)){stop("Please ensure both inputs have the same metadata columns.")}

  d <- as.data.frame(rbind(top, bottom))

  d$POS <- as.numeric(as.character(d$POS))
  d$CHR <- droplevels(factor(d$CHR, levels = as.character(chroms)))
  d <- d[d$CHR %in% chroms, ]
  d_order <- d[order(d$CHR, d$POS), ]
  d_order$pos_index <- seq.int(nrow(d_order))
  d_order_sub <- d_order[, c("SNP", "CHR", "POS", "pvalue", "pos_index")]

  #Set up dataframe with color and position info
  maxRows <- by(d_order_sub, d_order_sub$CHR, function(x) x[which.max(x$pos_index),])
  minRows <- by(d_order_sub, d_order_sub$CHR, function(x) x[which.min(x$pos_index),])
  milimits <- do.call(rbind, minRows)
  malimits <- do.call(rbind, maxRows)
  lims <- merge(milimits, malimits, by="CHR")
  names(lims) <- c("Color", "snpx", "px", "posx", "posmin", "snpy", "py", "posy", "posmax")
  lims$av <- (lims$posmin + lims$posmax)/2
  lims <- lims[order(lims$Color),]
  lims$shademap <- rep(c("shade_ffffff", "shade_ebebeb"), length.out=nrow(lims), each=1)

  #Set up colors
  nchrcolors <- nlevels(factor(lims$Color))

  #Color by CHR
  colnames(d_order)[2] <- "Color"
  newcols <-c(rep(x=c(chrcolor1, chrcolor2), length.out=nchrcolors, each=1), "#FFFFFF", "#EBEBEB")
  names(newcols) <-c(levels(factor(lims$Color)), "shade_ffffff", "shade_ebebeb")

  #Info for y-axis
  if(log10==TRUE){
    d_order$pval <- -log10(d_order$pvalue)
    yaxislab1 <- expression(paste("-log"[10], "(p-value)", sep=""))
    yaxislab2 <- expression(paste("-log"[10], "(p-value)", sep=""))
    if(!missing(tline)) {tredline <- -log10(tline)}
    if(!missing(bline)) {bredline <- -log10(bline)}
  } else {
    d_order$pval <- d_order$pvalue
    yaxislab1 <- yaxis[1]
    yaxislab2 <- yaxis[2]
    if(!missing(tline)) {tredline <- tline}
    if(!missing(bline)) {bredline <- bline}
  }
  yaxismax1 <- ifelse(freey==FALSE, max(d_order$pval[which(d_order$pval< Inf)]), max(d_order$pval[which(d_order$pval< Inf) & d_order$Location=="Top"]))
  yaxismax2 <- ifelse(freey==FALSE, max(d_order$pval[which(d_order$pval< Inf)]), max(d_order$pval[which(d_order$pval< Inf) & d_order$Location=="Bottom"]))
  yaxismin1 <- ifelse(freey==FALSE, 0, min(d_order$pval[d_order$Location=="Top"]))
  yaxismin2 <- ifelse(freey==FALSE, 0, min(d_order$pval[d_order$Location=="Bottom"]))

  #Theme options
  backpanel1 <- ifelse(background=="white", "NULL", "ggplot2::geom_rect(data = lims, ggplot2::aes(xmin = posmin-.5, xmax = posmax+.5, ymin = yaxismin1, ymax = Inf, fill=factor(shademap)), alpha = 0.5)" )
  backpanel2 <- ifelse(background=="white", "NULL", "ggplot2::geom_rect(data = lims, ggplot2::aes(xmin = posmin-.5, xmax = posmax+.5, ymin = yaxismin2, ymax = Inf, fill=factor(shademap)), alpha = 0.5)" )

  #Start plotting
  #TOP PLOT
  p1 <- ggplot2::ggplot() + eval(parse(text=backpanel1))
  #Add shape info if available
  if("Shape" %in% topn){
    p1 <- p1 + ggplot2::geom_point(data=d_order[d_order$Location=="Top",], ggplot2::aes(x=pos_index, y=pval, color=factor(Color), shape=factor(Shape)), alpha=opacity)
  } else {
    p1 <- p1 + ggplot2::geom_point(data=d_order[d_order$Location=="Top",], ggplot2::aes(x=pos_index, y=pval, color=factor(Color)), alpha=opacity)
  }
  p1 <- p1 + ggplot2::scale_x_continuous(breaks=lims$av, labels=lims$Color, expand=c(0,0))
  if(chrblocks==TRUE){
    p1 <- p1 + ggplot2::geom_rect(data = lims, ggplot2::aes(xmin = posmin-.5, xmax = posmax+.5, ymin = -Inf, ymax = min(d_order$pval), fill=as.factor(Color)), alpha = 1)
  }
  p1 <- p1 + ggplot2::scale_colour_manual(name = "Color", values = newcols) + ggplot2::scale_fill_manual(name = "Color", values = newcols)
  p1 <- p1 + ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(), axis.title.x=ggplot2::element_blank(), legend.position="top", legend.title = ggplot2::element_blank())

  #BOTTOM PLOT
  p2 <- ggplot2::ggplot() + eval(parse(text=backpanel2))
  #Add shape info if available
  if("Shape" %in% bottomn){
    p2 <- p2 + ggplot2::geom_point(data=d_order[d_order$Location=="Bottom",], ggplot2::aes(x=pos_index, y=pval, color=factor(Color), shape=factor(Shape)), alpha=opacity)
  } else {
    p2 <- p2 + ggplot2::geom_point(data=d_order[d_order$Location=="Bottom",], ggplot2::aes(x=pos_index, y=pval, color=factor(Color)), alpha=opacity)
  }
  p2 <- p2 + ggplot2::scale_x_continuous(breaks=lims$av, labels=lims$Color, expand=c(0,0))
  if(chrblocks==TRUE){
    p2 <- p2 + ggplot2::geom_rect(data = lims, ggplot2::aes(xmin = posmin-.5, xmax = posmax+.5, ymin = -Inf, ymax = min(d_order$pval), fill=as.factor(Color)), alpha = 1)
  }
  p2 <- p2 + ggplot2::scale_colour_manual(name = "Color", values = newcols) + ggplot2::scale_fill_manual(name = "Color", values = newcols)
  p2 <- p2 + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90), panel.grid.minor.x = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank(), legend.position="bottom", legend.title = ggplot2::element_blank())

  #Highlight if given
  if(!missing(highlight_snp)){
    if("Shape" %in% topn){
      p1 <- p1 + ggplot2::geom_point(data=d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Top", ], ggplot2::aes(x=pos_index, y=pval, shape=Shape), colour=highlighter)
      p1 <- p1 + ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(colour = "black")))
    } else {
      p1 <- p1 + ggplot2::geom_point(data=d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Top", ], ggplot2::aes(x=pos_index, y=pval), colour=highlighter)
    }
    if("Shape" %in% bottomn){
      p2 <- p2 + ggplot2::geom_point(data=d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Bottom", ], ggplot2::aes(x=pos_index, y=pval, shape=Shape), colour=highlighter)
      p2 <- p2 + ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(colour = "black")))
    } else {
      p2 <- p2 + ggplot2::geom_point(data=d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Bottom", ], ggplot2::aes(x=pos_index, y=pval), colour=highlighter)
    }
  }
  if(!missing(highlight_p)){
    if("Shape" %in% topn){
      p1 <- p1 + ggplot2::geom_point(data=d_order[d_order$pvalue < highlight_p[1] & d_order$Location=="Top", ], ggplot2::aes(x=pos_index, y=pval, shape=Shape), colour=highlighter)
      p1 <- p1 + ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(colour = "black")))
    } else {
      p1 <- p1 + ggplot2::geom_point(data=d_order[d_order$pvalue < highlight_p[1] & d_order$Location=="Top", ], ggplot2::aes(x=pos_index, y=pval), colour=highlighter)
    }
    if("Shape" %in% bottomn){
      p2 <- p2 + ggplot2::geom_point(data=d_order[d_order$pvalue < highlight_p[2] & d_order$Location=="Bottom", ], ggplot2::aes(x=pos_index, y=pval, shape=Shape), colour=highlighter)
      p2 <- p2 + ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(colour = "black")))
    } else {
      p2 <- p2 + ggplot2::geom_point(data=d_order[d_order$pvalue < highlight_p[2] & d_order$Location=="Bottom", ], ggplot2::aes(x=pos_index, y=pval), colour=highlighter)
    }
  }
  #Add pvalue threshold line
  if(!missing(tline)){
    for(i in 1:length(tline)){
      p1 <- p1 + ggplot2::geom_hline(yintercept = tredline[i], colour="red")
    }
  }
  if(!missing(bline)){
    for(i in 1:length(bline)){
      p2 <- p2 + ggplot2::geom_hline(yintercept = bredline[i], colour="red")
    }
  }
  #Annotate
  if(!missing(annotate_p)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE) {
      print("Consider installing 'ggrepel' for improved text annotation")
      p1 <- p1 + ggplot2::geom_text(data=d_order[d_order$pvalue < annotate_p[1] & d_order$Location=="Top",], ggplot2::aes(pos_index,pval,label=SNP))
      p2 <- p2 + ggplot2::geom_text(data=d_order[d_order$pvalue < annotate_p[2] & d_order$Location=="Bottom",], ggplot2::aes(pos_index,pval,label=SNP))
    } else {
      p1 <- p1 + ggrepel::geom_text_repel(data=d_order[d_order$pvalue < annotate_p[1] & d_order$Location=="Top",], ggplot2::aes(pos_index,pval,label=SNP))
      p2 <- p2 + ggrepel::geom_text_repel(data=d_order[d_order$pvalue < annotate_p[2] & d_order$Location=="Bottom",], ggplot2::aes(pos_index,pval,label=SNP))
    }
  }
  if(!missing(annotate_snp)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE){
      print("Consider installing 'ggrepel' for improved text annotation")
      p1 <- p1 + ggplot2::geom_text(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Top",], ggplot2::aes(pos_index,pval,label=SNP))
      p2 <- p2 + ggplot2::geom_text(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Bottom",], ggplot2::aes(pos_index,pval,label=SNP))
    } else {
      p1 <- p1 + ggrepel::geom_text_repel(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Top",], ggplot2::aes(pos_index,pval,label=SNP))
      p2 <- p2 + ggrepel::geom_text_repel(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Bottom",], ggplot2::aes(pos_index,pval,label=SNP))
    }
  }
  #Add title and y axis title
  p1 <- p1 + ggplot2::ylab(yaxislab1)
  p2 <- p2 + ggplot2::ylab(yaxislab2)

  #Format
  if(chrblocks==TRUE){
    if(freey==TRUE){
      print("Sorry, drawing chrblocks with freey=TRUE is currently unsupported and will be ignored.")
    } else {
      p1 <- p1+ggplot2::theme(axis.text.x = ggplot2::element_text(vjust=1),axis.ticks.x = ggplot2::element_blank()) + ggplot2::ylim(c(yaxismin1,yaxismax1))
      p2 <- p2+ggplot2::scale_y_reverse(limits=c(yaxismax2, yaxismin2)) + ggplot2::theme(axis.text.x = ggplot2::element_blank(),axis.ticks.x = ggplot2::element_blank())
    }
  } else {
    p1 <- p1+ggplot2::theme(axis.text.x = ggplot2::element_text(vjust=1),axis.ticks.x = ggplot2::element_blank())+ ggplot2::scale_y_continuous(limits=c(yaxismin1, yaxismax1),expand=expansion(mult=c(0,0.1)))
    p2 <- p2+ggplot2::scale_y_reverse(limits=c(yaxismax2,yaxismin2), expand=expansion(mult=c(0.1,0))) + ggplot2::theme(axis.text.x = ggplot2::element_blank(),axis.ticks.x = ggplot2::element_blank())
  }

  if(background=="white"){
    p1 <- p1 + ggplot2::theme(panel.background = element_rect(fill="white"))
    p2 <- p2 + ggplot2::theme(panel.background = element_rect(fill="white"))
  }
  p1 <- p1 + ggplot2::guides(fill="none", color="none")
  p2 <- p2 + ggplot2::guides(fill="none", color="none")
  #Save
  print(paste0("Saving plot to ", file, ".", type))
  p <- gridExtra::grid.arrange(gridExtra::arrangeGrob(p1, top=toptitle), gridExtra::arrangeGrob(p2, bottom=bottomtitle), padding=0, heights=c(hgtratio,1-hgtratio))
  ggplot2::ggsave(p, filename=paste0(file, ".", type), dpi=res, units="in", height=hgt, width=wi)
  return(p)
}
## Function 37
######### Added in 3.0
performLDClumping <- function(ldclump, DataDir, ResultDir, LDreference, summarystat, clump_p1, clump_p2, clump_r2, clump_kb, byCHR) {
  if (ldclump) {
    clumpedResult <- ClumpLD(DataDir = DataDir, ResultDir = ResultDir, finput = LDreference, SNPdata = list(summarystat),
                             clump_p1 = clump_p1, clump_p2 = clump_p2, clump_r2 = clump_r2, clump_kb = clump_kb, byCHR = byCHR)

    write.table(clumpedResult$SNP, file = paste0(ResultDir,"/","Valid.SNP"), quote = FALSE, row.names = FALSE)

    return(list(clumpExtract = "--extract", clumpSNP = paste0(ResultDir,"/","Valid.SNP")))
  } else {
    return(list(clumpExtract = NULL, clumpSNP = NULL))
  }
}

## Function 38
######### Added in 3.0
preparePhenotypeData <- function(phenofile, nPC, DataDir, ResultDir, finput, highLD_regions, ld_prunning, window_size, step_size, r2_threshold) {
  phenotype <- cbind(phenofile[,1:2], phenofile[,"Pheno1"])
  colnames(phenotype) <- c("FID","IID","Pheno1")

  if (nPC > 0) {
    GP <- ComputeGeneticPC(DataDir = DataDir, ResultDir = ResultDir, finput = finput, countPC = nPC, highLD_regions = highLD_regions,
                           ld_prunning = ld_prunning, window_size = window_size, step_size = step_size, r2_threshold = r2_threshold, plotPC = FALSE)
    colnames(GP) <- c("FID", "IID", paste0("PC",1:nPC))
    return(GP)
  } else {
    print("Parameter 'nPC' is either zero or negative. Genetic PC will not be computed.")
    return(phenotype[,1:2])
  }
}

## Function 39
######### Added in 3.0
computeNullModel <- function(pheno, pheno_type) {
  if (pheno_type == "binary") {
    # Convert binary phenotype to 0 and 1
    pheno$Pheno1 <- as.integer(pheno$Pheno1 == 2)
    # Fit logistic regression model
    null_model <- glm(Pheno1 ~ ., data = pheno[, !colnames(pheno) %in% c("FID", "IID")], family = binomial)
  } else {
    # Fit linear regression model
    null_model <- lm(Pheno1 ~ ., data = pheno[, !colnames(pheno) %in% c("FID", "IID")])
  }

  # Calculate R-squared value
  null_r2 <- summary(null_model)$r.squared

  return(list(model = null_model, r_squared = null_r2))
}

## Function 40
######### Added in 3.0
prsFun <- function(pthreshold, ResultDir, DataDir, finput, clumpExtract, clumpSNP, pheno, pheno_type, null_model) {

  print(pthreshold)
  print(paste0("Computing PRS for threshold ", pthreshold))

  pt <- data.table::as.data.table(cbind(pthreshold, 0, pthreshold))
  colnames(pt) <- c("Threshold", "Lowerbound", "UpperBound")
  write.table(pt, file = paste0(ResultDir, "/range_list"), quote = FALSE, row.names = FALSE)

  ######
  # x1 <-  read.table(paste0(DataDir, "/", finput,".bim"))
  # x2 <- read.table(paste0(ResultDir, "/","prssummarystat"),header = TRUE)
  # x3 <- merge(x2,x1,by.x = "SNP",by.y ="V2")
  # x4 <- x3[,1:ncol(x2)]
  # write.table(paste0(ResultDir, "/","prssummarystat"),quote=FALSE, row.names=FALSE)
  # SNP.pvalue <- unique(x4[,c("SNP","P")])
  # write.table(SNP.pvalue, file=paste0(ResultDir,"/","SNP.pvalue"), quote=FALSE, row.names=FALSE)
  ######

  invisible(sys::exec_wait(
    paste0(ResultDir, "/","./plink"),
    args = c(
      "--bfile", paste0(DataDir, "/", finput),
      "--score", paste0(ResultDir, "/","prssummarystat"), 1, 2, 3, "header",
      "--q-score-range", paste0(ResultDir, "/range_list"), paste0(ResultDir, "/","SNP.pvalue"),
      clumpExtract, clumpSNP,
      "--out",
      paste0(ResultDir, "/","PRS"),
      "--silent"
    ),
    std_out = FALSE,
    std_err = FALSE
  ))

  prs <- read.table(paste0(ResultDir, "/","PRS.", pthreshold, ".profile"), header = TRUE)
  pheno.prs <- merge(pheno, prs[, c("FID", "IID", "SCORE")], by = c("FID", "IID"))

  if (pheno_type == "binary") {
    phenoBin.prs <- pheno.prs
    phenoBin.prs$Pheno1 <- as.integer(phenoBin.prs$Pheno1 == 2)
    model <- glm(Pheno1 ~ ., data = phenoBin.prs[, !colnames(phenoBin.prs) %in% c("FID", "IID")], family = binomial)
    McFaddenR2 <- 1 - stats::logLik(model) / stats::logLik(null_model)
    prs.r2 <- McFaddenR2[1]
  } else {
    model <- lm(Pheno1 ~ ., data = pheno.prs[, !colnames(pheno.prs) %in% c("FID", "IID")])
    prs.r2 <- summary(model)$r.squared - summary(null_model)$r.squared
  }

  prs.coef <- summary(model)$coeff["SCORE", ]
  prs.result <- rbind(data.frame(pthreshold, R2 = prs.r2, P = prs.coef[4], BETA = prs.coef[1], SE = prs.coef[2]))
  colnames(prs.result) <- c("Threshold", "R2", "P", "BETA", "SE")

  return(prs.result)
}

## Function 41
######### Added in 3.0
createPRSPlot <- function(prsResult) {
  p1 <- ggplot2::ggplot(data = prsResult, ggplot2::aes(x = factor(prsResult$Threshold), y = prsResult$R2)) +
    ggplot2::geom_text(
      ggplot2::aes(label = paste(prsResult$WriteP)),
      vjust = -1.5,
      hjust = 0,
      angle = 45,
      cex = 2,
      parse = TRUE
    ) +
    ggplot2::scale_y_continuous(limits = c(0, max(prsResult$R2) * 1.25)) +
    ggplot2::xlab("P-value thresholds") +
    ggplot2::ylab("PRS model fit: R^2") +
    ggplot2::geom_bar(ggplot2::aes(fill = -log10(prsResult$P)), stat = "identity") +
    ggplot2::scale_fill_gradient2(
      low = "dodgerblue",
      high = "firebrick",
      mid = "dodgerblue",
      midpoint = 1e-4
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold", size = 7),
      axis.text = ggplot2::element_text(size = 7),
      legend.title = ggplot2::element_text(face = "bold", size = 7),
      legend.text = ggplot2::element_text(size = 7),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::ggtitle("P-value thresholds vs PRS model fit") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold"))

  return(p1)
}

## Function 42
######### Added in 3.0
createSexDistributionPlot <- function(dat) {
  dat <- as.data.frame(dat)
  # Determine the title based on the available data
  title <- if (nrow(dat[dat$SEX == 1,]) != 0 && nrow(dat[dat$SEX == 2,]) != 0) {
    "Best PRS distribution\n(males vs females)"
  } else {
    "Best PRS distribution"
  }

  # Create the plot
  p2 <- ggplot2::ggplot(dat, ggplot2::aes(x = SCORE, color = SEX)) +
    ggplot2::geom_density() +
    ggplot2::ggtitle(title) +
    ggplot2::xlab("PRS") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 10, face = "bold"),
      axis.title.x = ggplot2::element_text(size = 8)
    )

  return(p2)
}

## Function 43
######### Added in 3.0
createBinaryPhenotypePlots <- function(dat, p1, p2) {
  # Convert numeric phenotype to 'control' and 'cases'
  dat$Pheno1 <- factor(ifelse(dat$Pheno1 == 1, "control", "cases"))

  # Generate density plot for overall distribution
  p3 <- ggplot2::ggplot(dat, ggplot2::aes(x = `SCORE`, color = Pheno1)) +
    ggplot2::geom_density() +
    ggplot2::ggtitle("Best PRS distribution\n(cases vs controls)") +
    ggplot2::xlab("PRS") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold")) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = 8))

  # Generate density plots separated by sex
  mdat <- dat[dat$SEX == "Male", ]
  fdat <- dat[dat$SEX == "Female", ]

  p4 <- ggplot2::ggplot(mdat, ggplot2::aes(x = SCORE, color = Pheno1)) +
    ggplot2::geom_density() +
    ggplot2::ggtitle("Best PRS distribution in males\n(cases vs controls)") +
    ggplot2::xlab("PRS") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold")) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = 8))

  p5 <- ggplot2::ggplot(fdat, ggplot2::aes(x = SCORE, color = Pheno1)) +
    ggplot2::geom_density() +
    ggplot2::ggtitle("Best PRS distribution in females\n(cases vs controls)") +
    ggplot2::xlab("PRS") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold")) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = 8))

  # Arrange plots
  if (nrow(mdat) != 0 && nrow(fdat) != 0) {
    plotArrangement <- ggpubr::ggarrange(p1, p2, p3, p4, p5, ncol = 3, nrow = 2, labels = LETTERS[1:5], widths = c(1.5, 1.5, 1.5, 1.5))
  } else {
    plotArrangement <- ggpubr::ggarrange(p1, p2, p3, ncol = 2, nrow = 2, labels = LETTERS[1:3], widths = c(1.5, 1.5, 1.5, 1.5))
  }

  return(plotArrangement)
}

## Function 44
Run_newcovarfile <- function(DataDir,ResultDir,covarfile,covartest){

  #covarfile <- "covarfile.txt"
  covarfile1 <-
    as.data.frame(read.table(paste0(DataDir,"/",covarfile),
                             stringsAsFactors = FALSE,
                             header = TRUE))

  covarfile <- cbind(covarfile1[, 1:2,drop = FALSE], covarfile1[, covartest, drop = FALSE])

  write.table(
    covarfile,
    file = paste0(ResultDir,"/","newcovarfile.txt"),
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    eol = "\r\n",
    sep = " "
  )
}

## Function 45
## Updated in 3.0
FMsub <- function(ResultDir, plot.jpeg, plotname, snp_pval, annotateTopSnp, suggestiveline, genomewideline){


  if (file.exists(paste0(ResultDir,"/allsnpsresults.rda"))[1] == TRUE) {

    load(paste0(ResultDir,"/allsnpsresults.rda"))
    #globalVariables("allsnpsresults")
    XWAS <- data.table::as.data.table(allsnpsresults)
    gc(reset=TRUE)
    XWAS <- na.omit(XWAS)
    gc(reset=TRUE)
    XWAS_ADD <- na.omit(XWAS[XWAS$TEST == "ADD",c("SNP","CHR","BP","P")])
    gc(reset=TRUE)
    # From p-values, calculate chi-squared statistic
    chisq <- qchisq(1-na.omit(XWAS_ADD$P),1)
    lamdaGC <- median(chisq)/qchisq(0.5,1)

    XWAS_ADD_X <- na.omit(XWAS_ADD[XWAS_ADD$CHR == 23,])
    chisq1 <- qchisq(1-na.omit(XWAS_ADD_X$P),1)
    lamdaGC1 <- median(chisq1)/qchisq(0.5,1)

    if (plot.jpeg[1] == TRUE){
      options(bitmapType='cairo')
      grDevices::jpeg(paste0(ResultDir,"/",plotname,".jpeg"),  width = 20,
                      height = 10,
                      units = 'in',
                      res = 300)
      graphics::par(mfrow = c(2, 2))

      uplim <- -log10(min(XWAS_ADD$P))+1
      suppressWarnings(qqman::manhattan(XWAS_ADD, ylim = c(0,uplim),suggestiveline = suggestiveline, genomewideline = genomewideline, annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of GWAS"))
      suppressWarnings(qqman::qq(XWAS_ADD$P, main = paste0(("Q-Q plot of GWAS p-values with GIF = "), round(lamdaGC,3))))

      if (nrow(XWAS_ADD_X)!=0){
        uplim <- -log10(min(XWAS_ADD_X$P))+1
        suppressWarnings(qqman::manhattan(XWAS_ADD_X, ylim = c(0,uplim),suggestiveline = suggestiveline, genomewideline = genomewideline, annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of XWAS"))
        suppressWarnings(qqman::qq(XWAS_ADD_X$P, main = paste0(("Q-Q plot of XWAS p-values with GIF = "), round(lamdaGC1,3))))

      }else{
        print("X-chromosome may not be present.")
      }

      gc(reset=TRUE)
      dev.off()

      return(na.omit(XWAS))

    }else if (plot.jpeg[1] == FALSE){
      par(mar = c(1, 1, 1, 1))
      # graphics::par(mfrow = c(2, 2))
      uplim <- -log10(min(XWAS_ADD$P))+1
      suppressWarnings(qqman::manhattan(XWAS_ADD, ylim = c(0,uplim), suggestiveline = suggestiveline, genomewideline = genomewideline, annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of GWAS"))
      suppressWarnings(qqman::qq(XWAS_ADD$P, main = paste0(("Q-Q plot of GWAS p-values with GIF = "), round(lamdaGC,3))))

      if (nrow(XWAS_ADD_X)!=0){
        uplim <- -log10(min(XWAS_ADD_X$P))+1
        suppressWarnings(qqman::manhattan(XWAS_ADD_X, ylim = c(0,uplim), suggestiveline = suggestiveline, genomewideline = genomewideline, annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of XWAS"))
        suppressWarnings(qqman::qq(XWAS_ADD_X$P, main = paste0(("Q-Q plot of XWAS p-values with GIF = "), round(lamdaGC1,3))))


      }else{
        print("X-chromosome association betas may not be present.")
      }
      gc(reset=TRUE)
    }

    return(na.omit(XWAS))
  } else if (file.exists(paste0(ResultDir,"/","allsnpsresults.rda"))[1] == FALSE) {
    print(paste0("GWAS cannot be performed. Check the log file in ResultDir for checking the error."))
  }
}

## Function 46
## Added in 3.0
paraGwas <- function(chunks,chunk,ResultDir,DataDir,finput,trait, modelv,regress,standard_b,noxsexv,sexv
                     ,interactionv,parameterv,Inphenocovv,covar,covarv,
                     snpfile){

  ## Chunkfile create
  print(paste0("Chunk index processing: ",chunks))
  if (nrow(snpfile)>=(chunks+chunk)){
    snp_names <- snpfile$V2[chunks:(chunks+chunk)]
  }else{
    snp_names <- snpfile$V2[chunks:nrow(snpfile)]
  }
  ## Create snps with chunk
  write.table(snp_names,file = paste0(ResultDir,"/",chunks,"_snps"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  # Make plink bedfiles for this SNP in DataDir

  file_path <- paste0(ResultDir,"/sink_file.txt")

  # Create an empty file (if it doesn't exist)
  if (!file.exists(file_path)) {
    file.create(file_path)
  }

  # Redirect output to the specified file
  sink(file_path)

  invisible(sys::exec_wait(
    paste0(ResultDir,"/","./plink"),
    args = c(
      "--bfile",
      paste0(DataDir,"/",finput),
      "--extract",
      paste0(ResultDir,"/",chunks,"_snps"),
      "--xchr-model",
      modelv,
      "--freq",
      "--ci", 0.95,
      regress, sexv, noxsexv, "intercept", ## 3.0 adding sex
      standard_b,
      interactionv,
      parameterv,
      Inphenocovv,
      covar,
      covarv,
      "--out",
      paste0(ResultDir,"/",chunks,"_ss"),
      "--silent"
    ),
    std_out = FALSE,
    std_err = FALSE
  ))

  # Reset the sink to stop redirecting the output to the file
  if (sink.number() > 0) {
    sink()  # This line resets the output redirection
  }

  if (trait[1] == "binary"){
    single_snp_result <- read.table(paste0(ResultDir,"/",chunks,"_ss.assoc.logistic"),header = T)
  }else if (trait[1] == "quantitative"){
    single_snp_result <- read.table(paste0(ResultDir,"/",chunks,"_ss.assoc.linear"),header = T)
  }
  return(single_snp_result)
}


## Function 47
## Updated in 3.0
FMmain <- function(DataDir, ResultDir, finput, trait, standard_beta, xmodel,
                   sex, xsex, covarfile, interaction, covartest, Inphenocov, plot.jpeg, plotname, snp_pval, annotateTopSnp, suggestiveline, genomewideline, ncores){


  if (xmodel[1] == "FMcomb01"|xmodel[1] == "FMstatrified"){
    modelv <- 1
  }else if (xmodel[1] == "FMcomb02"){
    modelv <- 2
  }

  if (trait[1] == "binary"){
    regress <- "--logistic"
  }else if (trait[1] == "quantitative"){
    regress <- "--linear"}

  if (trait[1] == "quantitative" && standard_beta[1] == TRUE){
    standard_b <- "--standard-beta"
  }else{
    standard_b <- "beta"
  }

  if (xsex[1] == TRUE){
    noxsexv = NULL
  }else{
    noxsexv = "no-x-sex"
  }

  ####### Adding this in 3.0
  if (sex[1] == TRUE){
    noxsexv = NULL
    sexv = "sex"
  }else{
    sexv = NULL
  }

  if (is.null(covarfile)){
    covar <- NULL
    covarv <- NULL
    ##NEWLY
    covartest <- NULL
  }else{

    if (covartest == "ALL") {
      invisible(file.copy(paste0(DataDir,"/",list(covarfile)),ResultDir))

      covarfile <- paste0(ResultDir,"/",covarfile)

    } else if (covartest != "ALL" & !is.null(covarfile)) {
      Run_newcovarfile(DataDir = DataDir, ResultDir=ResultDir,covarfile = covarfile,covartest = covartest)
      covarfile <- paste0(ResultDir,"/newcovarfile.txt")
    }
    covar <- "--covar"
    covarv <- covarfile
  }

  if (interaction[1] == FALSE) {
    interactionv <- NULL
    Inphenocovv <- NULL
    parameterv <- NULL

  }else{
    interactionv <- interaction
    if (Inphenocov[1] == "ALL"){
      Inphenocovv <- NULL
      parameterv <- NULL

    }else if (Inphenocov[1] != "ALL"){
      Inphenocovv <- Inphenocov
      parameterv <- "--parameters"
    }
  }

  if (ncores == 0){

    print("If you want parallel computation, please provide non-zero value for argument ncores.")

    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--bfile",
        paste0(DataDir,"/",finput),
        "--xchr-model",
        modelv,
        "--freq",
        "--ci", 0.95,
        regress,
        sexv,
        noxsexv,
        "intercept", ## 3.0 adding sex
        standard_b,
        interactionv,
        parameterv,
        Inphenocovv,
        covar,
        covarv,
        "--out",
        paste0(ResultDir,"/",xmodel),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))

    if (trait[1] == "binary"){
      single_snp_result <- read.table(paste0(ResultDir,"/",xmodel,".assoc.logistic"),header = T)
    }else if (trait[1] == "quantitative"){
      single_snp_result <- read.table(paste0(ResultDir,"/",xmodel,".assoc.linear"),header = T)
    }
    allsnpsresults <- unique(single_snp_result)
    save(allsnpsresults, file=paste0(ResultDir,"/allsnpsresults.rda"))

  }else{

    snpfile = read.table(paste0(DataDir,"/",finput,".bim"))
    chunk <- round(nrow(snpfile)/ncores)+1
    chunks <- round(seq(1, nrow(snpfile), by = chunk),0)

    #### Parallel computation
    print("Parallel computation is in progress --------->")
    cl <- parallel::makeCluster(mc <- getOption("cl.cores", ncores),type = "FORK")

    invisible(parallel::clusterEvalQ(cl, library(data.table)))
    invisible(parallel::clusterEvalQ(cl, library(parallel)))
    parallel::clusterExport(cl=cl,NULL,envir=environment())

    allsnpsresults <- data.table::as.data.table(data.table::rbindlist(parallel::parLapply(cl,chunks,paraGwas,chunk=chunk,ResultDir=ResultDir,
                                                                                          DataDir=DataDir,finput=finput,trait = trait, modelv=modelv,regress=regress,
                                                                                          standard_b=standard_b,noxsexv=noxsexv,sexv = sexv,
                                                                                          interactionv=interactionv,parameterv=parameterv,
                                                                                          Inphenocovv=Inphenocovv,covar=covar,covarv=covarv,snpfile=snpfile)))
    allsnpsresults <- unique(allsnpsresults)
    save(allsnpsresults, file=paste0(ResultDir,"/allsnpsresults.rda"))

    parallel::stopCluster(cl)
  }

  x <- FMsub(ResultDir = ResultDir, plot.jpeg = plot.jpeg, plotname = plotname, snp_pval =snp_pval, annotateTopSnp = annotateTopSnp, suggestiveline = suggestiveline, genomewideline = genomewideline)

  gc(reset = TRUE)
  gc(reset = TRUE)
  return(x)
}


## Function 48
# Taken from metaseqR package source code
# Copied from ex-CRAN package MADAM and exported (https://rdrr.io/rforge/MADAM/man/fisher.method.html). The man pages are copied from
# the original package.
fisher.sum <- function(p, zero.sub = 0.00001, na.rm = FALSE) {
  if (any(p > 1, na.rm = TRUE) || any(p < 0, na.rm = TRUE))
    stop("You provided bad p-values")
  stopifnot(zero.sub >= 0 & zero.sub <= 1 || length(zero.sub) != 1)
  p[p == 0] <- zero.sub
  if (na.rm)
    p <- p[!is.na(p)]
  S = -2 * sum(log(p))
  res <- data.frame(S = S, num.p = length(p))
  return(res)
}

## Function 49
# Copied from ex-CRAN package MADAM and exported. The man pages are copied from the original package.
fisher.method <-
  function(pvals,
           method = c("fisher"),
           p.corr = c("bonferroni", "BH",
                      "none"),
           zero.sub = 0.00001,
           na.rm = FALSE,
           mc.cores = NULL) {
    stopifnot(method %in% c("fisher"))
    stopifnot(p.corr %in% c("none", "bonferroni", "BH"))
    stopifnot(all(pvals >= 0, na.rm = TRUE) &
                all(pvals <= 1, na.rm = TRUE))
    stopifnot(zero.sub >= 0 & zero.sub <= 1 || length(zero.sub) != 1)
    if (is.null(dim(pvals)))
      stop("pvals must have a dim attribute")
    p.corr <- ifelse(length(p.corr) != 1, "BH", p.corr)
    ##substitute p-values of 0
    pvals[pvals == 0] <- zero.sub
    if (is.null(mc.cores)) {
      fisher.sums <- data.frame(do.call(
        rbind,
        apply(
          pvals,
          1,
          fisher.sum,
          zero.sub = zero.sub,
          na.rm = na.rm
        )
      ))
    }
    else {
      fisher.sums <- parallel::mclapply(1:nrow(pvals), function(i) {
        fisher.sum(pvals[i, ], zero.sub = zero.sub, na.rm = na.rm)
      }, mc.cores = mc.cores)
      fisher.sums <- data.frame(do.call(rbind, fisher.sums))
    }

    rownames(fisher.sums) <- rownames(pvals)
    fisher.sums$p.value <-
      1 - stats::pchisq(fisher.sums$S, df = 2 * fisher.sums$num.p)
    fisher.sums$p.adj <- switch(
      p.corr,
      bonferroni = stats::p.adjust(fisher.sums$p.value, "bonferroni"),
      BH = stats::p.adjust(fisher.sums$p.value, "BH"),
      none = fisher.sums$p.value
    )
    return(fisher.sums)
  }

## Function 50
# Copied from ex-CRAN package MADAM and exported. The man pages are copied from the original package.
fisher.method.perm <-
  function(pvals,
           p.corr = c("bonferroni", "BH", "none"),
           zero.sub = 0.00001,
           B = 10000,
           mc.cores = NULL,
           blinker = 1000) {
    stopifnot(is.na(blinker) || blinker > 0)
    stopifnot(p.corr %in% c("none", "bonferroni", "BH"))
    stopifnot(all(pvals >= 0, na.rm = TRUE) &
                all(pvals <= 1, na.rm = TRUE))
    stopifnot(zero.sub >= 0 & zero.sub <= 1 || length(zero.sub) != 1)
    if (is.null(dim(pvals)))
      stop("pvals must have a dim attribute")
    p.corr <- ifelse(length(p.corr) != 1, "BH", p.corr)
    pvals[pvals == 0] <- zero.sub

    res.perm <- lapply(1:nrow(pvals), function(i) {
      if (!is.na(blinker) & i %% blinker == 0)
        message("=", appendLF = FALSE)
      ##which studies contribute to S (don't have a NA in row i)
      good.p <- which(!is.na(pvals[i, ]))
      S.obs = fisher.sum(pvals[i, good.p], na.rm = FALSE)
      if (is.null(mc.cores)) {
        S.rand <- unlist(lapply(1:B, function(b) {
          ##get non NA p-values from studies contributing to S
          myp <- sapply(good.p, function(pc) {
            sample(stats::na.exclude(pvals[, pc]), 1)
          })
          fisher.sum(myp)$S
        }))
      } else {
        S.rand <- unlist(parallel::mclapply(1:B, function(b) {
          ##get non NA p-values from studies contributing to S
          myp <- sapply(good.p, function(pc) {
            sample(stats::na.exclude(pvals[, pc]), 1)
          })
          fisher.sum(myp)$S
        }, mc.cores = mc.cores))
      }
      p.value <- sum(S.rand >= S.obs$S) / B
      data.frame(S = S.obs$S,
                 num.p = S.obs$num.p,
                 p.value = p.value)
    })
    res.perm <- data.frame(do.call(rbind, res.perm))

    if (!is.na(blinker) && blinker > 0)
      message()
    ## rownames(res.perm) <- rownames(pvals)
    res.perm$p.adj <- switch(
      p.corr,
      bonferroni = stats::p.adjust(res.perm$p.value, "bonferroni"),
      BH = stats::p.adjust(res.perm$p.value, "BH"),
      none = res.perm$p.value
    )
    return(res.perm)
  }


## Function 51
stouffer.method <-
  function(pvals,
           method = c("stouffer"),
           p.corr = c("bonferroni", "BH",
                      "none"),
           zero.sub = 0.00001,
           na.rm = FALSE,
           mc.cores = NULL) {

    stopifnot(method %in% c("stouffer"))
    stopifnot(p.corr %in% c("none", "bonferroni", "BH"))
    stopifnot(all(pvals >= 0, na.rm = TRUE) &
                all(pvals <= 1, na.rm = TRUE))
    stopifnot(zero.sub >= 0 & zero.sub <= 1 || length(zero.sub) != 1)
    if (is.null(dim(pvals)))
      stop("pvals must have a dim attribute")
    p.corr <- ifelse(length(p.corr) != 1, "BH", p.corr)
    ##substitute p-values of 0
    pvals[pvals == 0] <- zero.sub

    i <- 1:nrow(pvals)
    stpFun <- function(i){
      p1 <-
        poolr::stouffer(
          as.matrix(pvals[i, ]),
          adjust = "none",
          side = 2,
          nearpd = TRUE
        )
      p2 <- p1$p
      pvals <- as.data.frame(pvals)
      p3 <- as.data.frame(cbind(pvals[i, , drop = FALSE], p2))
      p3$p.adj <-
        switch(
          p.corr,
          bonferroni = stats::p.adjust(p3$p2, "bonferroni", n = length(as.matrix(pvals[i, ]))),
          BH = stats::p.adjust(p3$p2, "BH", n = length(as.matrix(pvals[i, ]))),
          none = p3$p2
        )
      return(p3)
    }

    stofP <- data.table::rbindlist(lapply(i,stpFun))

    return(stofP)
  }


## Function 52
#### Parallel Stouffer method
paraStouffer <- function(chunks,chunk,pvals,MF.p.corr,MF.zero.sub,MF.na.rm){

  if (nrow(pvals)>=(chunks+chunk)){
    pval_chunk <- pvals[chunks:(chunks+chunk),]
  }else{
    pval_chunk <- pvals[chunks:nrow(pvals),]
  }
  stouffer.method <-
    function(pvals,
             method = c("stouffer"),
             p.corr = c("bonferroni", "BH",
                        "none"),
             zero.sub = 0.00001,
             na.rm = FALSE,
             mc.cores = NULL) {

      stopifnot(method %in% c("stouffer"))
      stopifnot(p.corr %in% c("none", "bonferroni", "BH"))
      stopifnot(all(pvals >= 0, na.rm = TRUE) &
                  all(pvals <= 1, na.rm = TRUE))
      stopifnot(zero.sub >= 0 & zero.sub <= 1 || length(zero.sub) != 1)
      if (is.null(dim(pvals)))
        stop("pvals must have a dim attribute")
      p.corr <- ifelse(length(p.corr) != 1, "BH", p.corr)
      ##substitute p-values of 0
      pvals[pvals == 0] <- zero.sub

      i <- 1:nrow(pvals)
      stpFun <- function(i){
        p1 <-
          poolr::stouffer(
            as.matrix(pvals[i, ]),
            adjust = "none",
            side = 2,
            nearpd = TRUE
          )
        p2 <- p1$p
        pvals <- as.data.frame(pvals)
        p3 <- as.data.frame(cbind(pvals[i, , drop = FALSE], p2))
        p3$p.adj <-
          switch(
            p.corr,
            bonferroni = stats::p.adjust(p3$p2, "bonferroni", n = length(as.matrix(pvals[i, ]))),
            BH = stats::p.adjust(p3$p2, "BH", n = length(as.matrix(pvals[i, ]))),
            none = p3$p2
          )
        return(p3)
      }

      stofP <- data.table::rbindlist(lapply(i,stpFun))

      return(stofP)
    }
  Pnew1 <-
    stouffer.method(
      pvals = pval_chunk ,
      p.corr = MF.p.corr,
      zero.sub = MF.zero.sub,
      na.rm = MF.na.rm,
      mc.cores = NULL
    )
  return(Pnew1)
}

## Function 53
###### Added in 3.0
applyStoufferMethod <- function(pvals, MF.p.corr, MF.zero.sub, MF.na.rm, MF.mc.cores, ncores) {
  if (ncores == 0) {
    Pnew <- stouffer.method(
      pvals = pvals,
      p.corr = MF.p.corr,
      zero.sub = MF.zero.sub,
      na.rm = MF.na.rm,
      mc.cores = MF.mc.cores
    )
  } else {
    chunk <- round(nrow(pvals) / ncores) + 1
    chunks <- round(seq(1, nrow(pvals), by = chunk), 0)

    print("cluster making started")
    cl <- parallel::makeCluster(mc <- getOption("cl.cores", ncores), type = "FORK")

    invisible(parallel::clusterEvalQ(cl, library(data.table)))
    invisible(parallel::clusterEvalQ(cl, library(parallel)))
    print("cluster export started")
    parallel::clusterExport(cl = cl, NULL, envir = environment())

    print("cluster making done")
    Pnew <- data.table::as.data.table(data.table::rbindlist(parallel::parLapply(cl, chunks, paraStouffer, chunk = chunk, pvals = pvals, MF.p.corr = MF.p.corr, MF.zero.sub = MF.zero.sub, MF.na.rm = MF.na.rm)))

    parallel::stopCluster(cl)
    print("clusters stopped.")
  }
  return(Pnew)
}

## Function 54
###### Added in 3.0
generateGWASPlots <- function(plot.jpeg, plotname, FemaleWAS, MaleWAS, gwas.t2, gwas.b2, Result1, XWAS_ADD_X1, snp_pval, annotateTopSnp, suggestiveline, genomewideline, lamdaGC, lamdaGC1, ResultDir) {


  if (plot.jpeg[1] == TRUE) {

    options(bitmapType='cairo')
    grDevices::jpeg(paste0(plotname,".jpeg"),  width = 20,
                    height = 10,
                    units = 'in',
                    res = 300)
    graphics::par(mfrow = c(2, 2))

    file_path <- paste0(ResultDir,"/sink_file.txt")

    # Create an empty file (if it doesn't exist)
    if (!file.exists(file_path)) {
      file.create(file_path)
    }

    # Redirect output to the specified file
    sink(file_path)

    suppressWarnings(gmirror(top=FemaleWAS, bottom=MaleWAS, tline= snp_pval, bline=snp_pval,
                             toptitle="GWAS of females", bottomtitle = "GWAS of males",
                             highlight_p = c(snp_pval,snp_pval), highlighter="green", chrblocks = TRUE,file = "Stratified_GWAS"))

    # Reset the sink to stop redirecting the output to the file
    if (sink.number() > 0) {
      sink()  # This line resets the output redirection
    }

    print("Miami plot of stratified GWAS is saved in working directory.")
    gc(reset=TRUE)

    if (nrow(gwas.t2)!=0 && nrow(gwas.b2)!=0){

      file_path <- paste0(ResultDir,"/sink_file.txt")

      # Create an empty file (if it doesn't exist)
      if (!file.exists(file_path)) {
        file.create(file_path)
      }

      # Redirect output to the specified file
      sink(file_path)

      suppressWarnings(gmirror(top=gwas.t2, bottom=gwas.b2, tline=snp_pval, bline=snp_pval,
                               toptitle="XWAS of females", bottomtitle = "XWAS of males",
                               highlight_p = c(snp_pval,snp_pval), highlighter="green", chrblocks = TRUE, file = paste0(ResultDir,"/","Stratified_XWAS")))


      # Reset the sink to stop redirecting the output to the file
      if (sink.number() > 0) {
        sink()  # This line resets the output redirection
      }

      print("Miami plot of stratified XWAS is saved in working directory.")

    }else{
      print("Miami plot for stratified XWAS cannot be drawn.")
    }
    gc(reset=TRUE)

    ## Updating this in 3.0
    uplim <- -log(min(Result1$P))+1
    #suppressWarnings(qqman::manhattan(Result1, ylim = c(0,16), suggestiveline = suggestiveline, genomewideline = genomewideline ,annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of male-female combined GWAS"))
    suppressWarnings(qqman::manhattan(Result1, ylim = c(0,uplim), suggestiveline = suggestiveline, genomewideline = genomewideline ,annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of male-female combined GWAS"))

    gc(reset=TRUE)

    if (nrow(XWAS_ADD_X1)!=0){
      uplim <- -log10(min(XWAS_ADD_X1$P))+1
      suppressWarnings(qqman::manhattan(XWAS_ADD_X1, ylim = c(0,uplim), suggestiveline = suggestiveline, genomewideline = genomewideline, annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of male-female combined XWAS"))
    }else{
      print("There may not be any X chromosome in the data.")
    }
    gc(reset=TRUE)
    if (sum(Result1$P) == nrow(Result1)){
      print("All adjusted p values are 1, qq plot cannot be created.")
    }else{
      suppressWarnings(qqman::qq(Result1$P, main = paste0(("Q-Q plot of GWAS male-female combined with GIF = "), round(lamdaGC,3))))
      gc(reset=TRUE)
    }

    if (sum(XWAS_ADD_X1$P) == nrow(XWAS_ADD_X1)){
      print("All adjusted p values are 1, qq plot cannot be created.")
    }else{
      suppressWarnings(qqman::qq(XWAS_ADD_X1$P, main = paste0(("Q-Q plot of XWAS male-female combined p-values with GIF = "), lamdaGC1)))
    }

    dev.off()
    # This section includes the creation of JPEG files and plotting using gmirror, qqman::manhattan, qqman::qq, etc.
  } else {

    #graphics::par(mfrow = c(2, 2))
    gc(reset=TRUE)

    file_path <- paste0(ResultDir,"/sink_file.txt")

    # Create an empty file (if it doesn't exist)
    if (!file.exists(file_path)) {
      file.create(file_path)
    }

    # Redirect output to the specified file
    sink(file_path)

    print(suppressWarnings(gmirror(top=FemaleWAS, bottom=MaleWAS, tline= snp_pval, bline=snp_pval,
                                   toptitle="GWAS of females", bottomtitle = "GWAS of males",
                                   highlight_p = c(snp_pval,snp_pval), highlighter="green", chrblocks = TRUE, file = paste0(ResultDir,"/","Stratified_GWAS"))))

    # Reset the sink to stop redirecting the output to the file
    if (sink.number() > 0) {
      sink()  # This line resets the output redirection
    }

    gc(reset=TRUE)

    file_path <- paste0(ResultDir,"/sink_file.txt")

    # Create an empty file (if it doesn't exist)
    if (!file.exists(file_path)) {
      file.create(file_path)
    }

    # Redirect output to the specified file
    sink(file_path)

    print(suppressWarnings(gmirror(top=gwas.t2, bottom=gwas.b2, tline=snp_pval, bline=snp_pval,
                                   toptitle="XWAS of females", bottomtitle = "XWAS of males",
                                   highlight_p = c(snp_pval,snp_pval), highlighter="green", chrblocks = TRUE, file = paste0(ResultDir,"/","Stratified_XWAS"))))
    # Reset the sink to stop redirecting the output to the file
    if (sink.number() > 0) {
      sink()  # This line resets the output redirection
    }

    gc(reset=TRUE)
    par(mar = c(1, 1, 1, 1))
    uplim <- -log10(min(Result1$P))+1
    suppressWarnings(qqman::manhattan(Result1, ylim = c(0,uplim), suggestiveline = suggestiveline, genomewideline = genomewideline ,annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of male-female combined GWAS"))
    gc(reset=TRUE)
    uplim <- -log10(min(XWAS_ADD_X1$P))+1
    suppressWarnings(qqman::manhattan(XWAS_ADD_X1, ylim = c(0,uplim), suggestiveline = suggestiveline, genomewideline = genomewideline, annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of male-female combined XWAS"))
    gc(reset=TRUE)

    if (sum(Result1$P) == nrow(Result1)){
      print("All adjusted p values are 1, qq plot cannot be created.")
    }else{
      gc(reset=TRUE)
      suppressWarnings(qqman::qq(Result1$P, main = paste0(("Q-Q plot of male-female combined GWAS with GIF = "), round(lamdaGC,3))))
      gc(reset=TRUE)
    }

    if (sum(XWAS_ADD_X1$P) == nrow(XWAS_ADD_X1)){
      print("All adjusted p values are 1, qq plot cannot be created.")
    }else{
      gc(reset=TRUE)
      qqman::qq(XWAS_ADD_X1$P, main = paste0(("Q-Q plot of male-female combined XWAS with GIF = "), round(lamdaGC1,3)))
      gc(reset=TRUE)
    }

    # This section includes plotting using gmirror, qqman::manhattan, qqman::qq, etc.
  }
}

## Function 55
###### Added in 3.0
createPlots <- function(ResultDir, plotname, FemaleWAS, MaleWAS, gwas.t2, gwas.b2, Result, XWAS_ADD_X, snp_pval, suggestiveline, genomewideline, annotateTopSnp, lamdaGC, lamdaGC1, plot.jpeg) {

  if (plot.jpeg[1] == TRUE) {
    options(bitmapType='cairo')
    grDevices::jpeg(paste0(ResultDir,"/",plotname,".jpeg"), width = 20, height = 10, units = 'in', res = 300)
    graphics::par(mfrow = c(2, 2))
  }

  # Stratified GWAS plot
  file_path <- paste0(ResultDir,"/sink_file.txt")

  # Create an empty file (if it doesn't exist)
  if (!file.exists(file_path)) {
    file.create(file_path)
  }

  # Redirect output to the specified file
  sink(file_path)

  print(suppressWarnings(gmirror(top=FemaleWAS, bottom=MaleWAS, tline= snp_pval, bline=snp_pval,
                                 toptitle="GWAS of females", bottomtitle = "GWAS of males",
                                 highlight_p = c(snp_pval,snp_pval), highlighter="green", chrblocks = TRUE, file = paste0(ResultDir,"/","Stratified_GWAS"))))

  # Reset the sink to stop redirecting the output to the file
  if (sink.number() > 0) {
    sink()  # This line resets the output redirection
  }

  if (nrow(gwas.t2) != 0 && nrow(gwas.b2) != 0) {

    file_path <- paste0(ResultDir,"/sink_file.txt")

    # Create an empty file (if it doesn't exist)
    if (!file.exists(file_path)) {
      file.create(file_path)
    }

    # Redirect output to the specified file
    sink(file_path)

    print(suppressWarnings(gmirror(top=gwas.t2, bottom=gwas.b2, tline=snp_pval, bline=snp_pval,
                                   toptitle="XWAS of females", bottomtitle = "XWAS of males",
                                   highlight_p = c(snp_pval,snp_pval), highlighter="green", chrblocks = TRUE, file = paste0(ResultDir,"/","Stratified_XWAS"))))

    # Reset the sink to stop redirecting the output to the file
    if (sink.number() > 0) {
      sink()  # This line resets the output redirection
    }

  } else {
    print("Miami plot for stratified XWAS cannot be drawn.")
  }

  # Manhattan plot
  uplim <- -log10(min(Result$P, na.rm = TRUE)) + 1
  suppressWarnings(qqman::manhattan(Result, ylim = c(0, uplim), suggestiveline = suggestiveline, genomewideline = genomewideline, annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of male-female combined GWAS"))

  # Manhattan plot for XWAS
  if (nrow(XWAS_ADD_X) != 0) {
    uplim <- -log10(min(XWAS_ADD_X$P, na.rm = TRUE)) + 1
    suppressWarnings(qqman::manhattan(XWAS_ADD_X, ylim = c(0, uplim), suggestiveline = suggestiveline, genomewideline = genomewideline, annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of male-female combined XWAS"))
  } else {
    print("There may not be any X chromosome in the data.")
  }

  # QQ plot
  if (sum(Result$P, na.rm = TRUE) != nrow(Result)) {
    suppressWarnings(qqman::qq(Result$P, main = paste0("Q-Q plot of GWAS male-female combined with GIF = ", round(lamdaGC, 3))))
  } else {
    print("All adjusted p values are 1, qq plot cannot be created.")
  }

  # QQ plot for XWAS
  if (sum(XWAS_ADD_X$P, na.rm = TRUE) != nrow(XWAS_ADD_X)) {
    suppressWarnings(qqman::qq(XWAS_ADD_X$P, main = paste0("Q-Q plot of XWAS male-female combined p-values with GIF = ", round(lamdaGC1, 3))))
  } else {
    print("All adjusted p values are 1, qq plot cannot be created.")
  }

  if (plot.jpeg[1] == TRUE) {
    dev.off()
  }
}

## Function 56
###### Updated in 3.0
FMcomb_sub <- function(ResultDir, combtest, MF.p.corr,
                       MF.zero.sub,
                       MF.na.rm,
                       MF.mc.cores,
                       B,
                       plot.jpeg,
                       plotname,
                       snp_pval,
                       annotateTopSnp,
                       suggestiveline,
                       genomewideline,
                       ncores){


  load(paste0(ResultDir,"/MaleWAS.Rda"))
  MaleWAS <- data.table::as.data.table(MaleWAS)
  gc(reset=TRUE)
  load(paste0(ResultDir,"/FemaleWAS.Rda"))
  FemaleWAS <- data.table::as.data.table(FemaleWAS)
  gc(reset=TRUE)
  MFWAS <- merge(FemaleWAS, MaleWAS, by = c("SNP", "A1", "TEST"))
  gc(reset=TRUE)
  MFWAS1 <- MFWAS[MFWAS$TEST == "ADD",]
  pvals <- as.data.frame(MFWAS1[, c(12, 21)])

  if (combtest[1] == "stouffer.method"){

    if (ncores == 0){

      Pnew <-
        stouffer.method(
          pvals = pvals ,
          p.corr = MF.p.corr,
          zero.sub = MF.zero.sub,
          na.rm = MF.na.rm,
          mc.cores = MF.mc.cores
        )

    }else{

      chunk <- round(nrow(pvals)/ncores)+1
      chunks <- round(seq(1, nrow(pvals), by = chunk),0)

      #### Parallel computation
      print("Parallel computation is in progress --------->")
      cl <- parallel::makeCluster(mc <- getOption("cl.cores", ncores),type = "FORK")

      invisible(parallel::clusterEvalQ(cl, library(data.table)))
      invisible(parallel::clusterEvalQ(cl, library(parallel)))
      parallel::clusterExport(cl=cl,NULL,envir=environment())

      Pnew <- data.table::as.data.table(data.table::rbindlist(parallel::parLapply(cl,chunks,paraStouffer,chunk=chunk,pvals=pvals,MF.p.corr=MF.p.corr,MF.zero.sub=MF.zero.sub,MF.na.rm=MF.na.rm)))

      parallel::stopCluster(cl)
    }

    Result <- cbind(MFWAS1[, 1:5], Pnew[, 3:4])
    gc(reset=TRUE)
    Result <- Result[Result$TEST == "ADD",c("SNP","CHR.x","BP.x","p2")]# we could choose "p.adj" as well.
    gc(reset=TRUE)
    colnames(Result) <- c("SNP","CHR","BP","P")

    XWAS_ADD_X <- Result[Result$CHR == 23,]

  } else if (combtest[1] == "fisher.method") {

    Pnew <-
      fisher.method(
        pvals = pvals ,
        p.corr = MF.p.corr,
        zero.sub = MF.zero.sub,
        na.rm = MF.na.rm,
        mc.cores = MF.mc.cores
      )
    Result <- cbind(MFWAS[, 1:5], Pnew[, 3:4])
    gc(reset=TRUE)
    Result <- Result[Result$TEST == "ADD",c("SNP","CHR.x","BP.x","p.value")]# we could choose "p.adj" as well.
    colnames(Result) <- c("SNP","CHR","BP","P")
    gc(reset=TRUE)
    XWAS_ADD_X <- Result[Result$CHR == 23,]

  } else if (combtest[1] == "fisher.method.perm") {
    Pnew <-
      fisher.method.perm(
        pvals = pvals,
        p.corr = MF.p.corr,
        zero.sub =
          MF.zero.sub,
        B = B,
        mc.cores = MF.mc.cores,
        blinker = 1000
      )
    Result <- cbind(MFWAS[, 1:7], Pnew[, 3:4])
    gc(reset=TRUE)
    Result <- Result[Result$TEST == "ADD",c("SNP","CHR","BP","P")]
    gc(reset=TRUE)
    XWAS_ADD_X <- Result[Result$CHR == 23,]


  }

  # From p-values, calculate chi-squared statistic
  chisq <- qchisq(1-na.omit(Result$P),1)
  lamdaGC <- median(chisq)/qchisq(0.5,1)
  chisq1 <- qchisq(1-na.omit(XWAS_ADD_X$P),1)
  lamdaGC1 <- median(chisq1)/qchisq(0.5,1)

  # Stratified GWAS plot
  # Manhattan and QQ-plots will be produced using P values from additive effect only. For all other tests, please use the final output.
  FemaleWAS <- na.omit(FemaleWAS[FemaleWAS$TEST=="ADD",c("SNP","CHR","BP","P")])
  gc(reset=TRUE)
  MaleWAS <- na.omit(MaleWAS[MaleWAS$TEST=="ADD",c("SNP","CHR","BP","P")])
  gc(reset=TRUE)
  colnames(FemaleWAS) <- c("SNP","CHR","POS","pvalue")
  colnames(MaleWAS) <- c("SNP","CHR","POS","pvalue")
  FemaleWAS <- as.data.frame(FemaleWAS)
  FemaleWAS[FemaleWAS$CHR == "23","CHR"]<-"X"
  MaleWAS <- as.data.frame(MaleWAS)
  MaleWAS[MaleWAS$CHR == "23","CHR"]<-"X"

  # Stratified XWAS plot
  gwas.t2 <- FemaleWAS[FemaleWAS$CHR=="X",]
  gwas.b2 <- MaleWAS[MaleWAS$CHR=="X",]


  # Call to the helper function for plotting
  print("Plots are initiated.")
  createPlots(ResultDir, plotname, FemaleWAS, MaleWAS, gwas.t2, gwas.b2, Result, XWAS_ADD_X, snp_pval, suggestiveline, genomewideline, annotateTopSnp, lamdaGC, lamdaGC1, plot.jpeg)

  gc(reset=TRUE)
  return(na.omit(Result))
}


## Function 57
## Updated in 3.0
FMcomb <-
  function(DataDir,
           ResultDir,
           trait,
           standard_beta,
           xmodel,
           covarfile,
           covartest,
           interaction,
           Inphenocov,
           combtest,
           B,
           MF.p.corr,
           MF.zero.sub,
           MF.na.rm,
           MF.mc.cores,
           plot.jpeg,
           plotname,
           snp_pval,
           annotateTopSnp,
           suggestiveline,
           genomewideline,
           ncores) {



    gc(reset=TRUE)
    # Here DataDir becomes ResultDir, need to put the covarfile always in ResultDir by copying
    if (!is.null(covarfile)){
      invisible(file.copy(paste0(DataDir,"/",list(covarfile)),ResultDir))
    }else{
      covarfile = covarfile
    }
    ## Removing perm = perm, mperm = mperm, genedrop = genedrop
    MaleWAS <- FMmain(DataDir = ResultDir, ResultDir = ResultDir, finput = "finput.male", trait = trait, standard_beta = standard_beta, xmodel = xmodel,
                      sex = FALSE, xsex = FALSE, covarfile = covarfile, interaction = interaction, covartest = covartest, Inphenocov = Inphenocov, plot.jpeg = TRUE, plotname = paste0(xmodel,"_MaleWAS"), snp_pval = snp_pval, annotateTopSnp = annotateTopSnp, suggestiveline = suggestiveline, genomewideline = genomewideline, ncores = ncores)
    gc(reset=TRUE)
    MaleWAS <- na.omit(MaleWAS)
    save(MaleWAS, file = paste0(ResultDir,"/MaleWAS.Rda"))
    gc(reset=TRUE)
    FemaleWAS <- FMmain(DataDir = ResultDir, ResultDir = ResultDir, finput = "finput.female", trait = trait, standard_beta = standard_beta, xmodel = xmodel,
                        sex = FALSE, xsex = FALSE, covarfile = covarfile, interaction = interaction, covartest = covartest, Inphenocov = Inphenocov, plot.jpeg = TRUE, plotname = paste0(xmodel,"_FemaleWAS"), snp_pval = snp_pval, annotateTopSnp = annotateTopSnp, suggestiveline = suggestiveline, genomewideline = genomewideline, ncores = ncores)
    gc(reset=TRUE)
    FemaleWAS <- na.omit(FemaleWAS)
    save(FemaleWAS, file = paste0(ResultDir,"/FemaleWAS.Rda"))

    CombinedWAS <- FMcomb_sub(ResultDir = ResultDir, combtest = combtest, MF.p.corr = MF.p.corr,
                              MF.zero.sub = MF.zero.sub,
                              MF.na.rm = MF.na.rm,
                              MF.mc.cores = MF.mc.cores,
                              B = B,
                              plot.jpeg = plot.jpeg,
                              plotname = plotname,
                              snp_pval = snp_pval,
                              annotateTopSnp = annotateTopSnp,
                              suggestiveline = suggestiveline,
                              genomewideline = genomewideline,
                              ncores = ncores)
    gc(reset=TRUE)
    CombinedWAS <- na.omit(CombinedWAS)
    save(CombinedWAS,file = paste0(ResultDir,"/","CombinedWAS.Rda"))
    gc(reset=TRUE)
    print(paste0("Three dataframes such as, CombinedWAS, MaleWAS and FemaleWAS are produced in", ResultDir))
    return(list(CombinedWAS = CombinedWAS, MaleWAS = MaleWAS, FemaleWAS = FemaleWAS))
  }


## Function 58
## Added in 3.0
paraGwasAuto <- function(chunks,chunk,ResultDir,finput,regress,sexv, standard_b,interactionv,parameterv,Inphenocovv,covar,covarv,
                         snpfile){

  ## Chunkfile create
  if (nrow(snpfile)>=(chunks+chunk)){
    snp_names <- snpfile$V2[chunks:(chunks+chunk)]
  }else{
    snp_names <- snpfile$V2[chunks:nrow(snpfile)]
  }

  write.table(snp_names,file = paste0(ResultDir,"/",chunks,"_snps"), quote = FALSE, row.names = FALSE, col.names = FALSE)

  # make plink bedfiles for this SNP in DataDir
  invisible(sys::exec_wait(
    paste0(ResultDir,"/","./plink"),
    args = c(
      "--bfile",
      paste0(ResultDir,"/",finput),
      "--extract",
      paste0(ResultDir,"/",chunks,"_snps"),
      "--freq",
      "--ci", 0.95,
      regress,
      sexv,
      noxsexv,
      "intercept",
      standard_b,
      interactionv,
      parameterv,
      Inphenocovv,
      covar,
      covarv,
      "--out",
      paste0(ResultDir,"/",chunks,"_ss"),
      "--silent"
    ),
    std_out = FALSE,
    std_err = FALSE
  ))

  single_snp_result <- read.table(paste0(ResultDir,"/",chunks,"_ss.assoc.logistic"),header = T)
  single_snp_result <- unique(single_snp_result)
  return(single_snp_result)
}

## Function 59
## Added in 3.0
autoFun <- function(DataDir, ResultDir, finput, sex, standard_beta, covarfile, interaction, covartest, Inphenocov, ncores){

  regress <- "--logistic"
  standard_b <- "beta"

  ####### Adding this in 3.0
  if (sex[1] == TRUE){
    sexv = "sex"
  }else{
    sexv = NULL
  }

  if (is.null(covarfile)){
    covar <- NULL
    covarv <- NULL
  }else{

    if (covartest == "ALL") {
      ## Copy covarfile to ResultDir
      invisible(file.copy(paste0(DataDir,"/",list(covarfile)),ResultDir))

      covarfile <- paste0(ResultDir,"/",covarfile)

    } else if (covartest != "ALL" & !is.null(covarfile)) {
      Run_newcovarfile(DataDir = DataDir, ResultDir=ResultDir,covarfile = covarfile,covartest = covartest)
      covarfile <- paste0(ResultDir,"/newcovarfile.txt")
    }
    covar <- "--covar"
    covarv <- covarfile

  }

  if (interaction[1] == FALSE) {
    interactionv <- NULL
    Inphenocovv <- NULL
    parameterv <- NULL

  }else{
    interactionv <- interaction
    if (Inphenocov[1] == "ALL"){
      Inphenocovv <- NULL
      parameterv <- NULL

    }else if (Inphenocov[1] != "ALL"){
      Inphenocovv <- Inphenocov
      parameterv <- "--parameters"
    }
  }

  if (ncores == 0){
    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./plink"),
      args = c(
        "--bfile",
        paste0(ResultDir,"/",finput),
        "--freq",
        "--ci", 0.95,
        regress,
        sexv,
        noxsexv,
        "intercept",
        standard_b,
        interactionv,
        parameterv,
        Inphenocovv,
        covar,
        covarv,
        "--out",
        paste0(ResultDir,"/","AutoWAS"),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))

  }else{

    snpfile = read.table(paste0(ResultDir,"/",finput,".bim"))
    chunk <- round(nrow(snpfile)/ncores)+1
    chunks <- round(seq(1, nrow(snpfile), by = chunk),0)

    #### Parallel computation
    print("Parallel computation is in progress --------->")
    cl <- parallel::makeCluster(mc <- getOption("cl.cores", ncores),type = "FORK")

    invisible(parallel::clusterEvalQ(cl, library(data.table)))
    invisible(parallel::clusterEvalQ(cl, library(parallel)))
    parallel::clusterExport(cl=cl,NULL,envir=environment())

    all_snps_results <- data.table::as.data.table(data.table::rbindlist(parallel::parLapply(cl,chunks,paraGwasAuto,chunk=chunk,ResultDir=ResultDir,finput=finput,regress=regress,sexv = sexv,
                                                                                            standard_b=standard_b,interactionv=interactionv,parameterv=parameterv,
                                                                                            Inphenocovv=Inphenocovv,covar=covar,covarv=covarv,snpfile=snpfile)))
    all_snps_results <- unique(all_snps_results)

    parallel::stopCluster(cl)

    write.table(all_snps_results, file= paste0(ResultDir,"/","AutoWAS",".assoc.logistic"), quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\r\n", sep = " ")
    rm(all_snps_results)
    ftemp <- list.files(paste0(ResultDir,"/"),pattern = paste0(chunks,"_ss"))
    invisible(file.remove(paste0(ResultDir,"/",ftemp)))
    ftemp <- list.files(paste0(ResultDir,"/"),pattern = paste0(chunks,"_snps"))
    invisible(file.remove(paste0(ResultDir,"/",ftemp)))
  }
  fl <- paste0("AutoWAS",".assoc.logistic")

  if (file.exists(paste0(ResultDir,"/",fl))[1] == TRUE) {
    XWAS <-
      read.table(file = paste0(ResultDir,"/",fl),
                 stringsAsFactors = FALSE,
                 header = TRUE)
    return(XWAS)

  }else if (file.exists(paste0(ResultDir,"/",fl))[1] == FALSE) {
    print(paste0("AutosomeWAS cannot be performed. Check the ", stringr::str_sub(fl, 1,5),"log file in DataDir for checking the error."))
  }

}

## Function 60
## Added in 3.0
createPlotsXCMA <- function(XWAS_ADD, XWAS_ADD_X, ResultDir, plotname, snp_pval, lamdaGC, lamdaGC1, suggestiveline, genomewideline, annotateTopSnp, plot.jpeg) {
  if (plot.jpeg[1] == TRUE) {
    options(bitmapType='cairo')
    grDevices::jpeg(paste0(ResultDir,"/",plotname,".jpeg"), width = 20, height = 10, units = 'in', res = 300)
    graphics::par(mfrow = c(2, 2))
  } else {
    par(mar = c(1, 1, 1, 1))
  }

  uplim_GWAS <- -log10(min(XWAS_ADD$P, na.rm = TRUE)) + 1
  suppressWarnings(qqman::manhattan(XWAS_ADD, ylim = c(0, uplim_GWAS), suggestiveline = suggestiveline, genomewideline = genomewideline, annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of GWAS"))
  suppressWarnings(qqman::qq(XWAS_ADD$P, main = paste0("Q-Q plot of GWAS p-values with GIF = ", round(lamdaGC,3))))

  uplim_XWAS <- -log10(min(XWAS_ADD_X$P, na.rm = TRUE)) + 1
  suppressWarnings(qqman::manhattan(XWAS_ADD_X, ylim = c(0, uplim_XWAS), suggestiveline = suggestiveline, genomewideline = genomewideline, annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of XWAS"))
  suppressWarnings(qqman::qq(XWAS_ADD_X$P, main = paste0("Q-Q plot of XWAS p-values with GIF = ", round(lamdaGC1,3))))

  if (plot.jpeg[1] == TRUE) {
    dev.off()
  }
}


## Function 61
# This function is used to test the association between an X chromosomal marker and a binary trait. One SNP at a time. # Add sex as covariate ## Can't be for this function
#
# Arguments of the function
# data: The data should contain the information of phenotype, genotype, sex, and possible additional covariates . Each row represents the data of a specific individual.
# The first column records the phenotype information, with 1 representing the case and 0 representing the control.
# The second column provides the genotype information, with 0, 1 and 2 representing the number of risk alleles.
# The third column contains the sex information, with 0 being male and 1 being female.
# The remaining columns records the information of possible additional covariates.
#data <- read.table("/projects/b1137/BBose/ProjectGWAS/XWAS_QC/Panscan/DG_98846_1_H55_c1_panscan_v3/han.csv",sep = ",",header = TRUE)

XCMAX4 <- function(data){
  D= as.matrix(data[,1])
  G= as.matrix(data[,2])

  X=as.matrix(data[,-c(1,2)])
  gender <- as.matrix(data["gender"])
  n <- length(D)
  Ind_AA <- rep(0,n)
  Ind_Aa <- rep(0,n)
  Ind_AO <- rep(0,n)
  Ind_AA[gender==1 & G==2] <- 1
  Ind_Aa[gender==1 & G==1] <- 1
  Ind_AO[gender==0 & G==1] <- 1

  esta <- summary(stats::glm(D ~ X , family = stats::binomial(link = "logit")))$coefficients[,1]

  X <- cbind(rep(1,n),X)

  estf=function(esta,X){
    re=1/(1+exp(-X%*%esta))
    return(re)
  }
  estpen <- estf(esta,X)

  infora=function(estpen,X){
    l=dim(X)[2]
    Ia=matrix(0, nrow=l,ncol=l)
    for(i in 1:l){
      for(j in 1:l){
        Ia[i,j]=sum(X[,i]*X[,j]*(1-estpen)*estpen)
      }
    }
    return(Ia)
  }
  Ia=infora(estpen,X)

  inforb=function(estpen,z1,z2,X,Ind_AO,Ind_Aa,Ind_AA){
    G = 2*Ind_AA + z1*Ind_Aa + z2*Ind_AO
    Ib=sum(G*G*(1-estpen)*estpen)
    return(Ib)
  }
  Ib11=inforb(estpen,1,1,X,Ind_AO,Ind_Aa,Ind_AA)
  Ib02=inforb(estpen,0,2,X,Ind_AO,Ind_Aa,Ind_AA)
  Ib12=inforb(estpen,1,2,X,Ind_AO,Ind_Aa,Ind_AA)
  Ib22=inforb(estpen,2,2,X,Ind_AO,Ind_Aa,Ind_AA)


  inforba=function(estpen,z1,z2,X,Ind_AO,Ind_Aa,Ind_AA){
    l=dim(X)[2]
    G = 2*Ind_AA + z1*Ind_Aa + z2*Ind_AO
    Iba=NULL
    for(i in 1:l){
      Iba[i]=sum(X[,i]*G*(1-estpen)*estpen)
    }
    return(Iba)
  }
  Iba11=inforba(estpen,1,1,X,Ind_AO,Ind_Aa,Ind_AA)
  Iba02=inforba(estpen,0,2,X,Ind_AO,Ind_Aa,Ind_AA)
  Iba12=inforba(estpen,1,2,X,Ind_AO,Ind_Aa,Ind_AA)
  Iba22=inforba(estpen,2,2,X,Ind_AO,Ind_Aa,Ind_AA)

  stest=function(estpen,D,z1,z2,Ind_AO,Ind_Aa,Ind_AA,Ib,Iba,Ia){
    G <- 2*Ind_AA + z1*Ind_Aa + z2*Ind_AO
    score <- sum(G*(D-estpen))
    variance <- Ib-t(as.matrix(Iba))%*%solve(Ia)%*% as.matrix(Iba)
    re=score/sqrt(variance) ## test statistic
    #standerror <- sqrt(variance)
    return(re)
  }
  s11=stest(estpen,D,1,1,Ind_AO,Ind_Aa,Ind_AA,Ib11,Iba11,Ia)
  s02=stest(estpen,D,0,2,Ind_AO,Ind_Aa,Ind_AA,Ib02,Iba02,Ia)
  s12=stest(estpen,D,1,2,Ind_AO,Ind_Aa,Ind_AA,Ib12,Iba12,Ia)
  s22=stest(estpen,D,2,2,Ind_AO,Ind_Aa,Ind_AA,Ib22,Iba22,Ia)

  ####
  stestE=function(estpen,D,z1,z2,Ind_AO,Ind_Aa,Ind_AA,Ib,Iba,Ia){
    G <- 2*Ind_AA + z1*Ind_Aa + z2*Ind_AO
    score <- sum(G*(D-estpen))
    variance <- Ib-t(as.matrix(Iba))%*%solve(Ia)%*%as.matrix(Iba)
    #re=score/sqrt(variance) ## test statistic
    standerror <- sqrt(variance)/nrow(D) ## standard deviation/sample size
    return(standerror)
  }
  s11E=stestE(estpen,D,1,1,Ind_AO,Ind_Aa,Ind_AA,Ib11,Iba11,Ia)
  s02E=stestE(estpen,D,0,2,Ind_AO,Ind_Aa,Ind_AA,Ib02,Iba02,Ia)
  s12E=stestE(estpen,D,1,2,Ind_AO,Ind_Aa,Ind_AA,Ib12,Iba12,Ia)
  s22E=stestE(estpen,D,2,2,Ind_AO,Ind_Aa,Ind_AA,Ib22,Iba22,Ia)
  ####

  Iba11_02=rbind(Iba11,Iba02)
  Iba11_12=rbind(Iba11,Iba12)
  Iba11_22=rbind(Iba11,Iba22)
  Iba02_12=rbind(Iba02,Iba12)
  Iba02_22=rbind(Iba02,Iba22)
  Iba12_22=rbind(Iba12,Iba22)

  infor_dob=function(estpen,z11,z21,z12,z22,D,Ind_AO,Ind_Aa,Ind_AA){
    G1 = 2*Ind_AA + z11*Ind_Aa + z21*Ind_AO
    G2 = 2*Ind_AA + z12*Ind_Aa + z22*Ind_AO
    I1=sum(G1*G1*(1-estpen)*estpen)
    I12=sum(G1*G2*(1-estpen)*estpen)
    I2=sum(G2*G2*(1-estpen)*estpen)
    re=matrix(  c(I1,I12,I12,I2),nrow=2 )
    return(re)
  }
  inf11_02=infor_dob(estpen,1,1,0,2,D,Ind_AO,Ind_Aa,Ind_AA)
  inf11_12=infor_dob(estpen,1,1,1,2,D,Ind_AO,Ind_Aa,Ind_AA)
  inf11_22=infor_dob(estpen,1,1,2,2,D,Ind_AO,Ind_Aa,Ind_AA)
  inf02_12=infor_dob(estpen,0,2,1,2,D,Ind_AO,Ind_Aa,Ind_AA)
  inf02_22=infor_dob(estpen,0,2,2,2,D,Ind_AO,Ind_Aa,Ind_AA)
  inf12_22=infor_dob(estpen,1,2,2,2,D,Ind_AO,Ind_Aa,Ind_AA)

  C11_02=inf11_02 - as.matrix(Iba11_02)%*%solve(Ia)%*%t(as.matrix(Iba11_02))
  C11_12=inf11_12 - as.matrix(Iba11_12)%*%solve(Ia)%*%t(as.matrix(Iba11_12))
  C11_22=inf11_22 - as.matrix(Iba11_22)%*%solve(Ia)%*%t(as.matrix(Iba11_22))
  C02_12=inf02_12 - as.matrix(Iba02_12)%*%solve(Ia)%*%t(as.matrix(Iba02_12))
  C02_22=inf02_22 - as.matrix(Iba02_22)%*%solve(Ia)%*%t(as.matrix(Iba02_22))
  C12_22=inf12_22 - as.matrix(Iba12_22)%*%solve(Ia)%*%t(as.matrix(Iba12_22))

  RC11_02 <- C11_02[1,2] / sqrt( C11_02[1,1]* C11_02[2,2])
  RC11_12 <- C11_12[1,2] / sqrt( C11_12[1,1]* C11_12[2,2])
  RC11_22 <- C11_22[1,2] / sqrt( C11_22[1,1]* C11_22[2,2])
  RC02_12 <- C02_12[1,2] / sqrt( C02_12[1,1]* C02_12[2,2])
  RC02_22 <- C02_22[1,2] / sqrt( C02_22[1,1]* C02_22[2,2])
  RC12_22 <- C12_22[1,2] / sqrt( C12_22[1,1]* C12_22[2,2])

  vacov=matrix(c( 1,        RC11_02, RC11_12, RC11_22,
                  RC11_02,  1,       RC02_12, RC02_22,
                  RC11_12,  RC02_12, 1,       RC12_22,
                  RC11_22,  RC02_22, RC12_22, 1       ),
               ncol=4,
               byrow = TRUE)

  zmax1 <- max(abs(s11),abs(s02),abs(s12),abs(s22))
  zmax1E <- max(abs(s11E),abs(s02E),abs(s12E),abs(s22E))

  #rhombus formula
  p_rh <- function(zmax1,vacov,a,b,c,d){
    part1 <- 2*(stats::pnorm(zmax1)-pnorm(-zmax1)-1)
    l12 <- acos(vacov[a,b]);l23 <- acos(vacov[b,c]);l34<-acos(vacov[c,d])
    part2 <- stats::pnorm(l12*zmax1/2) + stats::pnorm((pi-l12)*zmax1/2) - 1
    part3 <- stats::pnorm(l23*zmax1/2) + stats::pnorm((pi-l23)*zmax1/2) - 1
    part4 <- stats::pnorm(l34*zmax1/2) + stats::pnorm((pi-l34)*zmax1/2) - 1
    p_rh <- part1 + 4*stats::dnorm(zmax1)/zmax1*(part2+part3+part4)
    return(p_rh)
  }
  p_series <- c(p_rh(zmax1,vacov,1,2,3,4),
                p_rh(zmax1,vacov,1,2,4,3),
                p_rh(zmax1,vacov,1,3,2,4),
                p_rh(zmax1,vacov,1,3,4,2),
                p_rh(zmax1,vacov,1,4,2,3),
                p_rh(zmax1,vacov,1,4,3,2),
                p_rh(zmax1,vacov,2,1,3,4),
                p_rh(zmax1,vacov,2,1,4,3),
                p_rh(zmax1,vacov,2,3,1,4),
                p_rh(zmax1,vacov,2,4,1,3),
                p_rh(zmax1,vacov,3,1,2,4),
                p_rh(zmax1,vacov,3,2,1,4))
  p_rh <- min(p_series,1)

  return(list("statictic"=zmax1,"standard-error"=zmax1E,
              "p-value"=p_rh))


}

## Function 62
## Updated in 3.0
XCMAX4_data1 <- function(Snp, genosnp, P, Samp){
  genosnp <- as.data.frame(genosnp)
  g <- genosnp[,Snp,drop = FALSE]
  colnames(g) <- c("Snp")
  g$IID <- Samp
  Pg <- merge(g,P,by = "IID")
  Data <- as.data.frame(Pg[,c(4,2,3)])
  Data$phenotype <- as.numeric(as.character(Data$phenotype))
  colnames(Data) <- c("D","genotype","gender")
  Data$gender <- as.numeric(as.character(Data$gender))
  # Replace the internal call to XCMAX4 with the standalone function
  invisible(Z <- XCMAX4(Data))
  invisible(Z <- XCMAX4(Data))
  pval <- Z$`p-value`
  statistic <- Z$statictic
  SE <- Z$`standard-error`
  result <- as.data.frame(cbind(Snp,statistic,pval,SE))
  colnames(result) <- c("SNP","STAT","P","SE")
  return(result)
}

## Function 63
## Added in 3.0
xcmaPara <- function(chunks, chunk, Snp, genosnp, P, Samp){
  ## Chunkfile create
  if (length(Snp)>=(chunks+chunk)){
    snp_names <- Snp[chunks:(chunks+chunk)]
  }else{
    snp_names <- Snp[chunks:length(Snp)]
  }
  # Replace the internal definition of XCMAX4_data1 with a call to the external function
  Result_pval1 <- data.table::rbindlist(lapply(snp_names, XCMAX4_data1, genosnp = genosnp, P = P, Samp = Samp))
  gc(reset = TRUE)
  return(Result_pval1)
}

## Function 64
## Updated in 3.0
## with covariate file
XCMAX4_data2 <- function(Snp, DataDir, genosnp, P, covarfile, Samp) {

  covarfile1 <- read.table(file = paste0(DataDir,"/",covarfile),
                           stringsAsFactors = FALSE,
                           header = TRUE)
  genosnp <- as.data.frame(genosnp)
  g <- genosnp[,Snp,drop = FALSE]
  colnames(g) <- c("Snp")
  g$IID <- Samp
  Pg <- merge(g,P,by = "IID")
  colnames(covarfile1[,-1])
  colnames(Pg)
  Data1 <- merge(Pg,covarfile1[,-1],by = "IID")[,-1]
  Data2 <- Data1[,c(3,1,2)]
  colnames(Data2) <- c("D","genotype","gender")
  Data3 <- Data1[,-c(1,2,3)]
  Data <- cbind(Data2,Data3)
  Datatest <- as.data.frame(Data)
  Z <- XCMAX4(Datatest)
  pval <- Z$`p-value`
  statistic <- Z$statictic
  SE <- Z$`standard-error`
  result <- as.data.frame(cbind(Snp,statistic,pval,SE))
  colnames(result) <- c("SNP","STAT","P","SE")
  return(result)
}

## Function 65
## Added in 3.0
xcmaParaCovar <- function(chunks, chunk, Snp, DataDir, genosnp, P, covarfile, Samp){

  #print(chunks)
  if (length(Snp)>=(chunks+chunk)){
    snp_names <- Snp[chunks:(chunks+chunk)]
  }else{
    snp_names <- Snp[chunks:length(Snp)]
  }

  ## with covariate file
  # Replace the internal definition of XCMAX4_data2 with a call to the external function
  Result_pval1 <- data.table::rbindlist(lapply(snp_names, XCMAX4_data2, DataDir = DataDir, genosnp = genosnp, P = P, covarfile = covarfile, Samp = Samp))
  #Result_pval1 <- data.table::rbindlist(lapply(snp_names,XCMAX4_data2,DataDir=DataDir,genosnp=genosnp,P=P,covarfile=covarfile,Samp=Samp))
  gc(reset = TRUE)
  return(Result_pval1)
}

## Function 66
## Added in 3.0
prepareData <- function(DataDir, ResultDir, finput) {
  ## Making separate plink files foor x chr and autosome and preparing input for XCMAX4().
  gdsfmt::showfile.gds(closeall=TRUE)
  # Separate X chromosome data
  invisible(sys::exec_wait(
    paste0(ResultDir, "/","./plink"),
    args = c("--bfile", paste0(DataDir, "/", finput), "--chr", 23, "--make-bed", "--out", paste0(ResultDir, "/","XPLINK"), "--silent"),
    std_out = FALSE, std_err = FALSE
  ))

  # Separate autosomal data
  invisible(sys::exec_wait(
    paste0(ResultDir, "/","./plink"),
    args = c("--bfile", paste0(DataDir, "/", finput), "--not-chr", 23, "--make-bed", "--out", paste0(ResultDir, "/","FilteredX"), "--silent"),
    std_out = FALSE, std_err = FALSE
  ))

  # File paths for the GDS conversion
  bedfile <- paste0(ResultDir, "/", "XPLINK.bed")
  famfile <- paste0(ResultDir, "/", "XPLINK.fam")
  bimfile <- paste0(ResultDir, "/", "XPLINK.bim")
  gdsfile <- paste0(ResultDir, "/", "IPCgeno.gdsX")

  # Convert PLINK files to GDS format
  gdsfileCreated <- SNPRelate::snpgdsBED2GDS(bedfile, famfile, bimfile, gdsfile)
  genofile <- SNPRelate::snpgdsOpen(gdsfileCreated)

  return(list(bedfile = bedfile, famfile = famfile, bimfile = bimfile, genofile = genofile))
}

## Function 67
## Added in 3.0
removePatternFiles <- function(ResultDir, patterns) {
  for (pattern in patterns) {
    files_to_remove <- list.files(ResultDir, pattern = pattern)
    full_paths <- file.path(ResultDir, files_to_remove)
    suppressWarnings(invisible(file.remove(full_paths)))
  }
}

## Function 68
## Updated in 3.0
XCMAFun <- function(DataDir, ResultDir, finput, standard_beta,
                    covarfile,sex, covartest,interaction,Inphenocov,plot.jpeg, plotname, snp_pval,
                    annotateTopSnp, suggestiveline = suggestiveline, genomewideline = genomewideline,
                    ncores = ncores){

  # Prepare data and get necessary files
  file_path <- paste0(ResultDir,"/sink_file.txt")

  # Create an empty file (if it doesn't exist)
  if (!file.exists(file_path)) {
    file.create(file_path)
  }

  # Redirect output to the specified file
  sink(file_path)

  dataFiles <- prepareData(DataDir, ResultDir, finput)

  # Reset the sink to stop redirecting the output to the file
  if (sink.number() > 0) {
    sink()  # This line resets the output redirection
  }

  # Extract additional data from the GDS file
  P <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(dataFiles$genofile, "sample.annot"))
  Snp <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(dataFiles$genofile, "snp.id"))
  Samp <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(dataFiles$genofile, "sample.id"))
  G <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(dataFiles$genofile, "genotype"))

  P[P$sex == "M",1] <- 0
  P[P$sex == "F",1] <- 1
  P$sex <- as.numeric(as.character(P$sex))
  P["phenotype"][P["phenotype"] == 1] <- 0
  P["phenotype"][P["phenotype"] == 2] <- 1
  P$IID <- Samp
  genosnp <- data.table::as.data.table(SNPRelate::snpgdsGetGeno(dataFiles$genofile, snp.id = Snp,verbose = FALSE))
  colnames(genosnp)<- Snp

  if (is.null(covarfile)) {

    if (ncores == 0){

      Result_pval <- data.table::rbindlist(lapply(Snp,XCMAX4_data1,genosnp=genosnp,P=P,Samp=Samp))

    }else{

      chunk <- (round(length(Snp)/ncores)+1)/3
      chunks <- round(seq(1, length(Snp), by = chunk),0)

      #### Parallel computation
      gc(reset = TRUE)
      print("Parallel computation is in progress --------->")
      cl <- parallel::makeCluster(mc <- getOption("cl.cores", ncores),type = "FORK")

      invisible(parallel::clusterEvalQ(cl, library(data.table)))
      invisible(parallel::clusterEvalQ(cl, library(parallel)))
      parallel::clusterExport(cl=cl,NULL,envir=environment())

      Result_pval <- data.table::rbindlist(parallel::parLapply(cl, chunks, xcmaPara, chunk = chunk, Snp = Snp, genosnp = genosnp, P = P, Samp = Samp))
      parallel::stopCluster(cl)
      gc(reset = TRUE)
    }

    Result_pval$CHR <- "23"
    Result_pval$TEST <- "GWAScxci"
    Result_pval$BETA <- "NA"
    Result_pval$L95 <- "NA"
    Result_pval$U95 <- "NA"

  } else{

    DataDir=DataDir
    covarfile=covarfile
    print("Running GWAScxci model for X chromosome.")

    if (ncores == 0){

      Result_pval <- data.table::rbindlist(lapply(Snp,XCMAX4_data2,DataDir=DataDir,genosnp=genosnp,P=P,Samp = Samp,covarfile=covarfile))

    }else{

      chunk <- (round(length(Snp)/ncores)+1)/3
      chunks <- round(seq(1, length(Snp), by = chunk),0)

      #### Parallel computation
      gc(reset = TRUE)
      print("Parallel computation is in progress --------->")
      cl <- parallel::makeCluster(mc <- getOption("cl.cores", ncores),type = "FORK")

      invisible(parallel::clusterEvalQ(cl, library(data.table)))
      invisible(parallel::clusterEvalQ(cl, library(parallel)))
      invisible(parallel::clusterEvalQ(cl, library(SNPRelate)))
      parallel::clusterExport(cl=cl,NULL,envir=environment())

      Result_pval <- data.table::rbindlist(parallel::parLapply(cl, chunks, xcmaParaCovar, chunk = chunk, Snp = Snp, DataDir = DataDir, genosnp = genosnp, P = P, covarfile = covarfile, Samp = Samp))
      parallel::stopCluster(cl)
      gc(reset = TRUE)

    }

    Result_pval$CHR <- "23"
    Result_pval$TEST <- "GWAScxci"
    Result_pval$BETA <- "NA"

    Result_pval$L95 <- "NA"
    Result_pval$U95 <- "NA"
    Xchr_Result_pval <- Result_pval
    save(Xchr_Result_pval, file = paste0(ResultDir,"/Xchr_Result_pval.Rda"))

    rm(Xchr_Result_pval)
    gc(reset = TRUE)
  }

  ## Performing Autosome WAS
  print("Running GWAScxci model for autosome.")
  XWAS_2 <- autoFun(DataDir = DataDir, ResultDir = ResultDir, finput = "FilteredX", sex = sex, standard_beta = standard_beta,
                    covarfile = covarfile, interaction =interaction,
                    covartest = covartest, Inphenocov = Inphenocov, ncores = ncores)

  ## Extra run for X related SNP for gathering other information
  invisible(sys::exec_wait(
    paste0(ResultDir,"/","./plink"),
    args = c(
      "--bfile",
      paste0(ResultDir,"/","XPLINK"),
      "--logistic",
      "beta",
      "intercept", "--ci", 0.95,
      "--out",
      paste0(ResultDir,"/","XChrRun"),
      "--silent"
    ),
    std_out = FALSE,
    std_err = FALSE
  ))

  if (file.exists(paste0(ResultDir,"/","XChrRun.assoc.logistic")) == TRUE) {
    XWAS_3 <-
      read.table(file = paste0(ResultDir,"/","XChrRun.assoc.logistic"),
                 stringsAsFactors = FALSE,
                 header = TRUE)

    XWAS_4 <- XWAS_3[,c(2:4,6)]

    XWAS_5 <- unique(merge(Result_pval,XWAS_4,by = "SNP"))
    XWAS_6 <- XWAS_5[,c(5,1,10,11,6,12,7,4,8,9,2,3)]
    XWAS <- rbind(XWAS_2,XWAS_6)
    XWAS$CHR <- as.integer(as.character(XWAS$CHR))
    XWAS$P <- as.numeric(as.character(XWAS$P ))
    ## Plots
    XWAS_ADD <- na.omit(XWAS[XWAS$TEST == "ADD" | XWAS$TEST == "GWAScxci",c("SNP","CHR","BP","P")])
    XWAS_ADD$P <- as.numeric(as.character(XWAS_ADD$P))
    # From p-values, calculate chi-squared statistic
    chisq <- qchisq(1-na.omit(XWAS_ADD$P),1)
    lamdaGC <- median(chisq)/qchisq(0.5,1)

    XWAS_ADD_X <- na.omit(XWAS_ADD[XWAS_ADD$CHR == 23,])
    chisq1 <- qchisq(1-na.omit(XWAS_ADD_X$P),1)
    lamdaGC1 <- median(chisq1)/qchisq(0.5,1)

    createPlotsXCMA(XWAS_ADD, XWAS_ADD_X, ResultDir, plotname, snp_pval, lamdaGC, lamdaGC1, suggestiveline, genomewideline, annotateTopSnp, plot.jpeg)

    GXWAS <- XWAS
    save(GXWAS, file = paste0(ResultDir,"/GXWAS_XCGA.Rda"))
    print(paste0("A dataframe named GXWAS_XCGA.Rda is saved in ",ResultDir))

  } else if (file.exists(paste0(ResultDir,"/","XChrRun.assoc.logistic")) == FALSE) {
    print("XWAS cannot be performed. Check the XChrRun.log file for checking the error.")

  }

  return(GXWAS)
}


## Function 69
ComputeLD <- function(DataDir,ResultDir,finput, ByCHR = FALSE, CHRnum = NULL, r2_LD ){


  if (ByCHR == FALSE){
    chr = NULL
    CHRnum = NULL
  }else{
    chr = "--chr"
    CHRnum = CHRnum
  }
  invisible(sys::exec_wait(
    paste0(ResultDir,"/","./plink"),
    args = c(
      "--bfile",
      paste0(DataDir,"/",finput),
      chr,CHRnum,
      "--r2",
      "--ld-window-r2",r2_LD,
      "--out",
      paste0(ResultDir,"/","snpcorr"),
      "--silent"
    ),
    std_out = FALSE,
    std_err = FALSE
  ))
  snpld <- read.table(paste0(ResultDir,"/snpcorr.ld"),header = T)
  return(snpld)
}

## Function 70
ComputeLDSC <- function(snpld, test.df_beta, ncores, LDSC_blocks){

  snps1 <- data.table::data.table(intersect(unique(snpld[["SNP_A"]]), unique(test.df_beta[["rsid"]])))
  colnames(snps1) <- "SNP"
  snps2 <- data.table::data.table(intersect(unique(snpld[["SNP_B"]]), unique(test.df_beta[["rsid"]])))
  colnames(snps2) <- "SNP"
  snps <- data.table::data.table(intersect(unique(snps1[["SNP"]]),unique(snps2[["SNP"]])))
  colnames(snps) <- "SNP"
  snpld1 <- as.data.frame(snpld)
  snpld1 <- merge(snpld,snps, by.x = "SNP_A",by.y="SNP")
  snpld1 <- merge(snpld1,snps, by.x = "SNP_B",by.y="SNP")

  snpld1 <- snpld1[,c(2,1,7)]
  finalsnp <- data.table::as.data.table(snpld1[["SNP_A"]])
  colnames(finalsnp) <- "SNP"
  test.df_beta <- data.table::as.data.table(test.df_beta)
  test.df_beta <- unique(merge(test.df_beta,finalsnp,by.x = "rsid",by.y = "SNP"))

  snpld1 = transform(snpld1, SNP_A= factor(SNP_A),SNP_B = factor(SNP_B))
  test.corr0 = Matrix::sparseMatrix(as.integer(snpld1$SNP_A), as.integer(snpld1$SNP_B), x = snpld1$R2)
  colnames(test.corr0) = levels(snpld1$SNP_B)
  rownames(test.corr0) = levels(snpld1$SNP_A)

  test.df_beta1 <- test.df_beta[,c("beta","beta_se","n_eff")]
  result <- bigsnpr::snp_ldsc2(test.corr0, test.df_beta1, blocks = LDSC_blocks, intercept = NULL, ncores = ncores)
  ## If there is negative h2, then make it 0. Discuss.
  return(result)

}

## Function 71
setupGCTA <- function(wdir) {
  # Helper function to download and unzip files
  downloadAndUnzip <- function(file_url, dest_file, wdir) {
    utils::download.file(destfile = dest_file, file_url, quiet = TRUE)
    utils::unzip(dest_file, exdir = wdir, junkpaths = TRUE)
  }

  # Helper function to remove files
  removeFiles <- function(files, wdir) {
    full_paths <- file.path(wdir, files)
    invisible(sapply(full_paths, file.remove))
  }

  # Get operating system information
  OS <- Sys.info()['sysname']
  os_specific_data <- list(
    Linux = list(
      file_url = "https://github.com/boseb/bose_binaries/raw/main/gcta-1.94.1-linux-kernel-3-x86_64.zip",
      file_name = "gcta-1.94.1-linux-kernel-3-x86_64.zip",
      exec = "gcta-1.94.1",
      remove = c("MIT_License.txt", "README.txt")
    ),
    Windows = list(
      file_url = "https://figshare.com/ndownloader/files/42229122",
      file_name = "gcta-1.94.1",
      exec = "gcta-1.94.1",
      remove = c()
      # Add additional setup steps for Windows if needed
    ),
    macOS = list(
      file_url = "https://github.com/boseb/bose_binaries/blob/main/gcta-1.94.1-macOS-x86_64.zip",
      file_name = "gcta-1.94.1-macOS-x86_64.zip",
      exec = "gcta-1.94.1",
      remove = c("MIT_License.txt", "README.txt")
    )
  )

  if (!OS %in% names(os_specific_data)) {
    stop("Unsupported operating system.")
  }

  # Perform the setup
  os_data <- os_specific_data[[OS]]
  dest_file <- file.path(wdir, os_data$file_name)

  tryCatch({
    downloadAndUnzip(os_data$file_url, dest_file, wdir)
    Sys.chmod(file.path(wdir, os_data$exec), mode = "0777", use_umask = TRUE)
    removeFiles(c(os_data$file_name, os_data$remove), wdir)
    print("GCTA setup complete. This program is currently only supported in Linux.")
  }, error = function(e) {
    cat("An error occurred during GCTA setup:", e$message, "\n")
  })
}

## Function 72
## Added in 3.0
executeGCTA <- function(ResultDir, args) {
  tryCatch({
    # Determine the destinations for standard output and standard error based on the operating system
    stdout_dest <- ifelse(.Platform$OS.type == "windows", "NUL", "/dev/null")
    stderr_dest <- stdout_dest  # Redirecting stderr to the same null device

    file_path <- paste0(ResultDir,"/sink_file.txt")

    # Create an empty file (if it doesn't exist)
    if (!file.exists(file_path)) {
      file.create(file_path)
    }

    # Redirect output to the specified file
    sink(file_path)

    # Execute GCTA command

    invisible(sys::exec_wait(
      file.path(ResultDir, "gcta-1.94.1"),  # Path to the GCTA executable
      args = args,                          # Arguments for the GCTA command
      std_out = stdout_dest,                # Standard output redirection
      std_err = stderr_dest                 # Standard error redirection
    ))
  }, error = function(e) {
    stop("An error occurred while executing GCTA: ", e$message)
  })

  # Reset the sink to stop redirecting the output to the file
  if (sink.number() > 0) {
    sink()  # This line resets the output redirection
  }
}

## Function 73
ComputeGRMauto <- function(DataDir, ResultDir, finput,partGRM, nGRM, cripticut,  minMAF = NULL, maxMAF = NULL, ByCHR = FALSE, CHRnum = NULL){

  if (ByCHR == FALSE){
    chr = NULL
    CHRnum = NULL
    autosome = "--autosome"
  }else{
    chr = "--chr"
    CHRnum = CHRnum
    autosome = NULL
  }

  if (is.null(maxMAF)){
    maxmaf <- NULL
    maxmafval <- NULL
  }else{
    maxmaf <- "--max-maf"
    maxmafval <- maxMAF
  }

  if (is.null(minMAF)){
    minmaf <- NULL
    minmafval <- NULL
    maxmaf <- NULL
    maxmafval <- NULL
  }else{
    minmaf <- "--maf"
    minmafval <- minMAF
  }

  if (cripticut == 0){
    grmcutoff <- NULL
    crip <- NULL

  }else{
    grmcutoff <- "--grm-cutoff"
    crip <- cripticut
  }

  if (partGRM == FALSE){
    args <- c(
      "--bfile", paste0(DataDir, "/", finput),
      autosome,
      chr, CHRnum,
      minmaf, minmafval,
      maxmaf, maxmafval,
      grmcutoff, crip,
      "--make-grm",
      "--out", paste0(ResultDir, "/test")
    )
    executeGCTA(ResultDir, args)

  }else {

    partGRMfun <- function(i){

      args <- c(
        "--bfile", paste0(DataDir, "/", finput),
        autosome,
        chr, CHRnum,
        minmaf, minmafval,
        maxmaf, maxmafval,
        grmcutoff,
        crip,
        "--make-grm-part",
        nGRM,
        i, # Part number
        "--thread-num", 5,
        "--out", paste0(ResultDir, "/test")
      )

      # Execute GCTA with the specified arguments
      executeGCTA(ResultDir, args)
    }

    i = 1:nGRM
    lapply(i,partGRMfun)

    OS <- Sys.info()['sysname']

    if (OS == "Linux" | OS == "Mac") {

      system(paste0("cat ",ResultDir,"/","test.part_",nGRM,"_*.grm.id > ",ResultDir,"/test.grm.id"))
      system(paste0("cat ",ResultDir,"/","test.part_",nGRM,"_*.grm.bin > ",ResultDir,"/test.grm.bin"))
      system(paste0("cat ",ResultDir,"/","test.part_",nGRM,"_*.grm.N.bin > ",ResultDir,"/test.grm.N.bin"))


    }else if (OS == "Windows"){

      # copy /b test.part_3_*.grm.id test.grm.id
      # copy /b test.part_3_*.grm.bin test.grm.bin
      # copy /b test.part_3_*.grm.N.bin test.grm.N.bin
    }

  }

  ## Compute GRM: GRM is calculated using the equation sum{[(xij - 2pi)*(xik - 2pi)] / [2pi(1-pi)]} as described in Yang et al. 2010 Nat Genet.


}

## Function 74
ComputeGRMX <- function(DataDir, ResultDir, finput, partGRM, nGRM, minMAF = NULL, maxMAF = NULL){

  ## Removing PAR regions and adding it to autosomes, PAR = chr 24 in .bim file.


  if (is.null(maxMAF)){
    maxmaf <- NULL
    maxmafval <- NULL
  }else{
    maxmaf <- "--max-maf"
    maxmafval <- maxMAF
  }

  if (is.null(minMAF)){
    minmaf <- NULL
    minmafval <- NULL
    maxmaf <- NULL
    maxmafval <- NULL
  }else{
    minmaf <- "--maf"
    minmafval <- minMAF
  }
  if (partGRM == FALSE){

    args <- c(
      "--bfile", paste0(DataDir, "/", finput),
      minmaf, minmafval,
      maxmaf, maxmafval,
      "--make-grm-xchr",
      "--out", paste0(ResultDir, "/xtest")
    )

    # Execute GCTA with the specified arguments for X chromosome GRM
    executeGCTA(ResultDir, args)


  }else {

    partGRM <- function(i){

      args <- c(
        "--bfile", paste0(DataDir, "/", finput),
        minmaf, minmafval,
        maxmaf, maxmafval,
        "--make-grm-xchr-part",
        nGRM,
        i,  # Part number for X chromosome GRM
        "--thread-num", 5,
        "--out", paste0(ResultDir, "/xtest")
      )

      # Execute GCTA with the specified arguments for X chromosome part
      executeGCTA(ResultDir, args)

    }
    i = 1:nGRM
    lapply(i,partGRM)

    OS <- Sys.info()['sysname']
    if (OS == "Linux" | OS == "Mac") {

      system(paste0("cat ",ResultDir,"/","xtest.part_",nGRM,"_*.grm.id > ",ResultDir,"/xtest.grm.id"))
      system(paste0("cat ",ResultDir,"/","xtest.part_",nGRM,"_*.grm.bin > ",ResultDir,"/xtest.grm.bin"))
      system(paste0("cat ",ResultDir,"/","xtest.part_",nGRM,"_*.grm.N.bin > ",ResultDir,"/xtest.grm.N.bin"))


    }else if (OS == "Windows"){

      # copy /b test.part_3_*.grm.id test.grm.id
      # copy /b test.part_3_*.grm.bin test.grm.bin
      # copy /b test.part_3_*.grm.N.bin test.grm.N.bin
    }
  }

  ## Compute GRM: GRM is calculated using the equation sum{[(xij - 2pi)*(xik - 2pi)] / [2pi(1-pi)]} as described in Yang et al. 2010 Nat Genet.


}


## Function 75
ComputeREMLone <- function(DataDir, ResultDir, REMLalgo = c(0,1,2), nitr = 100, phenofile = NULL, cat_covarfile = NULL, quant_covarfile = NULL,
                           prevalance = 0.01, chr ,grmfile, ncores){

  if (is.null(phenofile)){
    pheno = NULL
    phenofile = NULL
  }else{
    pheno = "--pheno"
    phenofile = phenofile
  }
  if (prevalance == 0){
    preval = NULL
    prevalval = NULL
  }else{
    preval = "--prevalence"
    prevalval = prevalance
  }

  if (is.null(cat_covarfile)){
    catcovar = NULL
    catcovarval = NULL
  }else{
    catcovar = "--ccovar"
    catcovarval = paste0(DataDir,"/",cat_covarfile)
  }

  if (is.null(quant_covarfile)){
    quantcovar = NULL
    quantcovarval = NULL
  }else{
    quantcovar = "--qcovar"
    quantcovarval = paste0(DataDir,"/",quant_covarfile)
  }

  # Constructing the arguments for the GCTA command
  args <- c(
    "--reml",
    "--reml-alg", REMLalgo,
    "--reml-maxit", nitr,
    "--grm", paste0(ResultDir, "/", grmfile),
    pheno, paste0(ResultDir, "/", phenofile),
    preval, prevalval,
    catcovar, catcovarval,
    quantcovar, quantcovarval,
    "--reml-no-lrt",
    "--reml-no-constrain",
    "--thread-num", ncores,
    "--out", paste0(ResultDir, "/", chr, "test_reml")
  )

  # Executing GCTA with the specified arguments
  executeGCTA(ResultDir, args)


  if (file.exists(paste0(ResultDir,"/",chr,"test_reml.hsq"))){

    resultREML <- na.omit(data.table::fread(paste0(ResultDir,"/",chr,"test_reml.hsq"),fill = TRUE))
    return(resultREML)
  }else{
    print(grep("Error", readLines(paste0(ResultDir,"/", chr,"test_reml.log")), value = TRUE))
    print("Convergence issue occurs, please check the models, use byCHR = TRUE, check different options, SNP partitioning, or quality of the  data")
    x <- data.frame(NA,NA,NA)
    colnames(x) <- c("Source", "Variance","SE")
    return()
  }

  #gcta64  --reml  --grm test  --pheno test_cc.phen  --prevalence 0.01  --qcovar test_10PCs.txt  --out test_cc



}

## Function 76
ComputeREMLmulti <- function(DataDir, ResultDir, REMLalgo = c(0,1,2), nitr = 100, phenofile = NULL, GE = FALSE, cat_covarfile = NULL, quant_covarfile = NULL,
                             prevalance = 0.01, grmfile = "multi_GRMs.txt", ncores){

  if (is.null(phenofile)){
    pheno = NULL
    phenofile = NULL
  }else{
    pheno = "--pheno"
    phenofile = phenofile
  }
  if (prevalance == 0){
    preval = NULL
    prevalval = NULL
  }else{
    preval = "--prevalence"
    prevalval = prevalance
  }

  if (is.null(cat_covarfile)){
    catcovar = NULL
    catcovarval = NULL
  }else{
    catcovar = "--ccovar"
    catcovarval = paste0(DataDir,"/",cat_covarfile)
  }

  if (is.null(quant_covarfile)){
    quantcovar = NULL
    quantcovarval = NULL
  }else{
    quantcovar = "--qcovar"
    quantcovarval = paste0(DataDir,"/",quant_covarfile)
  }

  ## Create multi_GRMs.txt in ResultDir
  fileConn<-file(paste0(ResultDir,"/multi_GRMs.txt"))
  writeLines(c(paste0(ResultDir,"/test"),paste0(ResultDir,"/xtest")), fileConn)
  close(fileConn)

  invisible(sys::exec_wait(
    paste0(ResultDir,"/","./gcta-1.94.1"),
    args = c(
      "--reml",
      "--reml-alg",REMLalgo,
      "--reml-maxit",nitr,
      "--mgrm",
      paste0(ResultDir,"/","multi_GRMs.txt"),
      pheno, paste0(ResultDir,"/",phenofile),
      preval, prevalval,
      catcovar, catcovarval,
      quantcovar, quantcovarval,
      "--reml-no-lrt",
      "--reml-no-constrain",
      "--thread-num",ncores,
      "--out",
      paste0(ResultDir,"/","test_reml")
    ),
    std_out = FALSE,
    std_err = FALSE
  ))

  if (file.exists(paste0(ResultDir,"/test_reml.hsq"))){

    resultREML <- na.omit(data.table::fread(paste0(ResultDir,"/test_reml.hsq"),fill = TRUE))
    return(resultREML)
  }else{
    print(grep("Error", readLines(paste0(ResultDir,"/test_reml.log")), value = TRUE))
    print("Convergence issue occurs, please check the models, use byCHR = TRUE, check different options, SNP partitioning or quality of the data")
    x <- data.frame(NA,NA,NA)
    colnames(x) <- c("Source", "Variance","SE")
    return()
  }

  #gcta64  --reml  --grm test  --pheno test_cc.phen  --prevalence 0.01  --qcovar test_10PCs.txt  --out test_cc

}


## Function 77
## Getting chromosome length
#hg = "hg19"
ChrLength <- function(hg){
  x <- GenomeInfoDb::getChromInfoFromUCSC(paste0(hg))
  x$size_mb <- x$size/1000000
  x_mb <- x[1:26,c(1,5)]
  #library(stringr)
  x_mb$chrom <- 1:26
  return(x_mb)

}

## Function 78
GettingGene <- function(gene_file, gene_range, SNP_bimfile, finput){
  genes <- read.table(gene_file)
  colnames(genes) <- c(c("gene_name", "X", "chr", "Y", "start", "end"))
  genes$up_Mb <- genes$start - gene_range
  genes$down_Mb <- genes$end + gene_range
  genes.gr <- GenomicRanges::makeGRangesFromDataFrame(genes, keep.extra.columns = T)

  suppressWarnings(SNPfile <-  read.table(
    file = paste0(SNP_bimfile, ".bim"),
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
}

## Function 79
ChrwiseLDprun <- function(DataDir,ResultDir,finput,chromosome, highLD_regions, IndepSNP_window_size,
                          IndepSNP_step_size,IndepSNP_r2_threshold){

  if (is.null(highLD_regions)){
    highLD_regions <- NULL
    excludev <- NULL
  }else{
    write.table(highLD_regions, file = paste0(ResultDir,"/","highLD_regions_temp"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    #highLD_regions <- paste0(DataDir,"/",highLD_regions)
    highLD_regions <- paste0(ResultDir,"/","highLD_regions_temp")
    excludev <- "--exclude"
  }

  invisible(sys::exec_wait(
    paste0(ResultDir,"/","./plink"),
    args = c("--bfile",
             paste0(DataDir,"/",finput),
             "--chr",chromosome,
             #"--exclude",
             excludev,
             #paste0(DataDir,"/",highLD_regions),
             highLD_regions,
             "--indep-pairwise",
             IndepSNP_window_size,
             IndepSNP_step_size,
             IndepSNP_r2_threshold,
             "--out",
             paste0(ResultDir,"/","LDsnp"),
             "--silent"
    ),
    std_out = FALSE,
    std_err = FALSE
  ))

  invisible(sys::exec_wait(
    paste0(ResultDir,"/","./plink"),
    args = c("--bfile",
             paste0(DataDir,"/",finput),
             "--extract",
             paste0(ResultDir,"/","LDsnp.prune.in"),
             "--allow-no-sex",
             "--make-bed",
             "--out",
             paste0(ResultDir,"/","LDfiltered"),
             "--silent"
    ),
    std_out = FALSE,
    std_err = FALSE
  ))

  ftemp <- list.files(paste0(ResultDir,"/"),pattern = "_temp")
  invisible(file.remove(paste0(ResultDir,"/",ftemp)))
}
## Function 80
GeneProtein <- function(ResultDir,hg,chromosome){ ## Automatically using HG data from extdata
  if (hg == "hg19"){
    #DataDir <- system.file("extdata", package = "GXwasR")

    #Downloading HG19 data from figshare to ResultDir
    utils::download.file(destfile = paste0(ResultDir,"/HumanGenome19info.txt"),
                         "https://figshare.com/ndownloader/files/42118098", quiet = TRUE,)

    genes <- unique(read.table(file = paste0(ResultDir,"/","HumanGenome19info.txt"))[,c(1,7,8)])#Using hgnc_name

    genes$Chrom <- stringr::str_sub(genes$Chrom,4)
    genes1 <- unique(genes[genes$Chrom==chromosome,2,drop=FALSE])
    no.of.genes <- nrow(genes1)
    proteins <- genes[genes$gene_biotype=="protein_coding",]
    proteins1 <- unique(proteins[proteins$Chrom==chromosome,2,drop=FALSE])
    no.of.proteins <- nrow(proteins1)
    GP <- data.table::as.data.table(cbind(no.of.genes,no.of.proteins))
    return(GP)
  }else{
    utils::download.file(destfile = paste0(ResultDir,"/HumanGenome38info.txt"),
                         "https://figshare.com/ndownloader/files/42118242", quiet = TRUE,)

    genes <- unique(read.table(file = paste0("inst/extdata/HumanGenome38info.txt"))[,c(1,7,8)])#Using hgnc_name
    genes$Chrom <- stringr::str_sub(genes$Chrom,4)
    genes1 <- unique(genes[genes$Chrom==chromosome,2,drop=FALSE])
    no.of.genes <- nrow(genes1)
    proteins <- genes[genes$gene_biotype=="protein_coding",]
    proteins1 <- unique(proteins[proteins$Chrom==chromosome,2,drop=FALSE])
    no.of.proteins <- nrow(proteins1)
    GP <- data.table::as.data.table(cbind(no.of.genes,no.of.proteins))
    return(GP)
  }
  ftemp <- list.files(paste0(ResultDir,"/"),pattern = "info.txt")
  invisible(file.remove(paste0(ResultDir,"/",ftemp)))
}

## Function 81
PlotHeritability <- function(Hdata,miMAF,maMAF){
  #create plot with regression line, regression equation, Pearson correlation and p-value.
  Hdata <- na.omit(Hdata)
  p1<- ggplot2::ggplot(data=Hdata, ggplot2::aes(x=Hdata$size_mb, y=Hdata$Variance)) +
    ggplot2::labs(title="Chromosome-wise heritability",
                  x = "Chromosome length (mb)", y = "Heritability") +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, level=0.95) +      ## For CI
    ggplot2::geom_point(size = 5, shape = 21, fill = "white") +
    ggplot2::geom_text(label = Hdata$chromosome,vjust = 0.4, size = 3)+
    ggpubr::stat_regline_equation(label.x=median(Hdata$size_mb), label.y= max(Hdata$Variance) + max(Hdata$Variance)/2) +
    ggpubr::stat_cor(method = "pearson", label.x = median(Hdata$size_mb), label.y = max(Hdata$Variance) + max(Hdata$Variance)/3)+
    ggplot2::xlim(min(Hdata$size_mb), max(Hdata$size_mb))

  p2 <- ggplot2::ggplot(data=Hdata, ggplot2::aes(x=Hdata$snp_proportion, y=Hdata$Variance)) +
    ggplot2::labs(title="SNP proportion per \nchromosome vs heritability",
                  x = "snp_proportion", y = "Heritability") +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x) +
    ggplot2::geom_point(size = 5, shape = 21, fill = "white") +
    ggplot2::geom_text(label = Hdata$chromosome,vjust = 0.4, size = 3)+
    ggpubr::stat_regline_equation(label.x= median(Hdata$snp_proportion), label.y= max(Hdata$Variance) + max(Hdata$Variance)/2) +
    ggpubr::stat_cor(method = "pearson", label.x = median(Hdata$snp_proportion), label.y = max(Hdata$Variance) + max(Hdata$Variance)/3)+
    ggplot2::xlim(min(Hdata$snp_proportion), max(Hdata$snp_proportion))


  p3 <- ggplot2::ggplot(data=Hdata, ggplot2::aes(x=Hdata$no.of.genes, y=Hdata$Variance)) +
    ggplot2::labs(title="Number of genes per \nchromosome vs heritability",
                  x = "no.of.genes", y = "Heritability") +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x) +
    ggplot2::geom_point(size = 5, shape = 21, fill = "white") +
    ggplot2::geom_text(label = Hdata$chromosome,vjust = 0.4, size = 3)+
    ggpubr::stat_regline_equation(label.x= median(Hdata$no.of.genes), label.y= max(Hdata$Variance) + max(Hdata$Variance)/2) +
    ggpubr::stat_cor(method = "pearson", label.x = median(Hdata$no.of.genes), label.y = max(Hdata$Variance) + max(Hdata$Variance)/3)+
    ggplot2::xlim(min(Hdata$no.of.genes), max(Hdata$no.of.genes))

  p4 <- ggplot2::ggplot(data=Hdata, ggplot2::aes(x=Hdata$no.of.proteins, y=Hdata$Variance)) +
    ggplot2::labs(title="Number of proteins per \nchromosome vs heritability",
                  x = "no.of.proteins", y = "Heritability") +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x) +
    ggplot2::geom_point(size = 5, shape = 21, fill = "white") +
    ggplot2::geom_text(label = Hdata$chromosome,vjust = 0.4, size = 3)+
    ggpubr::stat_regline_equation(label.x= median(Hdata$no.of.proteins), label.y= max(Hdata$Variance) + max(Hdata$Variance)/2) +
    ggpubr::stat_cor(method = "pearson", label.x = median(Hdata$no.of.proteins), label.y = max(Hdata$Variance) + max(Hdata$Variance)/3)+
    ggplot2::xlim(min(Hdata$no.of.proteins), max(Hdata$no.of.proteins))

  plot1 <- ggpubr::ggarrange(p1, p2,p3,p4,
                             labels = c("A", "B","C","D"),
                             ncol = 2, nrow = 2)

  print(ggpubr::annotate_figure(plot1, top = ggpubr::text_grob(paste0(miMAF,",",maMAF),
                                                               color = "red", face = "bold", size = 10)))
}

## Function 82
## Added in 3.0
computeMAFRange <- function(DataDir, ResultDir, finput, minMAF, maxMAF) {
  if (is.null(minMAF) | is.null(maxMAF)) {
    # Compute min and max MAF using Plink
    plink_exec_path <- paste0(ResultDir, "/plink")
    args <- c("--bfile", paste0(DataDir, "/", finput), "--freq", "--out", paste0(ResultDir, "/MAF"), "--silent")

    invisible(sys::exec_wait(plink_exec_path, args = args, std_out = FALSE, std_err = FALSE))

    MAF <- na.omit(read.table(paste0(ResultDir, "/MAF.frq"), header = TRUE))
    minmaf <- min(MAF$MAF)
    maxmaf <- max(MAF$MAF)
    miMAF <- paste0("MAF = [", minmaf)
    maMAF <- paste0(maxmaf, "]")
  } else {
    miMAF <- paste0("MAF < ", minMAF)
    maMAF <- paste0(maxMAF, " >")
  }

  return(list(miMAF = miMAF, maMAF = maMAF))
}

## Function 83
## Added in 3.0
processGWASData <- function(DataDir, finput, summarystat) {
  # Construct the file path for the GWAS bed file
  gxwas_bedfile <- paste0(DataDir, "/", finput, ".bed")

  # Read the bedfile and store the data in a temporary directory
  rds <- bigsnpr::snp_readBed(gxwas_bedfile, backingfile = tempfile())

  # Load the data from backing files
  test.bigSNP <- bigsnpr::snp_attach(rds)
  test.G <- test.bigSNP$genotypes
  y <- test.bigSNP$fam$affection - 1
  test.map <- setNames(test.bigSNP$map[-c(3)], c("chr", "rsid", "pos", "a1", "a0"))

  # Merge summary statistics with the test map
  summarystat1 <- merge(summarystat, test.map, by = c("chr", "rsid", "pos", "a1"))
  summarystat2 <- summarystat1[, c(1:4, 8, 5, 6, 7)]

  # Match summary statistics to test map and omit NA values
  test.df_beta <- bigsnpr::snp_match(summarystat2, test.map, join_by_pos = FALSE)
  test.df_beta <- na.omit(test.df_beta)

  # Return a list of processed data
  return(list(test.bigSNP = test.bigSNP, test.G = test.G, y = y, test.map = test.map, test.df_beta = test.df_beta))
}

## Function 84
## Added in 3.0
processLDSCModel <- function(DataDir, ResultDir, finput, summarystat, byCHR, r2_LD, LDSC_blocks, ncores, prevalance, PlotIndepSNP, highLD_regions, IndepSNP_window_size, IndepSNP_step_size, IndepSNP_r2_threshold, hg, miMAF, maMAF) {
  processed_data <- processGWASData(DataDir, finput, summarystat)
  test.df_beta <- processed_data$test.df_beta

  if (byCHR == FALSE) {
    snpld <- ComputeLD(DataDir = DataDir, ResultDir = ResultDir, finput = finput, ByCHR = FALSE, CHRnum = NULL, r2_LD = r2_LD)
    result <- ComputeLDSC(snpld = snpld, test.df_beta = test.df_beta, ncores = ncores, LDSC_blocks = LDSC_blocks)
    herit_result <- data.table::as.data.table(t(as.data.frame(result)))
    herit_result$Source <- "V(G)/Vp"
    h2_liab <- herit_result$h2 * bigsnpr::coef_to_liab(prevalance)
    herit_result2 <- herit_result
    herit_result2$h2 <- h2_liab
    herit_result2$Source <- "V(G)/Vp_L"
    herit_result3 <- rbind(herit_result, herit_result2)
    herit_result3 <- herit_result3[, c(5, 3, 4, 1, 2)]
    colnames(herit_result3) <- c("Source", "Variance", "SE", "Intercept", "Int_SE")
    return(herit_result3)
  } else {
    bimfile <- read.table(paste0(DataDir, "/", finput, ".bim"))
    chrnum <- 1:length(unique(bimfile$V1))

    chrwiseLD <- function(chrnum){

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
      GP <- GeneProtein(ResultDir = ResultDir, hg = hg,chromosome = chromosome)

      print(paste0("Processing chromosome ",chromosome))

      if (chrnum %in% unique(test.df_beta$chr)){
        snpld <- ComputeLD(DataDir = DataDir,ResultDir = ResultDir,finput = finput, ByCHR = TRUE, CHRnum = chromosome, r2_LD = r2_LD )

        result1 <- ComputeLDSC(snpld = snpld, test.df_beta = test.df_beta, ncores = ncores, LDSC_blocks = LDSC_blocks)
      }else{
        print("Summary statistics doesn't contain all the chromosomes.")
        result1 <- t(data.frame(NA,NA,NA,NA))
        colnames(result1) <- "result1"
      }
      result1 <- as.data.frame(result1)
      #print(result1)
      herit_result1 <- data.table::as.data.table(t(result1))
      herit_result <- data.table::as.data.table(cbind(chromosome,snp_proportion,GP,herit_result1))
      herit_result$Source <- "V(G)/Vp"
      #Scaling coefficient to convert e.g. heritability to the liability scale.
      h2_liab <- herit_result$h2 * bigsnpr::coef_to_liab(prevalance)# keeping the default value of K_gwas = 0.5, since n_eff is taken care of.
      herit_result2 <- herit_result
      herit_result2$h2 <- h2_liab
      herit_result2$Source <- "V(G)/Vp_L"
      #herit_result3 <- rbind(herit_result,herit_result2)
      return(herit_result2)

    }

    result <- data.table::rbindlist(lapply(chrnum, chrwiseLD), fill = TRUE)
    x <- na.omit(result)
    if (nrow(x)<3){
      print("Not enough data points for plots.")
      result1 <- result[,1:9]
      result2 <- na.omit(result1)
      colnames(result2) <- c("chromosome","snp_proportion","no.of.genes","no.of.proteins","Intercept","Int_SE","Variance","SE","Source")

      return(result2)

    }else{
      result1 <- result[result$Source == "V(G)/Vp",]##Check if needed
      chr_mb <- ChrLength(hg=hg)
      result2 <- merge(result1,chr_mb, by.y = "chrom", by.x = "chromosome")
      colnames(result2) <- c("chromosome","snp_proportion","no.of.genes","no.of.proteins","Intercept","Int_SE","Variance","SE","Source","size_mb")
      #create plot with regression line, regression equation, Pearson correlation and p-value.
      suppressWarnings(PlotHeritability(Hdata = result2,miMAF = miMAF, maMAF = maMAF))

      result3 <- result2
      colnames(result3) <- c("Chromosome","SNP_proportion","No.of.genes","No.of.proteins","Intercept","Int_SE","Variance","SE","Source","Size_mb")
      result3 <- result3[,c(1:4,10,9,7,8,5,6)]
      return(result3)
    }

  }
}

## Function 85
## Added in 3.0
processGREMLModel <- function(DataDir, ResultDir, finput, byCHR, autosome, Xsome, partGRM, nGRM, cripticut, minMAF, maxMAF, REMLalgo, nitr, cat_covarfile, quant_covarfile, prevalance, PlotIndepSNP, highLD_regions, IndepSNP_window_size, IndepSNP_step_size, IndepSNP_r2_threshold, hg, miMAF, maMAF, ncores) {

  if (!checkFiles(DataDir, finput)) {
    stop("Missing required Plink files in the specified DataDir.")
  }

  setupGCTA(ResultDir)

  famfile <- read.table(paste0(DataDir, "/", finput, ".fam"))[, c(1, 2, 6)]
  write.table(famfile, file = paste0(ResultDir, "/phenofile.phen"), sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)

  if (byCHR == FALSE){

    if (autosome == TRUE && Xsome == FALSE){

      ## Compute GRM
      ComputeGRMauto(DataDir = DataDir, ResultDir = ResultDir, finput = finput,
                     partGRM = partGRM, nGRM = nGRM, cripticut = cripticut,  minMAF = minMAF, maxMAF = maxMAF)

      ## Updated this part in V3.0

      herit_result <- ComputeREMLone(DataDir = DataDir, ResultDir = ResultDir, REMLalgo = REMLalgo, nitr = nitr, phenofile = "phenofile.phen", cat_covarfile = cat_covarfile,
                                     quant_covarfile = quant_covarfile, prevalance = prevalance, chr = "chromosome", grmfile = "test", ncores = ncores)

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
                                     quant_covarfile = quant_covarfile, prevalance = prevalance, chr = "chromosome", grmfile = "xtest", ncores = ncores)

      return(herit_result)

    }else{
      print("Autosome and Xsome cannot be set as FALSE together.")
    }

  }else{

    bimfile <- read.table(paste0(DataDir, "/", finput, ".bim"))

    chrnum <- 1:length(unique(bimfile$V1))
    #chrnum <- 1:3

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
      GP <- GeneProtein(ResultDir = ResultDir, hg = hg,chromosome = chromosome)

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

        herit_result <- ComputeREMLone(DataDir = DataDir, ResultDir = ResultDir, REMLalgo = REMLalgo, nitr = nitr, phenofile = "phenofile.phen", cat_covarfile = cat_covarfile,
                                       quant_covarfile = quant_covarfile, prevalance = prevalance, chr = chromosome, grmfile = "test", ncores = ncores)


        herit_result <- data.table::as.data.table(cbind(chromosome,snp_proportion,GP,herit_result))

        return(herit_result)
      }

    }

    result <- data.table::rbindlist(lapply(chrnum,chrwiseRELM),fill = TRUE)

    x <- na.omit(result)
    if (nrow(x)<3){
      print("Not enough data points for plots.")
      return(result)
    }else{
      result1 <- result[result$Source == "V(G)/Vp",]

      chr_mb <- ChrLength(hg=hg)
      result2 <- merge(result1,chr_mb, by.y = "chrom", by.x = "chromosome")

      #create plot with regression line, regression equation, Pearson correlation and p-value.
      suppressWarnings(PlotHeritability(Hdata = result2, miMAF = miMAF, maMAF = maMAF))

      result3 <- merge(result,chr_mb, by.y = "chrom", by.x = "chromosome")
      result3 <- result3[,c(1:4,8,5:7)]
      colnames(result3) <- c("Chromosome","SNP_proportion","No.of.genes","No.of.proteins","Size_mb","Source","Variance","SE")

      return(result3)
    }
  }
}

# Function 86
getCI = function(mn1, se1, method){
  remov = c(0, NA)
  mn    = mn1[! mn1 %in% remov]
  se    = se1[! mn1 %in% remov]
  vars  <- se^2
  vwts  <- 1/vars

  fixedsumm <- sum(vwts * mn)/sum(vwts)
  Q         <- sum(((mn - fixedsumm)^2)/vars)
  df        <- length(mn) - 1
  tau2      <- max(0, (Q - df)/(sum(vwts) - sum(vwts^2)/sum(vwts)) )

  if (method == "fixed"){ wt <- 1/vars } else { wt <- 1/(vars + tau2) }

  summ <- sum(wt * mn)/sum(wt)
  if (method == "fixed")
    varsum <- sum(wt * wt * vars)/(sum(wt)^2)
  else varsum <- sum(wt * wt * (vars + tau2))/(sum(wt)^2)

  summtest   <- summ/sqrt(varsum)
  df         <- length(vars) - 1
  se.summary <- sqrt(varsum)
  pval       = 1 - pchisq(summtest^2,1)
  pvalhet    = 1 - pchisq(Q, df)
  L95        = summ - 1.96*se.summary
  U95        = summ + 1.96*se.summary
  # out = c(round(c(summ,L95,U95),2), format(pval,scientific=TRUE), pvalhet)
  # c("OR","L95","U95","p","ph")
  # return(out)

  out = c(paste(round(summ,3), ' [', round(L95,3), ', ', round(U95,3), ']', sep=""),
          format(pval, scientific=TRUE), round(pvalhet,3))
  # c("OR","L95","U95","p","ph")
  return(out)
}

## Function 87
## Updated in 3.0
topForestplot <- function(i,MR2,Sbeta){

  SNPs <- MR2$SNP[i]
  Fixed_Effect <- MR2[i,c("BETA","CIfixedLL","CIfixedUL")]
  Random_Effect <- MR2[i,c("BETA.R.","CIrandomLL","CIrandomUL")]
  Weighted_Effect <- MR2[i,c("WEIGHTED_Z","CIweightedLL","CIweightedUL")]
  Study_EFFect <- Sbeta[Sbeta$SNP == SNPs,,drop = FALSE]
  Study_EFFect <- Study_EFFect[,-1]
  D1 <- rbind(Study_EFFect,Fixed_Effect,Random_Effect,Weighted_Effect,use.names = FALSE)
  #D1 <- rbind(Study_EFFect,Fixed_Effect,Random_Effect,use.names = FALSE)
  D1$V2 <- row.names(D1)
  D1$index <- as.integer(D1$V2)
  D1$study <- paste0("S",D1$index)
  D1 <- D1[,c(6,5,1:3)]
  colnames(D1) <- c("study","index","effect","lower","upper")
  D1[nrow(D1),"study"] <- "W"
  D1[nrow(D1)-1,"study"] <- "R"
  D1[nrow(D1)-2,"study"] <- "F"
  D1$CI <- paste0("(",D1$lower,",",D1$upper,")")
  df <- D1
  ############
  forest.plot <- function(x, intervals, labels = NULL, main = NULL, xlab = "Effect size",
                          pchs = rep(19,length(x)), cols = rep("black", length(x)),
                          cexs = rep(1,length(x))){
    K = length(x)
    stopifnot(nrow(intervals) == K)
    graphics::plot(0, col="white", xlim = c( min(c(intervals[,1],0) - 0.05), max(c(intervals[,2],0) + 0.05)),
                   ylim = c(0, K+1), xlab = xlab, ylab = "", yaxt = "n",main = main)
    graphics::axis(2, at = K:1, labels = labels, cex.axis = 0.8)
    graphics::arrows(intervals[,1], K:1, intervals[,2], K:1,
                     code = 3, angle = 90, length = 0.02, col = cols)
    graphics::points(x, K:1, pch = pchs, cex = cexs, col = cols)
    graphics::abline(v = 0,lty = 2)
  }

  suppressWarnings(forest.plot(D1$effect, intervals = as.matrix(D1[,4:5]), labels = D1$study, main = SNPs, xlab = "Effect size (beta, 95% CI)\nW:Weighted, R:Random, F:Fixed, S1, S2,..:Studies",
                               pchs = c(rep(19,length(D1$effect)-3),18,18,18), cexs = c(rep(.8,length(D1$effect)-3),1.3,1.3,1.3), cols = c(rep(1,length(D1$effect)-3),4,4,4)))

}

## Function 88
## Updated in 3.0
allForestplot <- function(i,MR2,Sbeta){

  SNPs <- MR2$SNP[i]
  Fixed_Effect <- MR2[i,c("BETA","CIfixedLL","CIfixedUL")]
  Random_Effect <- MR2[i,c("BETA.R.","CIrandomLL","CIrandomUL")]
  Weighted_Effect <- MR2[i,c("WEIGHTED_Z","CIweightedLL","CIweightedUL")]
  Study_EFFect <- Sbeta[Sbeta$SNP == SNPs,,drop = FALSE]
  Study_EFFect <- Study_EFFect[,-1]
  D1 <- rbind(Study_EFFect,Fixed_Effect,Random_Effect,Weighted_Effect,use.names = FALSE)
  #D1 <- rbind(Study_EFFect,Fixed_Effect,Random_Effect,use.names = FALSE)
  D1$V2 <- row.names(D1)
  D1$index <- as.integer(D1$V2)
  D1$study <- paste0("S",D1$index)
  D1 <- D1[,c(6,5,1:3)]
  colnames(D1) <- c("study","index","effect","lower","upper")
  D1[nrow(D1),"study"] <- "W"
  D1[nrow(D1)-1,"study"] <- "R"
  D1[nrow(D1)-2,"study"] <- "F"
  D1$CI <- paste0("(",round(D1$lower,2),",",round(D1$upper,2),")")
  df <- D1

  plot1 <- ggplot2::ggplot(data=df, ggplot2::aes(y=df$index, x=df$effect, xmin=df$lower, xmax=df$upper)) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbarh(height=.1) +
    ggplot2::scale_y_continuous(breaks=1:nrow(df), labels=df$study) +
    ggplot2::labs(title= paste0('Forest plot for ',SNPs), x='Effect Size (95% CI)', y = 'Studies and tests') +
    ggplot2::geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
    #theme_minimal()
    ggplot2::theme(plot.margin = grid::unit(c(.5, 0, .5, 1), "cm"))

  ## Create the table-base pallete
  table_base <- ggplot2::ggplot(dat = df, ggplot2::aes(y=df$study)) +
    ggplot2::ylab(NULL) + ggplot2::xlab("  ") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size=12),
                   axis.text.x = ggplot2::element_text(color="white", hjust = -3, size = 25), ## This is used to help with alignment
                   axis.line = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   legend.position = "none",
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   plot.margin = grid::unit(c(.6, 0, 0, .6), "cm"))

  ## OR point estimate table
  tab1 <- table_base +
    ggplot2::labs(title = "space") +
    ggplot2::geom_text(ggplot2::aes(y = rev(df$index), x = 1, label = sprintf("%0.1f", round(rev(df$effect), digits = 1))), size = 4) + ## decimal places
    ggplot2::ggtitle("BETA")

  ## 95% CI table
  tab2 <- table_base +
    ggplot2::geom_text(ggplot2::aes(y = rev(df$index), x = 1, label = rev(df$CI)), size = 4) +
    ggplot2::ggtitle("95% CI")
  ## Merge tables with plot
  lay <-  matrix(c(1,1,1,1,1,1,1,1,1,1,3,4,4), nrow = 1)

  gridExtra::grid.arrange(plot1, tab1, tab2, layout_matrix = lay)

}

## Function 89
## Added in 3.0
## Apply Genomic Control
getGCse <- function(SummData,ResultDir){
  print("Applying study-specific genomic control.")
  s1 <- read.table(paste0(ResultDir,"/",SummData),header = TRUE)

  # From p-values, calculate chi-squared statistic
  chisq <- qchisq(1-na.omit(s1$P),1)
  lamdaGC <- median(chisq)/qchisq(0.5,1)
  s1$SE <- s1$SE * sqrt(lamdaGC)
  write.table(s1,paste0(ResultDir,"/",SummData),col.names = TRUE, row.names = FALSE, quote = FALSE)
  return()
}

## Function 90
# Added in 3.0
# Compute confidence interval
getStudyCI <- function(SummData,MRsnps,ResultDir){
  s1 <- read.table(paste0(ResultDir,"/",SummData),header = TRUE)
  s11 <- merge(MRsnps,s1, by = "SNP")
  S1 <- s11[,c("SNP","BETA","L95","U95")]
  return(S1)
}

## Function 91
# Added in 3.0
adjustPvalThreshold <- function(top_snp_pval, MR, pval_filter) {
  min_pval <- switch(pval_filter,
                     "R" = min(MR$P, na.rm = TRUE),
                     "F" = min(MR$P.F., na.rm = TRUE),
                     "W" = min(MR$P.WZ., na.rm = TRUE))
  if (top_snp_pval < min_pval) {
    adjusted_pval <- min_pval * 100
    print(paste0("Minimum p-value is higher than provided threshold. Using ", adjusted_pval))
    return(adjusted_pval)
  } else {
    return(top_snp_pval)
  }
}

## Function 92
# Added in 3.0
filterSNPsForForestPlot <- function(MR, top_snp_pval, pval_filter) {
  pval_col <- switch(pval_filter,
                     "R" = "P",
                     "F" = "P.F.",
                     "W" = "P.WZ.")
  filtered_MR <- MR[MR[[pval_col]] <= top_snp_pval, , drop = FALSE]
  return(filtered_MR[order(filtered_MR[[pval_col]]), , drop = FALSE])
}

## Function 93
# Added in 3.0
# Meta-analysis Processing
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

## Function 94
# Added in 3.0
generatePlots <- function(MRfiltered, Sbeta, ResultDir, plotname, useSNPposition, pval_threshold_manplot, chosen_snps_file) {

  # Determine the number of SNPs to plot
  numSNPs <- min(length(unique(MRfiltered$SNP)), 10)
  if (length(unique(MRfiltered$SNP)) > 10) {
    print("Maximum 10 Forest plots of 10 chosen SNPs will be drawn in the plot window. For all other Forest plots, please check ResultDir.")
  }

  # Generate Forest plots
  invisible(suppressWarnings(lapply(1:numSNPs, topForestplot, MR2 = MRfiltered, Sbeta = Sbeta)))

  # Add mtext for all plots
  options(warn = -1)
  invisible(suppressWarnings(graphics::mtext(text = "Studies (black) and tests (blue)", side = 4, line = 0, outer = TRUE)))
  invisible(suppressWarnings(graphics::mtext(text = "Effect size (beta, 95% CI)", side = 1, line = 1, outer = TRUE)))

  if (is.null(chosen_snps_file)) {
    invisible(suppressWarnings(graphics::mtext(text = "Forest plots of a few chosen SNPs", side = 3, line = .5, outer = TRUE)))
  } else {
    invisible(suppressWarnings(graphics::mtext(text = "Forest plots of the top SNPs", side = 3, line = .5, outer = TRUE)))
  }


}

## Function 95
ComputeBivarREMLone <- function(DataDir, ResultDir, REMLalgo = c(0,1,2), nitr = 100, phenofile, cat_covarfile = NULL,
                                quant_covarfile = NULL, grmfile = "test", excludeResidual = c("FALSE","TRUE"), chr, ncores = bigparallelr::nb_cores()){

  if (excludeResidual=="FALSE"){
    ExResi = NULL
  }else{
    ExResi = "--reml-bivar-nocove"
  }


  if (is.null(cat_covarfile)){
    catcovar = NULL
    catcovarval = NULL
  }else{
    catcovar = "--ccovar"
    catcovarval = paste0(DataDir,"/",cat_covarfile)
  }

  if (is.null(quant_covarfile)){
    quantcovar = NULL
    quantcovarval = NULL
  }else{
    quantcovar = "--qcovar"
    quantcovarval = paste0(DataDir,"/",quant_covarfile)
  }

  invisible(sys::exec_wait(
    paste0(ResultDir,"/","./gcta-1.94.1"),
    args = c(
      "--reml-bivar",
      "--reml-alg",REMLalgo,
      "--reml-maxit",nitr,
      "--grm",
      paste0(ResultDir,"/",grmfile),
      "--pheno", paste0(ResultDir,"/",phenofile),
      catcovar, catcovarval,
      quantcovar, quantcovarval,
      ExResi,
      "--reml-no-constrain",
      "--thread-num",ncores,
      "--out",
      paste0(ResultDir,"/",chr,"test_bireml")
    ),
    std_out = FALSE,
    std_err = FALSE
  ))

  if (file.exists(paste0(ResultDir,"/",chr,"test_bireml.hsq"))){

    resultREML <- na.omit(data.table::fread(paste0(ResultDir,"/",chr,"test_bireml.hsq"),fill = TRUE))
    return(resultREML)
  }else{
    print(grep("Error", readLines(paste0(ResultDir,"/",chr,"test_bireml.log")), value = TRUE))
    print("An error occurs, please try byCHR = TRUE, check different options, SNP partitioning or data")
    x <- data.frame(NA,NA,NA)
    colnames(x) <- c("Source", "Variance","SE")
    return(x)
  }



}

## Function 96
ComputeBivarREMLmulti <- function(DataDir, ResultDir, REMLalgo = c(0,1,2), nitr = 100, phenofile, cat_covarfile = NULL,
                                  quant_covarfile = NULL, grmfile = "multi_GRMs.txt",excludeResidual = c("FALSE","TRUE"), ncores = bigparallelr::nb_cores()){

  if (excludeResidual=="FALSE"){
    ExResi = NULL
  }else{
    ExResi = "--reml-bivar-nocove"
  }
  if (is.null(cat_covarfile)){
    catcovar = NULL
    catcovarval = NULL
  }else{
    catcovar = "--ccovar"
    catcovarval = paste0(DataDir,"/",cat_covarfile)
  }

  if (is.null(quant_covarfile)){
    quantcovar = NULL
    quantcovarval = NULL
  }else{
    quantcovar = "--qcovar"
    quantcovarval = paste0(DataDir,"/",quant_covarfile)
  }

  ## Create multi_GRMs.txt in ResultDir
  fileConn<-file(paste0(ResultDir,"/multi_GRMs.txt"))
  writeLines(c(paste0(ResultDir,"/test"),paste0(ResultDir,"/xtest")), fileConn)
  close(fileConn)

  args <- c(
    "--reml-alg", REMLalgo,
    "--reml-bivar", 1, 2,
    "--reml-maxit", nitr,
    "--mgrm", paste0(ResultDir, "/multi_GRMs.txt"),
    "--pheno", paste0(ResultDir, "/", phenofile),
    catcovar, catcovarval,
    quantcovar, quantcovarval,
    ExResi,
    "--reml-no-constrain",
    "--thread-num", ncores,
    "--out", paste0(ResultDir, "/test_bireml")
  )

  # Executing GCTA with the specified arguments
  executeGCTA(ResultDir, args)

  if (file.exists(paste0(ResultDir,"/test_bireml.hsq"))){

    resultREML <- na.omit(data.table::fread(paste0(ResultDir,"/test_bireml.hsq"),fill = TRUE))

    return(resultREML)

  }else{

    x1 <- read.table(text = grep(paste0(nitr-1,"\t"), readLines(paste0(ResultDir,"/test_bireml.log")), value = TRUE))
    colnames(x1)<- NULL
    rownames(x1)<- NULL
    #print(x1)
    print(grep("Error", readLines(paste0(ResultDir,"/test_bireml.log")), value = TRUE))
    print(grep("Note: to constrain", readLines(paste0(ResultDir,"/test_bireml.log")), value = TRUE))
    print("Convergence problem occurs, please try byCHR = TRUE, check different options, SNP partitioning or ensure the quality of the input data.")
    print("The result will be provided for the last iteration.")
    x <- as.data.frame(t(x1))
    x <- x[-1,]
    colnames(x) <- c("Source", "Variance")
    return(x)
  }

}

## Function 97
# sumFREGAT (2017-2022) Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS
geneTestScoreFile <- function(ResultDir,data, reference = 'ref1KG.MAC5.EUR_AF.RData', output.file.prefix) {

  OS <- Sys.info()['sysname']
  if (OS == "Windows"){
    print("Currently this function maynot work on Windows as bgzip and tabix for windows are down. Pleas use linux environment.")

    utils::download.file(destfile = paste0(ResultDir,"/","bgzip_tabix.zip"),
                         "https://github.com/boseb/bose_binaries/raw/main/bgzip_tabix.zip", mode = "wb", quiet = TRUE,
    )
  }else{
    utils::download.file(destfile = paste0(ResultDir,"/","bgzip_tabix.zip"),
                         "https://github.com/boseb/bose_binaries/raw/main/bgzip_tabix.zip", quiet = TRUE,
    )
  }
  utils::unzip(paste0(ResultDir,"/","bgzip_tabix.zip"), exdir = ResultDir)

  Sys.chmod(paste0(ResultDir,"/bgzip"), mode = "0777", use_umask = TRUE)
  Sys.chmod(paste0(ResultDir,"/tabix"), mode = "0777", use_umask = TRUE)
  # 'CHROM', 'POS', 'ID', 'EA', 'P', 'BETA', 'EAF'
  if (length(data) == 1) {
    input.file <- data
    if (requireNamespace("data.table", quietly = TRUE)) {
      suppressWarnings(df <- data.table::fread(input.file, header = TRUE, data.table = FALSE))
    } else {
      df <- read.table(input.file, header = TRUE, as.is = TRUE)
    }
  } else if (length(data) > 1) {
    df <- data
    input.file <- 'scores'
  }

  cn <- toupper(colnames(df))
  v <- which(cn %in% c('CHR', 'CHROMOSOME', 'CHROM'))
  if (length(v) == 1) colnames(df)[v] <- 'CHROM'
  v <- which(cn %in% c('POSITION', 'POSITIONS', 'MAP', 'POS'))
  if (length(v) == 1) colnames(df)[v] <- 'POS'
  v <- which(cn %in% c('PVALUE', 'PV', 'PVAL', 'P.VALUE', 'P_VALUE', 'P'))
  if (length(v) == 1) colnames(df)[v] <- 'P'
  v <- which(cn %in% c('RSID', 'RS.ID', 'RS_ID', 'SNP.ID', 'SNP_ID', 'ID'))
  if (length(v) == 1) colnames(df)[v] <- 'ID'
  v <- which(cn == 'EA')
  if (length(v) == 1) {
    colnames(df)[v] <- 'EFFECT.ALLELE'
    df[, 'EFFECT.ALLELE'] <- toupper(df[, 'EFFECT.ALLELE'])
  }

  # ID and PVAL mandatory
  # others from user file or reference

  ColNames <- c('ID', 'P')
  v <- !ColNames %in% colnames(df)
  if (sum(v)) stop(paste("Mandatory column(s) missing:", paste(ColNames[v], collapse = ', ')))

  df <- df[!is.na(df$P) & !is.na(df$ID), ]
  if (dim(df)[1] == 0) stop("No values assigned for P or ID")

  ColNames <- c('CHROM', 'POS', 'EAF')
  v <- !ColNames %in% colnames(df)
  take <- ColNames[v]
  if (sum(v)) print(paste("Columns that are missing and will be looked for in reference data:", paste(take, collapse = ', ')))
  take[take == 'EAF'] <- 'AF'

  if ('BETA' %in% colnames(df)) {
    df$BETA[df$BETA == 0] <- 1e-16
    if ('EFFECT.ALLELE' %in% colnames(df)) {
      colnames(df)[which(colnames(df) == 'REF')] <- 'REF0'
      colnames(df)[which(colnames(df) == 'ALT')] <- 'ALT0'
      take <- c(take, 'REF', 'ALT')
    } else {
      print("Effect allele column not found, effect sizes cannot be linked")
    }
  } else {
    print("Effect sizes (beta) column not found")
  }
  if (length(take) > 0) {
    is.ref <- 0
    is.ref.object <- 0
    if (length(reference) == 1) {
      if (!is.na(reference)) {
        if (file.exists(reference)) {
          is.ref <- 1
        } else {
          if (reference != '') print ("Reference file not found! Please download it from https://mga.bionet.nsc.ru/sumFREGAT/ref1KG.MAC5.EUR_AF.RData to use 1000 Genome Reference correlation matrices")
        }
      }
    } else if (length(reference) > 1) is.ref <- is.ref.object <- 1

    if (is.ref) {
      if (is.ref.object) {
        ref <- reference
      } else {
        print('Loading reference file...')
        ref <- get(load(reference))
      }
      colnames(ref) <- toupper(colnames(ref))
      if ('CHROM' %in% take & !'CHROM' %in% colnames(ref)) stop ("No CHROM column in data and reference")
      if ('POS' %in% take & !'POS' %in% colnames(ref)) stop ("No POS column in data and reference")
      v <- match(df$ID, ref$ID)

      if (!sum(v, na.rm = TRUE)) {
        if (all(c('CHROM', 'POS') %in% colnames(df))) {
          df$ind <- paste(df$CHROM, df$POS, sep = ':')
          print('No IDs matching, trying to link through map data...')
          ref$ind <- paste(ref$CHROM, ref$POS, sep = ':')
          v <- match(df$ind, ref$ind)
          if (sum(!is.na(v)) < (length(v) / 2)) {
            print("Too few variants match between input file and reference data")
            v <- NA
          }
        }
      }
      if (sum(v, na.rm = TRUE)) {
        print(paste(sum(!is.na(v)), "of", length(v), "variants found in reference"))
        vv <- take %in% colnames(ref)
        if (sum(!vv)) {
          print(paste("Columns that are missing in reference data:", paste(take[!vv], collapse = ', ')))
          if ('REF' %in% take & !'REF' %in% colnames(ref)) {
            print ("Reference alleles not found, effect sizes cannot be linked")
            df$BETA <- df$EFFECT.ALLELE <- NULL
          }
          if ('AF' %in% take & !'AF' %in% colnames(ref)) print ("Allele frequencies not found, some weighted tests will be unavailable")
        }
        df <- cbind(df, ref[v, take[vv]])
      }
    } else {
      v <- NA
    }
    if (sum(v, na.rm = TRUE) == 0) { # fail to open or link reference data
      if (any(c('CHROM', 'POS') %in% take)) stop ("Cannot find map data (chromosome, position)")
      if ('BETA' %in% colnames(df)) {
        warning ("Reference unavailable, effect sizes not linked")
        df$BETA <- df$EFFECT.ALLELE <- NULL
      }
    }
  }

  if ('REF' %in% colnames(df) & 'EFFECT.ALLELE' %in% colnames(df)) {
    v <- c()
    if (all(c('REF', 'REF0', 'ALT', 'ALT0') %in% colnames(df))) {
      v <- which((df$REF0 != df$REF & df$REF0 != df$ALT) | (df$ALT0 != df$REF & df$ALT0 != df$ALT))
    }
    if ('ALT' %in% colnames(df)) {
      v <- unique(c(v, which(df$EFFECT.ALLELE != df$REF & df$EFFECT.ALLELE != df$ALT)))
    }
    if (sum(v, na.rm = T)) {
      print(paste("Effect alleles or REF/ALT alleles do not match reference data for", sum(v), "variant(s)"))
      df[v, 'BETA'] <- NA
    }
    df[is.na(df$EFFECT.ALLELE) | is.na(df$REF), 'BETA'] <- NA
    v <- which(df$EFFECT.ALLELE == df$REF)
    #here we go
    df$BETA[v] <- -df$BETA[v]
    if ('EAF' %in% colnames(df)) {
      df$EAF[v] <- 1 - df$EAF[v]
      colnames(df)[colnames(df) == 'EAF'] <- 'AF'
    }
    print(paste('Effect sizes recoded for', length(v), 'variant(s)'))
  }

  if (any(df$P == 0)) {
    print("Some P values equal zero, will be assigned to minimum value in the sample")
    df$P[df$P == 0] <- min(df$P[df$P > 0])
  }
  df$Z <- qnorm(df$P / 2, lower.tail = FALSE)
  if ('BETA' %in% colnames(df)) {
    df$Z <- df$Z * sign(df$BETA)
    df$SE.BETA <- df$BETA / df$Z
  }

  if (!missing(output.file.prefix)) {
    fn <- paste(output.file.prefix, 'vcf', sep = '.')
  } else {
    fn <- paste(input.file, 'vcf', sep = '.')
  }

  df <- df[order(df[, 'POS']), ]
  df <- df[order(df[, 'CHROM']), ]
  if (!'ALT' %in% colnames(df)) df$ALT <- NA
  if (!'REF' %in% colnames(df)) df$REF <- NA
  vcf <- df[, c('CHROM', 'POS', 'ID', 'REF', 'ALT')]
  colnames(vcf)[1] <- '#CHROM'
  vcf$POS <- format(vcf$POS, scientific = FALSE)
  vcf$POS <- gsub(' ', '', vcf$POS)
  vcf <- cbind(vcf, QUAL = '.', FILTER = '.')
  vcf$INFO <- paste0('Z=', df$Z)
  title <- c('##INFO=<ID=Z,Number=1,Type=Float,Description="Z statistics">')

  if ('BETA' %in% colnames(df)) {
    vcf$INFO <- paste0(vcf$INFO, ';SE.Beta=', df$SE.BETA)
    title <- c(title, '##INFO=<ID=SE.Beta,Number=1,Type=Float,Description="SE Beta">')
  }

  if ('EAF' %in% colnames(df)) colnames(df)[colnames(df) == 'EAF'] <- 'AF'
  if ('AF' %in% colnames(df)) {
    vcf$INFO <- paste0(vcf$INFO, ';AF=', df$AF)
    title <- c(title, '##INFO=<ID=AF,Number=1,Type=Float,Description="Frequency of alternative allele">')
    print(paste0('Allele frequencies found and linked'))
  }

  a <- grep('\\bW', colnames(df))
  if (length(a) == 1) {
    vcf$INFO <- paste0(vcf$INFO, ';W=', df[, a])
    title <- c(title, '##INFO=<ID=W,Number=1,Type=Float,Description="Weights">')
    print(paste0("User weights ('", colnames(df)[a], "') found and linked"))
  }

  a <- grep('\\bANNO', colnames(df), value = TRUE)
  if (length(a) == 1) {
    vcf$INFO <- paste0(vcf$INFO, ';ANNO=', df[, a])
    title <- c(title, '##INFO=<ID=ANNO,Number=1,Type=String,Description="Variants annotations">')
    print(paste0("Annotations ('", colnames(df)[a], "') found and linked"))
  }

  a <- grep('\\bPROB', colnames(df), value = TRUE)
  for (an in a) {
    vcf$INFO <- paste0(vcf$INFO, ';', an, '=', df[, as.character(an)])
    title <- c(title, paste0("##INFO=<ID=", an, ",Number=1,Type=Float,Description='", an, "'>"))
    print(paste0("Column '", an, "' linked"))
  }

  #write.table(title, fn, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
  write.table(title, paste0(ResultDir,"/",fn), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')

  # if (requireNamespace("data.table", quietly = TRUE)) {
  #   suppressWarnings(data.table::fwrite(vcf, fn, row.names = FALSE, quote = FALSE, append = TRUE, col.names = TRUE, sep = '\t', na = 'NA'))
  # } else {
  #   suppressWarnings(write.table(vcf, fn, row.names = FALSE, quote = FALSE, append = TRUE, sep = '\t'))
  # }

  if (requireNamespace("data.table", quietly = TRUE)) {
    suppressWarnings(data.table::fwrite(vcf, paste0(ResultDir,"/",fn), row.names = FALSE, quote = FALSE, append = TRUE, col.names = TRUE, sep = '\t', na = 'NA'))
  } else {
    suppressWarnings(write.table(vcf, paste0(ResultDir,"/",fn), row.names = FALSE, quote = FALSE, append = TRUE, sep = '\t'))
  }

  fn.gz <- paste(fn, 'gz', sep = '.')
  # if (file.exists(fn.gz)) system(paste('rm', fn.gz))
  # system(paste('./bgzip', fn))
  # system(paste('./tabix -p vcf', fn.gz))
  # print(paste('File', fn.gz, 'has been created'))
  if (file.exists(paste0(ResultDir,"/",fn.gz))) system(paste('rm', paste0(ResultDir,"/",fn.gz)))
  system(paste0(ResultDir,"/",'./bgzip ', ResultDir,"/",fn))
  system(paste0(ResultDir,"/",'./tabix -p vcf ', ResultDir,"/",fn.gz))
  print(paste('File', fn.gz, 'has been created'))

}

## Function 98
# Download_reference(refdata = "HapMapIII_NCBI36", wdir = DataDir)
# Download_reference <- function(refdata, wdir){
#   OS <- Sys.info()['sysname']
#   options(timeout=200)
#   if (refdata == "HapMapIII_NCBI36"){
#     if (OS == "Linux"){
#
#     utils::download.file(destfile = paste0(wdir,"/",refdata,".zip"),
#                        "https://figshare.com/ndownloader/files/40585145", quiet = TRUE,
#     )
#     }else if (OS == "Windows"){
#
#       utils::download.file(destfile = paste0(wdir,"/",refdata,".zip"),
#                            "https://figshare.com/ndownloader/files/40585145", quiet = TRUE,mode = "wb",
#       )
#     }else if (OS == "macOS"){
#       utils::download.file(destfile = paste0(wdir,"/",refdata,".zip"),
#                            "https://figshare.com/ndownloader/files/40585145", quiet = TRUE,
#       )
#     }
#     utils::unzip(paste0(wdir,"/",refdata,".zip"), exdir = wdir, junkpaths = TRUE)
#     utils::unzip(paste0(wdir,"/",refdata,".zip"), exdir = wdir)
#     invisible(file.remove(paste0(wdir,"/",refdata,".zip")))
#
#     }else if(refdata == "ThousandGenome"){
#       options(timeout=200)##Takes more time, that's why this line
#       if (OS == "Linux"){
#         utils::download.file(destfile = paste0(wdir,"/",refdata,".tar.gz"),
#                              "https://figshare.com/ndownloader/files/40728539", quiet = TRUE,
#         )
#       }else if (OS == "Windows"){
#         utils::download.file(destfile = paste0(wdir,"/",refdata,".tar.gz"),
#                              "https://figshare.com/ndownloader/files/40728539", quiet = TRUE,mode = "wb",
#         )
#       }else if (OS == "macOS"){
#         utils::download.file(destfile = paste0(wdir,"/",refdata,".tar.gz"),
#                              "https://figshare.com/ndownloader/files/40728539", quiet = TRUE,
#         )
#       }
#
#       untar(paste0(wdir,"/",refdata,".tar.gz"),exdir = wdir)
#       invisible(file.remove(paste0(wdir,"/",refdata,".tar.gz")))
#     }
#
#   print(paste0("Reference:", refdata, " downloaded."))
#
# }

## Function 99
## Added in 3.0
validateGXwasInputs <- function(DataDir, ResultDir, finput, trait, standard_beta, xmodel, sex, xsex, covarfile, interaction, covartest, Inphenocov, combtest, MF.zero.sub, B, MF.mc.cores, MF.na.rm, MF.p.corr, plot.jpeg, plotname, snp_pval, annotateTopSnp, suggestiveline, genomewideline, ncores) {
  # Validate existence of required Plink files
  if (!checkFiles(DataDir, finput)) return("Missing required Plink files in the specified DataDir.")

  # Validate trait, standard_beta, sex, xsex, interaction, plot.jpeg, annotateTopSnp
  if (!trait %in% c("binary", "quantitative")) return("Invalid trait. Choose 'binary' or 'quantitative'.")
  if (!is.logical(standard_beta)) return("standard_beta must be TRUE or FALSE.")
  if (!is.logical(sex) || !is.logical(xsex)) return("Sex and xsex must be TRUE or FALSE.")
  if (!is.logical(interaction)) return("interaction must be TRUE or FALSE.")
  if (!is.logical(plot.jpeg)) return("plot.jpeg must be TRUE or FALSE.")
  if (!is.logical(annotateTopSnp)) return("annotateTopSnp must be TRUE or FALSE.")

  # Validate xmodel
  if (!xmodel %in% c("FMcomb01", "FMcomb02", "FMstatrified", "GWAScxci")) return("Invalid xmodel value.")

  # Validate covarfile
  if (!is.null(covarfile) && !file.exists(file.path(DataDir, covarfile))) return("Specified covarfile does not exist.")

  # Validate covartest and Inphenocov against covarfile
  validateCovarParams <- function(param, covarfile) {
    if (param == "ALL"|| is.null(param)) return(TRUE)
    if (is.numeric(param)) {
      if (is.null(covarfile)) return("covarfile needed for numeric covartest or Inphenocov.")
      covarData <- tryCatch(read.table(file.path(DataDir, covarfile), header = TRUE), error = function(e) NULL)
      if (is.null(covarData)) return("Failed to read covarfile.")
      maxIndex <- ncol(covarData) - 2
      if (any(param < 1) || any(param > maxIndex)) return("Indices out of bounds in covartest or Inphenocov.")
      return(TRUE)
    }
    return(FALSE)
  }
  if (!validateCovarParams(covartest, covarfile)) return("Invalid covartest.")
  if (!validateCovarParams(Inphenocov, covarfile)) return("Invalid Inphenocov.")

  # Validate combtest
  if (!all(combtest %in% c("fisher.method", "fisher.method.perm", "stouffer.method"))) return("Invalid combtest methods.")

  # Validate numeric parameters
  if (!is.numeric(MF.zero.sub) || MF.zero.sub < 0) return("Invalid value for MF.zero.sub.")
  if (!is.numeric(B) || B <= 0) return("Invalid value for B.")
  if (!is.null(MF.mc.cores) && (!is.numeric(MF.mc.cores) || MF.mc.cores < 1 || floor(MF.mc.cores) != MF.mc.cores)) {
    return("Invalid value for MF.mc.cores.; must be >= 1")
  }
  if (!is.numeric(snp_pval) || snp_pval <= 0 || snp_pval > 1) return("Invalid value for snp_pval.")
  if (!is.numeric(suggestiveline) || suggestiveline <= 0) return("Invalid value for suggestiveline.")
  if (!is.numeric(genomewideline) || genomewideline <= 0) return("Invalid value for genomewideline.")
  if (!is.null(ncores) && (!is.numeric(ncores) || ncores < 0 || floor(ncores) != ncores)) {
    return("Invalid value for MF.mc.cores.")
  }

  # Return NULL if all validations pass
  return(NULL)
}

## Function 100
## Added in 3.0
# Input Validation Helper Function
validateAncestryCheckInputs <- function(DataDir, ResultDir, finput, reference, filterSNP, studyLD, studyLD_window_size, studyLD_step_size, studyLD_r2_threshold, referLD, referLD_window_size, referLD_step_size, referLD_r2_threshold, highLD_regions, study_pop, outlier, outlierOf, outlier_threshold) {

  # Validate directory paths
  if (!dir.exists(DataDir)) {
    stop("DataDir does not exist: ", DataDir)
  }
  if (!dir.exists(ResultDir)) {
    stop("ResultDir does not exist: ", ResultDir)
  }

  # Validate file input
  if (!file.exists(paste0(DataDir, "/", finput, ".bim"))) {
    stop("BIM file not found for input: ", finput)
  }

  # Validate reference
  if (!is.character(reference) || !(reference %in% c("HapMapIII_NCBI36", "ThousandGenome"))) {
    stop("Reference must be either 'HapMapIII_NCBI36' or 'ThousandGenome'")
  }

  # Validate boolean flags
  if (!is.logical(filterSNP) || !is.logical(studyLD) || !is.logical(referLD) || !is.logical(outlier)) {
    stop("filterSNP, studyLD, referLD, and outlier must be boolean values")
  }

  # Validate numeric inputs
  if (!is.numeric(studyLD_window_size) || !is.numeric(studyLD_step_size) || !is.numeric(studyLD_r2_threshold) ||
      !is.numeric(referLD_window_size) || !is.numeric(referLD_step_size) || !is.numeric(referLD_r2_threshold) ||
      !is.numeric(outlier_threshold)) {
    stop("Window sizes, step sizes, r2 thresholds, and outlier_threshold must be numeric")
  }

  if (!is.null(highLD_regions)) {
    if (!is.data.frame(highLD_regions) || ncol(highLD_regions) != 4) {
      stop("highLD_regions must be a dataframe with exactly four columns.")
    }
  }

  # Read IID from .fam file
  famFilePath <- paste0(DataDir, "/", finput, ".fam")
  if (!file.exists(famFilePath)) {
    stop(".fam file not found: ", famFilePath)
  }
  famData <- read.table(famFilePath, stringsAsFactors = FALSE)
  famIIDs <- famData[, 2]

  # Validate study_pop
  if (!is.null(study_pop)) {
    if (!is.data.frame(study_pop) || ncol(study_pop) != 2) {
      stop("study_pop must be a dataframe with exactly two columns.")
    }
    # if (!all(study_pop[, 1] %in% famIIDs)) {
    #   stop("All IIDs in study_pop must exist in the .fam file.")
    # }
    if (!all(famIIDs %in% study_pop[, 1])) {
      stop("All IIDs in input file  must exist in the study_pop.")
    }
  }

  # Validate outlierOf
  if (outlier && !is.character(outlierOf)) {
    stop("outlierOf must be a character value when outlier is TRUE")
  }

  # Return TRUE if all validations pass
  return(TRUE)
}

## Function 101
## Added in 3.0
processLDstudyData <- function(ResultDir, highLD_regions, studyLD, studyLD_window_size, studyLD_step_size, studyLD_r2_threshold) {
  if (studyLD) {
    # Using the provided helper function to execute PLINK with specific arguments
    executePlinkAd(ResultDir, c(
      "--bfile", paste0(ResultDir, "/filtered_study_temp1"),
      "--exclude", "range", highLD_regions,
      "--indep-pairwise", studyLD_window_size, studyLD_step_size, studyLD_r2_threshold,
      "--allow-no-sex",                                                     ## Adding in 4.0
      "--make-bed", "--out", paste0(ResultDir, "/filtered_study_temp2"),
      "--silent"
    ))
  } else {
    executePlinkAd(ResultDir, c(
      "--bfile", paste0(ResultDir, "/filtered_study_temp1"),
      "--make-bed", "--out", paste0(ResultDir, "/filtered_study_temp2"),
      "--silent"
    ))
  }

  # Return a message based on the LD pruning status
  if (studyLD) {
    "LD pruning was done for study dataset."
  } else {
    "LD pruning for study dataset is recommended. Set studyLD == TRUE."
  }
}

## Function 102
## Added in 3.0
processLDreferenceData <- function(ResultDir, highLD_regions, referLD, referLD_window_size, referLD_step_size, referLD_r2_threshold) {
  if (referLD) {
    executePlinkAd(ResultDir, c(
      "--bfile", paste0(ResultDir, "/filtered_ref_temp1"),
      "--exclude", "range", highLD_regions,
      "--indep-pairwise", referLD_window_size, referLD_step_size, referLD_r2_threshold,
      "--allow-no-sex",                                                     ## Adding in 4.0
      "--make-bed", "--out", paste0(ResultDir, "/filtered_ref_temp2"),
      "--silent"
    ))
  } else {
    executePlinkAd(ResultDir, c(
      "--bfile", paste0(ResultDir, "/filtered_ref_temp1"),
      "--make-bed", "--out", paste0(ResultDir, "/filtered_ref_temp2"),
      "--silent"
    ))
  }

  if (referLD) {
    "LD pruning was done for reference dataset."
  } else {
    "LD pruning for reference dataset is recommended. Set referLD == TRUE."
  }
}

## Function 103
## Added in 3.0
filterATGCSNPs <- function(DataDir, ResultDir, finput, reference) {
  filterSNPs <- function(DataDir, file, SNPPrefix) {
    bimData <- read.table(file = paste0(DataDir, "/", file, ".bim"), stringsAsFactors = FALSE)
    SNP_AT <- bimData[which(bimData[, 5] == "A" & bimData[, 6] == "T"), 2, drop = FALSE]
    SNP_TA <- bimData[which(bimData[, 5] == "T" & bimData[, 6] == "A"), 2, drop = FALSE]
    SNP_GC <- bimData[which(bimData[, 5] == "G" & bimData[, 6] == "C"), 2, drop = FALSE]
    SNP_CG <- bimData[which(bimData[, 5] == "C" & bimData[, 6] == "G"), 2, drop = FALSE]

    SNP <- rbind(SNP_AT, SNP_TA, SNP_GC, SNP_CG)
    write.table(SNP, file = paste0(ResultDir, "/", SNPPrefix, "_SNP"), quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\r\n", sep = " ")

    return(SNP)
  }

  study_SNP <- filterSNPs(DataDir, finput, "study")
  ref_SNP <- filterSNPs(ResultDir, reference, "ref")

  # Execute Plink Commands
  executePlinkAd(ResultDir, c("--bfile", paste0(DataDir, "/", finput), "--exclude", paste0(ResultDir, "/study_SNP"), "--allow-no-sex", "--make-bed", "--out", paste0(ResultDir, "/filtered_study_temp1"), "--silent")) ## Added "--allow-no-sex" in 4.0
  executePlinkAd(ResultDir, c("--bfile", paste0(ResultDir, "/", reference), "--exclude", paste0(ResultDir, "/ref_SNP"), "--allow-no-sex", "--make-bed", "--out", paste0(ResultDir, "/filtered_ref_temp1"), "--silent")) ## Added "--allow-no-sex" in 4.0

  # Print Messages
  if (nrow(study_SNP) == 0) {
    print("No SNP had 'A-T' and 'G-C' in study data.")
  } else {
    print(paste0(nrow(study_SNP), " SNPs had 'A-T' and 'G-C' in study data. These SNPs were removed."))
  }

  if (nrow(ref_SNP) == 0) {
    print("No SNP had 'A-T' and 'G-C' in reference data.")
  } else {
    print(paste0(nrow(ref_SNP), " SNPs were 'A-T' and 'G-C' in reference data. These SNPs were removed."))
  }
}

## Function 104
## Added in 3.0
executePlinkForUnfilteredData <- function(DataDir, ResultDir, finput, reference) {
  executePlinkAd(ResultDir, c("--bfile", paste0(DataDir, "/", finput), "--make-bed", "--out", paste0(ResultDir, "/filtered_study_temp2"), "--silent"))
  executePlinkAd(ResultDir, c("--bfile", paste0(ResultDir, "/", reference), "--make-bed", "--out", paste0(ResultDir, "/filtered_ref_temp2"), "--silent"))

  print("A-T and C-G SNPs recommended to remove from both the reference and study data set by setting filterSNP == TRUE")
}

## Function 105
## Added in 3.0
findCommonSNPs <- function(ResultDir) {
  pruned_study <- read.table(file = paste0(ResultDir, "/","filtered_study_temp2", ".bim"), stringsAsFactors = FALSE)
  pruned_ref <- read.table(file = paste0(ResultDir, "/","filtered_ref_temp2", ".bim"), stringsAsFactors = FALSE)
  common_snps <- intersect(pruned_study$V2, pruned_ref$V2)

  if (length(common_snps) == 0){
    stop("No common SNPs found between study and reference data. This analysis cannot be done.")
  }

  write.table(common_snps, file = paste0(ResultDir, "/","common_snps"), quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\r\n", sep = " ")

  return(list(common_snps = common_snps, pruned_study = pruned_study, pruned_ref = pruned_ref))
}

## Function 106
## Added in 3.0
processCommonSNPs <- function(ResultDir, common_snps) {
  # Write common SNPs to a file
  write.table(
    common_snps,
    file = paste0(ResultDir, "/","common_snps"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    eol = "\r\n",
    sep = " "
  )

  # Extract common SNPs from study and reference datasets
  executePlinkAd(ResultDir, c("--bfile", paste0(ResultDir, "/","filtered_study_temp2"), "--extract", paste0(ResultDir, "/","common_snps"), "--allow-no-sex", "--make-bed", "--out", paste0(ResultDir, "/","filtered_study_temp3"),"--silent")) ## Added "--allow-no-sex" in 4.0
  executePlinkAd(ResultDir, c("--bfile", paste0(ResultDir, "/","filtered_ref_temp2"), "--extract", paste0(ResultDir, "/","common_snps"), "--allow-no-sex", "--make-bed", "--out", paste0(ResultDir, "/","filtered_ref_temp3"),"--silent"))

  # Remove temporary files
  removeTempFiles(ResultDir, "filtered_study_temp2")
  removeTempFiles(ResultDir, "filtered_ref_temp2")
}

## Function 107
## Added in 3.0
handleSnpAlleleFlips <- function(ResultDir, snp_allele_flips) {
  snp_allele_flips_count <- length(unique(snp_allele_flips))

  # Write snp_allele_flips to file if any exist
  if (snp_allele_flips_count != 0) {
    write.table(
      snp_allele_flips,
      file = paste0(ResultDir, "/snp_allele_flips"),
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      eol = "\r\n",
      sep = " "
    )


    executePlinkAd(ResultDir, c("--bfile", paste0(ResultDir, "/filtered_ref_temp4"),
                                "--flip", paste0(ResultDir, "/snp_allele_flips"),"--allow-no-sex",          ## Added in 4.0
                                "--make-bed", "--out", paste0(ResultDir, "/filtered_ref_temp5"),"--silent"))



  } else {
    # Execute PLINK without allele flips
    executePlinkAd(
      ResultDir,
      c("--bfile", paste0(ResultDir, "/filtered_ref_temp4"),
        "--make-bed", "--out", paste0(ResultDir, "/filtered_ref_temp5"),
        "--silent")
    )
  }

  # Remove temporary files using the provided helper function
  removeTempFiles(ResultDir, "filtered_ref_temp4")
}

## Function 108
## Added in 3.0
# Helper function for correcting chromosome mismatches
correctChromosomeMismatches <- function(ResultDir, common_snps, pruned_study, pruned_ref) {
  S1 <- pruned_study[match(common_snps, pruned_study[, 2]), , drop = FALSE]
  S2 <- pruned_ref[match(common_snps, pruned_ref[, "V2"]), , drop = FALSE]

  whChrNotSame <- which(S1[, 1] != S2[, "V1"])
  whPosNotSame <- which(S1[, 4] != S2[, "V4"])
  whChrPosNotSame <- union(whChrNotSame, whPosNotSame)

  snpSameNameDiffPos <- S1[whChrPosNotSame, c(1, 4, 2), drop = FALSE]

  if (nrow(snpSameNameDiffPos) != 0) {
    write.table(snpSameNameDiffPos, file = paste0(ResultDir, "/","snpSameNameDiffPos"), quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\r\n", sep = " ")
    executePlinkAd(ResultDir, c("--bfile", paste0(ResultDir, "/","filtered_ref_temp3"), "--update-chr", paste0(ResultDir, "/","snpSameNameDiffPos"), 1, 3, "--update-map", paste0(ResultDir, "/","snpSameNameDiffPos"), 2, 3, "--allow-no-sex", "--make-bed", "--out", paste0(ResultDir, "/","filtered_ref_temp4"),"--silent")) ## Adding "--allow-no-sex" in 4.0.
  } else {
    executePlinkAd(ResultDir, c("--bfile", paste0(ResultDir, "/","filtered_ref_temp3"), "--make-bed", "--out", paste0(ResultDir, "/","filtered_ref_temp4"), "--silent"))
  }

  removeTempFiles(ResultDir, "filtered_ref_temp3")

  return(list(S1 = S1, S2 = S2, snpSameNameDiffPos = snpSameNameDiffPos))
}

## Function 109
## Added in 3.0
handleAlleleFlipsWrong <- function(ResultDir, allele_flips_wrong) {
  allele_flips_wrong_count <- length(unique(allele_flips_wrong))

  # Write allele_flips_wrong to file if any exist
  if (allele_flips_wrong_count != 0) {
    write.table(
      allele_flips_wrong,
      file = paste0(ResultDir, "/allele_flips_wrong"),
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      eol = "\r\n",
      sep = " "
    )

    # Execute PLINK for allele flips wrong
    executePlinkAd(
      ResultDir,
      c("--bfile", paste0(ResultDir, "/filtered_ref_temp5"),
        "--exclude", paste0(ResultDir, "/allele_flips_wrong"),
        "--allow-no-sex",                                       ## Adding in 4.0.
        "--make-bed", "--out", paste0(ResultDir, "/filtered_ref_temp6"), "--silent")
    )


  } else {
    # Execute PLINK without allele flips wrong
    executePlinkAd(
      ResultDir,
      c("--bfile", paste0(ResultDir, "/filtered_ref_temp4"),
        "--make-bed", "--out", paste0(ResultDir, "/filtered_ref_temp5"),
        "--silent")
    )


  }

  # Remove temporary files using the provided helper function
  removeTempFiles(ResultDir, "filtered_ref_temp5")
}

## Function 110
## Added in 3.0
mergeDatasetsAndPerformPCA <- function(ResultDir) {
  # Merge study and reference datasets
  executePlinkAd(
    ResultDir,
    c("--bfile", paste0(ResultDir, "/filtered_study_temp3"),
      "--bmerge", paste0(ResultDir, "/filtered_ref_temp6"),
      "--allow-no-sex",                                    ## Adding in 4.0
      "--make-bed", "--out", paste0(ResultDir, "/study_ref_merge"),
      "--silent")
  )


  # Check for strand inconsistency
  if (file.exists(paste0(ResultDir, "/study_ref_merge-merge.missnp"))) {
    # Re-merge study and reference datasets after excluding inconsistent SNPs
    executePlinkAd(
      ResultDir,
      c("--bfile", paste0(ResultDir, "/filtered_ref_temp6"),
        "--exclude", paste0(ResultDir, "/study_ref_merge-merge.missnp"),
        "--allow-no-sex",                                                ## Adding in 4.0
        "--make-bed", "--out", paste0(ResultDir, "/filtered_ref_temp7"), "--silent")
    )


    executePlinkAd(
      ResultDir,
      c("--bfile", paste0(ResultDir, "/filtered_study_temp3"),
        "--bmerge", paste0(ResultDir, "/filtered_ref_temp7"),
        "--allow-no-sex",                                     ## Adding in 4.0.
        "--make-bed", "--out", paste0(ResultDir, "/study_ref_merge"),
        "--silent")
    )

  }

  # Remove temporary files using the provided helper function
  removeTempFiles(ResultDir, "filtered_study_temp3")

  # Perform PCA for ancestry
  executePlinkAd(
    ResultDir,
    c("--bfile", paste0(ResultDir, "/study_ref_merge"),
      "--pca", "--out", paste0(ResultDir, "/study_ref_merge"),
      "--silent")
  )


  # Clean up temporary files
  removeTempFiles(ResultDir, "filtered_ref_temp6")
  removeTempFiles(ResultDir, "filtered_ref_temp7")

  print("PCA done.")
}

## Function 111
## Added in 3.0
loadAndProcessReferenceAncestry <- function(reference) {
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

  return(ref_ancestry_CEU_YRI_ASIAN)
}

## Function 112
## Added in 3.0
prepareAncestryData <- function(study_pop, ref_ancestry_CEU_YRI_ASIAN) {
  study_ancestry <- study_pop
  colnames(study_ancestry) <- c("ID", "Ancestry")
  study_ancestry$Ancestry <- paste0("Study_", study_ancestry$Ancestry)

  combined_pop <- rbind(ref_ancestry_CEU_YRI_ASIAN, study_ancestry)
  return(combined_pop)
}

## Function 113
## Added in 3.0
loadPCAData <- function(ResultDir, combined_pop) {
  pca <- read.table(file = paste0(ResultDir, "/study_ref_merge.eigenvec"), stringsAsFactors = FALSE, header = FALSE)
  tab <- data.frame(
    sample = pca$V2,
    pop = factor(combined_pop$Ancestry)[match(pca$V2, combined_pop$ID)],
    PC1 = pca[, 3], # the first eigenvector
    PC2 = pca[, 4], # the second eigenvector
    stringsAsFactors = FALSE
  )
  return(na.omit(tab))
}

## Function 114
## Added in 3.0
createPopulationTypeData <- function(tab) {
  tab <- na.omit(tab)
  #head(tab)
  pop_type <- as.data.frame(unique(tab$pop))
  colnames(pop_type) <- "type"
  pop_type$value <- 20
  pop_type[grepl("Study_.*", pop_type$type), "value"] <- 3

  # Similarly, for the reference population:
  pop_type_ref <- pop_type[grepl("Ref_.*", pop_type$type), ]
  pop_type_study <- pop_type[grepl("Study_.*", pop_type$type), ]

  # And then merging them
  pop_type <- rbind(pop_type_ref, pop_type_study)

  return(pop_type)
}

## Function 115
## Added in 3.0
plotPCA <- function(tab, pop_type) {
  p <- ggplot2::ggplot(data = tab, ggplot2::aes(
    x = tab$PC1, y = tab$PC2, color = tab$pop, shape = tab$pop
  )) +
    ggplot2::geom_hline(yintercept = 0, lty = 2) +
    ggplot2::geom_vline(xintercept = 0, lty = 2) +
    ggplot2::guides(color = ggplot2::guide_legend(title = "Ancestry"), shape = ggplot2::guide_legend(title = "Ancestry")) +
    ggplot2::scale_shape_manual(values = pop_type$value) +
    ggplot2::geom_point(alpha = 1, size = 4)

  print(p)
}

## Function 116
## Added in 3.0
reportAlleleFlips <- function(snp_allele_flips, ResultDir) {
  snp_allele_flips_count <- length(unique(snp_allele_flips))

  if (snp_allele_flips_count == 0) {
    print("No allele flips between study and reference data.")
  } else {
    print(paste0(snp_allele_flips_count, " allele flips identified between study and reference data."))
    mergebim <- read.table(paste0(ResultDir, "/study_ref_merge.bim"))
    print(paste0(length(unique(mergebim$V2)), " SNPs were finally retained in study and reference data after correcting for position mismatch and allele flips."))
  }
}

## Function 117
## Added in 3.0
detectOutliers <- function(tab, ResultDir, outlier, outlierOf, outlier_threshold) {
  if (!outlier) {
    return(data.frame())  # Return an empty data frame if outlier detection is not required
  }

  # Perform outlier detection
  ref_pop <- tab[tab$pop == paste0("Ref_", outlierOf), ]
  ref_pc1_median <- median(ref_pop$PC1)
  ref_pc2_median <- median(ref_pop$PC2)

  ref_pop$dis <- sqrt((ref_pop$PC1 - ref_pc1_median) ^ 2 + (ref_pop$PC2 - ref_pc2_median) ^ 2)
  ref_max_dis <- max(ref_pop$dis)

  study_pop <- tab[grepl("Study_.*", tab$pop), ]

  study_pop$dis <- sqrt((study_pop$PC1 - ref_pc1_median) ^ 2 + (study_pop$PC2 - ref_pc2_median) ^ 2)

  Non_Eu <- study_pop[study_pop$dis > ref_max_dis * outlier_threshold, ]
  Outlier_samples <- unique(Non_Eu$sample)

  # Report results
  if (length(Outlier_samples) != 0) {
    write.table(
      Outlier_samples,
      file = paste0(ResultDir, "/Outlier_ancestry"),
      quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\r\n", sep = " "
    )
    print(paste0(length(Outlier_samples), " Non-European samples as outlier."))
  } else {
    print("There is no non-European sample as outlier.")
  }

  # Prepare Outlier_samples1
  Outlier_samples1 <- Non_Eu[, 1:2]
  Outlier_samples1$pop <- as.character(Outlier_samples1$pop)
  Outlier_samples1 <- Outlier_samples1[grepl("Study_.*", Outlier_samples1$pop), ]
  return(Outlier_samples1)
}

## Function 118
## Added in 3.0
runSKAT <- function(score.file, gene.file, genes, cor.path, gene_approximation, anno.type, beta.par, weights.function, gen.var.weights, kernel_p_method, acc_devies, lim_devies, rho, skato_p_threshold) {
  sumFREGAT::SKAT(
    score.file = score.file,
    gene.file = gene.file,
    genes = genes,
    cor.path = cor.path,
    approximation = gene_approximation,
    anno.type = anno.type,
    beta.par = beta.par,
    weights.function = weights.function,
    user.weights = FALSE,
    gen.var.weights = gen.var.weights,
    method = kernel_p_method,
    acc = acc_devies,
    lim = lim_devies,
    rho = rho,
    p.threshold = skato_p_threshold,
    write.file = FALSE,
    quiet = TRUE
  )
}

## Function 119
## Added in 3.0
runSKATO <- function(score.file, gene.file, genes, cor.path, anno.type, gene_approximation, beta.par, weights.function, kernel_p_method, acc_devies, lim_devies, rho, skato_p_threshold) {
  sumFREGAT::SKATO(
    score.file = score.file,
    gene.file = gene.file,
    genes = genes,
    cor.path = cor.path,
    anno.type = anno.type,
    approximation = gene_approximation,
    beta.par = beta.par,
    weights.function = weights.function,
    user.weights = FALSE,
    method = kernel_p_method,
    acc = acc_devies,
    lim = lim_devies,
    rho = rho,
    p.threshold = skato_p_threshold,
    write.file = FALSE,
    quiet = TRUE
  )
}

## Function 120
## Added in 3.0
runSumChi <- function(score.file, gene.file, genes, cor.path, gene_approximation, anno.type, kernel_p_method, acc_devies, lim_devies) {
  sumFREGAT::sumchi(
    score.file = score.file,
    gene.file = gene.file,
    genes = genes,
    cor.path = cor.path,
    approximation = gene_approximation,
    anno.type = anno.type,
    method = kernel_p_method,
    acc = acc_devies,
    lim = lim_devies,
    write.file = FALSE,
    quiet = TRUE
  )
}

## Function 121
## Added in 3.0
runACAT <- function(score.file, gene.file, genes, anno.type, beta.par, weights.function, gen.var.weights, mac.threshold, sample_size) {
  sumFREGAT::ACAT(
    score.file = score.file,
    gene.file = gene.file,
    genes = genes,
    anno.type = anno.type,
    beta.par = beta.par,
    weights.function = weights.function,
    user.weights = FALSE,
    gen.var.weights = gen.var.weights,
    mac.threshold = mac.threshold,
    n = sample_size,
    write.file = FALSE,
    quiet = TRUE
  )
}

## Function 122
## Added in 3.0
runBT <- function(score.file, gene.file, genes, cor.path, anno.type, beta.par, weights.function) {
  sumFREGAT::BT(
    score.file = score.file,
    gene.file = gene.file,
    genes = genes,
    cor.path = cor.path,
    anno.type = anno.type,
    beta.par = beta.par,
    weights.function = weights.function,
    user.weights = FALSE,
    write.file = FALSE,
    quiet = TRUE
  )
}

## Function 123
## Added in 3.0
runPCA <- function(score.file, gene.file, genes, cor.path, gene_approximation, anno.type, sample_size, beta.par, weights.function, reference.matrix, regularize.fun, var.fraction) {
  sumFREGAT::PCA(
    score.file = score.file,
    gene.file = gene.file,
    genes = genes,
    cor.path = cor.path,
    approximation = gene_approximation,
    anno.type = anno.type,
    n = sample_size,
    beta.par = beta.par,
    weights.function = weights.function,
    user.weights = FALSE,
    reference.matrix = reference.matrix,
    fun = regularize.fun,
    var.fraction = var.fraction,
    write.file = FALSE,
    quiet = TRUE
  )
}

## Function 124
## Added in 3.0
runFLM <- function(score.file, gene.file, genes, cor.path, gene_approximation, anno_type, sample_size, beta.par, weights_function, flm_basis_function, flm_num_basis, flm_poly_order, flip_genotypes, Fan, reference_matrix_used, fun) {
  sumFREGAT::FLM(
    score.file = score.file,
    gene.file = gene.file,
    genes = genes,
    cor.path = cor.path,
    approximation = gene_approximation,
    anno.type = anno_type,
    n = sample_size,
    beta.par = beta.par,
    weights.function = weights_function,
    user.weights = FALSE,
    basis.function = flm_basis_function,
    k = flm_num_basis,
    order = flm_poly_order,
    flip.genotypes = flip_genotypes,
    Fan = Fan,
    reference.matrix = reference_matrix_used,
    fun = fun,
    write.file = FALSE,
    quiet = TRUE
  )
}

## Function 125
## Added in 3.0
runSimpleM <- function(score.file, gene.file, genes, anno_type, pca_var_fraction) {
  sumFREGAT::simpleM(
    score.file = score.file,
    gene.file = gene.file,
    genes = genes,
    anno.type = anno_type,
    var.fraction = pca_var_fraction,
    write.file = FALSE,
    quiet = TRUE
  )
}

## Function 126
## Added in 3.0
runMinP <- function(score.file, gene.file, genes, anno_type) {
  sumFREGAT::minp(
    score.file = score.file,
    gene.file = gene.file,
    genes = genes,
    anno.type = anno_type,
    write.file = FALSE,
    quiet = TRUE
  )
}

## Function 127
## Added in 3.0
validateInputDataSexDiff <- function(Mfile, Ffile) {
  requiredColumnsMfile <- c("SNP", "CHR", "BP", "A1", "BETA_M", "SE_M")
  requiredColumnsFfile <- c("SNP", "CHR", "BP", "A1", "BETA_F", "SE_F")

  # Check if Mfile and Ffile are dataframes
  if (!is.data.frame(Mfile)) {
    stop("Error: Mfile must be a dataframe.")
  }
  if (!is.data.frame(Ffile)) {
    stop("Error: Ffile must be a dataframe.")
  }

  # Check for existence of required columns in Mfile
  missingColsMfile <- setdiff(requiredColumnsMfile, names(Mfile))
  if (length(missingColsMfile) > 0) {
    stop(paste("Error: Mfile is missing these columns:", paste(missingColsMfile, collapse=", ")))
  }

  # Check for existence of required columns in Ffile
  missingColsFfile <- setdiff(requiredColumnsFfile, names(Ffile))
  if (length(missingColsFfile) > 0) {
    stop(paste("Error: Ffile is missing these columns:", paste(missingColsFfile, collapse=", ")))
  }

  return(TRUE)
}

## Function 128
## Added in 3.0
validateInputForQCsnp <- function(DataDir, ResultDir = tempdir(), finput, foutput, casecontrol = FALSE, hweCase = NULL, hweControl = NULL, hwe = NULL, maf = 0.05, geno = 0.1, monomorphicSNPs = FALSE, caldiffmiss = FALSE, diffmissFilter = FALSE, dmissX = FALSE, dmissAutoY = FALSE, highLD_regions, ld_prunning = FALSE, window_size = 50, step_size = 5, r2_threshold = 0.02) {

  # Validate directories
  if (!dir.exists(DataDir)) {
    stop("Error in DataDir: Directory does not exist.")
  }
  if (!dir.exists(ResultDir)) {
    stop("Error in ResultDir: Directory does not exist.")
  }

  # Validate file prefixes
  if (!is.character(finput)) {
    stop("Error in finput: Must be a character string.")
  }
  if (!is.character(foutput)) {
    stop("Error in foutput: Must be a character string.")
  }

  # Validate boolean parameters
  boolean_params <- list(casecontrol = casecontrol, monomorphicSNPs = monomorphicSNPs, caldiffmiss = caldiffmiss, diffmissFilter = diffmissFilter, dmissX = dmissX, dmissAutoY = dmissAutoY, ld_prunning = ld_prunning)
  for (param_name in names(boolean_params)) {
    param_value <- boolean_params[[param_name]]
    if (!is.logical(param_value)) {
      stop(paste("Error in", param_name, ": Must be a boolean value."))
    }
  }

  # Validate numeric parameters
  numeric_params <- list(hweCase = hweCase, hweControl = hweControl, hwe = hwe, maf = maf, geno = geno, r2_threshold = r2_threshold)
  for (param_name in names(numeric_params)) {
    param_value <- numeric_params[[param_name]]
    if (!is.null(param_value) && (!is.numeric(param_value) || param_value < 0 || param_value > 1)) {
      stop(paste("Error in", param_name, ": Must be a numeric value between 0 and 1."))
    }
  }

  # Validate highLD_regions if not NULL
  if (!is.null(highLD_regions) && !is.data.frame(highLD_regions)) {
    stop("Error in highLD_regions: Must be a dataframe.")
  }

  # Validate integer-like parameters
  if (!is.numeric(window_size) || window_size <= 0 || window_size != as.integer(window_size)) {
    stop("Error in window_size: Must be a positive whole number.")
  }
  if (!is.numeric(step_size) || step_size <= 0 || step_size != as.integer(step_size)) {
    stop("Error in step_size: Must be a positive whole number.")
  }

  return(TRUE)
}

## Function 129
## Added in 3.0
validateInputForEstimateHerit <- function(DataDir, ResultDir = tempdir(), finput, summarystat = NULL, ncores = parallel::detectCores(), model = c("LDSC", "GREML"), byCHR = FALSE, r2_LD = 0, LDSC_blocks = 20, REMLalgo = c(0, 1, 2), nitr = 100, cat_covarfile = NULL, quant_covarfile = NULL, prevalance = 0.01, partGRM = FALSE, autosome = TRUE, Xsome = TRUE, nGRM = 3, cripticut = 0.025, minMAF = NULL, maxMAF = NULL, hg = c("hg19", "hg38"), PlotIndepSNP = c(TRUE, FALSE), IndepSNP_window_size = 50, IndepSNP_step_size = 5, IndepSNP_r2_threshold = 0.02, highLD_regions = NULL) {

  # Validate directories
  if (!dir.exists(DataDir)) {
    stop("Error in DataDir: Directory does not exist.")
  }
  if (!dir.exists(ResultDir)) {
    stop("Error in ResultDir: Directory does not exist.")
  }

  # Validate file prefix
  if (!is.character(finput)) {
    stop("Error in finput: Must be a character string.")
  }

  # Validate summarystat if not NULL
  if (!is.null(summarystat) && !is.data.frame(summarystat)) {
    stop("Error in summarystat: Must be a dataframe.")
  }

  # Validate ncores
  if (!is.numeric(ncores) || ncores < 0) {
    stop("Error in ncores: Must be a non-negative integer.")
  }

  # Validate model
  if (!is.character(model) || !all(model %in% c("LDSC", "GREML"))) {
    stop("Error in model: Must be either 'LDSC' or 'GREML', or both.")
  }

  # Validate boolean parameters
  if (!is.logical(byCHR) || !is.logical(partGRM) || !is.logical(autosome) || !is.logical(Xsome) || !is.logical(PlotIndepSNP)) {
    stop("Error in Boolean parameters: Must be TRUE or FALSE.")
  }

  # Validate numerical parameters
  if (!is.numeric(r2_LD) || r2_LD < 0) {
    stop("Error in r2_LD: Must be a non-negative number.")
  }
  if (!is.numeric(prevalance) || prevalance < 0) {
    stop("Error in prevalance: Must be a non-negative number.")
  }
  if (!is.numeric(cripticut) || cripticut < 0) {
    stop("Error in cripticut: Must be a non-negative number.")
  }

  # Validate integer-like parameters
  if (!is.null(LDSC_blocks) && (!is.numeric(LDSC_blocks) || LDSC_blocks <= 0 || LDSC_blocks != as.integer(LDSC_blocks))) {
    stop("Error in LDSC_blocks: Must be a positive whole number.")
  }
  if (!is.null(nitr) && (!is.numeric(nitr) || nitr <= 0 || nitr != as.integer(nitr))) {
    stop("Error in nitr: Must be a positive whole number.")
  }
  if (!is.null(nGRM) && (!is.numeric(nGRM) || nGRM <= 0 || nGRM != as.integer(nGRM))) {
    stop("Error in nGRM: Must be a positive whole number.")
  }
  if (!is.numeric(IndepSNP_window_size) || IndepSNP_window_size <= 0 || IndepSNP_window_size != as.integer(IndepSNP_window_size)) {
    stop("Error in IndepSNP_window_size: Must be a positive whole number.")
  }
  if (!is.numeric(IndepSNP_step_size) || IndepSNP_step_size <= 0 || IndepSNP_step_size != as.integer(IndepSNP_step_size)) {
    stop("Error in IndepSNP_step_size: Must be a positive whole number.")
  }


  # Validate REMLalgo
  if (!is.numeric(REMLalgo) || !all(REMLalgo %in% c(0, 1, 2))) {
    stop("Error in REMLalgo: Must be 0, 1, or 2.")
  }

  # Validate MAF parameters
  if ((!is.null(minMAF) && (!is.numeric(minMAF) || minMAF < 0 || minMAF > 1)) || (!is.null(maxMAF) && (!is.numeric(maxMAF) || maxMAF < 0 || maxMAF > 1))) {
    stop("Error in MAF parameters: Must be within the range [0, 1].")
  }

  # Validate hg parameter
  if (!is.character(hg) || !all(hg %in% c("hg19", "hg38"))) {
    stop("Error in hg: Must be 'hg19', 'hg38', or both.")
  }

  # Validate highLD_regions if not NULL
  if (!is.null(highLD_regions) && !is.data.frame(highLD_regions)) {
    stop("Error in highLD_regions: Must be a dataframe.")
  }

  return(TRUE)
}

## Function 130
## Added in 3.0
validateInputForComputeGeneticPC <- function(DataDir, ResultDir = tempdir(), finput, countPC = 10, plotPC = TRUE, highLD_regions = NULL, ld_prunning = TRUE, window_size = 50, step_size = 5, r2_threshold = 0.02) {

  # Validate directories
  if (!dir.exists(DataDir)) {
    stop("Error in DataDir: Directory does not exist.")
  }
  if (!dir.exists(ResultDir)) {
    stop("Error in ResultDir: Directory does not exist.")
  }

  # Validate file prefix
  if (!is.character(finput)) {
    stop("Error in finput: Must be a character string.")
  }

  # Validate countPC
  if (!is.numeric(countPC) || countPC <= 0 || countPC != as.integer(countPC)) {
    stop("Error in countPC: Must be a positive whole number.")
  }

  # Validate plotPC
  if (!is.logical(plotPC)) {
    stop("Error in plotPC: Must be a boolean value.")
  }

  # Validate ld_prunning
  if (!is.logical(ld_prunning)) {
    stop("Error in ld_prunning: Must be a boolean value.")
  }

  # Validate highLD_regions if not NULL
  if (!is.null(highLD_regions) && !is.data.frame(highLD_regions)) {
    stop("Error in highLD_regions: Must be a dataframe.")
  }

  # Validate integer-like parameters for window_size and step_size
  if (!is.numeric(window_size) || window_size <= 0 || (window_size != as.integer(window_size) && window_size != 50)) {
    stop("Error in window_size: Must be a positive whole number or the default value of 50.")
  }
  if (!is.numeric(step_size) || step_size <= 0 || (step_size != as.integer(step_size) && step_size != 5)) {
    stop("Error in step_size: Must be a positive whole number or the default value of 5.")
  }

  # Validate r2_threshold
  if (!is.numeric(r2_threshold) || r2_threshold < 0 || r2_threshold > 1) {
    stop("Error in r2_threshold: Must be a numeric value between 0 and 1.")
  }

  return(TRUE)
}

## Function 131
## Added in 3.0
validateInputForComputePRS <- function(DataDir, ResultDir = tempdir(), finput, summarystat, phenofile, covarfile = NULL, effectsize = c("BETA", "OR"), ldclump = FALSE, LDreference, clump_p1 = 0.0001, clump_p2 = 0.01, clump_r2 = 0.50, clump_kb = 250, byCHR = TRUE, pthreshold = c(0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5), highLD_regions = "high-LD-regions-hg19-GRCh37.txt", ld_prunning = FALSE, window_size = 50, step_size = 5, r2_threshold = 0.02, nPC = 6, pheno_type = "binary") {

  # Validate directories
  if (!dir.exists(DataDir)) {
    stop("Error in DataDir: Directory does not exist.")
  }
  if (!dir.exists(ResultDir)) {
    stop("Error in ResultDir: Directory does not exist.")
  }

  # Validate file prefix
  if (!is.character(finput)) {
    stop("Error in finput: Must be a character string.")
  }

  if (!is.data.frame(summarystat)) {
    stop("Error in summarystat: Must be a dataframe.")
  }

  # Check for required column names
  requiredColumns <- c("SNP", "A1", c("BETA", "OR"))
  actualColumns <- colnames(summarystat)[1:3]

  # Check if the first two columns are "SNP" and "A1"
  if (!all(requiredColumns[1:2] %in% actualColumns[1:2])) {
    stop("Error in summarystat: The first two columns must be 'SNP' and 'A1'.")
  }

  # Check if the third column is either "BETA" or "OR"
  if (!actualColumns[3] %in% requiredColumns[3]) {
    stop("Error in summarystat: The third column must be either 'BETA' or 'OR'.")
  }

  # Validate phenofile
  if (!is.data.frame(phenofile)) {
    stop("Error in phenofile: Must be a dataframe..")
  }

  # Validate covarfile if not NULL
  if (!is.data.frame(covarfile) && !is.character(covarfile)) {
    stop("Error in covarfile: Must be a dataframe..")
  }

  # Validate effectsize
  if (!effectsize %in% c("BETA", "OR")) {
    stop("Error in effectsize: Must be 'BETA' or 'OR'.")
  }

  # Validate ldclump
  if (!is.logical(ldclump)) {
    stop("Error in ldclump: Must be a boolean value.")
  }

  # Validate LDreference if ldclump is TRUE
  if (ldclump && !is.character(LDreference)) {
    stop("Error in LDreference: Must be a character string.")
  }

  # Validate numeric parameters
  numeric_params <- list(clump_p1 = clump_p1, clump_p2 = clump_p2, clump_r2 = clump_r2, r2_threshold = r2_threshold)
  for (param_name in names(numeric_params)) {
    param_value <- numeric_params[[param_name]]
    if (!is.numeric(param_value) || param_value < 0 || param_value > 1) {
      stop(paste("Error in", param_name, ": Must be a numeric value between 0 and 1."))
    }
  }

  # Validate clump_kb
  if (!is.numeric(clump_kb) || clump_kb <= 0 || clump_kb != as.integer(clump_kb)) {
    stop("Error in clump_kb: Must be a positive whole number.")
  }

  # Validate byCHR
  if (!is.logical(byCHR)) {
    stop("Error in byCHR: Must be a boolean value.")
  }

  # Validate pthreshold
  if (!is.vector(pthreshold) || any(pthreshold < 0) || any(pthreshold > 1)) {
    stop("Error in pthreshold: Must be a numeric vector with values between 0 and 1.")
  }

  # Validate ld_prunning
  if (!is.logical(ld_prunning)) {
    stop("Error in ld_prunning: Must be a boolean value.")
  }

  # Validate integer-like parameters for window_size and step_size
  if (!is.numeric(window_size) || window_size <= 0 || (window_size != as.integer(window_size) && window_size != 50)) {
    stop("Error in window_size: Must be a positive whole number or the default value of 50.")
  }
  if (!is.numeric(step_size) || step_size <= 0 || (step_size != as.integer(step_size) && step_size != 5)) {
    stop("Error in step_size: Must be a positive whole number or the default value of 5.")
  }

  # Validate nPC
  if (!is.numeric(nPC) || nPC <= 0 || nPC != as.integer(nPC)) {
    stop("Error in nPC: Must be a positive whole number.")
  }

  # Validate pheno_type
  if (!pheno_type %in% c("binary", "quantitative")) {
    stop("Error in pheno_type: Must be 'binary' or 'quantitative'.")
  }

  return(TRUE)
}

## Function 132
## Added in 3.0
validateInputForMergeRegion <- function(DataDir, ResultDir = tempdir(), finput1, finput2, foutput, use_common_snps = TRUE) {

  # Validate directories
  if (!dir.exists(DataDir)) {
    stop("Error in DataDir: Directory does not exist.")
  }
  if (!dir.exists(ResultDir)) {
    stop("Error in ResultDir: Directory does not exist.")
  }

  # Validate file prefixes
  if (!is.character(finput1)) {
    stop("Error in finput1: Must be a character string.")
  }
  if (!is.character(finput2)) {
    stop("Error in finput2: Must be a character string.")
  }
  if (!is.character(foutput)) {
    stop("Error in foutput: Must be a character string.")
  }

  # Validate use_common_snps
  if (!is.logical(use_common_snps)) {
    stop("Error in use_common_snps: Must be a boolean value.")
  }

  return(TRUE)
}

## Function 133
## Added in 3.0
validateInputForPlinkVCF <- function(DataDir, ResultDir = tempdir(), finput, foutput, VtoP = TRUE, PtoV = TRUE, Famfile = NULL, PVbyCHR = TRUE) {

  # Validate directories
  if (!dir.exists(DataDir)) {
    stop("Error in DataDir: Directory does not exist.")
  }
  if (!dir.exists(ResultDir)) {
    stop("Error in ResultDir: Directory does not exist.")
  }

  # Validate file prefix
  if (!is.character(finput)) {
    stop("Error in finput: Must be a character string.")
  }
  if (!is.character(foutput)) {
    stop("Error in foutput: Must be a character string.")
  }

  # Validate VtoP and PtoV
  if (!is.logical(VtoP)) {
    stop("Error in VtoP: Must be a boolean value.")
  }
  if (!is.logical(PtoV)) {
    stop("Error in PtoV: Must be a boolean value.")
  }

  # Validate Famfile if not NULL
  if (!is.null(Famfile) && !is.character(Famfile)) {
    stop("Error in Famfile: Must be a character string.")
  }

  # Validate PVbyCHR
  if (!is.logical(PVbyCHR)) {
    stop("Error in PVbyCHR: Must be a boolean value.")
  }

  return(TRUE)
}

## Function 134
## Added in 3.0
validateInputForSexCheck <- function(DataDir, ResultDir = tempdir(), finput, impute_sex = FALSE, compute_freq = FALSE, LD = TRUE, LD_window_size = 50, LD_step_size = 5, LD_r2_threshold = 0.02, fmax_F = 0.2, mmin_F = 0.8) {

  # Validate directories
  if (!dir.exists(DataDir)) {
    stop("Error in DataDir: Directory does not exist.")
  }
  if (!dir.exists(ResultDir)) {
    stop("Error in ResultDir: Directory does not exist.")
  }

  # Validate file prefix
  if (!is.character(finput)) {
    stop("Error in finput: Must be a character string.")
  }

  # Validate boolean parameters
  boolean_params <- list(impute_sex = impute_sex, compute_freq = compute_freq, LD = LD)
  for (param_name in names(boolean_params)) {
    param_value <- boolean_params[[param_name]]
    if (!is.logical(param_value)) {
      stop(paste("Error in", param_name, ": Must be a boolean value."))
    }
  }

  # Validate integer-like parameters for LD_window_size and LD_step_size
  if (!is.numeric(LD_window_size) || LD_window_size <= 0 || (LD_window_size != as.integer(LD_window_size) && LD_window_size != 50)) {
    stop("Error in LD_window_size: Must be a positive whole number or the default value of 50.")
  }
  if (!is.numeric(LD_step_size) || LD_step_size <= 0 || (LD_step_size != as.integer(LD_step_size) && LD_step_size != 5)) {
    stop("Error in LD_step_size: Must be a positive whole number or the default value of 5.")
  }

  # Validate numeric parameters
  numeric_params <- list(LD_r2_threshold = LD_r2_threshold, fmax_F = fmax_F, mmin_F = mmin_F)
  for (param_name in names(numeric_params)) {
    param_value <- numeric_params[[param_name]]
    if (!is.numeric(param_value) || param_value < 0 || param_value > 1) {
      stop(paste("Error in", param_name, ": Must be a numeric value between 0 and 1."))
    }
  }

  return(TRUE)
}

## Function 135
## Added in 3.0
validateInputForFilterPlinkSample <- function(DataDir, ResultDir = tempdir(), finput, foutput, filter_sample = "cases", keep_remove_sample_file = NULL, keep = TRUE) {

  # Validate directories
  if (!dir.exists(DataDir)) {
    stop("Error in DataDir: Directory does not exist.")
  }
  if (!dir.exists(ResultDir)) {
    stop("Error in ResultDir: Directory does not exist.")
  }

  # Validate file prefixes
  if (!is.character(finput)) {
    stop("Error in finput: Must be a character string.")
  }
  if (!is.null(foutput) && !is.character(foutput)) {
    stop("Error in foutput: Must be a character string.")
  }

  # Validate filter_sample
  valid_filters <- c("cases", "controls", "males", "females")
  if (!is.character(filter_sample) || !filter_sample %in% valid_filters) {
    stop("Error in filter_sample: Must be one of 'cases', 'controls', 'males', or 'females'.")
  }

  # Validate filter_sample
  valid_filters <- c("cases", "controls", "males", "females")
  if (!is.character(filter_sample) || length(filter_sample) != 1 || !filter_sample %in% valid_filters) {
    stop("Error in filter_sample: Must be a single value, one of 'cases', 'controls', 'males', or 'females'.")
  }

  # Validate keep
  if (!is.logical(keep)) {
    stop("Error in keep: Must be a boolean value.")
  }

  return(TRUE)
}

## Function 136
## Added in 3.0
validateInputForGetMFPlink <- function(DataDir, ResultDir = tempdir(), finput, foutput, sex, xplink = FALSE, autoplink = FALSE) {

  # Validate directories
  if (!dir.exists(DataDir)) {
    stop("Error in DataDir: Directory does not exist.")
  }
  if (!dir.exists(ResultDir)) {
    stop("Error in ResultDir: Directory does not exist.")
  }

  # Validate file prefix
  if (!is.character(finput)) {
    stop("Error in finput: Must be a character string.")
  }
  if (!is.character(foutput)) {
    stop("Error in foutput: Must be a character string.")
  }

  # Validate sex
  valid_sexes <- c("males", "females")
  if (!is.character(sex) || !sex %in% valid_sexes) {
    stop("Error in sex: Must be 'males' or 'females'.")
  }

  # Validate xplink and autoplink
  if (!is.logical(xplink)) {
    stop("Error in xplink: Must be a boolean value.")
  }
  if (!is.logical(autoplink)) {
    stop("Error in autoplink: Must be a boolean value.")
  }

  return(TRUE)
}

## Function 136
## Added in 3.0
validateInputForXhwe <- function(DataDir, ResultDir = tempdir(), finput, foutput, filterSNP = TRUE) {

  # Validate directories
  if (!dir.exists(DataDir)) {
    stop("Error in DataDir: Directory does not exist.")
  }
  if (!dir.exists(ResultDir)) {
    stop("Error in ResultDir: Directory does not exist.")
  }

  # Validate file prefixes
  if (!is.character(finput)) {
    stop("Error in finput: Must be a character string.")
  }
  if (!is.character(foutput)) {
    stop("Error in foutput: Must be a character string.")
  }

  # Validate filterSNP
  if (!is.logical(filterSNP)) {
    stop("Error in filterSNP: Must be a boolean value.")
  }

  return(TRUE)
}

## Function 137
## Added in 3.0
validateInputForMAFdiffSexControl <- function(DataDir, ResultDir = tempdir(), finput, filterSNP = FALSE, foutput = NULL) {

  # Validate directories
  if (!dir.exists(DataDir)) {
    stop("Error in DataDir: Directory does not exist.")
  }
  if (!dir.exists(ResultDir)) {
    stop("Error in ResultDir: Directory does not exist.")
  }

  # Validate file prefix
  if (!is.character(finput)) {
    stop("Error in finput: Must be a character string.")
  }

  # Validate foutput if not NULL
  if (!is.null(foutput) && !is.character(foutput)) {
    stop("Error in foutput: Must be a character string.")
  }

  # Validate filterSNP
  if (!is.logical(filterSNP)) {
    stop("Error in filterSNP: Must be a boolean value.")
  }

  return(TRUE)
}

## Function 138
## Added in 3.0
validateInputForQCsample <- function(DataDir, ResultDir, finput, foutput, imiss, het, small_sample_mod, IBD, IBDmatrix, ambi_out, legend_text_size, legend_title_size, axis_text_size, axis_title_size, title_size, filterSample) {
  # Validation for directory paths (should be strings)
  if (!is.character(DataDir) || !is.character(ResultDir)) {
    stop("DataDir and ResultDir must be strings representing directory paths.")
  }

  # Validation for file names (should be strings)
  if (!is.character(finput) || (!is.null(foutput) && !is.character(foutput))) {
    stop("finput and foutput must be strings representing file names.")
  }

  # Validation for numeric parameters
  numeric_params <- list(imiss = imiss, het = het)

  # Validate 'imiss' to be between 0 and 1
  param_value_imiss <- numeric_params[['imiss']]
  if (!is.numeric(param_value_imiss) || param_value_imiss < 0 || param_value_imiss > 1) {
    stop("imiss must be a numeric value between 0 and 1.")
  }

  # Validate 'het' to be numeric
  param_value_het <- numeric_params[['het']]
  if (!is.null(param_value_het) && !is.numeric(param_value_het)) {
    stop("het must be a numeric value or NULL.")
  }

  # Validation for size parameters (should be positive integers)
  size_params <- list(legend_text_size = legend_text_size, legend_title_size = legend_title_size, axis_text_size = axis_text_size, axis_title_size = axis_title_size, title_size = title_size)
  for (param_name in names(size_params)) {
    if (!is.numeric(size_params[[param_name]]) || size_params[[param_name]] <= 0) {
      stop(paste0(param_name, " must be a positive integer."))
    }
  }
}

## Function 139
## Added in 3.0
validateInputForGXWASmiami <- function(ResultDir, FemaleWAS, MaleWAS, snp_pval, Xchr) {

  # Validate ResultDir
  if (!dir.exists(ResultDir)) {
    stop("Error in ResultDir: Directory does not exist.")
  }

  # Validate FemaleWAS and MaleWAS dataframes
  required_cols <- c("SNP", "CHR", "POS", "pvalue")
  if (!is.data.frame(FemaleWAS) || !all(required_cols %in% names(FemaleWAS))) {
    stop("Error in FemaleWAS: Must be a dataframe with columns SNP, CHR, POS, pvalue.")
  }
  if (!is.data.frame(MaleWAS) || !all(required_cols %in% names(MaleWAS))) {
    stop("Error in MaleWAS: Must be a dataframe with columns SNP, CHR, POS, pvalue.")
  }

  # Validate snp_pval
  if (!is.numeric(snp_pval) || length(snp_pval) != 1 || snp_pval < 0 || snp_pval > 1) {
    stop("Error in snp_pval: Must be a numeric value between 0 and 1.")
  }

  # Validate Xchr
  if (!is.logical(Xchr)) {
    stop("Error in Xchr: Must be a boolean value.")
  }

  return(TRUE)
}

## Function 140
## Added in 3.0
validateInputForMetaGWAS <- function(DataDir, ResultDir, SummData, SNPfile, useSNPposition, UseA1, GCse, plotname, pval_filter, top_snp_pval, max_top_snps, chosen_snps_file, byCHR, pval_threshold_manplot) {

  # Validate DataDir
  if (!dir.exists(DataDir)) {
    stop("Error in DataDir: Directory does not exist.")
  }

  # Validate ResultDir
  if (!dir.exists(ResultDir)) {
    stop("Error in ResultDir: Directory does not exist.")
  }


  # Validate SummData - must be a list of dataframes
  if (!is.list(SummData)) {
    stop("Error in SummData: Must be a list of dataframes.")
  }

  required_cols <- c("SNP", "BETA", "SE", "P", "NMISS")

  # Loop through each dataframe in SummData
  for (dataframe in SummData) {
    if (!is.data.frame(dataframe)) {
      stop("Error in SummData: Each element must be a dataframe.")
    }

    # Check if the dataframe contains required columns
    if (!all(required_cols %in% names(dataframe))) {
      stop("Error in SummData: Each dataframe must contain columns SNP, BETA, SE, P, NMISS.")
    }
  }

  # Validate SNPfile - should be NULL or a character string
  if (!is.null(SNPfile) && !is.character(SNPfile)) {
    stop("Error in SNPfile: Must be a character string.")
  }

  # Validate logical parameters
  if (!is.logical(useSNPposition) || !is.logical(UseA1) || !is.logical(GCse) || !is.logical(byCHR)) {
    stop("Logical parameters (useSNPposition, UseA1, GCse, byCHR) must be boolean values.")
  }

  # Validate plotname and pval_filter - must be character strings
  if (!is.character(plotname) || !is.character(pval_filter)) {
    stop("Error in plotname/pval_filter: Must be a character string.")
  }

  # Validate pval parameters - must be numeric between 0 and 1
  pval_params <- list(top_snp_pval = top_snp_pval, pval_threshold_manplot = pval_threshold_manplot)
  for (param_name in names(pval_params)) {
    param_value <- pval_params[[param_name]]
    if (!is.numeric(param_value) || param_value < 0 || param_value > 1) {
      stop(paste0("Error in ", param_name, ": Must be a numeric value between 0 and 1."))
    }
  }

  # Validate max_top_snps - must be a positive integer
  if (!is.numeric(max_top_snps) || max_top_snps != as.integer(max_top_snps) || max_top_snps <= 0) {
    stop("Error in max_top_snps: Must be a positive integer.")
  }

  # Validate chosen_snps_file - should be NULL or a character string
  if (!is.null(chosen_snps_file) && !is.character(chosen_snps_file)) {
    stop("Error in chosen_snps_file: Must be a character string.")
  }

  return(TRUE)
}

## Function 141
## Added in 3.0
validateInputForClumpLD <- function(DataDir, finput, SNPdata, ResultDir, clump_p1, clump_p2, clump_r2, clump_kb, byCHR) {

  # Validate DataDir
  if (!dir.exists(DataDir)) {
    stop("Error in DataDir: Directory does not exist.")
  }

  # Validate ResultDir
  if (!dir.exists(ResultDir)) {
    stop("Error in ResultDir: Directory does not exist.")
  }

  # Validate finput - should be a character string
  if (!is.character(finput)) {
    stop("Error in finput: Must be a character string.")
  }

  # Validate SNPdata - must be a list of dataframes
  if (!is.list(SNPdata) || !all(sapply(SNPdata, is.data.frame))) {
    stop("Error in SNPdata: Must be a list of dataframes.")
  }

  # Validate clump parameters - must be numeric
  clump_params <- list(clump_p1 = clump_p1, clump_p2 = clump_p2, clump_r2 = clump_r2, clump_kb = clump_kb)
  for (param_name in names(clump_params)) {
    param_value <- clump_params[[param_name]]
    if (!is.numeric(param_value)) {
      stop(paste0("Error in ", param_name, ": Must be a numeric value."))
    }
  }

  # Validate byCHR - should be a logical value
  if (!is.logical(byCHR)) {
    stop("Error in byCHR: Must be a boolean value.")
  }

  return(TRUE)
}

## Function 142
## Added in 3.0
validateInputForDiffZeroOne <- function(inputdata, diffzero, diffone) {

  # Validate inputdata - must be a dataframe
  if (!is.data.frame(inputdata)) {
    stop("Error in inputdata: Must be a dataframe.")
  }

  # Validate required columns in inputdata
  required_cols <- c("Stat", "SE")
  if (!all(required_cols %in% names(inputdata))) {
    stop("Error in inputdata: Dataframe must contain columns 'Stat' and 'SE'.")
  }

  # Validate diffzero - must be a single boolean value
  if (!is.logical(diffzero) || length(diffzero) != 1) {
    stop("Error in diffzero: Must be a single boolean value.")
  }

  # Validate diffone - must be a single boolean value
  if (!is.logical(diffone) || length(diffone) != 1) {
    stop("Error in diffone: Must be a single boolean value.")
  }

  return(TRUE)
}


## Function 143
## Added in 3.0
validateInputForSexDiffZcrore <- function(inputdata) {

  # Validate inputdata - must be a dataframe
  if (!is.data.frame(inputdata)) {
    stop("Error in inputdata: Must be a dataframe.")
  }

  # Validate required columns in inputdata
  required_cols <- c("Fstat", "Fse", "Mstat", "Mse")
  if (!all(required_cols %in% names(inputdata))) {
    stop("Error in inputdata: Dataframe must contain columns 'Fstat', 'Fse', 'Mstat', 'Mse'.")
  }

  return(TRUE)
}

## Function 144
## Added in 3.0
validateInputForGeneticCorrBT <- function(DataDir, ResultDir, finput, byCHR, REMLalgo, nitr, phenofile, cat_covarfile, quant_covarfile, partGRM, autosome, Xsome, nGRM, cripticut, minMAF, maxMAF, excludeResidual, ncores) {

  # Validate DataDir and ResultDir
  if (!dir.exists(DataDir)) {
    stop("Error in DataDir: Directory does not exist.")
  }
  if (!dir.exists(ResultDir)) {
    stop("Error in ResultDir: Directory does not exist.")
  }

  # Validate finput
  if (!is.character(finput)) {
    stop("Error in finput: Must be a character string.")
  }

  # Validate phenofile - must be a dataframe with exactly four columns
  if (!is.null(phenofile)) {
    if (!is.data.frame(phenofile)) {
      stop("Error in phenofile: Must be a dataframe.")
    }
    if (ncol(phenofile) != 4) {
      stop("Error in phenofile: Dataframe must contain exactly four columns.")
    }
  }


  # Validate categorical and quantitative covariate files
  if (!is.null(cat_covarfile) && !file.exists(paste0(DataDir, "/", cat_covarfile))) {
    stop("Error in cat_covarfile: Specified file does not exist in DataDir.")
  }
  if (!is.null(quant_covarfile) && !file.exists(paste0(DataDir, "/", quant_covarfile))) {
    stop("Error in quant_covarfile: Specified file does not exist in DataDir.")
  }

  # Validate boolean parameters
  if (!is.logical(byCHR) || !is.logical(partGRM) || !is.logical(autosome) || !is.logical(Xsome) || !is.logical(excludeResidual)) {
    stop("Error in boolean parameters: All must be TRUE or FALSE.")
  }

  # Validate REMLalgo
  if (!is.numeric(REMLalgo) || !all(REMLalgo %in% c(0, 1, 2))) {
    stop("Error in REMLalgo: Must be 0, 1, or 2.")
  }

  # Validate nitr
  if (!is.numeric(nitr) || nitr <= 0 || nitr != as.integer(nitr)) {
    stop("Error in nitr: Must be a positive whole number.")
  }

  # Validate nGRM
  if (!is.numeric(nGRM) || nGRM <= 0 || nGRM != as.integer(nGRM)) {
    stop("Error in nGRM: Must be a positive whole number.")
  }

  # Validate cripticut, minMAF, maxMAF
  if (!is.null(cripticut) && (!is.numeric(cripticut) || cripticut < 0 || cripticut > 1)) {
    stop("Error in cripticut: Must be within the range [0, 1].")
  }
  if (!is.null(minMAF) && (!is.numeric(minMAF) || minMAF < 0 || minMAF > 1)) {
    stop("Error in minMAF: Must be within the range [0, 1].")
  }
  if (!is.null(maxMAF) && (!is.numeric(maxMAF) || maxMAF < 0 || maxMAF > 1)) {
    stop("Error in maxMAF: Must be within the range [0, 1].")
  }

  # Validate ncores
  if (!is.numeric(ncores) || ncores <= 0 || ncores != as.integer(ncores)) {
    stop("Error in ncores: Must be a positive whole number.")
  }

  return(TRUE)
}

## Function 145
## Added in 3.0
validateInputForSexRegress <- function(fdata, regressor_index, response_index) {

  # Validate fdata - must be a dataframe
  if (!is.data.frame(fdata)) {
    stop("Error in fdata: Must be a dataframe.")
  }

  # Validate response_index and regressor_index - must be numeric and within the column index range of fdata
  if (!is.numeric(response_index) || any(response_index < 1) || any(response_index > ncol(fdata))) {
    stop("Error in response_index: Must be a numeric value within the column index range of fdata.")
  }

  if (!is.numeric(regressor_index) || any(regressor_index < 1) || any(regressor_index > ncol(fdata))) {
    stop("Error in regressor_index: Must be a numeric value within the column index range of fdata.")
  }

  # Validate that response_index and regressor_index are not overlapping
  if (any(regressor_index %in% response_index)) {
    stop("Error: regressor_index and response_index must refer to different columns.")
  }

  return(TRUE)
}

## Function 146
## Added in 3.0
validateFilterRegionParams <- function(DataDir, ResultDir, finput, foutput, CHRX, CHRY, filterPAR, filterXTR, filterAmpliconic, regionfile, filterCHR, Hg, exclude) {
  # Validate DataDir
  if (!is.character(DataDir) || !dir.exists(DataDir)) {
    stop("DataDir must be a valid directory path.")
  }

  # Validate ResultDir
  if (!is.character(ResultDir) || (!dir.exists(ResultDir) && ResultDir != tempdir())) {
    stop("ResultDir must be a valid directory path or the default tempdir().")
  }

  # Validate finput and foutput
  if (!is.character(finput) || !is.character(foutput)) {
    stop("finput and foutput must be character strings.")
  }

  # Validate boolean parameters
  boolean_params <- list(CHRX = CHRX, CHRY = CHRY, filterPAR = filterPAR, filterXTR = filterXTR, filterAmpliconic = filterAmpliconic, exclude = exclude)
  for (param in names(boolean_params)) {
    if (!is.logical(boolean_params[[param]]) || length(boolean_params[[param]]) != 1) {
      stop(sprintf("%s must be a single boolean value.", param))
    }
  }

  # Validate regionfile
  if (!is.logical(regionfile) && !is.character(regionfile)) {
    stop("regionfile must be either a boolean or a character string.")
  }

  # Validate filterCHR
  if (!is.null(filterCHR)) {
    if (is.numeric(filterCHR)) {
      if (any(filterCHR < 1) || any(filterCHR > 30) || any(filterCHR != as.integer(filterCHR))) {
        stop("filterCHR must contain positive integers between 1 and 24.")
      }
    } else {
      stop("filterCHR must be a numeric vector (single positive integer or array of positive integers between 1 and 24).")
    }
  }

  # Validate Hg
  if (!is.character(Hg) || !Hg %in% c("19", "38")) {
    stop("Hg must be '19' or '38'.")
  }

  # Parameters are valid
  return(TRUE)
}

## Function 147
## Added in 3.0
validatePvalCombInputs <- function(SumstatMale, SumstatFemale, combtest, MF.p.corr, MF.zero.sub, MF.na.rm, MF.mc.cores, B, plot.jpeg, plotname, snp_pval, annotateTopSnp, suggestiveline, genomewideline, ncores) {
  # Validate SumstatMale and SumstatFemale
  if (!is.data.frame(SumstatMale) || !is.data.frame(SumstatFemale)) {
    return("SumstatMale and SumstatFemale must be data frames.")
  }

  # Validate combtest
  if (!combtest %in% c("stouffer.method", "fisher.method", "fisher.method.perm")) {
    return("Invalid combtest. Choose 'stouffer.method', 'fisher.method', or 'fisher.method.perm'.")
  }

  # Validate numeric parameters
  numeric_params <- list(MF.zero.sub = MF.zero.sub, B = B, snp_pval = snp_pval, suggestiveline = suggestiveline, genomewideline = genomewideline, ncores = ncores)
  for (param in names(numeric_params)) {
    if (!is.numeric(numeric_params[[param]]) || numeric_params[[param]] <= 0) {
      return(sprintf("%s must be a positive numeric value.", param))
    }
  }

  # Validate boolean parameters
  boolean_params <- list(MF.na.rm = MF.na.rm, plot.jpeg = plot.jpeg, annotateTopSnp = annotateTopSnp)
  for (param in names(boolean_params)) {
    if (!is.logical(boolean_params[[param]]) || length(boolean_params[[param]]) != 1) {
      return(sprintf("%s must be a single boolean value.", param))
    }
  }

  # Validate MF.mc.cores
  if (!is.null(MF.mc.cores) && (!is.numeric(MF.mc.cores) || MF.mc.cores < 1)) {
    return("MF.mc.cores must be a positive integer.")
  }

  # Validate plotname
  if (!is.character(plotname) || nchar(plotname) == 0) {
    return("plotname must be a non-empty character string.")
  }

  # Return NULL if all validations pass
  return(NULL)
}

## Function 148
## Added in 3.0
validateFilterSNPInputs <- function(DataDir, finput, SNPvec, extract) {
  # Validate DataDir
  if (!is.character(DataDir) || !dir.exists(DataDir)) {
    return("DataDir must be a valid directory path.")
  }

  # Validate finput
  if (!is.character(finput) || !file.exists(file.path(DataDir, paste0(finput, ".bed")))) {
    return("finput must be a valid filename in DataDir.")
  }

  # Validate SNPvec
  if (length(SNPvec) == 0) {
    return("SNPvec must be a non-empty vector.")
  }

  # Validate extract
  if (!is.logical(extract) || length(extract) != 1) {
    return("extract must be a single boolean value.")
  }

  # Return NULL if all validations pass
  return(NULL)
}

## Function 149
## Added in 3.0
setFilterParameters <- function(CHRX, CHRY, filterCHR, regionfile, filterPAR, filterXTR, filterAmpliconic) {

  if (CHRX == FALSE && CHRY == FALSE){
    CHRX = TRUE
  }else{
    CHRX = CHRX
  }
  if (is.null(filterCHR)){
    fch <- NULL
  }else{
    fch <- filterCHR
    CHRX <- TRUE
    #CHR <-"chrX"
  }
  if (regionfile == FALSE){
    rf <- NULL
  }else{
    rf <- regionfile
    CHRX <- TRUE
    #CHR <-"chrX"
  }

  if (!is.null(fch) && !is.null(rf)){
    print("filterCHR and regionfile, cannot be in effect together.")
    return()
  }else if (!is.null(fch) && (filterPAR == TRUE| filterXTR == TRUE|filterAmpliconic == TRUE)){
    print("filterCHR cannot be in effect with other filters such as PAR, XTR, Ampliconic. Other filters were implicitly FALSE.")
  }else if (!is.null(rf) && (filterPAR == TRUE| filterXTR == TRUE|filterAmpliconic == TRUE)){
    print("regionfile cannot be in effect with other filters such as PAR, XTR, Ampliconic. Other filters were implicitly FALSE.")
  }

  return(list(CHRX = CHRX, fch = fch, rf = rf))
}

## Function 150
## Added in 3.0
processRegionFilter <- function(x, filterPAR, filterXTR, filterAmpliconic, ResultDir, DataDir, finput, foutput) {

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
        "--allow-no-sex",                      ## Adding in 4.0.
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
        "--allow-no-sex",                      ## Adding in 4.0.
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
        "--allow-no-sex",                             ## Adding in 4.0.
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


  return(list(par_snps, xtr_snps, ampliconic_snps, filterPAR, filterXTR, filterAmpliconic))
}

## Function 151
## Added in 3.0
readGenomicFeatures <- function(DataDir1, CHRX, CHRY, CHR, HG) {

  if (CHRX == TRUE && CHRY == TRUE) {
    x1 <- as.data.frame(
      read.table(
        file = paste0(DataDir1, "/", "X", "_genomic_features_", HG, ".bed.txt"),
        header = FALSE,
        sep = ""
      )
    )

    x2 <- as.data.frame(
      read.table(
        file = paste0(DataDir1, "/", "Y", "_genomic_features_", HG, ".bed.txt"),
        header = FALSE,
        sep = ""
      )
    )

    x <- rbind(x1, x2)
  } else {
    x <- as.data.frame(
      read.table(
        file = paste0(DataDir1, "/", CHR, "_genomic_features_", HG, ".bed.txt"),
        header = FALSE,
        sep = ""
      )
    )
  }
  return(x)
}

## Function 152
## Added in 3.0
filterGenomicFeatures <- function(x, filterPAR, filterXTR, filterAmpliconic) {

  if (filterPAR == TRUE && filterXTR == TRUE && filterAmpliconic == TRUE) {
    y <- x
  } else if (filterPAR == TRUE && filterXTR == TRUE && filterAmpliconic == FALSE) {
    y <- x[x$V4 == "PAR1" | x$V4 == "PAR2" | x$V4 == "XTR", ]
  } else if (filterPAR == TRUE && filterXTR == FALSE && filterAmpliconic == TRUE) {
    x1 <- x[grep("^Ampliconic", x$V4),]
    x2 <- x[x$V4 == "PAR1" | x$V4 == "PAR2", ]
    y <- rbind(x1, x2)
  } else if (filterPAR == FALSE && filterXTR == TRUE && filterAmpliconic == TRUE) {
    x1 <- x[grep("^Ampliconic", x$V4),]
    x2 <- x[x$V4 == "XTR", ]
    y <- rbind(x1, x2)
  } else if (filterPAR == TRUE && filterXTR == FALSE && filterAmpliconic == FALSE) {
    y <- x[x$V4 == "PAR1" | x$V4 == "PAR2", ]
  } else if (filterPAR == FALSE && filterXTR == TRUE && filterAmpliconic == FALSE) {
    y <- x[x$V4 == "XTR", ]
  } else if (filterPAR == FALSE && filterXTR == FALSE && filterAmpliconic == TRUE) {
    y <- x[grep("^Ampliconic", x$V4),]
  }
  return(y)
}

## Function 153
## Added in 3.0
executePlinkExcludeExtract <- function(ResultDir, DataDir, finput, rangefile, foutput) {

  executePlink <- function(DataDir, ResultDir, args) {

    tryCatch({
      stderr_dest <- ifelse(.Platform$OS.type == "windows", "NUL", "/dev/null")
      invisible(sys::exec_wait(file.path(ResultDir, "./plink"), args = args, std_err = stderr_dest))
    }, error = function(e) {
      stop("An error occurred while executing Plink: ", e$message)
    })
  }
  #globalVariables("DataDir","ResultDir")

  # Exclude region
  exclude_args <- c(
    "--bfile",
    paste0(DataDir, "/", finput),
    "--exclude",
    "range",
    rangefile,
    "--allow-no-sex",                ## Adding in 4.0
    "--make-bed",
    "--out",
    paste0(ResultDir, "/", foutput),
    "--silent"
  )

  executePlink(DataDir, ResultDir, exclude_args)

  # Extract region
  extract_args <- c(
    "--bfile",
    paste0(DataDir, "/", finput),
    "--extract",
    "range",
    rangefile,
    "--allow-no-sex",                                 ## Adding in 4.0
    "--make-bed",
    "--out",
    paste0(ResultDir, "/", foutput, "_snps_extracted"),
    "--silent"
  )

  executePlink(DataDir, ResultDir, extract_args)
}

## Function 154
## Added in 3.0
executePlinkChrFilter <- function(ResultDir, DataDir, finput, filterCHR, foutput) {

  executePlink <- function(ResultDir, DataDir, args) {
    # globalVariables("ResultDir")
    tryCatch({
      stderr_dest <- ifelse(.Platform$OS.type == "windows", "NUL", "/dev/null")
      invisible(sys::exec_wait(file.path(ResultDir, "./plink"), args = args, std_err = stderr_dest))
    }, error = function(e) {
      stop("An error occurred while executing Plink: ", e$message)
    })
  }

  # Filter SNPs not on specified chromosomes
  not_chr_args <- c(
    "--bfile",
    paste0(DataDir, "/", finput),
    "--not-chr",
    filterCHR,
    "--make-bed",
    "--out",
    paste0(ResultDir, "/", foutput),
    "--silent"
  )

  executePlink(ResultDir, DataDir, not_chr_args)

  # Filter SNPs on specified chromosomes
  chr_args <- c(
    "--bfile",
    paste0(DataDir, "/", finput),
    "--chr",
    filterCHR,
    "--make-bed",
    "--out",
    paste0(ResultDir, "/", foutput, "_snps_extracted"),
    "--silent"
  )

  executePlink(ResultDir, DataDir, chr_args)
}

## Function 155
## Added in 3.0
removeFiles <- function(fileNames, directory) {
  for (fileName in fileNames) {
    file.remove(file.path(directory, fileName))
  }
}




