# GXwasR helper functions
# @author Banabithi Bose
#' @importFrom stats na.omit median qchisq binomial dnorm glm na.exclude p.adjust pchisq pnorm lm pt qnorm setNames
#' @importFrom utils download.file read.table unzip write.table
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
  #Set working directory
  #xdir <- getwd()
  #setwd(wdir)
  #Specify operating system
  OS <- Sys.info()['sysname']
  if (OS == "Linux") {

    utils::download.file(destfile = paste0(wdir,"/","plink_linux_x86_64_20220402.zip"),
                         "https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20220402.zip", quiet = TRUE,
    )

    utils::unzip(paste0(wdir,"/","plink_linux_x86_64_20220402.zip"), exdir = wdir)
    Sys.chmod(paste0(wdir,"/plink"), mode = "0777", use_umask = TRUE)
    fileremove <- c("LICENSE","plink_linux_x86_64_20220402.zip","prettify","toy.map","toy.ped")
    f <- paste0(wdir,"/",fileremove)
    invisible(do.call(file.remove,list(f)))
    print("Program is continuing.")

  } else if (OS == "Windows") {

    utils::download.file(destfile = paste0(wdir,"/","plink_win64_20230116.zip"),
                         "https://s3.amazonaws.com/plink1-assets/plink_win64_20230116.zip", quiet = TRUE,
    )
    utils::unzip(paste0(wdir,"/","plink_win64_20230116.zip"), exdir = wdir)
    Sys.chmod(paste0(wdir,"/plink"), mode = "0777", use_umask = TRUE)
    fileremove <- c("LICENSE","plink_win64_20230116.zip","prettify.exe","toy.map","toy.ped")
    f <- paste0(wdir,"/",fileremove)
    invisible(do.call(file.remove,list(f)))
    print("Program starts.")

  }else if (OS == "macOS"){

    utils::download.file(destfile = paste0(wdir,"/","plink_mac_20230116.zip"),
                         "https://s3.amazonaws.com/plink1-assets/plink_mac_20230116.zip", quiet = TRUE,
    )
    utils::unzip(paste0(wdir,"/","plink_mac_20230116.zip"), exdir = wdir)
    Sys.chmod(paste0(wdir,"/plink"), mode = "0777", use_umask = TRUE)
    fileremove <- c("LICENSE","plink_mac_20230116.zip","prettify","toy.map","toy.ped")
    f <- paste0(wdir,"/",fileremove)
    invisible(do.call(file.remove,list(f)))
    print("Program starts.")

  } else {
    return(print("OS is not found by Sys.info()['sysname']."))
    print(
      "Plink excecutable needs to be in the 'wdir', Please check https://www.cog-genomics.org/plink2/"
    )
    print(
      "After unzipping, please type in R console: system('chmod 755 plink')' for OS = Linux and for OS = Windows"
    )
  }
}


###########
## this function is using datadir for plink...check
MFsplitPlink <- function(DataDir, ResultDir, finput, foutput, sex, xplink = FALSE, autoplink = FALSE){

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

  # cdir <- getwd()
  # setwd(DataDir)
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
  #return()
}


###########
## Function 2
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

#fl <- "FM01.assoc.logistic"
## Function 3
FMsub <- function(ResultDir, fl,plot.jpeg, plotname, snp_pval, annotateTopSnp, suggestiveline, genomewideline){

  if (file.exists(paste0(ResultDir,"/",fl))[1] == TRUE) {
    XWAS <-
      data.table::as.data.table(read.table(file = paste0(ResultDir,"/",fl),
                                           stringsAsFactors = FALSE,
                                           header = TRUE))
    gc(reset=TRUE)
    XWAS <- na.omit(XWAS)
    gc(reset=TRUE)
    XWAS_ADD <- na.omit(XWAS[XWAS$TEST == "ADD",c("SNP","CHR","BP","P")])
    gc(reset=TRUE)
    # From p-values, calculate chi-squared statistic
    chisq <- qchisq(1-XWAS_ADD$P,1)
    lamdaGC <- median(chisq)/qchisq(0.5,1)

    XWAS_ADD_X <- na.omit(XWAS_ADD[XWAS_ADD$CHR == 23,])
    chisq1 <- qchisq(1-XWAS_ADD_X$P,1)
    lamdaGC1 <- median(chisq1)/qchisq(0.5,1)

    if (plot.jpeg[1] == TRUE){
      options(bitmapType='cairo')
      grDevices::jpeg(paste0(ResultDir,"/",plotname,".jpeg"),  width = 20,
                      height = 10,
                      units = 'in',
                      res = 300)
      graphics::par(mfrow = c(2, 2))

      suppressWarnings(qqman::manhattan(XWAS_ADD, ylim = c(0,16),suggestiveline = suggestiveline, genomewideline = genomewideline, annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of GWAS"))
      suppressWarnings(qqman::qq(XWAS_ADD$P, main = paste0(("Q-Q plot of GWAS p-values with GIF = "), round(lamdaGC,3))))

      if (nrow(XWAS_ADD_X)!=0){
        suppressWarnings(qqman::manhattan(XWAS_ADD_X, ylim = c(0,16),suggestiveline = suggestiveline, genomewideline = genomewideline, annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of XWAS"))
        suppressWarnings(qqman::qq(XWAS_ADD_X$P, main = paste0(("Q-Q plot of XWAS p-values with GIF = "), round(lamdaGC1,3))))

      }else{
        print("X-chromosome may not be present.")
      }

      gc(reset=TRUE)
      dev.off()
    }else if (plot.jpeg[1] == FALSE){
      par(mar = c(1, 1, 1, 1))
      # graphics::par(mfrow = c(2, 2))
      suppressWarnings(qqman::manhattan(XWAS_ADD, ylim = c(0,16), suggestiveline = suggestiveline, genomewideline = genomewideline, annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of GWAS"))
      suppressWarnings(qqman::qq(XWAS_ADD$P, main = paste0(("Q-Q plot of GWAS p-values with GIF = "), round(lamdaGC,3))))

      if (nrow(XWAS_ADD_X)!=0){

        suppressWarnings(qqman::manhattan(XWAS_ADD_X, ylim = c(0,16), suggestiveline = suggestiveline, genomewideline = genomewideline, annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of XWAS"))
        suppressWarnings(qqman::qq(XWAS_ADD_X$P, main = paste0(("Q-Q plot of XWAS p-values with GIF = "), round(lamdaGC1,3))))
        # graphics::par(mfrow = c(1, 1))

      }else{
        print("X-chromosome is may not present.")
      }
      gc(reset=TRUE)
    }

    return(na.omit(XWAS))
  } else if (file.exists(paste0(ResultDir,"/",fl))[1] == FALSE) {
    print(paste0("XWAS cannot be performed. Check the ", stringr::str_sub(fl, 1,5),"log file in ResultDir for checking the error."))
  }
}

## Function 4
## Only DataDir
FMmain <- function(DataDir, ResultDir, finput, trait, standard_beta, perm, mperm, genedrop, xmodel,
                   noxsex, covarfile, interaction, covartest, Inphenocov, plot.jpeg, plotname, snp_pval, annotateTopSnp, suggestiveline, genomewideline, ncores){


  if (xmodel[1] == "FM01"){
    modelv <- 1
  }else if (xmodel[1] == "FM02"){
    modelv <- 2
  }

  if (xmodel[1] == "FM01comb"){
    modelv <- 1
  }else if (xmodel[1] == "FM02comb"){
    modelv <- 2}

  if (trait[1] == "binary"){
    regress <- "--logistic"
  }else if (trait[2] == "quantitative"){
    regress <- "--linear"}

  if (trait[1] == "quantitative" && standard_beta[1] == TRUE){
    standard_b <- "--standard-beta"
  }else{
    standard_b <- "beta"
  }

  if (mperm == 0){
    mpermv = NULL
    mperm = NULL
  }else{
    mpermv = "--mperm"
    mperm = mperm
  }

  # if (perm == FALSE | mperm != 0){
  #   permv = NULL
  # }else{
  #   permv = "perm"
  # }
  if (perm[1] == FALSE){
    permv = NULL
  }else{
    permv = "perm"
  }

  # if (mperm == 0){
  #   permv = NULL
  # }else{
  #   permv = "perm"
  # }

  if (noxsex[1] == FALSE){
    noxsexv = NULL
  }else{
    noxsexv = "no-x-sex"
  }

  if (genedrop[1] == FALSE){
    genedropv = NULL
  }else{
    genedropv = "genedrop"
  }

  if (is.null(covarfile)){
    covar <- NULL
    covarv <- NULL
    ##NEWLY
    covartest <- NULL
  }else{

    if (covartest == "ALL") {
      ## Copy covarfile to ResultDir
      #ftemp <- list.files(paste0(ResultDir,"/"),pattern = "PostimputeEX_QC2")
      invisible(file.copy(paste0(DataDir,"/",list(covarfile)),ResultDir))

      covarfile <- paste0(ResultDir,"/",covarfile)

    } else if (covartest != "ALL" & !is.null(covarfile)) {
      Run_newcovarfile(DataDir = DataDir, ResultDir=ResultDir,covarfile = covarfile,covartest = covartest)
      covarfile <- paste0(ResultDir,"/newcovarfile.txt")
    }
    covar <- "--covar"
    #covarv <- covarfile
    #covarv <- paste0(DataDir,"/",covarfile)
    covarv <- covarfile
  }

  #print("line 180")


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
        standard_b,
        genedropv,
        noxsexv,
        permv,
        "intercept",
        interactionv,
        parameterv,
        Inphenocovv,
        mpermv,
        mperm,
        covar,
        covarv,
        "--out",
        paste0(ResultDir,"/",xmodel),
        "--silent"
      ),
      std_out = FALSE,
      std_err = FALSE
    ))
  }else{

    snpfile = read.table(paste0(DataDir,"/",finput,".bim"))
    chunk <- round(nrow(snpfile)/ncores)+1
    chunks <- round(seq(1, length(snpfile), by = chunk),0)
    #chunks <- seq(1, nrow(snpfile), by = chunk)

    # for (i in 1:ncores){
    #   j <- i*chunk
    #   k <-
    #   print(i)
    #   snp_names <- snpfile$V2[i:j]
    #   write.table(snp_names, file = paste0(ResultDir,"/chunk_",i), quote = FALSE, row.names = FALSE, col.names = FALSE)
    #
    # }
    #
    # chunks <- c(chunk,chunk+1)

    ## Making chunkwise snpfiles
    #snpfiles = 1:nrow(snpfile)
    #snpfiles = 1:10
    paraGwas <- function(chunks,chunk,ResultDir,DataDir,finput,modelv,regress,standard_b,genedropv,noxsexv,
                         permv,interactionv,parameterv,Inphenocovv,mpermv,mperm,covar,covarv,
                         snpfile){

      ## Chunkfile create
      if (nrow(snpfile)>=(chunks+chunk)){
        snp_names <- snpfile$V2[chunks:(chunks+chunk)]
      }else{
        snp_names <- snpfile$V2[chunks:nrow(snpfile)]
      }
      #snp_names <- snpfile$V2
      ## Create snps with chunk
      write.table(snp_names,file = paste0(ResultDir,"/",chunks,"_snps"), quote = FALSE, row.names = FALSE, col.names = FALSE)
      #snp_names <- snpfile$V2[chunks:chunks+chunk]
      #snp_name <- snpfile$V2[snpfiles]
      # make plink bedfiles for this SNP in DataDir
      invisible(sys::exec_wait(
        paste0(ResultDir,"/","./plink"),
        args = c(
          "--bfile",
          paste0(DataDir,"/",finput),
          "--extract",
          paste0(ResultDir,"/",chunks,"_snps"),
          #"--snps",
          #snp_names,
          "--xchr-model",
          modelv,
          "--freq",
          "--ci", 0.95,
          regress,
          standard_b,
          genedropv,
          noxsexv,
          permv,
          "intercept",
          interactionv,
          parameterv,
          Inphenocovv,
          mpermv,
          mperm,
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
      # ftemp <- list.files(paste0(ResultDir,"/"),pattern = paste0(chunks,"_ss"))
      # invisible(file.remove(paste0(ResultDir,"/",ftemp)))
      # ftemp <- list.files(paste0(ResultDir,"/"),pattern = paste0(chunks,"_snps"))
      # invisible(file.remove(paste0(ResultDir,"/",ftemp)))
      # single_snp_result <- unique(single_snp_result)
      return(single_snp_result)
    }

    #### Parallel computation
    print("cluster making started")
    cl <- parallel::makeCluster(mc <- getOption("cl.cores", ncores),type = "FORK")

    invisible(parallel::clusterEvalQ(cl, library(data.table)))
    invisible(parallel::clusterEvalQ(cl, library(parallel)))
    print("cluster export started")
    ## Took time
    # parallel::clusterExport(cl=cl,c("ResultDir","DataDir","finput","modelv","regress","standard_b","genedropv","noxsexv",
    #                         "permv","interactionv","parameterv","Inphenocovv","mpermv","mperm","covar","covarv",
    #                         "snpfile"),envir=environment())
    parallel::clusterExport(cl=cl,NULL,envir=environment())


    print("cluster making done")
    all_snps_results <- data.table::as.data.table(data.table::rbindlist(parallel::parLapply(cl,chunks,paraGwas,chunk=chunk,ResultDir=ResultDir,
                                                                                            DataDir=DataDir,finput=finput,modelv=modelv,regress=regress,
                                                                                            standard_b=standard_b,genedropv=genedropv,noxsexv=noxsexv,
                                                                                            permv=permv,interactionv=interactionv,parameterv=parameterv,
                                                                                            Inphenocovv=Inphenocovv,mpermv=mpermv,mperm=mperm,covar=covar,covarv=covarv,snpfile=snpfile)))
    all_snps_results <- unique(all_snps_results)
    #all_snps_results <- data.table::as.data.table(data.table::rbindlist(parallel::parLapplyLB(cl,1:2,paraGwas)))

    #save(all_snps_results, file=paste0(ResultDir,"/all_snps_results.rda"))
    parallel::stopCluster(cl)
    print("clusters stopped.")
    ####
    write.table(all_snps_results, file= paste0(ResultDir,"/",xmodel,".assoc.logistic"), quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\r\n", sep = " ")
    rm(all_snps_results)
    # ftemp <- list.files(paste0(ResultDir,"/"),pattern = paste0(chunks,"_ss"))
    # invisible(file.remove(paste0(ResultDir,"/",ftemp)))
    # ftemp <- list.files(paste0(ResultDir,"/"),pattern = paste0(chunks,"_snps"))
    # invisible(file.remove(paste0(ResultDir,"/",ftemp)))
    # single_snp_result <- unique(single_snp_result)
  }
  fl <- paste0(xmodel,".assoc.logistic")
  x <- FMsub(ResultDir = ResultDir, fl,plot.jpeg = plot.jpeg, plotname = plotname, snp_pval =snp_pval, annotateTopSnp = annotateTopSnp, suggestiveline = suggestiveline, genomewideline = genomewideline)
  gc(reset = TRUE)
  ## Remove the files from DataDir
  #ftemp <- list.files(paste0(ResultDir,"/"),pattern = paste0(xmodel))
  #invisible(file.remove(paste0(ResultDir,"/",ftemp)))
  #ftemp <- list.files(paste0(ResultDir,"/"),pattern = "finput")
  #invisible(file.remove(paste0(DataDir,"/",ftemp)))
  gc(reset = TRUE)
  return(x)
}

## Function 5
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

## Function 6
# Copied from ex-CRAN package MADAM and exported. The man pages are copied from
# the original package.
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


## Function 8
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

## Function 9
# FMcomb_sub <- function(ResultDir, FemaleWAS, MaleWAS, combtest, MF.p.corr,
#                        MF.zero.sub,
#                        MF.na.rm,
#                        MF.mc.cores,
#                        MF.comb.na.rm,
#                        MF.comb.mc.cores,
#                        B,
#                        plot.jpeg,
#                        plotname,
#                        snp_pval,
#                        annotateTopSnp,
#                        suggestiveline,
#                        genomewideline)
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
                       genomewideline){

  #load(paste0("~/",ResultDir,"/MaleWAS.Rda"))
  load(paste0(ResultDir,"/MaleWAS.Rda"))
  gc(reset=TRUE)
  MaleWAS <- data.table::as.data.table(MaleWAS)
  gc(reset=TRUE)
  #load(paste0("~/",ResultDir,"/FemaleWAS.Rda"))
  load(paste0(ResultDir,"/FemaleWAS.Rda"))
  FemaleWAS <- data.table::as.data.table(FemaleWAS)
  gc(reset=TRUE)
  MFWAS <- merge(FemaleWAS, MaleWAS, by = c("SNP", "A1", "TEST"))
  #MFWAS <- merge(FemaleWAS, MaleWAS, by = c("SNP", "TEST"))
  gc(reset=TRUE)
  MFWAS1 <- MFWAS[MFWAS$TEST == "ADD",]
  pvals <- as.data.frame(MFWAS1[, c(12, 21)])
  #rm(MaleWAS)
  #rm(FemaleWAS)
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
      print("cluster making started")
      cl <- parallel::makeCluster(mc <- getOption("cl.cores", ncores),type = "FORK")

      invisible(parallel::clusterEvalQ(cl, library(data.table)))
      invisible(parallel::clusterEvalQ(cl, library(parallel)))
      print("cluster export started")
      parallel::clusterExport(cl=cl,NULL,envir=environment())


      print("cluster making done")

      Pnew <- data.table::as.data.table(data.table::rbindlist(parallel::parLapply(cl,chunks,paraStouffer,chunk=chunk,pvals=pvals,MF.p.corr=MF.p.corr,MF.zero.sub=MF.zero.sub,MF.na.rm=MF.na.rm)))

      parallel::stopCluster(cl)
      print("clusters stopped.")
    }
    #Result <- cbind(MFWAS[, 1:7], Pnew[, 3:4])
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
  chisq <- qchisq(1-Result$P,1)
  lamdaGC <- median(chisq)/qchisq(0.5,1)
  chisq1 <- qchisq(1-XWAS_ADD_X$P,1)
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


  if (plot.jpeg[1] == TRUE){
    options(bitmapType='cairo')
    grDevices::jpeg(paste0(ResultDir,"/",plotname,".jpeg"),  width = 20,
                    height = 10,
                    units = 'in',
                    res = 300)
    graphics::par(mfrow = c(2, 2))

    invisible(gmirror(top=FemaleWAS, bottom=MaleWAS, tline= snp_pval, bline=snp_pval,
                              toptitle="GWAS of females", bottomtitle = "GWAS of males",
                              highlight_p = c(snp_pval,snp_pval), highlighter="green", chrblocks = TRUE,file = paste0(ResultDir,"/","Stratified_GWAS")))

    print(paste0("Miami plot of stratified GWAS is saved in ", ResultDir))
    gc(reset=TRUE)

    if (nrow(gwas.t2)!=0 && nrow(gwas.b2)!=0){
      invisible(gmirror(top=gwas.t2, bottom=gwas.b2, tline=snp_pval, bline=snp_pval,
                                toptitle="XWAS of females", bottomtitle = "XWAS of males",
                                highlight_p = c(snp_pval,snp_pval), highlighter="green", chrblocks = TRUE, file = paste0(ResultDir,"/","Stratified_XWAS")))

      print(paste0("Miami plot of stratified XWAS is saved in ", ResultDir))

    }else{
      print("Miami plot for stratified XWAS cannot be drawn.")
    }
    gc(reset=TRUE)


    suppressWarnings(qqman::manhattan(Result, ylim = c(0,16), suggestiveline = suggestiveline, genomewideline = genomewideline ,annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of male-female combined GWAS"))
    gc(reset=TRUE)

    if (nrow(XWAS_ADD_X)!=0){
      suppressWarnings(qqman::manhattan(XWAS_ADD_X, ylim = c(0,16), suggestiveline = suggestiveline, genomewideline = genomewideline, annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of male-female combined XWAS"))
    }else{
      print("There may not be any X chromosome in the data.")
    }
    gc(reset=TRUE)
    if (sum(Result$P) == nrow(Result)){
      print("All adjusted p values are 1, qq plot cannot be created.")
    }else{
      suppressWarnings(qqman::qq(Result$P, main = paste0(("Q-Q plot of GWAS male-female combined with GIF = "), round(lamdaGC,3))))
      gc(reset=TRUE)
    }

    if (sum(XWAS_ADD_X$P) == nrow(XWAS_ADD_X)){
      print("All adjusted p values are 1, qq plot cannot be created.")
    }else{
      suppressWarnings(qqman::qq(XWAS_ADD_X$P, main = paste0(("Q-Q plot of XWAS male-female combined p-values with GIF = "), lamdaGC1)))
    }

    dev.off()
  }else{
    #graphics::par(mfrow = c(2, 2))
    gc(reset=TRUE)
    print(suppressWarnings(gmirror(top=FemaleWAS, bottom=MaleWAS, tline= snp_pval, bline=snp_pval,
                                           toptitle="GWAS of females", bottomtitle = "GWAS of males",
                                           highlight_p = c(snp_pval,snp_pval), highlighter="green", chrblocks = TRUE, file = paste0(ResultDir,"/","Stratified_GWAS"))))


    gc(reset=TRUE)
    print(suppressWarnings(gmirror(top=gwas.t2, bottom=gwas.b2, tline=snp_pval, bline=snp_pval,
                                           toptitle="XWAS of females", bottomtitle = "XWAS of males",
                                           highlight_p = c(snp_pval,snp_pval), highlighter="green", chrblocks = TRUE, file = paste0(ResultDir,"/","Stratified_XWAS"))))
    gc(reset=TRUE)
    par(mar = c(1, 1, 1, 1))
    suppressWarnings(qqman::manhattan(Result, ylim = c(0,16), suggestiveline = suggestiveline, genomewideline = genomewideline ,annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of male-female combined GWAS"))
    gc(reset=TRUE)
    suppressWarnings(qqman::manhattan(XWAS_ADD_X, ylim = c(0,16), suggestiveline = suggestiveline, genomewideline = genomewideline, annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of male-female combined XWAS"))
    gc(reset=TRUE)
    if (sum(Result$P) == nrow(Result)){
      print("All adjusted p values are 1, qq plot cannot be created.")
    }else{
      gc(reset=TRUE)
      suppressWarnings(qqman::qq(Result$P, main = paste0(("Q-Q plot of male-female combined GWAS with GIF = "), round(lamdaGC,3))))
      gc(reset=TRUE)
    }

    if (sum(XWAS_ADD_X$P) == nrow(XWAS_ADD_X)){
      print("All adjusted p values are 1, qq plot cannot be created.")
    }else{
      gc(reset=TRUE)
      qqman::qq(XWAS_ADD_X$P, main = paste0(("Q-Q plot of male-female combined XWAS with GIF = "), round(lamdaGC1,3)))
      gc(reset=TRUE)
    }
    #graphics::par(mfrow = c(1, 1))

  }
  gc(reset=TRUE)
  return(na.omit(Result))
}

## Function 10
FMcomb <-
  function(DataDir,ResultDir, perm,
           trait,
           mperm,
           standard_beta,
           genedrop,
           noxsex,
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
    invisible(file.copy(paste0(DataDir,"/",list(covarfile)),ResultDir))

    MaleWAS <- FMmain(DataDir = ResultDir, ResultDir = ResultDir, finput = "finput.male", trait = trait, standard_beta = standard_beta, perm = perm, mperm = mperm, genedrop = genedrop, xmodel = xmodel,
                      noxsex = noxsex, covarfile = covarfile, interaction = interaction, covartest = covartest, Inphenocov = Inphenocov, plot.jpeg = TRUE, plotname = paste0(xmodel,"_MaleWAS"), snp_pval = snp_pval, annotateTopSnp = annotateTopSnp, suggestiveline = suggestiveline, genomewideline = genomewideline, ncores = ncores)
    gc(reset=TRUE)
    MaleWAS <- na.omit(MaleWAS)
    save(MaleWAS, file = paste0(ResultDir,"/MaleWAS.Rda"))
    rm(MaleWAS)
    gc(reset=TRUE)
    FemaleWAS <- FMmain(DataDir = ResultDir, ResultDir = ResultDir, finput = "finput.female", trait = trait, standard_beta = standard_beta, perm = perm, mperm = mperm, genedrop = genedrop, xmodel = xmodel,
                        noxsex = noxsex, covarfile = covarfile, interaction = interaction, covartest = covartest, Inphenocov = Inphenocov, plot.jpeg = TRUE, plotname = paste0(xmodel,"_FemaleWAS"), snp_pval = snp_pval, annotateTopSnp = annotateTopSnp, suggestiveline = suggestiveline, genomewideline = genomewideline, ncores = ncores)
    gc(reset=TRUE)
    FemaleWAS <- na.omit(FemaleWAS)
    save(FemaleWAS, file = paste0(ResultDir,"/FemaleWAS.Rda"))
    rm(FemaleWAS)
    ## Remove the files from DataDir
    #ftemp <- list.files(paste0(ResultDir,"/"),pattern = paste0(xmodel))
    #invisible(file.remove(paste0(ResultDir,"/",ftemp)))
    #ftemp <- list.files(paste0(ResultDir,"/"),pattern = "finput")
    #invisible(file.remove(paste0(DataDir,"/",ftemp)))

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
                              genomewideline = genomewideline)
    gc(reset=TRUE)
    #Result <- na.omit(Result)
    CombinedWAS <- na.omit(CombinedWAS)
    save(CombinedWAS,file = paste0(ResultDir,"/","CombinedWAS.Rda"))
    rm(CombinedWAS)
    gc(reset=TRUE)
    print(paste0("Three dataframes such as, CombinedWAS, MaleWAS and FemaleWAS are produced in", ResultDir))
    return()
    #return(list(CombinedWAS=Result,MaleWAS=MaleWAS,FemaleWAS=FemaleWAS))
  }


#####################
#####################

###XCMAX4

# 1. Function: XCMAX4.
# This function is used to test the association between an X chromosomal marker and a binary trait. One SNP at a time. # Add sex as covariate ## Can't be for this function
#
# 2. Arguments of the function
# data: The data should contain the information of phenotype, genotype, sex, and possible additional covariates . Each row represents the data of a specific individual.
# The first column records the phenotype information, with 1 representing the case and 0 representing the control.
# The second column provides the genotype information, with 0, 1 and 2 representing the number of risk alleles.
# The third column contains the sex information, with 0 being male and 1 being female.
# The remaining columns records the information of possible additional covariates.
#data <- read.table("/projects/b1137/BBose/ProjectGWAS/XWAS_QC/Panscan/DG_98846_1_H55_c1_panscan_v3/han.csv",sep = ",",header = TRUE)

## Function 11
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

## Function 12
## without covariate file
XCMAX4_data1 <- function(Snp,genosnp,P,Samp){

  #g <- data.table::as.data.table(SNPRelate::snpgdsGetGeno(genofile, snp.id = Snp,verbose = FALSE))
  genosnp <- as.data.frame(genosnp)
  g <- genosnp[,Snp,drop = FALSE]
  colnames(g) <- c("Snp")
  g$IID <- Samp
  Pg <- merge(g,P,by = "IID")
  Data <- as.data.frame(Pg[,c(4,2,3)])
  Data$phenotype <- as.numeric(as.character(Data$phenotype))
  colnames(Data) <- c("D","genotype","gender")
  Data$gender <- as.numeric(as.character(Data$gender))

  # Add sex as covariate ## Can't be for this function
  #Data$Cov1 <- Data$gender
  #Data$Cov1 <- as.numeric(as.character(Data$Cov1))
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
  invisible(Z <- XCMAX4(Data))
  pval <- Z$`p-value`
  statistic <- Z$statictic
  SE <- Z$`standard-error`
  result <- as.data.frame(cbind(Snp,statistic,pval,SE))
  colnames(result) <- c("SNP","STAT","P","SE")
  return(result)
}

## Function 13
## with covariate file
XCMAX4_data2 <- function(Snp,DataDir,genosnp,P,covarfile,Samp){
  #XCMAX4_data2 <- function(Snp,DataDir,ResultDir,P,covarfile){

  #gdsfmt::showfile.gds(closeall=TRUE)
  #genofile <- SNPRelate::snpgdsOpen(f)
  # bedfile<- paste0(ResultDir,"/","XPLINK.bed")
  # famfile<- paste0(ResultDir,"/","XPLINK.fam")
  # bimfile<- paste0(ResultDir,"/","XPLINK.bim")
  #
  # ### Convert PLINK FILES INTO GDS FORMAT FOR R
  #
  # f<-SNPRelate::snpgdsBED2GDS(bedfile, famfile, bimfile, paste0(ResultDir,"/","IPCgeno.gdsX"))
  # print("f done")
  # # Open the GDS file and
  #genofile <- SNPRelate::snpgdsOpen(f)
  covarfile1 <- read.table(file = paste0(DataDir,"/",covarfile),
                           stringsAsFactors = FALSE,
                           header = TRUE)
  #g <- data.table::as.data.table(SNPRelate::snpgdsGetGeno(genofile, snp.id = Snp, verbose = FALSE))
  genosnp <- as.data.frame(genosnp)
  g <- genosnp[,Snp,drop = FALSE]
  colnames(g) <- c("Snp")
  #Samp <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id"))
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

  Z <- XCMAX4(Datatest)
  pval <- Z$`p-value`
  statistic <- Z$statictic
  SE <- Z$`standard-error`
  result <- as.data.frame(cbind(Snp,statistic,pval,SE))
  colnames(result) <- c("SNP","STAT","P","SE")
  return(result)
}

## Function 14
## This function as using DataDir for plink
autoFun <- function(DataDir, ResultDir, finput, trait, standard_beta, perm, mperm, genedrop,
                    noxsex, covarfile, interaction, covartest, Inphenocov,xmodel, ncores){

  #cat(paste0("Main arguments are: \n", "finput = ",finput,"\n trait = ",trait,"\n standard_beta = ",standard_beta,"\n perm = ",perm,"\n mperm = ",mperm,"\n genedrop = ",genedrop,"\n xmodel = ",xmodel,"\n noxsex = ",noxsex,"\n covarfile = ",covarfile,"\n interaction = ",interaction,"\n"))

  modelv <- 1
  regress <- "--logistic"
  standard_b <- "beta"

  if (mperm == 0){
    mpermv = NULL
    mperm = NULL
  }else{
    mpermv = "--mperm"
    mperm = mperm
  }

  if (perm[1] == FALSE){
    permv = NULL
  }else{
    permv = "perm"
  }

  # if (mperm != 0){
  #   permv = NULL
  # }else{
  #   permv = "perm"
  # }
  if (noxsex[1] == FALSE){
    noxsexv = NULL
  }else{
    noxsexv = "no-x-sex"
  }

  if (genedrop[1] == FALSE){
    genedropv = NULL
  }else{
    genedropv = "genedrop"
  }

  if (is.null(covarfile)){
    covar <- NULL
    covarv <- NULL
  }else{

    ########
    if (covartest == "ALL") {
      ## Copy covarfile to ResultDir
      #ftemp <- list.files(paste0(ResultDir,"/"),pattern = "PostimputeEX_QC2")
      invisible(file.copy(paste0(DataDir,"/",list(covarfile)),ResultDir))

      covarfile <- paste0(ResultDir,"/",covarfile)

    } else if (covartest != "ALL" & !is.null(covarfile)) {
      Run_newcovarfile(DataDir = DataDir, ResultDir=ResultDir,covarfile = covarfile,covartest = covartest)
      covarfile <- paste0(ResultDir,"/newcovarfile.txt")
    }
    covar <- "--covar"
    covarv <- covarfile

  }

  print("line 180")


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
        "--xchr-model",
        modelv,
        "--freq",
        "--ci", 0.95,
        regress,
        standard_b,
        genedropv,
        noxsexv,
        permv,
        "intercept",
        interactionv,
        parameterv,
        Inphenocovv,
        mpermv,
        mperm,
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

    paraGwas <- function(chunks,chunk,ResultDir,finput,modelv,regress,standard_b,genedropv,noxsexv,
                         permv,interactionv,parameterv,Inphenocovv,mpermv,mperm,covar,covarv,
                         snpfile){

      ## Chunkfile create
      if (nrow(snpfile)>=(chunks+chunk)){
        snp_names <- snpfile$V2[chunks:(chunks+chunk)]
      }else{
        snp_names <- snpfile$V2[chunks:nrow(snpfile)]
      }
      #snp_names <- snpfile$V2
      ## Create snps with chunk
      write.table(snp_names,file = paste0(ResultDir,"/",chunks,"_snps"), quote = FALSE, row.names = FALSE, col.names = FALSE)
      #snp_names <- snpfile$V2[chunks:chunks+chunk]
      #snp_name <- snpfile$V2[snpfiles]
      # make plink bedfiles for this SNP in DataDir
      invisible(sys::exec_wait(
        paste0(ResultDir,"/","./plink"),
        args = c(
          "--bfile",
          paste0(ResultDir,"/",finput),
          "--extract",
          paste0(ResultDir,"/",chunks,"_snps"),
          "--xchr-model",
          modelv,
          "--freq",
          "--ci", 0.95,
          regress,
          standard_b,
          genedropv,
          noxsexv,
          permv,
          "intercept",
          interactionv,
          parameterv,
          Inphenocovv,
          mpermv,
          mperm,
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
      # ftemp <- list.files(paste0(ResultDir,"/"),pattern = paste0(chunks,"_ss"))
      # invisible(file.remove(paste0(ResultDir,"/",ftemp)))
      # ftemp <- list.files(paste0(ResultDir,"/"),pattern = paste0(chunks,"_snps"))
      # invisible(file.remove(paste0(ResultDir,"/",ftemp)))
      single_snp_result <- unique(single_snp_result)
      return(single_snp_result)
    }

    #### Parallel computation
    print("cluster making started")
    cl <- parallel::makeCluster(mc <- getOption("cl.cores", ncores),type = "FORK")

    invisible(parallel::clusterEvalQ(cl, library(data.table)))
    invisible(parallel::clusterEvalQ(cl, library(parallel)))
    print("cluster export started")
    parallel::clusterExport(cl=cl,NULL,envir=environment())


    print("cluster making done")
    all_snps_results <- data.table::as.data.table(data.table::rbindlist(parallel::parLapply(cl,chunks,paraGwas,chunk=chunk,ResultDir=ResultDir,finput=finput,modelv=modelv,regress=regress,
                                                                                            standard_b=standard_b,genedropv=genedropv,noxsexv=noxsexv,
                                                                                            permv=permv,interactionv=interactionv,parameterv=parameterv,
                                                                                            Inphenocovv=Inphenocovv,mpermv=mpermv,mperm=mperm,covar=covar,covarv=covarv,snpfile=snpfile)))
    all_snps_results <- unique(all_snps_results)
    #all_snps_results <- data.table::as.data.table(data.table::rbindlist(parallel::parLapplyLB(cl,1:2,paraGwas)))

    #save(all_snps_results, file=paste0(ResultDir,"/all_snps_results.rda"))
    parallel::stopCluster(cl)
    print("clusters stopped.")
    ####
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


## Function 15
#This function is also using DataDir for plink
XCMAFun <- function(DataDir, ResultDir, finput,trait, standard_beta, perm, mperm, genedrop,noxsex,
                    covarfile,covartest,interaction,Inphenocov,plot.jpeg, plotname, snp_pval,
                    annotateTopSnp, suggestiveline = suggestiveline, genomewideline = genomewideline,
                    ncores = ncores){

  ## Making separate plink files foor x chr and autosome and preparing input for XCMAX4().
  gdsfmt::showfile.gds(closeall=TRUE)
  ## Making plinks x chr and other chr separately

  invisible(sys::exec_wait(
    paste0(ResultDir,"/","./plink"),
    args = c(
      "--bfile",
      paste0(DataDir,"/",finput),
      "--chr",
      23,
      "--make-bed",
      "--out",
      paste0(ResultDir,"/","XPLINK"),
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
      "--not-chr",
      23,
      "--make-bed",
      "--out",
      paste0(ResultDir,"/","FilteredX"),
      "--silent"
    ),
    std_out = FALSE,
    std_err = FALSE
  ))


  bedfile<- paste0(ResultDir,"/","XPLINK.bed")
  famfile<- paste0(ResultDir,"/","XPLINK.fam")
  bimfile<- paste0(ResultDir,"/","XPLINK.bim")

  ### Convert PLINK FILES INTO GDS FORMAT FOR R

  f<-SNPRelate::snpgdsBED2GDS(bedfile, famfile, bimfile, paste0(ResultDir,"/","IPCgeno.gdsX"))
  print("f done")
  # Open the GDS file and
  genofile <- SNPRelate::snpgdsOpen(f)
  print("g done")
  G <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "genotype"))
  print("G done")
  Samp <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id"))
  print("Samp done")
  ## can test with few SNPs
  Snp <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.id"))
  print("Snp done")
  #print(Snp)
  #gdsfmt::showfile.gds(closeall=TRUE)

  P <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.annot"))
  P[P$sex == "M",1] <- 0
  P[P$sex == "F",1] <- 1
  P$sex <- as.numeric(as.character(P$sex))
  P["phenotype"][P["phenotype"] == 1] <- 0
  P["phenotype"][P["phenotype"] == 2] <- 1
  P$IID <- Samp
  ########
  genosnp <- data.table::as.data.table(SNPRelate::snpgdsGetGeno(genofile, snp.id = Snp,verbose = FALSE))
  colnames(genosnp)<- Snp

  #ncores = 0
  ########


  if (is.null(covarfile)) {

    if (ncores == 0){

      Result_pval <- data.table::rbindlist(lapply(Snp,XCMAX4_data1,genosnp=genosnp,P=P,Samp=Samp))

    }else{
      #snpfile = read.table(paste0(ResultDir,"/",finput,".bim"))
      chunk <- (round(length(Snp)/ncores)+1)/3
      chunks <- round(seq(1, length(Snp), by = chunk),0)

      xcmaPara <- function(chunks,chunk,Snp,genosnp,P,Samp){
        ## Chunkfile create
        if (length(Snp)>=(chunks+chunk)){
          snp_names <- Snp[chunks:(chunks+chunk)]
        }else{
          snp_names <- Snp[chunks:length(Snp)]
        }
        ## without covariate file
        XCMAX4_data1 <- function(Snp,genosnp,P,Samp){

          #g <- data.table::as.data.table(SNPRelate::snpgdsGetGeno(genofile, snp.id = Snp,verbose = FALSE))
          genosnp <- as.data.frame(genosnp)
          g <- genosnp[,Snp,drop = FALSE]
          colnames(g) <- c("Snp")
          g$IID <- Samp
          Pg <- merge(g,P,by = "IID")
          Data <- as.data.frame(Pg[,c(4,2,3)])
          Data$phenotype <- as.numeric(as.character(Data$phenotype))
          colnames(Data) <- c("D","genotype","gender")
          Data$gender <- as.numeric(as.character(Data$gender))

          # Add sex as covariate ## Can't be for this function
          #Data$Cov1 <- Data$gender
          #Data$Cov1 <- as.numeric(as.character(Data$Cov1))
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
          invisible(Z <- XCMAX4(Data))
          pval <- Z$`p-value`
          statistic <- Z$statictic
          SE <- Z$`standard-error`
          result <- as.data.frame(cbind(Snp,statistic,pval,SE))
          colnames(result) <- c("SNP","STAT","P","SE")
          return(result)
        }
        Result_pval1 <- data.table::rbindlist(lapply(snp_names,XCMAX4_data1,genosnp=genosnp,P=P,Samp=Samp))
        gc(reset = TRUE)
        return(Result_pval1)
      }
      #### Parallel computation
      gc(reset = TRUE)
      print("cluster making started")
      cl <- parallel::makeCluster(mc <- getOption("cl.cores", ncores),type = "FORK")

      invisible(parallel::clusterEvalQ(cl, library(data.table)))
      invisible(parallel::clusterEvalQ(cl, library(parallel)))
      print("cluster export started")
      parallel::clusterExport(cl=cl,NULL,envir=environment())
      print("cluster making done")
      Result_pval <- data.table::rbindlist(parallel::parLapply(cl, chunks,xcmaPara,chunk=chunk,Snp=Snp,genosnp=genosnp,P=P,Samp=Samp))
      parallel::stopCluster(cl)
      print("clusters stopped.")
      gc(reset = TRUE)
    }

    Result_pval$CHR <- "23"
    Result_pval$TEST <- "XCGA"
    Result_pval$BETA <- "NA"
    Result_pval$L95 <- "NA"
    Result_pval$U95 <- "NA"

  } else{

    DataDir=DataDir
    covarfile=covarfile
    print("Running XCGA model for X chromosome.")

    if (ncores == 0){

      Result_pval <- data.table::rbindlist(lapply(Snp,XCMAX4_data2,DataDir=DataDir,genosnp=genosnp,P=P,Samp = Samp,covarfile=covarfile))

    }else{

      chunk <- (round(length(Snp)/ncores)+1)/3
      chunks <- round(seq(1, length(Snp), by = chunk),0)



      xcmaPara <- function(chunks,chunk,Snp,DataDir,genosnp,P,covarfile,Samp){

        ## Chunkfile create
        print(chunks)
        if (length(Snp)>=(chunks+chunk)){
          snp_names <- Snp[chunks:(chunks+chunk)]
        }else{
          snp_names <- Snp[chunks:length(Snp)]
        }
        ## without covariate file
        ## with covariate file
        XCMAX4_data2 <- function(Snp,DataDir,genosnp,P,covarfile,Samp){
          #XCMAX4_data2 <- function(Snp,DataDir,ResultDir,P,covarfile){

          #gdsfmt::showfile.gds(closeall=TRUE)
          #genofile <- SNPRelate::snpgdsOpen(f)
          # bedfile<- paste0(ResultDir,"/","XPLINK.bed")
          # famfile<- paste0(ResultDir,"/","XPLINK.fam")
          # bimfile<- paste0(ResultDir,"/","XPLINK.bim")
          #
          # ### Convert PLINK FILES INTO GDS FORMAT FOR R
          #
          # f<-SNPRelate::snpgdsBED2GDS(bedfile, famfile, bimfile, paste0(ResultDir,"/","IPCgeno.gdsX"))
          # print("f done")
          # # Open the GDS file and
          #genofile <- SNPRelate::snpgdsOpen(f)
          covarfile1 <- read.table(file = paste0(DataDir,"/",covarfile),
                                   stringsAsFactors = FALSE,
                                   header = TRUE)
          #g <- data.table::as.data.table(SNPRelate::snpgdsGetGeno(genofile, snp.id = Snp, verbose = FALSE))
          genosnp <- as.data.frame(genosnp)
          g <- genosnp[,Snp,drop = FALSE]
          colnames(g) <- c("Snp")
          #Samp <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id"))
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

          Z <- XCMAX4(Datatest)
          pval <- Z$`p-value`
          statistic <- Z$statictic
          SE <- Z$`standard-error`
          result <- as.data.frame(cbind(Snp,statistic,pval,SE))
          colnames(result) <- c("SNP","STAT","P","SE")
          return(result)
        }

        Result_pval1 <- data.table::rbindlist(lapply(snp_names,XCMAX4_data2,DataDir=DataDir,genosnp=genosnp,P=P,covarfile=covarfile,Samp=Samp))
        gc(reset = TRUE)
        return(Result_pval1)
      }

      #### Parallel computation
      gc(reset = TRUE)
      print("cluster making started")
      cl <- parallel::makeCluster(mc <- getOption("cl.cores", ncores),type = "FORK")

      invisible(parallel::clusterEvalQ(cl, library(data.table)))
      invisible(parallel::clusterEvalQ(cl, library(parallel)))
      invisible(parallel::clusterEvalQ(cl, library(SNPRelate)))
      print("cluster export started")
      parallel::clusterExport(cl=cl,NULL,envir=environment())
      print("cluster making done")

      Result_pval <- data.table::rbindlist(parallel::parLapply(cl,chunks,xcmaPara,chunk=chunk,Snp=Snp,DataDir=DataDir,genosnp=genosnp,P=P,covarfile=covarfile,Samp=Samp))
      parallel::stopCluster(cl)
      gc(reset = TRUE)
      print("clusters stopped.")

    }

    Result_pval$CHR <- "23"
    Result_pval$TEST <- "XCGA"
    Result_pval$BETA <- "NA"

    Result_pval$L95 <- "NA"
    Result_pval$U95 <- "NA"
    Xchr_Result_pval <- Result_pval
    save(Xchr_Result_pval, file = paste0(ResultDir,"/Xchr_Result_pval.Rda"))

    rm(Xchr_Result_pval)
    gc(reset = TRUE)
  }

  ## Performing Autosome WAS
  print("Running XCGA model for autosome.")
  XWAS_2 <- autoFun(DataDir = DataDir, ResultDir = ResultDir, finput = "FilteredX", trait = trait, standard_beta = standard_beta, perm = perm, mperm = mperm, genedrop = genedrop,
                    noxsex = noxsex, covarfile = covarfile, interaction =interaction,
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
  #######
  if (file.exists(paste0(ResultDir,"/","XChrRun.assoc.logistic")) == TRUE) {
    XWAS_3 <-
      read.table(file = paste0(ResultDir,"/","XChrRun.assoc.logistic"),
                 stringsAsFactors = FALSE,
                 header = TRUE)

    XWAS_4 <- XWAS_3[,c(2:4,6)]

    print("line 4836 reached")
    XWAS_5 <- unique(merge(Result_pval,XWAS_4,by = "SNP"))
    XWAS_6 <- XWAS_5[,c(5,1,10,11,6,12,7,4,8,9,2,3)]
    XWAS <- rbind(XWAS_2,XWAS_6)
    XWAS$CHR <- as.integer(as.character(XWAS$CHR))

    ## Plots
    XWAS_ADD <- na.omit(XWAS[XWAS$TEST == "ADD" | XWAS$TEST == "XCGA",c("SNP","CHR","BP","P")])
    XWAS_ADD$P <- as.numeric(as.character(XWAS_ADD$P))
    # From p-values, calculate chi-squared statistic
    chisq <- qchisq(1-XWAS_ADD$P,1)
    lamdaGC <- median(chisq)/qchisq(0.5,1)

    XWAS_ADD_X <- na.omit(XWAS_ADD[XWAS_ADD$CHR == 23,])
    chisq1 <- qchisq(1-XWAS_ADD_X$P,1)
    lamdaGC1 <- median(chisq1)/qchisq(0.5,1)

    if (plot.jpeg[1] == TRUE){
      options(bitmapType='cairo')
      grDevices::jpeg(paste0(ResultDir,"/",plotname,".jpeg"),  width = 20,
                      height = 10,
                      units = 'in',
                      res = 300)
      graphics::par(mfrow = c(2, 2))

      suppressWarnings(qqman::manhattan(XWAS_ADD, ylim = c(0,16),suggestiveline = suggestiveline, genomewideline = genomewideline, annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of GWAS"))
      suppressWarnings(qqman::qq(XWAS_ADD$P, main = paste0(("Q-Q plot of GWAS p-values with GIF = "), round(lamdaGC,3))))

      suppressWarnings(qqman::manhattan(XWAS_ADD_X, ylim = c(0,16),suggestiveline = suggestiveline, genomewideline = genomewideline, annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of XWAS"))
      suppressWarnings(qqman::qq(XWAS_ADD_X$P, main = paste0(("Q-Q plot of XWAS p-values with GIF = "), round(lamdaGC1,3))))

      dev.off()
    }else if (plot.jpeg[1] == FALSE){
      par(mar = c(1, 1, 1, 1))
      # graphics::par(mfrow = c(2, 2))
      suppressWarnings(qqman::manhattan(XWAS_ADD, ylim = c(0,16), suggestiveline = suggestiveline, genomewideline = genomewideline, annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of GWAS"))
      suppressWarnings(qqman::qq(XWAS_ADD$P, main = paste0(("Q-Q plot of GWAS p-values with GIF = "), round(lamdaGC,3))))

      suppressWarnings(qqman::manhattan(XWAS_ADD_X, ylim = c(0,16), suggestiveline = suggestiveline, genomewideline = genomewideline, annotatePval = snp_pval, annotateTop = annotateTopSnp, main = "Manhattan plot of XWAS"))
      suppressWarnings(qqman::qq(XWAS_ADD_X$P, main = paste0(("Q-Q plot of XWAS p-values with GIF = "), round(lamdaGC1,3))))
      #graphics::par(mfrow = c(1, 1))
    }
    #return(XWAS)
    GXWAS <- XWAS
    save(GXWAS, file = paste0(ResultDir,"/GXWAS_XCGA.Rda"))
    print(paste0("A dataframe named GXWAS_XCGA.Rda is saved in ",ResultDir))

  } else if (file.exists(paste0(ResultDir,"/","XChrRun.assoc.logistic")) == FALSE) {
    print("XWAS cannot be performed. Check the XChrRun.log file for checking the error.")
  }

}

############
## Function 16
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


## Function 17
ComputeLDSC <- function(snpld, test.df_beta, ncores, LDSC_blocks){

   #library(data.table)
  snps1 <- data.table::data.table(intersect(unique(snpld$SNP_A),unique(test.df_beta$rsid)))
  colnames(snps1) <- "SNP"
  snps2 <- data.table::data.table(intersect(unique(snpld$SNP_B),unique(test.df_beta$rsid)))
  colnames(snps2) <- "SNP"
  snps <- data.table::data.table(intersect(unique(snps1$SNP),unique(snps2$SNP)))
  colnames(snps) <- "SNP"
  snpld <- data.table::as.data.table(snpld)
  snpld1 <- merge(snpld,snps, by.x = "SNP_A",by.y="SNP")
  snpld1 <- merge(snpld1,snps, by.x = "SNP_B",by.y="SNP")
  snpld1 <- snpld1[,c(2,1,7)]
  finalsnp <- data.table::as.data.table(snpld1$SNP_A)
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


}


## Function 18

setupGCTA <- function(wdir) {
  #Set working directory

  #xdir <- getwd()
  #setwd(wdir)
  #Specify operating system
  OS <- Sys.info()['sysname']
  if (OS == "Linux") {

  utils::download.file(destfile = paste0(wdir,"/","gcta-1.94.1-linux-kernel-3-x86_64.zip"),
                         "https://github.com/boseb/bose_binaries/raw/main/gcta-1.94.1-linux-kernel-3-x86_64.zip", quiet = TRUE,
    )

    utils::unzip(paste0(wdir,"/","gcta-1.94.1-linux-kernel-3-x86_64.zip"), exdir = wdir, junkpaths = TRUE)
    Sys.chmod(paste0(wdir,"/gcta-1.94.1"), mode = "0777", use_umask = TRUE)
    print("GCTA program is succesfully setup.")
    invisible(file.remove(paste0(wdir,"/","MIT_License.txt")))
    invisible(file.remove(paste0(wdir,"/","README.txt")))

  } else if (OS == "Windows") {
print("Currently this function maynot work on Windows as GCTA for windows is down. Pleas use linux environment.")
    utils::download.file(destfile = paste0(wdir,"/","gcta-1.94.1"),
                         "https://figshare.com/ndownloader/files/42229122", quiet = TRUE,mode = "wb",
    )

    #utils::unzip(paste0(wdir,"/","gcta-1.94.1-Win-x86_64.zip"), exdir = wdir)
    Sys.chmod(paste0(wdir,"/gcta-1.94.1"), mode = "0777", use_umask = TRUE)
    print("GCTA program is succesfully setup.")

  } else if (OS == "macOS"){

    utils::download.file(destfile = paste0(wdir,"/","gcta-1.94.1-macOS-x86_64.zip"),
                         "https://github.com/boseb/bose_binaries/blob/main/gcta-1.94.1-macOS-x86_64.zip", quiet = TRUE,
    )

    utils::unzip(paste0(wdir,"/","gcta-1.94.1-macOS-x86_64.zip"), exdir = wdir, junkpaths = TRUE)
    Sys.chmod(paste0(wdir,"/gcta-1.94.1"), mode = "0777", use_umask = TRUE)
    print("GCTA program is succesfully setup.")
    invisible(file.remove(paste0(wdir,"/","MIT_License.txt")))
    invisible(file.remove(paste0(wdir,"/","README.txt")))

  }else {
    return(print("OS is not found by Sys.info()['sysname']."))
    print(
      "GCTA excecutable needs to be in the 'wdir', Please check https://www.cog-genomics.org/plink2/"
    )
    print(
      "After unzipping, please type in R console: system('chmod 755 plink')' for OS = Linux and for OS = Windows"
    )
  }
  OS <- Sys.info()['sysname']


}
## Function 19
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

    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./gcta-1.94.1"),
      args = c(
        "--bfile",
        paste0(DataDir,"/",finput),
        autosome,
        chr,CHRnum,
        minmaf,
        minmafval,
        maxmaf,
        maxmafval,
        grmcutoff,
        crip,
        "--make-grm",
        "--out",
        paste0(ResultDir,"/","test")
      ),
      std_out = FALSE,
      std_err = FALSE
    ))

  }else {

    partGRMfun <- function(i){

      invisible(sys::exec_wait(
        paste0(ResultDir,"/","./gcta-1.94.1"),
        args = c(
          "--bfile",
          paste0(DataDir,"/",finput),
          autosome,
          chr,CHRnum,
          minmaf,
          minmafval,
          maxmaf,
          maxmafval,
          grmcutoff,
          crip,
          "--make-grm-part",
          nGRM,
          i,
          "--thread-num",
          5,
          "--out",
          paste0(ResultDir,"/","test")
        ),
        std_out = FALSE,
        std_err = FALSE
      ))
    #gcta64 --bfile test --autosome --make-grm --make-grm-part nGRM i --thread-num 5 --out test
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

## Function 20
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

    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./gcta-1.94.1"),
      args = c(
        "--bfile",
        paste0(DataDir,"/",finput),
        minmaf,
        minmafval,
        maxmaf,
        maxmafval,
        "--make-grm-xchr",
        "--out",
        paste0(ResultDir,"/","xtest")
      ),
      std_out = FALSE,
      std_err = FALSE
    ))



    #system("./gcta-1.94.1  --bfile test  --autosome  --make-grm  --out test")
    #gcta64  --bfile test  --autosome  --make-grm  --out test
  }else {

    partGRM <- function(i){

      invisible(sys::exec_wait(
        paste0(ResultDir,"/","./gcta-1.94.1"),
        args = c(
          "--bfile",
          paste0(DataDir,"/",finput),
          minmaf,
          minmafval,
          maxmaf,
          maxmafval,
          "--make-grm-xchr-part",
          nGRM,
          i,
          "--thread-num",
          5,
          "--out",
          paste0(ResultDir,"/","xtest")
        ),
        std_out = FALSE,
        std_err = FALSE
      ))
      #gcta64 --bfile test --autosome --make-grm --make-grm-part nGRM i --thread-num 5 --out test
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


## Function 21
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

      invisible(sys::exec_wait(
        paste0(ResultDir,"/","./gcta-1.94.1"),
        args = c(
          "--reml",
          "--reml-alg",REMLalgo,
          "--reml-maxit",nitr,
          "--grm",
          paste0(ResultDir,"/",grmfile),
          pheno, paste0(ResultDir,"/",phenofile),
          preval, prevalval,
          catcovar, catcovarval,
          quantcovar, quantcovarval,
          "--reml-no-lrt",
          "--reml-no-constrain",
          "--thread-num",ncores,
          "--out",
          paste0(ResultDir,"/",chr,"test_reml")
        ),
        std_out = FALSE,
        std_err = FALSE
      ))

      if (file.exists(paste0(ResultDir,"/",chr,"test_reml.hsq"))){

        resultREML <- na.omit(data.table::fread(paste0(ResultDir,"/",chr,"test_reml.hsq"),fill = TRUE))
        return(resultREML)
      }else{
        print(grep("Error", readLines(paste0(ResultDir,"/", chr,"test_reml.log")), value = TRUE))
        print("An error occurs, please check the models, use byCHR = TRUE, check different options, SNP partitioning, or quality of the  data")
        x <- data.frame(NA,NA,NA)
        colnames(x) <- c("Source", "Variance","SE")
        return()
        }

    #gcta64  --reml  --grm test  --pheno test_cc.phen  --prevalence 0.01  --qcovar test_10PCs.txt  --out test_cc



}

## Function 22
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
        print("An error occurs, please check the models, use byCHR = TRUE, check different options, SNP partitioning or quality of the data")
        x <- data.frame(NA,NA,NA)
        colnames(x) <- c("Source", "Variance","SE")
        return()
        }

    #gcta64  --reml  --grm test  --pheno test_cc.phen  --prevalence 0.01  --qcovar test_10PCs.txt  --out test_cc

}


## Function 23
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
## Function 24

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

## Function 25
GeneProtein <- function(hg,chromosome){ ## Automatically using HG data from extdata
  if (hg == "hg19"){
    DataDir <- system.file("extdata", package = "GXwasR")

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

## Function 26
PlotHeritability <- function(Hdata,miMAF,maMAF){
  #create plot with regression line, regression equation, Pearson correlation and p-value.
  #data=result2
  #y=Variance
  p1<- ggplot2::ggplot(data=Hdata, ggplot2::aes(x=size_mb, y=Variance)) +
    ggplot2::labs(title="Chromosome-wise heritability",
         x = "Chromosome length (mb)", y = "Heritability") +
    ggplot2::geom_smooth(method="lm",level=0.95) +      ## For CI
    ggplot2::geom_point(size = 5, shape = 21, fill = "white") +
    ggplot2::geom_text(label = Hdata$chromosome,vjust = 0.4, size = 3)+
    ggpubr::stat_regline_equation(label.x=5, label.y= max(Hdata$Variance + 0.016)) +
    ggpubr::stat_cor(method = "pearson", label.x = 5, label.y = max(Hdata$Variance + 0.01))

  p2 <- ggplot2::ggplot(data=Hdata, ggplot2::aes(x=snp_proportion, y=Variance)) +
    ggplot2::labs(title="SNP proportion per chromosome vs heritability",
         x = "snp_proportion", y = "Heritability") +
    ggplot2::geom_smooth(method="lm") +
    ggplot2::geom_point(size = 5, shape = 21, fill = "white") +
    ggplot2::geom_text(label = Hdata$chromosome,vjust = 0.4, size = 3)+
    ggpubr::stat_regline_equation(label.x= median(Hdata$snp_proportion), label.y= max(Hdata$Variance + 0.016)) +
    ggpubr::stat_cor(method = "pearson", label.x = 5, label.y = max(Hdata$Variance + 0.01))

  p3 <- ggplot2::ggplot(data=Hdata, ggplot2::aes(x=no.of.genes, y=Variance)) +
    ggplot2::labs(title="Number of genes per chromosome vs heritability",
         x = "no.of.genes", y = "Heritability") +
    ggplot2::geom_smooth(method="lm") +
    ggplot2::geom_point(size = 5, shape = 21, fill = "white") +
    ggplot2::geom_text(label = Hdata$chromosome,vjust = 0.4, size = 3)+
    ggpubr::stat_regline_equation(label.x= median(Hdata$no.of.genes), label.y= max(Hdata$Variance + 0.016)) +
    ggpubr::stat_cor(method = "pearson", label.x = 5, label.y = max(Hdata$Variance + 0.01))

  p4 <- ggplot2::ggplot(data=Hdata, ggplot2::aes(x=no.of.proteins, y=Variance)) +
    ggplot2::labs(title="Number of proteins per chromosome vs heritability",
         x = "no.of.proteins", y = "Heritability") +
    ggplot2::geom_smooth(method="lm") +
    ggplot2::geom_point(size = 5, shape = 21, fill = "white") +
    ggplot2::geom_text(label = Hdata$chromosome,vjust = 0.4, size = 3)+
    ggpubr::stat_regline_equation(label.x= median(Hdata$no.of.proteins), label.y= max(Hdata$Variance + 0.016)) +
    ggpubr::stat_cor(method = "pearson", label.x = 5, label.y = max(Hdata$Variance + 0.01))

  plot1 <- ggpubr::ggarrange(p1, p2,p3,p4,
                     labels = c("A", "B","C","D"),
                     ncol = 2, nrow = 2)

  print(ggpubr::annotate_figure(plot1, top = ggpubr::text_grob(paste0(miMAF,",",maMAF),
                                               color = "red", face = "bold", size = 10)))
}


# Function 27
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

## Function 28
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

## Function 29

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
    #library(ggplot2)
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

## Function 30
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

## Function 31 ##PROBLEM
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

    invisible(sys::exec_wait(
      paste0(ResultDir,"/","./gcta-1.94.1"),
      args = c(
        "--reml-bivar",
        "--reml-alg",
        REMLalgo,
        "--reml-maxit",
        nitr,
        "--mgrm",
        paste0(ResultDir,"/","multi_GRMs.txt"),
        "--pheno",
        paste0(ResultDir,"/",phenofile),
        catcovar,
        catcovarval,
        quantcovar,
        quantcovarval,
        ExResi,
        "--reml-no-constrain",
        "--thread-num",
        ncores,
        "--out",
        paste0(ResultDir,"/","test_bireml")
      ),
      std_out = FALSE,
      std_err = FALSE
    ))

    if (file.exists(paste0(ResultDir,"/test_bireml.hsq"))){

      resultREML <- na.omit(data.table::fread(paste0(ResultDir,"/test_bireml.hsq"),fill = TRUE))
      return(resultREML)
    }else{
      print(grep("Error", readLines(paste0(ResultDir,"/test_bireml.log")), value = TRUE))
      print("An error occurs, please try byCHR = TRUE, check different options, SNP partitioning or quality of the data")
      x <- data.frame(NA,NA,NA)
      colnames(x) <- c("Source", "Variance","SE")
      return(x)
      }



}



## Function 32
## Function 121
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

## Function 33
outersect <- function(x, y, ...) {
  big.vec <- c(x, y, ...)
  duplicates <- big.vec[duplicated(big.vec)]
  setdiff(big.vec, unique(duplicates))
}

## Function 34
#' Title
#'
#' @param refdata dataset name (prefix)
#' @param wdir  directory name
#'
#' @return NULL
#' @export
#'
#' @examples
#' ## Not Run
#' # Download_reference(refdata = "HapMapIII_NCBI36", wdir = DataDir)
Download_reference <- function(refdata, wdir){
  OS <- Sys.info()['sysname']
  options(timeout=200)
  if (refdata == "HapMapIII_NCBI36"){
    if (OS == "Linux"){

    utils::download.file(destfile = paste0(wdir,"/",refdata,".zip"),
                       "https://figshare.com/ndownloader/files/40585145", quiet = TRUE,
    )
    }else if (OS == "Windows"){

      utils::download.file(destfile = paste0(wdir,"/",refdata,".zip"),
                           "https://figshare.com/ndownloader/files/40585145", quiet = TRUE,mode = "wb",
      )
    }else if (OS == "macOS"){
      utils::download.file(destfile = paste0(wdir,"/",refdata,".zip"),
                           "https://figshare.com/ndownloader/files/40585145", quiet = TRUE,
      )
    }
    utils::unzip(paste0(wdir,"/",refdata,".zip"), exdir = wdir, junkpaths = TRUE)
    utils::unzip(paste0(wdir,"/",refdata,".zip"), exdir = wdir)
    invisible(file.remove(paste0(wdir,"/",refdata,".zip")))

    }else if(refdata == "ThousandGenome"){
      options(timeout=200)##Takes more time, that's why this line
      if (OS == "Linux"){
        utils::download.file(destfile = paste0(wdir,"/",refdata,".tar.gz"),
                             "https://figshare.com/ndownloader/files/40728539", quiet = TRUE,
        )
      }else if (OS == "Windows"){
        utils::download.file(destfile = paste0(wdir,"/",refdata,".tar.gz"),
                             "https://figshare.com/ndownloader/files/40728539", quiet = TRUE,mode = "wb",
        )
      }else if (OS == "macOS"){
        utils::download.file(destfile = paste0(wdir,"/",refdata,".tar.gz"),
                             "https://figshare.com/ndownloader/files/40728539", quiet = TRUE,
        )
      }

      untar(paste0(wdir,"/",refdata,".tar.gz"),exdir = wdir)
      invisible(file.remove(paste0(wdir,"/",refdata,".tar.gz")))
    }

  print(paste0("Reference:", refdata, " downloaded."))

}

##gmirror function from hudson

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
  backpanel1 <- ifelse(background=="white", "NULL", "ggplot2::geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = yaxismin1, ymax = Inf, fill=factor(shademap)), alpha = 0.5)" )
  backpanel2 <- ifelse(background=="white", "NULL", "ggplot2::geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = yaxismin2, ymax = Inf, fill=factor(shademap)), alpha = 0.5)" )

  #Start plotting
  #TOP PLOT
  p1 <- ggplot() + eval(parse(text=backpanel1))
  #Add shape info if available
  if("Shape" %in% topn){
    p1 <- p1 + ggplot2::geom_point(data=d_order[d_order$Location=="Top",], aes(x=pos_index, y=pval, color=factor(Color), shape=factor(Shape)), alpha=opacity)
  } else {
    p1 <- p1 + ggplot2::geom_point(data=d_order[d_order$Location=="Top",], aes(x=pos_index, y=pval, color=factor(Color)), alpha=opacity)
  }
  p1 <- p1 + ggplot2::scale_x_continuous(breaks=lims$av, labels=lims$Color, expand=c(0,0))
  if(chrblocks==TRUE){
    p1 <- p1 + ggplot2::geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = -Inf, ymax = min(d_order$pval), fill=as.factor(Color)), alpha = 1)
  }
  p1 <- p1 + ggplot2::scale_colour_manual(name = "Color", values = newcols) + ggplot2::scale_fill_manual(name = "Color", values = newcols)
  p1 <- p1 + ggplot2::theme(panel.grid.minor.x = element_blank(), panel.grid.major.x=element_blank(), axis.title.x=element_blank(), legend.position="top", legend.title=element_blank())

  #BOTTOM PLOT
  p2 <- ggplot() + eval(parse(text=backpanel2))
  #Add shape info if available
  if("Shape" %in% bottomn){
    p2 <- p2 + ggplot2::geom_point(data=d_order[d_order$Location=="Bottom",], aes(x=pos_index, y=pval, color=factor(Color), shape=factor(Shape)), alpha=opacity)
  } else {
    p2 <- p2 + ggplot2::geom_point(data=d_order[d_order$Location=="Bottom",], aes(x=pos_index, y=pval, color=factor(Color)), alpha=opacity)
  }
  p2 <- p2 + ggplot2::scale_x_continuous(breaks=lims$av, labels=lims$Color, expand=c(0,0))
  if(chrblocks==TRUE){
    p2 <- p2 + ggplot2::geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = -Inf, ymax = min(d_order$pval), fill=as.factor(Color)), alpha = 1)
  }
  p2 <- p2 + ggplot2::scale_colour_manual(name = "Color", values = newcols) + ggplot2::scale_fill_manual(name = "Color", values = newcols)
  p2 <- p2 + ggplot2::theme(axis.text.x=element_text(angle=90), panel.grid.minor.x = element_blank(), panel.grid.major.x=element_blank(), axis.title.x=element_blank(), legend.position="bottom", legend.title=element_blank())

  #Highlight if given
  if(!missing(highlight_snp)){
    if("Shape" %in% topn){
      p1 <- p1 + ggplot2::geom_point(data=d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Top", ], aes(x=pos_index, y=pval, shape=Shape), colour=highlighter)
      p1 <- p1 + ggplot2::guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p1 <- p1 + ggplot2::geom_point(data=d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Top", ], aes(x=pos_index, y=pval), colour=highlighter)
    }
    if("Shape" %in% bottomn){
      p2 <- p2 + ggplot2::geom_point(data=d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Bottom", ], aes(x=pos_index, y=pval, shape=Shape), colour=highlighter)
      p2 <- p2 + ggplot2::guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p2 <- p2 + ggplot2::geom_point(data=d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Bottom", ], aes(x=pos_index, y=pval), colour=highlighter)
    }
  }
  if(!missing(highlight_p)){
    if("Shape" %in% topn){
      p1 <- p1 + ggplot2::geom_point(data=d_order[d_order$pvalue < highlight_p[1] & d_order$Location=="Top", ], aes(x=pos_index, y=pval, shape=Shape), colour=highlighter)
      p1 <- p1 + ggplot2::guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p1 <- p1 + ggplot2::geom_point(data=d_order[d_order$pvalue < highlight_p[1] & d_order$Location=="Top", ], aes(x=pos_index, y=pval), colour=highlighter)
    }
    if("Shape" %in% bottomn){
      p2 <- p2 + ggplot2::geom_point(data=d_order[d_order$pvalue < highlight_p[2] & d_order$Location=="Bottom", ], aes(x=pos_index, y=pval, shape=Shape), colour=highlighter)
      p2 <- p2 + ggplot2::guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p2 <- p2 + ggplot2::geom_point(data=d_order[d_order$pvalue < highlight_p[2] & d_order$Location=="Bottom", ], aes(x=pos_index, y=pval), colour=highlighter)
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
      p1 <- p1 + ggplot2::geom_text(data=d_order[d_order$pvalue < annotate_p[1] & d_order$Location=="Top",], aes(pos_index,pval,label=SNP))
      p2 <- p2 + ggplot2::geom_text(data=d_order[d_order$pvalue < annotate_p[2] & d_order$Location=="Bottom",], aes(pos_index,pval,label=SNP))
    } else {
      p1 <- p1 + ggrepel::geom_text_repel(data=d_order[d_order$pvalue < annotate_p[1] & d_order$Location=="Top",], aes(pos_index,pval,label=SNP))
      p2 <- p2 + ggrepel::geom_text_repel(data=d_order[d_order$pvalue < annotate_p[2] & d_order$Location=="Bottom",], aes(pos_index,pval,label=SNP))
    }
  }
  if(!missing(annotate_snp)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE){
      print("Consider installing 'ggrepel' for improved text annotation")
      p1 <- p1 + ggplot2::geom_text(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Top",], aes(pos_index,pval,label=SNP))
      p2 <- p2 + ggplot2::geom_text(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Bottom",], aes(pos_index,pval,label=SNP))
    } else {
      p1 <- p1 + ggrepel::geom_text_repel(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Top",], aes(pos_index,pval,label=SNP))
      p2 <- p2 + ggrepel::geom_text_repel(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Bottom",], aes(pos_index,pval,label=SNP))
    }
  }
  #Add title and y axis title
  p1 <- p1 + ylab(yaxislab1)
  p2 <- p2 + ylab(yaxislab2)

  #Format
  if(chrblocks==TRUE){
    if(freey==TRUE){
      print("Sorry, drawing chrblocks with freey=TRUE is currently unsupported and will be ignored.")
    } else {
      p1 <- p1+ggplot2::theme(axis.text.x = element_text(vjust=1),axis.ticks.x = element_blank())+ggplot2::ylim(c(yaxismin1,yaxismax1))
      p2 <- p2+ggplot2::scale_y_reverse(limits=c(yaxismax2, yaxismin2)) + ggplot2::theme(axis.text.x = ggplot2::element_blank(),axis.ticks.x = ggplot2::element_blank())
    }
  } else {
    p1 <- p1+ggplot2::theme(axis.text.x = element_text(vjust=1),axis.ticks.x = element_blank())+ scale_y_continuous(limits=c(yaxismin1, yaxismax1),expand=expansion(mult=c(0,0.1)))
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
  p <- gridExtra::grid.arrange(arrangeGrob(p1, top=toptitle), arrangeGrob(p2, bottom=bottomtitle), padding=0, heights=c(hgtratio,1-hgtratio))
  ggplot2::ggsave(p, filename=paste0(file, ".", type), dpi=res, units="in", height=hgt, width=wi)
  return(p)
}
