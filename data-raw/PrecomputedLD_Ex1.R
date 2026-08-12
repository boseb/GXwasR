## code to prepare `PrecomputedLD_Ex1` dataset
# load("PrecomputedLD_hg19_EUR.Rda")
data("Summary_Stat_Ex1", package = "GXwasR")

test.sumstats <- na.omit(
    Summary_Stat_Ex1[
        Summary_Stat_Ex1$TEST == "ADD",
        c(seq_len(4), 6:8)
    ]
)

colnames(test.sumstats) <- c(
    "chr", "rsid", "pos", "a1",
    "n_eff", "beta", "beta_se"
)

PrecomputedLD_Ex1 <- PrecomputedLD_hg19_EUR[
    PrecomputedLD_hg19_EUR$SNP %in% test.sumstats$rsid,
    c("CHR", "SNP", "ld_score")
]
usethis::use_data(PrecomputedLD_Ex1, overwrite = TRUE)
