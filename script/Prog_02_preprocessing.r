##------------------------------------------------------------
# DNA METHYLATION PREPROCESSING
# Brain IDH-mutant FFPE tissue
#
# Platforms:
#   - Illumina MethylationEPIC v1: 3 samples
#   - Illumina MethylationEPIC v2: 26 samples
#
# Current objective:
#   1. Load raw RGChannelSet objects
#   2. Verify and repair sample metadata
#   3. Calculate detection P-values
#   4. Perform sample-level QC
#   5. Normalize EPIC v1 and EPIC v2 separately using Noob
#   6. Filter unreliable probes
#   7. Generate Beta and M-value matrices
#   8. Save processed objects
#
# Important:
#   EPIC v1 and EPIC v2 are NOT combined in this script.
##------------------------------------------------------------
##################################################
## library
#################################################
library(minfi)
library(SummarizedExperiment)
library(S4Vectors)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

##################################################
## functions
##################################################
run_pca <- function(m_values, targets_data, n_probes = 10000) {

  common_samples <- intersect(
    colnames(m_values),
    targets_data$Sample_Name
  )

  m_use <- m_values[, common_samples, drop = FALSE]

  probe_var <- apply(
    m_use,
    1,
    var,
    na.rm = TRUE
  )

  probe_var <- probe_var[is.finite(probe_var)]

  top_probes <- names(
    sort(probe_var, decreasing = TRUE)
  )[seq_len(min(n_probes, length(probe_var)))]

  pca <- prcomp(
    t(m_use[top_probes, , drop = FALSE]),
    center = TRUE,
    scale. = FALSE
  )

  variance_explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)

  pca_df <- data.frame(
    Sample_Name = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  ) %>%
    left_join(
      targets_data,
      by = "Sample_Name"
    )

  list(
    pca = pca,
    data = pca_df,
    variance_explained = variance_explained
  )
}

#################################################
## setup directory
#################################################
dir_input  <- "data/proc"
dir_output <- "result/data"
dir_figures <- "result/qc"

#################################################
## load data
#################################################
load(file.path(dir_input, 'RGSets_targets.RData'))
ls()

############################################################
## Calculate detection P-values
############################################################
detP_EPICv1 <- detectionP(rgSet_EPICv1)
detP_EPICv2 <- detectionP(rgSet_EPICv2)

dim(detP_EPICv1)
dim(detP_EPICv2)

stopifnot(
  identical(
    colnames(detP_EPICv1),
    sampleNames(rgSet_EPICv1)
  )
)

stopifnot(
  identical(
    colnames(detP_EPICv2),
    sampleNames(rgSet_EPICv2)
  )
)

#################################################
## Summarize sample-level detection performance
#################################################
detection_cutoff <- 0.01

sample_qc_EPICv1 <- data.frame(
  Sample_Name = colnames(detP_EPICv1),
  n_probes = nrow(detP_EPICv1),
  n_failed = colSums(detP_EPICv1 > detection_cutoff, na.rm = TRUE),
  prop_failed = colMeans(detP_EPICv1 > detection_cutoff, na.rm = TRUE),
  stringsAsFactors = FALSE
)

sample_qc_EPICv2 <- data.frame(
  Sample_Name = colnames(detP_EPICv2),
  n_probes = nrow(detP_EPICv2),
  n_failed = colSums(detP_EPICv2 > detection_cutoff, na.rm = TRUE),
  prop_failed = colMeans(detP_EPICv2 > detection_cutoff, na.rm = TRUE),
  stringsAsFactors = FALSE
)

sample_qc_EPICv1 <- sample_qc_EPICv1 %>%
  arrange(desc(prop_failed))

sample_qc_EPICv2 <- sample_qc_EPICv2 %>%
  arrange(desc(prop_failed))

sample_qc_EPICv1
sample_qc_EPICv2
 
sample_qc_EPICv1 <- sample_qc_EPICv1 %>%
  mutate(
    percent_failed = 100 * prop_failed
  )

sample_qc_EPICv2 <- sample_qc_EPICv2 %>%
  mutate(
    percent_failed = 100 * prop_failed
  )

sample_qc_EPICv1 %>%
  select(
    Sample_Name,
    n_failed,
    percent_failed
  )

sample_qc_EPICv2 %>%
  select(
    Sample_Name,
    n_failed,
    percent_failed
  )

#################################################
## Add phenotypes
#################################################
sample_qc_EPICv1 <- sample_qc_EPICv1 %>%
  left_join(
    targets_EPICv1 %>%
      select(
        Sample_Name,
        Subject,
        Array,
        Genome_build,
        Sex,
        Path_dx,
        Grade
      ),
    by = "Sample_Name"
  )

sample_qc_EPICv2 <- sample_qc_EPICv2 %>%
  left_join(
    targets_EPICv2 %>%
      select(
        Sample_Name,
        Subject,
        Array,
        Genome_build,
        Sex,
        Path_dx,
        Grade
      ),
    by = "Sample_Name"
  )

#################################################
## Flag poor samples
#################################################
sample_qc_EPICv1 <- sample_qc_EPICv1 %>%
  mutate(
    detection_flag = case_when(
      prop_failed > 0.05 ~ "Review: >5% failed",
      prop_failed > 0.01 ~ "Review: 1-5% failed",
      TRUE ~ "Pass"
    )
  )

sample_qc_EPICv2 <- sample_qc_EPICv2 %>%
  mutate(
    detection_flag = case_when(
      prop_failed > 0.05 ~ "Review: >5% failed",
      prop_failed > 0.01 ~ "Review: 1-5% failed",
      TRUE ~ "Pass"
    )
  )

table(sample_qc_EPICv1$detection_flag)
table(sample_qc_EPICv2$detection_flag)

##################################################
## Figures --- sample-level detection performance
##################################################
pdf(
  file.path(dir_figures, "EPICv1_detection_failure_by_sample.pdf"),
  width = 8,
  height = 7
)

barplot(
  sample_qc_EPICv1$percent_failed,
  names.arg = sample_qc_EPICv1$Sample_Name,
  las = 2,
  ylab = "Failed probes (%)",
  main = "EPIC v1: detection P-value failures"
)

abline(h = 1, lty = 2)
abline(h = 5, lty = 3)

dev.off()

pdf(
  file.path(dir_figures, "EPICv2_detection_failure_by_sample.pdf"),
  width = 12,
  height = 7
)

barplot(
  sample_qc_EPICv2$percent_failed,
  names.arg = sample_qc_EPICv2$Sample_Name,
  las = 2,
  ylab = "Failed probes (%)",
  main = "EPIC v2: detection P-value failures"
)

abline(h = 1, lty = 2)
abline(h = 5, lty = 3)

dev.off()

############################################################
## Assess raw methylated and unmethylated signal
############################################################
mSet_raw_EPICv1 <- preprocessRaw(rgSet_EPICv1)
mSet_raw_EPICv2 <- preprocessRaw(rgSet_EPICv2)

qc_signal_EPICv1 <- getQC(mSet_raw_EPICv1)
qc_signal_EPICv2 <- getQC(mSet_raw_EPICv2)

qc_signal_EPICv1
qc_signal_EPICv2

qc_signal_EPICv1 <- as.data.frame(qc_signal_EPICv1) %>%
  tibble::rownames_to_column("Sample_Name")

qc_signal_EPICv2 <- as.data.frame(qc_signal_EPICv2) %>%
  tibble::rownames_to_column("Sample_Name")

sample_qc_EPICv1 <- sample_qc_EPICv1 %>%
  left_join(
    qc_signal_EPICv1,
    by = "Sample_Name"
  )

sample_qc_EPICv2 <- sample_qc_EPICv2 %>%
  left_join(
    qc_signal_EPICv2,
    by = "Sample_Name"
  )

pdf(
  file.path(dir_figures, "EPICv1_raw_signal_QC.pdf"),
  width = 7,
  height = 6
)

plotQC(qc_signal_EPICv1)

dev.off()

pdf(
  file.path(dir_figures, "EPICv2_raw_signal_QC.pdf"),
  width = 7,
  height = 6
)

plotQC(qc_signal_EPICv2)

dev.off()

rgSet_EPICv1_qc <- rgSet_EPICv1
rgSet_EPICv2_qc <- rgSet_EPICv2[
  ,
  sampleNames(rgSet_EPICv2) != "GSM9325977"
]

dim(rgSet_EPICv1_qc)
dim(rgSet_EPICv2_qc)
sampleNames(rgSet_EPICv2_qc)

sample_qc_EPICv2 <- sample_qc_EPICv2 %>%
  mutate(
    final_qc_status = case_when(
      Sample_Name == "GSM9325977" ~ "Exclude",
      Sample_Name %in% c(
        "GSM9325997",
        "GSM9325994"
      ) ~ "Retain and monitor",
      TRUE ~ "Retain"
    ),
    qc_reason = case_when(
      Sample_Name == "GSM9325977" ~
        paste(
          "Excluded due to >5% detection failures",
          "and markedly reduced methylated and",
          "unmethylated median signal intensities."
        ),
      Sample_Name %in% c(
        "GSM9325997",
        "GSM9325994"
      ) ~
        "Elevated detection failures; retain pending post-normalization QC.",
      TRUE ~
        "Acceptable detection and signal-intensity QC."
    )
  )

write.csv(
  sample_qc_EPICv1,
  file.path(dir_output, "EPICv1_sample_QC.csv"),
  row.names = FALSE
)

write.csv(
  sample_qc_EPICv2,
  file.path(dir_output, "EPICv2_sample_QC.csv"),
  row.names = FALSE
)

save(
  detP_EPICv1,
  detP_EPICv2,
  mSet_raw_EPICv1,
  mSet_raw_EPICv2,
  qc_signal_EPICv1,
  qc_signal_EPICv2,
  sample_qc_EPICv1,
  sample_qc_EPICv2,
  file = file.path(
    dir_output,
    "raw_QC_results.RData"
  )
)

############################################################
## Normalize each platform separately using Noob
############################################################
mSet_noob_EPICv1 <- preprocessNoob(
  rgSet_EPICv1_qc,
  dyeCorr = TRUE,
  verbose = TRUE
)

mSet_noob_EPICv2 <- preprocessNoob(
  rgSet_EPICv2_qc,
  dyeCorr = TRUE,
  verbose = TRUE
)

############################################################
## Post-normalization QC: basic checks
############################################################
dim(mSet_noob_EPICv1)
dim(mSet_noob_EPICv2)

sampleNames(mSet_noob_EPICv1)
sampleNames(mSet_noob_EPICv2)

stopifnot(
  identical(
    sampleNames(mSet_noob_EPICv1),
    sampleNames(rgSet_EPICv1_qc)
  )
)

stopifnot(
  identical(
    sampleNames(mSet_noob_EPICv2),
    sampleNames(rgSet_EPICv2_qc)
  )
)

############################################################
## Generate Noob-normalized Beta and M-value matrices
############################################################

beta_noob_EPICv1 <- getBeta(mSet_noob_EPICv1)
beta_noob_EPICv2 <- getBeta(mSet_noob_EPICv2)

mval_noob_EPICv1 <- getM(mSet_noob_EPICv1)
mval_noob_EPICv2 <- getM(mSet_noob_EPICv2)

range(beta_noob_EPICv1, na.rm = TRUE)
range(beta_noob_EPICv2, na.rm = TRUE)

sum(is.na(beta_noob_EPICv1))
sum(is.na(beta_noob_EPICv2))

sum(!is.finite(mval_noob_EPICv1))
sum(!is.finite(mval_noob_EPICv2))

dim(beta_noob_EPICv1)
dim(beta_noob_EPICv2)

summary(beta_noob_EPICv1[, 1])
summary(beta_noob_EPICv2[, 1])

pdf(
  file.path(dir_figures, "EPICv1_noob_beta_density.pdf"),
  width = 8,
  height = 6
)

densityPlot(
  beta_noob_EPICv1,
  main = "EPIC v1: Beta-value density after Noob",
  legend = FALSE
)

dev.off()

pdf(
  file.path(dir_figures, "EPICv2_noob_beta_density.pdf"),
  width = 8,
  height = 6
)

densityPlot(
  beta_noob_EPICv2,
  main = "EPIC v2: Beta-value density after Noob",
  legend = FALSE
)

dev.off()

###########################################
## QC: after Noob normalization
###########################################
qc_noob_EPICv1 <- getQC(mSet_noob_EPICv1)
qc_noob_EPICv2 <- getQC(mSet_noob_EPICv2)

qc_noob_EPICv1
qc_noob_EPICv2

pdf(
  file.path(dir_figures, "EPICv1_noob_signal_QC.pdf"),
  width = 7,
  height = 7
)

plotQC(qc_noob_EPICv1)

dev.off()

pdf(
  file.path(dir_figures, "EPICv2_noob_signal_QC.pdf"),
  width = 8,
  height = 7
)

plotQC(qc_noob_EPICv2)

dev.off()

#################################################
## Run PCA
#################################################
targets_EPICv1_qc <- targets_EPICv1 %>%
  filter(Sample_Name %in% colnames(mval_noob_EPICv1))

targets_EPICv2_qc <- targets_EPICv2 %>%
  filter(Sample_Name %in% colnames(mval_noob_EPICv2))

pca_EPICv1_noob <- run_pca(
  mval_noob_EPICv1,
  targets_EPICv1_qc
)

pca_EPICv2_noob <- run_pca(
  mval_noob_EPICv2,
  targets_EPICv2_qc
)


pdf(
  file.path(dir_figures, "EPICv2_noob_PCA.pdf"),
  width = 8,
  height = 7
)

plot(
  pca_EPICv2_noob$data$PC1,
  pca_EPICv2_noob$data$PC2,
  pch = 19,
  xlab = paste0(
    "PC1 (",
    round(pca_EPICv2_noob$variance_explained[1], 1),
    "%)"
  ),
  ylab = paste0(
    "PC2 (",
    round(pca_EPICv2_noob$variance_explained[2], 1),
    "%)"
  ),
  main = "EPIC v2: PCA after Noob normalization"
)

text(
  pca_EPICv2_noob$data$PC1,
  pca_EPICv2_noob$data$PC2,
  labels = pca_EPICv2_noob$data$Sample_Name,
  pos = 3,
  cex = 0.65
)

dev.off()
