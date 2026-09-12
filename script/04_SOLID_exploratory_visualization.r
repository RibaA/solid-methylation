############################################################
## 04_SOLID_exploratory_visualization.r
##
## SOLID TUMOR–PLASMA EXPLORATORY ANALYSIS
##
## Purpose:
##   1. Load the frozen downstream analysis-ready object
##   2. Apply the existing valid-pair mask
##   3. Summarize tumor/plasma Beta values by patient
##   4. Summarize ΔBeta = Plasma - Tumor
##   5. Review valid-region availability by patient
##   6. Generate descriptive figures
##   7. Incorporate selected clinical metadata
##
## Clinical variables available:
##   Grade
##   Age
##   Sex
##   ECOG
##   Response_RANO
##   pfs / PFS_status
##   os / OS_status
##
## IMPORTANT:
##   - No new filtering is performed
##   - No PCA is performed
##   - No clustering is performed
##   - No DMR testing is performed
############################################################


############################################################
## 0. CLEAN SESSION
############################################################

rm(list = ls())

options(
  stringsAsFactors = FALSE,
  width = 180,
  scipen = 999,
  warn = 1
)

############################################################
## 1. REQUIRED PACKAGE
############################################################

required_packages <- c(
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {

  stop(
    "Missing package(s): ",
    paste(
      missing_packages,
      collapse = ", "
    ),
    "\nInstall with:\n",
    "install.packages(c(",
    paste(
      sprintf(
        '"%s"',
        missing_packages
      ),
      collapse = ", "
    ),
    "))"
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
})

############################################################
## 2. PROJECT DIRECTORIES
############################################################

project_dir <- "C:/solid-methylation"

input_file <- file.path(
  project_dir,
  "result",
  "03_matched_tissue_plasma",
  "SOLID_downstream_analysis_ready.rds"
)

exploratory_dir <- file.path(
  project_dir,
  "result",
  "03_matched_tissue_plasma",
  "exploratory"
)

figure_dir <- file.path(
  exploratory_dir,
  "figures"
)

table_dir <- file.path(
  exploratory_dir,
  "tables"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(
  file.exists(input_file)
)

############################################################
## 3. LOAD ANALYSIS-READY OBJECT
############################################################

obj <- readRDS(
  input_file
)

cat(
  "\n============================================\n"
)

cat(
  "SOLID EXPLORATORY VISUALIZATION\n"
)

cat(
  "============================================\n"
)

############################################################
## 4. BASIC VALIDATION
############################################################

required_components <- c(
  "tissue_beta",
  "plasma_beta",
  "delta_beta",
  "valid_pair_mask",
  "region_annotation",
  "patient_metadata"
)

missing_components <- setdiff(
  required_components,
  names(obj)
)

if (length(missing_components) > 0L) {

  stop(
    "Missing object components: ",
    paste(
      missing_components,
      collapse = ", "
    )
  )
}

tissue_beta <- obj$tissue_beta
plasma_beta <- obj$plasma_beta
delta_beta <- obj$delta_beta
valid_mask <- obj$valid_pair_mask

clinical <- as.data.frame(
  obj$patient_metadata
)

stopifnot(
  identical(
    dim(tissue_beta),
    dim(plasma_beta)
  ),
  identical(
    dim(tissue_beta),
    dim(delta_beta)
  ),
  identical(
    dim(tissue_beta),
    dim(valid_mask)
  ),
  ncol(tissue_beta) == 13L,
  nrow(tissue_beta) == 124961L,
  nrow(clinical) == 13L
)

patient_ids <- colnames(
  tissue_beta
)

stopifnot(
  identical(
    patient_ids,
    clinical$patient_id
  )
)

cat(
  "\nInput validation: PASS\n"
)

cat(
  "Patients:",
  ncol(tissue_beta),
  "\n"
)

cat(
  "Regions:",
  nrow(tissue_beta),
  "\n"
)

############################################################
## 5. APPLY VALID-PAIR MASK
##
## No observations are removed permanently.
## Masking is used only so all descriptive comparisons use
## the same valid tumor-plasma measurements.
############################################################

tissue_beta_valid <- tissue_beta
plasma_beta_valid <- plasma_beta
delta_beta_valid <- delta_beta

tissue_beta_valid[
  !valid_mask
] <- NA_real_

plasma_beta_valid[
  !valid_mask
] <- NA_real_

delta_beta_valid[
  !valid_mask
] <- NA_real_

############################################################
## 6. PATIENT-LEVEL SUMMARY
############################################################

patient_summary <- data.frame(
  patient_id = patient_ids,
  stringsAsFactors = FALSE
)

############################################################
## Valid-region counts
############################################################

patient_summary$n_valid_regions <- colSums(
  valid_mask,
  na.rm = TRUE
)

patient_summary$pct_valid_regions <-
  100 *
  patient_summary$n_valid_regions /
  nrow(valid_mask)

############################################################
## Tumor Beta summaries
############################################################

patient_summary$mean_tumor_beta <- vapply(
  seq_along(patient_ids),
  function(j) {
    mean(
      tissue_beta_valid[, j],
      na.rm = TRUE
    )
  },
  numeric(1)
)

patient_summary$median_tumor_beta <- vapply(
  seq_along(patient_ids),
  function(j) {
    median(
      tissue_beta_valid[, j],
      na.rm = TRUE
    )
  },
  numeric(1)
)

############################################################
## Plasma Beta summaries
############################################################

patient_summary$mean_plasma_beta <- vapply(
  seq_along(patient_ids),
  function(j) {
    mean(
      plasma_beta_valid[, j],
      na.rm = TRUE
    )
  },
  numeric(1)
)

patient_summary$median_plasma_beta <- vapply(
  seq_along(patient_ids),
  function(j) {
    median(
      plasma_beta_valid[, j],
      na.rm = TRUE
    )
  },
  numeric(1)
)

############################################################
## Delta Beta summaries
############################################################

patient_summary$mean_delta_beta <- vapply(
  seq_along(patient_ids),
  function(j) {
    mean(
      delta_beta_valid[, j],
      na.rm = TRUE
    )
  },
  numeric(1)
)

patient_summary$median_delta_beta <- vapply(
  seq_along(patient_ids),
  function(j) {
    median(
      delta_beta_valid[, j],
      na.rm = TRUE
    )
  },
  numeric(1)
)

patient_summary$q25_delta_beta <- vapply(
  seq_along(patient_ids),
  function(j) {
    quantile(
      delta_beta_valid[, j],
      probs = 0.25,
      na.rm = TRUE,
      names = FALSE
    )
  },
  numeric(1)
)

patient_summary$q75_delta_beta <- vapply(
  seq_along(patient_ids),
  function(j) {
    quantile(
      delta_beta_valid[, j],
      probs = 0.75,
      na.rm = TRUE,
      names = FALSE
    )
  },
  numeric(1)
)

############################################################
## Direction summaries
############################################################

patient_summary$pct_plasma_higher <- vapply(
  seq_along(patient_ids),
  function(j) {

    x <- delta_beta_valid[, j]

    x <- x[
      is.finite(x)
    ]

    100 *
      mean(
        x > 0
      )
  },
  numeric(1)
)

patient_summary$pct_tumor_higher <- vapply(
  seq_along(patient_ids),
  function(j) {

    x <- delta_beta_valid[, j]

    x <- x[
      is.finite(x)
    ]

    100 *
      mean(
        x < 0
      )
  },
  numeric(1)
)

############################################################
## 7. ADD CLINICAL METADATA
############################################################

clinical_for_merge <- clinical[
  ,
  c(
    "patient_id",
    "Grade",
    "Age",
    "Sex",
    "ECOG",
    "Response_RANO",
    "pfs",
    "PFS_status",
    "os",
    "OS_status"
  ),
  drop = FALSE
]

patient_summary <- merge(
  patient_summary,
  clinical_for_merge,
  by = "patient_id",
  all.x = TRUE,
  sort = FALSE
)


patient_summary <- patient_summary[
  match(
    patient_ids,
    patient_summary$patient_id
  ),
  ,
  drop = FALSE
]


stopifnot(
  identical(
    patient_summary$patient_id,
    patient_ids
  )
)


############################################################
## 8. SAVE PATIENT SUMMARY
############################################################

patient_summary_file <- file.path(
  table_dir,
  "SOLID_patient_exploratory_summary.tsv"
)

write.table(
  patient_summary,
  patient_summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


cat(
  "\nPatient-level exploratory summary:\n"
)

print(
  patient_summary
)


############################################################
## 9. SAMPLE-TYPE DISTRIBUTION SUMMARY
############################################################

distribution_summary <- rbind(

  data.frame(
    sample_type = "Tumor",
    mean_beta = mean(
      tissue_beta_valid,
      na.rm = TRUE
    ),
    median_beta = median(
      tissue_beta_valid,
      na.rm = TRUE
    ),
    q25_beta = quantile(
      tissue_beta_valid,
      0.25,
      na.rm = TRUE
    ),
    q75_beta = quantile(
      tissue_beta_valid,
      0.75,
      na.rm = TRUE
    )
  ),

  data.frame(
    sample_type = "Plasma",
    mean_beta = mean(
      plasma_beta_valid,
      na.rm = TRUE
    ),
    median_beta = median(
      plasma_beta_valid,
      na.rm = TRUE
    ),
    q25_beta = quantile(
      plasma_beta_valid,
      0.25,
      na.rm = TRUE
    ),
    q75_beta = quantile(
      plasma_beta_valid,
      0.75,
      na.rm = TRUE
    )
  )
)


distribution_summary_file <- file.path(
  table_dir,
  "SOLID_tumor_plasma_beta_distribution_summary.tsv"
)

write.table(
  distribution_summary,
  distribution_summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


cat(
  "\nTumor/plasma Beta summary:\n"
)

print(
  distribution_summary
)


############################################################
## 10. REPRODUCIBLE REGION SUBSET FOR DENSITY PLOTS
##
## Plotting every matrix cell is unnecessary.
## Use the same randomly selected regions for all patients
## and both sample types.
############################################################

set.seed(
  20260912
)

n_plot_regions <- min(
  20000L,
  nrow(tissue_beta_valid)
)

plot_region_index <- sort(
  sample(
    seq_len(
      nrow(tissue_beta_valid)
    ),
    size = n_plot_regions,
    replace = FALSE
  )
)


tissue_plot_matrix <- tissue_beta_valid[
  plot_region_index,
  ,
  drop = FALSE
]

plasma_plot_matrix <- plasma_beta_valid[
  plot_region_index,
  ,
  drop = FALSE
]


beta_density_data <- rbind(

  data.frame(
    beta = as.vector(
      tissue_plot_matrix
    ),
    sample_type = "Tumor"
  ),

  data.frame(
    beta = as.vector(
      plasma_plot_matrix
    ),
    sample_type = "Plasma"
  )
)


beta_density_data <- beta_density_data[
  is.finite(
    beta_density_data$beta
  ),
  ,
  drop = FALSE
]


beta_density_data$sample_type <- factor(
  beta_density_data$sample_type,
  levels = c(
    "Tumor",
    "Plasma"
  )
)


############################################################
## 11. FIGURE 1
## Tumor vs plasma Beta distributions
############################################################

p1 <- ggplot(
  beta_density_data,
  aes(
    x = beta,
    color = sample_type,
    fill = sample_type
  )
) +
  geom_density(
    alpha = 0.20,
    linewidth = 0.9
  ) +
  labs(
    title = "SOLID tumor and plasma methylation distributions",
    subtitle = paste0(
      "Representative subset of ",
      format(
        n_plot_regions,
        big.mark = ","
      ),
      " retained 1-kb regions"
    ),
    x = expression(beta),
    y = "Density",
    color = "Sample type",
    fill = "Sample type"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )


ggsave(
  file.path(
    figure_dir,
    "SOLID_tumor_plasma_beta_density.png"
  ),
  p1,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(
    figure_dir,
    "SOLID_tumor_plasma_beta_density.pdf"
  ),
  p1,
  width = 7,
  height = 5
)


############################################################
## 12. PREPARE PAIRED MEDIAN-BETA DATA
############################################################

median_beta_long <- rbind(

  data.frame(
    patient_id =
      patient_summary$patient_id,

    sample_type =
      "Tumor",

    median_beta =
      patient_summary$median_tumor_beta,

    Grade =
      patient_summary$Grade,

    Sex =
      patient_summary$Sex,

    stringsAsFactors = FALSE
  ),

  data.frame(
    patient_id =
      patient_summary$patient_id,

    sample_type =
      "Plasma",

    median_beta =
      patient_summary$median_plasma_beta,

    Grade =
      patient_summary$Grade,

    Sex =
      patient_summary$Sex,

    stringsAsFactors = FALSE
  )
)


median_beta_long$sample_type <- factor(
  median_beta_long$sample_type,
  levels = c(
    "Tumor",
    "Plasma"
  )
)


median_beta_long$Grade <- factor(
  median_beta_long$Grade
)

median_beta_long$Sex <- factor(
  median_beta_long$Sex
)


############################################################
## 13. FIGURE 2
## Paired median Beta by patient
############################################################

p2 <- ggplot(
  median_beta_long,
  aes(
    x = sample_type,
    y = median_beta,
    group = patient_id
  )
) +
  geom_line(
    alpha = 0.5
  ) +
  geom_point(
    aes(
      color = Grade,
      shape = Sex
    ),
    size = 3
  ) +
  labs(
    title = "SOLID matched tumor–plasma median methylation",
    subtitle = "Each line represents one matched patient",
    x = NULL,
    y = expression("Median " * beta),
    color = "Grade",
    shape = "Sex"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )


ggsave(
  file.path(
    figure_dir,
    "SOLID_patient_paired_median_beta.png"
  ),
  p2,
  width = 7,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(
    figure_dir,
    "SOLID_patient_paired_median_beta.pdf"
  ),
  p2,
  width = 7,
  height = 5.5
)


############################################################
## 14. FIGURE 3
## Valid regions per patient
############################################################

patient_summary$patient_id_factor <- factor(
  patient_summary$patient_id,
  levels = rev(
    patient_summary$patient_id
  )
)


p3 <- ggplot(
  patient_summary,
  aes(
    x = patient_id_factor,
    y = pct_valid_regions,
    fill = factor(Grade)
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "SOLID valid paired regions by patient",
    subtitle = "Eligibility based on the frozen paired-analysis criteria",
    x = "Patient",
    y = "Valid retained regions (%)",
    fill = "Grade"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )


ggsave(
  file.path(
    figure_dir,
    "SOLID_valid_region_percentage_by_patient.png"
  ),
  p3,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(
    figure_dir,
    "SOLID_valid_region_percentage_by_patient.pdf"
  ),
  p3,
  width = 7,
  height = 6
)


############################################################
## 15. FIGURE 4
## Median ΔBeta by patient
##
## ΔBeta = Plasma - Tumor
############################################################

p4 <- ggplot(
  patient_summary,
  aes(
    x = patient_id_factor,
    y = median_delta_beta,
    fill = factor(Grade)
  )
) +
  geom_col() +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  coord_flip() +
  labs(
    title = expression(
      paste(
        "SOLID median ",
        Delta,
        beta,
        " by patient"
      )
    ),
    subtitle = expression(
      Delta * beta == "Plasma" - "Tumor"
    ),
    x = "Patient",
    y = expression(
      "Median " * Delta * beta
    ),
    fill = "Grade"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )


ggsave(
  file.path(
    figure_dir,
    "SOLID_patient_median_delta_beta.png"
  ),
  p4,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(
    figure_dir,
    "SOLID_patient_median_delta_beta.pdf"
  ),
  p4,
  width = 7,
  height = 6
)


############################################################
## 16. PREPARE SAMPLED ΔBeta DATA
############################################################

delta_plot_list <- vector(
  "list",
  length(patient_ids)
)


for (j in seq_along(patient_ids)) {

  x <- delta_beta_valid[
    plot_region_index,
    j
  ]

  delta_plot_list[[j]] <- data.frame(

    patient_id =
      patient_ids[j],

    delta_beta =
      x,

    Grade =
      clinical$Grade[j],

    Sex =
      clinical$Sex[j],

    stringsAsFactors = FALSE
  )
}


delta_plot_data <- do.call(
  rbind,
  delta_plot_list
)


delta_plot_data <- delta_plot_data[
  is.finite(
    delta_plot_data$delta_beta
  ),
  ,
  drop = FALSE
]


delta_plot_data$patient_id <- factor(
  delta_plot_data$patient_id,
  levels = patient_ids
)


delta_plot_data$Grade <- factor(
  delta_plot_data$Grade
)


############################################################
## 17. FIGURE 5
## Per-patient ΔBeta distribution
############################################################

p5 <- ggplot(
  delta_plot_data,
  aes(
    x = patient_id,
    y = delta_beta,
    fill = Grade
  )
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  coord_flip() +
  labs(
    title = expression(
      paste(
        "SOLID ",
        Delta,
        beta,
        " distributions by patient"
      )
    ),
    subtitle = paste0(
      "Plasma - Tumor; representative subset of ",
      format(
        n_plot_regions,
        big.mark = ","
      ),
      " regions"
    ),
    x = "Patient",
    y = expression(
      Delta * beta
    ),
    fill = "Grade"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )


ggsave(
  file.path(
    figure_dir,
    "SOLID_patient_delta_beta_distribution.png"
  ),
  p5,
  width = 7.5,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(
    figure_dir,
    "SOLID_patient_delta_beta_distribution.pdf"
  ),
  p5,
  width = 7.5,
  height = 7
)


############################################################
## 18. FIGURE 6
## Age vs median ΔBeta
##
## Exploratory clinical visualization only.
############################################################

p6 <- ggplot(
  patient_summary,
  aes(
    x = Age,
    y = median_delta_beta,
    color = factor(Grade),
    shape = Sex
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  geom_point(
    size = 3
  ) +
  geom_text(
    aes(
      label = patient_id
    ),
    nudge_y = 0.005,
    check_overlap = TRUE,
    size = 3
  ) +
  labs(
    title = expression(
      paste(
        "Age and median SOLID ",
        Delta,
        beta
      )
    ),
    subtitle = "Exploratory visualization; no association test performed",
    x = "Age",
    y = expression(
      "Median " * Delta * beta
    ),
    color = "Grade",
    shape = "Sex"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )


ggsave(
  file.path(
    figure_dir,
    "SOLID_age_vs_median_delta_beta.png"
  ),
  p6,
  width = 7,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(
    figure_dir,
    "SOLID_age_vs_median_delta_beta.pdf"
  ),
  p6,
  width = 7,
  height = 5.5
)


############################################################
## 19. GLOBAL ΔBETA SUMMARY
############################################################

all_delta <- as.vector(
  delta_beta_valid
)

all_delta <- all_delta[
  is.finite(
    all_delta
  )
]


delta_global_summary <- data.frame(

  statistic = c(
    "n_valid_pair_measurements",
    "mean_delta_beta",
    "median_delta_beta",
    "q25_delta_beta",
    "q75_delta_beta",
    "pct_delta_beta_positive",
    "pct_delta_beta_negative"
  ),

  value = c(

    length(
      all_delta
    ),

    mean(
      all_delta
    ),

    median(
      all_delta
    ),

    quantile(
      all_delta,
      0.25,
      names = FALSE
    ),

    quantile(
      all_delta,
      0.75,
      names = FALSE
    ),

    100 *
      mean(
        all_delta > 0
      ),

    100 *
      mean(
        all_delta < 0
      )
  ),

  stringsAsFactors = FALSE
)


delta_summary_file <- file.path(
  table_dir,
  "SOLID_global_delta_beta_summary.tsv"
)

write.table(
  delta_global_summary,
  delta_summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
## 20. SAVE PLOTTING REGION IDS
##
## Allows exact reproduction of sampled density/boxplot
## figures.
############################################################

plot_regions <- obj$region_annotation[
  plot_region_index,
  ,
  drop = FALSE
]


plot_regions_file <- file.path(
  table_dir,
  "SOLID_exploratory_plot_region_subset.tsv"
)

write.table(
  plot_regions,
  plot_regions_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
## 21. FINAL CONSOLE SUMMARY
############################################################

cat(
  "\n============================================\n"
)

cat(
  "SOLID EXPLORATORY ANALYSIS COMPLETE\n"
)

cat(
  "============================================\n"
)

cat(
  "\nPatients:",
  length(patient_ids),
  "\n"
)

cat(
  "Retained regions:",
  nrow(tissue_beta),
  "\n"
)

cat(
  "Regions sampled for visualization:",
  n_plot_regions,
  "\n"
)

cat(
  "\nValid-region percentage summary:\n"
)

print(
  summary(
    patient_summary$pct_valid_regions
  )
)

cat(
  "\nMedian tumor Beta summary:\n"
)

print(
  summary(
    patient_summary$median_tumor_beta
  )
)

cat(
  "\nMedian plasma Beta summary:\n"
)

print(
  summary(
    patient_summary$median_plasma_beta
  )
)

cat(
  "\nMedian patient ΔBeta summary:\n"
)

print(
  summary(
    patient_summary$median_delta_beta
  )
)

cat(
  "\nGlobal ΔBeta summary:\n"
)

print(
  delta_global_summary
)

cat(
  "\nPatient summary:\n",
  patient_summary_file,
  "\n",
  sep = ""
)

cat(
  "\nFigures:\n",
  figure_dir,
  "\n",
  sep = ""
)

cat(
  "\nTables:\n",
  table_dir,
  "\n",
  sep = ""
)


############################################################
## 22. SESSION INFORMATION
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    exploratory_dir,
    "sessionInfo_04_SOLID_exploratory_visualization.txt"
  )
)
