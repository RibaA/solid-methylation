############################################################
## 04_SOLID_exploratory_visualization.r
##
## SOLID TUMOR-PLASMA EXPLORATORY ANALYSIS
##
## Purpose:
##   1. Load frozen downstream analysis-ready object
##   2. Apply existing valid-pair mask
##   3. Summarize tumor/plasma Beta values by patient
##   4. Summarize Delta Beta = Plasma - Tumor
##   5. Review valid-region availability
##   6. Calculate matched tumor-plasma correlations
##   7. Generate descriptive figures
##   8. Retain selected clinical annotations only:
##        Grade
##        Sex
##        ECOG
##        Response_RANO
##
## IMPORTANT:
##   - Run from repository root
##   - No new filtering
##   - No PCA
##   - No clustering
##   - No DMR testing
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
## 1. REQUIRED PACKAGES
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
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
})


############################################################
## 2. SHARED PLOTTING STYLE
############################################################

plot_style_file <- file.path(
  "script",
  "00_plot_style_and_palettes.r"
)

if (!file.exists(plot_style_file)) {
  stop(
    "Shared plotting style file not found:\n",
    plot_style_file
  )
}

source(plot_style_file)


############################################################
## 3. INPUT / OUTPUT
############################################################

input_file <- file.path(
  "result",
  "03_matched_tissue_plasma",
  "SOLID_downstream_analysis_ready.rds"
)

exploratory_dir <- file.path(
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

if (!file.exists(input_file)) {
  stop(
    "Input object not found:\n",
    input_file
  )
}


############################################################
## 4. HELPERS
############################################################

safe_mean <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  mean(x)
}


safe_median <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  median(x)
}


safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  quantile(
    x,
    probs = probability,
    names = FALSE
  )
}


safe_cor <- function(
    x,
    y,
    method = "pearson"
) {
  keep <- is.finite(x) &
    is.finite(y)

  if (sum(keep) < 3L) {
    return(NA_real_)
  }

  suppressWarnings(
    cor(
      x[keep],
      y[keep],
      method = method
    )
  )
}


save_plot <- function(
    plot_object,
    filename_base,
    width,
    height
) {
  ggsave(
    filename = file.path(
      figure_dir,
      paste0(
        filename_base,
        ".png"
      )
    ),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300
  )

  ggsave(
    filename = file.path(
      figure_dir,
      paste0(
        filename_base,
        ".pdf"
      )
    ),
    plot = plot_object,
    width = width,
    height = height
  )
}


############################################################
## 5. LOAD OBJECT
############################################################

obj <- readRDS(
  input_file
)

cat(
  "\n============================================\n",
  "SOLID EXPLORATORY VISUALIZATION\n",
  "============================================\n",
  sep = ""
)


############################################################
## 6. BASIC VALIDATION
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
  "\nInput validation: PASS\n",
  "Patients: ",
  ncol(tissue_beta),
  "\nRegions: ",
  nrow(tissue_beta),
  "\n",
  sep = ""
)


############################################################
## 7. APPLY FROZEN VALID-PAIR MASK
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
## 8. PATIENT-LEVEL SUMMARY
############################################################

patient_summary <- data.frame(
  patient_id = patient_ids,
  stringsAsFactors = FALSE
)


############################################################
## Valid-region support
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
## Tumor summaries
############################################################

patient_summary$mean_tumor_beta <- vapply(
  seq_along(patient_ids),
  function(j) {
    safe_mean(
      tissue_beta_valid[, j]
    )
  },
  numeric(1)
)

patient_summary$median_tumor_beta <- vapply(
  seq_along(patient_ids),
  function(j) {
    safe_median(
      tissue_beta_valid[, j]
    )
  },
  numeric(1)
)


############################################################
## Plasma summaries
############################################################

patient_summary$mean_plasma_beta <- vapply(
  seq_along(patient_ids),
  function(j) {
    safe_mean(
      plasma_beta_valid[, j]
    )
  },
  numeric(1)
)

patient_summary$median_plasma_beta <- vapply(
  seq_along(patient_ids),
  function(j) {
    safe_median(
      plasma_beta_valid[, j]
    )
  },
  numeric(1)
)


############################################################
## Delta-Beta summaries
############################################################

patient_summary$mean_delta_beta <- vapply(
  seq_along(patient_ids),
  function(j) {
    safe_mean(
      delta_beta_valid[, j]
    )
  },
  numeric(1)
)

patient_summary$median_delta_beta <- vapply(
  seq_along(patient_ids),
  function(j) {
    safe_median(
      delta_beta_valid[, j]
    )
  },
  numeric(1)
)

patient_summary$q25_delta_beta <- vapply(
  seq_along(patient_ids),
  function(j) {
    safe_quantile(
      delta_beta_valid[, j],
      0.25
    )
  },
  numeric(1)
)

patient_summary$q75_delta_beta <- vapply(
  seq_along(patient_ids),
  function(j) {
    safe_quantile(
      delta_beta_valid[, j],
      0.75
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
    x <- x[is.finite(x)]

    if (length(x) == 0L) {
      return(NA_real_)
    }

    100 * mean(x > 0)
  },
  numeric(1)
)


patient_summary$pct_tumor_higher <- vapply(
  seq_along(patient_ids),
  function(j) {
    x <- delta_beta_valid[, j]
    x <- x[is.finite(x)]

    if (length(x) == 0L) {
      return(NA_real_)
    }

    100 * mean(x < 0)
  },
  numeric(1)
)


############################################################
## Matched correlations
############################################################

patient_summary$tumor_plasma_pearson <- vapply(
  seq_along(patient_ids),
  function(j) {
    safe_cor(
      tissue_beta_valid[, j],
      plasma_beta_valid[, j],
      method = "pearson"
    )
  },
  numeric(1)
)


patient_summary$tumor_plasma_spearman <- vapply(
  seq_along(patient_ids),
  function(j) {
    safe_cor(
      tissue_beta_valid[, j],
      plasma_beta_valid[, j],
      method = "spearman"
    )
  },
  numeric(1)
)


############################################################
## 9. SELECT CLINICAL ANNOTATIONS
############################################################

annotation_variables <- c(
  "Grade",
  "Sex",
  "ECOG",
  "Response_RANO"
)


missing_annotation_variables <- setdiff(
  annotation_variables,
  names(clinical)
)


if (length(missing_annotation_variables) > 0L) {
  stop(
    "Missing clinical annotation(s): ",
    paste(
      missing_annotation_variables,
      collapse = ", "
    )
  )
}


clinical_for_merge <- clinical[
  ,
  c(
    "patient_id",
    annotation_variables
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


for (v in annotation_variables) {
  patient_summary[[v]] <- factor(
    patient_summary[[v]]
  )
}


############################################################
## 10. SAVE PATIENT SUMMARY
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
## 11. SAMPLE-TYPE DISTRIBUTION SUMMARY
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
      na.rm = TRUE,
      names = FALSE
    ),
    q75_beta = quantile(
      tissue_beta_valid,
      0.75,
      na.rm = TRUE,
      names = FALSE
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
      na.rm = TRUE,
      names = FALSE
    ),
    q75_beta = quantile(
      plasma_beta_valid,
      0.75,
      na.rm = TRUE,
      names = FALSE
    )
  )
)


write.table(
  distribution_summary,
  file.path(
    table_dir,
    "SOLID_tumor_plasma_beta_distribution_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
## 12. REPRODUCIBLE REGION SUBSET
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


############################################################
## 13. FIGURE 1
## Tumor vs Plasma Beta density
############################################################

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


p1 <- ggplot(
  beta_density_data,
  aes(
    x = beta,
    color = sample_type,
    fill = sample_type
  )
) +

  geom_density(
    alpha = 0.18,
    linewidth = 0.9
  ) +

  scale_color_manual(
    values = sample_cols[
      c(
        "Tumor",
        "Plasma"
      )
    ]
  ) +

  scale_fill_manual(
    values = sample_cols[
      c(
        "Tumor",
        "Plasma"
      )
    ]
  ) +

  labs(
    title =
      "SOLID tumor and plasma methylation distributions",
    subtitle =
      paste0(
        "Representative subset of ",
        format(
          n_plot_regions,
          big.mark = ","
        ),
        " retained 1-kb regions"
      ),
    x = "Beta value",
    y = "Density",
    color = "Sample type",
    fill = "Sample type"
  ) +

  theme_project()


save_plot(
  p1,
  "SOLID_tumor_plasma_beta_density",
  width = 7,
  height = 5
)


############################################################
## 14. PAIRED MEDIAN-BETA DATA
############################################################

median_beta_long <- rbind(

  data.frame(
    patient_id =
      patient_summary$patient_id,
    sample_type =
      "Tumor",
    median_beta =
      patient_summary$median_tumor_beta,
    Sex =
      patient_summary$Sex
  ),

  data.frame(
    patient_id =
      patient_summary$patient_id,
    sample_type =
      "Plasma",
    median_beta =
      patient_summary$median_plasma_beta,
    Sex =
      patient_summary$Sex
  )
)


median_beta_long$sample_type <- factor(
  median_beta_long$sample_type,
  levels = c(
    "Tumor",
    "Plasma"
  )
)


############################################################
## 15. FIGURE 2
## Paired median Beta
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
    color = "grey70",
    alpha = 0.65,
    linewidth = 0.65
  ) +

  geom_point(
    aes(
      color = sample_type,
      shape = Sex
    ),
    size = 3.2
  ) +

  scale_color_manual(
    values = sample_cols[
      c(
        "Tumor",
        "Plasma"
      )
    ]
  ) +

  labs(
    title =
      "SOLID matched tumor-plasma median methylation",
    subtitle =
      "Each line represents one matched patient",
    x = NULL,
    y = "Median Beta",
    color = "Sample type",
    shape = "Sex"
  ) +

  theme_project()


save_plot(
  p2,
  "SOLID_patient_paired_median_beta",
  width = 7,
  height = 5.5
)


############################################################
## 16. FIGURE 3
## Valid paired regions
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
    fill = Grade
  )
) +

  geom_col(
    width = 0.72
  ) +

  coord_flip() +

  scale_fill_manual(
    values = grade_cols
  ) +

  labs(
    title =
      "SOLID valid paired regions by patient",
    subtitle =
      "Eligibility based on frozen paired-analysis criteria",
    x = "Patient",
    y = "Valid retained regions (%)",
    fill = "Grade"
  ) +

  theme_project()


save_plot(
  p3,
  "SOLID_valid_region_percentage_by_patient",
  width = 7,
  height = 6
)


############################################################
## 17. FIGURE 4
## Median Delta Beta
############################################################

p4 <- ggplot(
  patient_summary,
  aes(
    x = patient_id_factor,
    y = median_delta_beta,
    fill = Grade
  )
) +

  geom_col(
    width = 0.72
  ) +

  geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.5,
    color = "grey45"
  ) +

  coord_flip() +

  scale_fill_manual(
    values = grade_cols
  ) +

  labs(
    title =
      "SOLID median delta Beta by patient",
    subtitle =
      "Delta Beta = Plasma - Tumor",
    x = "Patient",
    y = "Median delta Beta",
    fill = "Grade"
  ) +

  theme_project()


save_plot(
  p4,
  "SOLID_patient_median_delta_beta",
  width = 7,
  height = 6
)


############################################################
## 18. SAMPLED DELTA-BETA DATA
############################################################

delta_plot_list <- vector(
  "list",
  length(patient_ids)
)


for (j in seq_along(patient_ids)) {

  delta_plot_list[[j]] <- data.frame(

    patient_id =
      patient_ids[j],

    delta_beta =
      delta_beta_valid[
        plot_region_index,
        j
      ],

    Grade =
      patient_summary$Grade[j]
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
## 19. FIGURE 5
## Per-patient Delta-Beta distributions
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
    outlier.shape = NA,
    linewidth = 0.45
  ) +

  geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.5,
    color = "grey45"
  ) +

  coord_flip() +

  scale_fill_manual(
    values = grade_cols
  ) +

  labs(
    title =
      "SOLID delta Beta distributions by patient",
    subtitle =
      paste0(
        "Plasma - Tumor; representative subset of ",
        format(
          n_plot_regions,
          big.mark = ","
        ),
        " regions"
      ),
    x = "Patient",
    y = "Delta Beta",
    fill = "Grade"
  ) +

  theme_project()


save_plot(
  p5,
  "SOLID_patient_delta_beta_distribution",
  width = 7.5,
  height = 7
)


############################################################
## 20. FIGURE 6
## Patient-level tumor-plasma correlation
############################################################

correlation_plot_data <- patient_summary


correlation_plot_data$patient_id_cor <- factor(
  correlation_plot_data$patient_id,
  levels =
    correlation_plot_data$patient_id[
      order(
        correlation_plot_data$tumor_plasma_pearson
      )
    ]
)


p6 <- ggplot(
  correlation_plot_data,
  aes(
    x = tumor_plasma_pearson,
    y = patient_id_cor,
    color = Grade
  )
) +

  geom_point(
    size = 3.2
  ) +

  scale_color_manual(
    values = grade_cols
  ) +

  labs(
    title =
      "SOLID tumor-plasma methylation correlation",
    subtitle =
      "Matched patient-level Pearson correlation",
    x =
      "Pearson correlation",
    y =
      "Patient",
    color =
      "Grade"
  ) +

  theme_project()


save_plot(
  p6,
  "SOLID_patient_tumor_plasma_correlation",
  width = 7,
  height = 6
)


############################################################
## 21. GLOBAL DELTA-BETA SUMMARY
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
  )
)


write.table(
  delta_global_summary,
  file.path(
    table_dir,
    "SOLID_global_delta_beta_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
## 22. CORRELATION SUMMARY
############################################################

correlation_summary <- data.frame(

  statistic = c(
    "minimum",
    "Q1",
    "median",
    "mean",
    "Q3",
    "maximum"
  ),

  pearson = c(
    min(
      patient_summary$tumor_plasma_pearson,
      na.rm = TRUE
    ),
    quantile(
      patient_summary$tumor_plasma_pearson,
      0.25,
      na.rm = TRUE,
      names = FALSE
    ),
    median(
      patient_summary$tumor_plasma_pearson,
      na.rm = TRUE
    ),
    mean(
      patient_summary$tumor_plasma_pearson,
      na.rm = TRUE
    ),
    quantile(
      patient_summary$tumor_plasma_pearson,
      0.75,
      na.rm = TRUE,
      names = FALSE
    ),
    max(
      patient_summary$tumor_plasma_pearson,
      na.rm = TRUE
    )
  ),

  spearman = c(
    min(
      patient_summary$tumor_plasma_spearman,
      na.rm = TRUE
    ),
    quantile(
      patient_summary$tumor_plasma_spearman,
      0.25,
      na.rm = TRUE,
      names = FALSE
    ),
    median(
      patient_summary$tumor_plasma_spearman,
      na.rm = TRUE
    ),
    mean(
      patient_summary$tumor_plasma_spearman,
      na.rm = TRUE
    ),
    quantile(
      patient_summary$tumor_plasma_spearman,
      0.75,
      na.rm = TRUE,
      names = FALSE
    ),
    max(
      patient_summary$tumor_plasma_spearman,
      na.rm = TRUE
    )
  )
)


write.table(
  correlation_summary,
  file.path(
    table_dir,
    "SOLID_tumor_plasma_correlation_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
## 23. SAVE PLOTTING REGION IDS
############################################################

plot_regions <- obj$region_annotation[
  plot_region_index,
  ,
  drop = FALSE
]


write.table(
  plot_regions,
  file.path(
    table_dir,
    "SOLID_exploratory_plot_region_subset.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
## 24. FINAL SUMMARY
############################################################

cat(
  "\n============================================\n",
  "SOLID EXPLORATORY ANALYSIS COMPLETE\n",
  "============================================\n",
  sep = ""
)


cat(
  "\nPatients: ",
  length(patient_ids),
  "\n",
  sep = ""
)


cat(
  "Retained regions: ",
  nrow(tissue_beta),
  "\n",
  sep = ""
)


cat(
  "Regions sampled for visualization: ",
  n_plot_regions,
  "\n",
  sep = ""
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
  "\nMedian patient Delta Beta summary:\n"
)

print(
  summary(
    patient_summary$median_delta_beta
  )
)


cat(
  "\nTumor-plasma correlation summary:\n"
)

print(
  correlation_summary
)


cat(
  "\nGlobal Delta Beta summary:\n"
)

print(
  delta_global_summary
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
## 25. SESSION INFORMATION
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    exploratory_dir,
    "sessionInfo_04_SOLID_exploratory_visualization.txt"
  )
)