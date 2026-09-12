############################################################
## 06_SOLID_matched_concordance.r
##
## SOLID MATCHED TUMOR–PLASMA CONCORDANCE
##
## Main questions:
##
##   1. How strongly does each tumor correlate with its
##      matched plasma sample?
##
##   2. Is the matched tumor-plasma correlation stronger
##      than correlations with unmatched plasma samples?
##
##   3. Which 1-kb regions show recurrent tumor-plasma
##      concordance or discordance across patients?
##
## Primary scale:
##   Beta values
##
## Delta Beta:
##   Plasma - Tumor
##
## IMPORTANT:
##   - Uses frozen 124,961-region paired dataset
##   - No new QC filtering
##   - No DMR testing
##   - No final candidate cutoff imposed here
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
  "ggplot2",
  "pheatmap"
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
    "\nInstall with:\ninstall.packages(c(",
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
  library(pheatmap)
})

############################################################
## 2. DIRECTORIES
############################################################

project_dir <- "C:/solid-methylation"

input_file <- file.path(
  project_dir,
  "result",
  "03_matched_tissue_plasma",
  "SOLID_downstream_analysis_ready.rds"
)

output_dir <- file.path(
  project_dir,
  "result",
  "03_matched_tissue_plasma",
  "concordance"
)

figure_dir <- file.path(
  output_dir,
  "figures"
)

table_dir <- file.path(
  output_dir,
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
  "SOLID MATCHED CONCORDANCE ANALYSIS\n"
)

cat(
  "============================================\n"
)

############################################################
## 4. LOAD COMPONENTS
############################################################

tissue_beta <- obj$tissue_beta
plasma_beta <- obj$plasma_beta
delta_beta <- obj$delta_beta
valid_mask <- obj$valid_pair_mask

annotation <- as.data.frame(
  obj$region_annotation
)

clinical <- as.data.frame(
  obj$patient_metadata
)

############################################################
## 5. BASIC VALIDATION
############################################################

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

  nrow(annotation) == 124961L,

  nrow(clinical) == 13L
)

patient_ids <- colnames(
  tissue_beta
)

stopifnot(
  identical(
    patient_ids,
    colnames(plasma_beta)
  ),

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
  length(patient_ids),
  "\n"
)

cat(
  "Regions:",
  nrow(tissue_beta),
  "\n"
)

############################################################
## 6. APPLY FROZEN VALID-PAIR MASK
############################################################

tissue_valid <- tissue_beta
plasma_valid <- plasma_beta
delta_valid <- delta_beta

tissue_valid[
  !valid_mask
] <- NA_real_

plasma_valid[
  !valid_mask
] <- NA_real_

delta_valid[
  !valid_mask
] <- NA_real_

############################################################
## 7. PRIMARY MATCHED-PAIR CORRELATIONS
##
## Each patient uses all regions valid for that patient's
## matched tumor-plasma pair.
############################################################

matched_rho <- numeric(
  length(patient_ids)
)

matched_n_regions <- integer(
  length(patient_ids)
)

for (j in seq_along(patient_ids)) {

  good <- is.finite(
    tissue_valid[, j]
  ) &
    is.finite(
      plasma_valid[, j]
    )

  matched_n_regions[j] <- sum(
    good
  )

  matched_rho[j] <- cor(
    tissue_valid[good, j],
    plasma_valid[good, j],
    method = "spearman"
  )
}

matched_summary <- data.frame(
  patient_id = patient_ids,
  n_valid_regions = matched_n_regions,
  matched_spearman_rho = matched_rho,
  stringsAsFactors = FALSE
)

############################################################
## 8. ADD CLINICAL VARIABLES
############################################################

matched_summary <- cbind(
  matched_summary,
  clinical[
    ,
    setdiff(
      names(clinical),
      "patient_id"
    ),
    drop = FALSE
  ]
)

############################################################
## 9. COMMON COMPLETE REGION SET FOR FAIR
##    MATCHED-vs-UNMATCHED COMPARISON
##
## Require a region to have a valid paired measurement
## in every one of the 13 patients.
##
## This gives the same features for every cross-patient
## tumor-plasma correlation.
############################################################

common_complete <- apply(
  valid_mask,
  1,
  all
)

n_common_complete <- sum(
  common_complete
)

cat(
  "\nRegions valid in all 13 matched pairs:",
  n_common_complete,
  "\n"
)

cat(
  "Percentage of retained regions:",
  round(
    100 *
      n_common_complete /
      nrow(valid_mask),
    3
  ),
  "%\n"
)

if (n_common_complete < 1000L) {

  stop(
    "Too few common complete regions for matched/unmatched comparison."
  )
}

tissue_common <- tissue_beta[
  common_complete,
  ,
  drop = FALSE
]

plasma_common <- plasma_beta[
  common_complete,
  ,
  drop = FALSE
]

stopifnot(
  all(
    is.finite(
      tissue_common
    )
  ),

  all(
    is.finite(
      plasma_common
    )
  )
)

############################################################
## 10. ALL 13 x 13 TUMOR-PLASMA CORRELATIONS
############################################################

cor_matrix <- matrix(
  NA_real_,
  nrow = length(patient_ids),
  ncol = length(patient_ids),
  dimnames = list(
    paste0(
      "Tumor_",
      patient_ids
    ),
    paste0(
      "Plasma_",
      patient_ids
    )
  )
)

for (i in seq_along(patient_ids)) {

  for (j in seq_along(patient_ids)) {

    cor_matrix[i, j] <- cor(
      tissue_common[, i],
      plasma_common[, j],
      method = "spearman"
    )
  }
}

############################################################
## 11. MATCHED AND UNMATCHED CORRELATION SUMMARY
############################################################

matched_common_rho <- diag(
  cor_matrix
)

median_unmatched_rho <- numeric(
  length(patient_ids)
)

mean_unmatched_rho <- numeric(
  length(patient_ids)
)

matched_rank <- integer(
  length(patient_ids)
)

matched_percentile <- numeric(
  length(patient_ids)
)

for (i in seq_along(patient_ids)) {

  all_rho_i <- cor_matrix[
    i,
  ]

  unmatched_i <- all_rho_i[
    -i
  ]

  median_unmatched_rho[i] <- median(
    unmatched_i
  )

  mean_unmatched_rho[i] <- mean(
    unmatched_i
  )


  ##########################################################
  ## Rank matched plasma among all 13 plasma samples
  ## rank 1 = highest tumor-plasma correlation
  ##########################################################

  matched_rank[i] <- rank(
    -all_rho_i,
    ties.method = "min"
  )[i]


  ##########################################################
  ## Percent of unmatched correlations lower than matched
  ##########################################################

  matched_percentile[i] <-
    100 *
    mean(
      unmatched_i <
        all_rho_i[i]
    )
}

matched_summary$matched_common_rho <-
  matched_common_rho

matched_summary$median_unmatched_rho <-
  median_unmatched_rho

matched_summary$mean_unmatched_rho <-
  mean_unmatched_rho

matched_summary$matched_minus_median_unmatched <-
  matched_common_rho -
  median_unmatched_rho

matched_summary$matched_rank_among_13_plasma <-
  matched_rank

matched_summary$matched_percentile_vs_unmatched <-
  matched_percentile

############################################################
## 12. PATIENT-LEVEL MATCHED vs UNMATCHED TEST
##
## Each patient contributes:
##
##   matched correlation
##   median of its 12 unmatched correlations
##
## Paired Wilcoxon test avoids treating all 156 unmatched
## correlations as independent observations.
############################################################

matched_unmatched_test <- wilcox.test(
  matched_common_rho,
  median_unmatched_rho,
  paired = TRUE,
  exact = FALSE
)

test_summary <- data.frame(
  test = "paired Wilcoxon: matched rho vs median unmatched rho",
  n_patients = length(patient_ids),
  median_matched_rho = median(
    matched_common_rho
  ),
  median_unmatched_rho = median(
    median_unmatched_rho
  ),
  median_difference = median(
    matched_common_rho -
      median_unmatched_rho
  ),
  statistic = unname(
    matched_unmatched_test$statistic
  ),
  p_value = matched_unmatched_test$p.value,
  stringsAsFactors = FALSE
)

############################################################
## 13. SAVE CORRELATION TABLES
############################################################

write.table(
  matched_summary,
  file.path(
    table_dir,
    "SOLID_matched_tumor_plasma_correlation_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  cor_matrix,
  file.path(
    table_dir,
    "SOLID_tumor_plasma_all_pair_correlations.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)


write.table(
  test_summary,
  file.path(
    table_dir,
    "SOLID_matched_vs_unmatched_correlation_test.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
## 14. FIGURE 1
## Matched Spearman correlation by patient
############################################################

matched_summary$patient_factor <- factor(
  matched_summary$patient_id,
  levels = rev(
    matched_summary$patient_id
  )
)


p1 <- ggplot(
  matched_summary,
  aes(
    x = patient_factor,
    y = matched_spearman_rho,
    fill = factor(Grade)
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "SOLID matched tumor–plasma concordance",
    subtitle = "Spearman correlation across eligible 1-kb regions",
    x = "Patient",
    y = "Spearman correlation",
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
    "SOLID_matched_tumor_plasma_spearman_by_patient.png"
  ),
  p1,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(
    figure_dir,
    "SOLID_matched_tumor_plasma_spearman_by_patient.pdf"
  ),
  p1,
  width = 7,
  height = 6
)


############################################################
## 15. FIGURE 2
## Matched vs median unmatched correlation
############################################################

comparison_long <- rbind(

  data.frame(
    patient_id = patient_ids,
    comparison = "Matched",
    rho = matched_common_rho,
    stringsAsFactors = FALSE
  ),

  data.frame(
    patient_id = patient_ids,
    comparison = "Median unmatched",
    rho = median_unmatched_rho,
    stringsAsFactors = FALSE
  )
)


comparison_long$comparison <- factor(
  comparison_long$comparison,
  levels = c(
    "Median unmatched",
    "Matched"
  )
)


p2 <- ggplot(
  comparison_long,
  aes(
    x = comparison,
    y = rho,
    group = patient_id
  )
) +
  geom_line(
    alpha = 0.5
  ) +
  geom_point(
    size = 3
  ) +
  labs(
    title = "SOLID matched vs unmatched tumor–plasma correlation",
    subtitle = paste0(
      "Common complete region set; paired Wilcoxon P = ",
      signif(
        matched_unmatched_test$p.value,
        3
      )
    ),
    x = NULL,
    y = "Spearman correlation"
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
    "SOLID_matched_vs_unmatched_correlation.png"
  ),
  p2,
  width = 7,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(
    figure_dir,
    "SOLID_matched_vs_unmatched_correlation.pdf"
  ),
  p2,
  width = 7,
  height = 5.5
)


############################################################
## 16. FIGURE 3
## 13 x 13 tumor-plasma correlation heatmap
##
## Diagonal = true matched pairs
############################################################

pdf(
  file.path(
    figure_dir,
    "SOLID_tumor_plasma_all_pair_correlation_heatmap.pdf"
  ),
  width = 9,
  height = 8
)

pheatmap(
  cor_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = FALSE,
  main =
    "SOLID tumor–plasma Spearman correlations"
)

dev.off()


png(
  file.path(
    figure_dir,
    "SOLID_tumor_plasma_all_pair_correlation_heatmap.png"
  ),
  width = 2600,
  height = 2400,
  res = 300
)

pheatmap(
  cor_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = FALSE,
  main =
    "SOLID tumor–plasma Spearman correlations"
)

dev.off()


############################################################
## 17. REGION-LEVEL CONCORDANCE SUMMARY
##
## No final threshold is imposed here.
##
## We calculate several descriptive measures so that
## concordance cutoffs can be chosen transparently later.
############################################################

region_summary <- annotation


region_summary$mean_delta_beta <- apply(
  delta_valid,
  1,
  mean,
  na.rm = TRUE
)


region_summary$median_delta_beta <- apply(
  delta_valid,
  1,
  median,
  na.rm = TRUE
)


region_summary$mean_abs_delta_beta <- apply(
  abs(
    delta_valid
  ),
  1,
  mean,
  na.rm = TRUE
)


region_summary$median_abs_delta_beta <- apply(
  abs(
    delta_valid
  ),
  1,
  median,
  na.rm = TRUE
)


############################################################
## Percentage of valid patients with small |Delta Beta|
############################################################

pct_within_threshold <- function(
    x,
    threshold
) {

  x <- x[
    is.finite(x)
  ]

  if (length(x) == 0L) {

    return(
      NA_real_
    )
  }

  100 *
    mean(
      abs(x) <= threshold
    )
}


region_summary$pct_abs_delta_le_0.05 <- apply(
  delta_valid,
  1,
  pct_within_threshold,
  threshold = 0.05
)


region_summary$pct_abs_delta_le_0.10 <- apply(
  delta_valid,
  1,
  pct_within_threshold,
  threshold = 0.10
)


############################################################
## Recurrent directional differences
############################################################

region_summary$pct_plasma_higher_0.10 <- apply(
  delta_valid,
  1,
  function(x) {

    x <- x[
      is.finite(x)
    ]

    100 *
      mean(
        x >= 0.10
      )
  }
)


region_summary$pct_tumor_higher_0.10 <- apply(
  delta_valid,
  1,
  function(x) {

    x <- x[
      is.finite(x)
    ]

    100 *
      mean(
        x <= -0.10
      )
  }
)


region_summary$pct_plasma_higher_0.20 <- apply(
  delta_valid,
  1,
  function(x) {

    x <- x[
      is.finite(x)
    ]

    100 *
      mean(
        x >= 0.20
      )
  }
)


region_summary$pct_tumor_higher_0.20 <- apply(
  delta_valid,
  1,
  function(x) {

    x <- x[
      is.finite(x)
    ]

    100 *
      mean(
        x <= -0.20
      )
  }
)


############################################################
## 18. SAVE FULL REGION CONCORDANCE TABLE
############################################################

region_summary_file <- file.path(
  table_dir,
  "SOLID_region_tumor_plasma_concordance_summary.tsv.gz"
)


con <- gzfile(
  region_summary_file,
  open = "wt"
)

write.table(
  region_summary,
  con,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

close(
  con
)


############################################################
## 19. TOP CONCORDANT REGIONS
##
## Ranking only.
## This is NOT a final biological candidate definition.
############################################################

concordant_rank <- order(
  region_summary$median_abs_delta_beta,
  region_summary$mean_abs_delta_beta
)


top_concordant <- region_summary[
  concordant_rank[
    seq_len(
      min(
        1000L,
        nrow(region_summary)
      )
    )
  ],
  ,
  drop = FALSE
]


write.table(
  top_concordant,
  file.path(
    table_dir,
    "SOLID_top1000_tumor_plasma_concordant_regions.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
## 20. TOP DIRECTIONAL DISCORDANT REGIONS
##
## Exploratory ranking only.
############################################################

plasma_higher_rank <- order(
  region_summary$median_delta_beta,
  decreasing = TRUE
)

tumor_higher_rank <- order(
  region_summary$median_delta_beta,
  decreasing = FALSE
)


top_plasma_higher <- region_summary[
  plasma_higher_rank[
    seq_len(
      min(
        1000L,
        nrow(region_summary)
      )
    )
  ],
  ,
  drop = FALSE
]


top_tumor_higher <- region_summary[
  tumor_higher_rank[
    seq_len(
      min(
        1000L,
        nrow(region_summary)
      )
    )
  ],
  ,
  drop = FALSE
]


write.table(
  top_plasma_higher,
  file.path(
    table_dir,
    "SOLID_top1000_plasma_higher_regions_exploratory.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  top_tumor_higher,
  file.path(
    table_dir,
    "SOLID_top1000_tumor_higher_regions_exploratory.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
## 21. CONCORDANCE THRESHOLD SENSITIVITY
##
## Count regions meeting different recurrence thresholds.
##
## This helps us decide later what "concordant" should mean.
############################################################

recurrence_levels <- c(
  50,
  70,
  80,
  90,
  100
)


threshold_summary <- do.call(
  rbind,
  lapply(
    recurrence_levels,
    function(pct) {

      data.frame(

        recurrence_required_pct =
          pct,

        n_regions_abs_delta_le_0.05 =
          sum(
            region_summary$pct_abs_delta_le_0.05 >= pct,
            na.rm = TRUE
          ),

        n_regions_abs_delta_le_0.10 =
          sum(
            region_summary$pct_abs_delta_le_0.10 >= pct,
            na.rm = TRUE
          ),

        n_regions_plasma_higher_0.10 =
          sum(
            region_summary$pct_plasma_higher_0.10 >= pct,
            na.rm = TRUE
          ),

        n_regions_tumor_higher_0.10 =
          sum(
            region_summary$pct_tumor_higher_0.10 >= pct,
            na.rm = TRUE
          ),

        stringsAsFactors = FALSE
      )
    }
  )
)


write.table(
  threshold_summary,
  file.path(
    table_dir,
    "SOLID_region_concordance_threshold_sensitivity.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
## 22. FIGURE 4
## Distribution of regional median absolute Delta Beta
############################################################

p4 <- ggplot(
  region_summary,
  aes(
    x = median_abs_delta_beta
  )
) +
  geom_histogram(
    bins = 80
  ) +
  geom_vline(
    xintercept = c(
      0.05,
      0.10
    ),
    linetype = 2
  ) +
  labs(
    title = "SOLID regional tumor–plasma concordance",
    subtitle = "Smaller median |Δβ| indicates greater tumor–plasma similarity",
    x = expression(
      "Median |" * Delta * beta * "|"
    ),
    y = "Number of regions"
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
    "SOLID_region_median_absolute_delta_beta_distribution.png"
  ),
  p4,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(
    figure_dir,
    "SOLID_region_median_absolute_delta_beta_distribution.pdf"
  ),
  p4,
  width = 7,
  height = 5
)


############################################################
## 23. FIGURE 5
## Regional direction/effect summary
############################################################

p5 <- ggplot(
  region_summary,
  aes(
    x = median_delta_beta,
    y = median_abs_delta_beta
  )
) +
  geom_point(
    alpha = 0.15,
    size = 0.6
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 2
  ) +
  labs(
    title = "SOLID regional tumor–plasma effect structure",
    subtitle = expression(
      Delta * beta == "Plasma - Tumor"
    ),
    x = expression(
      "Median " * Delta * beta
    ),
    y = expression(
      "Median |" * Delta * beta * "|"
    )
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
    "SOLID_region_delta_beta_effect_structure.png"
  ),
  p5,
  width = 7,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(
    figure_dir,
    "SOLID_region_delta_beta_effect_structure.pdf"
  ),
  p5,
  width = 7,
  height = 5.5
)


############################################################
## 24. FINAL CONSOLE SUMMARY
############################################################

cat(
  "\n============================================\n"
)

cat(
  "SOLID MATCHED CONCORDANCE COMPLETE\n"
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
  "Common complete regions used for matched/unmatched:",
  n_common_complete,
  "\n"
)


cat(
  "\nMatched Spearman correlations:\n"
)

print(
  summary(
    matched_summary$matched_spearman_rho
  )
)


cat(
  "\nMatched correlations on common complete regions:\n"
)

print(
  summary(
    matched_common_rho
  )
)


cat(
  "\nMedian unmatched correlations:\n"
)

print(
  summary(
    median_unmatched_rho
  )
)


cat(
  "\nMatched - median unmatched difference:\n"
)

print(
  summary(
    matched_common_rho -
      median_unmatched_rho
  )
)


cat(
  "\nMatched pair rank among 13 plasma samples:\n"
)

print(
  table(
    matched_rank
  )
)


cat(
  "\nPaired matched-vs-unmatched test:\n"
)

print(
  test_summary
)


cat(
  "\nRegional concordance threshold sensitivity:\n"
)

print(
  threshold_summary
)


cat(
  "\nMatched correlation table:\n",
  file.path(
    table_dir,
    "SOLID_matched_tumor_plasma_correlation_summary.tsv"
  ),
  "\n",
  sep = ""
)


cat(
  "\nRegional concordance table:\n",
  region_summary_file,
  "\n",
  sep = ""
)


cat(
  "\nFigures:\n",
  figure_dir,
  "\n",
  sep = ""
)


############################################################
## 25. SESSION INFORMATION
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "sessionInfo_06_SOLID_matched_concordance.txt"
  )
)
