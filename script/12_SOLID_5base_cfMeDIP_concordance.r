##------------------------------------------------------------
# SOLID PLASMA 5-BASE <-> cfMeDIP CONCORDANCE
#
# Script:
#   12_SOLID_5base_cfMeDIP_concordance.r
#
# Goal:
#   Compare plasma methylation measured by:
#
#   1. 5-base sequencing
#      - regional Beta values
#
#   2. cfMeDIP
#      - regional methylated-fragment enrichment
#      - aggregated to matched 1-kb regions in Script 11
#
# Cohort:
#   13 matched SOLID patients
#
# Main analyses:
#   1. Confirm identical patients
#   2. Confirm common 1-kb regions
#   3. Patient-wise Spearman correlation
#   4. Cross-patient correlation matrix
#   5. Matched vs unmatched patient comparison
#   6. True-match ranking
#   7. Region-wise Spearman correlation across patients
#
# Important:
#   - 5-base Beta and cfMeDIP log2CPM are different scales.
#   - Spearman correlation is therefore used.
#   - Beta and log2CPM are NOT directly subtracted.
#
# Region IDs:
#   Both assays use:
#
#      chr:start:end
#
# Heatmap:
#   Rows and columns are displayed in the SAME patient order.
#   Therefore true matched samples always lie on the diagonal.
##------------------------------------------------------------


##------------------------------------------------------------
# 0. Setup
##------------------------------------------------------------

rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  scipen = 999
)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})


##------------------------------------------------------------
# 1. Directories
##------------------------------------------------------------

dir_paired <- file.path(
  "result",
  "03_matched_tissue_plasma",
  "paired_filtered"
)

dir_cfmedip_object <- file.path(
  "result",
  "05_plasma_cfMeDIP",
  "harmonized_1kb",
  "objects"
)


##------------------------------------------------------------
# Output directories
##------------------------------------------------------------

dir_result <- file.path(
  "result",
  "06_multiassay_integration",
  "5base_cfMeDIP_concordance"
)

dir_table <- file.path(
  dir_result,
  "tables"
)

dir_fig <- file.path(
  dir_result,
  "figures"
)

dir_object <- file.path(
  dir_result,
  "objects"
)


dir.create(
  dir_table,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dir_fig,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dir_object,
  recursive = TRUE,
  showWarnings = FALSE
)


##------------------------------------------------------------
# 2. Input files
##------------------------------------------------------------

file_5base_beta <- file.path(
  dir_paired,
  "SOLID_paired_plasma_beta_filtered.rds"
)

file_5base_annotation <- file.path(
  dir_paired,
  "SOLID_paired_retained_region_annotation.tsv.gz"
)

file_medip <- file.path(
  dir_cfmedip_object,
  "SOLID_cfMeDIP_matched13_1kb_log2CPM.rds"
)

file_medip_annotation <- file.path(
  dir_cfmedip_object,
  "SOLID_cfMeDIP_matched13_1kb_region_annotation.rds"
)


##------------------------------------------------------------
# Input check
##------------------------------------------------------------

input_files <- c(
  plasma_5base_beta =
    file_5base_beta,

  plasma_5base_annotation =
    file_5base_annotation,

  cfMeDIP_log2CPM =
    file_medip,

  cfMeDIP_annotation =
    file_medip_annotation
)

input_check <- data.frame(
  Input = names(input_files),
  Path = unname(input_files),
  Exists = file.exists(input_files),
  stringsAsFactors = FALSE
)


cat("\n")
cat("====================================================\n")
cat("SCRIPT 12 INPUT CHECK\n")
cat("====================================================\n")

print(
  input_check
)


if (
  !all(
    input_check$Exists
  )
) {

  stop(
    "One or more Script 12 input files are missing."
  )
}


##------------------------------------------------------------
# 3. Load matched plasma 5-base Beta matrix
##------------------------------------------------------------

beta_5base <- readRDS(
  file_5base_beta
)

beta_5base <- as.matrix(
  beta_5base
)

storage.mode(
  beta_5base
) <- "numeric"


cat("\n")
cat("====================================================\n")
cat("MATCHED PLASMA 5-BASE MATRIX\n")
cat("====================================================\n")

cat(
  "Regions:",
  format(
    nrow(beta_5base),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Patients:",
  ncol(beta_5base),
  "\n"
)

cat(
  "\nPatient IDs:\n"
)

print(
  colnames(beta_5base)
)


if (
  nrow(beta_5base) == 0
) {

  stop(
    "5-base matrix contains zero regions."
  )
}


if (
  ncol(beta_5base) != 13
) {

  stop(
    paste0(
      "Expected 13 5-base patients, but found ",
      ncol(beta_5base),
      "."
    )
  )
}


##------------------------------------------------------------
# 4. Load retained 5-base region annotation
##------------------------------------------------------------

annotation_5base <- read.delim(
  file_5base_annotation,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required_5base_cols <- c(
  "region_id",
  "chr",
  "start",
  "end"
)

missing_5base_cols <- setdiff(
  required_5base_cols,
  colnames(annotation_5base)
)


if (
  length(
    missing_5base_cols
  ) > 0
) {

  stop(
    paste(
      "Missing 5-base annotation columns:",
      paste(
        missing_5base_cols,
        collapse = ", "
      )
    )
  )
}


if (
  nrow(annotation_5base) !=
    nrow(beta_5base)
) {

  stop(
    "5-base matrix and annotation have different row counts."
  )
}


##------------------------------------------------------------
# Canonical region-ID helper
##------------------------------------------------------------

canonical_region_id <- function(
    chr,
    start,
    end
) {

  chr <- as.character(chr)

  chr <- sub(
    "^chr",
    "",
    chr,
    ignore.case = TRUE
  )

  chr <- paste0(
    "chr",
    chr
  )

  paste(
    chr,
    as.integer(start),
    as.integer(end),
    sep = ":"
  )
}


##------------------------------------------------------------
# Rebuild canonical 5-base IDs
##------------------------------------------------------------

annotation_5base$canonical_region_id <- canonical_region_id(
  annotation_5base$chr,
  annotation_5base$start,
  annotation_5base$end
)


if (
  anyDuplicated(
    annotation_5base$canonical_region_id
  )
) {

  stop(
    "Duplicated 5-base region IDs detected."
  )
}


rownames(
  beta_5base
) <- annotation_5base$canonical_region_id


##------------------------------------------------------------
# 5. Load corrected cfMeDIP 1-kb matrix
##------------------------------------------------------------

medip <- readRDS(
  file_medip
)

medip <- as.matrix(
  medip
)

storage.mode(
  medip
) <- "numeric"


medip_annotation <- readRDS(
  file_medip_annotation
)

medip_annotation <- as.data.frame(
  medip_annotation
)


cat("\n")
cat("====================================================\n")
cat("MATCHED cfMeDIP 1-kb MATRIX\n")
cat("====================================================\n")

cat(
  "Regions:",
  format(
    nrow(medip),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Patients:",
  ncol(medip),
  "\n"
)

cat(
  "Annotation rows:",
  format(
    nrow(medip_annotation),
    big.mark = ","
  ),
  "\n"
)

cat(
  "\nPatient IDs:\n"
)

print(
  colnames(medip)
)


##------------------------------------------------------------
# Critical cfMeDIP checks
##------------------------------------------------------------

if (
  nrow(medip) == 0
) {

  stop(
    "cfMeDIP matrix contains zero regions."
  )
}


if (
  ncol(medip) != 13
) {

  stop(
    paste0(
      "Expected 13 cfMeDIP patients, but found ",
      ncol(medip),
      "."
    )
  )
}


if (
  nrow(medip_annotation) !=
    nrow(medip)
) {

  stop(
    "cfMeDIP matrix and annotation have different row counts."
  )
}


if (
  is.null(
    rownames(medip)
  )
) {

  stop(
    "cfMeDIP matrix has no region IDs."
  )
}


##------------------------------------------------------------
# Confirm canonical cfMeDIP IDs
##------------------------------------------------------------

valid_medip_id <- grepl(
  "^chr[^:]+:[0-9]+:[0-9]+$",
  rownames(medip)
)


cat(
  "\nValid cfMeDIP canonical IDs:",
  sum(valid_medip_id),
  "/",
  length(valid_medip_id),
  "\n"
)


if (
  !all(
    valid_medip_id
  )
) {

  cat(
    "\nUnexpected cfMeDIP region IDs:\n"
  )

  print(
    head(
      rownames(medip)[
        !valid_medip_id
      ],
      20
    )
  )

  stop(
    "Some cfMeDIP region IDs are not in chr:start:end format."
  )
}


if (
  anyDuplicated(
    rownames(medip)
  )
) {

  stop(
    "Duplicated cfMeDIP region IDs detected."
  )
}


##------------------------------------------------------------
# Annotation correspondence
##------------------------------------------------------------

if (
  !"region_id" %in%
    colnames(medip_annotation)
) {

  stop(
    "cfMeDIP annotation does not contain region_id."
  )
}


if (
  !identical(
    rownames(medip),
    as.character(
      medip_annotation$region_id
    )
  )
) {

  stop(
    "cfMeDIP matrix row names do not match annotation region IDs."
  )
}


##------------------------------------------------------------
# 6. Match patients
##------------------------------------------------------------

patients_5base <- colnames(
  beta_5base
)

patients_medip <- colnames(
  medip
)


common_patients <- intersect(
  patients_5base,
  patients_medip
)


cat("\n")
cat("====================================================\n")
cat("PATIENT MATCHING\n")
cat("====================================================\n")

cat(
  "5-base patients:",
  length(patients_5base),
  "\n"
)

cat(
  "cfMeDIP patients:",
  length(patients_medip),
  "\n"
)

cat(
  "Common patients:",
  length(common_patients),
  "\n"
)

print(
  common_patients
)


if (
  length(
    common_patients
  ) != 13
) {

  stop(
    paste0(
      "Expected 13 common patients; found ",
      length(common_patients),
      "."
    )
  )
}


##------------------------------------------------------------
# Preserve 5-base patient order
##------------------------------------------------------------

patient_order <- patients_5base[
  patients_5base %in%
    common_patients
]


beta_5base <- beta_5base[
  ,
  patient_order,
  drop = FALSE
]

medip <- medip[
  ,
  patient_order,
  drop = FALSE
]


stopifnot(
  identical(
    colnames(beta_5base),
    colnames(medip)
  )
)


##------------------------------------------------------------
# 7. Match genomic regions
##------------------------------------------------------------

common_regions <- intersect(
  rownames(beta_5base),
  rownames(medip)
)


cat("\n")
cat("====================================================\n")
cat("REGION MATCHING\n")
cat("====================================================\n")

cat(
  "Filtered 5-base regions:",
  format(
    nrow(beta_5base),
    big.mark = ","
  ),
  "\n"
)

cat(
  "cfMeDIP-supported regions:",
  format(
    nrow(medip),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Common regions:",
  format(
    length(common_regions),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Percent cfMeDIP regions matched:",
  round(
    100 *
      length(common_regions) /
      nrow(medip),
    2
  ),
  "%\n"
)


if (
  length(
    common_regions
  ) == 0
) {

  stop(
    "No common regions found between 5-base and cfMeDIP."
  )
}


if (
  length(
    common_regions
  ) !=
    nrow(
      medip
    )
) {

  warning(
    paste0(
      "Not all cfMeDIP regions matched 5-base. ",
      "cfMeDIP = ",
      nrow(medip),
      "; common = ",
      length(common_regions)
    )
  )
}


##------------------------------------------------------------
# Preserve 5-base region order
##------------------------------------------------------------

region_order <- rownames(
  beta_5base
)

region_order <- region_order[
  region_order %in%
    common_regions
]


beta_common <- beta_5base[
  region_order,
  ,
  drop = FALSE
]

medip_common <- medip[
  region_order,
  ,
  drop = FALSE
]


stopifnot(
  identical(
    rownames(beta_common),
    rownames(medip_common)
  ),

  identical(
    colnames(beta_common),
    colnames(medip_common)
  )
)


cat("\n")
cat("Final concordance matrices:\n")

cat(
  "5-base:",
  nrow(beta_common),
  "x",
  ncol(beta_common),
  "\n"
)

cat(
  "cfMeDIP:",
  nrow(medip_common),
  "x",
  ncol(medip_common),
  "\n"
)


##------------------------------------------------------------
# 8. Save common matrices
##------------------------------------------------------------

saveRDS(
  beta_common,
  file.path(
    dir_object,
    "SOLID_5base_common_with_cfMeDIP_beta.rds"
  )
)

saveRDS(
  medip_common,
  file.path(
    dir_object,
    "SOLID_cfMeDIP_common_with_5base_log2CPM.rds"
  )
)


##------------------------------------------------------------
# 9. Patient-wise concordance
##------------------------------------------------------------

patient_results <- lapply(
  patient_order,
  function(patient) {

    x <- beta_common[
      ,
      patient
    ]

    y <- medip_common[
      ,
      patient
    ]

    keep <- (
      is.finite(x) &
        is.finite(y)
    )

    n_valid <- sum(
      keep
    )


    if (
      n_valid < 3
    ) {

      return(
        data.frame(
          Patient_ID = patient,
          N_regions = n_valid,
          Spearman_rho = NA_real_,
          P_value = NA_real_,
          stringsAsFactors = FALSE
        )
      )
    }


    ct <- suppressWarnings(
      cor.test(
        x[
          keep
        ],
        y[
          keep
        ],
        method = "spearman",
        exact = FALSE
      )
    )


    data.frame(
      Patient_ID =
        patient,

      N_regions =
        n_valid,

      Spearman_rho =
        unname(
          ct$estimate
        ),

      P_value =
        ct$p.value,

      stringsAsFactors = FALSE
    )
  }
)


patient_concordance <- do.call(
  rbind,
  patient_results
)


patient_concordance$FDR <- p.adjust(
  patient_concordance$P_value,
  method = "BH"
)


cat("\n")
cat("====================================================\n")
cat("PATIENT-WISE CONCORDANCE\n")
cat("====================================================\n")

print(
  patient_concordance
)


write.csv(
  patient_concordance,
  file.path(
    dir_table,
    "SOLID_5base_cfMeDIP_patient_concordance.csv"
  ),
  row.names = FALSE
)


##------------------------------------------------------------
# 10. Patient-wise summary
##------------------------------------------------------------

patient_summary <- data.frame(

  Metric = c(
    "Patients",
    "Median Spearman rho",
    "Mean Spearman rho",
    "Minimum Spearman rho",
    "Maximum Spearman rho",
    "Positive rho patients",
    "Negative rho patients"
  ),

  Value = c(

    nrow(
      patient_concordance
    ),

    median(
      patient_concordance$Spearman_rho,
      na.rm = TRUE
    ),

    mean(
      patient_concordance$Spearman_rho,
      na.rm = TRUE
    ),

    min(
      patient_concordance$Spearman_rho,
      na.rm = TRUE
    ),

    max(
      patient_concordance$Spearman_rho,
      na.rm = TRUE
    ),

    sum(
      patient_concordance$Spearman_rho > 0,
      na.rm = TRUE
    ),

    sum(
      patient_concordance$Spearman_rho < 0,
      na.rm = TRUE
    )
  ),

  stringsAsFactors = FALSE
)


print(
  patient_summary
)


write.csv(
  patient_summary,
  file.path(
    dir_table,
    "SOLID_5base_cfMeDIP_patient_concordance_summary.csv"
  ),
  row.names = FALSE
)


##------------------------------------------------------------
# 11. Patient-wise Spearman plot
##------------------------------------------------------------

patient_plot_df <- patient_concordance

patient_plot_df$Patient_ID <- factor(
  patient_plot_df$Patient_ID,
  levels =
    patient_plot_df$Patient_ID[
      order(
        patient_plot_df$Spearman_rho
      )
    ]
)


p_patient <- ggplot(
  patient_plot_df,
  aes(
    x = Patient_ID,
    y = Spearman_rho
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.4
  ) +
  geom_point(
    size = 3
  ) +
  coord_flip() +
  labs(
    title =
      "Plasma 5-base vs cfMeDIP concordance",
    subtitle =
      "Spearman correlation across common 1-kb regions",
    x = NULL,
    y = "Spearman rho"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )


ggsave(
  file.path(
    dir_fig,
    "SOLID_5base_cfMeDIP_patient_Spearman.pdf"
  ),
  p_patient,
  width = 7,
  height = 5.5
)

ggsave(
  file.path(
    dir_fig,
    "SOLID_5base_cfMeDIP_patient_Spearman.png"
  ),
  p_patient,
  width = 7,
  height = 5.5,
  dpi = 300
)


##------------------------------------------------------------
# 12. Cross-patient correlation matrix
#
# Rows:
#   5-base patients
#
# Columns:
#   cfMeDIP patients
##------------------------------------------------------------

cross_patient_cor <- cor(
  beta_common,
  medip_common,
  method = "spearman",
  use = "pairwise.complete.obs"
)


##------------------------------------------------------------
# Force exact matched patient order
##------------------------------------------------------------

cross_patient_cor <- cross_patient_cor[
  patient_order,
  patient_order,
  drop = FALSE
]


stopifnot(
  identical(
    rownames(cross_patient_cor),
    patient_order
  ),

  identical(
    colnames(cross_patient_cor),
    patient_order
  ),

  identical(
    rownames(cross_patient_cor),
    colnames(cross_patient_cor)
  )
)


write.csv(
  cross_patient_cor,
  file.path(
    dir_table,
    "SOLID_5base_cfMeDIP_cross_patient_correlation_matrix.csv"
  )
)


##------------------------------------------------------------
# 13. Matched vs unmatched correlations
##------------------------------------------------------------

cross_long <- expand.grid(

  Patient_5base =
    rownames(
      cross_patient_cor
    ),

  Patient_cfMeDIP =
    colnames(
      cross_patient_cor
    ),

  stringsAsFactors = FALSE
)


cross_long$Spearman_rho <- as.vector(
  cross_patient_cor
)


cross_long$Pair_type <- ifelse(
  cross_long$Patient_5base ==
    cross_long$Patient_cfMeDIP,
  "Matched",
  "Unmatched"
)


write.csv(
  cross_long,
  file.path(
    dir_table,
    "SOLID_5base_cfMeDIP_matched_vs_unmatched_correlations.csv"
  ),
  row.names = FALSE
)


matched_unmatched_summary <- aggregate(
  Spearman_rho ~ Pair_type,
  data = cross_long,
  FUN = function(x) {

    c(
      N =
        length(x),

      Median =
        median(
          x,
          na.rm = TRUE
        ),

      Mean =
        mean(
          x,
          na.rm = TRUE
        ),

      SD =
        sd(
          x,
          na.rm = TRUE
        ),

      Min =
        min(
          x,
          na.rm = TRUE
        ),

      Max =
        max(
          x,
          na.rm = TRUE
        )
    )
  }
)


cat("\n")
cat("====================================================\n")
cat("MATCHED VS UNMATCHED CORRELATIONS\n")
cat("====================================================\n")

print(
  matched_unmatched_summary
)


write.csv(
  cross_long,
  file.path(
    dir_table,
    "SOLID_5base_cfMeDIP_matched_vs_unmatched_correlations.csv"
  ),
  row.names = FALSE
)


##------------------------------------------------------------
# Exploratory Wilcoxon comparison
##------------------------------------------------------------

matched_cor <- cross_long$Spearman_rho[
  cross_long$Pair_type == "Matched"
]

unmatched_cor <- cross_long$Spearman_rho[
  cross_long$Pair_type == "Unmatched"
]


wilcox_matched <- wilcox.test(
  matched_cor,
  unmatched_cor,
  alternative = "greater"
)


cat(
  "\nExploratory Wilcoxon test:\n"
)

print(
  wilcox_matched
)


##------------------------------------------------------------
# 14. True-match ranking
#
# For each 5-base patient:
# rank all 13 cfMeDIP samples from highest to lowest
# correlation.
#
# Rank 1 = best match.
##------------------------------------------------------------

match_rank <- sapply(
  patient_order,
  function(patient) {

    x <- cross_patient_cor[
      patient,
    ]

    ranks <- rank(
      -x,
      ties.method = "average"
    )

    as.numeric(
      ranks[
        patient
      ]
    )
  }
)


match_rank_table <- data.frame(

  Patient =
    patient_order,

  True_match_rank =
    as.numeric(
      match_rank[
        patient_order
      ]
    ),

  Matched_Spearman_rho =
    diag(
      cross_patient_cor
    ),

  stringsAsFactors = FALSE
)


cat("\n")
cat("====================================================\n")
cat("TRUE-MATCH RANKING\n")
cat("====================================================\n")

print(
  match_rank_table
)


cat(
  "\nMedian true-match rank:",
  median(
    match_rank_table$True_match_rank
  ),
  "\n"
)

cat(
  "Rank #1 matches:",
  sum(
    match_rank_table$True_match_rank == 1
  ),
  "/13\n"
)

cat(
  "Top-3 matches:",
  sum(
    match_rank_table$True_match_rank <= 3
  ),
  "/13\n"
)


write.csv(
  match_rank_table,
  file.path(
    dir_table,
    "SOLID_5base_cfMeDIP_true_match_ranking.csv"
  ),
  row.names = FALSE
)


##------------------------------------------------------------
# 15. Cross-patient heatmap
#
# IMPORTANT:
# Both axes use identical patient order.
# Matched samples therefore appear on the diagonal.
#
# No independent clustering.
##------------------------------------------------------------

cross_patient_plot <- cross_patient_cor[
  patient_order,
  patient_order,
  drop = FALSE
]


stopifnot(
  identical(
    rownames(
      cross_patient_plot
    ),
    colnames(
      cross_patient_plot
    )
  )
)


##------------------------------------------------------------
# Symmetric color scale centered at zero
##------------------------------------------------------------

cor_min <- min(
  cross_patient_plot,
  na.rm = TRUE
)

cor_max <- max(
  cross_patient_plot,
  na.rm = TRUE
)

cat(
  "\nHeatmap correlation range:",
  round(
    cor_min,
    4
  ),
  "to",
  round(
    cor_max,
    4
  ),
  "\n"
)


if (
  cor_min < 0 &&
    cor_max > 0
) {

  ## Correlations span zero

  cor_abs <- max(
    abs(
      c(
        cor_min,
        cor_max
      )
    )
  )

  color_breaks <- c(
    -cor_abs,
    0,
    cor_abs
  )

} else {

  ## Correlations are entirely positive
  ## or entirely negative

  cor_mid <- (
    cor_min +
      cor_max
  ) / 2

  color_breaks <- c(
    cor_min,
    cor_mid,
    cor_max
  )
}


cat(
  "Heatmap color breaks:",
  paste(
    round(
      color_breaks,
      4
    ),
    collapse = ", "
  ),
  "\n"
)


col_fun <- colorRamp2(
  color_breaks,
  c(
    "#2166AC",
    "#F7F7F7",
    "#B2182B"
  )
)

ht <- Heatmap(

  cross_patient_plot,

  name = "Spearman\nrho",

  col = col_fun,

  cluster_rows = FALSE,
  cluster_columns = FALSE,

  show_row_names = TRUE,
  show_column_names = TRUE,

  row_names_gp = gpar(
    fontsize = 8
  ),

  column_names_gp = gpar(
    fontsize = 8
  ),

  row_title =
    "5-base plasma",

  column_title =
    "cfMeDIP plasma",

  column_title_gp = gpar(
    fontface = "bold"
  ),

  heatmap_legend_param = list(
    title = "Spearman\nrho"
  )
)


pdf(
  file.path(
    dir_fig,
    "SOLID_5base_cfMeDIP_cross_patient_heatmap_matched_order.pdf"
  ),
  width = 8,
  height = 8
)

draw(
  ht,
  heatmap_legend_side = "right"
)

dev.off()


png(
  file.path(
    dir_fig,
    "SOLID_5base_cfMeDIP_cross_patient_heatmap_matched_order.png"
  ),
  width = 2400,
  height = 2400,
  res = 300
)

draw(
  ht,
  heatmap_legend_side = "right"
)

dev.off()


##------------------------------------------------------------
# 16. Region-wise concordance
#
# For each common 1-kb region:
# correlate 5-base Beta vs cfMeDIP across patients.
##------------------------------------------------------------

minimum_valid_patients <- 10L


region_results <- lapply(
  seq_len(
    nrow(
      beta_common
    )
  ),
  function(i) {

    x <- beta_common[
      i,
    ]

    y <- medip_common[
      i,
    ]

    keep <- (
      is.finite(x) &
        is.finite(y)
    )

    n_valid <- sum(
      keep
    )


    if (
      n_valid <
        minimum_valid_patients
    ) {

      return(
        c(
          N_valid =
            n_valid,

          Spearman_rho =
            NA_real_
        )
      )
    }


    rho <- suppressWarnings(
      cor(
        x[
          keep
        ],
        y[
          keep
        ],
        method = "spearman"
      )
    )


    c(
      N_valid =
        n_valid,

      Spearman_rho =
        rho
    )
  }
)


region_results <- do.call(
  rbind,
  region_results
)


region_concordance <- data.frame(

  region_id =
    rownames(
      beta_common
    ),

  N_valid =
    as.integer(
      region_results[
        ,
        "N_valid"
      ]
    ),

  Spearman_rho =
    as.numeric(
      region_results[
        ,
        "Spearman_rho"
      ]
    ),

  stringsAsFactors = FALSE
)


##------------------------------------------------------------
# 17. Add genomic coordinates
##------------------------------------------------------------

region_split <- strsplit(
  region_concordance$region_id,
  ":",
  fixed = TRUE
)


region_concordance$chr <- vapply(
  region_split,
  function(x) x[1],
  character(1)
)


region_concordance$start <- as.integer(
  vapply(
    region_split,
    function(x) x[2],
    character(1)
  )
)


region_concordance$end <- as.integer(
  vapply(
    region_split,
    function(x) x[3],
    character(1)
  )
)


write.table(
  region_concordance,
  file.path(
    dir_table,
    "SOLID_5base_cfMeDIP_region_concordance.tsv.gz"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


##------------------------------------------------------------
# 18. Region-wise summary
##------------------------------------------------------------

valid_region_rho <- region_concordance$Spearman_rho[
  is.finite(
    region_concordance$Spearman_rho
  )
]


region_summary <- data.frame(

  Metric = c(
    "Common regions",
    "Regions with >=10 valid patients",
    "Median region Spearman rho",
    "Mean region Spearman rho",
    "Positive-rho regions",
    "Negative-rho regions"
  ),

  Value = c(

    nrow(
      region_concordance
    ),

    length(
      valid_region_rho
    ),

    median(
      valid_region_rho,
      na.rm = TRUE
    ),

    mean(
      valid_region_rho,
      na.rm = TRUE
    ),

    sum(
      valid_region_rho > 0,
      na.rm = TRUE
    ),

    sum(
      valid_region_rho < 0,
      na.rm = TRUE
    )
  ),

  stringsAsFactors = FALSE
)


cat("\n")
cat("====================================================\n")
cat("REGION-WISE CONCORDANCE SUMMARY\n")
cat("====================================================\n")

print(
  region_summary
)


write.csv(
  region_summary,
  file.path(
    dir_table,
    "SOLID_5base_cfMeDIP_region_concordance_summary.csv"
  ),
  row.names = FALSE
)


##------------------------------------------------------------
# 19. Region-wise rho distribution
##------------------------------------------------------------

p_region <- ggplot(
  region_concordance[
    is.finite(
      region_concordance$Spearman_rho
    ),
    ,
    drop = FALSE
  ],
  aes(
    x = Spearman_rho
  )
) +
  geom_histogram(
    bins = 60
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.4
  ) +
  labs(
    title =
      "Region-wise plasma assay concordance",
    subtitle =
      "5-base Beta vs cfMeDIP log2CPM across matched patients",
    x =
      "Region-wise Spearman rho",
    y =
      "Number of 1-kb regions"
  ) +
  theme_bw() +
  theme(
    plot.title =
      element_text(
        face = "bold"
      )
  )


ggsave(
  file.path(
    dir_fig,
    "SOLID_5base_cfMeDIP_region_Spearman_distribution.pdf"
  ),
  p_region,
  width = 7,
  height = 5
)

ggsave(
  file.path(
    dir_fig,
    "SOLID_5base_cfMeDIP_region_Spearman_distribution.png"
  ),
  p_region,
  width = 7,
  height = 5,
  dpi = 300
)


##------------------------------------------------------------
# 20. Save checkpoint
##------------------------------------------------------------

saveRDS(
  list(

    beta_5base =
      beta_common,

    cfMeDIP_log2CPM =
      medip_common,

    common_patients =
      patient_order,

    common_regions =
      region_order,

    patient_concordance =
      patient_concordance,

    cross_patient_correlation =
      cross_patient_cor,

    matched_unmatched =
      cross_long,

    matched_unmatched_summary =
      matched_unmatched_summary,

    wilcox_matched =
      wilcox_matched,

    match_rank =
      match_rank_table,

    region_concordance =
      region_concordance,

    patient_summary =
      patient_summary,

    region_summary =
      region_summary

  ),

  file.path(
    dir_object,
    "SOLID_5base_cfMeDIP_concordance_checkpoint.rds"
  )
)


##------------------------------------------------------------
# 21. Final validation
##------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("FINAL SCRIPT 12 VALIDATION\n")
cat("====================================================\n")

cat(
  "Common patients:",
  length(
    patient_order
  ),
  "\n"
)

cat(
  "Common regions:",
  format(
    length(
      region_order
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "5-base final matrix:",
  nrow(beta_common),
  "x",
  ncol(beta_common),
  "\n"
)

cat(
  "cfMeDIP final matrix:",
  nrow(medip_common),
  "x",
  ncol(medip_common),
  "\n"
)

cat(
  "Row IDs identical:",
  identical(
    rownames(beta_common),
    rownames(medip_common)
  ),
  "\n"
)

cat(
  "Patient IDs identical:",
  identical(
    colnames(beta_common),
    colnames(medip_common)
  ),
  "\n"
)

cat(
  "Heatmap rows/columns identical:",
  identical(
    rownames(cross_patient_plot),
    colnames(cross_patient_plot)
  ),
  "\n"
)


stopifnot(

  length(
    patient_order
  ) == 13,

  nrow(
    beta_common
  ) > 0,

  identical(
    dim(
      beta_common
    ),
    dim(
      medip_common
    )
  ),

  identical(
    rownames(
      beta_common
    ),
    rownames(
      medip_common
    )
  ),

  identical(
    colnames(
      beta_common
    ),
    colnames(
      medip_common
    )
  ),

  identical(
    rownames(
      cross_patient_plot
    ),
    colnames(
      cross_patient_plot
    )
  )
)


##------------------------------------------------------------
# 22. Final summary
##------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("SOLID 5-BASE <-> cfMeDIP CONCORDANCE COMPLETE\n")
cat("====================================================\n")

cat(
  "Matched patients:",
  length(
    patient_order
  ),
  "\n"
)

cat(
  "Common 1-kb regions:",
  format(
    length(
      region_order
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Median patient-wise Spearman rho:",
  round(
    median(
      patient_concordance$Spearman_rho,
      na.rm = TRUE
    ),
    4
  ),
  "\n"
)

cat(
  "Matched-pair median rho:",
  round(
    median(
      matched_cor,
      na.rm = TRUE
    ),
    4
  ),
  "\n"
)

cat(
  "Unmatched-pair median rho:",
  round(
    median(
      unmatched_cor,
      na.rm = TRUE
    ),
    4
  ),
  "\n"
)

cat(
  "Median true-match rank:",
  median(
    match_rank_table$True_match_rank
  ),
  "\n"
)

cat(
  "Top-3 true matches:",
  sum(
    match_rank_table$True_match_rank <= 3
  ),
  "/13\n"
)

cat(
  "Median region-wise Spearman rho:",
  round(
    median(
      valid_region_rho,
      na.rm = TRUE
    ),
    4
  ),
  "\n"
)

cat(
  "\nOutput directory:\n",
  dir_result,
  "\n"
)

sessionInfo()