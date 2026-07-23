##------------------------------------------------------------
# DNA METHYLATION PROBE FILTERING AND PLATFORM HARMONIZATION
# Brain IDH-mutant FFPE tissue
#
# Platforms:
#   - Illumina MethylationEPIC v1: 3 retained samples
#   - Illumina MethylationEPIC v2: 25 retained samples
#
# Input data:
#   - Platform-specific Noob-normalized MethylSet objects
#   - Noob-normalized Beta-value and M-value matrices
#   - Raw detection P-value matrices
#   - Retained sample metadata and QC annotations
#
# Current objective:
#   1. Load the post-Noob normalization checkpoint
#   2. Restrict detection P-values to retained samples
#   3. Filter probes using raw detection P-values
#   4. Remove probes with missing or non-finite measurements
#   5. Standardize EPIC v2 probe identifiers
#   6. Identify EPIC v1/EPIC v2 manifest-level overlap
#   7. Restrict EPIC v2 probes to CpGs represented on EPIC v1
#   8. Collapse replicated EPIC v2 probes to canonical CpG IDs
#   9. Select CpGs retained on both platforms
#  10. Combine harmonized Beta-value and M-value matrices
#  11. Create and validate combined sample metadata
#  12. Assess residual platform effects using PCA
#  13. Save harmonized matrices, metadata, and QC summaries
#
# Important:
#   - Probe filtering is performed separately for each platform.
#   - EPIC v1 and EPIC v2 are combined only after normalization,
#     filtering, and probe-ID harmonization.
#   - EPIC v2 suffixes are removed to recover canonical CpG IDs.
#   - Replicated EPIC v2 features are summarized before merging.
#   - Beta and M-values are harmonized using canonical CpG IDs.
#   - hg19 and hg38 genomic coordinates are not mixed directly.
#   - Unified genomic annotation will be added after the final
#     common CpG set has been created.
##------------------------------------------------------------
############################################################
## Load libraries
############################################################
suppressPackageStartupMessages({
  library(minfi)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(GenomicRanges)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(IlluminaHumanMethylationEPICmanifest)
  library(IlluminaHumanMethylationEPICv2manifest)
})

############################################################
## functions
############################################################
collapse_EPICv2_replicates <- function(
    matrix_data,
    canonical_ids,
    summary_function = c("median", "mean")
) {
  summary_function <- match.arg(summary_function)

  if (!is.matrix(matrix_data)) {
    matrix_data <- as.matrix(matrix_data)
  }

  if (length(canonical_ids) != nrow(matrix_data)) {
    stop(
      "The number of canonical IDs (", length(canonical_ids),
      ") must equal the number of matrix rows (", nrow(matrix_data), ")."
    )
  }

  if (anyNA(canonical_ids)) {
    stop("canonical_ids contains missing values.")
  }

  # Split matrix row indices by canonical CpG ID
  probe_groups <- split(seq_len(nrow(matrix_data)), canonical_ids)

  collapsed_matrix <- vapply(
    probe_groups,
    FUN = function(idx) {
      if (length(idx) == 1L) {
        return(as.numeric(matrix_data[idx, ]))
      }

      if (summary_function == "median") {
        apply(
          matrix_data[idx, , drop = FALSE],
          2,
          median,
          na.rm = TRUE
        )
      } else {
        colMeans(
          matrix_data[idx, , drop = FALSE],
          na.rm = TRUE
        )
      }
    },
    FUN.VALUE = numeric(ncol(matrix_data))
  )

  # vapply returns samples × CpGs, so transpose
  collapsed_matrix <- t(collapsed_matrix)

  rownames(collapsed_matrix) <- names(probe_groups)
  colnames(collapsed_matrix) <- colnames(matrix_data)

  collapsed_matrix
}

############################################################
## Define directories
############################################################
dir_input   <- "result/data"
dir_output  <- "result/data"
dir_figures <- "result/qc"

############################################################
## Load post-Noob checkpoint
############################################################
load(file.path(dir_input, "post_noob_checkpoint.RData"))

############################################################
## Restrict detection P-values to retained samples
############################################################
detP_EPICv1_retained <- detP_EPICv1[
  ,
  colnames(beta_noob_EPICv1),
  drop = FALSE
]

detP_EPICv2_retained <- detP_EPICv2[
  ,
  colnames(beta_noob_EPICv2),
  drop = FALSE
]

stopifnot(
  identical(
    colnames(detP_EPICv1_retained),
    colnames(beta_noob_EPICv1)
  )
)

stopifnot(
  identical(
    colnames(detP_EPICv2_retained),
    colnames(beta_noob_EPICv2)
  )
)

stopifnot(
  identical(
    rownames(detP_EPICv1_retained),
    rownames(beta_noob_EPICv1)
  )
)

stopifnot(
  identical(
    rownames(detP_EPICv2_retained),
    rownames(beta_noob_EPICv2)
  )
)

############################################################
## Filter unreliable probes using detection P-values
############################################################
detection_cutoff <- 0.01

# EPIC v1:
# retain probes passing detection in all 3 retained samples
keep_probe_EPICv1 <- rowSums(
  detP_EPICv1_retained >
    detection_cutoff,
  na.rm = TRUE
) == 0

# EPIC v2:
# retain probes failing detection in no more than 1 of 25 samples
max_failed_EPICv2 <- 1

keep_probe_EPICv2 <- rowSums(
  detP_EPICv2_retained >
    detection_cutoff,
  na.rm = TRUE
) <= max_failed_EPICv2

print(
  table(keep_probe_EPICv1)
)

print(
  table(keep_probe_EPICv2)
)

beta_filtered_EPICv1 <- beta_noob_EPICv1[
  keep_probe_EPICv1,
  ,
  drop = FALSE
]

mval_filtered_EPICv1 <- mval_noob_EPICv1[
  keep_probe_EPICv1,
  ,
  drop = FALSE
]

beta_filtered_EPICv2 <- beta_noob_EPICv2[
  keep_probe_EPICv2,
  ,
  drop = FALSE
]

mval_filtered_EPICv2 <- mval_noob_EPICv2[
  keep_probe_EPICv2,
  ,
  drop = FALSE
]

############################################################
## Remove probes with missing or non-finite measurements
############################################################

keep_valid_EPICv1 <- (
  rowSums(
    is.na(beta_filtered_EPICv1)
  ) == 0 &
    rowSums(
      !is.finite(mval_filtered_EPICv1)
    ) == 0
)

keep_valid_EPICv2 <- (
  rowSums(
    is.na(beta_filtered_EPICv2)
  ) == 0 &
    rowSums(
      !is.finite(mval_filtered_EPICv2)
    ) == 0
)

print(
  table(keep_valid_EPICv1)
)

print(
  table(keep_valid_EPICv2)
)

beta_filtered_EPICv1 <- beta_filtered_EPICv1[
  keep_valid_EPICv1,
  ,
  drop = FALSE
]

mval_filtered_EPICv1 <- mval_filtered_EPICv1[
  keep_valid_EPICv1,
  ,
  drop = FALSE
]

beta_filtered_EPICv2 <- beta_filtered_EPICv2[
  keep_valid_EPICv2,
  ,
  drop = FALSE
]

mval_filtered_EPICv2 <- mval_filtered_EPICv2[
  keep_valid_EPICv2,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(beta_filtered_EPICv1),
    rownames(mval_filtered_EPICv1)
  )
)

stopifnot(
  identical(
    rownames(beta_filtered_EPICv2),
    rownames(mval_filtered_EPICv2)
  )
)

############################################################
## Summarize initial probe filtering
############################################################

probe_filter_summary <- tibble(
  Platform = c(
    "EPICv1",
    "EPICv2"
  ),

  Probes_before = c(
    nrow(beta_noob_EPICv1),
    nrow(beta_noob_EPICv2)
  ),

  Removed_detection = c(
    sum(!keep_probe_EPICv1),
    sum(!keep_probe_EPICv2)
  ),

  Retained_after_detection = c(
    sum(keep_probe_EPICv1),
    sum(keep_probe_EPICv2)
  ),

  Removed_invalid = c(
    sum(!keep_valid_EPICv1),
    sum(!keep_valid_EPICv2)
  ),

  Final_after_initial_QC = c(
    nrow(beta_filtered_EPICv1),
    nrow(beta_filtered_EPICv2)
  )
)

print(probe_filter_summary)

write.csv(
  probe_filter_summary,
  file.path(
    dir_output,
    "probe_filter_summary.csv"
  ),
  row.names = FALSE
)

############################################################
## Inspect platform-specific probe identifiers
############################################################

print(
  head(
    rownames(beta_filtered_EPICv1)
  )
)

print(
  head(
    rownames(beta_filtered_EPICv2)
  )
)

cat(
  "Duplicated EPIC v1 row names:",
  anyDuplicated(
    rownames(beta_filtered_EPICv1)
  ),
  "\n"
)

cat(
  "Duplicated EPIC v2 row names:",
  anyDuplicated(
    rownames(beta_filtered_EPICv2)
  ),
  "\n"
)

cat(
  "EPIC v1 IDs containing suffixes:",
  sum(
    grepl(
      "_",
      rownames(beta_filtered_EPICv1)
    )
  ),
  "\n"
)

cat(
  "EPIC v2 IDs containing suffixes:",
  sum(
    grepl(
      "_",
      rownames(beta_filtered_EPICv2)
    )
  ),
  "\n"
)

############################################################
## Standardize EPIC v2 probe identifiers
############################################################
probe_map_EPICv2 <- tibble(
  EPICv2_ID = rownames(
    beta_filtered_EPICv2
  ),
  Canonical_CpG = sub(
    "_.*$",
    "",
    rownames(beta_filtered_EPICv2)
  )
)

canonical_id_EPICv2 <- probe_map_EPICv2$Canonical_CpG

print(
  head(probe_map_EPICv2)
)

############################################################
## Summarize replicated EPIC v2 CpG IDs
############################################################
replicate_counts_EPICv2 <- table(
  canonical_id_EPICv2
)

replication_summary_EPICv2 <- as.data.frame(
  table(
    replicate_counts_EPICv2
  )
)

colnames(
  replication_summary_EPICv2
) <- c(
  "Number_of_features_per_CpG",
  "Number_of_CpGs"
)

n_unique_EPICv2 <- length(
  unique(
    canonical_id_EPICv2
  )
)

n_duplicated_rows_EPICv2 <- sum(
  duplicated(
    canonical_id_EPICv2
  )
)

n_replicated_CpGs_EPICv2 <- sum(
  replicate_counts_EPICv2 > 1
)

cat(
  "Filtered EPIC v2 array features:",
  nrow(beta_filtered_EPICv2),
  "\n"
)

cat(
  "Unique canonical EPIC v2 CpGs:",
  n_unique_EPICv2,
  "\n"
)

cat(
  "Duplicated EPIC v2 array-feature rows after suffix removal:",
  n_duplicated_rows_EPICv2,
  "\n"
)

cat(
  "Canonical CpGs represented by multiple EPIC v2 features:",
  n_replicated_CpGs_EPICv2,
  "\n"
)

print(
  replication_summary_EPICv2
)

duplicated_cpgs_EPICv2 <- names(
  replicate_counts_EPICv2[
    replicate_counts_EPICv2 > 1
  ]
)

duplicated_probe_examples_EPICv2 <- probe_map_EPICv2 %>%
  filter(
    Canonical_CpG %in%
      head(
        duplicated_cpgs_EPICv2,
        5
      )
  )

print(
  duplicated_probe_examples_EPICv2
)

############################################################
## Identify manifest-level overlap between EPIC v1 and v2
############################################################

manifest_probes_EPICv1 <- getManifestInfo(
  IlluminaHumanMethylationEPICmanifest,
  "locusNames"
)

manifest_probes_EPICv2 <- getManifestInfo(
  IlluminaHumanMethylationEPICv2manifest,
  "locusNames"
)

manifest_probes_EPICv1 <- unique(
  manifest_probes_EPICv1
)

manifest_probes_EPICv2_canonical <- unique(
  sub(
    "_.*$",
    "",
    manifest_probes_EPICv2
  )
)

common_manifest_probes <- intersect(
  manifest_probes_EPICv1,
  manifest_probes_EPICv2_canonical
)

cat(
  "Unique EPIC v1 manifest CpGs:",
  length(manifest_probes_EPICv1),
  "\n"
)

cat(
  "Unique EPIC v2 canonical manifest CpGs:",
  length(manifest_probes_EPICv2_canonical),
  "\n"
)

cat(
  "Common manifest CpGs:",
  length(common_manifest_probes),
  "\n"
)

############################################################
## Restrict EPIC v2 features to CpGs retained on EPIC v1
############################################################

keep_EPICv2_overlap <- (
  canonical_id_EPICv2 %in%
    rownames(beta_filtered_EPICv1)
)

print(
  table(keep_EPICv2_overlap)
)

beta_EPICv2_overlap <- beta_filtered_EPICv2[
  keep_EPICv2_overlap,
  ,
  drop = FALSE
]

mval_EPICv2_overlap <- mval_filtered_EPICv2[
  keep_EPICv2_overlap,
  ,
  drop = FALSE
]

canonical_overlap_EPICv2 <- canonical_id_EPICv2[
  keep_EPICv2_overlap
]

stopifnot(
  nrow(beta_EPICv2_overlap) ==
    length(canonical_overlap_EPICv2)
)


############################################################
## Collapse replicated EPIC v2 features by canonical CpG ID
############################################################

# Use mean aggregation to remain consistent with the reference
# workflow. Change to "median" for a sensitivity analysis if needed.
replicate_summary_method <- "mean"

beta_EPICv2_canonical <- collapse_EPICv2_replicates(
  matrix_data = beta_EPICv2_overlap,
  canonical_ids = canonical_overlap_EPICv2,
  summary_function = replicate_summary_method
)

mval_EPICv2_canonical <- collapse_EPICv2_replicates(
  matrix_data = mval_EPICv2_overlap,
  canonical_ids = canonical_overlap_EPICv2,
  summary_function = replicate_summary_method
)

stopifnot(
  !anyDuplicated(
    rownames(beta_EPICv2_canonical)
  )
)

stopifnot(
  !anyDuplicated(
    rownames(mval_EPICv2_canonical)
  )
)

stopifnot(
  identical(
    rownames(beta_EPICv2_canonical),
    rownames(mval_EPICv2_canonical)
  )
)

stopifnot(
  identical(
    colnames(beta_EPICv2_canonical),
    colnames(beta_filtered_EPICv2)
  )
)

cat(
  "Canonical EPIC v2 CpGs after replicate collapse:",
  nrow(beta_EPICv2_canonical),
  "\n"
)


############################################################
## Identify CpGs retained on both platforms
############################################################

common_filtered_probes <- rownames(
  beta_filtered_EPICv1
)[
  rownames(beta_filtered_EPICv1) %in%
    rownames(beta_EPICv2_canonical)
]

cat(
  "Common filtered CpGs:",
  length(common_filtered_probes),
  "\n"
)


############################################################
## Subset both platforms to the same CpGs and order
############################################################

beta_common_EPICv1 <- beta_filtered_EPICv1[
  common_filtered_probes,
  ,
  drop = FALSE
]

beta_common_EPICv2 <- beta_EPICv2_canonical[
  common_filtered_probes,
  ,
  drop = FALSE
]

mval_common_EPICv1 <- mval_filtered_EPICv1[
  common_filtered_probes,
  ,
  drop = FALSE
]

mval_common_EPICv2 <- mval_EPICv2_canonical[
  common_filtered_probes,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(beta_common_EPICv1),
    rownames(beta_common_EPICv2)
  )
)

stopifnot(
  identical(
    rownames(mval_common_EPICv1),
    rownames(mval_common_EPICv2)
  )
)

stopifnot(
  identical(
    rownames(beta_common_EPICv1),
    rownames(mval_common_EPICv1)
  )
)


############################################################
## Combine harmonized Beta-value and M-value matrices
############################################################

beta_combined <- cbind(
  beta_common_EPICv1,
  beta_common_EPICv2
)

mval_combined <- cbind(
  mval_common_EPICv1,
  mval_common_EPICv2
)

stopifnot(
  ncol(beta_combined) == 28
)

stopifnot(
  ncol(mval_combined) == 28
)

stopifnot(
  identical(
    colnames(beta_combined),
    colnames(mval_combined)
  )
)

stopifnot(
  identical(
    rownames(beta_combined),
    rownames(mval_combined)
  )
)

stopifnot(
  !anyDuplicated(
    colnames(beta_combined)
  )
)

stopifnot(
  all(
    is.finite(beta_combined)
  )
)

stopifnot(
  all(
    is.finite(mval_combined)
  )
)

cat(
  "Combined Beta matrix dimensions:",
  paste(
    dim(beta_combined),
    collapse = " x "
  ),
  "\n"
)

cat(
  "Combined M-value matrix dimensions:",
  paste(
    dim(mval_combined),
    collapse = " x "
  ),
  "\n"
)


############################################################
## Create and align combined sample metadata
############################################################

targets_combined <- bind_rows(
  targets_EPICv1_qc,
  targets_EPICv2_qc
)

targets_combined <- targets_combined[
  match(
    colnames(beta_combined),
    targets_combined$Sample_Name
  ),
  ,
  drop = FALSE
]

stopifnot(
  !anyNA(
    targets_combined$Sample_Name
  )
)

stopifnot(
  identical(
    targets_combined$Sample_Name,
    colnames(beta_combined)
  )
)

print(
  table(
    targets_combined$Array
  )
)


############################################################
## Assess residual platform effects using PCA
############################################################

pca_combined <- run_pca(
  m_values = mval_combined,
  targets_data = targets_combined,
  n_probes = 10000
)

p_pca_platform <- ggplot(
  pca_combined$data,
  aes(
    x = PC1,
    y = PC2,
    color = Array
  )
) +
  geom_point(
    size = 3
  ) +
  geom_text(
    aes(
      label = Sample_Name
    ),
    size = 2.5,
    vjust = -0.7,
    show.legend = FALSE
  ) +
  labs(
    title = "Combined EPIC v1 and EPIC v2 PCA",
    subtitle = "Noob-normalized, detection-filtered, common CpG probes",
    x = paste0(
      "PC1 (",
      round(
        pca_combined$variance_explained[1],
        1
      ),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(
        pca_combined$variance_explained[2],
        1
      ),
      "%)"
    ),
    color = "Platform"
  ) +
  theme_bw()

ggsave(
  filename = file.path(
    dir_figures,
    "combined_EPICv1_EPICv2_PCA_by_platform.pdf"
  ),
  plot = p_pca_platform,
  width = 8,
  height = 7
)

write.csv(
  pca_combined$data,
  file.path(
    dir_output,
    "combined_EPICv1_EPICv2_PCA_coordinates.csv"
  ),
  row.names = FALSE
)


############################################################
## Summarize platform harmonization
############################################################

harmonization_summary <- tibble(
  Measure = c(
    "EPICv1 probes before filtering",
    "EPICv2 features before filtering",
    "EPICv1 probes after QC",
    "EPICv2 features after QC",
    "Unique canonical EPICv2 CpGs after QC",
    "Common manifest CpGs",
    "Common filtered CpGs",
    "Combined samples"
  ),

  Value = c(
    nrow(beta_noob_EPICv1),
    nrow(beta_noob_EPICv2),
    nrow(beta_filtered_EPICv1),
    nrow(beta_filtered_EPICv2),
    nrow(beta_EPICv2_canonical),
    length(common_manifest_probes),
    length(common_filtered_probes),
    ncol(beta_combined)
  )
)

print(
  harmonization_summary
)


############################################################
## Save harmonized matrices as RDS files
############################################################

saveRDS(
  beta_combined,
  file.path(
    dir_output,
    "beta_combined_EPICv1_EPICv2.rds"
  )
)

saveRDS(
  mval_combined,
  file.path(
    dir_output,
    "mval_combined_EPICv1_EPICv2.rds"
  )
)

saveRDS(
  targets_combined,
  file.path(
    dir_output,
    "targets_combined_EPICv1_EPICv2.rds"
  )
)


############################################################
## Save harmonized data checkpoint
############################################################

save(
  detection_cutoff,
  max_failed_EPICv2,
  replicate_summary_method,
  probe_filter_summary,
  harmonization_summary,
  probe_map_EPICv2,
  replicate_counts_EPICv2,
  replication_summary_EPICv2,
  common_manifest_probes,
  common_filtered_probes,
  beta_filtered_EPICv1,
  mval_filtered_EPICv1,
  beta_filtered_EPICv2,
  mval_filtered_EPICv2,
  beta_EPICv2_canonical,
  mval_EPICv2_canonical,
  beta_common_EPICv1,
  beta_common_EPICv2,
  mval_common_EPICv1,
  mval_common_EPICv2,
  beta_combined,
  mval_combined,
  targets_combined,
  pca_combined,
  file = file.path(
    dir_output,
    "harmonized_EPICv1_EPICv2.RData"
  )
)


############################################################
## Save summary tables and probe lists
############################################################

write.csv(
  harmonization_summary,
  file.path(
    dir_output,
    "harmonization_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  targets_combined,
  file.path(
    dir_output,
    "harmonized_sample_metadata.csv"
  ),
  row.names = FALSE
)

write.csv(
  probe_map_EPICv2,
  file.path(
    dir_output,
    "EPICv2_probe_ID_mapping.csv"
  ),
  row.names = FALSE
)

write.csv(
  replication_summary_EPICv2,
  file.path(
    dir_output,
    "EPICv2_replicate_summary.csv"
  ),
  row.names = FALSE
)

write.table(
  common_filtered_probes,
  file.path(
    dir_output,
    "common_filtered_EPICv1_EPICv2_CpGs.txt"
  ),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


############################################################
## Final validation
############################################################

harmonized_checkpoint <- file.path(
  dir_output,
  "harmonized_EPICv1_EPICv2.RData"
)

stopifnot(
  file.exists(
    harmonized_checkpoint
  )
)

stopifnot(
  identical(
    targets_combined$Sample_Name,
    colnames(beta_combined)
  )
)

stopifnot(
  identical(
    rownames(beta_combined),
    common_filtered_probes
  )
)

message(
  "Platform harmonization completed successfully."
)

message(
  "Saved checkpoint: ",
  harmonized_checkpoint
)

message(
  "Final dimensions: ",
  nrow(beta_combined),
  " CpGs x ",
  ncol(beta_combined),
  " samples."
)