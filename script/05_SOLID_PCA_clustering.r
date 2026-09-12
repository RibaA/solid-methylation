############################################################
## 05_SOLID_PCA_clustering.r
##
## SOLID MATCHED TUMOR–PLASMA
## PCA + CLUSTERING
##
## Purpose:
##   1. Load frozen downstream analysis-ready object
##   2. Use M-values for multivariate structure
##   3. Apply existing valid-pair mask
##   4. Select complete regions across all 26 samples
##   5. Rank regions by variance
##   6. PCA on top 10,000 variable regions
##   7. Sample correlation + hierarchical clustering
##   8. Heatmap of top 500 variable regions
##   9. Add clinical annotations
##
## IMPORTANT:
##   - No new biological filtering
##   - No DMR testing
##   - Variable-region selection is ONLY for
##     PCA/clustering/visualization
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
## 2. PROJECT DIRECTORIES
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
  "PCA_clustering"
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
  "SOLID PCA + CLUSTERING\n"
)

cat(
  "============================================\n"
)


############################################################
## 4. BASIC INPUT VALIDATION
############################################################

required_components <- c(
  "tissue_M",
  "plasma_M",
  "valid_pair_mask",
  "region_annotation",
  "patient_metadata",
  "sample_metadata"
)

missing_components <- setdiff(
  required_components,
  names(obj)
)

if (length(missing_components) > 0L) {

  stop(
    "Missing required object component(s): ",
    paste(
      missing_components,
      collapse = ", "
    )
  )
}


tissue_M <- obj$tissue_M
plasma_M <- obj$plasma_M
valid_mask <- obj$valid_pair_mask

region_annotation <- as.data.frame(
  obj$region_annotation
)

patient_metadata <- as.data.frame(
  obj$patient_metadata
)


stopifnot(
  identical(
    dim(tissue_M),
    dim(plasma_M)
  ),
  identical(
    dim(tissue_M),
    dim(valid_mask)
  ),
  ncol(tissue_M) == 13L,
  nrow(tissue_M) == 124961L
)


patient_ids <- colnames(
  tissue_M
)

stopifnot(
  identical(
    patient_ids,
    colnames(plasma_M)
  ),
  identical(
    patient_ids,
    patient_metadata$patient_id
  )
)

cat(
  "\nInput validation: PASS\n"
)

############################################################
## 5. APPLY EXISTING VALID-PAIR MASK
############################################################

tissue_M_valid <- tissue_M
plasma_M_valid <- plasma_M

tissue_M_valid[
  !valid_mask
] <- NA_real_

plasma_M_valid[
  !valid_mask
] <- NA_real_


############################################################
## 6. COMBINE TUMOR + PLASMA MATRIX
##
## Rows    = regions
## Columns = 26 samples
############################################################

combined_M <- cbind(
  tissue_M_valid,
  plasma_M_valid
)

colnames(
  combined_M
) <- c(
  paste0(
    "Tumor_",
    patient_ids
  ),
  paste0(
    "Plasma_",
    patient_ids
  )
)


############################################################
## 7. BUILD MATCHING SAMPLE METADATA
############################################################

sample_metadata <- rbind(

  data.frame(
    sample_key =
      paste0(
        "Tumor_",
        patient_ids
      ),

    patient_id =
      patient_ids,

    sample_type =
      "Tumor",

    patient_metadata[
      ,
      setdiff(
        names(patient_metadata),
        "patient_id"
      ),
      drop = FALSE
    ],

    stringsAsFactors = FALSE
  ),

  data.frame(
    sample_key =
      paste0(
        "Plasma_",
        patient_ids
      ),

    patient_id =
      patient_ids,

    sample_type =
      "Plasma",

    patient_metadata[
      ,
      setdiff(
        names(patient_metadata),
        "patient_id"
      ),
      drop = FALSE
    ],

    stringsAsFactors = FALSE
  )
)


stopifnot(
  identical(
    sample_metadata$sample_key,
    colnames(combined_M)
  )
)


############################################################
## 8. COMPLETE-REGION REQUIREMENT
##
## PCA/clustering should not depend on imputation.
## Keep regions finite across all 26 samples.
############################################################

complete_region <- apply(
  combined_M,
  1,
  function(x) {
    all(
      is.finite(x)
    )
  }
)

n_complete_regions <- sum(
  complete_region
)

cat(
  "\nComplete regions across all 26 samples:",
  n_complete_regions,
  "of",
  nrow(combined_M),
  "\n"
)

cat(
  "Complete-region percentage:",
  100 *
    n_complete_regions /
    nrow(combined_M),
  "\n"
)


if (n_complete_regions < 10000L) {

  stop(
    "Fewer than 10,000 complete regions available."
  )
}


combined_complete <- combined_M[
  complete_region,
  ,
  drop = FALSE
]

annotation_complete <- region_annotation[
  complete_region,
  ,
  drop = FALSE
]


############################################################
## 9. REGION VARIANCE
############################################################

region_variance <- apply(
  combined_complete,
  1,
  var
)

stopifnot(
  all(
    is.finite(
      region_variance
    )
  )
)


variance_order <- order(
  region_variance,
  decreasing = TRUE
)


############################################################
## 10. SELECT VARIABLE REGIONS
############################################################

n_pca_regions <- min(
  10000L,
  length(
    variance_order
  )
)

n_heatmap_regions <- min(
  500L,
  length(
    variance_order
  )
)


pca_region_index <- variance_order[
  seq_len(
    n_pca_regions
  )
]

heatmap_region_index <- variance_order[
  seq_len(
    n_heatmap_regions
  )
]


M_pca <- combined_complete[
  pca_region_index,
  ,
  drop = FALSE
]

M_heatmap <- combined_complete[
  heatmap_region_index,
  ,
  drop = FALSE
]


############################################################
## 11. SAVE VARIABLE-REGION TABLE
############################################################

variable_region_table <- annotation_complete[
  variance_order,
  ,
  drop = FALSE
]

variable_region_table$variance_M <-
  region_variance[
    variance_order
  ]

variable_region_table$variance_rank <-
  seq_len(
    nrow(
      variable_region_table
    )
  )

variable_region_table$used_for_PCA <-
  variable_region_table$variance_rank <=
  n_pca_regions

variable_region_table$used_for_heatmap <-
  variable_region_table$variance_rank <=
  n_heatmap_regions


write.table(
  variable_region_table,
  file.path(
    table_dir,
    "SOLID_M_value_region_variance_ranking.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
## 12. PCA
##
## Samples are observations.
############################################################

pca_fit <- prcomp(
  t(M_pca),
  center = TRUE,
  scale. = FALSE
)


variance_explained <-
  100 *
  pca_fit$sdev^2 /
  sum(
    pca_fit$sdev^2
  )


pca_coordinates <- as.data.frame(
  pca_fit$x
)

pca_coordinates$sample_key <-
  rownames(
    pca_coordinates
  )


pca_coordinates <- merge(
  pca_coordinates,
  sample_metadata,
  by = "sample_key",
  all.x = TRUE,
  sort = FALSE
)


pca_coordinates <- pca_coordinates[
  match(
    rownames(
      pca_fit$x
    ),
    pca_coordinates$sample_key
  ),
  ,
  drop = FALSE
]


stopifnot(
  identical(
    pca_coordinates$sample_key,
    rownames(
      pca_fit$x
    )
  )
)


############################################################
## 13. SAVE PCA COORDINATES
############################################################

write.table(
  pca_coordinates,
  file.path(
    table_dir,
    "SOLID_PCA_sample_coordinates.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


pca_variance_table <- data.frame(
  PC = paste0(
    "PC",
    seq_along(
      variance_explained
    )
  ),
  variance_explained_pct =
    variance_explained,
  stringsAsFactors = FALSE
)


write.table(
  pca_variance_table,
  file.path(
    table_dir,
    "SOLID_PCA_variance_explained.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
## 14. PCA FIGURE 1
## Main view: sample type
############################################################

p1 <- ggplot(
  pca_coordinates,
  aes(
    x = PC1,
    y = PC2,
    color = sample_type
  )
) +
  geom_point(
    size = 3.5
  ) +
  geom_text(
    aes(
      label = patient_id
    ),
    nudge_y = 0.2,
    check_overlap = TRUE,
    size = 3
  ) +
  labs(
    title = "SOLID tumor–plasma PCA",
    subtitle = paste0(
      "Top ",
      format(
        n_pca_regions,
        big.mark = ","
      ),
      " variable 1-kb regions; M-values"
    ),
    x = paste0(
      "PC1 (",
      round(
        variance_explained[1],
        1
      ),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(
        variance_explained[2],
        1
      ),
      "%)"
    ),
    color = "Sample type"
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
    "SOLID_PCA_by_sample_type.png"
  ),
  p1,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(
    figure_dir,
    "SOLID_PCA_by_sample_type.pdf"
  ),
  p1,
  width = 7,
  height = 6
)


############################################################
## 15. PCA FIGURE 2
## Paired patient trajectories
############################################################

p2 <- ggplot(
  pca_coordinates,
  aes(
    x = PC1,
    y = PC2
  )
) +
  geom_line(
    aes(
      group = patient_id
    ),
    alpha = 0.45
  ) +
  geom_point(
    aes(
      color = sample_type,
      shape = factor(Grade)
    ),
    size = 3.5
  ) +
  labs(
    title = "SOLID paired tumor–plasma PCA",
    subtitle = "Lines connect tumor and plasma samples from the same patient",
    x = paste0(
      "PC1 (",
      round(
        variance_explained[1],
        1
      ),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(
        variance_explained[2],
        1
      ),
      "%)"
    ),
    color = "Sample type",
    shape = "Grade"
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
    "SOLID_PCA_paired_patients.png"
  ),
  p2,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(
    figure_dir,
    "SOLID_PCA_paired_patients.pdf"
  ),
  p2,
  width = 7,
  height = 6
)


############################################################
## 16. PCA FIGURE 3
## Clinical view: Grade + Sex
############################################################

p3 <- ggplot(
  pca_coordinates,
  aes(
    x = PC1,
    y = PC2,
    color = factor(Grade),
    shape = Sex
  )
) +
  geom_point(
    size = 3.5
  ) +
  facet_wrap(
    ~ sample_type
  ) +
  labs(
    title = "SOLID PCA with clinical annotation",
    subtitle = "Tumor and plasma shown separately",
    x = paste0(
      "PC1 (",
      round(
        variance_explained[1],
        1
      ),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(
        variance_explained[2],
        1
      ),
      "%)"
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
    "SOLID_PCA_grade_sex.png"
  ),
  p3,
  width = 9,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(
    figure_dir,
    "SOLID_PCA_grade_sex.pdf"
  ),
  p3,
  width = 9,
  height = 5.5
)


############################################################
## 17. SAMPLE CORRELATION
##
## Spearman correlation over top variable regions.
############################################################

sample_cor <- cor(
  M_pca,
  method = "spearman",
  use = "pairwise.complete.obs"
)


write.table(
  sample_cor,
  file.path(
    table_dir,
    "SOLID_sample_spearman_correlation.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)


############################################################
## 18. HEATMAP ANNOTATION
############################################################

heatmap_annotation <- sample_metadata[
  ,
  c(
    "sample_key",
    "sample_type",
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


rownames(
  heatmap_annotation
) <- heatmap_annotation$sample_key


heatmap_annotation$sample_key <- NULL


############################################################
## Ensure categorical variables are factors
############################################################

factor_fields <- intersect(
  c(
    "sample_type",
    "Grade",
    "Sex",
    "ECOG",
    "Response_RANO",
    "PFS_status",
    "OS_status"
  ),
  names(
    heatmap_annotation
  )
)

for (v in factor_fields) {

  heatmap_annotation[[v]] <- factor(
    heatmap_annotation[[v]]
  )
}


############################################################
## 19. SAMPLE CORRELATION HEATMAP
############################################################

pdf(
  file.path(
    figure_dir,
    "SOLID_sample_correlation_heatmap.pdf"
  ),
  width = 12,
  height = 11
)

pheatmap(
  sample_cor,
  annotation_col =
    heatmap_annotation,
  annotation_row =
    heatmap_annotation,
  clustering_distance_rows =
    as.dist(
      1 - sample_cor
    ),
  clustering_distance_cols =
    as.dist(
      1 - sample_cor
    ),
  clustering_method =
    "complete",
  main =
    "SOLID sample Spearman correlation"
)

dev.off()


png(
  file.path(
    figure_dir,
    "SOLID_sample_correlation_heatmap.png"
  ),
  width = 3200,
  height = 3000,
  res = 300
)

pheatmap(
  sample_cor,
  annotation_col =
    heatmap_annotation,
  annotation_row =
    heatmap_annotation,
  clustering_distance_rows =
    as.dist(
      1 - sample_cor
    ),
  clustering_distance_cols =
    as.dist(
      1 - sample_cor
    ),
  clustering_method =
    "complete",
  main =
    "SOLID sample Spearman correlation"
)

dev.off()


############################################################
## 20. VARIABLE-REGION HEATMAP
##
## Row-standardize M-values so the figure emphasizes
## relative methylation patterns rather than absolute scale.
############################################################

row_mean <- rowMeans(
  M_heatmap
)

row_sd <- apply(
  M_heatmap,
  1,
  sd
)


heatmap_z <- sweep(
  M_heatmap,
  1,
  row_mean,
  "-"
)

heatmap_z <- sweep(
  heatmap_z,
  1,
  row_sd,
  "/"
)


finite_rows <- apply(
  heatmap_z,
  1,
  function(x) {
    all(
      is.finite(x)
    )
  }
)


heatmap_z <- heatmap_z[
  finite_rows,
  ,
  drop = FALSE
]


############################################################
## 21. TOP-VARIABLE-REGION HEATMAP
############################################################

pdf(
  file.path(
    figure_dir,
    "SOLID_top500_variable_regions_heatmap.pdf"
  ),
  width = 11,
  height = 12
)

pheatmap(
  heatmap_z,
  annotation_col =
    heatmap_annotation,
  show_rownames =
    FALSE,
  show_colnames =
    TRUE,
  cluster_rows =
    TRUE,
  cluster_cols =
    TRUE,
  clustering_method =
    "complete",
  main =
    "SOLID top variable 1-kb methylation regions"
)

dev.off()


png(
  file.path(
    figure_dir,
    "SOLID_top500_variable_regions_heatmap.png"
  ),
  width = 3000,
  height = 3400,
  res = 300
)

pheatmap(
  heatmap_z,
  annotation_col =
    heatmap_annotation,
  show_rownames =
    FALSE,
  show_colnames =
    TRUE,
  cluster_rows =
    TRUE,
  cluster_cols =
    TRUE,
  clustering_method =
    "complete",
  main =
    "SOLID top variable 1-kb methylation regions"
)

dev.off()


############################################################
## 22. HIERARCHICAL CLUSTER MEMBERSHIP
##
## Save a simple 2-cluster and 3-cluster solution.
## These are exploratory only.
############################################################

sample_distance <- as.dist(
  1 - sample_cor
)

sample_hclust <- hclust(
  sample_distance,
  method = "complete"
)


cluster_table <- data.frame(
  sample_key =
    names(
      sample_hclust$order
    ),
  stringsAsFactors = FALSE
)


cluster_table <- data.frame(
  sample_key =
    colnames(
      sample_cor
    ),

  cluster_k2 =
    cutree(
      sample_hclust,
      k = 2
    )[
      colnames(
        sample_cor
      )
    ],

  cluster_k3 =
    cutree(
      sample_hclust,
      k = 3
    )[
      colnames(
        sample_cor
      )
    ],

  stringsAsFactors = FALSE
)


cluster_table <- merge(
  cluster_table,
  sample_metadata,
  by = "sample_key",
  all.x = TRUE,
  sort = FALSE
)


cluster_table <- cluster_table[
  match(
    colnames(sample_cor),
    cluster_table$sample_key
  ),
  ,
  drop = FALSE
]


write.table(
  cluster_table,
  file.path(
    table_dir,
    "SOLID_sample_hierarchical_clusters.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
## 23. SAVE PCA OBJECT
############################################################

saveRDS(
  pca_fit,
  file.path(
    output_dir,
    "SOLID_PCA_top10000_variable_M_regions.rds"
  )
)


############################################################
## 24. FINAL SUMMARY
############################################################

cat(
  "\n============================================\n"
)

cat(
  "SOLID PCA + CLUSTERING COMPLETE\n"
)

cat(
  "============================================\n"
)

cat(
  "\nTotal retained regions:",
  nrow(combined_M),
  "\n"
)

cat(
  "Complete regions across all 26 samples:",
  n_complete_regions,
  "\n"
)

cat(
  "Regions used for PCA/correlation:",
  n_pca_regions,
  "\n"
)

cat(
  "Regions used for heatmap:",
  nrow(heatmap_z),
  "\n"
)

cat(
  "\nVariance explained:\n"
)

cat(
  "PC1:",
  round(
    variance_explained[1],
    2
  ),
  "%\n"
)

cat(
  "PC2:",
  round(
    variance_explained[2],
    2
  ),
  "%\n"
)

cat(
  "PC3:",
  round(
    variance_explained[3],
    2
  ),
  "%\n"
)

cat(
  "\nPCA sample coordinates:\n",
  file.path(
    table_dir,
    "SOLID_PCA_sample_coordinates.tsv"
  ),
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
## 25. SESSION INFORMATION
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "sessionInfo_05_SOLID_PCA_clustering.txt"
  )
)
