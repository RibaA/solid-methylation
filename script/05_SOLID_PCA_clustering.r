############################################################
## 05_SOLID_PCA_clustering.r
##
## SOLID MATCHED TUMOR-PLASMA
## PCA + CLUSTERING
##
## Purpose:
##   1. Load frozen downstream analysis-ready object
##   2. Use M-values for multivariate structure
##   3. Apply the existing valid-pair mask
##   4. Select complete regions across all 26 samples
##   5. Rank complete regions by variance
##   6. PCA on top 10,000 variable regions
##   7. Sample correlation + hierarchical clustering
##   8. Heatmap of top 500 variable regions
##   9. Annotate with:
##        Grade
##        Sex
##        ECOG
##        Response_RANO
##
## IMPORTANT:
##   - Run from repository root:
##       C:/solid-methylation
##
##   - No new biological filtering
##   - No imputation
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
    )
  )
}


suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
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


source(
  plot_style_file
)


############################################################
## 3. INPUT / OUTPUT DIRECTORIES
############################################################

input_file <- file.path(
  "result",
  "03_matched_tissue_plasma",
  "SOLID_downstream_analysis_ready.rds"
)


output_dir <- file.path(
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


if (!file.exists(input_file)) {

  stop(
    "Input object not found:\n",
    input_file
  )
}


############################################################
## 4. HELPER FUNCTIONS
############################################################

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


make_annotation_palette <- function(
    values,
    available_colors
) {

  values <- unique(
    as.character(
      values[
        !is.na(values)
      ]
    )
  )

  values <- sort(
    values
  )

  if (length(values) == 0L) {
    return(NULL)
  }


  if (length(values) > length(available_colors)) {

    available_colors <- grDevices::colorRampPalette(
      available_colors
    )(
      length(values)
    )
  }


  setNames(
    available_colors[
      seq_along(values)
    ],
    values
  )
}


############################################################
## 5. LOAD ANALYSIS-READY OBJECT
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
## 6. BASIC INPUT VALIDATION
############################################################

required_components <- c(
  "tissue_M",
  "plasma_M",
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
## 7. KEEP SELECTED CLINICAL ANNOTATIONS ONLY
############################################################

annotation_variables <- c(
  "Grade",
  "Sex",
  "ECOG",
  "Response_RANO"
)


missing_annotation_variables <- setdiff(
  annotation_variables,
  names(patient_metadata)
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


patient_annotation <- patient_metadata[
  ,
  c(
    "patient_id",
    annotation_variables
  ),
  drop = FALSE
]


############################################################
## Convert annotation variables to factors
############################################################

for (v in annotation_variables) {

  patient_annotation[[v]] <- factor(
    patient_annotation[[v]]
  )
}


############################################################
## Display observed annotation levels
############################################################

cat(
  "\nClinical annotation levels:\n"
)


for (v in annotation_variables) {

  cat(
    "\n",
    v,
    ":\n",
    sep = ""
  )

  print(
    levels(
      patient_annotation[[v]]
    )
  )
}


############################################################
## 8. APPLY EXISTING VALID-PAIR MASK
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
## 9. COMBINE TUMOR + PLASMA MATRIX
##
## Rows    = 1-kb regions
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
## 10. BUILD SAMPLE METADATA
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

    patient_annotation[
      ,
      annotation_variables,
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

    patient_annotation[
      ,
      annotation_variables,
      drop = FALSE
    ],

    stringsAsFactors = FALSE
  )
)


sample_metadata$sample_type <- factor(
  sample_metadata$sample_type,
  levels = c(
    "Tumor",
    "Plasma"
  )
)


for (v in annotation_variables) {

  sample_metadata[[v]] <- factor(
    sample_metadata[[v]]
  )
}


stopifnot(
  identical(
    sample_metadata$sample_key,
    colnames(combined_M)
  )
)


############################################################
## Save analysis sample metadata
############################################################

write.table(
  sample_metadata,
  file.path(
    table_dir,
    "SOLID_PCA_sample_metadata.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
## 11. COMPLETE-REGION REQUIREMENT
##
## No imputation.
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
  round(
    100 *
      n_complete_regions /
      nrow(combined_M),
    2
  ),
  "%\n"
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
## 12. REGION VARIANCE
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
## 13. SELECT VARIABLE REGIONS
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
## 14. SAVE VARIABLE-REGION TABLE
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
## 15. PCA
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
## 16. SAVE PCA COORDINATES
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
## Common PCA axis labels
############################################################

pc1_label <- paste0(
  "PC1 (",
  round(
    variance_explained[1],
    1
  ),
  "%)"
)


pc2_label <- paste0(
  "PC2 (",
  round(
    variance_explained[2],
    1
  ),
  "%)"
)


############################################################
## 17. PCA FIGURE 1
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
    size = 3.6,
    alpha = 0.9
  ) +

  scale_color_manual(
    values =
      sample_cols[
        c(
          "Tumor",
          "Plasma"
        )
      ]
  ) +

  labs(
    title =
      "SOLID tumor-plasma PCA",

    subtitle =
      paste0(
        "Top ",
        format(
          n_pca_regions,
          big.mark = ","
        ),
        " variable 1-kb regions; M-values"
      ),

    x =
      pc1_label,

    y =
      pc2_label,

    color =
      "Sample type"
  ) +

  theme_project()


save_plot(
  p1,
  "SOLID_PCA_by_sample_type",
  width = 7,
  height = 6
)


############################################################
## 18. PCA FIGURE 2
## Matched patient trajectories
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
    color = "grey70",
    alpha = 0.60,
    linewidth = 0.6
  ) +

  geom_point(
    aes(
      color = sample_type,
      shape = Grade
    ),
    size = 3.6
  ) +

  scale_color_manual(
    values =
      sample_cols[
        c(
          "Tumor",
          "Plasma"
        )
      ]
  ) +

  labs(
    title =
      "SOLID paired tumor-plasma PCA",

    subtitle =
      "Lines connect matched tumor and plasma samples",

    x =
      pc1_label,

    y =
      pc2_label,

    color =
      "Sample type",

    shape =
      "Grade"
  ) +

  theme_project()


save_plot(
  p2,
  "SOLID_PCA_paired_patients",
  width = 7,
  height = 6
)


############################################################
## 19. CLINICAL ANNOTATION COLORS
##
## Grade is fixed by the shared project palette.
##
## Sex, ECOG and Response_RANO are generated from observed
## values here. Once their exact levels are confirmed, these
## mappings can be moved into 00_plot_style_and_palettes.r.
############################################################

sex_cols <- make_annotation_palette(
  sample_metadata$Sex,
  c(
    "#99B6BD",
    "#ECC9A0",
    "#BEB59C"
  )
)


ecog_cols <- make_annotation_palette(
  sample_metadata$ECOG,
  c(
    "#D2C396",
    "#A5A596",
    "#697878",
    "#4B5A69",
    "#5A4B3C"
  )
)


response_cols <- make_annotation_palette(
  sample_metadata$Response_RANO,
  c(
    "#99B6BD",
    "#B3A86A",
    "#ECC9A0",
    "#D0937D",
    "#78847F",
    "#9A9391"
  )
)


############################################################
## 20. PCA FIGURE 3
## Grade
############################################################

p3 <- ggplot(
  pca_coordinates,
  aes(
    x = PC1,
    y = PC2,
    color = Grade
  )
) +

  geom_point(
    size = 3.6,
    alpha = 0.9
  ) +

  facet_wrap(
    ~ sample_type
  ) +

  scale_color_manual(
    values = grade_cols
  ) +

  labs(
    title =
      "SOLID PCA by grade",

    subtitle =
      "Tumor and plasma shown separately",

    x =
      pc1_label,

    y =
      pc2_label,

    color =
      "Grade"
  ) +

  theme_project()


save_plot(
  p3,
  "SOLID_PCA_by_grade",
  width = 9,
  height = 5.5
)


############################################################
## 21. PCA FIGURE 4
## Sex
############################################################

p4 <- ggplot(
  pca_coordinates,
  aes(
    x = PC1,
    y = PC2,
    color = Sex
  )
) +

  geom_point(
    size = 3.6,
    alpha = 0.9
  ) +

  facet_wrap(
    ~ sample_type
  ) +

  scale_color_manual(
    values = sex_cols
  ) +

  labs(
    title =
      "SOLID PCA by sex",

    subtitle =
      "Tumor and plasma shown separately",

    x =
      pc1_label,

    y =
      pc2_label,

    color =
      "Sex"
  ) +

  theme_project()


save_plot(
  p4,
  "SOLID_PCA_by_sex",
  width = 9,
  height = 5.5
)


############################################################
## 22. PCA FIGURE 5
## ECOG
############################################################

p5 <- ggplot(
  pca_coordinates,
  aes(
    x = PC1,
    y = PC2,
    color = ECOG
  )
) +

  geom_point(
    size = 3.6,
    alpha = 0.9
  ) +

  facet_wrap(
    ~ sample_type
  ) +

  scale_color_manual(
    values = ecog_cols
  ) +

  labs(
    title =
      "SOLID PCA by ECOG",

    subtitle =
      "Tumor and plasma shown separately",

    x =
      pc1_label,

    y =
      pc2_label,

    color =
      "ECOG"
  ) +

  theme_project()


save_plot(
  p5,
  "SOLID_PCA_by_ECOG",
  width = 9,
  height = 5.5
)


############################################################
## 23. PCA FIGURE 6
## RANO response
############################################################

p6 <- ggplot(
  pca_coordinates,
  aes(
    x = PC1,
    y = PC2,
    color = Response_RANO
  )
) +

  geom_point(
    size = 3.6,
    alpha = 0.9
  ) +

  facet_wrap(
    ~ sample_type
  ) +

  scale_color_manual(
    values = response_cols
  ) +

  labs(
    title =
      "SOLID PCA by RANO response",

    subtitle =
      "Tumor and plasma shown separately",

    x =
      pc1_label,

    y =
      pc2_label,

    color =
      "RANO response"
  ) +

  theme_project()


save_plot(
  p6,
  "SOLID_PCA_by_RANO_response",
  width = 9,
  height = 5.5
)


############################################################
## 24. SAMPLE CORRELATION
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
## 25. HEATMAP ANNOTATION
############################################################

heatmap_annotation <- sample_metadata[
  ,
  c(
    "sample_key",
    "sample_type",
    "Grade",
    "Sex",
    "ECOG",
    "Response_RANO"
  ),
  drop = FALSE
]


rownames(
  heatmap_annotation
) <- heatmap_annotation$sample_key


heatmap_annotation$sample_key <- NULL


############################################################
## Drop patient_id intentionally.
## Only selected annotations are shown.
############################################################

for (v in names(heatmap_annotation)) {

  heatmap_annotation[[v]] <- factor(
    heatmap_annotation[[v]]
  )
}


############################################################
## 26. HEATMAP ANNOTATION COLORS
############################################################

grade_heatmap_cols <- grade_cols[
  intersect(
    names(grade_cols),
    levels(
      heatmap_annotation$Grade
    )
  )
]


heatmap_annotation_colors <- list(

  sample_type = c(
    Tumor =
      unname(
        sample_cols["Tumor"]
      ),

    Plasma =
      unname(
        sample_cols["Plasma"]
      )
  ),

  Grade =
    grade_heatmap_cols,

  Sex =
    sex_cols,

  ECOG =
    ecog_cols,

  Response_RANO =
    response_cols
)


############################################################
## Shared heatmap color gradient
############################################################

cor_heatmap_colors <- grDevices::colorRampPalette(
  c(
    "#F4F1EA",
    "#BEB59C",
    "#697878",
    "#4B5A69"
  )
)(
  100
)


methylation_heatmap_colors <- grDevices::colorRampPalette(
  c(
    "#526A83",
    "#BABAAF",
    "#F4F1EA",
    "#D9AF6B",
    "#BF816B"
  )
)(
  101
)


############################################################
## 27. SAMPLE CORRELATION HEATMAP
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

  color =
    cor_heatmap_colors,

  annotation_col =
    heatmap_annotation,

  annotation_row =
    heatmap_annotation,

  annotation_colors =
    heatmap_annotation_colors,

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

  border_color =
    NA,

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

  color =
    cor_heatmap_colors,

  annotation_col =
    heatmap_annotation,

  annotation_row =
    heatmap_annotation,

  annotation_colors =
    heatmap_annotation_colors,

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

  border_color =
    NA,

  main =
    "SOLID sample Spearman correlation"
)


dev.off()


############################################################
## 28. VARIABLE-REGION HEATMAP
##
## Row-standardize M-values to emphasize relative
## methylation patterns rather than absolute M-value scale.
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
## 29. TOP-VARIABLE-REGION HEATMAP
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

  color =
    methylation_heatmap_colors,

  annotation_col =
    heatmap_annotation,

  annotation_colors =
    heatmap_annotation_colors,

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

  border_color =
    NA,

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

  color =
    methylation_heatmap_colors,

  annotation_col =
    heatmap_annotation,

  annotation_colors =
    heatmap_annotation_colors,

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

  border_color =
    NA,

  main =
    "SOLID top variable 1-kb methylation regions"
)


dev.off()


############################################################
## 30. HIERARCHICAL CLUSTER MEMBERSHIP
##
## Save exploratory 2-cluster and 3-cluster solutions.
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
    colnames(
      sample_cor
    ),
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
## 31. SAVE PCA OBJECT
############################################################

saveRDS(
  pca_fit,
  file.path(
    output_dir,
    "SOLID_PCA_top10000_variable_M_regions.rds"
  )
)


############################################################
## 32. FINAL SUMMARY
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
  "\nClinical annotations retained:\n"
)


cat(
  paste(
    annotation_variables,
    collapse = ", "
  ),
  "\n"
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
## 33. SESSION INFORMATION
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "sessionInfo_05_SOLID_PCA_clustering.txt"
  )
)