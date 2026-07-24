##------------------------------------------------------------
## SOLID PLASMA cfDNA 5-BASE METHYLATION
## Exploratory analysis and visualization
##
## Objectives:
##   1. Load and validate processed SOLID methylation data
##   2. Summarize sample- and region-level QC
##   3. Visualize methylation distributions
##   4. Perform PCA and sample correlation analysis
##   5. Prepare SOLID data for tumour DMR projection
##------------------------------------------------------------
##############################################################
## load library
##############################################################
options(
  stringsAsFactors = FALSE,
  warn = 1
)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(matrixStats)
  library(patchwork)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
})

##############################################################
# set up directory
##############################################################
dir_input <- 'data/solid-5base'
dir_output <- 'result/solid-5base'

dir_figures <- file.path(
  dir_output,
  "figures"
)

dir_tables <- file.path(
  dir_output,
  "tables"
)

dir_objects <- file.path(
  dir_output,
  "objects"
)

dir.create(
  dir_figures,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dir_tables,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dir_objects,
  recursive = TRUE,
  showWarnings = FALSE
)

##############################################################
# Define SOLID input files
##############################################################
file_beta <- file.path(
  dir_input,
  "SOLID_1kb_beta_cov5_n7_standard_chr.rds"
)

file_M <- file.path(
  dir_input,
  "SOLID_1kb_M_cov5_n7_standard_chr.rds"
)

file_coverage <- file.path(
  dir_input,
  "SOLID_1kb_coverage_cov5_n7_standard_chr.rds"
)

file_covered_CpGs <- file.path(
  dir_input,
  "SOLID_1kb_covered_CpGs_cov5_n7_standard_chr.rds"
)

file_methylated_reads <- file.path(
  dir_input,
  "SOLID_1kb_methylated_reads_cov5_n7_standard_chr.rds"
)

file_unmethylated_reads <- file.path(
  dir_input,
  "SOLID_1kb_unmethylated_reads_cov5_n7_standard_chr.rds"
)

file_annotation <- file.path(
  dir_input,
  "SOLID_1kb_annotation_cov5_n7_standard_chr.tsv"
)

file_regions <- file.path(
  dir_input,
  "SOLID_1kb_regions_cov5_n7_standard_chr.tsv"
)

file_processing_metadata <- file.path(
  dir_input,
  "SOLID_1kb_processing_metadata.rds"
)

file_sample_qc <- file.path(
  dir_input,
  "SOLID_5base_sample_QC.tsv"
)

files_required <- c(
  beta = file_beta,
  M = file_M,
  coverage = file_coverage,
  covered_CpGs = file_covered_CpGs,
  methylated_reads = file_methylated_reads,
  unmethylated_reads = file_unmethylated_reads,
  annotation = file_annotation,
  regions = file_regions,
  processing_metadata = file_processing_metadata,
  sample_qc = file_sample_qc
)

file_check <- data.frame(
  object = names(files_required),
  file = basename(files_required),
  exists = file.exists(files_required),
  size_MB = round(file.info(files_required)$size / 1024^2, 2),
  row.names = NULL
)

print(file_check)

stopifnot(all(file_check$exists))

##------------------------------------------------------------
## Load SOLID annotation, regions, metadata, and sample QC
##------------------------------------------------------------

annotation_solid <- fread(
  file_annotation
)

regions_solid <- fread(
  file_regions
)

processing_metadata_solid <- readRDS(
  file_processing_metadata
)

qc_solid <- fread(
  file_sample_qc
)

cat("\nAnnotation dimensions:\n")
print(dim(annotation_solid))

cat("\nRegions dimensions:\n")
print(dim(regions_solid))

cat("\nAnnotation column names:\n")
print(names(annotation_solid))

cat("\nRegions column names:\n")
print(names(regions_solid))

cat("\nSample QC dimensions:\n")
print(dim(qc_solid))

cat("\nSample QC column names:\n")
print(names(qc_solid))

cat("\nProcessing metadata:\n")
print(processing_metadata_solid)

anyDuplicated(annotation_solid)
anyDuplicated(regions_solid)

colSums(is.na(annotation_solid))
colSums(is.na(regions_solid))
colSums(is.na(qc_solid))

dim(annotation_solid)
names(annotation_solid)

dim(regions_solid)
names(regions_solid)

dim(qc_solid)
names(qc_solid)

processing_metadata_solid

##------------------------------------------------------------
## 4. Validate regional annotation and sample metadata
##------------------------------------------------------------

setnames(
  regions_solid,
  old = c("V1", "V2", "V3", "V4"),
  new = c("chr", "start", "end", "n_samples_cov5")
)

## Confirm coordinate agreement between the two files
region_coordinate_check <- c(
  chr_identical = identical(
    regions_solid$chr,
    annotation_solid$chr
  ),
  start_identical = identical(
    regions_solid$start,
    annotation_solid$start
  ),
  end_identical = identical(
    regions_solid$end,
    annotation_solid$end
  ),
  n_samples_cov5_identical = identical(
    regions_solid$n_samples_cov5,
    annotation_solid$n_samples_cov5
  )
)

print(region_coordinate_check)

stopifnot(all(region_coordinate_check))

## Confirm QC sample IDs match processing metadata
sample_id_check <- identical(
  qc_solid$sample_id,
  processing_metadata_solid$sample_ids
)

cat(
  "QC sample order matches metadata:",
  sample_id_check,
  "\n"
)

stopifnot(sample_id_check)

summary(annotation_solid$n_samples_cov5)
summary(annotation_solid$n_samples_cov5_matrix)
summary(annotation_solid$n_samples_observed)

table(
  annotation_solid$n_samples_cov5 ==
    annotation_solid$n_samples_cov5_matrix
)

table(
  annotation_solid$n_samples_cov5_matrix ==
    annotation_solid$n_samples_observed
)

##------------------------------------------------------------
## 5. Load SOLID methylation matrices
##------------------------------------------------------------
beta_solid <- readRDS(
  file_beta
)

M_solid <- readRDS(
  file_M
)

coverage_solid <- readRDS(
  file_coverage
)

cat("\nBeta matrix dimensions:\n")
print(dim(beta_solid))

cat("\nM-value matrix dimensions:\n")
print(dim(M_solid))

cat("\nCoverage matrix dimensions:\n")
print(dim(coverage_solid))

##------------------------------------------------------------
## 6. Validate matrix dimensions and sample IDs
##------------------------------------------------------------
stopifnot(
  nrow(beta_solid) == nrow(annotation_solid),
  nrow(M_solid) == nrow(annotation_solid),
  nrow(coverage_solid) == nrow(annotation_solid),
  ncol(beta_solid) == processing_metadata_solid$total_samples,
  ncol(M_solid) == processing_metadata_solid$total_samples,
  ncol(coverage_solid) == processing_metadata_solid$total_samples
)

matrix_column_check <- c(
  beta_vs_M = identical(
    colnames(beta_solid),
    colnames(M_solid)
  ),
  beta_vs_coverage = identical(
    colnames(beta_solid),
    colnames(coverage_solid)
  ),
  beta_vs_metadata = identical(
    colnames(beta_solid),
    processing_metadata_solid$sample_ids
  ),
  beta_vs_qc = identical(
    colnames(beta_solid),
    qc_solid$sample_id
  )
)

print(matrix_column_check)

stopifnot(all(matrix_column_check))

##------------------------------------------------------------
## 7. Check missingness and value ranges
##------------------------------------------------------------

beta_missing_by_sample <- colSums(
  is.na(beta_solid)
)

M_missing_by_sample <- colSums(
  is.na(M_solid)
)

coverage_missing_by_sample <- colSums(
  is.na(coverage_solid)
)

missingness_summary <- data.table(
  sample_id = colnames(beta_solid),
  beta_missing_n = beta_missing_by_sample,
  beta_missing_pct = 100 * beta_missing_by_sample / nrow(beta_solid),
  M_missing_n = M_missing_by_sample,
  M_missing_pct = 100 * M_missing_by_sample / nrow(M_solid),
  coverage_missing_n = coverage_missing_by_sample,
  coverage_missing_pct =
    100 * coverage_missing_by_sample / nrow(coverage_solid)
)

print(missingness_summary)

beta_range <- range(
  beta_solid,
  na.rm = TRUE
)

M_range <- range(
  M_solid,
  na.rm = TRUE
)

coverage_range <- range(
  coverage_solid,
  na.rm = TRUE
)

cat("\nBeta range:\n")
print(beta_range)

cat("\nM-value range:\n")
print(M_range)

cat("\nCoverage range:\n")
print(coverage_range)

summary(qc_solid$weighted_global_beta)
summary(qc_solid$mean_cpg_coverage)

missingness_summary
matrix_column_check

dim(beta_solid)
dim(M_solid)
dim(coverage_solid)

##------------------------------------------------------------
## 8. Create the SOLID sample-level QC summary
##------------------------------------------------------------

qc_summary_solid <- merge(
  qc_solid,
  missingness_summary[
    ,
    .(
      sample_id,
      beta_missing_n,
      beta_missing_pct
    )
  ],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

## Restore matrix/sample order after merging
qc_summary_solid <- qc_summary_solid[
  match(
    colnames(beta_solid),
    sample_id
  )
]

qc_summary_solid[
  ,
  total_reads :=
    total_methylated_reads +
    total_unmethylated_reads
]

qc_summary_solid[
  ,
  methylated_read_pct :=
    100 * total_methylated_reads / total_reads
]

qc_summary_solid[
  ,
  observed_regions :=
    nrow(beta_solid) - beta_missing_n
]

qc_summary_solid[
  ,
  observed_regions_pct :=
    100 * observed_regions / nrow(beta_solid)
]

## Shorter labels for figures
qc_summary_solid[
  ,
  sample_label := sub(
    "_01_LB01-01$",
    "",
    sample_id
  )
]

print(qc_summary_solid)

stopifnot(
  identical(
    qc_summary_solid$sample_id,
    colnames(beta_solid)
  ),
  !anyNA(qc_summary_solid)
)

fwrite(
  qc_summary_solid,
  file.path(
    dir_tables,
    "SOLID_sample_QC_summary.tsv"
  ),
  sep = "\t"
)

##------------------------------------------------------------
## 9. SOLID sample-level QC visualization
##------------------------------------------------------------

qc_summary_solid[
  ,
  sample_label := factor(
    sample_label,
    levels = sample_label
  )
]

theme_solid <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(
      face = "bold",
      size = 11
    )
  )

p_total_reads <- ggplot(
  qc_summary_solid,
  aes(
    x = sample_label,
    y = total_reads
  )
) +
  geom_col() +
  scale_y_continuous(
    labels = scales::label_number(
      scale_cut = scales::cut_short_scale()
    ),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "A. Total CpG read coverage",
    x = NULL,
    y = "Methylated + unmethylated reads"
  ) +
  theme_solid

p_mean_coverage <- ggplot(
  qc_summary_solid,
  aes(
    x = sample_label,
    y = mean_cpg_coverage
  )
) +
  geom_col() +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "B. Mean CpG coverage",
    x = NULL,
    y = "Mean coverage"
  ) +
  theme_solid

p_global_beta <- ggplot(
  qc_summary_solid,
  aes(
    x = sample_label,
    y = weighted_global_beta
  )
) +
  geom_col() +
  coord_cartesian(
    ylim = c(
      min(qc_summary_solid$weighted_global_beta) - 0.02,
      max(qc_summary_solid$weighted_global_beta) + 0.02
    )
  ) +
  labs(
    title = "C. Weighted global methylation",
    x = NULL,
    y = "Weighted global Beta"
  ) +
  theme_solid

p_observed_regions <- ggplot(
  qc_summary_solid,
  aes(
    x = sample_label,
    y = observed_regions_pct
  )
) +
  geom_col() +
  coord_cartesian(
    ylim = c(
      min(qc_summary_solid$observed_regions_pct) - 1,
      100
    )
  ) +
  labs(
    title = "D. Regional data completeness",
    x = NULL,
    y = "Observed regions (%)"
  ) +
  theme_solid

p_qc_solid <- (
  p_total_reads |
    p_mean_coverage
) / (
  p_global_beta |
    p_observed_regions
)

p_qc_solid

ggsave(
  filename = file.path(
    dir_figures,
    "SOLID_sample_QC_overview.png"
  ),
  plot = p_qc_solid,
  width = 12,
  height = 8,
  dpi = 300
)

ggsave(
  filename = file.path(
    dir_figures,
    "SOLID_sample_QC_overview.pdf"
  ),
  plot = p_qc_solid,
  width = 12,
  height = 8
)

##------------------------------------------------------------
## 10. SOLID regional Beta-value distributions
##------------------------------------------------------------

set.seed(12345)

n_plot_regions <- 100000L

plot_region_index <- sample(
  seq_len(nrow(beta_solid)),
  size = min(n_plot_regions, nrow(beta_solid)),
  replace = FALSE
)

beta_plot_matrix <- beta_solid[
  plot_region_index,
  ,
  drop = FALSE
]

beta_plot_long <- as.data.table(
  as.table(beta_plot_matrix)
)

setnames(
  beta_plot_long,
  c("region_index", "sample_id", "beta")
)

beta_plot_long[
  ,
  sample_label := sub(
    "_01_LB01-01$",
    "",
    sample_id
  )
]

beta_plot_long[
  ,
  sample_label := factor(
    sample_label,
    levels = qc_summary_solid$sample_label
  )
]

beta_plot_long <- beta_plot_long[
  !is.na(beta)
]

dim(beta_plot_long)
summary(beta_plot_long$beta)

p_beta_density <- ggplot(
  beta_plot_long,
  aes(
    x = beta,
    group = sample_label
  )
) +
  geom_density(
    linewidth = 0.6,
    alpha = 0.15
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  labs(
    title = "SOLID plasma cfDNA regional methylation distributions",
    subtitle = paste0(
      "Random subset of ",
      format(n_plot_regions, big.mark = ","),
      " retained 1-kb regions"
    ),
    x = "Regional Beta value",
    y = "Density"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

p_beta_density

p_beta_density_facet <- ggplot(
  beta_plot_long,
  aes(x = beta)
) +
  geom_density(
    linewidth = 0.5,
    fill = "grey70"
  ) +
  facet_wrap(
    ~ sample_label,
    ncol = 4
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1)
  ) +
  labs(
    title = "SOLID plasma cfDNA Beta-value distributions",
    subtitle = paste0(
      "Random subset of ",
      format(n_plot_regions, big.mark = ","),
      " retained 1-kb regions"
    ),
    x = "Regional Beta value",
    y = "Density"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    strip.text = element_text(size = 9)
  )

p_beta_density_facet

beta_summary_solid <- data.table(
  sample_id = colnames(beta_solid),
  median_beta = matrixStats::colMedians(
    beta_solid,
    na.rm = TRUE
  ),
  mean_beta = colMeans(
    beta_solid,
    na.rm = TRUE
  ),
  proportion_hypomethylated = colMeans(
    beta_solid < 0.20,
    na.rm = TRUE
  ),
  proportion_intermediate = colMeans(
    beta_solid >= 0.20 & beta_solid <= 0.80,
    na.rm = TRUE
  ),
  proportion_hypermethylated = colMeans(
    beta_solid > 0.80,
    na.rm = TRUE
  )
)

beta_summary_solid[
  ,
  sample_label := sub(
    "_01_LB01-01$",
    "",
    sample_id
  )
]

print(beta_summary_solid)


ggsave(
  file.path(
    dir_figures,
    "SOLID_beta_density_overlay.png"
  ),
  p_beta_density,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(
    dir_figures,
    "SOLID_beta_density_faceted.png"
  ),
  p_beta_density_facet,
  width = 11,
  height = 9,
  dpi = 300
)

fwrite(
  beta_summary_solid,
  file.path(
    dir_tables,
    "SOLID_beta_distribution_summary.tsv"
  ),
  sep = "\t"
)

##------------------------------------------------------------
## 11. SOLID plasma cfDNA PCA
## Complete regions and top-variable M values
##------------------------------------------------------------
## Identify regions measured in all 14 samples
complete_region_index <- which(
  rowSums(is.na(M_solid)) == 0L
)

cat(
  "Complete regions:",
  format(length(complete_region_index), big.mark = ","),
  "\n"
)

cat(
  "Percentage complete:",
  round(
    100 * length(complete_region_index) / nrow(M_solid),
    2
  ),
  "%\n"
)

M_complete <- M_solid[
  complete_region_index,
  ,
  drop = FALSE
]

region_variance <- matrixStats::rowVars(
  M_complete
)

summary(region_variance)

stopifnot(
  length(region_variance) == nrow(M_complete),
  all(is.finite(region_variance))
)

n_variable_regions <- 10000L

top_variable_index_complete <- order(
  region_variance,
  decreasing = TRUE
)[seq_len(
  min(n_variable_regions, length(region_variance))
)]

top_variable_index_original <-
  complete_region_index[top_variable_index_complete]

M_pca <- M_solid[
  top_variable_index_original,
  ,
  drop = FALSE
]

dim(M_pca)
anyNA(M_pca)

pca_solid <- prcomp(
  t(M_pca),
  center = TRUE,
  scale. = FALSE
)

variance_explained <- 100 *
  pca_solid$sdev^2 /
  sum(pca_solid$sdev^2)

pca_solid_df <- data.table(
  sample_id = rownames(pca_solid$x),
  PC1 = pca_solid$x[, "PC1"],
  PC2 = pca_solid$x[, "PC2"],
  PC3 = pca_solid$x[, "PC3"]
)

pca_solid_df[
  ,
  sample_label := sub(
    "_01_LB01-01$",
    "",
    sample_id
  )
]

## Add QC metrics for interpretation
pca_solid_df <- merge(
  pca_solid_df,
  qc_summary_solid[
    ,
    .(
      sample_id,
      mean_cpg_coverage,
      weighted_global_beta,
      beta_missing_pct
    )
  ],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

pca_solid_df <- pca_solid_df[
  match(
    rownames(pca_solid$x),
    sample_id
  )
]

print(pca_solid_df)
print(round(variance_explained[1:5], 2))

p_pca_solid <- ggplot(
  pca_solid_df,
  aes(
    x = PC1,
    y = PC2,
    label = sample_label
  )
) +
  geom_point(
    size = 3
  ) +
  geom_text_repel(
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  labs(
    title = "SOLID plasma cfDNA methylation PCA",
    subtitle = paste0(
      "Top ",
      format(nrow(M_pca), big.mark = ","),
      " most variable complete 1-kb regions"
    ),
    x = paste0(
      "PC1 (",
      round(variance_explained[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(variance_explained[2], 1),
      "%)"
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

p_pca_solid

ggsave(
  file.path(
    dir_figures,
    "SOLID_PCA_top10000_variable_regions.png"
  ),
  p_pca_solid,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(
    dir_figures,
    "SOLID_PCA_top10000_variable_regions.pdf"
  ),
  p_pca_solid,
  width = 8,
  height = 6
)

fwrite(
  pca_solid_df,
  file.path(
    dir_tables,
    "SOLID_PCA_sample_coordinates.tsv"
  ),
  sep = "\t"
)

##------------------------------------------------------------
## 12. Assess technical associations with principal components
##------------------------------------------------------------

pc_qc_correlations <- rbindlist(
  lapply(
    c("PC1", "PC2", "PC3"),
    function(pc_name) {
      rbindlist(
        lapply(
          c(
            "mean_cpg_coverage",
            "weighted_global_beta",
            "beta_missing_pct"
          ),
          function(qc_name) {
            test_result <- cor.test(
              pca_solid_df[[pc_name]],
              pca_solid_df[[qc_name]],
              method = "spearman",
              exact = FALSE
            )

            data.table(
              principal_component = pc_name,
              qc_metric = qc_name,
              spearman_rho = unname(test_result$estimate),
              p_value = test_result$p.value
            )
          }
        )
      )
    }
  )
)

pc_qc_correlations[
  ,
  FDR := p.adjust(
    p_value,
    method = "BH"
  )
]

print(pc_qc_correlations)

fwrite(
  pc_qc_correlations,
  file.path(
    dir_tables,
    "SOLID_PCA_QC_correlations.tsv"
  ),
  sep = "\t"
)

##------------------------------------------------------------
## 13. PCA coloured by technical metrics
##------------------------------------------------------------

p_pca_coverage <- ggplot(
  pca_solid_df,
  aes(
    x = PC1,
    y = PC2,
    colour = mean_cpg_coverage,
    label = sample_label
  )
) +
  geom_point(size = 3.5) +
  geom_text_repel(
    size = 3.2,
    max.overlaps = Inf
  ) +
  labs(
    title = "A. PCA coloured by mean CpG coverage",
    x = paste0(
      "PC1 (",
      round(variance_explained[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(variance_explained[2], 1),
      "%)"
    ),
    colour = "Mean coverage"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

p_pca_missingness <- ggplot(
  pca_solid_df,
  aes(
    x = PC1,
    y = PC2,
    colour = beta_missing_pct,
    label = sample_label
  )
) +
  geom_point(size = 3.5) +
  geom_text_repel(
    size = 3.2,
    max.overlaps = Inf
  ) +
  labs(
    title = "B. PCA coloured by regional missingness",
    x = paste0(
      "PC1 (",
      round(variance_explained[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(variance_explained[2], 1),
      "%)"
    ),
    colour = "Missing regions (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

p_pca_global_beta <- ggplot(
  pca_solid_df,
  aes(
    x = PC1,
    y = PC2,
    colour = weighted_global_beta,
    label = sample_label
  )
) +
  geom_point(size = 3.5) +
  geom_text_repel(
    size = 3.2,
    max.overlaps = Inf
  ) +
  labs(
    title = "C. PCA coloured by global methylation",
    x = paste0(
      "PC1 (",
      round(variance_explained[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(variance_explained[2], 1),
      "%)"
    ),
    colour = "Weighted global Beta"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

p_pca_qc <- (
  p_pca_coverage |
    p_pca_missingness
) / p_pca_global_beta

p_pca_qc

ggsave(
  file.path(
    dir_figures,
    "SOLID_PCA_technical_assessment.png"
  ),
  p_pca_qc,
  width = 12,
  height = 10,
  dpi = 300
)

##------------------------------------------------------------
## 14. SOLID sample correlation heatmap
##------------------------------------------------------------

cor_solid <- cor(
  M_pca,
  method = "spearman",
  use = "pairwise.complete.obs"
)

stopifnot(
  all(dim(cor_solid) == c(14, 14)),
  all(is.finite(cor_solid))
)

short_sample_names <- sub(
  "_01_LB01-01$",
  "",
  colnames(cor_solid)
)

rownames(cor_solid) <- short_sample_names
colnames(cor_solid) <- short_sample_names

summary(cor_solid[upper.tri(cor_solid)])

annotation_heatmap <- HeatmapAnnotation(
  mean_coverage = qc_summary_solid$mean_cpg_coverage,
  missing_pct = qc_summary_solid$beta_missing_pct,
  global_beta = qc_summary_solid$weighted_global_beta,
  annotation_name_side = "left"
)

p_cor_solid <- Heatmap(
  cor_solid,
  name = "Spearman\ncorrelation",
  top_annotation = annotation_heatmap,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  row_names_gp = grid::gpar(fontsize = 9),
  column_names_gp = grid::gpar(fontsize = 9),
  column_names_rot = 45,
  border = TRUE,
  column_title = paste0(
    "SOLID plasma cfDNA sample correlation\n",
    "Top ",
    format(nrow(M_pca), big.mark = ","),
    " variable complete 1-kb regions"
  )
)

draw(
  p_cor_solid,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

pdf(
  file.path(
    dir_figures,
    "SOLID_sample_correlation_heatmap.pdf"
  ),
  width = 10,
  height = 9
)

draw(
  p_cor_solid,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

png(
  file.path(
    dir_figures,
    "SOLID_sample_correlation_heatmap.png"
  ),
  width = 3000,
  height = 2700,
  res = 300
)

draw(
  p_cor_solid,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

fwrite(
  as.data.table(
    cor_solid,
    keep.rownames = "sample_id"
  ),
  file.path(
    dir_tables,
    "SOLID_sample_spearman_correlation.tsv"
  ),
  sep = "\t"
)
