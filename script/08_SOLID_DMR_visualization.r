############################################################
## 08_SOLID_DMR_visualization.r
##
## SOLID DMR VISUALIZATION + IMPORTANT-REGION PRIORITIZATION
##
## Input:
##   result/03_matched_tissue_plasma/DMR/
##     SOLID_tissue_vs_plasma_DMR_object.rds
##
## Main goals:
##
##   1. Visualize paired tumor-plasma DMR results
##
##   2. Prioritize biologically compelling regions using:
##        - FDR
##        - median Delta Beta
##        - direction consistency across patients
##        - recurrent effect consistency
##        - valid-pair support
##        - EPIC / plasma technical support when available
##
##   3. Generate:
##        - Volcano plot with labeled important regions
##        - Effect-size sensitivity summary
##        - Delta-Beta distribution
##        - Direction summary
##        - Effect vs consistency scatter
##        - Delta-Beta heatmap
##        - Tumor + Plasma Beta heatmap
##        - Paired Beta plots for representative regions
##        - Genome-wide effect-size plot
##        - High-confidence adjacent DMR blocks
##
## IMPORTANT:
##
## Script 07 primary DMR:
##   FDR < 0.05
##   |median Delta Beta| >= 0.05
##
## Script 08 visualization/prioritization tiers:
##
## STRONG CONSISTENT:
##   FDR < 0.05
##   |median Delta Beta| >= 0.20
##   direction consistency >= 80%
##
## MODERATE CONSISTENT:
##   FDR < 0.05
##   |median Delta Beta| >= 0.10
##   direction consistency >= 80%
##
## PRIMARY DMR:
##   FDR < 0.05
##   |median Delta Beta| >= 0.05
##
## Delta Beta:
##   Plasma - Tumor
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
## 1. USER SETTINGS
############################################################

project_dir <- "C:/solid-methylation"

fdr_threshold <- 0.05

primary_delta_beta_threshold <- 0.05

moderate_delta_beta_threshold <- 0.10

strong_delta_beta_threshold <- 0.20

minimum_direction_consistency <- 80

############################################################
## Display settings
############################################################

n_volcano_labels_per_direction <- 8L

n_heatmap_regions_per_direction <- 20L

n_paired_plot_regions_per_direction <- 4L

n_top_table_per_direction <- 1000L

n_top_important_regions <- 200L


############################################################
## Heatmap display
############################################################

heatmap_clip_delta_beta <- 0.50


############################################################
## Genomic block definition
##
## Fixed 1-kb windows are 0-based half-open.
##
## gap = 0 means only directly adjacent significant
## windows are merged.
############################################################

maximum_gap_between_windows <- 0L

minimum_windows_per_block <- 2L


############################################################
## Volcano display
############################################################

maximum_minus_log10_FDR_display <- 50


############################################################
## Reproducibility
############################################################

random_seed <- 20260912L

set.seed(
  random_seed
)


############################################################
## 2. REQUIRED PACKAGES
############################################################

cran_packages <- c(
  "data.table",
  "ggplot2",
  "matrixStats",
  "scales",
  "ggrepel",
  "circlize"
)

bioc_packages <- c(
  "ComplexHeatmap"
)


missing_cran <- cran_packages[
  !vapply(
    cran_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

missing_bioc <- bioc_packages[
  !vapply(
    bioc_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]


if (length(missing_cran) > 0L) {

  stop(
    "Missing CRAN package(s): ",
    paste(
      missing_cran,
      collapse = ", "
    ),
    "\nInstall with:\ninstall.packages(c(",
    paste(
      sprintf(
        '"%s"',
        missing_cran
      ),
      collapse = ", "
    ),
    "))"
  )
}


if (length(missing_bioc) > 0L) {

  stop(
    "Missing Bioconductor package(s): ",
    paste(
      missing_bioc,
      collapse = ", "
    ),
    "\nInstall using BiocManager::install()."
  )
}


suppressPackageStartupMessages({

  library(data.table)

  library(ggplot2)

  library(matrixStats)

  library(scales)

  library(ggrepel)

  library(circlize)

  library(ComplexHeatmap)
})


############################################################
## 3. DIRECTORIES
############################################################

dmr_dir <- file.path(
  project_dir,
  "result",
  "03_matched_tissue_plasma",
  "DMR"
)


figure_dir <- file.path(
  dmr_dir,
  "figures"
)


priority_dir <- file.path(
  dmr_dir,
  "prioritization"
)


dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  priority_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dmr_object_file <- file.path(
  dmr_dir,
  "SOLID_tissue_vs_plasma_DMR_object.rds"
)


if (!file.exists(dmr_object_file)) {

  stop(
    "DMR object not found:\n",
    dmr_object_file,
    "\nRun Script 07 first."
  )
}


############################################################
## 4. HELPER FUNCTIONS
############################################################

message_header <- function(text) {

  cat(
    "\n",
    paste(
      rep(
        "=",
        72
      ),
      collapse = ""
    ),
    "\n",
    text,
    "\n",
    paste(
      rep(
        "=",
        72
      ),
      collapse = ""
    ),
    "\n",
    sep = ""
  )
}


save_plot <- function(
    plot_object,
    filename,
    width,
    height,
    dpi = 300
) {

  ggsave(
    filename = filename,
    plot = plot_object,
    width = width,
    height = height,
    dpi = dpi,
    limitsize = FALSE
  )

  invisible(
    filename
  )
}


standardize_chr <- function(x) {

  x <- as.character(x)

  x <- sub(
    "^chr",
    "",
    x,
    ignore.case = TRUE
  )


  factor(
    x,
    levels = c(
      as.character(
        1:22
      ),
      "X",
      "Y"
    ),
    ordered = TRUE
  )
}


############################################################
## 5. LOAD DMR OBJECT
############################################################

message_header(
  "LOADING SOLID DMR OBJECT"
)


dmr_object <- readRDS(
  dmr_object_file
)


required_elements <- c(
  "all_results",
  "delta_beta",
  "tissue_beta",
  "plasma_beta",
  "valid_pair_mask",
  "patient_metadata"
)


missing_elements <- setdiff(
  required_elements,
  names(
    dmr_object
  )
)


if (length(missing_elements) > 0L) {

  stop(
    "DMR object missing element(s): ",
    paste(
      missing_elements,
      collapse = ", "
    )
  )
}


results <- as.data.table(
  dmr_object$all_results
)


delta_beta <- dmr_object$delta_beta

tissue_beta <- dmr_object$tissue_beta

plasma_beta <- dmr_object$plasma_beta

valid_pair_mask <- dmr_object$valid_pair_mask


patient_metadata <- as.data.table(
  dmr_object$patient_metadata
)


############################################################
## 6. VALIDATE INPUT
############################################################

stopifnot(

  nrow(results) > 0L,

  is.matrix(
    delta_beta
  ),

  is.matrix(
    tissue_beta
  ),

  is.matrix(
    plasma_beta
  ),

  is.matrix(
    valid_pair_mask
  ),

  identical(
    dim(delta_beta),
    dim(tissue_beta)
  ),

  identical(
    dim(delta_beta),
    dim(plasma_beta)
  ),

  identical(
    dim(delta_beta),
    dim(valid_pair_mask)
  ),

  ncol(
    delta_beta
  ) == 13L,

  nrow(
    patient_metadata
  ) == 13L,

  setequal(
    results$region_id,
    rownames(
      delta_beta
    )
  )
)


############################################################
## Restore matrix order
############################################################

results <- results[
  match(
    rownames(
      delta_beta
    ),
    region_id
  )
]


stopifnot(
  identical(
    results$region_id,
    rownames(
      delta_beta
    )
  )
)


patient_ids <- colnames(
  delta_beta
)


stopifnot(
  identical(
    patient_ids,
    patient_metadata$patient_id
  )
)


cat(
  "Regions loaded:",
  format(
    nrow(results),
    big.mark = ","
  ),
  "\n"
)


cat(
  "Matched patients:",
  length(
    patient_ids
  ),
  "\n"
)


############################################################
## 7. CALCULATE PATIENT-LEVEL DIRECTION CONSISTENCY
############################################################

message_header(
  "CALCULATING REGION-LEVEL DIRECTION CONSISTENCY"
)


############################################################
## Valid observations
############################################################

valid_delta <- is.finite(
  delta_beta
)


n_valid <- rowSums(
  valid_delta
)


############################################################
## Sign consistency
############################################################

n_positive <- rowSums(
  delta_beta > 0,
  na.rm = TRUE
)


n_negative <- rowSums(
  delta_beta < 0,
  na.rm = TRUE
)


############################################################
## Strong-effect recurrence
############################################################

n_positive_010 <- rowSums(
  delta_beta >=
    moderate_delta_beta_threshold,
  na.rm = TRUE
)


n_negative_010 <- rowSums(
  delta_beta <=
    -moderate_delta_beta_threshold,
  na.rm = TRUE
)


n_positive_020 <- rowSums(
  delta_beta >=
    strong_delta_beta_threshold,
  na.rm = TRUE
)


n_negative_020 <- rowSums(
  delta_beta <=
    -strong_delta_beta_threshold,
  na.rm = TRUE
)


############################################################
## Direction-specific consistency
############################################################

results[
  ,
  direction_consistency_pct :=
    fifelse(
      median_delta_beta > 0,

      100 *
        n_positive /
        n_valid,

      fifelse(
        median_delta_beta < 0,

        100 *
          n_negative /
          n_valid,

        NA_real_
      )
    )
]


results[
  ,
  effect_consistency_0.10_pct :=
    fifelse(
      median_delta_beta > 0,

      100 *
        n_positive_010 /
        n_valid,

      fifelse(
        median_delta_beta < 0,

        100 *
          n_negative_010 /
          n_valid,

        NA_real_
      )
    )
]


results[
  ,
  effect_consistency_0.20_pct :=
    fifelse(
      median_delta_beta > 0,

      100 *
        n_positive_020 /
        n_valid,

      fifelse(
        median_delta_beta < 0,

        100 *
          n_negative_020 /
          n_valid,

        NA_real_
      )
    )
]


############################################################
## Ensure valid-pair count remains consistent
############################################################

results[
  ,
  n_valid_pairs_visualization :=
    n_valid
]


############################################################
## 8. DEFINE REGION DIRECTION
############################################################

results[
  ,
  DMR_direction :=
    fifelse(
      median_delta_beta > 0,
      "Plasma_higher",
      fifelse(
        median_delta_beta < 0,
        "Tumor_higher",
        "No_change"
      )
    )
]


############################################################
## 9. DEFINE PRIORITIZATION TIERS
##
## Do NOT replace Script 07's primary statistical definition.
## These tiers are for visualization and prioritization.
############################################################

message_header(
  "DEFINING IMPORTANT-REGION TIERS"
)


results[
  ,
  importance_tier :=
    fifelse(

      FDR < fdr_threshold &

        abs_median_delta_beta >=
          strong_delta_beta_threshold &

        direction_consistency_pct >=
          minimum_direction_consistency,

      "Strong_consistent",

      fifelse(

        FDR < fdr_threshold &

          abs_median_delta_beta >=
            moderate_delta_beta_threshold &

          direction_consistency_pct >=
            minimum_direction_consistency,

        "Moderate_consistent",

        fifelse(

          FDR < fdr_threshold &

            abs_median_delta_beta >=
              primary_delta_beta_threshold,

          "Primary_DMR",

          fifelse(
            FDR <
              fdr_threshold,
            "FDR_only",
            "Not_significant"
          )
        )
      )
    )
]


results[
  ,
  importance_tier :=
    factor(
      importance_tier,
      levels = c(
        "Strong_consistent",
        "Moderate_consistent",
        "Primary_DMR",
        "FDR_only",
        "Not_significant"
      ),
      ordered = TRUE
    )
]


############################################################
## 10. IMPORTANT REGION TABLE
############################################################

strong_regions <- results[
  importance_tier ==
    "Strong_consistent"
]


moderate_regions <- results[
  importance_tier %in%
    c(
      "Strong_consistent",
      "Moderate_consistent"
    )
]


############################################################
## Transparent ranking:
##
## 1. stronger recurrent effect across patients
## 2. larger absolute median Delta Beta
## 3. stronger direction consistency
## 4. smaller FDR
## 5. more valid pairs
############################################################

setorder(
  strong_regions,
  -effect_consistency_0.10_pct,
  -abs_median_delta_beta,
  -direction_consistency_pct,
  FDR,
  -n_valid_pairs
)


important_regions <- head(
  strong_regions,
  n_top_important_regions
)


############################################################
## Split direction
############################################################

important_tumor <- strong_regions[
  DMR_direction ==
    "Tumor_higher"
]


important_plasma <- strong_regions[
  DMR_direction ==
    "Plasma_higher"
]


############################################################
## 11. SUMMARY OF PRIORITIZATION TIERS
############################################################

tier_summary <- results[
  ,
  .(
    n_regions = .N
  ),
  by = .(
    importance_tier,
    DMR_direction
  )
]


tier_summary[
  ,
  percentage_tested :=
    100 *
    n_regions /
    nrow(
      results
    )
]


print(
  tier_summary
)


fwrite(
  tier_summary,
  file.path(
    priority_dir,
    "SOLID_DMR_importance_tier_summary.tsv"
  ),
  sep = "\t"
)


############################################################
## 12. SAVE PRIORITIZED REGION TABLES
############################################################

fwrite(
  results,
  file.path(
    priority_dir,
    "SOLID_all_DMR_results_with_consistency.tsv.gz"
  ),
  sep = "\t",
  compress = "gzip"
)


fwrite(
  strong_regions,
  file.path(
    priority_dir,
    "SOLID_strong_consistent_DMRs.tsv.gz"
  ),
  sep = "\t",
  compress = "gzip"
)


fwrite(
  important_regions,
  file.path(
    priority_dir,
    "SOLID_top200_important_DMR_regions.tsv"
  ),
  sep = "\t"
)


fwrite(
  head(
    important_tumor,
    n_top_table_per_direction
  ),
  file.path(
    priority_dir,
    "SOLID_top_tumor_higher_DMRs.tsv"
  ),
  sep = "\t"
)


fwrite(
  head(
    important_plasma,
    n_top_table_per_direction
  ),
  file.path(
    priority_dir,
    "SOLID_top_plasma_higher_DMRs.tsv"
  ),
  sep = "\t"
)


############################################################
## 13. EFFECT-SIZE / CONSISTENCY SENSITIVITY TABLE
############################################################

effect_thresholds <- c(
  0.05,
  0.10,
  0.20
)


consistency_thresholds <- c(
  50,
  70,
  80,
  90,
  100
)


sensitivity_table <- rbindlist(
  lapply(
    effect_thresholds,
    function(effect_cutoff) {

      rbindlist(
        lapply(
          consistency_thresholds,
          function(consistency_cutoff) {

            keep <-
              results$FDR <
                fdr_threshold &

              results$abs_median_delta_beta >=
                effect_cutoff &

              results$direction_consistency_pct >=
                consistency_cutoff


            data.table(

              abs_median_delta_beta =
                effect_cutoff,

              direction_consistency_pct =
                consistency_cutoff,

              n_regions =
                sum(
                  keep,
                  na.rm = TRUE
                ),

              tumor_higher =
                sum(
                  keep &
                    results$DMR_direction ==
                    "Tumor_higher",
                  na.rm = TRUE
                ),

              plasma_higher =
                sum(
                  keep &
                    results$DMR_direction ==
                    "Plasma_higher",
                  na.rm = TRUE
                )
            )
          }
        )
      )
    }
  )
)


fwrite(
  sensitivity_table,
  file.path(
    priority_dir,
    "SOLID_DMR_effect_consistency_sensitivity.tsv"
  ),
  sep = "\t"
)


cat(
  "\nEffect + consistency sensitivity:\n"
)

print(
  sensitivity_table
)


############################################################
## 14. VOLCANO PLOT WITH TOP IMPORTANT REGIONS LABELED
############################################################

message_header(
  "CREATING LABELED VOLCANO PLOT"
)


results[
  ,
  minus_log10_FDR :=
    -log10(
      pmax(
        FDR,
        .Machine$double.xmin
      )
    )
]


results[
  ,
  minus_log10_FDR_display :=
    pmin(
      minus_log10_FDR,
      maximum_minus_log10_FDR_display
    )
]


############################################################
## Select labels separately by direction
############################################################

volcano_labels <- rbindlist(
  list(

    head(
      important_tumor,
      n_volcano_labels_per_direction
    ),

    head(
      important_plasma,
      n_volcano_labels_per_direction
    )
  ),
  use.names = TRUE,
  fill = TRUE
)


volcano_labels[
  ,
  minus_log10_FDR_display :=
    pmin(
      -log10(
        pmax(
          FDR,
          .Machine$double.xmin
        )
      ),
      maximum_minus_log10_FDR_display
    )
]


volcano_plot <- ggplot(
  results,
  aes(
    x = median_delta_beta,
    y = minus_log10_FDR_display
  )
) +

  geom_point(
    data = results[
      importance_tier %in%
        c(
          "Not_significant",
          "FDR_only"
        )
    ],
    alpha = 0.15,
    size = 0.35,
    colour = "grey70"
  ) +

  geom_point(
    data = results[
      importance_tier ==
        "Primary_DMR"
    ],
    alpha = 0.20,
    size = 0.40,
    colour = "grey45"
  ) +

  geom_point(
    data = results[
      importance_tier ==
        "Moderate_consistent"
    ],
    aes(
      colour = DMR_direction
    ),
    alpha = 0.40,
    size = 0.50
  ) +

  geom_point(
    data = results[
      importance_tier ==
        "Strong_consistent"
    ],
    aes(
      colour = DMR_direction
    ),
    alpha = 0.65,
    size = 0.70
  ) +

  geom_vline(
    xintercept = c(
      -strong_delta_beta_threshold,
      strong_delta_beta_threshold
    ),
    linetype = 2
  ) +

  geom_hline(
    yintercept =
      min(
        -log10(
          fdr_threshold
        ),
        maximum_minus_log10_FDR_display
      ),
    linetype = 2
  ) +

  ggrepel::geom_text_repel(
    data = volcano_labels,
    aes(
      label = region_id
    ),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.20,
    min.segment.length = 0
  ) +

  scale_colour_manual(
    values = c(
      "Tumor_higher" =
        "#2166AC",
      "Plasma_higher" =
        "#B2182B"
    )
  ) +

  labs(
    title =
      "SOLID tumor–plasma differential methylation",
    subtitle = paste0(
      "Strong consistent regions: FDR < 0.05, |median Δβ| ≥ ",
      strong_delta_beta_threshold,
      ", direction consistency ≥ ",
      minimum_direction_consistency,
      "%\nΔβ = Plasma - Tumor"
    ),
    x =
      "Median Δβ (Plasma - Tumor)",
    y =
      paste0(
        "-log10(FDR), capped at ",
        maximum_minus_log10_FDR_display
      ),
    colour =
      "Direction"
  ) +

  theme_bw(
    base_size = 12
  ) +

  theme(
    legend.position = "top",
    plot.title =
      element_text(
        face = "bold"
      )
  )


save_plot(
  volcano_plot,
  file.path(
    figure_dir,
    "SOLID_DMR_volcano_important_regions.png"
  ),
  width = 12,
  height = 8
)


ggsave(
  file.path(
    figure_dir,
    "SOLID_DMR_volcano_important_regions.pdf"
  ),
  volcano_plot,
  width = 12,
  height = 8
)


############################################################
## 15. DELTA-BETA DISTRIBUTION
############################################################

delta_distribution <- ggplot(
  results,
  aes(
    x = median_delta_beta
  )
) +

  geom_histogram(
    bins = 150,
    fill = "grey55",
    colour = "white",
    linewidth = 0.1
  ) +

  geom_vline(
    xintercept = c(
      -0.05,
      0.05
    ),
    linetype = 3
  ) +

  geom_vline(
    xintercept = c(
      -0.10,
      0.10
    ),
    linetype = 2
  ) +

  geom_vline(
    xintercept = c(
      -0.20,
      0.20
    ),
    linewidth = 0.7
  ) +

  labs(
    title =
      "Distribution of SOLID paired methylation effects",
    subtitle =
      "Reference lines at |Δβ| = 0.05, 0.10 and 0.20",
    x =
      "Median Δβ (Plasma - Tumor)",
    y =
      "Number of regions"
  ) +

  theme_bw(
    base_size = 12
  )


save_plot(
  delta_distribution,
  file.path(
    figure_dir,
    "SOLID_DMR_delta_beta_distribution.png"
  ),
  width = 10,
  height = 7
)


############################################################
## 16. STRONG-CANDIDATE DIRECTION SUMMARY
############################################################

direction_summary <- strong_regions[
  ,
  .(
    n_regions = .N
  ),
  by = DMR_direction
]


direction_summary[
  ,
  percentage :=
    100 *
    n_regions /
    sum(
      n_regions
    )
]


direction_plot <- ggplot(
  direction_summary,
  aes(
    x = DMR_direction,
    y = n_regions,
    fill = DMR_direction
  )
) +

  geom_col(
    width = 0.65
  ) +

  geom_text(
    aes(
      label = paste0(
        comma(
          n_regions
        ),
        "\n",
        sprintf(
          "%.1f%%",
          percentage
        )
      )
    ),
    vjust = -0.3
  ) +

  scale_fill_manual(
    values = c(
      "Tumor_higher" =
        "#2166AC",
      "Plasma_higher" =
        "#B2182B"
    )
  ) +

  scale_y_continuous(
    labels = comma,
    expand = expansion(
      mult = c(
        0,
        0.15
      )
    )
  ) +

  labs(
    title =
      "Direction of strong consistent SOLID DMRs",
    subtitle = paste0(
      "|median Δβ| ≥ ",
      strong_delta_beta_threshold,
      " and consistency ≥ ",
      minimum_direction_consistency,
      "%"
    ),
    x = NULL,
    y =
      "Number of regions"
  ) +

  theme_bw(
    base_size = 12
  ) +

  theme(
    legend.position = "none"
  )


save_plot(
  direction_plot,
  file.path(
    figure_dir,
    "SOLID_DMR_strong_candidate_direction.png"
  ),
  width = 8,
  height = 6
)


############################################################
## 17. EFFECT SIZE vs DIRECTION CONSISTENCY
##
## This is one of the most useful prioritization figures.
############################################################

effect_consistency_plot <- ggplot(
  results[
    FDR <
      fdr_threshold &
      abs_median_delta_beta >=
      primary_delta_beta_threshold
  ],
  aes(
    x = median_delta_beta,
    y = direction_consistency_pct,
    colour = DMR_direction
  )
) +

  geom_point(
    alpha = 0.25,
    size = 0.7
  ) +

  geom_hline(
    yintercept =
      minimum_direction_consistency,
    linetype = 2
  ) +

  geom_vline(
    xintercept = c(
      -strong_delta_beta_threshold,
      strong_delta_beta_threshold
    ),
    linetype = 2
  ) +

  ggrepel::geom_text_repel(
    data = volcano_labels,
    aes(
      label = region_id
    ),
    size = 3,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +

  scale_colour_manual(
    values = c(
      "Tumor_higher" =
        "#2166AC",
      "Plasma_higher" =
        "#B2182B"
    )
  ) +

  coord_cartesian(
    ylim = c(
      0,
      100
    )
  ) +

  labs(
    title =
      "SOLID DMR effect size and patient consistency",
    subtitle =
      "Upper outer quadrants identify large and recurrent effects",
    x =
      "Median Δβ (Plasma - Tumor)",
    y =
      "Direction consistency across matched patients (%)",
    colour =
      "Direction"
  ) +

  theme_bw(
    base_size = 12
  ) +

  theme(
    legend.position = "top"
  )


save_plot(
  effect_consistency_plot,
  file.path(
    figure_dir,
    "SOLID_DMR_effect_vs_direction_consistency.png"
  ),
  width = 11,
  height = 8
)


############################################################
## 18. SELECT REPRESENTATIVE IMPORTANT REGIONS
############################################################

message_header(
  "SELECTING REPRESENTATIVE IMPORTANT REGIONS"
)


top_tumor_regions <- head(
  important_tumor,
  n_heatmap_regions_per_direction
)


top_plasma_regions <- head(
  important_plasma,
  n_heatmap_regions_per_direction
)


top_heatmap_regions <- rbindlist(
  list(
    top_tumor_regions,
    top_plasma_regions
  ),
  use.names = TRUE,
  fill = TRUE
)


top_heatmap_regions <- top_heatmap_regions[
  !duplicated(
    region_id
  )
]


if (nrow(top_heatmap_regions) < 2L) {

  stop(
    "Too few important regions for heatmap generation."
  )
}


fwrite(
  top_heatmap_regions,
  file.path(
    priority_dir,
    "SOLID_representative_important_DMRs.tsv"
  ),
  sep = "\t"
)


############################################################
## 19. DELTA-BETA HEATMAP
############################################################

message_header(
  "CREATING DELTA-BETA HEATMAP"
)


heatmap_idx <- match(
  top_heatmap_regions$region_id,
  rownames(
    delta_beta
  )
)


stopifnot(
  !anyNA(
    heatmap_idx
  )
)


delta_heatmap <- delta_beta[
  heatmap_idx,
  ,
  drop = FALSE
]


delta_heatmap_display <- pmax(
  pmin(
    delta_heatmap,
    heatmap_clip_delta_beta
  ),
  -heatmap_clip_delta_beta
)


rownames(
  delta_heatmap_display
) <- top_heatmap_regions$region_id


############################################################
## Patient annotation
############################################################

patient_annotation <- as.data.frame(
  patient_metadata[
    ,
    intersect(
      c(
        "Grade",
        "Sex",
        "ECOG",
        "Response_RANO"
      ),
      names(
        patient_metadata
      )
    ),
    with = FALSE
  ]
)


rownames(
  patient_annotation
) <- patient_metadata$patient_id


patient_annotation <- patient_annotation[
  patient_ids,
  ,
  drop = FALSE
]


heatmap_colours <- circlize::colorRamp2(
  c(
    -heatmap_clip_delta_beta,
    0,
    heatmap_clip_delta_beta
  ),
  c(
    "#2166AC",
    "white",
    "#B2182B"
  )
)


row_split <- factor(
  top_heatmap_regions$DMR_direction,
  levels = c(
    "Tumor_higher",
    "Plasma_higher"
  )
)


delta_heatmap_object <- ComplexHeatmap::Heatmap(

  delta_heatmap_display,

  name =
    "Δβ",

  col =
    heatmap_colours,

  row_split =
    row_split,

  top_annotation =
    ComplexHeatmap::HeatmapAnnotation(
      df =
        patient_annotation
    ),

  cluster_rows =
    TRUE,

  cluster_columns =
    TRUE,

  show_row_names =
    TRUE,

  show_column_names =
    TRUE,

  row_names_gp =
    grid::gpar(
      fontsize = 6
    ),

  column_names_rot =
    45,

  column_title =
    "SOLID within-patient Δβ",

  na_col =
    "grey90",

  heatmap_legend_param =
    list(
      title =
        "Plasma - Tumor"
    )
)


png(
  file.path(
    figure_dir,
    "SOLID_DMR_important_regions_delta_beta_heatmap.png"
  ),
  width = 3000,
  height = 2600,
  res = 300
)


ComplexHeatmap::draw(
  delta_heatmap_object,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)


dev.off()


pdf(
  file.path(
    figure_dir,
    "SOLID_DMR_important_regions_delta_beta_heatmap.pdf"
  ),
  width = 11,
  height = 10
)


ComplexHeatmap::draw(
  delta_heatmap_object,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)


dev.off()


############################################################
## 20. ACTUAL TUMOR + PLASMA BETA HEATMAP
##
## This complements the Delta-Beta heatmap and shows the
## actual methylation state in each sample.
############################################################

message_header(
  "CREATING TUMOR + PLASMA BETA HEATMAP"
)


tissue_heatmap <- tissue_beta[
  heatmap_idx,
  ,
  drop = FALSE
]


plasma_heatmap <- plasma_beta[
  heatmap_idx,
  ,
  drop = FALSE
]


beta_heatmap <- cbind(
  tissue_heatmap,
  plasma_heatmap
)


colnames(
  beta_heatmap
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


rownames(
  beta_heatmap
) <- top_heatmap_regions$region_id


############################################################
## 26-column annotation
############################################################

sample_annotation <- data.frame(

  Sample_Type = c(
    rep(
      "Tumor",
      length(
        patient_ids
      )
    ),
    rep(
      "Plasma",
      length(
        patient_ids
      )
    )
  ),

  Grade = rep(
    patient_metadata$Grade,
    2
  ),

  Sex = rep(
    patient_metadata$Sex,
    2
  ),

  ECOG = rep(
    patient_metadata$ECOG,
    2
  ),

  stringsAsFactors = FALSE
)


rownames(
  sample_annotation
) <- colnames(
  beta_heatmap
)


beta_colours <- circlize::colorRamp2(
  c(
    0,
    0.5,
    1
  ),
  c(
    "#2166AC",
    "white",
    "#B2182B"
  )
)


beta_heatmap_object <- ComplexHeatmap::Heatmap(

  beta_heatmap,

  name =
    "Beta",

  col =
    beta_colours,

  top_annotation =
    ComplexHeatmap::HeatmapAnnotation(
      df =
        sample_annotation
    ),

  column_split =
    factor(
      sample_annotation$Sample_Type,
      levels = c(
        "Tumor",
        "Plasma"
      )
    ),

  cluster_rows =
    TRUE,

  cluster_columns =
    TRUE,

  show_row_names =
    TRUE,

  show_column_names =
    TRUE,

  row_names_gp =
    grid::gpar(
      fontsize = 6
    ),

  column_names_gp =
    grid::gpar(
      fontsize = 7
    ),

  column_names_rot =
    45,

  column_title =
    "SOLID important DMRs: actual regional Beta values"
)


png(
  file.path(
    figure_dir,
    "SOLID_DMR_important_regions_tumor_plasma_beta_heatmap.png"
  ),
  width = 3600,
  height = 2800,
  res = 300
)


ComplexHeatmap::draw(
  beta_heatmap_object,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)


dev.off()


pdf(
  file.path(
    figure_dir,
    "SOLID_DMR_important_regions_tumor_plasma_beta_heatmap.pdf"
  ),
  width = 14,
  height = 10
)


ComplexHeatmap::draw(
  beta_heatmap_object,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)


dev.off()


############################################################
## 21. TOP-REGION PAIRED BETA PLOTS
############################################################

message_header(
  "CREATING TOP-REGION PAIRED BETA PLOTS"
)


top_paired_regions <- rbindlist(
  list(

    head(
      important_tumor,
      n_paired_plot_regions_per_direction
    ),

    head(
      important_plasma,
      n_paired_plot_regions_per_direction
    )
  ),
  use.names = TRUE,
  fill = TRUE
)


top_paired_regions <- top_paired_regions[
  !duplicated(
    region_id
  )
]


paired_plot_data <- rbindlist(
  lapply(
    seq_len(
      nrow(
        top_paired_regions
      )
    ),
    function(i) {

      region_id_i <-
        top_paired_regions$region_id[i]


      matrix_row <- match(
        region_id_i,
        rownames(
          tissue_beta
        )
      )


      data.table(

        region_id =
          rep(
            region_id_i,
            2 *
              length(
                patient_ids
              )
          ),

        patient_id =
          rep(
            patient_ids,
            2
          ),

        sample_type =
          rep(
            c(
              "Tumor",
              "Plasma"
            ),
            each =
              length(
                patient_ids
              )
          ),

        beta =
          c(
            as.numeric(
              tissue_beta[
                matrix_row,
              ]
            ),
            as.numeric(
              plasma_beta[
                matrix_row,
              ]
            )
          )
      )
    }
  )
)


paired_plot_data[
  ,
  sample_type :=
    factor(
      sample_type,
      levels = c(
        "Tumor",
        "Plasma"
      )
    )
]


paired_plot <- ggplot(
  paired_plot_data,
  aes(
    x = sample_type,
    y = beta,
    group = patient_id
  )
) +

  geom_line(
    alpha = 0.45,
    colour = "grey50"
  ) +

  geom_point(
    size = 1.8
  ) +

  facet_wrap(
    ~region_id,
    ncol = 2
  ) +

  coord_cartesian(
    ylim = c(
      0,
      1
    )
  ) +

  labs(
    title =
      "Representative high-confidence SOLID DMRs",
    subtitle =
      "Each line connects matched tumor and plasma from one patient",
    x = NULL,
    y =
      "Regional Beta value"
  ) +

  theme_bw(
    base_size = 10
  ) +

  theme(
    strip.text =
      element_text(
        face = "bold",
        size = 8
      )
  )


save_plot(
  paired_plot,
  file.path(
    figure_dir,
    "SOLID_DMR_important_regions_paired_beta.png"
  ),
  width = 12,
  height = 13
)


fwrite(
  paired_plot_data,
  file.path(
    priority_dir,
    "SOLID_DMR_important_regions_paired_beta_values.tsv"
  ),
  sep = "\t"
)


############################################################
## 22. GENOME-WIDE EFFECT-SIZE PLOT
##
## More informative here than significance alone because
## a very large fraction of regions are FDR significant.
############################################################

message_header(
  "CREATING GENOME-WIDE EFFECT-SIZE PLOT"
)


required_coordinate_columns <- c(
  "chr",
  "start",
  "end"
)


if (
  all(
    required_coordinate_columns %in%
      names(
        results
      )
  )
) {


  genome_data <- copy(
    results
  )


  genome_data[
    ,
    chr_plot :=
      standardize_chr(
        chr
      )
  ]


  genome_data <- genome_data[
    !is.na(
      chr_plot
    ) &
      !is.na(
        start
      )
  ]


  genome_data[
    ,
    chr_numeric :=
      as.integer(
        chr_plot
      )
  ]


  chromosome_sizes <- genome_data[
    ,
    .(
      chromosome_end =
        max(
          end,
          na.rm = TRUE
        )
    ),
    by = .(
      chr_plot,
      chr_numeric
    )
  ]


  setorder(
    chromosome_sizes,
    chr_numeric
  )


  chromosome_sizes[
    ,
    cumulative_offset :=
      c(
        0,
        head(
          cumsum(
            chromosome_end
          ),
          -1L
        )
      )
  ]


  genome_data <- merge(
    genome_data,
    chromosome_sizes[
      ,
      .(
        chr_plot,
        cumulative_offset
      )
    ],
    by = "chr_plot",
    all.x = TRUE,
    sort = FALSE
  )


  genome_data[
    ,
    genomic_position :=
      start +
      cumulative_offset
  ]


  chromosome_centres <- genome_data[
    ,
    .(
      centre =
        (
          min(
            genomic_position
          ) +
            max(
              genomic_position
            )
        ) /
        2
    ),
    by = chr_plot
  ]


  chromosome_centres[
    ,
    chr_numeric :=
      as.integer(
        chr_plot
      )
  ]


  setorder(
    chromosome_centres,
    chr_numeric
  )


  genome_effect_plot <- ggplot(
    genome_data,
    aes(
      x = genomic_position,
      y = median_delta_beta
    )
  ) +

    geom_point(
      data = genome_data[
        importance_tier !=
          "Strong_consistent"
      ],
      alpha = 0.15,
      size = 0.30,
      colour = "grey65"
    ) +

    geom_point(
      data = genome_data[
        importance_tier ==
          "Strong_consistent"
      ],
      aes(
        colour = DMR_direction
      ),
      alpha = 0.70,
      size = 0.55
    ) +

    geom_hline(
      yintercept = 0,
      linetype = 2
    ) +

    geom_hline(
      yintercept = c(
        -strong_delta_beta_threshold,
        strong_delta_beta_threshold
      ),
      linetype = 3
    ) +

    scale_colour_manual(
      values = c(
        "Tumor_higher" =
          "#2166AC",
        "Plasma_higher" =
          "#B2182B"
      )
    ) +

    scale_x_continuous(
      breaks =
        chromosome_centres$centre,
      labels =
        as.character(
          chromosome_centres$chr_plot
        )
    ) +

    labs(
      title =
        "Genome-wide SOLID tumor–plasma methylation effects",
      subtitle =
        "Strong consistent DMRs highlighted",
      x =
        "Chromosome",
      y =
        "Median Δβ (Plasma - Tumor)",
      colour =
        "Direction"
    ) +

    theme_bw(
      base_size = 11
    ) +

    theme(
      panel.grid.major.x =
        element_blank(),
      legend.position =
        "top"
    )


  save_plot(
    genome_effect_plot,
    file.path(
      figure_dir,
      "SOLID_DMR_genomewide_effect_size.png"
    ),
    width = 16,
    height = 7
  )


  ############################################################
  ## 23. MERGE ADJACENT STRONG CONSISTENT WINDOWS
  ##
  ## Only same-direction high-confidence windows are merged.
  ############################################################

  message_header(
    "MERGING ADJACENT STRONG CONSISTENT WINDOWS"
  )


  block_input <- strong_regions[
    !is.na(
      chr
    ) &
      !is.na(
        start
      ) &
      !is.na(
        end
      ),
    .(
      region_id,
      chr =
        as.character(
          chr
        ),
      start =
        as.integer(
          start
        ),
      end =
        as.integer(
          end
        ),
      FDR,
      median_delta_beta,
      mean_delta_beta,
      direction_consistency_pct,
      effect_consistency_0.10_pct,
      DMR_direction
    )
  ]


  setorder(
    block_input,
    chr,
    DMR_direction,
    start,
    end
  )


  block_input[
    ,
    previous_end :=
      shift(
        end
      ),
    by = .(
      chr,
      DMR_direction
    )
  ]


  block_input[
    ,
    new_block :=
      is.na(
        previous_end
      ) |
        (
          start -
            previous_end
        ) >
        maximum_gap_between_windows,
    by = .(
      chr,
      DMR_direction
    )
  ]


  block_input[
    ,
    block_number :=
      cumsum(
        new_block
      ),
    by = .(
      chr,
      DMR_direction
    )
  ]


  dmr_blocks <- block_input[
    ,
    .(

      block_start =
        min(
          start
        ),

      block_end =
        max(
          end
        ),

      ########################################################
      ## Coordinates are half-open:
      ## width = end - start
      ########################################################

      block_width =
        max(
          end
        ) -
        min(
          start
        ),

      n_windows =
        .N,

      minimum_FDR =
        min(
          FDR,
          na.rm = TRUE
        ),

      median_region_delta_beta =
        median(
          median_delta_beta,
          na.rm = TRUE
        ),

      maximum_absolute_region_delta_beta =
        max(
          abs(
            median_delta_beta
          ),
          na.rm = TRUE
        ),

      median_direction_consistency_pct =
        median(
          direction_consistency_pct,
          na.rm = TRUE
        ),

      median_effect_consistency_0.10_pct =
        median(
          effect_consistency_0.10_pct,
          na.rm = TRUE
        ),

      first_region_id =
        first(
          region_id
        ),

      last_region_id =
        last(
          region_id
        )

    ),
    by = .(
      chr,
      DMR_direction,
      block_number
    )
  ]


  dmr_blocks[
    ,
    retained_block :=
      n_windows >=
      minimum_windows_per_block
  ]


  retained_blocks <- dmr_blocks[
    retained_block ==
      TRUE
  ]


  if (nrow(retained_blocks) > 0L) {

    retained_blocks[
      ,
      block_id :=
        paste0(
          "SOLID_DMR_BLOCK_",
          sprintf(
            "%05d",
            seq_len(
              .N
            )
          )
        )
    ]


    setcolorder(
      retained_blocks,
      c(
        "block_id",
        setdiff(
          names(
            retained_blocks
          ),
          "block_id"
        )
      )
    )


    setorder(
      retained_blocks,
      -n_windows,
      -maximum_absolute_region_delta_beta,
      minimum_FDR
    )
  }


  fwrite(
    dmr_blocks,
    file.path(
      priority_dir,
      "SOLID_all_strong_consistent_DMR_blocks.tsv.gz"
    ),
    sep = "\t",
    compress = "gzip"
  )


  fwrite(
    retained_blocks,
    file.path(
      priority_dir,
      "SOLID_retained_multiwindow_DMR_blocks.tsv"
    ),
    sep = "\t"
  )


  cat(
    "\nStrong candidate regions:",
    format(
      nrow(
        strong_regions
      ),
      big.mark = ","
    ),
    "\n"
  )


  cat(
    "Merged blocks:",
    format(
      nrow(
        dmr_blocks
      ),
      big.mark = ","
    ),
    "\n"
  )


  cat(
    "Multi-window retained blocks:",
    format(
      nrow(
        retained_blocks
      ),
      big.mark = ","
    ),
    "\n"
  )

} else {

  warning(
    "chr/start/end columns missing: genome-wide plot and block merging skipped."
  )

  retained_blocks <- data.table()
}


############################################################
## 24. FINAL IMPORTANT-REGION SUMMARY
############################################################

message_header(
  "IMPORTANT REGION SUMMARY"
)


important_summary <- data.table(

  item = c(

    "regions_tested",

    "FDR_significant",

    "primary_DMR_FDR_and_absDeltaBeta_0.05",

    "moderate_consistent_absDeltaBeta_0.10_consistency_80",

    "strong_consistent_absDeltaBeta_0.20_consistency_80",

    "strong_consistent_tumor_higher",

    "strong_consistent_plasma_higher",

    "retained_multiwindow_blocks"

  ),

  value = c(

    nrow(
      results
    ),

    sum(
      results$FDR <
        fdr_threshold,
      na.rm = TRUE
    ),

    sum(
      results$FDR <
        fdr_threshold &
        results$abs_median_delta_beta >=
        primary_delta_beta_threshold,
      na.rm = TRUE
    ),

    sum(
      results$FDR <
        fdr_threshold &
        results$abs_median_delta_beta >=
        moderate_delta_beta_threshold &
        results$direction_consistency_pct >=
        minimum_direction_consistency,
      na.rm = TRUE
    ),

    nrow(
      strong_regions
    ),

    nrow(
      important_tumor
    ),

    nrow(
      important_plasma
    ),

    nrow(
      retained_blocks
    )
  )
)


print(
  important_summary
)


fwrite(
  important_summary,
  file.path(
    priority_dir,
    "SOLID_DMR_important_region_summary.tsv"
  ),
  sep = "\t"
)


############################################################
## 25. FINAL VALIDATION
############################################################

required_figures <- c(

  "SOLID_DMR_volcano_important_regions.png",

  "SOLID_DMR_delta_beta_distribution.png",

  "SOLID_DMR_strong_candidate_direction.png",

  "SOLID_DMR_effect_vs_direction_consistency.png",

  "SOLID_DMR_important_regions_delta_beta_heatmap.png",

  "SOLID_DMR_important_regions_tumor_plasma_beta_heatmap.png",

  "SOLID_DMR_important_regions_paired_beta.png"
)


missing_figures <- required_figures[
  !file.exists(
    file.path(
      figure_dir,
      required_figures
    )
  )
]


if (length(missing_figures) > 0L) {

  stop(
    "Missing expected figure(s):\n",
    paste(
      missing_figures,
      collapse = "\n"
    )
  )
}


############################################################
## 26. FINAL CONSOLE SUMMARY
############################################################

cat(
  "\n============================================\n"
)

cat(
  "SOLID DMR VISUALIZATION COMPLETE\n"
)

cat(
  "============================================\n"
)


cat(
  "\nRegions tested:",
  format(
    nrow(
      results
    ),
    big.mark = ","
  ),
  "\n"
)


cat(
  "FDR-significant:",
  format(
    sum(
      results$FDR <
        fdr_threshold,
      na.rm = TRUE
    ),
    big.mark = ","
  ),
  "\n"
)


cat(
  "\nPrimary DMRs |Δβ| >= 0.05:",
  format(
    sum(
      results$FDR <
        fdr_threshold &
        results$abs_median_delta_beta >=
        0.05,
      na.rm = TRUE
    ),
    big.mark = ","
  ),
  "\n"
)


cat(
  "Moderate consistent DMRs:",
  format(
    sum(
      results$FDR <
        fdr_threshold &
        results$abs_median_delta_beta >=
        0.10 &
        results$direction_consistency_pct >=
        minimum_direction_consistency,
      na.rm = TRUE
    ),
    big.mark = ","
  ),
  "\n"
)


cat(
  "Strong consistent DMRs:",
  format(
    nrow(
      strong_regions
    ),
    big.mark = ","
  ),
  "\n"
)


cat(
  "  Tumor higher:",
  format(
    nrow(
      important_tumor
    ),
    big.mark = ","
  ),
  "\n"
)


cat(
  "  Plasma higher:",
  format(
    nrow(
      important_plasma
    ),
    big.mark = ","
  ),
  "\n"
)


cat(
  "\nDirection consistency among strong DMRs:\n"
)


if (nrow(strong_regions) > 0L) {

  print(
    summary(
      strong_regions$direction_consistency_pct
    )
  )
}


cat(
  "\n0.10-effect consistency among strong DMRs:\n"
)


if (nrow(strong_regions) > 0L) {

  print(
    summary(
      strong_regions$effect_consistency_0.10_pct
    )
  )
}


cat(
  "\nFigures:\n",
  figure_dir,
  "\n",
  sep = ""
)


cat(
  "\nPrioritization tables:\n",
  priority_dir,
  "\n",
  sep = ""
)


############################################################
## 27. SESSION INFORMATION
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    dmr_dir,
    "sessionInfo_08_SOLID_DMR_visualization.txt"
  )
)