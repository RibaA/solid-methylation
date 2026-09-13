############################################################
## 09_SOLID_DMR_visualization.r
##
## SOLID FINAL DMR VISUALIZATION
##
## Input:
##   result/03_matched_tissue_plasma/DMR/final/
##     SOLID_DMR_final_object.rds
##
## Purpose:
##   1. Visualize the frozen Script-08 final SOLID results
##   2. Use Tier-A as the high-confidence SOLID DMR set
##   3. Generate figures for presentation and integration
##
## IMPORTANT:
##   - No statistical testing
##   - No new candidate definition
##   - No block reconstruction
##   - No changes to Tier-A membership
############################################################

rm(list = ls())

options(
  stringsAsFactors = FALSE,
  width = 180,
  scipen = 999,
  warn = 1
)

fdr_threshold <- 0.05
tierA_abs_delta_beta <- 0.30
maximum_minus_log10_FDR_display <- 50
heatmap_clip_delta_beta <- 0.50
n_volcano_labels_per_direction <- 8L
n_heatmap_regions_per_direction <- 20L
n_paired_plot_regions_per_direction <- 4L

cran_packages <- c(
  "data.table",
  "ggplot2",
  "ggrepel",
  "scales",
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
    paste(missing_cran, collapse = ", ")
  )
}

if (length(missing_bioc) > 0L) {
  stop(
    "Missing Bioconductor package(s): ",
    paste(missing_bioc, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(circlize)
  library(ComplexHeatmap)
})

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

required_plot_objects <- c(
  "sample_cols",
  "grade_cols",
  "neutral_cols",
  "theme_project"
)

missing_plot_objects <- required_plot_objects[
  !vapply(
    required_plot_objects,
    exists,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_plot_objects) > 0L) {
  stop(
    "Missing shared plotting object(s): ",
    paste(missing_plot_objects, collapse = ", ")
  )
}

dmr_dir <- file.path(
  "result",
  "03_matched_tissue_plasma",
  "DMR"
)

final_dir <- file.path(
  dmr_dir,
  "final"
)

figure_dir <- file.path(
  final_dir,
  "figures"
)

table_dir <- file.path(
  final_dir,
  "visualization_tables"
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

final_object_file <- file.path(
  final_dir,
  "SOLID_DMR_final_object.rds"
)

if (!file.exists(final_object_file)) {
  stop(
    "Final SOLID DMR object not found:\n",
    final_object_file,
    "\nRun Script 08 first."
  )
}

message_header <- function(text) {
  cat(
    "\n",
    paste(rep("=", 72), collapse = ""),
    "\n",
    text,
    "\n",
    paste(rep("=", 72), collapse = ""),
    "\n",
    sep = ""
  )
}

save_plot_both <- function(
    plot_object,
    basename,
    width,
    height
) {

  ggsave(
    filename =
      paste0(
        basename,
        ".png"
      ),
    plot =
      plot_object,
    width =
      width,
    height =
      height,
    dpi =
      300,
    limitsize =
      FALSE
  )

  ggsave(
    filename =
      paste0(
        basename,
        ".pdf"
      ),
    plot =
      plot_object,
    width =
      width,
    height =
      height,
    limitsize =
      FALSE
  )
}

make_annotation_palette <- function(
    values,
    available_colors
) {

  values <- sort(
    unique(
      as.character(
        values[
          !is.na(values)
        ]
      )
    )
  )

  if (length(values) == 0L) {
    return(NULL)
  }

  if (
    length(values) >
      length(available_colors)
  ) {

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

standardize_chr_character <- function(x) {
  sub(
    "^chr",
    "",
    as.character(x),
    ignore.case = TRUE
  )
}

message_header(
  "LOADING FINAL SOLID DMR OBJECT"
)

obj <- readRDS(
  final_object_file
)

required_elements <- c(
  "all_results",
  "TierA_regions",
  "TierA_tumor_higher",
  "TierA_plasma_higher",
  "TierA_retained_blocks",
  "TierA_block_members",
  "region_order",
  "delta_beta",
  "tissue_beta",
  "plasma_beta",
  "valid_pair_mask",
  "patient_metadata",
  "final_summary",
  "final_settings"
)

missing_elements <- setdiff(
  required_elements,
  names(obj)
)

if (length(missing_elements) > 0L) {
  stop(
    "Final object missing required element(s): ",
    paste(missing_elements, collapse = ", ")
  )
}

results <- as.data.table(
  obj$all_results
)

tierA <- as.data.table(
  obj$TierA_regions
)

tierA_tumor <- as.data.table(
  obj$TierA_tumor_higher
)

tierA_plasma <- as.data.table(
  obj$TierA_plasma_higher
)

blocks <- as.data.table(
  obj$TierA_retained_blocks
)

block_members <- as.data.table(
  obj$TierA_block_members
)

delta_beta <- obj$delta_beta
tissue_beta <- obj$tissue_beta
plasma_beta <- obj$plasma_beta

patient_metadata <- as.data.table(
  obj$patient_metadata
)

patient_ids <- colnames(
  delta_beta
)

stopifnot(
  identical(
    obj$region_order,
    results$region_id
  ),

  identical(
    obj$region_order,
    rownames(delta_beta)
  ),

  identical(
    patient_ids,
    patient_metadata$patient_id
  ),

  nrow(tierA) == 15322L,

  nrow(blocks) == 2920L,

  nrow(block_members) == 4492L
)

cat(
  "Final object validation: PASS\n"
)

direction_cols <- c(
  "Tumor_higher" =
    unname(
      sample_cols["Tumor"]
    ),

  "Plasma_higher" =
    unname(
      sample_cols["Plasma"]
    )
)

sex_cols <- make_annotation_palette(
  patient_metadata$Sex,
  c(
    "#99B6BD",
    "#ECC9A0",
    "#BEB59C"
  )
)

ecog_cols <- make_annotation_palette(
  patient_metadata$ECOG,
  c(
    "#D2C396",
    "#A5A596",
    "#697878",
    "#4B5A69",
    "#5A4B3C"
  )
)

response_cols <- make_annotation_palette(
  patient_metadata$Response_RANO,
  c(
    "#99B6BD",
    "#B3A86A",
    "#ECC9A0",
    "#D0937D",
    "#78847F",
    "#9A9391"
  )
)

message_header(
  "CREATING FINAL TIER-A VOLCANO PLOT"
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

volcano_labels <- rbindlist(
  list(
    head(
      tierA_tumor,
      n_volcano_labels_per_direction
    ),
    head(
      tierA_plasma,
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
    alpha = 0.12,
    size = 0.35,
    color =
      unname(
        neutral_cols["medium"]
      )
  ) +

  geom_point(
    data = results[
      TierA == TRUE
    ],
    aes(
      color = DMR_direction
    ),
    alpha = 0.70,
    size = 0.75
  ) +

  geom_vline(
    xintercept = c(
      -tierA_abs_delta_beta,
      tierA_abs_delta_beta
    ),
    linetype = 2,
    linewidth = 0.5
  ) +

  geom_hline(
    yintercept =
      -log10(
        fdr_threshold
      ),
    linetype = 2,
    linewidth = 0.5
  ) +

  ggrepel::geom_text_repel(
    data =
      volcano_labels,
    aes(
      label =
        region_id,
      color =
        DMR_direction
    ),
    size =
      3,
    max.overlaps =
      Inf,
    box.padding =
      0.35,
    point.padding =
      0.20,
    min.segment.length =
      0,
    show.legend =
      FALSE
  ) +

  scale_color_manual(
    values =
      direction_cols
  ) +

  labs(
    title =
      "SOLID tumor-plasma differential methylation",

    subtitle =
      "Tier-A highlighted; Delta Beta = Plasma - Tumor",

    x =
      "Median Delta Beta (Plasma - Tumor)",

    y =
      paste0(
        "-log10(FDR), capped at ",
        maximum_minus_log10_FDR_display
      ),

    color =
      "Direction"
  ) +

  theme_project() +

  theme(
    legend.position =
      "top"
  )

save_plot_both(
  volcano_plot,
  file.path(
    figure_dir,
    "SOLID_TierA_volcano"
  ),
  width = 12,
  height = 8
)

tierA_effect_plot <- ggplot(
  tierA,
  aes(
    x =
      median_delta_beta,
    fill =
      DMR_direction
  )
) +

  geom_histogram(
    bins =
      80,
    alpha =
      0.72,
    position =
      "identity"
  ) +

  geom_vline(
    xintercept =
      c(
        -tierA_abs_delta_beta,
        tierA_abs_delta_beta
      ),
    linetype =
      2,
    linewidth =
      0.5
  ) +

  scale_fill_manual(
    values =
      direction_cols
  ) +

  labs(
    title =
      "SOLID Tier-A methylation effects",

    subtitle =
      "FDR < 0.05; |median Delta Beta| >= 0.30; 100% direction consistency",

    x =
      "Median Delta Beta (Plasma - Tumor)",

    y =
      "Number of Tier-A regions",

    fill =
      "Direction"
  ) +

  theme_project() +

  theme(
    legend.position =
      "top"
  )

save_plot_both(
  tierA_effect_plot,
  file.path(
    figure_dir,
    "SOLID_TierA_delta_beta_distribution"
  ),
  width = 9,
  height = 6
)

direction_summary <- tierA[
  ,
  .(
    n_regions =
      .N
  ),
  by =
    DMR_direction
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
    x =
      DMR_direction,
    y =
      n_regions,
    fill =
      DMR_direction
  )
) +

  geom_col(
    width =
      0.65
  ) +

  geom_text(
    aes(
      label =
        paste0(
          scales::comma(
            n_regions
          ),
          "\n",
          sprintf(
            "%.1f%%",
            percentage
          )
        )
    ),
    vjust =
      -0.3
  ) +

  scale_fill_manual(
    values =
      direction_cols
  ) +

  scale_y_continuous(
    labels =
      scales::comma,
    expand =
      expansion(
        mult =
          c(
            0,
            0.15
          )
      )
  ) +

  labs(
    title =
      "Direction of SOLID Tier-A DMRs",

    x =
      NULL,

    y =
      "Number of Tier-A regions"
  ) +

  theme_project() +

  theme(
    legend.position =
      "none"
  )

save_plot_both(
  direction_plot,
  file.path(
    figure_dir,
    "SOLID_TierA_direction"
  ),
  width = 8,
  height = 6
)

effect_consistency_plot <- ggplot(
  results[
    FDR <
      fdr_threshold
  ],
  aes(
    x =
      median_delta_beta,
    y =
      effect_consistency_0.10_pct
  )
) +

  geom_point(
    alpha =
      0.12,
    size =
      0.45,
    color =
      unname(
        neutral_cols["medium"]
      )
  ) +

  geom_point(
    data =
      results[
        TierA ==
          TRUE
      ],
    aes(
      color =
        DMR_direction
    ),
    alpha =
      0.75,
    size =
      0.80
  ) +

  geom_hline(
    yintercept =
      80,
    linetype =
      2,
    linewidth =
      0.5
  ) +

  geom_vline(
    xintercept =
      c(
        -tierA_abs_delta_beta,
        tierA_abs_delta_beta
      ),
    linetype =
      2,
    linewidth =
      0.5
  ) +

  scale_color_manual(
    values =
      direction_cols
  ) +

  coord_cartesian(
    ylim =
      c(
        0,
        100
      )
  ) +

  labs(
    title =
      "SOLID DMR effect size and recurrent effect consistency",

    subtitle =
      "Tier-A highlighted",

    x =
      "Median Delta Beta (Plasma - Tumor)",

    y =
      "Same-direction |Delta Beta| >= 0.10 (%)",

    color =
      "Direction"
  ) +

  theme_project() +

  theme(
    legend.position =
      "top"
  )

save_plot_both(
  effect_consistency_plot,
  file.path(
    figure_dir,
    "SOLID_TierA_effect_vs_recurrent_effect_consistency"
  ),
  width = 11,
  height = 8
)

top_tumor_regions <- head(
  tierA_tumor,
  n_heatmap_regions_per_direction
)

top_plasma_regions <- head(
  tierA_plasma,
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

fwrite(
  top_heatmap_regions,
  file.path(
    table_dir,
    "SOLID_TierA_representative_heatmap_regions.tsv"
  ),
  sep = "\t"
)

message_header(
  "CREATING TIER-A DELTA-BETA HEATMAP"
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

patient_annotation <- as.data.frame(
  patient_metadata[
    ,
    c(
      "Grade",
      "Sex",
      "ECOG",
      "Response_RANO"
    ),
    with = FALSE
  ]
)

rownames(
  patient_annotation
) <- patient_ids

for (v in names(patient_annotation)) {
  patient_annotation[[v]] <- factor(
    patient_annotation[[v]]
  )
}

annotation_colors <- list(
  Grade =
    grade_cols[
      intersect(
        names(
          grade_cols
        ),
        levels(
          patient_annotation$Grade
        )
      )
    ],

  Sex =
    sex_cols,

  ECOG =
    ecog_cols,

  Response_RANO =
    response_cols
)

delta_heatmap_colours <- circlize::colorRamp2(
  c(
    -heatmap_clip_delta_beta,
    0,
    heatmap_clip_delta_beta
  ),
  c(
    "#526A83",
    "#F4F1EA",
    "#BF816B"
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
    "Delta Beta",

  col =
    delta_heatmap_colours,

  row_split =
    row_split,

  top_annotation =
    ComplexHeatmap::HeatmapAnnotation(
      df =
        patient_annotation,
      col =
        annotation_colors
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
    "SOLID Tier-A within-patient Delta Beta",

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
    "SOLID_TierA_delta_beta_heatmap.png"
  ),
  width = 3000,
  height = 2600,
  res = 300
)

ComplexHeatmap::draw(
  delta_heatmap_object,
  heatmap_legend_side =
    "right",
  annotation_legend_side =
    "right"
)

dev.off()

pdf(
  file.path(
    figure_dir,
    "SOLID_TierA_delta_beta_heatmap.pdf"
  ),
  width = 11,
  height = 10
)

ComplexHeatmap::draw(
  delta_heatmap_object,
  heatmap_legend_side =
    "right",
  annotation_legend_side =
    "right"
)

dev.off()

message_header(
  "CREATING TIER-A TUMOR + PLASMA BETA HEATMAP"
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

sample_annotation <- data.frame(
  Sample_Type =
    c(
      rep(
        "Tumor",
        length(patient_ids)
      ),
      rep(
        "Plasma",
        length(patient_ids)
      )
    ),

  Grade =
    rep(
      patient_metadata$Grade,
      2
    ),

  Sex =
    rep(
      patient_metadata$Sex,
      2
    ),

  ECOG =
    rep(
      patient_metadata$ECOG,
      2
    ),

  Response_RANO =
    rep(
      patient_metadata$Response_RANO,
      2
    ),

  stringsAsFactors =
    FALSE
)

rownames(
  sample_annotation
) <- colnames(
  beta_heatmap
)

for (v in names(sample_annotation)) {
  sample_annotation[[v]] <- factor(
    sample_annotation[[v]]
  )
}

sample_annotation_colors <- list(
  Sample_Type = c(
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
    grade_cols[
      intersect(
        names(
          grade_cols
        ),
        levels(
          sample_annotation$Grade
        )
      )
    ],

  Sex =
    sex_cols,

  ECOG =
    ecog_cols,

  Response_RANO =
    response_cols
)

beta_colours <- circlize::colorRamp2(
  c(
    0,
    0.5,
    1
  ),
  c(
    "#F4F1EA",
    "#AADCE0",
    "#526A83"
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
        sample_annotation,
      col =
        sample_annotation_colors
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
    "SOLID Tier-A regional Beta values"
)

png(
  file.path(
    figure_dir,
    "SOLID_TierA_tumor_plasma_beta_heatmap.png"
  ),
  width = 3600,
  height = 2800,
  res = 300
)

ComplexHeatmap::draw(
  beta_heatmap_object,
  heatmap_legend_side =
    "right",
  annotation_legend_side =
    "right"
)

dev.off()

pdf(
  file.path(
    figure_dir,
    "SOLID_TierA_tumor_plasma_beta_heatmap.pdf"
  ),
  width = 14,
  height = 10
)

ComplexHeatmap::draw(
  beta_heatmap_object,
  heatmap_legend_side =
    "right",
  annotation_legend_side =
    "right"
)

dev.off()

top_paired_regions <- rbindlist(
  list(
    head(
      tierA_tumor,
      n_paired_plot_regions_per_direction
    ),
    head(
      tierA_plasma,
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

      region_i <-
        top_paired_regions$region_id[i]

      idx <- match(
        region_i,
        rownames(
          tissue_beta
        )
      )

      data.table(
        region_id =
          rep(
            region_i,
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
                idx,
              ]
            ),
            as.numeric(
              plasma_beta[
                idx,
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

paired_plot_data[
  ,
  region_id :=
    factor(
      region_id,
      levels =
        top_paired_regions$region_id
    )
]

paired_plot <- ggplot(
  paired_plot_data,
  aes(
    x =
      sample_type,
    y =
      beta,
    group =
      patient_id
  )
) +

  geom_line(
    alpha =
      0.45,
    color =
      "grey65",
    linewidth =
      0.5
  ) +

  geom_point(
    aes(
      color =
        sample_type
    ),
    size =
      1.9
  ) +

  facet_wrap(
    ~ region_id,
    ncol =
      2
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

  coord_cartesian(
    ylim =
      c(
        0,
        1
      )
  ) +

  labs(
    title =
      "Representative SOLID Tier-A DMRs",

    subtitle =
      "Each line connects matched tumor and plasma from one patient",

    x =
      NULL,

    y =
      "Regional Beta value",

    color =
      "Sample type"
  ) +

  theme_project(
    base_size = 10
  ) +

  theme(
    legend.position =
      "top"
  )

save_plot_both(
  paired_plot,
  file.path(
    figure_dir,
    "SOLID_TierA_representative_paired_beta"
  ),
  width = 12,
  height = 13
)

fwrite(
  paired_plot_data,
  file.path(
    table_dir,
    "SOLID_TierA_representative_paired_beta_values.tsv"
  ),
  sep = "\t"
)

standard_chr_levels <- c(
  as.character(
    1:22
  ),
  "X",
  "Y"
)

tierA[
  ,
  chr_clean :=
    standardize_chr_character(
      chr
    )
]

tierA_chr <- tierA[
  chr_clean %in%
    standard_chr_levels,
  .(
    n_regions =
      .N
  ),
  by =
    chr_clean
]

chr_template <- data.table(
  chr_clean =
    standard_chr_levels,
  chr_order =
    seq_along(
      standard_chr_levels
    )
)

tierA_chr <- merge(
  chr_template,
  tierA_chr,
  by =
    "chr_clean",
  all.x =
    TRUE,
  sort =
    FALSE
)

tierA_chr[
  is.na(
    n_regions
  ),
  n_regions :=
    0L
]

setorder(
  tierA_chr,
  chr_order
)

tierA_chr[
  ,
  chr_factor :=
    factor(
      chr_clean,
      levels =
        standard_chr_levels
    )
]

chr_plot <- ggplot(
  tierA_chr,
  aes(
    x =
      chr_factor,
    y =
      n_regions
  )
) +

  geom_col(
    fill =
      unname(
        neutral_cols["dark"]
      )
  ) +

  scale_y_continuous(
    labels =
      scales::comma
  ) +

  labs(
    title =
      "SOLID Tier-A regions by chromosome",

    x =
      "Chromosome",

    y =
      "Number of Tier-A regions"
  ) +

  theme_project()

save_plot_both(
  chr_plot,
  file.path(
    figure_dir,
    "SOLID_TierA_regions_by_chromosome"
  ),
  width = 12,
  height = 6
)

block_size_plot <- ggplot(
  blocks,
  aes(
    x =
      n_windows,
    fill =
      DMR_direction
  )
) +

  geom_histogram(
    binwidth =
      1,
    alpha =
      0.70,
    position =
      "identity",
    boundary =
      0.5
  ) +

  scale_fill_manual(
    values =
      direction_cols
  ) +

  scale_y_continuous(
    labels =
      scales::comma
  ) +

  labs(
    title =
      "SOLID Tier-A multi-window block sizes",

    x =
      "Number of adjacent 1-kb windows",

    y =
      "Number of retained blocks",

    fill =
      "Direction"
  ) +

  theme_project() +

  theme(
    legend.position =
      "top"
  )

save_plot_both(
  block_size_plot,
  file.path(
    figure_dir,
    "SOLID_TierA_block_size_distribution"
  ),
  width = 9,
  height = 6
)

hg38_chr_lengths <- data.table(
  chr =
    c(
      as.character(
        1:22
      ),
      "X",
      "Y"
    ),

  chr_length =
    c(
      248956422,
      242193529,
      198295559,
      190214555,
      181538259,
      170805979,
      159345973,
      145138636,
      138394717,
      133797422,
      135086622,
      133275309,
      114364328,
      107043718,
      101991189,
      90338345,
      83257441,
      80373285,
      58617616,
      64444167,
      46709983,
      50818468,
      156040895,
      57227415
    )
)

hg38_chr_lengths[
  ,
  chr_order :=
    seq_len(.N)
]

hg38_chr_lengths[
  ,
  cumulative_offset :=
    c(
      0,
      head(
        cumsum(
          chr_length
        ),
        -1L
      )
    )
]

genome_data <- copy(
  results
)

genome_data[
  ,
  chr_clean :=
    standardize_chr_character(
      chr
    )
]

genome_data <- merge(
  genome_data,
  hg38_chr_lengths[
    ,
    .(
      chr,
      chr_order,
      cumulative_offset
    )
  ],
  by.x =
    "chr_clean",
  by.y =
    "chr",
  all.x =
    TRUE,
  sort =
    FALSE
)

genome_data <- genome_data[
  !is.na(
    chr_order
  ) &
    !is.na(
      start
    )
]

genome_data[
  ,
  genomic_position :=
    as.numeric(
      start
    ) +
    cumulative_offset
]

chromosome_centres <- hg38_chr_lengths[
  ,
  .(
    chr,
    centre =
      cumulative_offset +
      chr_length / 2
  )
]

genome_effect_plot <- ggplot(
  genome_data,
  aes(
    x =
      genomic_position,
    y =
      median_delta_beta
  )
) +

  geom_point(
    alpha =
      0.10,
    size =
      0.28,
    color =
      unname(
        neutral_cols["medium"]
      )
  ) +

  geom_point(
    data =
      genome_data[
        TierA ==
          TRUE
      ],
    aes(
      color =
        DMR_direction
    ),
    alpha =
      0.75,
    size =
      0.58
  ) +

  geom_hline(
    yintercept =
      0,
    linetype =
      2,
    linewidth =
      0.4
  ) +

  geom_hline(
    yintercept =
      c(
        -tierA_abs_delta_beta,
        tierA_abs_delta_beta
      ),
    linetype =
      3,
    linewidth =
      0.4
  ) +

  scale_color_manual(
    values =
      direction_cols
  ) +

  scale_x_continuous(
    breaks =
      chromosome_centres$centre,
    labels =
      chromosome_centres$chr
  ) +

  labs(
    title =
      "Genome-wide SOLID tumor-plasma methylation effects",

    subtitle =
      "Tier-A DMRs highlighted",

    x =
      "Chromosome",

    y =
      "Median Delta Beta (Plasma - Tumor)",

    color =
      "Direction"
  ) +

  theme_project(
    base_size = 11
  ) +

  theme(
    panel.grid.major.x =
      element_blank(),
    legend.position =
      "top"
  )

save_plot_both(
  genome_effect_plot,
  file.path(
    figure_dir,
    "SOLID_TierA_genomewide_effect_size"
  ),
  width = 16,
  height = 7
)

fwrite(
  direction_summary,
  file.path(
    table_dir,
    "SOLID_TierA_direction_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  tierA_chr,
  file.path(
    table_dir,
    "SOLID_TierA_chromosome_summary.tsv"
  ),
  sep = "\t"
)

block_summary <- blocks[
  ,
  .(
    n_blocks =
      .N,

    median_windows =
      median(
        n_windows,
        na.rm = TRUE
      ),

    maximum_windows =
      max(
        n_windows,
        na.rm = TRUE
      ),

    median_block_width =
      median(
        block_width,
        na.rm = TRUE
      ),

    maximum_block_width =
      max(
        block_width,
        na.rm = TRUE
      )
  ),
  by =
    DMR_direction
]

fwrite(
  block_summary,
  file.path(
    table_dir,
    "SOLID_TierA_block_summary.tsv"
  ),
  sep = "\t"
)

required_figures <- c(
  "SOLID_TierA_volcano.png",
  "SOLID_TierA_delta_beta_distribution.png",
  "SOLID_TierA_direction.png",
  "SOLID_TierA_effect_vs_recurrent_effect_consistency.png",
  "SOLID_TierA_delta_beta_heatmap.png",
  "SOLID_TierA_tumor_plasma_beta_heatmap.png",
  "SOLID_TierA_representative_paired_beta.png",
  "SOLID_TierA_regions_by_chromosome.png",
  "SOLID_TierA_block_size_distribution.png",
  "SOLID_TierA_genomewide_effect_size.png"
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

message_header(
  "SOLID FINAL DMR VISUALIZATION COMPLETE"
)

cat(
  "Regions tested: ",
  format(nrow(results), big.mark = ","),
  "\n",
  sep = ""
)

cat(
  "Tier-A regions: ",
  format(nrow(tierA), big.mark = ","),
  "\n",
  sep = ""
)

cat(
  "  Tumor higher: ",
  format(nrow(tierA_tumor), big.mark = ","),
  "\n",
  sep = ""
)

cat(
  "  Plasma higher: ",
  format(nrow(tierA_plasma), big.mark = ","),
  "\n",
  sep = ""
)

cat(
  "Tier-A windows in retained blocks: ",
  format(nrow(block_members), big.mark = ","),
  "\n",
  sep = ""
)

cat(
  "Tier-A retained blocks: ",
  format(nrow(blocks), big.mark = ","),
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
  "\nVisualization tables:\n",
  table_dir,
  "\n",
  sep = ""
)

capture.output(
  sessionInfo(),
  file = file.path(
    final_dir,
    "sessionInfo_09_SOLID_DMR_visualization.txt"
  )
)
