##------------------------------------------------------------
## 06. SOLID DMR FIGURES, PRIORITIZATION, AND BLOCK MERGING
##
## Purpose:
##   1. Load paired SOLID tissue-versus-plasma DMR results
##   2. apply the primary DMR-candidate definition
##   3. generate publication-style summary figures
##   4. generate a volcano plot aligned with the OCTANE workflow
##   5. create delta-Beta distribution and direction summaries
##   6. create a heatmap of representative top DMRs
##   7. create paired tissue-versus-plasma plots for top regions
##   8. create a chromosome-level Manhattan-style plot
##   9. merge adjacent significant 1-kb windows into broader DMR blocks
##  10. save prioritized region and block tables
##
## Primary high-confidence DMR definition:
##   - FDR < 0.05
##   - absolute median delta Beta >= 0.20
##
## Sensitivity thresholds reported:
##   - 0.05, 0.10, 0.15, and 0.20
##
## Delta Beta:
##   plasma minus tissue
##
## Therefore:
##   positive delta Beta = plasma higher
##   negative delta Beta = tissue higher
##
## Run after:
##   05_paired_DMR_SOLID_tissue_vs_plasma.R
##
## Inputs:
##   result/tissue-plasma/paired_DMR/
##     SOLID_tissue_vs_plasma_DMR_object.rds
##
## Outputs:
##   result/tissue-plasma/DMR_figures/
##   result/tissue-plasma/DMR_annotation/
##------------------------------------------------------------

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  warn = 1
)

##------------------------------------------------------------
## 1. USER SETTINGS
##------------------------------------------------------------

project_dir <- "C:/solid-methylation"

fdr_threshold <- 0.05
minimum_absolute_median_delta_beta <- 0.20

## Representative-region display settings
n_heatmap_regions_per_direction <- 50L
n_paired_plot_regions_per_direction <- 5L
n_top_table_per_direction <- 1000L

## Volcano and Manhattan background down-sampling
maximum_background_points <- 300000L
random_seed <- 123L

## DMR block-merging settings
maximum_gap_between_windows <- 1000L
minimum_windows_per_block <- 2L

## Heatmap settings
heatmap_clip_delta_beta <- 0.50

set.seed(
  random_seed
)

##------------------------------------------------------------
## 2. REQUIRED PACKAGES
##------------------------------------------------------------

cran_packages <- c(
  "data.table",
  "ggplot2",
  "matrixStats",
  "scales",
  "patchwork",
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
    paste(missing_cran, collapse = ", "),
    "\nInstall with:\n",
    "install.packages(c(",
    paste(
      sprintf('"%s"', missing_cran),
      collapse = ", "
    ),
    "))"
  )
}

if (length(missing_bioc) > 0L) {
  stop(
    "Missing Bioconductor package(s): ",
    paste(missing_bioc, collapse = ", "),
    "\nInstall with:\n",
    "if (!requireNamespace(\"BiocManager\", quietly = TRUE)) ",
    "install.packages(\"BiocManager\")\n",
    "BiocManager::install(c(",
    paste(
      sprintf('"%s"', missing_bioc),
      collapse = ", "
    ),
    "))"
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(matrixStats)
  library(scales)
  library(patchwork)
  library(circlize)
  library(ComplexHeatmap)
})

##------------------------------------------------------------
## 3. DIRECTORIES AND INPUT FILE
##------------------------------------------------------------

dmr_dir <- file.path(
  project_dir,
  "result",
  "tissue-plasma",
  "paired_DMR"
)

figure_dir <- file.path(
  project_dir,
  "result",
  "tissue-plasma",
  "DMR_figures"
)

annotation_dir <- file.path(
  project_dir,
  "result",
  "tissue-plasma",
  "DMR_annotation"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  annotation_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dmr_object_file <- file.path(
  dmr_dir,
  "SOLID_tissue_vs_plasma_DMR_object.rds"
)

if (!file.exists(dmr_object_file)) {
  stop(
    "Input DMR object not found:\n",
    dmr_object_file,
    "\nRun script 05 first."
  )
}

##------------------------------------------------------------
## 4. HELPER FUNCTIONS
##------------------------------------------------------------

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

  invisible(filename)
}

sample_background_points <- function(
    dat,
    maximum_points
) {
  if (nrow(dat) <= maximum_points) {
    return(copy(dat))
  }

  dat[
    sample(
      .N,
      maximum_points,
      replace = FALSE
    )
  ]
}

safe_filename <- function(x) {
  x <- gsub(
    "[^A-Za-z0-9_\\-]+",
    "_",
    x
  )

  gsub(
    "_+",
    "_",
    x
  )
}

standardize_chr <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)

  factor(
    x,
    levels = c(
      as.character(1:22),
      "X",
      "Y",
      "M",
      "MT"
    ),
    ordered = TRUE
  )
}

##------------------------------------------------------------
## 5. LOAD AND VALIDATE DMR OBJECT
##------------------------------------------------------------

message_header(
  "LOADING SOLID DMR OBJECT"
)

dmr_object <- readRDS(
  dmr_object_file
)

required_elements <- c(
  "all_results",
  "candidate_DMRs",
  "delta_M",
  "delta_beta",
  "tissue_beta",
  "plasma_beta",
  "matched_metadata"
)

missing_elements <- setdiff(
  required_elements,
  names(dmr_object)
)

if (length(missing_elements) > 0L) {
  stop(
    "The DMR object is missing required element(s): ",
    paste(missing_elements, collapse = ", ")
  )
}

results <- as.data.table(
  dmr_object$all_results
)

delta_beta <- dmr_object$delta_beta
tissue_beta <- dmr_object$tissue_beta
plasma_beta <- dmr_object$plasma_beta

matched_metadata <- as.data.table(
  dmr_object$matched_metadata
)

stopifnot(
  nrow(results) > 0L,
  is.matrix(delta_beta),
  is.matrix(tissue_beta),
  is.matrix(plasma_beta),
  nrow(delta_beta) == nrow(results),
  ncol(delta_beta) == 13L,
  identical(
    rownames(delta_beta),
    rownames(tissue_beta)
  ),
  identical(
    rownames(delta_beta),
    rownames(plasma_beta)
  ),
  identical(
    colnames(delta_beta),
    colnames(tissue_beta)
  ),
  identical(
    colnames(delta_beta),
    colnames(plasma_beta)
  ),
  setequal(
    results$region_id,
    rownames(delta_beta)
  )
)

## Restore matrix order in results
results <- results[
  match(
    rownames(delta_beta),
    region_id
  )
]

stopifnot(
  identical(
    results$region_id,
    rownames(delta_beta)
  )
)

cat(
  "Regions loaded:",
  format(nrow(results), big.mark = ","),
  "\n"
)

cat(
  "Matched patients:",
  ncol(delta_beta),
  "\n"
)

##------------------------------------------------------------
## 6. DEFINE PRIMARY CANDIDATES AND DIRECTIONS
##------------------------------------------------------------

message_header(
  "DEFINING PRIMARY DMR CANDIDATES"
)

results[
  ,
  primary_candidate :=
    FDR < fdr_threshold &
      abs(median_delta_beta) >=
        minimum_absolute_median_delta_beta
]

results[
  ,
  candidate_direction :=
    fifelse(
      primary_candidate &
        median_delta_beta > 0,
      "plasma_higher",
      fifelse(
        primary_candidate &
          median_delta_beta < 0,
        "tissue_higher",
        "not_candidate"
      )
    )
]

candidate_results <- results[
  primary_candidate == TRUE
]


## Create an effect-size sensitivity summary without rerunning limma
effect_cutoffs <- c(
  0.05,
  0.10,
  0.15,
  0.20
)

effect_cutoff_summary <- rbindlist(
  lapply(
    effect_cutoffs,
    function(cutoff) {
      candidate_mask <- (
        results$FDR < fdr_threshold &
          abs(results$median_delta_beta) >= cutoff
      )

      data.table(
        FDR_threshold =
          fdr_threshold,
        absolute_median_delta_beta_cutoff =
          cutoff,
        candidate_regions =
          sum(candidate_mask, na.rm = TRUE),
        tissue_higher_regions =
          sum(
            candidate_mask &
              results$median_delta_beta < 0,
            na.rm = TRUE
          ),
        plasma_higher_regions =
          sum(
            candidate_mask &
              results$median_delta_beta > 0,
            na.rm = TRUE
          ),
        percentage_of_tested_regions =
          100 *
          sum(candidate_mask, na.rm = TRUE) /
          nrow(results)
      )
    }
  )
)

print(
  effect_cutoff_summary
)

fwrite(
  effect_cutoff_summary,
  file.path(
    annotation_dir,
    "SOLID_DMR_effect_size_sensitivity.tsv"
  ),
  sep = "\t"
)

plasma_higher_results <- candidate_results[
  candidate_direction ==
    "plasma_higher"
]

tissue_higher_results <- candidate_results[
  candidate_direction ==
    "tissue_higher"
]

cat(
  "Primary DMR candidates:",
  format(nrow(candidate_results), big.mark = ","),
  "\n"
)

cat(
  "Plasma-higher candidates:",
  format(nrow(plasma_higher_results), big.mark = ","),
  "\n"
)

cat(
  "Tissue-higher candidates:",
  format(nrow(tissue_higher_results), big.mark = ","),
  "\n"
)

if (nrow(candidate_results) == 0L) {
  stop(
    "No regions satisfy the primary DMR definition."
  )
}

##------------------------------------------------------------
## 7. SAVE PRIORITIZED REGION TABLES
##------------------------------------------------------------

message_header(
  "SAVING PRIORITIZED REGION TABLES"
)

setorder(
  candidate_results,
  -abs_median_delta_beta,
  FDR
)

setorder(
  plasma_higher_results,
  -abs_median_delta_beta,
  FDR
)

setorder(
  tissue_higher_results,
  -abs_median_delta_beta,
  FDR
)

fwrite(
  candidate_results,
  file.path(
    annotation_dir,
    "SOLID_primary_DMR_candidates.tsv.gz"
  ),
  sep = "\t",
  compress = "gzip"
)

fwrite(
  head(
    plasma_higher_results,
    n_top_table_per_direction
  ),
  file.path(
    annotation_dir,
    "SOLID_top_plasma_higher_DMRs.tsv"
  ),
  sep = "\t"
)

fwrite(
  head(
    tissue_higher_results,
    n_top_table_per_direction
  ),
  file.path(
    annotation_dir,
    "SOLID_top_tissue_higher_DMRs.tsv"
  ),
  sep = "\t"
)

##------------------------------------------------------------
## 8. VOLCANO PLOT
##------------------------------------------------------------

message_header(
  "CREATING VOLCANO PLOT"
)

background_results <- results[
  primary_candidate == FALSE
]

background_display <- sample_background_points(
  background_results,
  maximum_background_points
)

volcano_data <- rbindlist(
  list(
    background_display[
      ,
      display_class :=
        "not_candidate"
    ],
    copy(candidate_results)[
      ,
      display_class :=
        candidate_direction
    ]
  ),
  use.names = TRUE,
  fill = TRUE
)

volcano_data[
  ,
  minus_log10_FDR :=
    -log10(
      pmax(
        FDR,
        .Machine$double.xmin
      )
    )
]

volcano_plot <- ggplot(
  volcano_data,
  aes(
    x = median_delta_beta,
    y = minus_log10_FDR
  )
) +
  geom_point(
    data = volcano_data[
      display_class ==
        "not_candidate"
    ],
    alpha = 0.15,
    size = 0.35,
    colour = "grey60"
  ) +
  geom_point(
    data = volcano_data[
      display_class ==
        "tissue_higher"
    ],
    alpha = 0.35,
    size = 0.45,
    colour = "#2b8cbe"
  ) +
  geom_point(
    data = volcano_data[
      display_class ==
        "plasma_higher"
    ],
    alpha = 0.45,
    size = 0.50,
    colour = "#d7301f"
  ) +
  geom_vline(
    xintercept = c(
      -minimum_absolute_median_delta_beta,
      minimum_absolute_median_delta_beta
    ),
    linetype = 2
  ) +
  geom_hline(
    yintercept =
      -log10(
        fdr_threshold
      ),
    linetype = 2
  ) +
  labs(
    title =
      "SOLID tissue versus plasma differential methylation",
    subtitle = paste0(
      "Primary candidates: FDR < ",
      fdr_threshold,
      " and |median delta Beta| >= ",
      minimum_absolute_median_delta_beta,
      "\nDelta Beta = plasma - tissue"
    ),
    x =
      "Median delta Beta (plasma - tissue)",
    y =
      expression(-log[10](FDR))
  ) +
  theme_bw(
    base_size = 12
  )

save_plot(
  volcano_plot,
  file.path(
    figure_dir,
    "SOLID_DMR_volcano.png"
  ),
  width = 11,
  height = 8
)

##------------------------------------------------------------
## 9. DELTA-BETA DISTRIBUTION AND DIRECTION SUMMARY
##------------------------------------------------------------

message_header(
  "CREATING EFFECT-DISTRIBUTION FIGURES"
)

distribution_sample <- sample_background_points(
  results,
  maximum_background_points
)

delta_distribution <- ggplot(
  distribution_sample,
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
      -minimum_absolute_median_delta_beta,
      minimum_absolute_median_delta_beta
    ),
    linetype = 2
  ) +
  labs(
    title =
      "Distribution of SOLID paired methylation effects",
    subtitle =
      "Median delta Beta is plasma minus tissue",
    x =
      "Median delta Beta",
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

candidate_direction_summary <- candidate_results[
  ,
  .(
    n_regions = .N
  ),
  by = candidate_direction
]

candidate_direction_summary[
  ,
  percentage :=
    100 *
      n_regions /
      sum(n_regions)
]

candidate_direction_summary[
  ,
  display_direction :=
    factor(
      candidate_direction,
      levels = c(
        "tissue_higher",
        "plasma_higher"
      ),
      labels = c(
        "Tissue higher",
        "Plasma higher"
      )
    )
]

candidate_direction_summary[
  ,
  display_label := paste0(
    comma(n_regions),
    "\n(",
    sprintf(
      "%.1f%%",
      percentage
    ),
    ")"
  )
]

direction_plot <- ggplot(
  candidate_direction_summary,
  aes(
    x = display_direction,
    y = percentage
  )
) +
  geom_col(
    width = 0.65,
    fill = "grey40"
  ) +
  geom_text(
    aes(
      label = display_label
    ),
    vjust = -0.30,
    size = 4
  ) +
  scale_y_continuous(
    labels = function(x) {
      paste0(x, "%")
    },
    expand = expansion(
      mult = c(0, 0.15)
    )
  ) +
  labs(
    title =
      "Direction of SOLID differential methylation candidates",
    subtitle = paste0(
      "FDR < ",
      fdr_threshold,
      " and |median delta Beta| >= ",
      minimum_absolute_median_delta_beta
    ),
    x = NULL,
    y =
      "Percentage of candidate regions"
  ) +
  theme_bw(
    base_size = 12
  )

save_plot(
  direction_plot,
  file.path(
    figure_dir,
    "SOLID_DMR_candidate_direction_percentage.png"
  ),
  width = 8,
  height = 7
)

fwrite(
  candidate_direction_summary,
  file.path(
    annotation_dir,
    "SOLID_DMR_candidate_direction_summary.tsv"
  ),
  sep = "\t"
)

##------------------------------------------------------------
## 10. SELECT REPRESENTATIVE TOP REGIONS
##------------------------------------------------------------

message_header(
  "SELECTING REPRESENTATIVE TOP DMR REGIONS"
)

top_plasma_regions <- head(
  plasma_higher_results,
  n_heatmap_regions_per_direction
)

top_tissue_regions <- head(
  tissue_higher_results,
  n_heatmap_regions_per_direction
)

top_heatmap_regions <- rbindlist(
  list(
    top_tissue_regions,
    top_plasma_regions
  ),
  use.names = TRUE,
  fill = TRUE
)

top_heatmap_regions <- top_heatmap_regions[
  !duplicated(region_id)
]

if (nrow(top_heatmap_regions) < 2L) {
  stop(
    "Too few candidate regions for heatmap generation."
  )
}

fwrite(
  top_heatmap_regions,
  file.path(
    annotation_dir,
    "SOLID_representative_heatmap_DMRs.tsv"
  ),
  sep = "\t"
)

##------------------------------------------------------------
## 11. HEATMAP OF WITHIN-PATIENT DELTA BETA
##------------------------------------------------------------

message_header(
  "CREATING DELTA-BETA HEATMAP"
)

heatmap_indices <- match(
  top_heatmap_regions$region_id,
  rownames(delta_beta)
)

stopifnot(
  !anyNA(heatmap_indices)
)

heatmap_matrix <- delta_beta[
  heatmap_indices,
  ,
  drop = FALSE
]

## Clip extreme values for display only
heatmap_matrix_display <- pmax(
  pmin(
    heatmap_matrix,
    heatmap_clip_delta_beta
  ),
  -heatmap_clip_delta_beta
)

row_labels <- paste0(
  top_heatmap_regions$region_id,
  " | ",
  ifelse(
    top_heatmap_regions$median_delta_beta > 0,
    "Plasma higher",
    "Tissue higher"
  )
)

rownames(heatmap_matrix_display) <-
  row_labels

heatmap_colours <- circlize::colorRamp2(
  c(
    -heatmap_clip_delta_beta,
    0,
    heatmap_clip_delta_beta
  ),
  c(
    "#2b8cbe",
    "white",
    "#d7301f"
  )
)

row_split <- factor(
  ifelse(
    top_heatmap_regions$median_delta_beta > 0,
    "Plasma higher",
    "Tissue higher"
  ),
  levels = c(
    "Tissue higher",
    "Plasma higher"
  )
)

heatmap_object <- ComplexHeatmap::Heatmap(
  heatmap_matrix_display,
  name = "Delta Beta",
  col = heatmap_colours,
  row_split = row_split,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names =
    nrow(heatmap_matrix_display) <= 80L,
  show_column_names = TRUE,
  column_title =
    "SOLID within-patient methylation differences",
  column_names_rot = 45,
  na_col = "grey90",
  heatmap_legend_param = list(
    title =
      "Plasma - tissue"
  )
)

png(
  filename = file.path(
    figure_dir,
    "SOLID_DMR_delta_beta_heatmap.png"
  ),
  width = 2600,
  height = 2400,
  res = 300
)

ComplexHeatmap::draw(
  heatmap_object,
  heatmap_legend_side = "right"
)

dev.off()

pdf(
  file = file.path(
    figure_dir,
    "SOLID_DMR_delta_beta_heatmap.pdf"
  ),
  width = 11,
  height = 10
)

ComplexHeatmap::draw(
  heatmap_object,
  heatmap_legend_side = "right"
)

dev.off()

##------------------------------------------------------------
## 12. PAIRED TISSUE-VERSUS-PLASMA PLOTS FOR TOP REGIONS
##------------------------------------------------------------

message_header(
  "CREATING TOP-REGION PAIRED PLOTS"
)

top_paired_regions <- rbindlist(
  list(
    head(
      plasma_higher_results,
      n_paired_plot_regions_per_direction
    ),
    head(
      tissue_higher_results,
      n_paired_plot_regions_per_direction
    )
  ),
  use.names = TRUE,
  fill = TRUE
)

top_paired_regions <- top_paired_regions[
  !duplicated(region_id)
]

paired_plot_data <- rbindlist(
  lapply(
    seq_len(
      nrow(top_paired_regions)
    ),
    function(i) {
      region_id_i <- top_paired_regions$
        region_id[i]

      matrix_row <- match(
        region_id_i,
        rownames(tissue_beta)
      )

      tissue_values <- tissue_beta[
        matrix_row,
      ]

      plasma_values <- plasma_beta[
        matrix_row,
      ]

      rbindlist(
        list(
          data.table(
            region_id =
              region_id_i,
            patient_id =
              colnames(tissue_beta),
            sample_type =
              "Tissue",
            beta =
              as.numeric(tissue_values)
          ),
          data.table(
            region_id =
              region_id_i,
            patient_id =
              colnames(plasma_beta),
            sample_type =
              "Plasma",
            beta =
              as.numeric(plasma_values)
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
        "Tissue",
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
    scales = "free_y",
    ncol = 2
  ) +
  coord_cartesian(
    ylim = c(0, 1)
  ) +
  labs(
    title =
      "Top SOLID paired tissue-plasma DMRs",
    subtitle =
      "Each line connects tissue and plasma from the same patient",
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
    "SOLID_DMR_top_regions_paired_beta.png"
  ),
  width = 12,
  height = 16
)

fwrite(
  paired_plot_data,
  file.path(
    annotation_dir,
    "SOLID_DMR_top_regions_paired_beta_values.tsv"
  ),
  sep = "\t"
)

##------------------------------------------------------------
## 13. MANHATTAN-STYLE PLOT
##------------------------------------------------------------

message_header(
  "CREATING MANHATTAN-STYLE PLOT"
)

required_coordinate_columns <- c(
  "chr",
  "start",
  "end"
)

missing_coordinate_columns <- setdiff(
  required_coordinate_columns,
  names(results)
)

if (length(missing_coordinate_columns) > 0L) {
  warning(
    "Skipping Manhattan plot and DMR block merging because ",
    "these coordinate columns are missing: ",
    paste(missing_coordinate_columns, collapse = ", ")
  )
} else {

  manhattan_data <- copy(
    results
  )

  manhattan_data[
    ,
    chr_plot :=
      standardize_chr(chr)
  ]

  manhattan_data <- manhattan_data[
    !is.na(chr_plot) &
      !is.na(start)
  ]

  manhattan_data[
    ,
    chr_numeric :=
      as.integer(chr_plot)
  ]

  chromosome_offsets <- manhattan_data[
    ,
    .(
      chromosome_length =
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
    chromosome_offsets,
    chr_numeric
  )

  chromosome_offsets[
    ,
    cumulative_offset :=
      c(
        0,
        head(
          cumsum(chromosome_length),
          -1L
        )
      )
  ]

  manhattan_data <- merge(
    manhattan_data,
    chromosome_offsets[
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

  manhattan_data[
    ,
    genomic_position :=
      start +
      cumulative_offset
  ]

  chromosome_centres <- manhattan_data[
    ,
    .(
      centre =
        (
          min(genomic_position) +
            max(genomic_position)
        ) /
        2
    ),
    by = chr_plot
  ]

  chromosome_centres[
    ,
    chr_numeric :=
      as.integer(chr_plot)
  ]

  setorder(
    chromosome_centres,
    chr_numeric
  )

  manhattan_background <- sample_background_points(
    manhattan_data[
      primary_candidate == FALSE
    ],
    maximum_background_points
  )

  manhattan_candidates <- manhattan_data[
    primary_candidate == TRUE
  ]

  manhattan_display <- rbindlist(
    list(
      manhattan_background[
        ,
        display_class :=
          "not_candidate"
      ],
      manhattan_candidates[
        ,
        display_class :=
          candidate_direction
      ]
    ),
    use.names = TRUE,
    fill = TRUE
  )

  manhattan_display[
    ,
    minus_log10_FDR :=
      -log10(
        pmax(
          FDR,
          .Machine$double.xmin
        )
      )
  ]

  manhattan_plot <- ggplot(
    manhattan_display,
    aes(
      x = genomic_position,
      y = minus_log10_FDR
    )
  ) +
    geom_point(
      data = manhattan_display[
        display_class ==
          "not_candidate"
      ],
      alpha = 0.15,
      size = 0.30,
      colour = "grey60"
    ) +
    geom_point(
      data = manhattan_display[
        display_class ==
          "tissue_higher"
      ],
      alpha = 0.35,
      size = 0.40,
      colour = "#2b8cbe"
    ) +
    geom_point(
      data = manhattan_display[
        display_class ==
          "plasma_higher"
      ],
      alpha = 0.45,
      size = 0.45,
      colour = "#d7301f"
    ) +
    geom_hline(
      yintercept =
        -log10(
          fdr_threshold
        ),
      linetype = 2
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
        "Genome-wide SOLID tissue-plasma differential methylation",
      subtitle =
        "Delta Beta = plasma - tissue",
      x =
        "Chromosome",
      y =
        expression(-log[10](FDR))
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.major.x =
        element_blank()
    )

  save_plot(
    manhattan_plot,
    file.path(
      figure_dir,
      "SOLID_DMR_manhattan.png"
    ),
    width = 16,
    height = 7
  )

  ##----------------------------------------------------------
  ## 14. MERGE ADJACENT SIGNIFICANT WINDOWS
  ##----------------------------------------------------------

  message_header(
    "MERGING ADJACENT SIGNIFICANT WINDOWS"
  )

  block_input <- candidate_results[
    !is.na(chr) &
      !is.na(start) &
      !is.na(end),
    .(
      region_id,
      chr =
        as.character(chr),
      start =
        as.integer(start),
      end =
        as.integer(end),
      FDR,
      p_value,
      median_delta_beta,
      mean_delta_beta,
      candidate_direction
    )
  ]

  setorder(
    block_input,
    chr,
    candidate_direction,
    start,
    end
  )

  block_input[
    ,
    previous_end :=
      shift(end),
    by = .(
      chr,
      candidate_direction
    )
  ]

  block_input[
    ,
    new_block :=
      is.na(previous_end) |
        start >
          previous_end +
            maximum_gap_between_windows,
    by = .(
      chr,
      candidate_direction
    )
  ]

  block_input[
    ,
    block_number :=
      cumsum(new_block),
    by = .(
      chr,
      candidate_direction
    )
  ]

  dmr_blocks <- block_input[
    ,
    .(
      block_start =
        min(start),
      block_end =
        max(end),
      block_width =
        max(end) -
        min(start) +
        1L,
      n_windows =
        .N,
      minimum_FDR =
        min(
          FDR,
          na.rm = TRUE
        ),
      median_window_FDR =
        median(
          FDR,
          na.rm = TRUE
        ),
      median_delta_beta =
        median(
          median_delta_beta,
          na.rm = TRUE
        ),
      mean_delta_beta =
        mean(
          mean_delta_beta,
          na.rm = TRUE
        ),
      maximum_absolute_delta_beta =
        max(
          abs(
            median_delta_beta
          ),
          na.rm = TRUE
        ),
      first_region_id =
        first(region_id),
      last_region_id =
        last(region_id)
    ),
    by = .(
      chr,
      candidate_direction,
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
    retained_block == TRUE
  ]

  retained_blocks[
    ,
    block_id := paste0(
      "SOLID_DMR_",
      sprintf(
        "%07d",
        seq_len(.N)
      )
    )
  ]

  setcolorder(
    retained_blocks,
    c(
      "block_id",
      setdiff(
        names(retained_blocks),
        "block_id"
      )
    )
  )

  fwrite(
    dmr_blocks,
    file.path(
      annotation_dir,
      "SOLID_all_merged_DMR_blocks.tsv.gz"
    ),
    sep = "\t",
    compress = "gzip"
  )

  fwrite(
    retained_blocks,
    file.path(
      annotation_dir,
      "SOLID_retained_merged_DMR_blocks.tsv"
    ),
    sep = "\t"
  )

  cat(
    "Merged candidate blocks:",
    format(nrow(dmr_blocks), big.mark = ","),
    "\n"
  )

  cat(
    "Retained multi-window blocks:",
    format(nrow(retained_blocks), big.mark = ","),
    "\n"
  )
}

##------------------------------------------------------------
## 15. CREATE FIGURE AND ANALYSIS SUMMARY
##------------------------------------------------------------

message_header(
  "CREATING FIGURE SUMMARY"
)

figure_summary <- data.table(
  item = c(
    "regions_tested",
    "primary_DMR_candidates",
    "plasma_higher_candidates",
    "tissue_higher_candidates",
    "heatmap_regions",
    "paired_plot_regions",
    "FDR_threshold",
    "absolute_median_delta_beta_threshold",
    "matched_patients",
    "primary_effect_size_cutoff"
  ),
  value = c(
    nrow(results),
    nrow(candidate_results),
    nrow(plasma_higher_results),
    nrow(tissue_higher_results),
    nrow(top_heatmap_regions),
    nrow(top_paired_regions),
    fdr_threshold,
    minimum_absolute_median_delta_beta,
    ncol(delta_beta),
    minimum_absolute_median_delta_beta
  )
)

print(
  figure_summary
)

fwrite(
  figure_summary,
  file.path(
    annotation_dir,
    "SOLID_DMR_figure_summary.tsv"
  ),
  sep = "\t"
)

##------------------------------------------------------------
## 16. FINAL VALIDATION
##------------------------------------------------------------

required_figure_files <- c(
  "SOLID_DMR_volcano.png",
  "SOLID_DMR_delta_beta_distribution.png",
  "SOLID_DMR_candidate_direction_percentage.png",
  "SOLID_DMR_delta_beta_heatmap.png",
  "SOLID_DMR_delta_beta_heatmap.pdf",
  "SOLID_DMR_top_regions_paired_beta.png"
)

missing_figure_files <- required_figure_files[
  !file.exists(
    file.path(
      figure_dir,
      required_figure_files
    )
  )
]

if (length(missing_figure_files) > 0L) {
  stop(
    "Expected figure file(s) were not created:\n",
    paste(
      missing_figure_files,
      collapse = "\n"
    )
  )
}

cat(
  "\nScript 06 completed successfully.\n"
)

cat(
  "Primary DMR candidates:",
  format(nrow(candidate_results), big.mark = ","),
  "\n"
)

cat(
  "Figures saved to:\n",
  figure_dir,
  "\n"
)

cat(
  "Prioritized tables saved to:\n",
  annotation_dir,
  "\n"
)