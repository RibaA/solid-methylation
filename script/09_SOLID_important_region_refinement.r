############################################################
## 09_SOLID_important_region_refinement.r
##
## SOLID IMPORTANT-REGION REFINEMENT
##
## Purpose:
##
##   1. Preserve all Script 07/08 discovery results
##
##   2. Create a stricter Tier-A set of highly compelling
##      SOLID tumor-plasma DMRs
##
##   3. Add technical-support information:
##        - EPIC probe support
##        - plasma coverage
##        - plasma CpG support
##
##   4. Add multi-window DMR-block support
##
##   5. Check chromosome representation explicitly
##      (including chr14-22/X/Y)
##
##   6. Rank important regions transparently for:
##        - later biological annotation
##        - OCTANE integration
##        - external validation
##
##
## Tier-A definition:
##
##   FDR < 0.05
##
##   |median Delta Beta| >= 0.30
##
##   direction consistency = 100%
##
##   same-direction |Delta Beta| >= 0.10
##      in >=80% of valid patients
##
##   valid matched pairs = 13/13
##
##
## Delta Beta:
##
##   Plasma - Tumor
##
## Positive:
##   Plasma higher
##
## Negative:
##   Tumor higher
##
##
## IMPORTANT:
##
##   - This script DOES NOT replace Script 08.
##
##   - It creates a stricter refinement layer.
##
##   - No new statistical testing is performed.
##
##   - Script 07 limma results remain authoritative.
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


############################################################
## Tier-A thresholds
############################################################

fdr_threshold <- 0.05

tierA_abs_delta_beta <- 0.30

tierA_direction_consistency <- 100

tierA_effect_consistency_010 <- 80

tierA_required_valid_pairs <- 13L


############################################################
## Output/display settings
############################################################

top_n_regions <- 100L

top_n_plot_regions_per_direction <- 10L


############################################################
## Standard chromosomes
############################################################

standard_chr_levels <- c(
  as.character(1:22),
  "X",
  "Y"
)


############################################################
## 2. REQUIRED PACKAGES
############################################################

required_packages <- c(
  "data.table",
  "ggplot2",
  "scales"
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

  library(data.table)

  library(ggplot2)

  library(scales)
})


############################################################
## 3. HELPER FUNCTIONS
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


############################################################
## Windows-safe TSV reader
##
## data.table::fread() occasionally treats .gz paths as
## shell commands under some Windows configurations.
##
## For .gz:
##   use base R gzfile() + read.delim()
##
## For regular TSV:
##   use fread()
############################################################

read_tsv_safe <- function(file) {

  if (!file.exists(file)) {

    stop(
      "Input file does not exist:\n",
      file
    )
  }


  if (grepl(
    "\\.gz$",
    file,
    ignore.case = TRUE
  )) {

    con <- gzfile(
      file,
      open = "rt"
    )


    on.exit(
      close(con),
      add = TRUE
    )


    out <- read.delim(
      con,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )


    out <- data.table::as.data.table(
      out
    )

  } else {

    out <- data.table::fread(
      file
    )
  }


  out
}


############################################################
## Save ggplot as both PNG and PDF
############################################################

save_plot_both <- function(
    plot_object,
    basename,
    width,
    height
) {

  ggsave(
    filename = paste0(
      basename,
      ".png"
    ),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300,
    limitsize = FALSE
  )


  ggsave(
    filename = paste0(
      basename,
      ".pdf"
    ),
    plot = plot_object,
    width = width,
    height = height,
    limitsize = FALSE
  )
}


############################################################
## 4. DIRECTORIES
############################################################

dmr_dir <- file.path(
  project_dir,
  "result",
  "03_matched_tissue_plasma",
  "DMR"
)


priority_dir <- file.path(
  dmr_dir,
  "prioritization"
)


refinement_dir <- file.path(
  dmr_dir,
  "refinement"
)


figure_dir <- file.path(
  refinement_dir,
  "figures"
)


table_dir <- file.path(
  refinement_dir,
  "tables"
)


dir.create(
  refinement_dir,
  recursive = TRUE,
  showWarnings = FALSE
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


############################################################
## 5. INPUT FILES
############################################################

all_results_file <- file.path(
  priority_dir,
  "SOLID_all_DMR_results_with_consistency.tsv.gz"
)


dmr_object_file <- file.path(
  dmr_dir,
  "SOLID_tissue_vs_plasma_DMR_object.rds"
)


block_file <- file.path(
  priority_dir,
  "SOLID_retained_multiwindow_DMR_blocks.tsv"
)


############################################################
## Validate mandatory input files
############################################################

if (!file.exists(all_results_file)) {

  stop(
    "Script-08 consistency results file not found:\n",
    all_results_file
  )
}


if (!file.exists(dmr_object_file)) {

  stop(
    "Script-07 DMR object not found:\n",
    dmr_object_file
  )
}


############################################################
## 6. LOAD SCRIPT-08 RESULTS
############################################################

message_header(
  "SOLID IMPORTANT-REGION REFINEMENT"
)


cat(
  "Loading Script-08 consistency results...\n"
)


results <- read_tsv_safe(
  all_results_file
)


cat(
  "Consistency results loaded:",
  format(
    nrow(results),
    big.mark = ","
  ),
  "regions\n"
)


############################################################
## 7. VALIDATE SCRIPT-08 TABLE
############################################################

required_script08_columns <- c(
  "region_id",
  "chr",
  "start",
  "end",
  "FDR",
  "median_delta_beta",
  "abs_median_delta_beta",
  "direction_consistency_pct",
  "effect_consistency_0.10_pct",
  "importance_tier",
  "n_valid_pairs",
  "DMR_direction"
)


missing_script08_columns <- setdiff(
  required_script08_columns,
  names(results)
)


if (length(missing_script08_columns) > 0L) {

  stop(
    "Script-08 table is missing required column(s): ",
    paste(
      missing_script08_columns,
      collapse = ", "
    )
  )
}


stopifnot(
  nrow(results) == 124961L,
  !anyNA(results$region_id),
  !anyDuplicated(results$region_id)
)


cat(
  "Script-08 consistency table validation: PASS\n"
)


############################################################
## 8. LOAD SCRIPT-07 DMR OBJECT
############################################################

cat(
  "\nLoading Script-07 DMR object...\n"
)


dmr_object <- readRDS(
  dmr_object_file
)


required_dmr_components <- c(
  "tissue_beta",
  "plasma_beta",
  "delta_beta",
  "valid_pair_mask",
  "patient_metadata"
)


missing_dmr_components <- setdiff(
  required_dmr_components,
  names(dmr_object)
)


if (length(missing_dmr_components) > 0L) {

  stop(
    "DMR object missing component(s): ",
    paste(
      missing_dmr_components,
      collapse = ", "
    )
  )
}


tissue_beta <- dmr_object$tissue_beta

plasma_beta <- dmr_object$plasma_beta

delta_beta <- dmr_object$delta_beta

valid_pair_mask <- dmr_object$valid_pair_mask

patient_metadata <- as.data.table(
  dmr_object$patient_metadata
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
    dim(valid_pair_mask)
  ),

  nrow(tissue_beta) ==
    nrow(results),

  ncol(tissue_beta) ==
    13L,

  setequal(
    rownames(tissue_beta),
    results$region_id
  )
)


############################################################
## Restore exact matrix order
############################################################

results <- results[
  match(
    rownames(tissue_beta),
    region_id
  )
]


stopifnot(
  identical(
    results$region_id,
    rownames(tissue_beta)
  )
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
    patient_metadata$patient_id
  )
)


cat(
  "DMR object validation: PASS\n"
)


############################################################
## 9. LOAD MULTI-WINDOW BLOCKS
############################################################

if (file.exists(block_file)) {

  blocks <- read_tsv_safe(
    block_file
  )


  cat(
    "\nMulti-window blocks loaded:",
    format(
      nrow(blocks),
      big.mark = ","
    ),
    "\n"
  )

} else {

  warning(
    "No retained multi-window block file found:\n",
    block_file
  )


  blocks <- data.table()
}


############################################################
## 10. STANDARDIZE CHROMOSOME LABELS
############################################################

results[
  ,
  chr_clean :=
    sub(
      "^chr",
      "",
      as.character(chr),
      ignore.case = TRUE
    )
]


results[
  ,
  chr_order :=
    match(
      chr_clean,
      standard_chr_levels
    )
]


############################################################
## 11. DEFINE TIER-A REGIONS
############################################################

message_header(
  "DEFINING TIER-A IMPORTANT REGIONS"
)


############################################################
## Use >=100 rather than ==100 to protect against floating
## representation, although expected values are exactly 100.
############################################################

results[
  ,
  TierA :=
    FDR <
      fdr_threshold &

    abs_median_delta_beta >=
      tierA_abs_delta_beta &

    direction_consistency_pct >=
      tierA_direction_consistency &

    effect_consistency_0.10_pct >=
      tierA_effect_consistency_010 &

    n_valid_pairs >=
      tierA_required_valid_pairs
]


tierA <- copy(
  results[
    TierA == TRUE
  ]
)


if (nrow(tierA) == 0L) {

  stop(
    "No regions satisfy the Tier-A definition."
  )
}


cat(
  "Tier-A regions:",
  format(
    nrow(tierA),
    big.mark = ","
  ),
  "\n"
)


############################################################
## 12. DIRECTION-SPECIFIC TIER-A SETS
############################################################

tierA_tumor <- copy(
  tierA[
    DMR_direction ==
      "Tumor_higher"
  ]
)


tierA_plasma <- copy(
  tierA[
    DMR_direction ==
      "Plasma_higher"
  ]
)


cat(
  "  Tumor higher:",
  format(
    nrow(tierA_tumor),
    big.mark = ","
  ),
  "\n"
)


cat(
  "  Plasma higher:",
  format(
    nrow(tierA_plasma),
    big.mark = ","
  ),
  "\n"
)


############################################################
## 13. TECHNICAL SUPPORT COLUMNS
############################################################

technical_columns <- intersect(
  c(
    "n_EPIC_probes",
    "median_plasma_coverage",
    "median_plasma_covered_cpgs",
    "n_tissue_samples_observed",
    "tissue_samples_observed_pct"
  ),
  names(tierA)
)


cat(
  "\nTechnical-support columns available:\n"
)


if (length(technical_columns) > 0L) {

  print(
    technical_columns
  )

} else {

  cat(
    "None of the optional technical-support columns were found.\n"
  )
}


############################################################
## 14. INITIALIZE BLOCK-SUPPORT COLUMNS
############################################################

tierA[
  ,
  in_multiwindow_block :=
    FALSE
]


tierA[
  ,
  block_id :=
    NA_character_
]


tierA[
  ,
  block_n_windows :=
    NA_integer_
]


tierA[
  ,
  block_width :=
    NA_integer_
]


############################################################
## 15. MAP TIER-A REGIONS TO MULTI-WINDOW BLOCKS
##
## Uses interval overlap rather than looping over every block.
############################################################

if (nrow(blocks) > 0L) {

  required_block_columns <- c(
    "block_id",
    "chr",
    "block_start",
    "block_end",
    "n_windows",
    "DMR_direction"
  )


  missing_block_columns <- setdiff(
    required_block_columns,
    names(blocks)
  )


  if (length(missing_block_columns) > 0L) {

    stop(
      "Block table missing required column(s): ",
      paste(
        missing_block_columns,
        collapse = ", "
      )
    )
  }


  ##########################################################
  ## Prepare Tier-A intervals
  ##########################################################

  tier_intervals <- tierA[
    ,
    .(
      tier_row =
        .I,

      chr =
        as.character(chr),

      start =
        as.integer(start),

      end =
        as.integer(end),

      region_direction =
        as.character(
          DMR_direction
        )
    )
  ]


  ##########################################################
  ## Prepare block intervals
  ##########################################################

  block_intervals <- blocks[
    ,
    .(
      block_id =
        as.character(
          block_id
        ),

      chr =
        as.character(chr),

      start =
        as.integer(
          block_start
        ),

      end =
        as.integer(
          block_end
        ),

      block_n_windows =
        as.integer(
          n_windows
        ),

      block_direction =
        as.character(
          DMR_direction
        ),

      block_width =
        as.integer(
          block_width
        )
    )
  ]


  ##########################################################
  ## foverlaps requires interval keys
  ##########################################################

  setkey(
    block_intervals,
    chr,
    start,
    end
  )


  overlaps <- foverlaps(
    tier_intervals,
    block_intervals,
    by.x = c(
      "chr",
      "start",
      "end"
    ),
    by.y = c(
      "chr",
      "start",
      "end"
    ),
    type = "within",
    nomatch = 0L
  )


  ##########################################################
  ## Require same DMR direction
  ##########################################################

  overlaps <- overlaps[
    region_direction ==
      block_direction
  ]


  ##########################################################
  ## A region should normally belong to one block.
  ## If more than one overlap exists, prioritize the
  ## block containing the largest number of windows.
  ##########################################################

  if (nrow(overlaps) > 0L) {

    setorder(
      overlaps,
      tier_row,
      -block_n_windows
    )


    overlap_best <- overlaps[
      ,
      .SD[1],
      by = tier_row
    ]


    tierA[
      overlap_best$tier_row,
      `:=`(
        in_multiwindow_block =
          TRUE,

        block_id =
          overlap_best$block_id,

        block_n_windows =
          overlap_best$block_n_windows,

        block_width =
          overlap_best$block_width
      )
    ]
  }
}


tierA_block_supported <- copy(
  tierA[
    in_multiwindow_block ==
      TRUE
  ]
)


cat(
  "\nTier-A regions in multi-window blocks:",
  format(
    nrow(
      tierA_block_supported
    ),
    big.mark = ","
  ),
  "\n"
)


############################################################
## 16. CREATE TRANSPARENT RANKING COMPONENTS
############################################################

############################################################
## Block support
############################################################

tierA[
  ,
  block_supported_rank :=
    as.integer(
      in_multiwindow_block
    )
]


############################################################
## EPIC probe support
############################################################

if ("n_EPIC_probes" %in% names(tierA)) {

  tierA[
    ,
    EPIC_probe_rank :=
      fifelse(
        is.na(
          n_EPIC_probes
        ),
        0,
        pmin(
          n_EPIC_probes,
          10
        )
      )
  ]

} else {

  tierA[
    ,
    EPIC_probe_rank :=
      0
  ]
}


############################################################
## Plasma coverage support
############################################################

if ("median_plasma_coverage" %in% names(tierA)) {

  tierA[
    ,
    plasma_coverage_rank :=
      fifelse(
        is.na(
          median_plasma_coverage
        ),
        0,
        median_plasma_coverage
      )
  ]

} else {

  tierA[
    ,
    plasma_coverage_rank :=
      0
  ]
}


############################################################
## Plasma covered CpG support
############################################################

if ("median_plasma_covered_cpgs" %in% names(tierA)) {

  tierA[
    ,
    plasma_cpg_rank :=
      fifelse(
        is.na(
          median_plasma_covered_cpgs
        ),
        0,
        median_plasma_covered_cpgs
      )
  ]

} else {

  tierA[
    ,
    plasma_cpg_rank :=
      0
  ]
}


############################################################
## 17. RANK TIER-A REGIONS
##
## Priority hierarchy:
##
##   1. multi-window block support
##   2. effect recurrence >=0.10
##   3. absolute median Delta Beta
##   4. direction consistency
##   5. EPIC probe support
##   6. plasma CpG support
##   7. plasma coverage
##   8. FDR
############################################################

setorder(
  tierA,

  -block_supported_rank,

  -effect_consistency_0.10_pct,

  -abs_median_delta_beta,

  -direction_consistency_pct,

  -EPIC_probe_rank,

  -plasma_cpg_rank,

  -plasma_coverage_rank,

  FDR
)


tierA[
  ,
  TierA_rank :=
    seq_len(.N)
]


############################################################
## Recreate direction-specific tables AFTER ranking
############################################################

tierA_tumor <- copy(
  tierA[
    DMR_direction ==
      "Tumor_higher"
  ]
)


tierA_plasma <- copy(
  tierA[
    DMR_direction ==
      "Plasma_higher"
  ]
)


tierA_block_supported <- copy(
  tierA[
    in_multiwindow_block ==
      TRUE
  ]
)


############################################################
## 18. SAVE TIER-A TABLES
############################################################

message_header(
  "SAVING TIER-A TABLES"
)


fwrite(
  tierA,
  file.path(
    table_dir,
    "SOLID_TierA_important_regions.tsv.gz"
  ),
  sep = "\t",
  compress = "gzip"
)


fwrite(
  head(
    tierA,
    top_n_regions
  ),
  file.path(
    table_dir,
    "SOLID_TierA_top100_regions.tsv"
  ),
  sep = "\t"
)


fwrite(
  tierA_block_supported,
  file.path(
    table_dir,
    "SOLID_TierA_block_supported_regions.tsv"
  ),
  sep = "\t"
)


fwrite(
  tierA_tumor,
  file.path(
    table_dir,
    "SOLID_TierA_tumor_higher_regions.tsv"
  ),
  sep = "\t"
)


fwrite(
  tierA_plasma,
  file.path(
    table_dir,
    "SOLID_TierA_plasma_higher_regions.tsv"
  ),
  sep = "\t"
)


############################################################
## 19. CHROMOSOME REPRESENTATION CHECK
############################################################

message_header(
  "CHECKING CHROMOSOME REPRESENTATION"
)


############################################################
## Build complete table so missing chromosomes appear as 0
############################################################

chromosome_template <- data.table(
  chr_clean =
    standard_chr_levels,
  chr_order =
    seq_along(
      standard_chr_levels
    )
)


chromosome_observed <- results[
  chr_clean %in%
    standard_chr_levels,
  .(

    regions_tested =
      .N,

    FDR_significant =
      sum(
        FDR <
          fdr_threshold,
        na.rm = TRUE
      ),

    strong_consistent =
      sum(
        as.character(
          importance_tier
        ) ==
          "Strong_consistent",
        na.rm = TRUE
      ),

    TierA_regions =
      sum(
        TierA,
        na.rm = TRUE
      )

  ),
  by = chr_clean
]


chromosome_check <- merge(
  chromosome_template,
  chromosome_observed,
  by = "chr_clean",
  all.x = TRUE,
  sort = FALSE
)


for (
  v in c(
    "regions_tested",
    "FDR_significant",
    "strong_consistent",
    "TierA_regions"
  )
) {

  set(
    chromosome_check,
    which(
      is.na(
        chromosome_check[[v]]
      )
    ),
    v,
    0
  )
}


setorder(
  chromosome_check,
  chr_order
)


fwrite(
  chromosome_check,
  file.path(
    table_dir,
    "SOLID_chromosome_representation_check.tsv"
  ),
  sep = "\t"
)


print(
  chromosome_check
)


############################################################
## Specifically inspect chromosomes that seemed absent
## from the previous genome-wide figure
############################################################

later_chr <- chromosome_check[
  chr_clean %in%
    c(
      as.character(
        14:22
      ),
      "X",
      "Y"
    )
]


cat(
  "\nChromosomes 14-22/X/Y:\n"
)


print(
  later_chr
)


############################################################
## 20. MULTI-WINDOW BLOCK SUMMARY
############################################################

message_header(
  "SUMMARIZING MULTI-WINDOW BLOCKS"
)


if (nrow(blocks) > 0L) {

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

      max_windows =
        max(
          n_windows,
          na.rm = TRUE
        ),

      median_block_width =
        median(
          block_width,
          na.rm = TRUE
        ),

      max_block_width =
        max(
          block_width,
          na.rm = TRUE
        )

    ),
    by = DMR_direction
  ]


  setorder(
    block_summary,
    DMR_direction
  )


  fwrite(
    block_summary,
    file.path(
      table_dir,
      "SOLID_multiwindow_block_summary.tsv"
    ),
    sep = "\t"
  )


  print(
    block_summary
  )


  ##########################################################
  ## Largest blocks
  ##########################################################

  largest_blocks <- copy(
    blocks
  )


  setorder(
    largest_blocks,
    -n_windows,
    -block_width
  )


  fwrite(
    head(
      largest_blocks,
      100L
    ),
    file.path(
      table_dir,
      "SOLID_top100_multiwindow_DMR_blocks.tsv"
    ),
    sep = "\t"
  )

} else {

  block_summary <- data.table()
}


############################################################
## 21. TIER-A EFFECT-SIZE DISTRIBUTION
############################################################

message_header(
  "CREATING TIER-A FIGURES"
)


p1 <- ggplot(
  tierA,
  aes(
    x = median_delta_beta,
    fill = DMR_direction
  )
) +

  geom_histogram(
    bins = 80,
    alpha = 0.70,
    position = "identity"
  ) +

  geom_vline(
    xintercept = c(
      -tierA_abs_delta_beta,
      tierA_abs_delta_beta
    ),
    linetype = 2
  ) +

  labs(
    title =
      "SOLID Tier-A important DMRs",

    subtitle =
      paste0(
        "FDR < ",
        fdr_threshold,
        "; |median Δβ| ≥ ",
        tierA_abs_delta_beta,
        "; 100% directional consistency"
      ),

    x =
      "Median Δβ (Plasma - Tumor)",

    y =
      "Number of Tier-A regions",

    fill =
      "Direction"
  ) +

  theme_bw(
    base_size = 12
  ) +

  theme(
    plot.title =
      element_text(
        face = "bold"
      ),
    legend.position =
      "top"
  )


save_plot_both(
  p1,
  file.path(
    figure_dir,
    "SOLID_TierA_delta_beta_distribution"
  ),
  width = 9,
  height = 6
)


############################################################
## 22. TIER-A REGIONS BY CHROMOSOME
############################################################

tierA_chr <- tierA[
  chr_clean %in%
    standard_chr_levels,
  .(
    n_regions = .N
  ),
  by = chr_clean
]


tierA_chr <- merge(
  chromosome_template,
  tierA_chr,
  by = "chr_clean",
  all.x = TRUE,
  sort = FALSE
)


tierA_chr[
  is.na(
    n_regions
  ),
  n_regions := 0L
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


p2 <- ggplot(
  tierA_chr,
  aes(
    x = chr_factor,
    y = n_regions
  )
) +

  geom_col() +

  labs(
    title =
      "SOLID Tier-A regions by chromosome",

    x =
      "Chromosome",

    y =
      "Number of Tier-A regions"
  ) +

  theme_bw(
    base_size = 12
  ) +

  theme(
    plot.title =
      element_text(
        face = "bold"
      )
  )


save_plot_both(
  p2,
  file.path(
    figure_dir,
    "SOLID_TierA_regions_by_chromosome"
  ),
  width = 12,
  height = 6
)


############################################################
## 23. BLOCK-SUPPORT SUMMARY FIGURE
############################################################

block_support_summary <- tierA[
  ,
  .(
    n_regions = .N
  ),
  by = .(
    DMR_direction,
    in_multiwindow_block
  )
]


block_support_summary[
  ,
  block_support :=
    fifelse(
      in_multiwindow_block,
      "Multi-window block",
      "Isolated 1-kb region"
    )
]


p3 <- ggplot(
  block_support_summary,
  aes(
    x = DMR_direction,
    y = n_regions,
    fill = block_support
  )
) +

  geom_col(
    position = "stack"
  ) +

  scale_y_continuous(
    labels = comma
  ) +

  labs(
    title =
      "Block support among SOLID Tier-A regions",

    x = NULL,

    y =
      "Number of Tier-A regions",

    fill =
      "Regional support"
  ) +

  theme_bw(
    base_size = 12
  ) +

  theme(
    legend.position =
      "top"
  )


save_plot_both(
  p3,
  file.path(
    figure_dir,
    "SOLID_TierA_block_support"
  ),
  width = 9,
  height = 6
)


############################################################
## 24. BLOCK SIZE DISTRIBUTION
############################################################

if (nrow(blocks) > 0L) {

  p4 <- ggplot(
    blocks,
    aes(
      x = n_windows,
      fill = DMR_direction
    )
  ) +

    geom_histogram(
      binwidth = 1,
      alpha = 0.65,
      position = "identity"
    ) +

    scale_y_continuous(
      labels = comma
    ) +

    labs(
      title =
        "SOLID multi-window DMR block sizes",

      x =
        "Number of adjacent 1-kb windows",

      y =
        "Number of DMR blocks",

      fill =
        "Direction"
    ) +

    theme_bw(
      base_size = 12
    ) +

    theme(
      legend.position =
        "top"
    )


  save_plot_both(
    p4,
    file.path(
      figure_dir,
      "SOLID_multiwindow_block_size_distribution"
    ),
    width = 9,
    height = 6
  )
}


############################################################
## 25. SELECT TOP REPRESENTATIVE TIER-A REGIONS
##
## Use both directions so one direction does not dominate.
############################################################

top_tumor_plot <- head(
  tierA_tumor,
  top_n_plot_regions_per_direction
)


top_plasma_plot <- head(
  tierA_plasma,
  top_n_plot_regions_per_direction
)


top_plot_regions <- rbindlist(
  list(
    top_tumor_plot,
    top_plasma_plot
  ),
  use.names = TRUE,
  fill = TRUE
)


top_plot_regions <- top_plot_regions[
  !duplicated(
    region_id
  )
]


fwrite(
  top_plot_regions,
  file.path(
    table_dir,
    "SOLID_TierA_representative_top_regions.tsv"
  ),
  sep = "\t"
)


############################################################
## 26. TOP TIER-A PAIRED BETA PLOTS
############################################################

if (nrow(top_plot_regions) > 0L) {

  paired_data <- rbindlist(
    lapply(
      seq_len(
        nrow(
          top_plot_regions
        )
      ),
      function(i) {

        region_i <-
          top_plot_regions$region_id[i]


        idx <- match(
          region_i,
          rownames(
            tissue_beta
          )
        )


        if (is.na(idx)) {

          stop(
            "Region not found in Beta matrix: ",
            region_i
          )
        }


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


  paired_data[
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


  ##########################################################
  ## Preserve region order from ranking
  ##########################################################

  paired_data[
    ,
    region_id :=
      factor(
        region_id,
        levels =
          top_plot_regions$region_id
      )
  ]


  p5 <- ggplot(
    paired_data,
    aes(
      x = sample_type,
      y = beta,
      group = patient_id
    )
  ) +

    geom_line(
      alpha = 0.40,
      colour = "grey55"
    ) +

    geom_point(
      size = 1.5
    ) +

    facet_wrap(
      ~region_id,
      ncol = 4
    ) +

    coord_cartesian(
      ylim = c(
        0,
        1
      )
    ) +

    labs(
      title =
        "Representative SOLID Tier-A DMRs",

      subtitle =
        "Ten highest-ranked regions per direction",

      x = NULL,

      y =
        "Regional Beta"
    ) +

    theme_bw(
      base_size = 9
    ) +

    theme(
      strip.text =
        element_text(
          size = 7,
          face = "bold"
        )
    )


  save_plot_both(
    p5,
    file.path(
      figure_dir,
      "SOLID_TierA_top20_paired_beta"
    ),
    width = 14,
    height = 12
  )


  fwrite(
    paired_data,
    file.path(
      table_dir,
      "SOLID_TierA_top20_paired_beta_values.tsv"
    ),
    sep = "\t"
  )
}


############################################################
## 27. FINAL REFINEMENT SUMMARY
############################################################

summary_table <- data.table(

  item = c(
    "regions_tested",
    "FDR_significant_regions",
    "strong_consistent_script08_regions",
    "TierA_regions",
    "TierA_tumor_higher",
    "TierA_plasma_higher",
    "TierA_block_supported",
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
      as.character(
        results$importance_tier
      ) ==
        "Strong_consistent",
      na.rm = TRUE
    ),

    nrow(
      tierA
    ),

    nrow(
      tierA_tumor
    ),

    nrow(
      tierA_plasma
    ),

    nrow(
      tierA_block_supported
    ),

    nrow(
      blocks
    )
  )
)


fwrite(
  summary_table,
  file.path(
    table_dir,
    "SOLID_TierA_refinement_summary.tsv"
  ),
  sep = "\t"
)


############################################################
## 28. FINAL VALIDATION
############################################################

message_header(
  "FINAL VALIDATION"
)


required_final_objects <- c(
  "results",
  "tierA",
  "tierA_tumor",
  "tierA_plasma",
  "tierA_block_supported",
  "blocks",
  "chromosome_check",
  "summary_table"
)


missing_final_objects <- required_final_objects[
  !vapply(
    required_final_objects,
    exists,
    logical(1),
    inherits = FALSE
  )
]


if (length(missing_final_objects) > 0L) {

  stop(
    "Script 09 did not complete correctly.\n",
    "Missing object(s): ",
    paste(
      missing_final_objects,
      collapse = ", "
    )
  )
}


required_output_files <- c(

  file.path(
    table_dir,
    "SOLID_TierA_important_regions.tsv.gz"
  ),

  file.path(
    table_dir,
    "SOLID_TierA_top100_regions.tsv"
  ),

  file.path(
    table_dir,
    "SOLID_TierA_block_supported_regions.tsv"
  ),

  file.path(
    table_dir,
    "SOLID_chromosome_representation_check.tsv"
  ),

  file.path(
    table_dir,
    "SOLID_TierA_refinement_summary.tsv"
  )
)


missing_output_files <- required_output_files[
  !file.exists(
    required_output_files
  )
]


if (length(missing_output_files) > 0L) {

  stop(
    "Expected output file(s) missing:\n",
    paste(
      missing_output_files,
      collapse = "\n"
    )
  )
}


cat(
  "Final validation: PASS\n"
)


############################################################
## 29. FINAL CONSOLE SUMMARY
############################################################

cat(
  "\n============================================\n"
)

cat(
  "SOLID IMPORTANT-REGION REFINEMENT COMPLETE\n"
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
  "Tier-A regions:",
  format(
    nrow(
      tierA
    ),
    big.mark = ","
  ),
  "\n"
)


cat(
  "  Tumor higher:",
  format(
    nrow(
      tierA_tumor
    ),
    big.mark = ","
  ),
  "\n"
)


cat(
  "  Plasma higher:",
  format(
    nrow(
      tierA_plasma
    ),
    big.mark = ","
  ),
  "\n"
)


cat(
  "Tier-A regions in multi-window blocks:",
  format(
    nrow(
      tierA_block_supported
    ),
    big.mark = ","
  ),
  "\n"
)


cat(
  "Retained multi-window blocks:",
  format(
    nrow(
      blocks
    ),
    big.mark = ","
  ),
  "\n"
)


cat(
  "\nTier-A absolute Delta Beta summary:\n"
)


print(
  summary(
    tierA$abs_median_delta_beta
  )
)


cat(
  "\nTier-A direction consistency summary:\n"
)


print(
  summary(
    tierA$direction_consistency_pct
  )
)


cat(
  "\nTier-A effect consistency >=0.10 summary:\n"
)


print(
  summary(
    tierA$effect_consistency_0.10_pct
  )
)


cat(
  "\nChromosomes 14-22/X/Y:\n"
)


print(
  later_chr
)


cat(
  "\nResults saved to:\n",
  refinement_dir,
  "\n",
  sep = ""
)


############################################################
## 30. SESSION INFORMATION
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    refinement_dir,
    "sessionInfo_09_SOLID_important_region_refinement.txt"
  )
)