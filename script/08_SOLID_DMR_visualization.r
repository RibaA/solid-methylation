############################################################
## 08_SOLID_DMR_finalization.r
##
## SOLID DMR FINALIZATION / DATA LOCK
##
## Purpose:
##   1. Load Script-07 paired DMR object
##   2. Recalculate patient-level consistency metrics
##   3. Define the final SOLID Tier-A region set
##   4. Build Tier-A multi-window blocks from Tier-A regions ONLY
##   5. Validate the expected frozen SOLID counts
##   6. Save canonical final tables/object for downstream integration
##
## Tier-A definition:
##   FDR < 0.05
##   |median Delta Beta| >= 0.30
##   direction consistency = 100%
##   same-direction |Delta Beta| >= 0.10 in >=80% of patients
##   valid matched pairs = 13/13
##
## Delta Beta = Plasma - Tumor
############################################################

rm(list = ls())

options(
  stringsAsFactors = FALSE,
  width = 180,
  scipen = 999,
  warn = 1
)

expected_matched_patients <- 13L
fdr_threshold <- 0.05
tierA_abs_delta_beta <- 0.30
tierA_direction_consistency <- 100
tierA_effect_consistency_010 <- 80
tierA_required_valid_pairs <- 13L
maximum_gap_between_windows <- 0L
minimum_windows_per_block <- 2L

expected_tierA_regions <- 15322L
expected_tierA_tumor_higher <- 11350L
expected_tierA_plasma_higher <- 3972L
expected_tierA_windows_in_multiwindow_blocks <- 4492L
expected_tierA_multiwindow_blocks <- 2920L

required_packages <- c(
  "data.table",
  "matrixStats"
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
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(matrixStats)
})

dmr_dir <- file.path(
  "result",
  "03_matched_tissue_plasma",
  "DMR"
)

final_dir <- file.path(
  dmr_dir,
  "final"
)

dir.create(
  final_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dmr_object_file <- file.path(
  dmr_dir,
  "SOLID_tissue_vs_plasma_DMR_object.rds"
)

if (!file.exists(dmr_object_file)) {
  stop(
    "Script-07 DMR object not found:\n",
    dmr_object_file,
    "\nRun Script 07 first."
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

standardize_chr_character <- function(x) {
  sub(
    "^chr",
    "",
    as.character(x),
    ignore.case = TRUE
  )
}

chromosome_order <- function(x) {
  match(
    standardize_chr_character(x),
    c(
      as.character(1:22),
      "X",
      "Y"
    )
  )
}

message_header(
  "LOADING SCRIPT-07 DMR OBJECT"
)

dmr_object <- readRDS(
  dmr_object_file
)

required_elements <- c(
  "all_results",
  "region_order",
  "delta_beta",
  "tissue_beta",
  "plasma_beta",
  "valid_pair_mask",
  "patient_metadata",
  "analysis_settings"
)

missing_elements <- setdiff(
  required_elements,
  names(dmr_object)
)

if (length(missing_elements) > 0L) {
  stop(
    "DMR object missing required element(s): ",
    paste(missing_elements, collapse = ", ")
  )
}

results <- as.data.table(
  dmr_object$all_results
)

canonical_region_order <- as.character(
  dmr_object$region_order
)

delta_beta <- dmr_object$delta_beta
tissue_beta <- dmr_object$tissue_beta
plasma_beta <- dmr_object$plasma_beta
valid_pair_mask <- dmr_object$valid_pair_mask

patient_metadata <- as.data.table(
  dmr_object$patient_metadata
)

message_header(
  "VALIDATING CANONICAL REGION ORDER"
)

stopifnot(
  nrow(results) == length(canonical_region_order),
  nrow(delta_beta) == length(canonical_region_order),
  ncol(delta_beta) == expected_matched_patients,

  identical(
    results$region_id,
    canonical_region_order
  ),

  identical(
    rownames(delta_beta),
    canonical_region_order
  ),

  identical(
    rownames(tissue_beta),
    canonical_region_order
  ),

  identical(
    rownames(plasma_beta),
    canonical_region_order
  ),

  identical(
    rownames(valid_pair_mask),
    canonical_region_order
  ),

  identical(
    colnames(delta_beta),
    colnames(tissue_beta)
  ),

  identical(
    colnames(delta_beta),
    colnames(plasma_beta)
  ),

  identical(
    colnames(delta_beta),
    colnames(valid_pair_mask)
  ),

  identical(
    colnames(delta_beta),
    patient_metadata$patient_id
  )
)

required_result_columns <- c(
  "region_id",
  "FDR",
  "median_delta_beta",
  "abs_median_delta_beta"
)

missing_result_columns <- setdiff(
  required_result_columns,
  names(results)
)

if (length(missing_result_columns) > 0L) {
  stop(
    "Script-07 results missing required column(s): ",
    paste(missing_result_columns, collapse = ", ")
  )
}

cat(
  "Canonical region-order validation: PASS\n"
)

cat(
  "Regions tested: ",
  format(nrow(results), big.mark = ","),
  "\n",
  sep = ""
)

message_header(
  "RECALCULATING CONSISTENCY METRICS"
)

valid_delta <-
  valid_pair_mask &
  is.finite(delta_beta)

n_valid_pairs_final <- rowSums(
  valid_delta
)

n_positive <- rowSums(
  (delta_beta > 0) &
    valid_delta,
  na.rm = TRUE
)

n_negative <- rowSums(
  (delta_beta < 0) &
    valid_delta,
  na.rm = TRUE
)

n_positive_010 <- rowSums(
  (delta_beta >= 0.10) &
    valid_delta,
  na.rm = TRUE
)

n_negative_010 <- rowSums(
  (delta_beta <= -0.10) &
    valid_delta,
  na.rm = TRUE
)

results[
  ,
  n_valid_pairs_final :=
    n_valid_pairs_final
]

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

results[
  ,
  direction_consistency_pct :=
    fifelse(
      median_delta_beta > 0,
      100 *
        n_positive /
        n_valid_pairs_final,
      fifelse(
        median_delta_beta < 0,
        100 *
          n_negative /
          n_valid_pairs_final,
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
        n_valid_pairs_final,
      fifelse(
        median_delta_beta < 0,
        100 *
          n_negative_010 /
          n_valid_pairs_final,
        NA_real_
      )
    )
]

if ("n_valid_pairs" %in% names(results)) {
  if (
    !identical(
      as.integer(results$n_valid_pairs),
      as.integer(results$n_valid_pairs_final)
    )
  ) {
    stop(
      "Recalculated n_valid_pairs does not match Script-07 n_valid_pairs."
    )
  }
}

message_header(
  "DEFINING FINAL SOLID TIER-A REGIONS"
)

results[
  ,
  TierA :=
    FDR < fdr_threshold &
    abs_median_delta_beta >=
      tierA_abs_delta_beta &
    direction_consistency_pct >=
      tierA_direction_consistency &
    effect_consistency_0.10_pct >=
      tierA_effect_consistency_010 &
    n_valid_pairs_final ==
      tierA_required_valid_pairs
]

tierA <- copy(
  results[
    TierA == TRUE
  ]
)

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

observed_tierA_regions <- nrow(
  tierA
)

observed_tierA_tumor <- nrow(
  tierA_tumor
)

observed_tierA_plasma <- nrow(
  tierA_plasma
)

cat(
  "Tier-A regions: ",
  format(observed_tierA_regions, big.mark = ","),
  "\n",
  sep = ""
)

cat(
  "  Tumor higher: ",
  format(observed_tierA_tumor, big.mark = ","),
  "\n",
  sep = ""
)

cat(
  "  Plasma higher: ",
  format(observed_tierA_plasma, big.mark = ","),
  "\n",
  sep = ""
)

if (observed_tierA_regions != expected_tierA_regions) {
  stop(
    "Tier-A region count mismatch. Expected ",
    expected_tierA_regions,
    "; observed ",
    observed_tierA_regions,
    "."
  )
}

if (observed_tierA_tumor != expected_tierA_tumor_higher) {
  stop(
    "Tier-A tumor-higher count mismatch. Expected ",
    expected_tierA_tumor_higher,
    "; observed ",
    observed_tierA_tumor,
    "."
  )
}

if (observed_tierA_plasma != expected_tierA_plasma_higher) {
  stop(
    "Tier-A plasma-higher count mismatch. Expected ",
    expected_tierA_plasma_higher,
    "; observed ",
    observed_tierA_plasma,
    "."
  )
}

cat(
  "Tier-A count validation: PASS\n"
)

message_header(
  "BUILDING TIER-A MULTI-WINDOW BLOCKS"
)

required_coordinate_columns <- c(
  "chr",
  "start",
  "end"
)

missing_coordinate_columns <- setdiff(
  required_coordinate_columns,
  names(tierA)
)

if (length(missing_coordinate_columns) > 0L) {
  stop(
    "Tier-A table missing genomic coordinate column(s): ",
    paste(missing_coordinate_columns, collapse = ", ")
  )
}

block_input <- tierA[
  !is.na(chr) &
    !is.na(start) &
    !is.na(end),
  .(
    region_id,
    chr =
      standardize_chr_character(chr),
    start =
      as.integer(start),
    end =
      as.integer(end),
    DMR_direction =
      as.character(DMR_direction),
    FDR,
    median_delta_beta,
    abs_median_delta_beta,
    direction_consistency_pct,
    effect_consistency_0.10_pct,
    n_valid_pairs_final
  )
]

block_input[
  ,
  chr_order :=
    chromosome_order(chr)
]

block_input <- block_input[
  !is.na(chr_order)
]

setorder(
  block_input,
  chr_order,
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
    is.na(previous_end) |
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

block_key <- unique(
  block_input[
    ,
    .(
      chr_order,
      chr,
      DMR_direction,
      block_number
    )
  ]
)

setorder(
  block_key,
  chr_order,
  DMR_direction,
  block_number
)

block_key[
  ,
  block_id :=
    paste0(
      "SOLID_TIERA_BLOCK_",
      sprintf(
        "%05d",
        seq_len(.N)
      )
    )
]

block_input <- merge(
  block_input,
  block_key,
  by = c(
    "chr_order",
    "chr",
    "DMR_direction",
    "block_number"
  ),
  all.x = TRUE,
  sort = FALSE
)

setorder(
  block_input,
  chr_order,
  DMR_direction,
  start,
  end
)

tierA_blocks <- block_input[
  ,
  .(
    block_start =
      min(start),

    block_end =
      max(end),

    block_width =
      max(end) -
      min(start),

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
        abs_median_delta_beta,
        na.rm = TRUE
      ),

    median_effect_consistency_0.10_pct =
      median(
        effect_consistency_0.10_pct,
        na.rm = TRUE
      ),

    first_region_id =
      first(region_id),

    last_region_id =
      last(region_id)
  ),
  by = .(
    block_id,
    chr_order,
    chr,
    DMR_direction
  )
]

tierA_blocks[
  ,
  retained_multiwindow :=
    n_windows >=
    minimum_windows_per_block
]

tierA_retained_blocks <- copy(
  tierA_blocks[
    retained_multiwindow == TRUE
  ]
)

setorder(
  tierA_blocks,
  chr_order,
  block_start,
  DMR_direction
)

setorder(
  tierA_retained_blocks,
  chr_order,
  block_start,
  DMR_direction
)

retained_block_ids <- tierA_retained_blocks$block_id

tierA_block_members <- block_input[
  block_id %in%
    retained_block_ids,
  .(
    block_id,
    region_id,
    chr,
    start,
    end,
    DMR_direction,
    median_delta_beta,
    FDR,
    direction_consistency_pct,
    effect_consistency_0.10_pct,
    n_valid_pairs_final
  )
]

setorder(
  tierA_block_members,
  chr,
  start,
  DMR_direction
)

observed_member_windows <- nrow(
  tierA_block_members
)

observed_retained_blocks <- nrow(
  tierA_retained_blocks
)

cat(
  "Tier-A windows in retained multi-window blocks: ",
  format(observed_member_windows, big.mark = ","),
  "\n",
  sep = ""
)

cat(
  "Tier-A retained multi-window blocks: ",
  format(observed_retained_blocks, big.mark = ","),
  "\n",
  sep = ""
)

if (
  observed_member_windows !=
    expected_tierA_windows_in_multiwindow_blocks
) {
  stop(
    "Tier-A block-member count mismatch. Expected ",
    expected_tierA_windows_in_multiwindow_blocks,
    "; observed ",
    observed_member_windows,
    "."
  )
}

if (
  observed_retained_blocks !=
    expected_tierA_multiwindow_blocks
) {
  stop(
    "Tier-A retained-block count mismatch. Expected ",
    expected_tierA_multiwindow_blocks,
    "; observed ",
    observed_retained_blocks,
    "."
  )
}

cat(
  "Tier-A block validation: PASS\n"
)

tierA[
  ,
  `:=`(
    in_multiwindow_block =
      region_id %in%
      tierA_block_members$region_id,

    block_id =
      NA_character_,

    block_n_windows =
      NA_integer_,

    block_width =
      NA_integer_
  )
]

if (nrow(tierA_block_members) > 0L) {

  member_map <- merge(
    tierA_block_members[
      ,
      .(
        block_id,
        region_id
      )
    ],
    tierA_retained_blocks[
      ,
      .(
        block_id,
        block_n_windows =
          n_windows,
        block_width
      )
    ],
    by = "block_id",
    all.x = TRUE,
    sort = FALSE
  )

  tierA[
    member_map,
    on = "region_id",
    `:=`(
      block_id =
        i.block_id,
      block_n_windows =
        i.block_n_windows,
      block_width =
        i.block_width
    )
  ]
}

tierA[
  ,
  block_supported_rank :=
    as.integer(
      in_multiwindow_block
    )
]

technical_columns <- intersect(
  c(
    "n_EPIC_probes",
    "median_plasma_coverage",
    "median_plasma_covered_cpgs"
  ),
  names(tierA)
)

tierA[
  ,
  EPIC_probe_rank :=
    if (
      "n_EPIC_probes" %in%
        names(tierA)
    ) {
      fifelse(
        is.na(n_EPIC_probes),
        0,
        pmin(
          n_EPIC_probes,
          10
        )
      )
    } else {
      0
    }
]

tierA[
  ,
  plasma_cpg_rank :=
    if (
      "median_plasma_covered_cpgs" %in%
        names(tierA)
    ) {
      fifelse(
        is.na(median_plasma_covered_cpgs),
        0,
        median_plasma_covered_cpgs
      )
    } else {
      0
    }
]

tierA[
  ,
  plasma_coverage_rank :=
    if (
      "median_plasma_coverage" %in%
        names(tierA)
    ) {
      fifelse(
        is.na(median_plasma_coverage),
        0,
        median_plasma_coverage
      )
    } else {
      0
    }
]

setorder(
  tierA,
  -block_supported_rank,
  -effect_consistency_0.10_pct,
  -abs_median_delta_beta,
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

final_summary <- data.table(
  item = c(
    "regions_tested",
    "FDR_significant_regions",
    "TierA_regions",
    "TierA_tumor_higher",
    "TierA_plasma_higher",
    "TierA_windows_in_multiwindow_blocks",
    "TierA_retained_multiwindow_blocks"
  ),

  value = c(
    nrow(results),

    sum(
      results$FDR <
        fdr_threshold,
      na.rm = TRUE
    ),

    nrow(tierA),

    nrow(tierA_tumor),

    nrow(tierA_plasma),

    nrow(tierA_block_members),

    nrow(tierA_retained_blocks)
  )
)

genome_build_value <- dmr_object$analysis_settings[
  setting == "genome_build",
  value
][1]

final_settings <- data.table(
  setting = c(
    "source_script",
    "source_DMR_object",
    "delta_beta_definition",
    "FDR_threshold",
    "TierA_abs_median_delta_beta_threshold",
    "TierA_direction_consistency_pct",
    "TierA_effect_consistency_0.10_pct",
    "TierA_required_valid_pairs",
    "coordinate_convention",
    "maximum_gap_between_windows",
    "minimum_windows_per_block",
    "block_input_population",
    "genome_build"
  ),

  value = c(
    "08_SOLID_DMR_finalization.r",
    dmr_object_file,
    "plasma_minus_tumor",
    fdr_threshold,
    tierA_abs_delta_beta,
    tierA_direction_consistency,
    tierA_effect_consistency_010,
    tierA_required_valid_pairs,
    "0_based_half_open",
    maximum_gap_between_windows,
    minimum_windows_per_block,
    "TierA_regions_only",
    as.character(genome_build_value)
  )
)

message_header(
  "SAVING FINAL SOLID DMR DATA LOCK"
)

fwrite(
  results,
  file.path(
    final_dir,
    "SOLID_DMR_all_regions_final.tsv.gz"
  ),
  sep = "\t",
  compress = "gzip"
)

fwrite(
  tierA,
  file.path(
    final_dir,
    "SOLID_TierA_regions.tsv.gz"
  ),
  sep = "\t",
  compress = "gzip"
)

fwrite(
  tierA_tumor,
  file.path(
    final_dir,
    "SOLID_TierA_tumor_higher.tsv.gz"
  ),
  sep = "\t",
  compress = "gzip"
)

fwrite(
  tierA_plasma,
  file.path(
    final_dir,
    "SOLID_TierA_plasma_higher.tsv.gz"
  ),
  sep = "\t",
  compress = "gzip"
)

fwrite(
  tierA_blocks,
  file.path(
    final_dir,
    "SOLID_TierA_all_blocks.tsv.gz"
  ),
  sep = "\t",
  compress = "gzip"
)

fwrite(
  tierA_retained_blocks,
  file.path(
    final_dir,
    "SOLID_TierA_retained_multiwindow_blocks.tsv"
  ),
  sep = "\t"
)

fwrite(
  tierA_block_members,
  file.path(
    final_dir,
    "SOLID_TierA_block_members.tsv.gz"
  ),
  sep = "\t",
  compress = "gzip"
)

fwrite(
  final_summary,
  file.path(
    final_dir,
    "SOLID_DMR_final_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  final_settings,
  file.path(
    final_dir,
    "SOLID_DMR_final_settings.tsv"
  ),
  sep = "\t"
)

final_object <- list(
  all_results =
    results,

  TierA_regions =
    tierA,

  TierA_tumor_higher =
    tierA_tumor,

  TierA_plasma_higher =
    tierA_plasma,

  TierA_all_blocks =
    tierA_blocks,

  TierA_retained_blocks =
    tierA_retained_blocks,

  TierA_block_members =
    tierA_block_members,

  region_order =
    canonical_region_order,

  delta_beta =
    delta_beta,

  tissue_beta =
    tissue_beta,

  plasma_beta =
    plasma_beta,

  valid_pair_mask =
    valid_pair_mask,

  patient_metadata =
    patient_metadata,

  final_summary =
    final_summary,

  final_settings =
    final_settings,

  provenance = list(
    source_DMR_object =
      dmr_object_file,

    source_script =
      "08_SOLID_DMR_finalization.r",

    source_statistical_script =
      "07_SOLID_paired_DMR.r",

    TierA_definition =
      paste0(
        "FDR<",
        fdr_threshold,
        "; abs_median_delta_beta>=",
        tierA_abs_delta_beta,
        "; direction_consistency_pct>=",
        tierA_direction_consistency,
        "; effect_consistency_0.10_pct>=",
        tierA_effect_consistency_010,
        "; n_valid_pairs==",
        tierA_required_valid_pairs
      ),

    block_definition =
      paste0(
        "TierA only; same chromosome; same direction; ",
        "0-based half-open windows; maximum gap=",
        maximum_gap_between_windows,
        "; minimum windows=",
        minimum_windows_per_block
      ),

    creation_date =
      as.character(
        Sys.Date()
      )
  )
)

final_object_file <- file.path(
  final_dir,
  "SOLID_DMR_final_object.rds"
)

saveRDS(
  final_object,
  final_object_file
)

validation_object <- readRDS(
  final_object_file
)

stopifnot(
  identical(
    validation_object$region_order,
    validation_object$all_results$region_id
  ),

  identical(
    validation_object$region_order,
    rownames(
      validation_object$delta_beta
    )
  ),

  nrow(
    validation_object$TierA_regions
  ) ==
    expected_tierA_regions,

  nrow(
    validation_object$TierA_tumor_higher
  ) ==
    expected_tierA_tumor_higher,

  nrow(
    validation_object$TierA_plasma_higher
  ) ==
    expected_tierA_plasma_higher,

  nrow(
    validation_object$TierA_block_members
  ) ==
    expected_tierA_windows_in_multiwindow_blocks,

  nrow(
    validation_object$TierA_retained_blocks
  ) ==
    expected_tierA_multiwindow_blocks
)

rm(
  validation_object
)

message_header(
  "SOLID DMR DATA LOCK PASSED"
)

print(
  final_summary
)

cat(
  "\nTechnical-support columns available in Tier-A table:\n"
)

if (length(technical_columns) > 0L) {
  print(
    technical_columns
  )
} else {
  cat(
    "None of the optional technical-support columns were present.\n"
  )
}

cat(
  "\nCanonical final object:\n",
  final_object_file,
  "\n",
  sep = ""
)

cat(
  "\nFinal output directory:\n",
  final_dir,
  "\n",
  sep = ""
)

capture.output(
  sessionInfo(),
  file = file.path(
    final_dir,
    "sessionInfo_08_SOLID_DMR_finalization.txt"
  )
)
