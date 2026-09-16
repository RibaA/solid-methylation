##------------------------------------------------------------
# SOLID cfMeDIP HARMONIZATION TO EXISTING 1-kb REGIONS
#
# Script:
#   11_SOLID_cfMeDIP_harmonization_1kb.r
#
# Goal:
#   Map the matched 13-patient cfMeDIP data from ~300-bp
#   regions to the SAME 1-kb genomic regions already used
#   in the SOLID tumor EPIC + plasma 5-base analysis.
#
# Input:
#   1. Matched 13-patient cfMeDIP log2CPM matrix
#   2. cfMeDIP feature annotation (~300-bp regions)
#   3. Existing SOLID tissue/plasma 1-kb region annotation
#
# Output:
#   - cfMeDIP 1-kb CPM-scale matrix
#   - cfMeDIP 1-kb log2CPM matrix
#   - harmonized 1-kb region annotation
#   - 300-bp -> 1-kb mapping table
#   - harmonization summary
#
# Important:
#   cfMeDIP input values are processed log2CPM enrichment values.
#   They are NOT methylation beta values.
#
# Aggregation:
#   1. Back-transform log2CPM to linear CPM-scale signal
#   2. Calculate an overlap-width-weighted mean across
#      cfMeDIP windows overlapping each existing 1-kb region
#   3. Transform the aggregated value back to log2 scale
#
# Note:
#   Because the exact upstream log2CPM transformation details
#   are not available here, 2^log2CPM is treated as a
#   back-transformed CPM-scale signal rather than assumed to
#   be the original raw CPM.
##------------------------------------------------------------


##------------------------------------------------------------
# 0. Setup
##------------------------------------------------------------

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
})


##------------------------------------------------------------
# 1. Directories
##------------------------------------------------------------

dir_cfmedip <- file.path(
  "result",
  "05_plasma_cfMeDIP"
)

dir_cfmedip_object <- file.path(
  dir_cfmedip,
  "objects"
)

dir_existing_matched <- file.path(
  "result",
  "03_matched_tissue_plasma"
)

dir_existing_preparation <- file.path(
  dir_existing_matched,
  "preparation"
)


##------------------------------------------------------------
# Output directories
##------------------------------------------------------------

dir_harmonized <- file.path(
  dir_cfmedip,
  "harmonized_1kb"
)

dir_harmonized_object <- file.path(
  dir_harmonized,
  "objects"
)

dir_harmonized_table <- file.path(
  dir_harmonized,
  "tables"
)

dir.create(
  dir_harmonized,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dir_harmonized_object,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dir_harmonized_table,
  recursive = TRUE,
  showWarnings = FALSE
)


##------------------------------------------------------------
# 2. Input files
##------------------------------------------------------------

file_medip_matched <- file.path(
  dir_cfmedip_object,
  "SOLID_cfMeDIP_matched13_log2CPM_matrix.rds"
)

file_medip_annotation <- file.path(
  dir_cfmedip,
  "tables",
  "cfMeDIP_feature_annotation.csv"
)

file_target_1kb <- file.path(
  dir_existing_preparation,
  "SOLID_tissue_1kb_region_annotation.tsv.gz"
)


##------------------------------------------------------------
# Check input files
##------------------------------------------------------------

input_files <- c(
  matched_cfMeDIP_matrix =
    file_medip_matched,

  cfMeDIP_feature_annotation =
    file_medip_annotation,

  existing_1kb_annotation =
    file_target_1kb
)

input_check <- data.frame(
  Input = names(
    input_files
  ),

  Path = unname(
    input_files
  ),

  Exists = file.exists(
    input_files
  ),

  stringsAsFactors = FALSE
)

cat("\n")
cat("====================================================\n")
cat("SCRIPT 11 INPUT FILE CHECK\n")
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
    paste0(
      "One or more required input files are missing.\n",
      "Do not continue until all three input paths are correct."
    )
  )
}


##------------------------------------------------------------
# 3. Load matched 13-patient cfMeDIP matrix
##------------------------------------------------------------

medip_matched <- readRDS(
  file_medip_matched
)

if (
  !is.matrix(
    medip_matched
  )
) {

  medip_matched <- as.matrix(
    medip_matched
  )
}

storage.mode(
  medip_matched
) <- "numeric"


cat("\n")
cat("====================================================\n")
cat("MATCHED cfMeDIP MATRIX\n")
cat("====================================================\n")

cat(
  "Regions:",
  format(
    nrow(
      medip_matched
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Patients:",
  ncol(
    medip_matched
  ),
  "\n"
)

cat(
  "\nPatient IDs:\n"
)

print(
  colnames(
    medip_matched
  )
)


##------------------------------------------------------------
# Safety check
##------------------------------------------------------------

if (
  ncol(
    medip_matched
  ) != 13
) {

  stop(
    paste0(
      "Expected 13 matched patients, but found ",
      ncol(
        medip_matched
      ),
      "."
    )
  )
}


##------------------------------------------------------------
# 4. Load cfMeDIP feature annotation
##------------------------------------------------------------

medip_annotation <- read.csv(
  file_medip_annotation,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

medip_annotation <- as.data.frame(
  medip_annotation
)


cat("\n")
cat("====================================================\n")
cat("cfMeDIP FEATURE ANNOTATION\n")
cat("====================================================\n")

cat(
  "Rows:",
  format(
    nrow(
      medip_annotation
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Columns:",
  paste(
    colnames(
      medip_annotation
    ),
    collapse = ", "
  ),
  "\n"
)


##------------------------------------------------------------
# Required annotation columns
##------------------------------------------------------------

required_medip_cols <- c(
  "chr",
  "start",
  "end",
  "region_id"
)

missing_medip_cols <- setdiff(
  required_medip_cols,
  colnames(
    medip_annotation
  )
)

if (
  length(
    missing_medip_cols
  ) > 0
) {

  stop(
    paste(
      "Missing cfMeDIP annotation columns:",
      paste(
        missing_medip_cols,
        collapse = ", "
      )
    )
  )
}


##------------------------------------------------------------
# Check matrix ↔ annotation correspondence
##------------------------------------------------------------

if (
  nrow(
    medip_annotation
  ) !=
    nrow(
      medip_matched
    )
) {

  stop(
    "cfMeDIP matrix and feature annotation have different numbers of regions."
  )
}


##------------------------------------------------------------
# Reorder annotation to match matrix if necessary
##------------------------------------------------------------

if (
  !identical(
    medip_annotation$region_id,
    rownames(
      medip_matched
    )
  )
) {

  medip_match <- match(
    rownames(
      medip_matched
    ),
    medip_annotation$region_id
  )

  if (
    any(
      is.na(
        medip_match
      )
    )
  ) {

    stop(
      "Could not match all cfMeDIP matrix regions to feature annotation."
    )
  }

  medip_annotation <- medip_annotation[
    medip_match,
    ,
    drop = FALSE
  ]
}


stopifnot(
  identical(
    medip_annotation$region_id,
    rownames(
      medip_matched
    )
  )
)


##------------------------------------------------------------
# 5. Load existing SOLID 1-kb annotation
##------------------------------------------------------------

target_1kb <- read.delim(
  file_target_1kb,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


cat("\n")
cat("====================================================\n")
cat("EXISTING SOLID 1-kb REGION ANNOTATION\n")
cat("====================================================\n")

cat(
  "Regions:",
  format(
    nrow(
      target_1kb
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "\nColumns:\n"
)

print(
  colnames(
    target_1kb
  )
)

cat(
  "\nFirst rows:\n"
)

print(
  head(
    target_1kb
  )
)


##------------------------------------------------------------
# 6. Identify genomic-coordinate columns
##------------------------------------------------------------

chr_candidates <- c(
  "chr",
  "chrom",
  "chromosome",
  "seqnames"
)

start_candidates <- c(
  "start",
  "Start"
)

end_candidates <- c(
  "end",
  "End"
)


chr_col <- chr_candidates[
  chr_candidates %in%
    colnames(
      target_1kb
    )
][1]

start_col <- start_candidates[
  start_candidates %in%
    colnames(
      target_1kb
    )
][1]

end_col <- end_candidates[
  end_candidates %in%
    colnames(
      target_1kb
    )
][1]


if (
  is.na(
    chr_col
  ) ||
    is.na(
      start_col
    ) ||
    is.na(
      end_col
    )
) {

  stop(
    paste0(
      "Could not identify chr/start/end columns in ",
      basename(
        file_target_1kb
      ),
      "."
    )
  )
}


cat(
  "\nUsing target-region coordinate columns:\n"
)

cat(
  "Chromosome:",
  chr_col,
  "\n"
)

cat(
  "Start:",
  start_col,
  "\n"
)

cat(
  "End:",
  end_col,
  "\n"
)


##------------------------------------------------------------
# Standardize target coordinates
##------------------------------------------------------------

target_1kb$chr_harmonized <- as.character(
  target_1kb[
    [
      chr_col
    ]
  ]
)

target_1kb$start_harmonized <- as.numeric(
  target_1kb[
    [
      start_col
    ]
  ]
)

target_1kb$end_harmonized <- as.numeric(
  target_1kb[
    [
      end_col
    ]
  ]
)


##------------------------------------------------------------
# Create target-region ID
##------------------------------------------------------------

target_1kb$region_1kb_id <- paste0(
  target_1kb$chr_harmonized,
  ":",
  target_1kb$start_harmonized,
  "-",
  target_1kb$end_harmonized
)


##------------------------------------------------------------
# Target-region widths
##------------------------------------------------------------

target_1kb$width_harmonized <- (
  target_1kb$end_harmonized -
    target_1kb$start_harmonized +
    1
)

cat(
  "\nMost common target-region widths:\n"
)

print(
  head(
    sort(
      table(
        target_1kb$width_harmonized
      ),
      decreasing = TRUE
    ),
    10
  )
)


##------------------------------------------------------------
# Remove invalid target coordinates
##------------------------------------------------------------

valid_target <- (
  !is.na(
    target_1kb$chr_harmonized
  ) &
    is.finite(
      target_1kb$start_harmonized
    ) &
    is.finite(
      target_1kb$end_harmonized
    ) &
    target_1kb$start_harmonized <=
      target_1kb$end_harmonized
)

if (
  !all(
    valid_target
  )
) {

  warning(
    paste(
      sum(
        !valid_target
      ),
      "target regions have invalid coordinates and will be removed."
    )
  )

  target_1kb <- target_1kb[
    valid_target,
    ,
    drop = FALSE
  ]
}


##------------------------------------------------------------
# 7. Harmonize chromosome naming
##------------------------------------------------------------

add_chr <- function(x) {

  x <- as.character(
    x
  )

  ifelse(
    grepl(
      "^chr",
      x,
      ignore.case = TRUE
    ),
    x,
    paste0(
      "chr",
      x
    )
  )
}


medip_annotation$chr_harmonized <- add_chr(
  medip_annotation$chr
)

target_1kb$chr_harmonized <- add_chr(
  target_1kb$chr_harmonized
)


##------------------------------------------------------------
# Keep standard chromosomes
##------------------------------------------------------------

standard_chr <- c(
  paste0(
    "chr",
    1:22
  ),
  "chrX",
  "chrY"
)

keep_medip_chr <- (
  medip_annotation$chr_harmonized %in%
    standard_chr
)

keep_target_chr <- (
  target_1kb$chr_harmonized %in%
    standard_chr
)


cat("\n")
cat("====================================================\n")
cat("STANDARD-CHROMOSOME FILTER\n")
cat("====================================================\n")

cat(
  "cfMeDIP regions retained:",
  sum(
    keep_medip_chr
  ),
  "/",
  length(
    keep_medip_chr
  ),
  "\n"
)

cat(
  "1-kb target regions retained:",
  sum(
    keep_target_chr
  ),
  "/",
  length(
    keep_target_chr
  ),
  "\n"
)


medip_annotation_use <- medip_annotation[
  keep_medip_chr,
  ,
  drop = FALSE
]

medip_use <- medip_matched[
  keep_medip_chr,
  ,
  drop = FALSE
]

target_1kb_use <- target_1kb[
  keep_target_chr,
  ,
  drop = FALSE
]


##------------------------------------------------------------
# 8. Build GRanges objects
##------------------------------------------------------------

gr_medip <- GRanges(
  seqnames =
    medip_annotation_use$chr_harmonized,

  ranges = IRanges(
    start =
      medip_annotation_use$start,

    end =
      medip_annotation_use$end
  )
)


gr_target <- GRanges(
  seqnames =
    target_1kb_use$chr_harmonized,

  ranges = IRanges(
    start =
      target_1kb_use$start_harmonized,

    end =
      target_1kb_use$end_harmonized
  )
)


##------------------------------------------------------------
# 9. Find overlaps
##------------------------------------------------------------

cat(
  "\nFinding cfMeDIP -> 1-kb overlaps...\n"
)

hits <- findOverlaps(
  query = gr_medip,
  subject = gr_target,
  ignore.strand = TRUE
)


if (
  length(
    hits
  ) == 0
) {

  stop(
    "No overlaps found between cfMeDIP regions and target 1-kb regions."
  )
}


q_hit <- queryHits(
  hits
)

s_hit <- subjectHits(
  hits
)


cat(
  "Total overlaps:",
  format(
    length(
      hits
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "cfMeDIP regions with >=1 target overlap:",
  format(
    length(
      unique(
        q_hit
      )
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Target 1-kb regions with >=1 cfMeDIP overlap:",
  format(
    length(
      unique(
        s_hit
      )
    ),
    big.mark = ","
  ),
  "\n"
)


##------------------------------------------------------------
# 10. Calculate bp overlap
##------------------------------------------------------------

overlap_start <- pmax(
  start(
    gr_medip
  )[
    q_hit
  ],

  start(
    gr_target
  )[
    s_hit
  ]
)

overlap_end <- pmin(
  end(
    gr_medip
  )[
    q_hit
  ],

  end(
    gr_target
  )[
    s_hit
  ]
)

overlap_bp <- (
  overlap_end -
    overlap_start +
    1
)


if (
  any(
    overlap_bp <= 0
  )
) {

  stop(
    "Invalid genomic overlap widths detected."
  )
}


##------------------------------------------------------------
# 11. Mapping summary per 1-kb region
##------------------------------------------------------------

n_windows_by_target <- tabulate(
  s_hit,
  nbins = length(
    gr_target
  )
)


overlap_bp_tmp <- rowsum(
  overlap_bp,
  group = s_hit,
  reorder = FALSE
)

overlap_bp_by_target <- numeric(
  length(
    gr_target
  )
)

overlap_bp_by_target[
  as.integer(
    rownames(
      overlap_bp_tmp
    )
  )
] <- overlap_bp_tmp[
  ,
  1
]


target_1kb_use$N_cfMeDIP_windows <- (
  n_windows_by_target
)

target_1kb_use$Total_cfMeDIP_overlap_bp <- (
  overlap_bp_by_target
)

target_1kb_use$Has_cfMeDIP <- (
  target_1kb_use$N_cfMeDIP_windows > 0
)


##------------------------------------------------------------
# Approximate fraction of 1-kb region supported by cfMeDIP
#
# pmin(1, ...) prevents values above 100% when overlapping
# cfMeDIP windows contribute duplicated genomic coverage.
##------------------------------------------------------------

target_1kb_use$cfMeDIP_coverage_fraction <- pmin(
  1,
  target_1kb_use$Total_cfMeDIP_overlap_bp /
    target_1kb_use$width_harmonized
)

target_1kb_use$cfMeDIP_coverage_percent <- (
  100 *
    target_1kb_use$cfMeDIP_coverage_fraction
)


cat("\n")
cat("====================================================\n")
cat("cfMeDIP -> 1-kb MAPPING SUMMARY\n")
cat("====================================================\n")

cat(
  "Target regions:",
  format(
    nrow(
      target_1kb_use
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Target regions with cfMeDIP support:",
  format(
    sum(
      target_1kb_use$Has_cfMeDIP
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Percent target regions with cfMeDIP support:",
  round(
    100 *
      mean(
        target_1kb_use$Has_cfMeDIP
      ),
    2
  ),
  "%\n"
)

cat(
  "\nNumber of cfMeDIP windows per supported 1-kb region:\n"
)

print(
  summary(
    target_1kb_use$N_cfMeDIP_windows[
      target_1kb_use$Has_cfMeDIP
    ]
  )
)

cat(
  "\nApproximate cfMeDIP coverage of supported 1-kb regions (%):\n"
)

print(
  summary(
    target_1kb_use$cfMeDIP_coverage_percent[
      target_1kb_use$Has_cfMeDIP
    ]
  )
)


##------------------------------------------------------------
# 12. Aggregate on linear CPM scale
#
# Step 1:
#   log2CPM -> linear CPM-scale signal
#
# Step 2:
#   overlap-width-weighted mean CPM-scale signal
#
# Step 3:
#   log2 transform aggregated signal
##------------------------------------------------------------


##------------------------------------------------------------
# Initialize linear-scale 1-kb matrix
##------------------------------------------------------------

medip_1kb_cpm <- matrix(
  NA_real_,
  nrow = length(
    gr_target
  ),
  ncol = ncol(
    medip_use
  )
)

rownames(
  medip_1kb_cpm
) <- target_1kb_use$region_1kb_id

colnames(
  medip_1kb_cpm
) <- colnames(
  medip_use
)


##------------------------------------------------------------
# Initialize log2-scale 1-kb matrix
##------------------------------------------------------------

medip_1kb_log2cpm <- matrix(
  NA_real_,
  nrow = length(
    gr_target
  ),
  ncol = ncol(
    medip_use
  )
)

rownames(
  medip_1kb_log2cpm
) <- target_1kb_use$region_1kb_id

colnames(
  medip_1kb_log2cpm
) <- colnames(
  medip_use
)


##------------------------------------------------------------
# Total overlap bp per target region
##------------------------------------------------------------

weight_sum <- rowsum(
  overlap_bp,
  group = s_hit,
  reorder = FALSE
)

weight_target_ids <- as.integer(
  rownames(
    weight_sum
  )
)


##------------------------------------------------------------
# Aggregate patient by patient
##------------------------------------------------------------

cat(
  "\nAggregating cfMeDIP on linear CPM scale to 1-kb regions...\n"
)


for (
  j in seq_len(
    ncol(
      medip_use
    )
  )
) {

  ##----------------------------------------------------------
  # Back-transform log2CPM
  ##----------------------------------------------------------

  cpm_linear <- 2^(
    medip_use[
      q_hit,
      j
    ]
  )


  ##----------------------------------------------------------
  # Weight by number of overlapping base pairs
  ##----------------------------------------------------------

  weighted_cpm <- (
    cpm_linear *
      overlap_bp
  )


  ##----------------------------------------------------------
  # Sum weighted signal within each target region
  ##----------------------------------------------------------

  weighted_sum <- rowsum(
    weighted_cpm,
    group = s_hit,
    reorder = FALSE
  )

  weighted_target_ids <- as.integer(
    rownames(
      weighted_sum
    )
  )


  ##----------------------------------------------------------
  # Safety check
  ##----------------------------------------------------------

  if (
    !identical(
      weighted_target_ids,
      weight_target_ids
    )
  ) {

    stop(
      "Internal overlap aggregation ordering mismatch."
    )
  }


  ##----------------------------------------------------------
  # Weighted mean CPM-scale signal
  ##----------------------------------------------------------

  mean_cpm <- (
    weighted_sum[
      ,
      1
    ] /
      weight_sum[
        ,
        1
      ]
  )


  ##----------------------------------------------------------
  # Save linear-scale matrix
  ##----------------------------------------------------------

  medip_1kb_cpm[
    weighted_target_ids,
    j
  ] <- mean_cpm


  ##----------------------------------------------------------
  # Return aggregated signal to log2 scale
  ##----------------------------------------------------------

  medip_1kb_log2cpm[
    weighted_target_ids,
    j
  ] <- log2(
    mean_cpm
  )
}


##------------------------------------------------------------
# 13. Keep 1-kb regions with cfMeDIP support
##------------------------------------------------------------

keep_1kb <- (
  target_1kb_use$Has_cfMeDIP
)

medip_1kb_cpm_mapped <- medip_1kb_cpm[
  keep_1kb,
  ,
  drop = FALSE
]

medip_1kb_mapped <- medip_1kb_log2cpm[
  keep_1kb,
  ,
  drop = FALSE
]

annotation_1kb_mapped <- target_1kb_use[
  keep_1kb,
  ,
  drop = FALSE
]


##------------------------------------------------------------
# Check row correspondence
##------------------------------------------------------------

stopifnot(
  identical(
    rownames(
      medip_1kb_cpm_mapped
    ),
    annotation_1kb_mapped$region_1kb_id
  )
)

stopifnot(
  identical(
    rownames(
      medip_1kb_mapped
    ),
    annotation_1kb_mapped$region_1kb_id
  )
)


cat("\n")
cat("====================================================\n")
cat("HARMONIZED cfMeDIP 1-kb MATRICES\n")
cat("====================================================\n")

cat(
  "Regions:",
  format(
    nrow(
      medip_1kb_mapped
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Patients:",
  ncol(
    medip_1kb_mapped
  ),
  "\n"
)

cat(
  "Missing linear-scale values:",
  format(
    sum(
      is.na(
        medip_1kb_cpm_mapped
      )
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Missing log2-scale values:",
  format(
    sum(
      is.na(
        medip_1kb_mapped
      )
    ),
    big.mark = ","
  ),
  "\n"
)


##------------------------------------------------------------
# 14. Harmonization summary
##------------------------------------------------------------

harmonization_summary <- data.frame(

  Metric = c(
    "Original cfMeDIP regions",
    "Matched patients",
    "Target 1-kb regions",
    "Target regions with cfMeDIP support",
    "Target regions without cfMeDIP support",
    "Percent target regions with cfMeDIP support",
    "Median supported-region cfMeDIP coverage percent",
    "Total cfMeDIP-target overlaps",
    "Final harmonized 1-kb regions",
    "Final harmonized patients"
  ),

  Value = c(

    nrow(
      medip_matched
    ),

    ncol(
      medip_matched
    ),

    nrow(
      target_1kb_use
    ),

    sum(
      target_1kb_use$Has_cfMeDIP
    ),

    sum(
      !target_1kb_use$Has_cfMeDIP
    ),

    round(
      100 *
        mean(
          target_1kb_use$Has_cfMeDIP
        ),
      4
    ),

    round(
      median(
        target_1kb_use$cfMeDIP_coverage_percent[
          target_1kb_use$Has_cfMeDIP
        ],
        na.rm = TRUE
      ),
      4
    ),

    length(
      hits
    ),

    nrow(
      medip_1kb_mapped
    ),

    ncol(
      medip_1kb_mapped
    )
  ),

  stringsAsFactors = FALSE
)


print(
  harmonization_summary
)


##------------------------------------------------------------
# 15. Save harmonized outputs
##------------------------------------------------------------

saveRDS(
  medip_1kb_cpm_mapped,
  file.path(
    dir_harmonized_object,
    "SOLID_cfMeDIP_matched13_1kb_CPMscale.rds"
  )
)


saveRDS(
  medip_1kb_mapped,
  file.path(
    dir_harmonized_object,
    "SOLID_cfMeDIP_matched13_1kb_log2CPM.rds"
  )
)


saveRDS(
  annotation_1kb_mapped,
  file.path(
    dir_harmonized_object,
    "SOLID_cfMeDIP_matched13_1kb_region_annotation.rds"
  )
)


write.table(
  annotation_1kb_mapped,
  file.path(
    dir_harmonized_table,
    "SOLID_cfMeDIP_matched13_1kb_region_annotation.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.csv(
  harmonization_summary,
  file.path(
    dir_harmonized_table,
    "SOLID_cfMeDIP_1kb_harmonization_summary.csv"
  ),
  row.names = FALSE
)


##------------------------------------------------------------
# 16. Save cfMeDIP -> 1-kb overlap mapping
##------------------------------------------------------------

overlap_mapping <- data.frame(

  cfMeDIP_region =
    medip_annotation_use$region_id[
      q_hit
    ],

  cfMeDIP_chr =
    medip_annotation_use$chr_harmonized[
      q_hit
    ],

  cfMeDIP_start =
    medip_annotation_use$start[
      q_hit
    ],

  cfMeDIP_end =
    medip_annotation_use$end[
      q_hit
    ],

  target_1kb_region =
    target_1kb_use$region_1kb_id[
      s_hit
    ],

  overlap_bp =
    overlap_bp,

  stringsAsFactors = FALSE
)


write.table(
  overlap_mapping,
  file.path(
    dir_harmonized_table,
    "SOLID_cfMeDIP_300bp_to_1kb_mapping.tsv.gz"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


##------------------------------------------------------------
# 17. Save Script 11 checkpoint
##------------------------------------------------------------

saveRDS(
  list(

    medip_1kb_CPMscale =
      medip_1kb_cpm_mapped,

    medip_1kb_log2CPM =
      medip_1kb_mapped,

    region_annotation_1kb =
      annotation_1kb_mapped,

    harmonization_summary =
      harmonization_summary,

    patient_ids =
      colnames(
        medip_1kb_mapped
      )

  ),

  file.path(
    dir_harmonized_object,
    "SOLID_cfMeDIP_1kb_harmonization_checkpoint.rds"
  )
)


##------------------------------------------------------------
# 18. Final console summary
##------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("SOLID cfMeDIP 1-kb HARMONIZATION COMPLETE\n")
cat("====================================================\n")

cat(
  "Input cfMeDIP regions:",
  format(
    nrow(
      medip_matched
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Input patients:",
  ncol(
    medip_matched
  ),
  "\n"
)

cat(
  "Existing target 1-kb regions:",
  format(
    nrow(
      target_1kb_use
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Final cfMeDIP-supported 1-kb regions:",
  format(
    nrow(
      medip_1kb_mapped
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Final patients:",
  ncol(
    medip_1kb_mapped
  ),
  "\n"
)

cat(
  "Median supported-region cfMeDIP coverage:",
  round(
    median(
      annotation_1kb_mapped$cfMeDIP_coverage_percent,
      na.rm = TRUE
    ),
    2
  ),
  "%\n"
)

cat(
  "\nSaved linear-scale matrix:\n",
  file.path(
    dir_harmonized_object,
    "SOLID_cfMeDIP_matched13_1kb_CPMscale.rds"
  ),
  "\n"
)

cat(
  "\nSaved log2-scale matrix:\n",
  file.path(
    dir_harmonized_object,
    "SOLID_cfMeDIP_matched13_1kb_log2CPM.rds"
  ),
  "\n"
)

cat(
  "\nOutput directory:\n",
  dir_harmonized,
  "\n"
)

sessionInfo()