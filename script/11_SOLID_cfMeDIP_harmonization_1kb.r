##------------------------------------------------------------
# SOLID cfMeDIP HARMONIZATION TO MATCHED 1-kb REGIONS
#
# Script:
#   11_SOLID_cfMeDIP_harmonization_1kb.r
#
# Goal:
#   Map matched 13-patient cfMeDIP data (~300-bp regions)
#   to the SAME retained 1-kb regions used in the matched
#   SOLID tumor EPIC + plasma 5-base analysis.
#
# Input:
#   1. Matched 13-patient cfMeDIP log2CPM matrix
#   2. cfMeDIP feature annotation
#   3. FINAL matched/filtered SOLID 1-kb region annotation
#
# Aggregation:
#   log2CPM
#      -> 2^x CPM-scale signal
#      -> overlap-width-weighted mean
#      -> log2 aggregated signal
#
# Important:
#   cfMeDIP is an enrichment assay.
#   These values are NOT methylation Beta values.
##------------------------------------------------------------


##------------------------------------------------------------
# 0. Setup
##------------------------------------------------------------

rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  scipen = 999
)

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

dir_paired <- file.path(
  "result",
  "03_matched_tissue_plasma",
  "paired_filtered"
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

file_medip <- file.path(
  dir_cfmedip_object,
  "SOLID_cfMeDIP_matched13_log2CPM_matrix.rds"
)

file_medip_annotation <- file.path(
  dir_cfmedip,
  "tables",
  "cfMeDIP_feature_annotation.csv"
)

file_target_1kb <- file.path(
  dir_paired,
  "SOLID_paired_retained_region_annotation.tsv.gz"
)

##------------------------------------------------------------
# Input check
##------------------------------------------------------------

input_files <- c(
  matched_cfMeDIP_matrix =
    file_medip,

  cfMeDIP_feature_annotation =
    file_medip_annotation,

  matched_1kb_annotation =
    file_target_1kb
)

input_check <- data.frame(
  Input = names(input_files),
  Path = unname(input_files),
  Exists = file.exists(input_files),
  stringsAsFactors = FALSE
)

cat("\n")
cat("====================================================\n")
cat("SCRIPT 11 INPUT CHECK\n")
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
    "One or more required Script 11 inputs are missing."
  )
}

##------------------------------------------------------------
# 3. Load matched cfMeDIP matrix
##------------------------------------------------------------

medip <- readRDS(
  file_medip
)

medip <- as.matrix(
  medip
)

storage.mode(
  medip
) <- "numeric"


cat("\n")
cat("====================================================\n")
cat("MATCHED cfMeDIP DATA\n")
cat("====================================================\n")

cat(
  "Regions:",
  format(
    nrow(medip),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Patients:",
  ncol(medip),
  "\n"
)

print(
  colnames(medip)
)


##------------------------------------------------------------
# Critical input checks
##------------------------------------------------------------

if (
  nrow(medip) == 0
) {

  stop(
    "Input matched cfMeDIP matrix has zero regions."
  )
}

if (
  ncol(medip) != 13
) {

  stop(
    paste0(
      "Expected 13 matched cfMeDIP patients; found ",
      ncol(medip),
      "."
    )
  )
}

if (
  is.null(
    rownames(medip)
  )
) {

  stop(
    "Original matched cfMeDIP matrix must have ~300-bp region IDs."
  )
}


##------------------------------------------------------------
# 4. Load cfMeDIP ~300-bp feature annotation
##------------------------------------------------------------

medip_annotation <- read.csv(
  file_medip_annotation,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_medip_cols <- c(
  "region_id",
  "chr",
  "start",
  "end"
)

missing_medip_cols <- setdiff(
  required_medip_cols,
  colnames(medip_annotation)
)

if (
  length(missing_medip_cols) > 0
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
# Matrix/annotation correspondence
##------------------------------------------------------------

if (
  nrow(medip_annotation) !=
    nrow(medip)
) {

  stop(
    "cfMeDIP matrix and feature annotation have different row counts."
  )
}


##------------------------------------------------------------
# Put annotation in exact matrix order
##------------------------------------------------------------

idx_medip <- match(
  rownames(medip),
  medip_annotation$region_id
)

if (
  anyNA(idx_medip)
) {

  stop(
    "Some cfMeDIP matrix regions are absent from feature annotation."
  )
}

medip_annotation <- medip_annotation[
  idx_medip,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    as.character(
      medip_annotation$region_id
    ),
    rownames(medip)
  )
)


##------------------------------------------------------------
# 5. Load FINAL matched 1-kb region annotation
##------------------------------------------------------------

target_1kb <- read.delim(
  file_target_1kb,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat("\n")
cat("====================================================\n")
cat("MATCHED/FILTERED 1-kb TARGET REGIONS\n")
cat("====================================================\n")

cat(
  "Regions:",
  format(
    nrow(target_1kb),
    big.mark = ","
  ),
  "\n"
)

print(
  colnames(target_1kb)
)


required_target_cols <- c(
  "region_id",
  "chr",
  "start",
  "end"
)

missing_target_cols <- setdiff(
  required_target_cols,
  colnames(target_1kb)
)

if (
  length(missing_target_cols) > 0
) {

  stop(
    paste(
      "Missing target annotation columns:",
      paste(
        missing_target_cols,
        collapse = ", "
      )
    )
  )
}


if (
  nrow(target_1kb) == 0
) {

  stop(
    "Matched 1-kb target annotation contains zero regions."
  )
}


##------------------------------------------------------------
# 6. Standardize chromosomes and coordinates
##------------------------------------------------------------

normalize_chr <- function(x) {

  x <- as.character(x)

  x <- sub(
    "^chr",
    "",
    x,
    ignore.case = TRUE
  )

  paste0(
    "chr",
    x
  )
}


medip_annotation$chr_std <- normalize_chr(
  medip_annotation$chr
)

target_1kb$chr_std <- normalize_chr(
  target_1kb$chr
)


medip_annotation$start <- as.integer(
  medip_annotation$start
)

medip_annotation$end <- as.integer(
  medip_annotation$end
)

target_1kb$start <- as.integer(
  target_1kb$start
)

target_1kb$end <- as.integer(
  target_1kb$end
)


##------------------------------------------------------------
# Standard chromosomes
##------------------------------------------------------------

standard_chr <- c(
  paste0(
    "chr",
    1:22
  ),
  "chrX",
  "chrY"
)


keep_medip <- (
  medip_annotation$chr_std %in%
    standard_chr
)

keep_target <- (
  target_1kb$chr_std %in%
    standard_chr
)


medip_annotation <- medip_annotation[
  keep_medip,
  ,
  drop = FALSE
]

medip <- medip[
  keep_medip,
  ,
  drop = FALSE
]

target_1kb <- target_1kb[
  keep_target,
  ,
  drop = FALSE
]


cat("\n")
cat("Standard chromosome cfMeDIP regions:",
    format(nrow(medip), big.mark = ","),
    "\n")

cat("Standard chromosome target regions:",
    format(nrow(target_1kb), big.mark = ","),
    "\n")


##------------------------------------------------------------
# Critical post-filter checks
##------------------------------------------------------------

if (
  nrow(medip) == 0
) {

  stop(
    "No cfMeDIP regions remain after chromosome filtering."
  )
}

if (
  nrow(target_1kb) == 0
) {

  stop(
    "No matched 1-kb target regions remain after chromosome filtering."
  )
}


##------------------------------------------------------------
# 7. Build GRanges
##------------------------------------------------------------

gr_medip <- GRanges(
  seqnames =
    medip_annotation$chr_std,

  ranges = IRanges(
    start =
      medip_annotation$start,

    end =
      medip_annotation$end
  )
)


gr_target <- GRanges(
  seqnames =
    target_1kb$chr_std,

  ranges = IRanges(
    start =
      target_1kb$start,

    end =
      target_1kb$end
  )
)


##------------------------------------------------------------
# 8. Find cfMeDIP -> matched 1-kb overlaps
##------------------------------------------------------------

hits <- findOverlaps(
  gr_medip,
  gr_target,
  ignore.strand = TRUE
)


if (
  length(hits) == 0
) {

  stop(
    "No cfMeDIP overlaps found with matched 1-kb regions."
  )
}


q_hit <- queryHits(
  hits
)

s_hit <- subjectHits(
  hits
)


##------------------------------------------------------------
# Exact matched target regions with cfMeDIP support
##------------------------------------------------------------

supported_target_idx <- sort(
  unique(
    s_hit
  )
)

n_supported <- length(
  supported_target_idx
)


cat("\n")
cat("====================================================\n")
cat("OVERLAP SUMMARY\n")
cat("====================================================\n")

cat(
  "Total cfMeDIP -> target overlaps:",
  format(
    length(hits),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Unique cfMeDIP regions overlapping targets:",
  format(
    length(unique(q_hit)),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Matched 1-kb regions with cfMeDIP support:",
  format(
    n_supported,
    big.mark = ","
  ),
  "\n"
)


if (
  n_supported == 0
) {

  stop(
    "No matched target regions have cfMeDIP support."
  )
}


##------------------------------------------------------------
# 9. Calculate overlap width
##------------------------------------------------------------

overlap_start <- pmax(
  start(gr_medip)[q_hit],
  start(gr_target)[s_hit]
)

overlap_end <- pmin(
  end(gr_medip)[q_hit],
  end(gr_target)[s_hit]
)

overlap_bp <- (
  overlap_end -
    overlap_start +
    1L
)


if (
  any(
    overlap_bp <= 0
  )
) {

  stop(
    "Invalid overlap widths detected."
  )
}


##------------------------------------------------------------
# 10. Convert target indices into FINAL output row indices
#
# supported_target_idx:
#   original row index in target_1kb
#
# output_target_idx:
#   row index in final supported matrix
##------------------------------------------------------------

output_target_idx <- match(
  s_hit,
  supported_target_idx
)


if (
  anyNA(
    output_target_idx
  )
) {

  stop(
    "Internal target-index mapping failed."
  )
}


##------------------------------------------------------------
# 11. Final supported target annotation
#
# IMPORTANT:
# This annotation uses EXACTLY the same row indices
# that define the output matrices.
##------------------------------------------------------------

annotation_1kb_mapped <- target_1kb[
  supported_target_idx,
  ,
  drop = FALSE
]


##------------------------------------------------------------
# Create canonical final region IDs
#
# Use the EXISTING matched 5-base convention:
#
#   chr:start:end
##------------------------------------------------------------

annotation_1kb_mapped$region_id <- paste(
  annotation_1kb_mapped$chr_std,
  annotation_1kb_mapped$start,
  annotation_1kb_mapped$end,
  sep = ":"
)


if (
  anyDuplicated(
    annotation_1kb_mapped$region_id
  )
) {

  stop(
    "Duplicated final 1-kb region IDs detected."
  )
}


##------------------------------------------------------------
# 12. Mapping QC per supported target
##------------------------------------------------------------

n_windows <- tabulate(
  output_target_idx,
  nbins = n_supported
)


total_overlap_bp <- rowsum(
  overlap_bp,
  group = output_target_idx,
  reorder = TRUE
)[
  ,
  1
]


annotation_1kb_mapped$N_cfMeDIP_windows <- (
  n_windows
)

annotation_1kb_mapped$Total_cfMeDIP_overlap_bp <- (
  total_overlap_bp
)


target_width <- (
  annotation_1kb_mapped$end -
    annotation_1kb_mapped$start +
    1L
)


annotation_1kb_mapped$cfMeDIP_coverage_fraction <- pmin(
  1,
  annotation_1kb_mapped$Total_cfMeDIP_overlap_bp /
    target_width
)

annotation_1kb_mapped$cfMeDIP_coverage_percent <- (
  100 *
    annotation_1kb_mapped$cfMeDIP_coverage_fraction
)


cat(
  "\ncfMeDIP windows per supported 1-kb region:\n"
)

print(
  summary(
    annotation_1kb_mapped$N_cfMeDIP_windows
  )
)

cat(
  "\nApproximate cfMeDIP region coverage (%):\n"
)

print(
  summary(
    annotation_1kb_mapped$cfMeDIP_coverage_percent
  )
)


##------------------------------------------------------------
# 13. Aggregate on LINEAR CPM scale
##------------------------------------------------------------

medip_1kb_cpm <- matrix(
  NA_real_,
  nrow = n_supported,
  ncol = ncol(medip),
  dimnames = list(
    annotation_1kb_mapped$region_id,
    colnames(medip)
  )
)


medip_1kb_log2cpm <- matrix(
  NA_real_,
  nrow = n_supported,
  ncol = ncol(medip),
  dimnames = list(
    annotation_1kb_mapped$region_id,
    colnames(medip)
  )
)


##------------------------------------------------------------
# Denominator for weighted mean
##------------------------------------------------------------

weight_sum <- rowsum(
  overlap_bp,
  group = output_target_idx,
  reorder = TRUE
)[
  ,
  1
]


if (
  length(weight_sum) != n_supported
) {

  stop(
    "Weight vector length does not match supported target count."
  )
}


##------------------------------------------------------------
# Aggregate each patient
##------------------------------------------------------------

cat(
  "\nAggregating cfMeDIP signal to matched 1-kb regions...\n"
)


for (
  j in seq_len(
    ncol(medip)
  )
) {

  ## Back-transform processed log2CPM
  ## to linear CPM-scale signal

  cpm_linear <- 2^(
    medip[
      q_hit,
      j
    ]
  )


  ## Weighted by overlap in bp

  weighted_signal <- (
    cpm_linear *
      overlap_bp
  )


  weighted_sum <- rowsum(
    weighted_signal,
    group = output_target_idx,
    reorder = TRUE
  )[
    ,
    1
  ]


  if (
    length(weighted_sum) != n_supported
  ) {

    stop(
      paste(
        "Aggregation failed for patient",
        colnames(medip)[j]
      )
    )
  }


  mean_cpm_scale <- (
    weighted_sum /
      weight_sum
  )


  medip_1kb_cpm[
    ,
    j
  ] <- mean_cpm_scale


  medip_1kb_log2cpm[
    ,
    j
  ] <- log2(
    mean_cpm_scale
  )
}


##------------------------------------------------------------
# 14. CRITICAL final validation
##------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("FINAL HARMONIZED OUTPUT VALIDATION\n")
cat("====================================================\n")


cat(
  "Annotation:",
  nrow(annotation_1kb_mapped),
  "rows\n"
)

cat(
  "CPM-scale matrix:",
  nrow(medip_1kb_cpm),
  "x",
  ncol(medip_1kb_cpm),
  "\n"
)

cat(
  "log2CPM matrix:",
  nrow(medip_1kb_log2cpm),
  "x",
  ncol(medip_1kb_log2cpm),
  "\n"
)


## ZERO rows are explicitly forbidden

if (
  nrow(
    annotation_1kb_mapped
  ) == 0
) {

  stop(
    "FINAL annotation has zero rows. Nothing will be saved."
  )
}


if (
  nrow(
    medip_1kb_cpm
  ) == 0
) {

  stop(
    "FINAL CPM matrix has zero rows. Nothing will be saved."
  )
}


if (
  nrow(
    medip_1kb_log2cpm
  ) == 0
) {

  stop(
    "FINAL log2CPM matrix has zero rows. Nothing will be saved."
  )
}


stopifnot(

  nrow(
    medip_1kb_cpm
  ) ==
    nrow(
      annotation_1kb_mapped
    ),

  nrow(
    medip_1kb_log2cpm
  ) ==
    nrow(
      annotation_1kb_mapped
    ),

  ncol(
    medip_1kb_log2cpm
  ) == 13,

  identical(
    rownames(
      medip_1kb_cpm
    ),
    annotation_1kb_mapped$region_id
  ),

  identical(
    rownames(
      medip_1kb_log2cpm
    ),
    annotation_1kb_mapped$region_id
  ),

  identical(
    colnames(
      medip_1kb_cpm
    ),
    colnames(
      medip
    )
  ),

  identical(
    colnames(
      medip_1kb_log2cpm
    ),
    colnames(
      medip
    )
  )
)


##------------------------------------------------------------
# Check non-finite values
##------------------------------------------------------------

cat(
  "Non-finite CPM-scale values:",
  sum(
    !is.finite(
      medip_1kb_cpm
    )
  ),
  "\n"
)

cat(
  "Non-finite log2CPM values:",
  sum(
    !is.finite(
      medip_1kb_log2cpm
    )
  ),
  "\n"
)


if (
  any(
    !is.finite(
      medip_1kb_cpm
    )
  )
) {

  stop(
    "Non-finite CPM-scale values detected."
  )
}


if (
  any(
    !is.finite(
      medip_1kb_log2cpm
    )
  )
) {

  stop(
    "Non-finite log2CPM values detected."
  )
}


##------------------------------------------------------------
# 15. Harmonization summary
##------------------------------------------------------------

harmonization_summary <- data.frame(

  Metric = c(
    "Input cfMeDIP regions",
    "Matched patients",
    "Matched 1-kb target regions",
    "cfMeDIP-supported matched 1-kb regions",
    "Percent matched regions with cfMeDIP support",
    "Total cfMeDIP-target overlaps",
    "Median cfMeDIP windows per supported region",
    "Median approximate cfMeDIP coverage percent"
  ),

  Value = c(

    nrow(
      medip
    ),

    ncol(
      medip
    ),

    nrow(
      target_1kb
    ),

    n_supported,

    round(
      100 *
        n_supported /
        nrow(target_1kb),
      4
    ),

    length(
      hits
    ),

    median(
      annotation_1kb_mapped$N_cfMeDIP_windows
    ),

    round(
      median(
        annotation_1kb_mapped$cfMeDIP_coverage_percent,
        na.rm = TRUE
      ),
      4
    )
  ),

  stringsAsFactors = FALSE
)


print(
  harmonization_summary
)


##------------------------------------------------------------
# 16. Save FINAL matched 1-kb objects
##------------------------------------------------------------

saveRDS(
  medip_1kb_cpm,
  file.path(
    dir_harmonized_object,
    "SOLID_cfMeDIP_matched13_1kb_CPMscale.rds"
  )
)


saveRDS(
  medip_1kb_log2cpm,
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
# 17. Save 300-bp -> matched 1-kb mapping
##------------------------------------------------------------

overlap_mapping <- data.frame(

  cfMeDIP_region =
    medip_annotation$region_id[
      q_hit
    ],

  cfMeDIP_chr =
    medip_annotation$chr_std[
      q_hit
    ],

  cfMeDIP_start =
    medip_annotation$start[
      q_hit
    ],

  cfMeDIP_end =
    medip_annotation$end[
      q_hit
    ],

  target_region_id =
    annotation_1kb_mapped$region_id[
      output_target_idx
    ],

  overlap_bp =
    overlap_bp,

  stringsAsFactors = FALSE
)


write.table(
  overlap_mapping,
  file.path(
    dir_harmonized_table,
    "SOLID_cfMeDIP_300bp_to_matched_1kb_mapping.tsv.gz"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


##------------------------------------------------------------
# 18. Save checkpoint
##------------------------------------------------------------

saveRDS(
  list(

    cfMeDIP_CPMscale =
      medip_1kb_cpm,

    cfMeDIP_log2CPM =
      medip_1kb_log2cpm,

    region_annotation =
      annotation_1kb_mapped,

    patient_ids =
      colnames(
        medip_1kb_log2cpm
      ),

    harmonization_summary =
      harmonization_summary

  ),

  file.path(
    dir_harmonized_object,
    "SOLID_cfMeDIP_matched13_1kb_harmonization_checkpoint.rds"
  )
)


##------------------------------------------------------------
# 19. Reload saved matrices and verify
##------------------------------------------------------------

check_log2 <- readRDS(
  file.path(
    dir_harmonized_object,
    "SOLID_cfMeDIP_matched13_1kb_log2CPM.rds"
  )
)

check_cpm <- readRDS(
  file.path(
    dir_harmonized_object,
    "SOLID_cfMeDIP_matched13_1kb_CPMscale.rds"
  )
)

check_annotation <- readRDS(
  file.path(
    dir_harmonized_object,
    "SOLID_cfMeDIP_matched13_1kb_region_annotation.rds"
  )
)


stopifnot(

  nrow(
    check_log2
  ) > 0,

  nrow(
    check_cpm
  ) > 0,

  nrow(
    check_annotation
  ) > 0,

  identical(
    dim(
      check_log2
    ),
    dim(
      medip_1kb_log2cpm
    )
  ),

  identical(
    dim(
      check_cpm
    ),
    dim(
      medip_1kb_cpm
    )
  ),

  identical(
    rownames(
      check_log2
    ),
    check_annotation$region_id
  )
)


##------------------------------------------------------------
# 20. Final summary
##------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("SCRIPT 11 COMPLETE\n")
cat("====================================================\n")

cat(
  "Matched patients:",
  ncol(
    check_log2
  ),
  "\n"
)

cat(
  "Matched 1-kb regions with cfMeDIP support:",
  format(
    nrow(
      check_log2
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Saved matrix dimensions:",
  nrow(
    check_log2
  ),
  "x",
  ncol(
    check_log2
  ),
  "\n"
)

cat(
  "Saved annotation rows:",
  nrow(
    check_annotation
  ),
  "\n"
)

cat(
  "\nOutput directory:\n",
  dir_harmonized,
  "\n"
)

sessionInfo()
