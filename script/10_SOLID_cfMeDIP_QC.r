##------------------------------------------------------------
# SOLID cfMeDIP PROCESSED-DATA QC
#
# Main goal:
#   QC cfMeDIP data and create the cfMeDIP subset corresponding
#   to the existing matched tumor EPIC + plasma 5-base cohort.
#
# Main analysis cohort:
#
#     Tumor EPIC
#         ∩
#     Plasma 5-base
#         ∩
#     Plasma cfMeDIP
#
# Expected final cohort:
#   13 matched SOLID patients
#
# cfMeDIP:
#   rows    = genomic regions (~300 bp)
#   columns = cfMeDIP samples
#   values  = processed log2CPM
#
# NOTE:
#   This script does NOT yet map cfMeDIP 300-bp regions to
#   the common 1-kb EPIC/5-base regions.
##------------------------------------------------------------


##------------------------------------------------------------
# 0. Setup
##------------------------------------------------------------

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(readxl)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})


##------------------------------------------------------------
# 1. Directories
##------------------------------------------------------------

dir_data <- "data"

dir_result <- file.path(
  "result",
  "05_plasma_cfMeDIP"
)

dir_fig <- file.path(
  dir_result,
  "figures"
)

dir_table <- file.path(
  dir_result,
  "tables"
)

dir_object <- file.path(
  dir_result,
  "objects"
)

dir.create(
  dir_result,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dir_fig,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dir_table,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dir_object,
  recursive = TRUE,
  showWarnings = FALSE
)


##------------------------------------------------------------
# 2. Input files
##------------------------------------------------------------

file_medip <- file.path(
  dir_data,
  "MeDIP_Solid_all_log2CPM.rds"
)

file_medip_map <- file.path(
  dir_data,
  "solid_cfDNA_leftover_with_medip_FINAL.xlsx"
)

file_research_meta <- file.path(
  dir_data,
  "meta_data_research.csv"
)

file_matched_meta <- file.path(
  "result",
  "03_matched_tissue_plasma",
  "preparation",
  "SOLID_tissue_harmonized_metadata.tsv"
)


##------------------------------------------------------------
# Confirm input files
##------------------------------------------------------------

input_files <- c(
  cfMeDIP_matrix = file_medip,
  cfMeDIP_mapping = file_medip_map,
  research_metadata = file_research_meta,
  matched_EPIC_metadata = file_matched_meta
)

cat("\n")
cat("====================================================\n")
cat("Input-file check\n")
cat("====================================================\n")

print(
  data.frame(
    File = names(input_files),
    Path = unname(input_files),
    Exists = file.exists(input_files),
    stringsAsFactors = FALSE
  )
)

if (
  !all(
    file.exists(
      input_files
    )
  )
) {

  stop(
    "One or more required input files are missing."
  )
}


##------------------------------------------------------------
# 3. Load cfMeDIP data
##------------------------------------------------------------

medip_raw <- readRDS(
  file_medip
)

cat("\n")
cat("====================================================\n")
cat("Original cfMeDIP object\n")
cat("====================================================\n")

cat(
  "Class:",
  paste(
    class(medip_raw),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Dimensions:",
  paste(
    dim(medip_raw),
    collapse = " x "
  ),
  "\n"
)


##------------------------------------------------------------
# 4. Convert cfMeDIP to numeric matrix
##------------------------------------------------------------

if (
  is.matrix(
    medip_raw
  )
) {

  medip <- medip_raw

} else if (
  is.data.frame(
    medip_raw
  )
) {

  numeric_cols <- vapply(
    medip_raw,
    is.numeric,
    logical(1)
  )

  if (
    sum(
      numeric_cols
    ) == 0
  ) {

    stop(
      "No numeric columns found in cfMeDIP object."
    )
  }

  medip <- as.matrix(
    medip_raw[
      ,
      numeric_cols,
      drop = FALSE
    ]
  )

} else {

  stop(
    paste(
      "Unexpected cfMeDIP object class:",
      paste(
        class(medip_raw),
        collapse = ", "
      )
    )
  )
}

storage.mode(
  medip
) <- "numeric"

n_features <- nrow(
  medip
)

n_samples <- ncol(
  medip
)

cat(
  "\ncfMeDIP regions:",
  format(
    n_features,
    big.mark = ","
  ),
  "\n"
)

cat(
  "cfMeDIP samples:",
  n_samples,
  "\n"
)


##------------------------------------------------------------
# 5. Parse cfMeDIP genomic coordinates
##------------------------------------------------------------

feature_ids <- rownames(
  medip
)

if (
  is.null(
    feature_ids
  )
) {

  stop(
    "cfMeDIP matrix does not have genomic region row names."
  )
}

feature_split <- strsplit(
  feature_ids,
  "\\."
)

valid_feature_format <- vapply(
  feature_split,
  length,
  integer(1)
) == 3

cat(
  "\nRegions matching chr.start.end format:",
  sum(
    valid_feature_format
  ),
  "/",
  length(
    feature_ids
  ),
  "\n"
)

feature_annotation <- do.call(
  rbind,
  lapply(
    feature_split,
    function(x) {

      if (
        length(
          x
        ) == 3
      ) {

        data.frame(
          chr = x[1],
          start = suppressWarnings(
            as.numeric(
              x[2]
            )
          ),
          end = suppressWarnings(
            as.numeric(
              x[3]
            )
          ),
          stringsAsFactors = FALSE
        )

      } else {

        data.frame(
          chr = NA_character_,
          start = NA_real_,
          end = NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  )
)

feature_annotation$region_id <- feature_ids

feature_annotation$width <- (
  feature_annotation$end -
    feature_annotation$start +
    1
)

cat(
  "\nMost common region widths:\n"
)

print(
  head(
    sort(
      table(
        feature_annotation$width
      ),
      decreasing = TRUE
    ),
    10
  )
)

write.csv(
  feature_annotation,
  file.path(
    dir_table,
    "cfMeDIP_feature_annotation.csv"
  ),
  row.names = FALSE
)


##------------------------------------------------------------
# 6. Load cfMeDIP mapping spreadsheet
##------------------------------------------------------------

medip_map <- read_excel(
  file_medip_map
)

medip_map <- as.data.frame(
  medip_map
)

cat("\n")
cat("====================================================\n")
cat("cfMeDIP mapping spreadsheet\n")
cat("====================================================\n")

cat(
  "Dimensions:",
  paste(
    dim(
      medip_map
    ),
    collapse = " x "
  ),
  "\n"
)

print(
  colnames(
    medip_map
  )
)


##------------------------------------------------------------
# 7. Load research metadata
##------------------------------------------------------------

research_meta <- read.csv(
  file_research_meta,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n")
cat("====================================================\n")
cat("SOLID research metadata\n")
cat("====================================================\n")

cat(
  "Dimensions:",
  paste(
    dim(
      research_meta
    ),
    collapse = " x "
  ),
  "\n"
)


##------------------------------------------------------------
# 8. Load matched tumor EPIC metadata
##------------------------------------------------------------

matched_meta <- read.delim(
  file_matched_meta,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n")
cat("====================================================\n")
cat("Matched tumor EPIC metadata\n")
cat("====================================================\n")

cat(
  "Dimensions:",
  paste(
    dim(
      matched_meta
    ),
    collapse = " x "
  ),
  "\n"
)

cat(
  "\nColumns:\n"
)

print(
  colnames(
    matched_meta
  )
)


##------------------------------------------------------------
# Confirm required tumor metadata variables
##------------------------------------------------------------

required_tumor_cols <- c(
  "Sample_Name",
  "Subject",
  "Array",
  "Genome_build"
)

missing_tumor_cols <- setdiff(
  required_tumor_cols,
  colnames(
    matched_meta
  )
)

if (
  length(
    missing_tumor_cols
  ) > 0
) {

  stop(
    paste(
      "Missing required columns from matched tumor metadata:",
      paste(
        missing_tumor_cols,
        collapse = ", "
      )
    )
  )
}

##------------------------------------------------------------
# 9. Basic QC of full cfMeDIP matrix
#
# Global integrity checks only.
# Biological analyses below use the final matched cohort.
##------------------------------------------------------------

identifier_summary <- data.frame(
  Metric = c(
    "Number of regions",
    "Number of cfMeDIP samples",
    "Duplicated sample IDs",
    "Duplicated region IDs",
    "Missing values",
    "Non-finite values"
  ),
  Value = c(
    nrow(medip),
    ncol(medip),
    sum(duplicated(colnames(medip))),
    sum(duplicated(rownames(medip))),
    sum(is.na(medip)),
    sum(!is.finite(medip), na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

cat("\n")
cat("====================================================\n")
cat("Global cfMeDIP integrity checks\n")
cat("====================================================\n")

print(identifier_summary)

write.csv(
  identifier_summary,
  file.path(
    dir_table,
    "cfMeDIP_global_identifier_QC.csv"
  ),
  row.names = FALSE
)


##------------------------------------------------------------
# 10. Global log2CPM summary
##------------------------------------------------------------

finite_values <- medip[
  is.finite(medip)
]

global_summary <- data.frame(
  Metric = c(
    "Minimum",
    "Q1",
    "Median",
    "Mean",
    "Q3",
    "Maximum",
    "SD"
  ),
  Value = c(
    min(finite_values),
    quantile(finite_values, 0.25),
    median(finite_values),
    mean(finite_values),
    quantile(finite_values, 0.75),
    max(finite_values),
    sd(finite_values)
  ),
  stringsAsFactors = FALSE
)

cat("\n")
cat("====================================================\n")
cat("Global cfMeDIP log2CPM summary\n")
cat("====================================================\n")

print(global_summary)

write.csv(
  global_summary,
  file.path(
    dir_table,
    "cfMeDIP_global_log2CPM_summary.csv"
  ),
  row.names = FALSE
)


##------------------------------------------------------------
# 11. Exact plasma 5-base <-> cfMeDIP mapping
#
# This mapping comes from the existing SOLID 5-base analysis.
##------------------------------------------------------------

matched_5base_medip <- read.delim(
  file.path(
    'result',
    "03_matched_tissue_plasma",
    "paired_filtered",
    "SOLID_plasma_sample_to_patient_mapping.tsv"
  ),
  stringsAsFactors = FALSE
)

colnames(matched_5base_medip) <- c("batch_number", "plasma_sample_id", 
                                   "Patient_ID", "MeDIP_ID")

##------------------------------------------------------------
# 12. Validate 5-base <-> cfMeDIP mapping
##------------------------------------------------------------

matched_5base_medip$MeDIP_in_matrix <- (
  matched_5base_medip$MeDIP_ID %in%
    colnames(medip)
)

cat("\n")
cat("====================================================\n")
cat("5-base <-> cfMeDIP mapping validation\n")
cat("====================================================\n")

print(
  matched_5base_medip
)

cat(
  "\nMapped cfMeDIP IDs found:",
  sum(
    matched_5base_medip$MeDIP_in_matrix
  ),
  "/",
  nrow(
    matched_5base_medip
  ),
  "\n"
)

if (
  !all(
    matched_5base_medip$MeDIP_in_matrix
  )
) {

  stop(
    "One or more mapped cfMeDIP IDs are absent from the cfMeDIP matrix."
  )
}

##------------------------------------------------------------
# 13. Define tumor EPIC patient IDs
##------------------------------------------------------------

epic_patient_ids <- unique(
  trimws(
    as.character(
      matched_meta$Subject
    )
  )
)

epic_patient_ids <- epic_patient_ids[
  !is.na(epic_patient_ids) &
    epic_patient_ids != ""
]

cat("\n")
cat("====================================================\n")
cat("Tumor EPIC cohort\n")
cat("====================================================\n")

cat(
  "Tumor EPIC patients:",
  length(
    epic_patient_ids
  ),
  "\n"
)

print(
  epic_patient_ids
)


##------------------------------------------------------------
# 14. Define final three-layer cohort
#
# Tumor EPIC
#     intersect
# Plasma 5-base
#     intersect
# Plasma cfMeDIP
##------------------------------------------------------------

three_layer_map <- matched_5base_medip[
  matched_5base_medip$Patient_ID %in%
    epic_patient_ids &
    matched_5base_medip$MeDIP_in_matrix,
  ,
  drop = FALSE
]

rownames(
  three_layer_map
) <- NULL


##------------------------------------------------------------
# Add tumor EPIC information
##------------------------------------------------------------

tumor_match <- match(
  three_layer_map$Patient_ID,
  matched_meta$Subject
)

three_layer_map$Tumor_Sample_Name <- matched_meta$Sample_Name[
  tumor_match
]

three_layer_map$Tumor_Array <- matched_meta$Array[
  tumor_match
]

three_layer_map$Tumor_Genome_build <- matched_meta$Genome_build[
  tumor_match
]


cat("\n")
cat("====================================================\n")
cat("FINAL EPIC + 5-BASE + cfMeDIP COHORT\n")
cat("====================================================\n")

cat(
  "Final matched patients:",
  nrow(
    three_layer_map
  ),
  "\n"
)

print(
  three_layer_map
)


##------------------------------------------------------------
# 15. Identify excluded 5-base/cfMeDIP samples
##------------------------------------------------------------

excluded_from_three_layer <- matched_5base_medip[
  !matched_5base_medip$Patient_ID %in%
    three_layer_map$Patient_ID,
  ,
  drop = FALSE
]

cat("\n")
cat("Excluded from three-layer cohort:\n")

print(
  excluded_from_three_layer
)


##------------------------------------------------------------
# Safety checks
##------------------------------------------------------------

if (
  nrow(
    three_layer_map
  ) != 13
) {

  stop(
    paste0(
      "Expected 13 final matched patients, but found ",
      nrow(
        three_layer_map
      ),
      ". Review the sample mapping before proceeding."
    )
  )
}

if (
  any(
    duplicated(
      three_layer_map$Patient_ID
    )
  )
) {

  stop(
    "Duplicated Patient_ID detected in final cohort."
  )
}

if (
  any(
    duplicated(
      three_layer_map$MeDIP_ID
    )
  )
) {

  stop(
    "Duplicated MeDIP_ID detected in final cohort."
  )
}


##------------------------------------------------------------
# Save matched sample information
##------------------------------------------------------------

write.csv(
  three_layer_map,
  file.path(
    dir_table,
    "SOLID_EPIC_5base_cfMeDIP_matched_samples.csv"
  ),
  row.names = FALSE
)

write.csv(
  excluded_from_three_layer,
  file.path(
    dir_table,
    "SOLID_three_layer_excluded_samples.csv"
  ),
  row.names = FALSE
)


##------------------------------------------------------------
# 16. Create final matched cfMeDIP matrix
##------------------------------------------------------------

medip_matched <- medip[
  ,
  three_layer_map$MeDIP_ID,
  drop = FALSE
]

stopifnot(
  ncol(
    medip_matched
  ) == 13
)

stopifnot(
  identical(
    colnames(
      medip_matched
    ),
    three_layer_map$MeDIP_ID
  )
)


##------------------------------------------------------------
# Rename cfMeDIP columns to patient IDs for downstream QC
##------------------------------------------------------------

colnames(
  medip_matched
) <- three_layer_map$Patient_ID


cat("\n")
cat("====================================================\n")
cat("FINAL MATCHED cfMeDIP MATRIX\n")
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
# Save matched cfMeDIP matrix now
##------------------------------------------------------------

saveRDS(
  medip_matched,
  file.path(
    dir_object,
    "SOLID_cfMeDIP_matched13_log2CPM_matrix.rds"
  )
)
