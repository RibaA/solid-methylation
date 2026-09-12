############################################################
## 03_prepare_SOLID_downstream_object.r
##
## SOLID DOWNSTREAM ANALYSIS SETUP
##
## Purpose:
##   1. Load frozen paired tumor-plasma object
##   2. Validate 13 matched patients and matrix alignment
##   3. Validate filtering settings
##   4. Load latest clinical metadata
##   5. Retain selected main clinical variables only
##   6. Validate clinical metadata against embedded metadata
##   7. Save one analysis-ready downstream object
##
## Clinical variables retained:
##   - Grade
##   - Age
##   - Sex
##   - ECOG
##   - Response_RANO
##   - PFS in months (pfs)
##   - PFS status
##   - OS in months (os)
##   - OS status
##
## Variable types:
##   Continuous:
##     Age, pfs, os
##
##   Categorical / event status:
##     Grade, Sex, ECOG, Response_RANO,
##     PFS_status, OS_status
##
## NO PCA / clustering / DMR testing is performed here.
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
## 1. PROJECT DIRECTORIES
############################################################

project_dir <- "C:/solid-methylation"

paired_dir <- file.path(
  project_dir,
  "result",
  "03_matched_tissue_plasma",
  "paired_filtered"
)

clinical_dir <- file.path(
  project_dir,
  "result",
  "04_clinical"
)

clinical_metadata_dir <- file.path(
  clinical_dir,
  "metadata"
)

clinical_table_dir <- file.path(
  clinical_dir,
  "tables"
)

downstream_dir <- file.path(
  project_dir,
  "result",
  "03_matched_tissue_plasma"
)

dir.create(
  clinical_metadata_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  clinical_table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
## 2. INPUT FILES
############################################################

paired_file <- file.path(
  paired_dir,
  "SOLID_paired_filtered_object.rds"
)

clinical_file <- file.path(
  project_dir,
  "data",
  "meta_data_research.csv"
)

stopifnot(
  file.exists(paired_file),
  file.exists(clinical_file)
)


############################################################
## 3. LOAD PAIRED OBJECT
############################################################

paired_obj <- readRDS(
  paired_file
)

cat(
  "\n============================================\n"
)

cat(
  "SOLID DOWNSTREAM SETUP\n"
)

cat(
  "============================================\n"
)


############################################################
## 4. REQUIRED OBJECT COMPONENTS
############################################################

required_components <- c(
  "tissue_beta",
  "plasma_beta",
  "delta_beta",
  "tissue_M",
  "plasma_M",
  "delta_M",
  "valid_pair_mask",
  "retained_annotation",
  "matched_metadata",
  "settings"
)

missing_components <- setdiff(
  required_components,
  names(paired_obj)
)

if (length(missing_components) > 0L) {

  stop(
    "Missing paired-object components: ",
    paste(
      missing_components,
      collapse = ", "
    )
  )
}

cat(
  "\nPaired-object components: PASS\n"
)


############################################################
## 5. MATRIX DIMENSIONS
############################################################

matrix_dims <- rbind(
  tissue_beta = dim(paired_obj$tissue_beta),
  plasma_beta = dim(paired_obj$plasma_beta),
  delta_beta = dim(paired_obj$delta_beta),
  tissue_M = dim(paired_obj$tissue_M),
  plasma_M = dim(paired_obj$plasma_M),
  delta_M = dim(paired_obj$delta_M),
  valid_pair_mask = dim(paired_obj$valid_pair_mask)
)

colnames(matrix_dims) <- c(
  "regions",
  "patients"
)

cat(
  "\nMatrix dimensions:\n"
)

print(
  matrix_dims
)

stopifnot(
  length(
    unique(
      matrix_dims[, "regions"]
    )
  ) == 1L,

  length(
    unique(
      matrix_dims[, "patients"]
    )
  ) == 1L,

  ncol(
    paired_obj$tissue_beta
  ) == 13L,

  nrow(
    paired_obj$tissue_beta
  ) == 124961L
)

cat(
  "\nMatrix dimension validation: PASS\n"
)


############################################################
## 6. PATIENT ORDER VALIDATION
############################################################

patient_ids <- colnames(
  paired_obj$tissue_beta
)

patient_order_pass <- all(
  identical(
    patient_ids,
    colnames(
      paired_obj$plasma_beta
    )
  ),

  identical(
    patient_ids,
    colnames(
      paired_obj$delta_beta
    )
  ),

  identical(
    patient_ids,
    colnames(
      paired_obj$tissue_M
    )
  ),

  identical(
    patient_ids,
    colnames(
      paired_obj$plasma_M
    )
  ),

  identical(
    patient_ids,
    colnames(
      paired_obj$delta_M
    )
  ),

  identical(
    patient_ids,
    colnames(
      paired_obj$valid_pair_mask
    )
  )
)

cat(
  "\nPatient order identical across matrices:",
  patient_order_pass,
  "\n"
)

stopifnot(
  patient_order_pass
)

cat(
  "\nMatched patients:\n"
)

print(
  patient_ids
)


############################################################
## 7. VERIFY DELTA DEFINITIONS
############################################################

beta_difference <- abs(
  paired_obj$delta_beta -
    (
      paired_obj$plasma_beta -
        paired_obj$tissue_beta
    )
)

max_beta_delta_error <- max(
  beta_difference,
  na.rm = TRUE
)


M_difference <- abs(
  paired_obj$delta_M -
    (
      paired_obj$plasma_M -
        paired_obj$tissue_M
    )
)

max_M_delta_error <- max(
  M_difference,
  na.rm = TRUE
)


cat(
  "\nMaximum ΔBeta reconstruction error:",
  max_beta_delta_error,
  "\n"
)

cat(
  "Maximum ΔM reconstruction error:",
  max_M_delta_error,
  "\n"
)

stopifnot(
  max_beta_delta_error < 1e-12,
  max_M_delta_error < 1e-12
)

cat(
  "Delta definitions: PASS\n"
)


############################################################
## 8. VERIFY PAIRED FILTERING SETTINGS
############################################################

cat(
  "\nPaired filtering settings:\n"
)

print(
  paired_obj$settings
)

stopifnot(

  as.numeric(
    paired_obj$settings$minimum_plasma_coverage
  ) == 10,

  as.numeric(
    paired_obj$settings$minimum_plasma_covered_cpgs
  ) == 3,

  as.numeric(
    paired_obj$settings$minimum_EPIC_probes
  ) == 2,

  as.numeric(
    paired_obj$settings$minimum_valid_pairs
  ) == 10,

  paired_obj$settings$delta_beta_definition ==
    "plasma_minus_tissue",

  paired_obj$settings$genome_build ==
    "hg38"
)

cat(
  "\nFiltering settings validation: PASS\n"
)


############################################################
## 9. VALIDATE RETAINED REGION ANNOTATION
############################################################

annotation <- as.data.frame(
  paired_obj$retained_annotation
)

stopifnot(

  nrow(annotation) ==
    nrow(
      paired_obj$tissue_beta
    ),

  !anyDuplicated(
    annotation$region_id
  ),

  all(
    annotation$n_EPIC_probes >= 2
  ),

  all(
    annotation$n_valid_pairs >= 10
  )
)

cat(
  "\nRetained region annotation: PASS\n"
)

cat(
  "Regions:",
  nrow(annotation),
  "\n"
)

cat(
  "Valid-pairs/region summary:\n"
)

print(
  summary(
    annotation$n_valid_pairs
  )
)


############################################################
## 10. READ LATEST CLINICAL METADATA
##
## check.names = TRUE makes duplicate column names unique.
############################################################

clinical_raw <- read.csv(
  clinical_file,
  stringsAsFactors = FALSE,
  check.names = TRUE
)

cat(
  "\nClinical metadata dimensions:",
  nrow(clinical_raw),
  "x",
  ncol(clinical_raw),
  "\n"
)

cat(
  "\nClinical columns:\n"
)

print(
  names(clinical_raw)
)


############################################################
## 11. VALIDATE PRIMARY SUBJECT COLUMN
############################################################

if (!"Subject" %in% names(clinical_raw)) {

  stop(
    "Primary Subject column not found in clinical metadata."
  )
}

clinical_raw$Subject <- trimws(
  as.character(
    clinical_raw$Subject
  )
)

stopifnot(
  !anyNA(
    clinical_raw$Subject
  ),

  !anyDuplicated(
    clinical_raw$Subject
  )
)


############################################################
## 12. MAIN CLINICAL VARIABLES TO RETAIN
############################################################

clinical_variables <- c(
  "Subject",
  "Grade",
  "Age",
  "Sex",
  "ECOG",
  "Response_RANO",
  "pfs",
  "PFS_status",
  "os",
  "OS_status"
)


available_clinical_variables <- intersect(
  clinical_variables,
  names(clinical_raw)
)

missing_clinical_variables <- setdiff(
  clinical_variables,
  names(clinical_raw)
)


cat(
  "\nClinical variables retained:\n"
)

print(
  available_clinical_variables
)


if (length(missing_clinical_variables) > 0L) {

  cat(
    "\nClinical variables not available:\n"
  )

  print(
    missing_clinical_variables
  )
}


############################################################
## 13. CREATE CLEAN CLINICAL TABLE
############################################################

clinical_clean <- clinical_raw[
  ,
  available_clinical_variables,
  drop = FALSE
]

names(clinical_clean)[
  names(clinical_clean) == "Subject"
] <- "patient_id"


############################################################
## 14. SUBSET TO THE 13 MATCHED PATIENTS
############################################################

clinical_matched <- clinical_clean[
  match(
    patient_ids,
    clinical_clean$patient_id
  ),
  ,
  drop = FALSE
]

rownames(
  clinical_matched
) <- NULL


############################################################
## 15. CONFIRM ALL 13 PATIENTS FOUND
############################################################

if (
  anyNA(
    clinical_matched$patient_id
  )
) {

  missing_patients <- patient_ids[
    is.na(
      clinical_matched$patient_id
    )
  ]

  stop(
    "Matched patients missing from clinical metadata: ",
    paste(
      missing_patients,
      collapse = ", "
    )
  )
}


stopifnot(
  identical(
    clinical_matched$patient_id,
    patient_ids
  )
)


cat(
  "\nClinical matching: PASS\n"
)

cat(
  "Matched clinical patients:",
  nrow(clinical_matched),
  "\n"
)


############################################################
## 16. VARIABLE TYPE CLEANING
############################################################

############################################################
## Continuous variables
############################################################

continuous_fields <- intersect(
  c(
    "Age",
    "pfs",
    "os"
  ),
  names(clinical_matched)
)


for (v in continuous_fields) {

  clinical_matched[[v]] <- suppressWarnings(
    as.numeric(
      clinical_matched[[v]]
    )
  )
}


############################################################
## Categorical / event-status variables
############################################################

categorical_fields <- intersect(
  c(
    "Grade",
    "Sex",
    "ECOG",
    "Response_RANO",
    "PFS_status",
    "OS_status"
  ),
  names(clinical_matched)
)


for (v in categorical_fields) {

  clinical_matched[[v]] <- trimws(
    as.character(
      clinical_matched[[v]]
    )
  )

  clinical_matched[[v]][
    clinical_matched[[v]] == ""
  ] <- NA_character_
}


############################################################
## 17. CLINICAL COMPLETENESS SUMMARY
############################################################

clinical_completeness <- data.frame(

  variable =
    names(
      clinical_matched
    ),

  n_nonmissing =
    vapply(
      clinical_matched,
      function(x) {

        sum(
          !is.na(x) &
            trimws(
              as.character(x)
            ) != ""
        )
      },
      integer(1)
    ),

  n_missing =
    vapply(
      clinical_matched,
      function(x) {

        sum(
          is.na(x) |
            trimws(
              as.character(x)
            ) == ""
        )
      },
      integer(1)
    ),

  stringsAsFactors = FALSE
)


clinical_completeness$pct_complete <-
  100 *
  clinical_completeness$n_nonmissing /
  nrow(
    clinical_matched
  )


cat(
  "\nClinical metadata completeness:\n"
)

print(
  clinical_completeness
)


############################################################
## 18. CLINICAL VALUE SUMMARIES
############################################################

cat(
  "\nAge summary:\n"
)

if ("Age" %in% names(clinical_matched)) {

  print(
    summary(
      clinical_matched$Age
    )
  )
}


cat(
  "\nSex distribution:\n"
)

if ("Sex" %in% names(clinical_matched)) {

  print(
    table(
      clinical_matched$Sex,
      useNA = "ifany"
    )
  )
}


cat(
  "\nGrade distribution:\n"
)

if ("Grade" %in% names(clinical_matched)) {

  print(
    table(
      clinical_matched$Grade,
      useNA = "ifany"
    )
  )
}


cat(
  "\nECOG distribution:\n"
)

if ("ECOG" %in% names(clinical_matched)) {

  print(
    table(
      clinical_matched$ECOG,
      useNA = "ifany"
    )
  )
}


cat(
  "\nResponse RANO distribution:\n"
)

if ("Response_RANO" %in% names(clinical_matched)) {

  print(
    table(
      clinical_matched$Response_RANO,
      useNA = "ifany"
    )
  )
}


cat(
  "\nPFS (months) summary:\n"
)

if ("pfs" %in% names(clinical_matched)) {

  print(
    summary(
      clinical_matched$pfs
    )
  )
}


cat(
  "\nPFS status distribution:\n"
)

if ("PFS_status" %in% names(clinical_matched)) {

  print(
    table(
      clinical_matched$PFS_status,
      useNA = "ifany"
    )
  )
}


cat(
  "\nOS (months) summary:\n"
)

if ("os" %in% names(clinical_matched)) {

  print(
    summary(
      clinical_matched$os
    )
  )
}


cat(
  "\nOS status distribution:\n"
)

if ("OS_status" %in% names(clinical_matched)) {

  print(
    table(
      clinical_matched$OS_status,
      useNA = "ifany"
    )
  )
}


############################################################
## 19. COMPARE WITH CLINICAL DATA ALREADY EMBEDDED
##
## Fresh meta_data_research.csv is authoritative.
## Embedded metadata are used only as a consistency check.
############################################################

embedded_meta <- as.data.frame(
  paired_obj$matched_metadata
)


embedded_meta <- embedded_meta[
  match(
    patient_ids,
    embedded_meta$patient_id
  ),
  ,
  drop = FALSE
]


comparison_variables <- intersect(
  c(
    "Grade",
    "Age",
    "Sex",
    "ECOG",
    "Response_RANO",
    "pfs",
    "PFS_status",
    "os",
    "OS_status"
  ),
  intersect(
    names(clinical_matched),
    names(embedded_meta)
  )
)


clinical_comparison <- data.frame()


for (v in comparison_variables) {

  new_value <- as.character(
    clinical_matched[[v]]
  )

  old_value <- as.character(
    embedded_meta[[v]]
  )


  new_missing <-
    is.na(new_value) |
    trimws(new_value) == ""

  old_missing <-
    is.na(old_value) |
    trimws(old_value) == ""


  same <- (
    new_value == old_value
  ) |
    (
      new_missing &
        old_missing
    )


  same[
    is.na(same)
  ] <- FALSE


  temp <- data.frame(

    patient_id =
      patient_ids,

    variable =
      v,

    current_clinical_value =
      new_value,

    embedded_value =
      old_value,

    identical =
      same,

    stringsAsFactors = FALSE
  )


  clinical_comparison <- rbind(
    clinical_comparison,
    temp
  )
}


cat(
  "\nClinical comparison with paired object:\n"
)

print(
  table(
    clinical_comparison$identical
  )
)


############################################################
## 20. SHOW ANY CLINICAL DISCREPANCIES
############################################################

clinical_discrepancies <- clinical_comparison[
  !clinical_comparison$identical,
  ,
  drop = FALSE
]


cat(
  "\nClinical discrepancies:\n"
)

if (nrow(clinical_discrepancies) == 0L) {

  cat(
    "None\n"
  )

} else {

  print(
    clinical_discrepancies
  )
}


############################################################
## 21. SAVE CLEAN CLINICAL TABLE
############################################################

clinical_matched_file <- file.path(
  clinical_metadata_dir,
  "SOLID_13_matched_clinical_metadata.tsv"
)


write.table(
  clinical_matched,
  clinical_matched_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
## 22. SAVE CLINICAL COMPLETENESS
############################################################

clinical_completeness_file <- file.path(
  clinical_table_dir,
  "SOLID_13_matched_clinical_metadata_completeness.tsv"
)


write.table(
  clinical_completeness,
  clinical_completeness_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
## 23. SAVE CLINICAL VALIDATION TABLE
############################################################

clinical_comparison_file <- file.path(
  clinical_table_dir,
  "SOLID_clinical_metadata_validation_vs_paired_object.tsv"
)


write.table(
  clinical_comparison,
  clinical_comparison_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
## 24. CREATE SAMPLE-LEVEL METADATA
##
## 26 rows:
##   13 tumor samples
##   13 plasma samples
##
## This will be used for:
##   - PCA
##   - clustering
##   - heatmap annotations
##   - exploratory visualizations
############################################################

matched_metadata <- as.data.frame(
  paired_obj$matched_metadata
)


matched_metadata <- matched_metadata[
  match(
    patient_ids,
    matched_metadata$patient_id
  ),
  ,
  drop = FALSE
]


stopifnot(
  identical(
    matched_metadata$patient_id,
    patient_ids
  )
)


############################################################
## Tumor samples
############################################################

tissue_sample_metadata <- clinical_matched

tissue_sample_metadata$sample_id <-
  matched_metadata$tissue_sample_id

tissue_sample_metadata$sample_type <-
  "Tumor"


############################################################
## Plasma samples
############################################################

plasma_sample_metadata <- clinical_matched

plasma_sample_metadata$sample_id <-
  matched_metadata$plasma_sample_id

plasma_sample_metadata$sample_type <-
  "Plasma"


############################################################
## Combine
############################################################

sample_metadata <- rbind(
  tissue_sample_metadata,
  plasma_sample_metadata
)


sample_metadata <- sample_metadata[
  ,
  c(
    "patient_id",
    "sample_id",
    "sample_type",
    setdiff(
      names(sample_metadata),
      c(
        "patient_id",
        "sample_id",
        "sample_type"
      )
    )
  ),
  drop = FALSE
]


rownames(
  sample_metadata
) <- NULL


stopifnot(
  nrow(
    sample_metadata
  ) == 26L
)


cat(
  "\nSample-level metadata:\n"
)

print(
  sample_metadata
)


############################################################
## 25. CREATE DOWNSTREAM ANALYSIS OBJECT
############################################################

analysis_ready <- list(

  tissue_beta =
    paired_obj$tissue_beta,

  plasma_beta =
    paired_obj$plasma_beta,

  delta_beta =
    paired_obj$delta_beta,

  tissue_M =
    paired_obj$tissue_M,

  plasma_M =
    paired_obj$plasma_M,

  delta_M =
    paired_obj$delta_M,

  valid_pair_mask =
    paired_obj$valid_pair_mask,

  region_annotation =
    annotation,

  patient_metadata =
    clinical_matched,

  sample_metadata =
    sample_metadata,

  matched_metadata =
    matched_metadata,

  settings =
    paired_obj$settings,

  provenance = list(

    paired_object =
      paired_file,

    clinical_metadata =
      clinical_file,

    n_patients =
      length(
        patient_ids
      ),

    n_regions =
      nrow(
        annotation
      ),

    delta_beta_definition =
      "plasma_minus_tissue",

    creation_date =
      as.character(
        Sys.Date()
      )
  )
)


############################################################
## 26. SAVE ANALYSIS-READY OBJECT
############################################################

analysis_ready_file <- file.path(
  downstream_dir,
  "SOLID_downstream_analysis_ready.rds"
)


saveRDS(
  analysis_ready,
  analysis_ready_file
)


############################################################
## 27. FINAL VALIDATION
############################################################

check_obj <- readRDS(
  analysis_ready_file
)


stopifnot(

  ncol(
    check_obj$tissue_beta
  ) == 13L,

  ncol(
    check_obj$plasma_beta
  ) == 13L,

  nrow(
    check_obj$tissue_beta
  ) == 124961L,

  nrow(
    check_obj$patient_metadata
  ) == 13L,

  nrow(
    check_obj$sample_metadata
  ) == 26L,

  identical(
    colnames(
      check_obj$tissue_beta
    ),
    check_obj$patient_metadata$patient_id
  )
)


cat(
  "\nFinal downstream object validation: PASS\n"
)


############################################################
## 28. FINAL SUMMARY
############################################################

cat(
  "\n============================================\n"
)

cat(
  "SOLID DOWNSTREAM SETUP COMPLETE\n"
)

cat(
  "============================================\n"
)

cat(
  "\nMatched patients:",
  length(
    patient_ids
  ),
  "\n"
)

cat(
  "Retained regions:",
  nrow(
    annotation
  ),
  "\n"
)

cat(
  "Clinical patients:",
  nrow(
    clinical_matched
  ),
  "\n"
)

cat(
  "Sample-level metadata rows:",
  nrow(
    sample_metadata
  ),
  "\n"
)

cat(
  "\nClinical variables retained:\n"
)

print(
  names(
    clinical_matched
  )
)

cat(
  "\nClinical completeness:\n"
)

print(
  clinical_completeness
)

cat(
  "\nClinical comparison with existing paired object:\n"
)

print(
  table(
    clinical_comparison$identical
  )
)

cat(
  "\nAnalysis-ready object:\n",
  analysis_ready_file,
  "\n",
  sep = ""
)

cat(
  "\nClean clinical metadata:\n",
  clinical_matched_file,
  "\n",
  sep = ""
)

cat(
  "\nClinical validation table:\n",
  clinical_comparison_file,
  "\n",
  sep = ""
)