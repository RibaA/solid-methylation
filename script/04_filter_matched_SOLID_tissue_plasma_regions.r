##------------------------------------------------------------
## 04. SOLID FILTER MATCHED TISSUE–PLASMA 1-KB REGIONS
##
## Purpose:
##   1. Load projected tumour EPIC 1-kb values from script 03
##   2. Load plasma 5-base Beta, coverage, and covered-CpG matrices
##   3. Match SOLID tumour–plasma patient pairs
##   4. merge clinical metadata for downstream figures
##   5. apply paired technical and regional filtering
##   6. save paired tissue, plasma, delta-Beta, and delta-M matrices
##
## Delta Beta = plasma - tissue
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

expected_matched_patients <- 13L

minimum_plasma_coverage <- 10L
minimum_plasma_covered_cpgs <- 3L
minimum_EPIC_probes <- 2L
minimum_valid_pairs <- 10L

beta_epsilon <- 1e-6
save_valid_pair_mask <- TRUE

##------------------------------------------------------------
## 2. PACKAGES
##------------------------------------------------------------

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

##------------------------------------------------------------
## 3. DIRECTORIES
##------------------------------------------------------------

projection_dir <- file.path(
  project_dir,
  "result",
  "tissue-plasma",
  "matched_preparation"
)

plasma_dir <- file.path(
  project_dir,
  "data",
  "solid-5base"
)

paired_dir <- file.path(
  project_dir,
  "result",
  "tissue-plasma",
  "paired_filtered"
)

clinical_metadata_file <- file.path(
  project_dir,
  "data",
  "meta_data_research.csv"
)

dir.create(
  paired_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(
  dir.exists(projection_dir),
  dir.exists(plasma_dir),
  dir.exists(paired_dir)
)

##------------------------------------------------------------
## 4. INPUT FILES
##------------------------------------------------------------

tissue_beta_file <- file.path(
  projection_dir,
  "SOLID_tissue_1kb_beta.rds"
)

projected_annotation_file <- file.path(
  projection_dir,
  "SOLID_tissue_1kb_region_annotation.tsv.gz"
)

tissue_metadata_file <- file.path(
  projection_dir,
  "SOLID_tissue_harmonized_metadata.tsv"
)

plasma_beta_file <- file.path(
  plasma_dir,
  "SOLID_1kb_beta_cov5_n7_standard_chr.rds"
)

plasma_coverage_file <- file.path(
  plasma_dir,
  "SOLID_1kb_coverage_cov5_n7_standard_chr.rds"
)

plasma_cpg_file <- file.path(
  plasma_dir,
  "SOLID_1kb_covered_CpGs_cov5_n7_standard_chr.rds"
)

plasma_processing_metadata_file <- file.path(
  plasma_dir,
  "SOLID_1kb_processing_metadata.rds"
)

plasma_qc_file <- file.path(
  plasma_dir,
  "SOLID_5base_sample_QC.tsv"
)

required_files <- c(
  tissue_beta_file,
  projected_annotation_file,
  tissue_metadata_file,
  plasma_beta_file,
  plasma_coverage_file,
  plasma_cpg_file,
  plasma_processing_metadata_file,
  plasma_qc_file
)

file_check <- data.table(
  object = c(
    "tissue_beta",
    "projected_annotation",
    "tissue_metadata",
    "plasma_beta",
    "plasma_coverage",
    "plasma_covered_cpgs",
    "plasma_processing_metadata",
    "plasma_qc"
  ),
  path = required_files,
  exists = file.exists(required_files)
)

print(file_check)

if (any(!file_check$exists)) {
  stop(
    "Missing required file(s):\n",
    paste(
      file_check[exists == FALSE, path],
      collapse = "\n"
    )
  )
}

##------------------------------------------------------------
## 5. HELPERS
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

detect_column <- function(
    dat,
    candidates,
    table_name,
    required = TRUE
) {
  detected <- candidates[
    candidates %in% names(dat)
  ]

  if (length(detected) == 0L) {
    if (required) {
      stop(
        "Could not identify a required column in ",
        table_name,
        ".\nTried: ",
        paste(candidates, collapse = ", "),
        "\nAvailable columns:\n",
        paste(names(dat), collapse = ", ")
      )
    }

    return(NULL)
  }

  detected[[1]]
}

extract_SOLID_id <- function(x) {
  x <- toupper(as.character(x))
  x <- gsub("_", "-", x)

  out <- rep(
    NA_character_,
    length(x)
  )

  has_id <- grepl(
    "SOLID-[0-9]+",
    x
  )

  out[has_id] <- sub(
    ".*?(SOLID-[0-9]+).*",
    "\\1",
    x[has_id]
  )

  out
}

safe_beta_to_M <- function(
    beta,
    epsilon = 1e-6
) {
  beta <- pmin(
    pmax(beta, epsilon),
    1 - epsilon
  )

  log2(
    beta /
      (1 - beta)
  )
}

safe_row_medians <- function(x) {
  out <- matrixStats::rowMedians(
    x,
    na.rm = TRUE
  )

  out[!is.finite(out)] <- NA_real_
  out
}

safe_row_means <- function(x) {
  out <- rowMeans(
    x,
    na.rm = TRUE
  )

  out[!is.finite(out)] <- NA_real_
  out
}

##------------------------------------------------------------
## 6. LOAD PROJECTED TISSUE DATA
##------------------------------------------------------------

message_header(
  "LOADING PROJECTED TISSUE DATA"
)

tissue_beta <- readRDS(
  tissue_beta_file
)

projected_annotation <- fread(
  projected_annotation_file
)

tissue_metadata <- fread(
  tissue_metadata_file
)

stopifnot(
  is.matrix(tissue_beta),
  nrow(projected_annotation) ==
    nrow(tissue_beta),
  identical(
    projected_annotation$region_id,
    rownames(tissue_beta)
  ),
  "plasma_row_index" %in%
    names(projected_annotation),
  "n_EPIC_probes" %in%
    names(projected_annotation)
)

tissue_sample_column <- detect_column(
  tissue_metadata,
  c(
    "Sample_Name",
    "sample_id",
    "SampleID",
    "sample"
  ),
  "tissue metadata"
)

tissue_patient_column <- detect_column(
  tissue_metadata,
  c(
    "Subject",
    "patient_id",
    "SolidID",
    "SOLID_ID",
    "subject_id",
    "Study"
  ),
  "tissue metadata"
)

tissue_metadata[
  ,
  tissue_sample_id :=
    as.character(
      get(tissue_sample_column)
    )
]

tissue_metadata[
  ,
  patient_id :=
    extract_SOLID_id(
      get(tissue_patient_column)
    )
]

tissue_metadata[
  is.na(patient_id),
  patient_id :=
    extract_SOLID_id(
      tissue_sample_id
    )
]

if (anyNA(tissue_metadata$patient_id)) {
  stop(
    "Some tissue samples lack a recognizable SOLID ID:\n",
    paste(
      tissue_metadata[
        is.na(patient_id),
        tissue_sample_id
      ],
      collapse = ", "
    )
  )
}

tissue_metadata <- tissue_metadata[
  match(
    colnames(tissue_beta),
    tissue_sample_id
  )
]

if (anyNA(tissue_metadata$tissue_sample_id)) {
  stop(
    "Tissue metadata could not be matched to all tissue matrix columns."
  )
}

if (anyDuplicated(tissue_metadata$patient_id)) {
  duplicated_ids <- unique(
    tissue_metadata[
      duplicated(patient_id) |
        duplicated(patient_id, fromLast = TRUE),
      patient_id
    ]
  )

  stop(
    "Multiple tissue samples found for patient(s): ",
    paste(duplicated_ids, collapse = ", "),
    "\nResolve duplicate patient samples before paired analysis."
  )
}

cat(
  "Tissue matrix:",
  paste(dim(tissue_beta), collapse = " x "),
  "\n"
)

cat(
  "Unique tissue patients:",
  uniqueN(tissue_metadata$patient_id),
  "\n"
)

##------------------------------------------------------------
## 7. LOAD PLASMA DATA
##------------------------------------------------------------

message_header(
  "LOADING PLASMA DATA"
)

plasma_beta <- readRDS(
  plasma_beta_file
)

plasma_coverage <- readRDS(
  plasma_coverage_file
)

plasma_covered_cpgs <- readRDS(
  plasma_cpg_file
)

plasma_processing_metadata <- readRDS(
  plasma_processing_metadata_file
)

plasma_qc <- fread(
  plasma_qc_file
)

stopifnot(
  is.matrix(plasma_beta),
  is.matrix(plasma_coverage),
  is.matrix(plasma_covered_cpgs),
  identical(
    dim(plasma_beta),
    dim(plasma_coverage)
  ),
  identical(
    dim(plasma_beta),
    dim(plasma_covered_cpgs)
  ),
  identical(
    colnames(plasma_beta),
    colnames(plasma_coverage)
  ),
  identical(
    colnames(plasma_beta),
    colnames(plasma_covered_cpgs)
  ),
  max(projected_annotation$plasma_row_index) <=
    nrow(plasma_beta)
)

plasma_metadata <- data.table(
  plasma_sample_id =
    colnames(plasma_beta),
  patient_id =
    extract_SOLID_id(
      colnames(plasma_beta)
    )
)

## Try to recover patient IDs from the plasma QC table if needed
qc_sample_column <- detect_column(
  plasma_qc,
  c(
    "sample_id",
    "Sample_ID",
    "SampleID",
    "sample",
    "Sample",
    "sample_name",
    "Sample_Name"
  ),
  "plasma QC",
  required = FALSE
)

qc_patient_column <- detect_column(
  plasma_qc,
  c(
    "patient_id",
    "Patient_ID",
    "Subject",
    "SolidID",
    "SOLID_ID",
    "subject_id",
    "Study"
  ),
  "plasma QC",
  required = FALSE
)

if (!is.null(qc_sample_column)) {
  plasma_qc[
    ,
    plasma_sample_id :=
      as.character(
        get(qc_sample_column)
      )
  ]

  if (!is.null(qc_patient_column)) {
    plasma_qc[
      ,
      qc_patient_id :=
        extract_SOLID_id(
          get(qc_patient_column)
        )
    ]
  } else {
    plasma_qc[
      ,
      qc_patient_id :=
        extract_SOLID_id(
          plasma_sample_id
        )
    ]

  }

  qc_id_table <- unique(
    plasma_qc[
      ,
      .(
        plasma_sample_id,
        qc_patient_id
      )
    ]
  )

  plasma_metadata <- merge(
    plasma_metadata,
    qc_id_table,
    by = "plasma_sample_id",
    all.x = TRUE,
    sort = FALSE
  )

  plasma_metadata[
    is.na(patient_id),
    patient_id :=
      qc_patient_id
  ]

  plasma_metadata[
    ,
    qc_patient_id := NULL
  ]
}

plasma_metadata <- plasma_metadata[
  match(
    colnames(plasma_beta),
    plasma_sample_id
  )
]

if (anyNA(plasma_metadata$patient_id)) {
  stop(
    "Could not identify SOLID patient IDs for plasma samples:\n",
    paste(
      plasma_metadata[
        is.na(patient_id),
        plasma_sample_id
      ],
      collapse = ", "
    ),
    "\nInspect colnames(plasma_beta) and the plasma QC table."
  )
}

if (anyDuplicated(plasma_metadata$patient_id)) {
  duplicated_ids <- unique(
    plasma_metadata[
      duplicated(patient_id) |
        duplicated(patient_id, fromLast = TRUE),
      patient_id
    ]
  )

  stop(
    "Multiple plasma samples found for patient(s): ",
    paste(duplicated_ids, collapse = ", "),
    "\nResolve duplicate patient samples before paired analysis."
  )
}

cat(
  "Plasma matrix:",
  paste(dim(plasma_beta), collapse = " x "),
  "\n"
)

cat(
  "Unique plasma patients:",
  uniqueN(plasma_metadata$patient_id),
  "\n"
)

##------------------------------------------------------------
## 8. IDENTIFY MATCHED PATIENTS
##------------------------------------------------------------

message_header(
  "MATCHING TISSUE AND PLASMA PATIENTS"
)

matched_patients <- sort(
  intersect(
    tissue_metadata$patient_id,
    plasma_metadata$patient_id
  )
)

tissue_only_patients <- sort(
  setdiff(
    tissue_metadata$patient_id,
    plasma_metadata$patient_id
  )
)

plasma_only_patients <- sort(
  setdiff(
    plasma_metadata$patient_id,
    tissue_metadata$patient_id
  )
)

cat(
  "Matched patients:",
  length(matched_patients),
  "\n"
)

cat(
  "Tissue-only patients:",
  length(tissue_only_patients),
  "\n"
)

cat(
  "Plasma-only patients:",
  length(plasma_only_patients),
  "\n"
)

if (
  length(matched_patients) !=
    expected_matched_patients
) {
  stop(
    "Expected ",
    expected_matched_patients,
    " matched patients but found ",
    length(matched_patients),
    ".\nMatched IDs:\n",
    paste(matched_patients, collapse = ", "),
    "\nTissue-only IDs:\n",
    paste(tissue_only_patients, collapse = ", "),
    "\nPlasma-only IDs:\n",
    paste(plasma_only_patients, collapse = ", ")
  )
}

matched_samples <- merge(
  tissue_metadata[
    patient_id %in% matched_patients,
    .(
      patient_id,
      tissue_sample_id
    )
  ],
  plasma_metadata[
    patient_id %in% matched_patients,
    .(
      patient_id,
      plasma_sample_id
    )
  ],
  by = "patient_id",
  all = FALSE,
  sort = TRUE
)

stopifnot(
  nrow(matched_samples) ==
    expected_matched_patients,
  !anyDuplicated(matched_samples$patient_id)
)

## Add the complete tissue metadata
matched_metadata <- merge(
  matched_samples,
  tissue_metadata,
  by = c(
    "patient_id",
    "tissue_sample_id"
  ),
  all.x = TRUE,
  sort = FALSE
)

## Add external clinical metadata if available
if (file.exists(clinical_metadata_file)) {
  clinical_metadata <- fread(
    clinical_metadata_file
  )

  clinical_patient_column <- detect_column(
    clinical_metadata,
    c(
      "Subject",
      "Study",
      "patient_id",
      "SolidID",
      "SOLID_ID",
      "subject_id"
    ),
    "clinical metadata",
    required = FALSE
  )

  if (!is.null(clinical_patient_column)) {
    clinical_metadata[
      ,
      patient_id :=
        extract_SOLID_id(
          get(clinical_patient_column)
        )
    ]

    clinical_metadata <- clinical_metadata[
      !is.na(patient_id)
    ]

    clinical_metadata <- clinical_metadata[
      !duplicated(patient_id)
    ]

    overlapping_columns <- intersect(
      setdiff(
        names(clinical_metadata),
        "patient_id"
      ),
      names(matched_metadata)
    )

    if (length(overlapping_columns) > 0L) {
      setnames(
        clinical_metadata,
        overlapping_columns,
        paste0(
          overlapping_columns,
          "_clinical"
        )
      )
    }

    matched_metadata <- merge(
      matched_metadata,
      clinical_metadata,
      by = "patient_id",
      all.x = TRUE,
      sort = FALSE
    )
  }
}

matched_metadata <- matched_metadata[
  match(
    matched_patients,
    patient_id
  )
]

stopifnot(
  identical(
    matched_metadata$patient_id,
    matched_patients
  )
)

fwrite(
  matched_metadata,
  file.path(
    paired_dir,
    "SOLID_matched_tissue_plasma_metadata.tsv"
  ),
  sep = "\t"
)

unmatched_table <- rbindlist(
  list(
    data.table(
      dataset = "tissue_only",
      patient_id = tissue_only_patients
    ),
    data.table(
      dataset = "plasma_only",
      patient_id = plasma_only_patients
    )
  ),
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  unmatched_table,
  file.path(
    paired_dir,
    "SOLID_unmatched_patients.tsv"
  ),
  sep = "\t"
)

##------------------------------------------------------------
## 9. MATCH SAMPLE COLUMNS
##------------------------------------------------------------

message_header(
  "CREATING MATCHED SAMPLE INDICES"
)

tissue_column_indices <- match(
  matched_metadata$tissue_sample_id,
  colnames(tissue_beta)
)

plasma_column_indices <- match(
  matched_metadata$plasma_sample_id,
  colnames(plasma_beta)
)

stopifnot(
  !anyNA(tissue_column_indices),
  !anyNA(plasma_column_indices)
)

patient_order <- matched_metadata$patient_id
n_pairs <- length(patient_order)

##------------------------------------------------------------
## 10. SUBSET TO PROJECTED REGIONS AND MATCHED PATIENTS
##------------------------------------------------------------

message_header(
  "SUBSETTING MATCHED MATRICES"
)

plasma_row_indices <- projected_annotation$
  plasma_row_index

matched_tissue_beta <- tissue_beta[
  ,
  tissue_column_indices,
  drop = FALSE
]

matched_plasma_beta <- plasma_beta[
  plasma_row_indices,
  plasma_column_indices,
  drop = FALSE
]

matched_plasma_coverage <- plasma_coverage[
  plasma_row_indices,
  plasma_column_indices,
  drop = FALSE
]

matched_plasma_cpgs <- plasma_covered_cpgs[
  plasma_row_indices,
  plasma_column_indices,
  drop = FALSE
]

rownames(matched_plasma_beta) <-
  projected_annotation$region_id

rownames(matched_plasma_coverage) <-
  projected_annotation$region_id

rownames(matched_plasma_cpgs) <-
  projected_annotation$region_id

colnames(matched_tissue_beta) <-
  patient_order

colnames(matched_plasma_beta) <-
  patient_order

colnames(matched_plasma_coverage) <-
  patient_order

colnames(matched_plasma_cpgs) <-
  patient_order

stopifnot(
  identical(
    dim(matched_tissue_beta),
    dim(matched_plasma_beta)
  ),
  identical(
    rownames(matched_tissue_beta),
    rownames(matched_plasma_beta)
  ),
  identical(
    colnames(matched_tissue_beta),
    colnames(matched_plasma_beta)
  )
)

cat(
  "Matched matrices:",
  paste(dim(matched_tissue_beta), collapse = " x "),
  "\n"
)

## Release full matrices
rm(
  tissue_beta,
  plasma_beta,
  plasma_coverage,
  plasma_covered_cpgs
)

invisible(
  gc(verbose = FALSE)
)

##------------------------------------------------------------
## 11. APPLY PAIR-SPECIFIC FILTERING
##------------------------------------------------------------

message_header(
  "APPLYING PAIRED TECHNICAL FILTERS"
)

valid_pair_mask <- (
  matched_plasma_coverage >=
    minimum_plasma_coverage
) &
  (
    matched_plasma_cpgs >=
      minimum_plasma_covered_cpgs
  ) &
  is.finite(
    matched_plasma_beta
  ) &
  is.finite(
    matched_tissue_beta
  )

n_valid_pairs <- rowSums(
  valid_pair_mask
)

region_pass_EPIC <- projected_annotation$
  n_EPIC_probes >=
  minimum_EPIC_probes

region_pass_pairs <- n_valid_pairs >=
  minimum_valid_pairs

retain_region <- region_pass_EPIC &
  region_pass_pairs

cat(
  "Projected regions evaluated:",
  format(length(retain_region), big.mark = ","),
  "\n"
)

cat(
  "Regions with >=",
  minimum_EPIC_probes,
  "EPIC probes:",
  format(sum(region_pass_EPIC), big.mark = ","),
  "\n"
)

cat(
  "Regions with >=",
  minimum_valid_pairs,
  "valid pairs:",
  format(sum(region_pass_pairs), big.mark = ","),
  "\n"
)

cat(
  "Regions passing both criteria:",
  format(sum(retain_region), big.mark = ","),
  "\n"
)

if (!any(retain_region)) {
  stop(
    "No regions passed the selected filtering criteria."
  )
}

##------------------------------------------------------------
## 12. CREATE FILTERED PAIRED MATRICES
##------------------------------------------------------------

message_header(
  "CREATING FILTERED PAIRED MATRICES"
)

paired_tissue_beta <- matched_tissue_beta[
  retain_region,
  ,
  drop = FALSE
]

paired_plasma_beta <- matched_plasma_beta[
  retain_region,
  ,
  drop = FALSE
]

paired_valid_pair_mask <- valid_pair_mask[
  retain_region,
  ,
  drop = FALSE
]

paired_tissue_beta[
  !paired_valid_pair_mask
] <- NA_real_

paired_plasma_beta[
  !paired_valid_pair_mask
] <- NA_real_

paired_delta_beta <- paired_plasma_beta -
  paired_tissue_beta

paired_tissue_M <- safe_beta_to_M(
  paired_tissue_beta,
  epsilon = beta_epsilon
)

paired_plasma_M <- safe_beta_to_M(
  paired_plasma_beta,
  epsilon = beta_epsilon
)

paired_delta_M <- paired_plasma_M -
  paired_tissue_M

stopifnot(
  identical(
    dim(paired_tissue_beta),
    dim(paired_plasma_beta)
  ),
  identical(
    dim(paired_tissue_beta),
    dim(paired_delta_beta)
  ),
  identical(
    rownames(paired_tissue_beta),
    rownames(paired_plasma_beta)
  ),
  identical(
    colnames(paired_tissue_beta),
    colnames(paired_plasma_beta)
  ),
  all(
    is.na(paired_tissue_beta) ==
      is.na(paired_plasma_beta)
  )
)

##------------------------------------------------------------
## 13. CREATE REGION SUMMARY
##------------------------------------------------------------

message_header(
  "CREATING REGION-LEVEL SUMMARY"
)

retained_annotation <- copy(
  projected_annotation[
    retain_region
  ]
)

retained_annotation[
  ,
  n_valid_pairs :=
    n_valid_pairs[
      retain_region
    ]
]

retained_annotation[
  ,
  median_plasma_coverage :=
    safe_row_medians(
      matched_plasma_coverage[
        retain_region,
        ,
        drop = FALSE
      ]
    )
]

retained_annotation[
  ,
  median_plasma_covered_cpgs :=
    safe_row_medians(
      matched_plasma_cpgs[
        retain_region,
        ,
        drop = FALSE
      ]
    )
]

retained_annotation[
  ,
  mean_tissue_beta :=
    safe_row_means(
      paired_tissue_beta
    )
]

retained_annotation[
  ,
  mean_plasma_beta :=
    safe_row_means(
      paired_plasma_beta
    )
]

retained_annotation[
  ,
  median_tissue_beta :=
    safe_row_medians(
      paired_tissue_beta
    )
]

retained_annotation[
  ,
  median_plasma_beta :=
    safe_row_medians(
      paired_plasma_beta
    )
]

retained_annotation[
  ,
  mean_delta_beta :=
    safe_row_means(
      paired_delta_beta
    )
]

retained_annotation[
  ,
  median_delta_beta :=
    safe_row_medians(
      paired_delta_beta
    )
]

stopifnot(
  identical(
    retained_annotation$region_id,
    rownames(paired_tissue_beta)
  )
)

##------------------------------------------------------------
## 14. FILTER SENSITIVITY SUMMARY
##------------------------------------------------------------

message_header(
  "CREATING FILTER SENSITIVITY SUMMARY"
)

EPIC_probe_thresholds <- c(
  1L,
  2L,
  3L,
  5L
)

valid_pair_thresholds <- sort(
  unique(
    c(
      8L,
      9L,
      10L,
      11L,
      12L,
      n_pairs
    )
  )
)

valid_pair_thresholds <- valid_pair_thresholds[
  valid_pair_thresholds <= n_pairs
]

filter_sensitivity <- rbindlist(
  lapply(
    EPIC_probe_thresholds,
    function(probe_threshold) {
      rbindlist(
        lapply(
          valid_pair_thresholds,
          function(pair_threshold) {
            data.table(
              minimum_EPIC_probes =
                probe_threshold,
              minimum_valid_pairs =
                pair_threshold,
              regions_retained =
                sum(
                  projected_annotation$
                    n_EPIC_probes >=
                    probe_threshold &
                    n_valid_pairs >=
                    pair_threshold
                )
            )
          }
        )
      )
    }
  )
)

print(
  filter_sensitivity
)

##------------------------------------------------------------
## 15. CREATE FILTER SUMMARY
##------------------------------------------------------------

filter_summary <- data.table(
  item = c(
    "tissue_samples_input",
    "plasma_samples_input",
    "matched_patients",
    "projected_regions_input",
    "minimum_plasma_coverage",
    "minimum_plasma_covered_cpgs",
    "minimum_EPIC_probes",
    "minimum_valid_pairs",
    "regions_passing_EPIC_probe_threshold",
    "regions_passing_valid_pair_threshold",
    "regions_retained"
  ),
  value = c(
    nrow(tissue_metadata),
    nrow(plasma_metadata),
    n_pairs,
    nrow(projected_annotation),
    minimum_plasma_coverage,
    minimum_plasma_covered_cpgs,
    minimum_EPIC_probes,
    minimum_valid_pairs,
    sum(region_pass_EPIC),
    sum(region_pass_pairs),
    sum(retain_region)
  )
)

print(
  filter_summary
)

filter_settings <- data.table(
  setting = c(
    "delta_beta_definition",
    "minimum_plasma_coverage",
    "minimum_plasma_covered_cpgs",
    "minimum_EPIC_probes",
    "minimum_valid_pairs",
    "beta_epsilon",
    "tissue_aggregation_method",
    "genome_build"
  ),
  value = c(
    "plasma_minus_tissue",
    minimum_plasma_coverage,
    minimum_plasma_covered_cpgs,
    minimum_EPIC_probes,
    minimum_valid_pairs,
    beta_epsilon,
    "median_EPIC_beta_within_1kb_region",
    "hg38"
  )
)

##------------------------------------------------------------
## 16. SAVE OUTPUTS
##------------------------------------------------------------

message_header(
  "SAVING FILTERED PAIRED OBJECTS"
)

saveRDS(
  paired_tissue_beta,
  file.path(
    paired_dir,
    "SOLID_paired_tissue_beta_filtered.rds"
  )
)

saveRDS(
  paired_plasma_beta,
  file.path(
    paired_dir,
    "SOLID_paired_plasma_beta_filtered.rds"
  )
)

saveRDS(
  paired_delta_beta,
  file.path(
    paired_dir,
    "SOLID_paired_delta_beta_filtered.rds"
  )
)

saveRDS(
  paired_tissue_M,
  file.path(
    paired_dir,
    "SOLID_paired_tissue_M_filtered.rds"
  )
)

saveRDS(
  paired_plasma_M,
  file.path(
    paired_dir,
    "SOLID_paired_plasma_M_filtered.rds"
  )
)

saveRDS(
  paired_delta_M,
  file.path(
    paired_dir,
    "SOLID_paired_delta_M_filtered.rds"
  )
)

if (save_valid_pair_mask) {
  saveRDS(
    paired_valid_pair_mask,
    file.path(
      paired_dir,
      "SOLID_paired_valid_pair_mask.rds"
    )
  )
}

fwrite(
  retained_annotation,
  file.path(
    paired_dir,
    "SOLID_paired_retained_region_annotation.tsv.gz"
  ),
  sep = "\t",
  compress = "gzip"
)

fwrite(
  filter_summary,
  file.path(
    paired_dir,
    "SOLID_paired_filter_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  filter_sensitivity,
  file.path(
    paired_dir,
    "SOLID_paired_filter_sensitivity.tsv"
  ),
  sep = "\t"
)

fwrite(
  filter_settings,
  file.path(
    paired_dir,
    "SOLID_paired_filter_settings.tsv"
  ),
  sep = "\t"
)

matched_indices <- list(
  matched_metadata =
    matched_metadata,
  matched_patients =
    patient_order,
  tissue_sample_ids =
    matched_metadata$tissue_sample_id,
  plasma_sample_ids =
    matched_metadata$plasma_sample_id,
  tissue_column_indices =
    tissue_column_indices,
  plasma_column_indices =
    plasma_column_indices,
  projected_plasma_row_indices =
    plasma_row_indices,
  retained_region_ids =
    retained_annotation$region_id
)

saveRDS(
  matched_indices,
  file.path(
    paired_dir,
    "SOLID_matched_indices.rds"
  )
)

paired_filtered_object <- list(
  tissue_beta =
    paired_tissue_beta,
  plasma_beta =
    paired_plasma_beta,
  delta_beta =
    paired_delta_beta,
  tissue_M =
    paired_tissue_M,
  plasma_M =
    paired_plasma_M,
  delta_M =
    paired_delta_M,
  valid_pair_mask =
    paired_valid_pair_mask,
  retained_annotation =
    retained_annotation,
  matched_metadata =
    matched_metadata,
  settings =
    as.list(
      setNames(
        filter_settings$value,
        filter_settings$setting
      )
    )
)

saveRDS(
  paired_filtered_object,
  file.path(
    paired_dir,
    "SOLID_paired_filtered_object.rds"
  )
)

cat(
  "\nScript 04 completed successfully.\n"
)

cat(
  "Matched patients:",
  n_pairs,
  "\n"
)

cat(
  "Retained regions:",
  format(
    nrow(retained_annotation),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Outputs saved to:\n",
  paired_dir,
  "\n"
)

