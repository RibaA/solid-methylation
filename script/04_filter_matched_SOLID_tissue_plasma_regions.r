##------------------------------------------------------------
## 04. SOLID FILTER MATCHED TISSUE–PLASMA 1-KB REGIONS
##
## Purpose:
##   1. Load projected tumour EPIC 1-kb values from script 03
##   2. Load plasma 5-base Beta, coverage, and covered-CpG matrices
##   3. Match 13 SOLID tumour–plasma patient pairs
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

  has_id <- grepl(
    "SOLID-[0-9]+",
    x
  )

  out <- rep(
    NA_character_,
    length(x)
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

## Find a sample-level data.frame inside a list object.
find_metadata_table <- function(x) {
  if (is.data.frame(x)) {
    return(as.data.table(x))
  }

  if (is.list(x)) {
    for (i in seq_along(x)) {
      candidate <- find_metadata_table(x[[i]])

      if (!is.null(candidate)) {
        return(candidate)
      }
    }
  }

  NULL
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
    "SampleID"
  ),
  "tissue metadata"
)

tissue_patient_column <- detect_column(
  tissue_metadata,
  c(
    "Subject",
    "patient_id",
    "SolidID",
    "SOLID_ID"
  ))
