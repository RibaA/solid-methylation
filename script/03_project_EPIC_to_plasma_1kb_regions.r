##------------------------------------------------------------
## 03. SOLID PROJECT TUMOUR EPIC METHYLATION
##     INTO PLASMA 1-KB REGIONS
##
## Purpose:
##   1. Load and validate harmonized tumour EPIC v1/v2 data
##   2. Obtain hg38 coordinates for harmonized canonical CpGs
##   3. Load and validate SOLID plasma 5-base 1-kb data
##   4. Map tumour EPIC CpGs into plasma 1-kb regions
##   5. Aggregate tumour Beta values within each plasma region
##   6. Calculate regional tumour M values from regional Beta
##   7. Save projected tissue-region objects for paired analysis
##
## Project structure used:
##
## C:/solid-methylation/
## ├── result/data/
## │   ├── beta_combined_EPICv1_EPICv2.rds
## │   ├── mval_combined_EPICv1_EPICv2.rds
## │   ├── harmonized_sample_metadata.csv
## │   └── common_filtered_EPICv1_EPICv2_CpGs.txt
## ├── data/solid-5base/
## │   ├── SOLID_1kb_beta_cov5_n7_standard_chr.rds
## │   ├── SOLID_1kb_M_cov5_n7_standard_chr.rds
## │   ├── SOLID_1kb_coverage_cov5_n7_standard_chr.rds
## │   ├── SOLID_1kb_covered_CpGs_cov5_n7_standard_chr.rds
## │   ├── SOLID_1kb_annotation_cov5_n7_standard_chr.tsv
## │   ├── SOLID_1kb_processing_metadata.rds
## │   └── SOLID_5base_sample_QC.tsv
## └── result/tissue-plasma/matched_preparation/
##
## Important:
##   - Tissue measurements originate from EPIC arrays.
##   - Plasma measurements originate from 5-base sequencing.
##   - Plasma 1-kb windows define the common genomic grid.
##   - Tissue regional Beta is the median Beta of retained EPIC
##     CpGs falling inside each plasma 1-kb region.
##   - This script performs projection only.
##   - Patient matching and paired filtering will occur in script 04.
##------------------------------------------------------------

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  warn = 1
)

##------------------------------------------------------------
## 1. USER SETTINGS AND DIRECTORIES
##------------------------------------------------------------

project_dir <- "C:/solid-methylation"

## Harmonized tumour EPIC outputs
tissue_dir <- file.path(
  project_dir,
  "result",
  "data"
)

## SOLID plasma 5-base outputs
plasma_dir <- file.path(
  project_dir,
  "data",
  "solid-5base"
)

## Output directories
result_dir <- file.path(
  project_dir,
  "result",
  "tissue-plasma"
)

projection_dir <- file.path(
  result_dir,
  "matched_preparation"
)

dir.create(
  projection_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(
  dir.exists(project_dir),
  dir.exists(tissue_dir),
  dir.exists(plasma_dir),
  dir.exists(projection_dir)
)

cat("Project directory:\n", project_dir, "\n\n")
cat("Tissue directory:\n", tissue_dir, "\n\n")
cat("Plasma directory:\n", plasma_dir, "\n\n")
cat("Projection output directory:\n", projection_dir, "\n\n")

## Save every plasma region containing at least this many EPIC probes.
## Stricter thresholds such as >=2 or >=3 will be assessed in script 04.
minimum_EPIC_probes_to_save <- 1L

## Small value used when converting regional Beta to M values
beta_epsilon <- 1e-6

## Number of projected regions handled in each aggregation chunk
chunk_size <- 100000L

##------------------------------------------------------------
## 2. REQUIRED PACKAGES
##------------------------------------------------------------

cran_packages <- c(
  "data.table",
  "matrixStats"
)

bioc_packages <- c(
  "minfi",
  "GenomicRanges",
  "IRanges",
  "S4Vectors",
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38"
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
  library(matrixStats)
  library(minfi)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
})

##------------------------------------------------------------
## 3. HELPER FUNCTIONS
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

make_region_id <- function(chr, start, end) {
  paste(chr, start, end, sep = ":")
}

normalize_chr <- function(chr) {
  chr <- as.character(chr)
  chr <- sub("^chr", "", chr, ignore.case = TRUE)
  paste0("chr", chr)
}

safe_beta_to_M <- function(beta, epsilon = 1e-6) {
  beta <- pmin(
    pmax(beta, epsilon),
    1 - epsilon
  )

  log2(
    beta /
      (1 - beta)
  )
}

detect_column <- function(
    dat,
    candidates,
    table_name,
    required = TRUE
) {
  matching_column <- candidates[
    candidates %in% names(dat)
  ]

  if (length(matching_column) == 0L) {
    if (required) {
      stop(
        "Could not identify the required column in ",
        table_name,
        ".\nTried: ",
        paste(candidates, collapse = ", "),
        "\nAvailable columns:\n",
        paste(names(dat), collapse = ", ")
      )
    }

    return(NULL)
  }

  matching_column[[1]]
}

##------------------------------------------------------------
## 4. DEFINE INPUT FILES
##------------------------------------------------------------

tissue_beta_file <- file.path(
  tissue_dir,
  "beta_combined_EPICv1_EPICv2.rds"
)

tissue_M_file <- file.path(
  tissue_dir,
  "mval_combined_EPICv1_EPICv2.rds"
)

tissue_metadata_file <- file.path(
  tissue_dir,
  "harmonized_sample_metadata.csv"
)

common_cpg_file <- file.path(
  tissue_dir,
  "common_filtered_EPICv1_EPICv2_CpGs.txt"
)

plasma_beta_file <- file.path(
  plasma_dir,
  "SOLID_1kb_beta_cov5_n7_standard_chr.rds"
)

plasma_M_file <- file.path(
  plasma_dir,
  "SOLID_1kb_M_cov5_n7_standard_chr.rds"
)

plasma_coverage_file <- file.path(
  plasma_dir,
  "SOLID_1kb_coverage_cov5_n7_standard_chr.rds"
)

plasma_cpg_file <- file.path(
  plasma_dir,
  "SOLID_1kb_covered_CpGs_cov5_n7_standard_chr.rds"
)

plasma_annotation_file <- file.path(
  plasma_dir,
  "SOLID_1kb_annotation_cov5_n7_standard_chr.tsv"
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
  tissue_M_file,
  tissue_metadata_file,
  common_cpg_file,
  plasma_beta_file,
  plasma_M_file,
  plasma_coverage_file,
  plasma_cpg_file,
  plasma_annotation_file,
  plasma_processing_metadata_file,
  plasma_qc_file
)

file_check <- data.table(
  file = basename(required_files),
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

cat("\nAll required tissue and plasma input files were found.\n")

##------------------------------------------------------------
## 5. LOAD AND VALIDATE HARMONIZED TISSUE EPIC DATA
##------------------------------------------------------------

message_header(
  "LOADING HARMONIZED TUMOUR EPIC DATA"
)

tissue_beta <- readRDS(
  tissue_beta_file
)

tissue_M <- readRDS(
  tissue_M_file
)

tissue_metadata <- fread(
  tissue_metadata_file
)

common_cpgs_table <- fread(
  common_cpg_file,
  header = FALSE
)

setnames(
  common_cpgs_table,
  names(common_cpgs_table)[1],
  "CpG"
)

stopifnot(
  is.matrix(tissue_beta),
  is.matrix(tissue_M),
  identical(
    dim(tissue_beta),
    dim(tissue_M)
  ),
  identical(
    rownames(tissue_beta),
    rownames(tissue_M)
  ),
  identical(
    colnames(tissue_beta),
    colnames(tissue_M)
  ),
  !is.null(rownames(tissue_beta)),
  !is.null(colnames(tissue_beta)),
  !anyDuplicated(rownames(tissue_beta)),
  !anyDuplicated(colnames(tissue_beta)),
  all(
    tissue_beta >= 0 &
      tissue_beta <= 1,
    na.rm = TRUE
  )
)

if (
  !setequal(
    common_cpgs_table$CpG,
    rownames(tissue_beta)
  )
) {
  stop(
    "The common-CpG file does not match the rows of the ",
    "harmonized tissue Beta matrix."
  )
}

common_cpgs_table <- common_cpgs_table[
  match(
    rownames(tissue_beta),
    CpG
  )
]

stopifnot(
  identical(
    common_cpgs_table$CpG,
    rownames(tissue_beta)
  )
)

tissue_sample_column <- detect_column(
  tissue_metadata,
  candidates = c(
    "Sample_Name",
    "sample_id",
    "SampleID",
    "sample"
  ),
  table_name = "tissue metadata"
)

tissue_patient_column <- detect_column(
  tissue_metadata,
  candidates = c(
    "Subject",
    "patient_id",
    "SolidID",
    "SOLID_ID",
    "subject_id"
  ),
  table_name = "tissue metadata"
)

if (
  !setequal(
    tissue_metadata[[tissue_sample_column]],
    colnames(tissue_beta)
  )
) {
  stop(
    "Tissue metadata sample IDs do not match ",
    "the tissue Beta matrix column names."
  )
}

tissue_metadata <- tissue_metadata[
  match(
    colnames(tissue_beta),
    get(tissue_sample_column)
  )
]

stopifnot(
  identical(
    tissue_metadata[[tissue_sample_column]],
    colnames(tissue_beta)
  )
)

cat(
  "Tissue Beta dimensions:",
  paste(dim(tissue_beta), collapse = " x "),
  "\n"
)

cat(
  "Tissue metadata rows:",
  nrow(tissue_metadata),
  "\n"
)

cat(
  "Unique tissue patients:",
  uniqueN(
    tissue_metadata[[tissue_patient_column]]
  ),
  "\n"
)

##------------------------------------------------------------
## 6. CREATE hg38 ANNOTATION FOR HARMONIZED CpGs
##------------------------------------------------------------

message_header(
  "CREATING hg38 ANNOTATION FOR HARMONIZED CpGs"
)

## Load EPIC v2 hg38 annotation
epicv2_manifest <- minfi::getAnnotation(
  IlluminaHumanMethylationEPICv2anno.20a1.hg38
)

## Convert DFrame annotation to data.table
epicv2_manifest_dt <- as.data.table(
  as.data.frame(epicv2_manifest),
  keep.rownames = "annotation_rowname"
)

cat(
  "EPIC v2 annotation dimensions:",
  paste(
    dim(epicv2_manifest_dt),
    collapse = " x "
  ),
  "\n"
)

cat("EPIC v2 annotation columns:\n")
print(names(epicv2_manifest_dt))

##------------------------------------------------------------
## 6A. IDENTIFY REQUIRED ANNOTATION COLUMNS
##------------------------------------------------------------

chr_column <- intersect(
  c(
    "chr",
    "CHR",
    "Chromosome",
    "chromosome"
  ),
  names(epicv2_manifest_dt)
)

position_column <- intersect(
  c(
    "pos",
    "MAPINFO",
    "position",
    "Position"
  ),
  names(epicv2_manifest_dt)
)

if (length(chr_column) == 0L) {
  stop(
    "Could not identify chromosome column.\n",
    "Available columns:\n",
    paste(
      names(epicv2_manifest_dt),
      collapse = ", "
    )
  )
}

if (length(position_column) == 0L) {
  stop(
    "Could not identify genomic-position column.\n",
    "Available columns:\n",
    paste(
      names(epicv2_manifest_dt),
      collapse = ", "
    )
  )
}

chr_column <- chr_column[[1]]
position_column <- position_column[[1]]

##------------------------------------------------------------
## 6B. CANONICALIZE EPIC v2 PROBE IDENTIFIERS
##------------------------------------------------------------

## EPIC v2 may contain replicated probe IDs such as:
## cg00000029_BC11
## cg00000029_TC11
##
## The harmonized tissue matrix uses canonical IDs:
## cg00000029

extract_canonical_CpG <- function(x) {

  x <- as.character(x)

  canonical_id <- sub(
    "^(cg[0-9]+).*$",
    "\\1",
    x
  )

  canonical_id[
    !grepl(
      "^cg[0-9]+$",
      canonical_id
    )
  ] <- NA_character_

  canonical_id
}

## Use Name when present; otherwise use row names
raw_probe_id <- if (
  "Name" %in% names(epicv2_manifest_dt)
) {
  as.character(
    epicv2_manifest_dt$Name
  )
} else {
  as.character(
    epicv2_manifest_dt$annotation_rowname
  )
}

cat("First raw EPIC v2 probe IDs:\n")
print(head(raw_probe_id, 10))

cat("First harmonized tissue CpG IDs:\n")
print(head(rownames(tissue_beta), 10))

##------------------------------------------------------------
## 6C. BUILD CANONICAL CpG COORDINATE TABLE
##------------------------------------------------------------

epicv2_annotation <- data.table(
  raw_probe_id = raw_probe_id,
  CpG = extract_canonical_CpG(
    raw_probe_id
  ),
  chr = normalize_chr(
    epicv2_manifest_dt[[chr_column]]
  ),
  position = as.integer(
    epicv2_manifest_dt[[position_column]]
  )
)

standard_chromosomes <- paste0(
  "chr",
  c(
    1:22,
    "X",
    "Y"
  )
)

epicv2_annotation <- epicv2_annotation[
  !is.na(CpG) &
    CpG != "" &
    !is.na(chr) &
    chr %in% standard_chromosomes &
    !is.na(position)
]

cat(
  "EPIC v2 annotation rows after canonicalization:",
  format(
    nrow(epicv2_annotation),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Unique canonical CpGs:",
  format(
    uniqueN(epicv2_annotation$CpG),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Tissue CpGs matching EPIC v2 annotation:",
  format(
    sum(
      rownames(tissue_beta) %in%
        epicv2_annotation$CpG
    ),
    big.mark = ","
  ),
  "of",
  format(
    nrow(tissue_beta),
    big.mark = ","
  ),
  "\n"
)

##------------------------------------------------------------
## 6D. RESOLVE REPLICATED EPIC v2 PROBES
##------------------------------------------------------------

coordinate_conflicts <- epicv2_annotation[
  ,
  .(
    n_probe_designs = .N,
    n_coordinates = uniqueN(
      paste(
        chr,
        position,
        sep = ":"
      )
    )
  ),
  by = CpG
][
  n_coordinates > 1L
]

cat(
  "Canonical CpGs with multiple hg38 coordinates:",
  format(
    nrow(coordinate_conflicts),
    big.mark = ","
  ),
  "\n"
)

if (nrow(coordinate_conflicts) > 0L) {
  warning(
    nrow(coordinate_conflicts),
    " canonical CpGs have multiple hg38 coordinates. ",
    "The first coordinate will be retained."
  )
}

setorder(
  epicv2_annotation,
  CpG,
  chr,
  position,
  raw_probe_id
)

## Retain one coordinate per canonical CpG
epicv2_annotation <- epicv2_annotation[
  !duplicated(CpG)
]

stopifnot(
  !anyNA(epicv2_annotation$CpG),
  !anyNA(epicv2_annotation$chr),
  !anyNA(epicv2_annotation$position),
  !anyDuplicated(epicv2_annotation$CpG)
)

##------------------------------------------------------------
## 6E. MATCH HARMONIZED TISSUE CpGs
##------------------------------------------------------------

annotation_match <- match(
  rownames(tissue_beta),
  epicv2_annotation$CpG
)

tissue_cpg_annotation <- epicv2_annotation[
  annotation_match
]

missing_annotation <- which(
  is.na(annotation_match)
)

cat(
  "Harmonized CpGs before annotation filtering:",
  format(
    nrow(tissue_beta),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Harmonized CpGs matching hg38 annotation:",
  format(
    sum(!is.na(annotation_match)),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Harmonized CpGs missing hg38 coordinates:",
  format(
    length(missing_annotation),
    big.mark = ","
  ),
  "\n"
)

keep_annotated_cpgs <- which(
  !is.na(annotation_match)
)

if (length(keep_annotated_cpgs) == 0L) {
  stop(
    "No harmonized tissue CpGs matched after ",
    "canonicalizing EPIC v2 probe identifiers."
  )
}

##------------------------------------------------------------
## 6F. FILTER TISSUE MATRICES TO ANNOTATED CpGs
##------------------------------------------------------------

tissue_beta <- tissue_beta[
  keep_annotated_cpgs,
  ,
  drop = FALSE
]

tissue_M <- tissue_M[
  keep_annotated_cpgs,
  ,
  drop = FALSE
]

tissue_cpg_annotation <- tissue_cpg_annotation[
  keep_annotated_cpgs
]

stopifnot(
  nrow(tissue_cpg_annotation) ==
    nrow(tissue_beta),
  identical(
    tissue_cpg_annotation$CpG,
    rownames(tissue_beta)
  ),
  identical(
    rownames(tissue_beta),
    rownames(tissue_M)
  )
)

##------------------------------------------------------------
## 6G. SAVE hg38 CpG ANNOTATION
##------------------------------------------------------------

fwrite(
  tissue_cpg_annotation,
  file.path(
    projection_dir,
    "SOLID_tissue_EPIC_CpG_annotation_hg38.tsv.gz"
  ),
  sep = "\t",
  compress = "gzip"
)

cat(
  "Harmonized CpGs retained with hg38 coordinates:",
  format(
    nrow(tissue_cpg_annotation),
    big.mark = ","
  ),
  "\n"
)

##------------------------------------------------------------
## 7. LOAD AND VALIDATE SOLID PLASMA DATA
##------------------------------------------------------------

message_header(
  "LOADING SOLID PLASMA 1-KB DATA"
)

plasma_beta <- readRDS(
  plasma_beta_file
)

plasma_M <- readRDS(
  plasma_M_file
)

plasma_coverage <- readRDS(
  plasma_coverage_file
)

plasma_covered_cpgs <- readRDS(
  plasma_cpg_file
)

plasma_annotation <- fread(
  plasma_annotation_file
)

plasma_processing_metadata <- readRDS(
  plasma_processing_metadata_file
)

plasma_qc <- fread(
  plasma_qc_file
)

stopifnot(
  is.matrix(plasma_beta),
  is.matrix(plasma_M),
  is.matrix(plasma_coverage),
  is.matrix(plasma_covered_cpgs),
  identical(
    dim(plasma_beta),
    dim(plasma_M)
  ),
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
    colnames(plasma_M)
  ),
  identical(
    colnames(plasma_beta),
    colnames(plasma_coverage)
  ),
  identical(
    colnames(plasma_beta),
    colnames(plasma_covered_cpgs)
  ),
  nrow(plasma_annotation) ==
    nrow(plasma_beta),
  all(
    plasma_beta >= 0 &
      plasma_beta <= 1,
    na.rm = TRUE
  )
)

if (is.data.frame(plasma_processing_metadata)) {
  plasma_processing_metadata <- as.data.table(
    plasma_processing_metadata
  )
}

cat(
  "Plasma Beta dimensions:",
  paste(dim(plasma_beta), collapse = " x "),
  "\n"
)

cat(
  "Plasma annotation rows:",
  format(
    nrow(plasma_annotation),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Plasma processing metadata class:",
  paste(
    class(plasma_processing_metadata),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Plasma QC dimensions:",
  paste(dim(plasma_qc), collapse = " x "),
  "\n"
)

##------------------------------------------------------------
## 8. STANDARDIZE PLASMA REGION ANNOTATION
##------------------------------------------------------------

message_header(
  "STANDARDIZING PLASMA REGION ANNOTATION"
)

required_plasma_annotation_columns <- c(
  "chr",
  "start",
  "end"
)

missing_plasma_columns <- setdiff(
  required_plasma_annotation_columns,
  names(plasma_annotation)
)

if (length(missing_plasma_columns) > 0L) {
  stop(
    "Plasma annotation is missing required columns: ",
    paste(
      missing_plasma_columns,
      collapse = ", "
    )
  )
}

plasma_annotation[
  ,
  chr := normalize_chr(chr)
]

plasma_annotation[
  ,
  start := as.integer(start)
]

plasma_annotation[
  ,
  end := as.integer(end)
]

if (!"region_id" %in% names(plasma_annotation)) {
  plasma_annotation[
    ,
    region_id := make_region_id(
      chr,
      start,
      end
    )
  ]
} else {
  plasma_annotation[
    ,
    region_id := as.character(region_id)
  ]
}

if (anyDuplicated(plasma_annotation$region_id)) {
  stop(
    "Plasma annotation contains duplicated region IDs."
  )
}

if (!is.null(rownames(plasma_beta))) {
  if (
    setequal(
      rownames(plasma_beta),
      plasma_annotation$region_id
    )
  ) {
    plasma_annotation <- plasma_annotation[
      match(
        rownames(plasma_beta),
        region_id
      )
    ]
  } else {
    warning(
      "Plasma matrix row names do not match region_id. ",
      "Assuming annotation rows are already in matrix row order."
    )
  }
}

stopifnot(
  nrow(plasma_annotation) ==
    nrow(plasma_beta)
)

##------------------------------------------------------------
## 9. BUILD GENOMIC RANGES AND FIND OVERLAPS
##------------------------------------------------------------

message_header(
  "MAPPING EPIC CpGs TO PLASMA 1-KB REGIONS"
)

tissue_cpg_gr <- GRanges(
  seqnames = tissue_cpg_annotation$chr,
  ranges = IRanges(
    start = tissue_cpg_annotation$position,
    end = tissue_cpg_annotation$position
  ),
  CpG = tissue_cpg_annotation$CpG,
  tissue_row_index = seq_len(
    nrow(tissue_cpg_annotation)
  )
)

plasma_region_gr <- GRanges(
  seqnames = plasma_annotation$chr,
  ranges = IRanges(
    start = plasma_annotation$start,
    end = plasma_annotation$end
  ),
  region_id = plasma_annotation$region_id,
  plasma_row_index = seq_len(
    nrow(plasma_annotation)
  )
)

overlap_hits <- findOverlaps(
  query = plasma_region_gr,
  subject = tissue_cpg_gr,
  ignore.strand = TRUE
)

if (length(overlap_hits) == 0L) {
  stop(
    "No overlaps were found between tissue EPIC CpGs ",
    "and plasma 1-kb regions. Confirm that both datasets ",
    "use hg38 coordinates and compatible chromosome names."
  )
}

mapping_table <- data.table(
  plasma_row_index = queryHits(
    overlap_hits
  ),
  tissue_row_index = subjectHits(
    overlap_hits
  )
)

mapping_table[
  ,
  region_id :=
    plasma_annotation$region_id[
      plasma_row_index
    ]
]

mapping_table[
  ,
  CpG :=
    tissue_cpg_annotation$CpG[
      tissue_row_index
    ]
]

mapping_table[
  ,
  chr :=
    plasma_annotation$chr[
      plasma_row_index
    ]
]

mapping_table[
  ,
  region_start :=
    plasma_annotation$start[
      plasma_row_index
    ]
]

mapping_table[
  ,
  region_end :=
    plasma_annotation$end[
      plasma_row_index
    ]
]

setorder(
  mapping_table,
  plasma_row_index,
  tissue_row_index
)

fwrite(
  mapping_table,
  file.path(
    projection_dir,
    "SOLID_EPIC_to_plasma_1kb_mapping.tsv.gz"
  ),
  sep = "\t",
  compress = "gzip"
)

cat(
  "CpG-to-region overlap records:",
  format(
    nrow(mapping_table),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Plasma regions containing at least one EPIC CpG:",
  format(
    uniqueN(mapping_table$plasma_row_index),
    big.mark = ","
  ),
  "\n"
)

##------------------------------------------------------------
## 10. CREATE REGIONAL EPIC PROBE COUNTS
##------------------------------------------------------------

message_header(
  "CALCULATING EPIC PROBE COUNTS PER REGION"
)

region_probe_counts <- mapping_table[
  ,
  .(
    n_EPIC_probes = uniqueN(CpG)
  ),
  by = .(
    plasma_row_index,
    region_id,
    chr,
    region_start,
    region_end
  )
]

region_probe_counts <- region_probe_counts[
  n_EPIC_probes >=
    minimum_EPIC_probes_to_save
]

setorder(
  region_probe_counts,
  plasma_row_index
)

retained_plasma_rows <- region_probe_counts$
  plasma_row_index

n_projected_regions <- length(
  retained_plasma_rows
)

if (n_projected_regions == 0L) {
  stop(
    "No regions passed the minimum EPIC-probe requirement."
  )
}

cat(
  "Regions retained for tissue projection:",
  format(
    n_projected_regions,
    big.mark = ","
  ),
  "\n"
)

##------------------------------------------------------------
## 11. AGGREGATE TISSUE EPIC BETA BY 1-KB REGION
##------------------------------------------------------------

message_header(
  "AGGREGATING TISSUE EPIC BETA VALUES"
)

## Initialize output matrix:
## rows = projected plasma 1-kb regions
## columns = tissue EPIC samples
tissue_region_beta <- matrix(
  NA_real_,
  nrow = n_projected_regions,
  ncol = ncol(tissue_beta),
  dimnames = list(
    region_probe_counts$region_id,
    colnames(tissue_beta)
  )
)

## Map original plasma row indices to output matrix row indices
output_row_lookup <- setNames(
  seq_len(n_projected_regions),
  as.character(retained_plasma_rows)
)

## Define chunks of projected regions
chunk_starts <- seq(
  from = 1L,
  to = n_projected_regions,
  by = chunk_size
)

for (chunk_number in seq_along(chunk_starts)) {

  start_index <- chunk_starts[chunk_number]

  end_index <- min(
    start_index + chunk_size - 1L,
    n_projected_regions
  )

  output_rows <- start_index:end_index

  plasma_rows_chunk <- retained_plasma_rows[
    output_rows
  ]

  cat(
    sprintf(
      "[%03d/%03d] Projected regions %s-%s\n",
      chunk_number,
      length(chunk_starts),
      format(start_index, big.mark = ","),
      format(end_index, big.mark = ",")
    )
  )

  chunk_mapping <- mapping_table[
    plasma_row_index %in% plasma_rows_chunk
  ]

  split_tissue_rows <- split(
    chunk_mapping$tissue_row_index,
    chunk_mapping$plasma_row_index
  )

  region_names <- names(split_tissue_rows)

  for (i in seq_along(split_tissue_rows)) {

    plasma_row_character <- region_names[i]

    output_row <- unname(
      output_row_lookup[plasma_row_character]
    )

    cpg_rows <- unique(
      split_tissue_rows[[i]]
    )

    if (
      length(output_row) != 1L ||
        is.na(output_row)
    ) {
      stop(
        "Could not identify output row for plasma row: ",
        plasma_row_character
      )
    }

    region_values <- tissue_beta[
      cpg_rows,
      ,
      drop = FALSE
    ]

    if (length(cpg_rows) == 1L) {

      tissue_region_beta[
        output_row,
      ] <- as.numeric(
        region_values[1, ]
      )

    } else {

      tissue_region_beta[
        output_row,
      ] <- matrixStats::colMedians(
        region_values,
        na.rm = TRUE
      )
    }
  }
}

tissue_region_beta[
  !is.finite(tissue_region_beta)
] <- NA_real_

stopifnot(
  is.matrix(tissue_region_beta),
  nrow(tissue_region_beta) == n_projected_regions,
  ncol(tissue_region_beta) == ncol(tissue_beta),
  identical(
    rownames(tissue_region_beta),
    region_probe_counts$region_id
  ),
  identical(
    colnames(tissue_region_beta),
    colnames(tissue_beta)
  ),
  all(
    tissue_region_beta >= 0 &
      tissue_region_beta <= 1,
    na.rm = TRUE
  )
)

cat(
  "Projected tissue Beta dimensions:",
  paste(dim(tissue_region_beta), collapse = " x "),
  "\n"
)
##------------------------------------------------------------
## 12. CREATE TISSUE REGIONAL M VALUES
##------------------------------------------------------------

message_header(
  "CREATING TISSUE REGIONAL M VALUES"
)

tissue_region_M <- safe_beta_to_M(
  tissue_region_beta,
  epsilon = beta_epsilon
)

rownames(tissue_region_M) <- rownames(
  tissue_region_beta
)

colnames(tissue_region_M) <- colnames(
  tissue_region_beta
)

##------------------------------------------------------------
## 13. CREATE PROJECTED REGION ANNOTATION
##------------------------------------------------------------

message_header(
  "CREATING PROJECTED REGION ANNOTATION"
)

projected_region_annotation <- copy(
  region_probe_counts
)

setnames(
  projected_region_annotation,
  old = c(
    "region_start",
    "region_end"
  ),
  new = c(
    "start",
    "end"
  )
)

projected_region_annotation[
  ,
  n_tissue_samples_observed :=
    rowSums(
      !is.na(tissue_region_beta)
    )
]

projected_region_annotation[
  ,
  tissue_samples_observed_pct :=
    100 *
      n_tissue_samples_observed /
      ncol(tissue_region_beta)
]

setcolorder(
  projected_region_annotation,
  c(
    "region_id",
    "chr",
    "start",
    "end",
    "plasma_row_index",
    "n_EPIC_probes",
    "n_tissue_samples_observed",
    "tissue_samples_observed_pct"
  )
)

stopifnot(
  identical(
    projected_region_annotation$region_id,
    rownames(tissue_region_beta)
  )
)

##------------------------------------------------------------
## 14. VALIDATE PROJECTED MATRICES
##------------------------------------------------------------

message_header(
  "VALIDATING PROJECTED TISSUE MATRICES"
)

stopifnot(
  is.matrix(tissue_region_beta),
  is.matrix(tissue_region_M),
  identical(
    dim(tissue_region_beta),
    dim(tissue_region_M)
  ),
  identical(
    rownames(tissue_region_beta),
    rownames(tissue_region_M)
  ),
  identical(
    colnames(tissue_region_beta),
    colnames(tissue_region_M)
  ),
  nrow(tissue_region_beta) ==
    nrow(projected_region_annotation),
  identical(
    rownames(tissue_region_beta),
    projected_region_annotation$region_id
  )
)

cat(
  "Projected tissue Beta dimensions:",
  paste(
    dim(tissue_region_beta),
    collapse = " x "
  ),
  "\n"
)

cat(
  "Projected tissue M dimensions:",
  paste(
    dim(tissue_region_M),
    collapse = " x "
  ),
  "\n"
)

cat(
  "Median EPIC probes per projected region:",
  median(
    projected_region_annotation$
      n_EPIC_probes
  ),
  "\n"
)

cat(
  "Regions with at least 3 EPIC probes:",
  format(
    sum(
      projected_region_annotation$
        n_EPIC_probes >= 3L
    ),
    big.mark = ","
  ),
  "\n"
)

##------------------------------------------------------------
## 15. CREATE PROJECTION SUMMARY
##------------------------------------------------------------

projection_summary <- data.table(
  item = c(
    "harmonized_tissue_CpGs_input",
    "harmonized_tissue_CpGs_with_hg38_coordinates",
    "plasma_1kb_regions_input",
    "CpG_region_overlap_records",
    "plasma_regions_with_EPIC_probes",
    "regions_saved",
    "regions_with_at_least_2_EPIC_probes",
    "regions_with_at_least_3_EPIC_probes",
    "regions_with_at_least_5_EPIC_probes",
    "tissue_samples",
    "plasma_samples"
  ),
  value = c(
    nrow(common_cpgs_table),
    nrow(tissue_cpg_annotation),
    nrow(plasma_annotation),
    nrow(mapping_table),
    uniqueN(mapping_table$plasma_row_index),
    nrow(projected_region_annotation),
    sum(
      projected_region_annotation$
        n_EPIC_probes >= 2L
    ),
    sum(
      projected_region_annotation$
        n_EPIC_probes >= 3L
    ),
    sum(
      projected_region_annotation$
        n_EPIC_probes >= 5L
    ),
    ncol(tissue_region_beta),
    ncol(plasma_beta)
  )
)

print(
  projection_summary
)

##------------------------------------------------------------
## 16. SAVE OUTPUTS
##------------------------------------------------------------

message_header(
  "SAVING PROJECTED TISSUE OBJECTS"
)

saveRDS(
  tissue_region_beta,
  file.path(
    projection_dir,
    "SOLID_tissue_1kb_beta.rds"
  )
)

saveRDS(
  tissue_region_M,
  file.path(
    projection_dir,
    "SOLID_tissue_1kb_M.rds"
  )
)

fwrite(
  projected_region_annotation,
  file.path(
    projection_dir,
    "SOLID_tissue_1kb_region_annotation.tsv.gz"
  ),
  sep = "\t",
  compress = "gzip"
)

fwrite(
  projection_summary,
  file.path(
    projection_dir,
    "SOLID_tissue_1kb_projection_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  tissue_metadata,
  file.path(
    projection_dir,
    "SOLID_tissue_harmonized_metadata.tsv"
  ),
  sep = "\t"
)

projection_object <- list(
  tissue_region_beta =
    tissue_region_beta,
  tissue_region_M =
    tissue_region_M,
  projected_region_annotation =
    projected_region_annotation,
  tissue_cpg_annotation =
    tissue_cpg_annotation,
  tissue_metadata =
    tissue_metadata,
  plasma_annotation =
    plasma_annotation,
  plasma_sample_ids =
    colnames(plasma_beta),
  plasma_processing_metadata =
    plasma_processing_metadata,
  plasma_qc =
    plasma_qc,
  settings = list(
    project_dir =
      project_dir,
    tissue_dir =
      tissue_dir,
    plasma_dir =
      plasma_dir,
    projection_dir =
      projection_dir,
    minimum_EPIC_probes_to_save =
      minimum_EPIC_probes_to_save,
    aggregation_method =
      "median",
    genome_build =
      "hg38",
    beta_epsilon =
      beta_epsilon,
    chunk_size =
      chunk_size
  )
)

saveRDS(
  projection_object,
  file.path(
    projection_dir,
    "SOLID_tissue_1kb_projection_object.rds"
  )
)

cat(
  "\nScript 03 completed successfully.\n"
)

cat(
  "Outputs saved to:\n",
  projection_dir,
  "\n"
)
