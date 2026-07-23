##------------------------------------------------------------
# CURATE PUBLIC HEALTHY PLASMA cfMeDIP-seq DATA
#
# Study:
# Nassiri et al., Nature Medicine 2020
# Detection and discrimination of intracranial tumors using
# plasma cell-free DNA methylomes
#
# Source:
# BrainTumourMeDIP_Reproducibility Zenodo archive
#
# Input:
# CombinedRPKM_Updated.RData
#
# Objects expected in the input file:
#   - Combined: processed cfMeDIP-seq RPKM matrix
#               rows = genomic windows
#               columns = plasma samples
#
#   - Input2: sample metadata
#             barcode, Class, Annotations
#
# Current objective:
#   1. Load the processed cfMeDIP data safely
#   2. Confirm matrix-to-metadata correspondence
#   3. Extract healthy plasma samples only
#   4. Parse genomic-window coordinates
#   5. Perform basic data integrity checks
#   6. Save a curated checkpoint
#
# Important:
# These healthy cfMeDIP samples will be used as an external
# healthy-plasma methylation reference. They should not be used
# as direct statistical controls for patient 5-base data.
##------------------------------------------------------------
##############################################################
# set up directory
##############################################################
dir_input <- 'data/raw'
dir_output <- 'data/proc'

##############################################################
# load data
##############################################################
cfmedip_env <- new.env(
  parent = emptyenv()
)

loaded_objects <- load(
  file.path(
  dir_input,
  "CombinedRPKM_Updated.RData"
  ),
  envir = cfmedip_env
)

message(
  "Objects loaded: ",
  paste(loaded_objects, collapse = ", ")
)

expected_objects <- c(
  "Combined",
  "Input2"
)

missing_objects <- setdiff(
  expected_objects,
  loaded_objects
)

if (length(missing_objects) > 0L) {
  stop(
    "The following expected objects are missing: ",
    paste(missing_objects, collapse = ", ")
  )
}

##############################################################
# extract data object
##############################################################
cfmedip_rpkm_all <- cfmedip_env$Combined
cfmedip_metadata_all <- cfmedip_env$Input2

rm(cfmedip_env)

##############################################################
# confirm object type 
##############################################################

if (!is.matrix(cfmedip_rpkm_all)) {
  stop(
    "`Combined` is expected to be a matrix, but its class is: ",
    paste(class(cfmedip_rpkm_all), collapse = ", ")
  )
}

if (!is.data.frame(cfmedip_metadata_all)) {
  stop(
    "`Input2` is expected to be a data.frame, but its class is: ",
    paste(class(cfmedip_metadata_all), collapse = ", ")
  )
}

required_metadata_columns <- c(
  "barcode",
  "Class",
  "Annotations"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(cfmedip_metadata_all)
)

if (length(missing_metadata_columns) > 0L) {
  stop(
    "Metadata columns are missing: ",
    paste(missing_metadata_columns, collapse = ", ")
  )
}

cat("\nFull processed cfMeDIP matrix:\n")
print(dim(cfmedip_rpkm_all))

cat("\nFull cfMeDIP metadata:\n")
print(dim(cfmedip_metadata_all))

##############################################################
# confirm matrix data and metadata
##############################################################
cfmedip_metadata_all$barcode <- trimws(
  as.character(cfmedip_metadata_all$barcode)
)

if (anyDuplicated(colnames(cfmedip_rpkm_all)) > 0L) {
  stop(
    "Duplicated sample IDs were found in the RPKM matrix."
  )
}

if (anyDuplicated(cfmedip_metadata_all$barcode) > 0L) {
  stop(
    "Duplicated sample IDs were found in the metadata."
  )
}

if (!setequal(
  colnames(cfmedip_rpkm_all),
  cfmedip_metadata_all$barcode
)) {
  stop(
    "Matrix sample IDs and metadata sample IDs do not match."
  )
}

if (!identical(
  colnames(cfmedip_rpkm_all),
  cfmedip_metadata_all$barcode
)) {
  message(
    "Metadata sample order differs from the matrix. ",
    "Reordering metadata to match matrix columns."
  )

  cfmedip_metadata_all <- cfmedip_metadata_all[
    match(
      colnames(cfmedip_rpkm_all),
      cfmedip_metadata_all$barcode
    ),
    ,
    drop = FALSE
  ]
}

stopifnot(
  identical(
    colnames(cfmedip_rpkm_all),
    cfmedip_metadata_all$barcode
  )
)

##############################################################
# review available samples
##############################################################

cat("\nSample counts by class:\n")

print(
  table(
    cfmedip_metadata_all$Class,
    useNA = "ifany"
  )
)

cat("\nSample counts by class and annotation:\n")

print(
  with(
    cfmedip_metadata_all,
    table(
      Class,
      Annotations,
      useNA = "ifany"
    )
  )
)

##############################################################
# extract healthy plasma samples
# The archive labels healthy plasma controls as:
# Class == "Normal"
##############################################################

healthy_index <- (
  trimws(cfmedip_metadata_all$Class) == "Normal"
)

if (!any(healthy_index)) {
  stop(
    "No samples with Class == 'Normal' were found."
  )
}

cfmedip_metadata_healthy <- cfmedip_metadata_all[
  healthy_index,
  ,
  drop = FALSE
]

cfmedip_rpkm_healthy <- cfmedip_rpkm_all[
  ,
  cfmedip_metadata_healthy$barcode,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(cfmedip_rpkm_healthy),
    cfmedip_metadata_healthy$barcode
  )
)

cat("\nHealthy cfMeDIP matrix dimensions:\n")
print(dim(cfmedip_rpkm_healthy))

cat("\nHealthy sample count:\n")
print(nrow(cfmedip_metadata_healthy))

if (ncol(cfmedip_rpkm_healthy) != 24L) {
  warning(
    "Expected 24 healthy samples, but found ",
    ncol(cfmedip_rpkm_healthy),
    "."
  )
}

##############################################################
# standardize healthy sample metadata
##############################################################

cfmedip_metadata_healthy <- data.frame(
  sample_id = cfmedip_metadata_healthy$barcode,
  assay = "cfMeDIP_seq",
  sample_type = "Healthy_plasma",
  disease_status = "Healthy",
  original_class = cfmedip_metadata_healthy$Class,
  original_annotation = cfmedip_metadata_healthy$Annotations,
  source_study = "Nassiri_2020",
  processed_measure = "RPKM",
  stringsAsFactors = FALSE
)

rownames(cfmedip_metadata_healthy) <-
  cfmedip_metadata_healthy$sample_id

##############################################################
# parse genomic window coordinates
# Example region ID:
# chr1.941701.942000
#
# Interpreted as:
# chromosome = chr1
# start      = 941701
# end        = 942000
##############################################################

region_ids <- rownames(cfmedip_rpkm_healthy)

if (is.null(region_ids)) {
  stop(
    "The processed cfMeDIP matrix has no row names."
  )
}

region_parts <- do.call(
  rbind,
  strsplit(
    region_ids,
    split = "\\."
  )
)

if (ncol(region_parts) != 3L) {
  stop(
    "Unexpected genomic-window row-name format."
  )
}

cfmedip_regions <- data.frame(
  region_id = region_ids,
  chr = region_parts[, 1],
  start = suppressWarnings(
    as.integer(region_parts[, 2])
  ),
  end = suppressWarnings(
    as.integer(region_parts[, 3])
  ),
  stringsAsFactors = FALSE
)

cfmedip_regions$width <- (
  cfmedip_regions$end -
    cfmedip_regions$start +
    1L
)

if (anyNA(cfmedip_regions$start) ||
    anyNA(cfmedip_regions$end)) {
  stop(
    "Some genomic coordinates could not be parsed."
  )
}

if (any(
  cfmedip_regions$start >
    cfmedip_regions$end
)) {
  stop(
    "Some genomic windows have start coordinates greater than end."
  )
}

if (anyDuplicated(cfmedip_regions$region_id) > 0L) {
  stop(
    "Duplicated genomic-window identifiers were found."
  )
}

stopifnot(
  identical(
    cfmedip_regions$region_id,
    rownames(cfmedip_rpkm_healthy)
  )
)

cat("\nGenomic-window width summary:\n")
print(
  summary(cfmedip_regions$width)
)

cat("\nMost frequent genomic-window widths:\n")
print(
  head(
    sort(
      table(cfmedip_regions$width),
      decreasing = TRUE
    ),
    10
  )
)

##############################################################
# QC
##############################################################

if (any(!is.finite(cfmedip_rpkm_healthy))) {
  stop(
    "The healthy RPKM matrix contains non-finite values."
  )
}

if (any(cfmedip_rpkm_healthy < 0, na.rm = TRUE)) {
  stop(
    "Negative values were found in the RPKM matrix."
  )
}

cfmedip_sample_qc <- data.frame(
  sample_id = colnames(cfmedip_rpkm_healthy),

  minimum_rpkm = apply(
    cfmedip_rpkm_healthy,
    2,
    min,
    na.rm = TRUE
  ),

  median_rpkm = apply(
    cfmedip_rpkm_healthy,
    2,
    median,
    na.rm = TRUE
  ),

  mean_rpkm = colMeans(
    cfmedip_rpkm_healthy,
    na.rm = TRUE
  ),

  maximum_rpkm = apply(
    cfmedip_rpkm_healthy,
    2,
    max,
    na.rm = TRUE
  ),

  zero_fraction = colMeans(
    cfmedip_rpkm_healthy == 0,
    na.rm = TRUE
  ),

  missing_fraction = colMeans(
    is.na(cfmedip_rpkm_healthy)
  ),

  stringsAsFactors = FALSE
)

cat("\nHealthy sample-level QC summary:\n")
print(
  summary(
    cfmedip_sample_qc[
      ,
      setdiff(
        names(cfmedip_sample_qc),
        "sample_id"
      ),
      drop = FALSE
    ]
  )
)

##############################################################
# create healthy reference
##############################################################

cfmedip_healthy_reference <- data.frame(
  region_id = rownames(cfmedip_rpkm_healthy),

  healthy_mean_rpkm = rowMeans(
    cfmedip_rpkm_healthy,
    na.rm = TRUE
  ),

  healthy_median_rpkm = apply(
    cfmedip_rpkm_healthy,
    1,
    median,
    na.rm = TRUE
  ),

  healthy_sd_rpkm = apply(
    cfmedip_rpkm_healthy,
    1,
    sd,
    na.rm = TRUE
  ),

  healthy_zero_fraction = rowMeans(
    cfmedip_rpkm_healthy == 0,
    na.rm = TRUE
  ),

  healthy_detected_fraction = rowMeans(
    cfmedip_rpkm_healthy > 0,
    na.rm = TRUE
  ),

  stringsAsFactors = FALSE
)

cfmedip_healthy_reference <- cbind(
  cfmedip_regions,
  cfmedip_healthy_reference[
    ,
    setdiff(
      names(cfmedip_healthy_reference),
      "region_id"
    ),
    drop = FALSE
  ]
)


##############################################################
# log-scale data
##############################################################

cfmedip_log2rpkm_healthy <- log2(
  cfmedip_rpkm_healthy + 1
)

##############################################################
# save curated data
##############################################################

file_checkpoint <- file.path(
  dir_output,
  "healthy_plasma_cfMeDIP_curated.RData"
)

save(
  cfmedip_rpkm_healthy,
  cfmedip_log2rpkm_healthy,
  cfmedip_metadata_healthy,
  cfmedip_regions,
  cfmedip_sample_qc,
  cfmedip_healthy_reference,
  file = file_checkpoint
)

message(
  "\nSaved curated checkpoint:\n",
  file_checkpoint
)

##############################################################
# EXPORT HUMAN-READABLE TABLES
##############################################################

write.csv(
  cfmedip_metadata_healthy,
  file.path(
    dir_output,
    "healthy_plasma_cfMeDIP_metadata.csv"
  ),
  row.names = FALSE
)

write.csv(
  cfmedip_sample_qc,
  file.path(
    dir_output,
    "healthy_plasma_cfMeDIP_sample_QC.csv"
  ),
  row.names = FALSE
)

write.csv(
  cfmedip_healthy_reference,
  file.path(
    dir_output,
    "healthy_plasma_cfMeDIP_region_reference.csv"
  ),
  row.names = FALSE
)

##############################################################
# FINAL SUMMARY
##############################################################

cat("\n========================================\n")
cat("HEALTHY cfMeDIP CURATION COMPLETE\n")
cat("========================================\n")

cat(
  "Processed measure: RPKM\n",
  "Healthy samples:   ",
  ncol(cfmedip_rpkm_healthy),
  "\n",
  "Genomic windows:   ",
  nrow(cfmedip_rpkm_healthy),
  "\n",
  "Matrix-metadata match: ",
  identical(
    colnames(cfmedip_rpkm_healthy),
    cfmedip_metadata_healthy$sample_id
  ),
  "\n",
  sep = ""
)

