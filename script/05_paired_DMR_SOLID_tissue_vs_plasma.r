##------------------------------------------------------------
## 05. PAIRED DMR ANALYSIS:
##     SOLID TISSUE EPIC VS PLASMA 5-BASE
##
## Statistical strategy:
##   - Input regions were filtered in script 04.
##   - For each matched patient:
##       delta M = plasma M - tissue M
##   - An intercept-only limma model tests whether mean delta M
##     differs from zero across matched patients.
##
## Effect-size reporting:
##   - median delta Beta = median(plasma Beta - tissue Beta)
##   - positive values: higher methylation in plasma
##   - negative values: higher methylation in tissue
##
## Candidate DMR definition:
##   - FDR < 0.05
##   - absolute median delta Beta >= 0.05
##
## Input:
##   result/tissue-plasma/paired_filtered/
##     SOLID_paired_filtered_object.rds
##
## Output:
##   result/tissue-plasma/paired_DMR/
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
minimum_valid_pairs <- 10L

fdr_threshold <- 0.05
delta_beta_threshold <- 0.05

use_robust_ebayes <- TRUE
use_trend_ebayes <- FALSE

##------------------------------------------------------------
## 2. REQUIRED PACKAGES
##------------------------------------------------------------

required_packages <- c(
  "data.table",
  "limma",
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
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall missing CRAN packages with install.packages(), ",
    "and Bioconductor packages with BiocManager::install()."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(matrixStats)
})

##------------------------------------------------------------
## 3. DIRECTORIES AND INPUT FILE
##------------------------------------------------------------

paired_dir <- file.path(
  project_dir,
  "result",
  "tissue-plasma",
  "paired_filtered"
)

dmr_dir <- file.path(
  project_dir,
  "result",
  "tissue-plasma",
  "paired_DMR"
)

dir.create(
  dmr_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

paired_object_file <- file.path(
  paired_dir,
  "SOLID_paired_filtered_object.rds"
)

if (!file.exists(paired_object_file)) {
  stop(
    "Input file not found:\n",
    paired_object_file,
    "\nRun script 04 first."
  )
}

##------------------------------------------------------------
## 4. HELPER FUNCTIONS
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
## 5. LOAD PAIRED FILTERED OBJECT
##------------------------------------------------------------

message_header(
  "LOADING FILTERED MATCHED TISSUE-PLASMA DATA"
)

paired_object <- readRDS(
  paired_object_file
)

required_elements <- c(
  "tissue_beta",
  "plasma_beta",
  "delta_beta",
  "tissue_M",
  "plasma_M",
  "delta_M",
  "valid_pair_mask",
  "retained_annotation",
  "matched_metadata"
)

missing_elements <- setdiff(
  required_elements,
  names(paired_object)
)

if (length(missing_elements) > 0L) {
  stop(
    "The paired object is missing required element(s): ",
    paste(missing_elements, collapse = ", ")
  )
}

tissue_beta <- paired_object$tissue_beta
plasma_beta <- paired_object$plasma_beta
delta_beta <- paired_object$delta_beta

tissue_M <- paired_object$tissue_M
plasma_M <- paired_object$plasma_M
delta_M <- paired_object$delta_M

valid_pair_mask <- paired_object$valid_pair_mask

region_annotation <- as.data.table(
  paired_object$retained_annotation
)

matched_metadata <- as.data.table(
  paired_object$matched_metadata
)

stopifnot(
  is.matrix(tissue_beta),
  is.matrix(plasma_beta),
  is.matrix(delta_beta),
  is.matrix(tissue_M),
  is.matrix(plasma_M),
  is.matrix(delta_M),

  nrow(delta_M) > 0L,
  ncol(delta_M) == expected_matched_patients,
  nrow(matched_metadata) == expected_matched_patients,

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
    dim(tissue_M)
  ),
  identical(
    dim(tissue_beta),
    dim(plasma_M)
  ),
  identical(
    dim(tissue_beta),
    dim(delta_M)
  ),

  identical(
    rownames(tissue_beta),
    rownames(delta_M)
  ),
  identical(
    colnames(tissue_beta),
    colnames(delta_M)
  ),

  nrow(region_annotation) ==
    nrow(delta_M),

  identical(
    as.character(region_annotation$region_id),
    rownames(delta_M)
  )
)

cat(
  "Regions loaded:",
  format(nrow(delta_M), big.mark = ","),
  "\n"
)

cat(
  "Matched patients:",
  ncol(delta_M),
  "\n"
)

##------------------------------------------------------------
## 6. FINAL VALIDITY CHECK
##------------------------------------------------------------

message_header(
  "CHECKING REGION-LEVEL VALID PAIR COUNTS"
)

n_valid_pairs <- rowSums(
  is.finite(delta_M) &
    is.finite(delta_beta)
)

keep_analysis <- n_valid_pairs >=
  minimum_valid_pairs

cat(
  "Regions with >=",
  minimum_valid_pairs,
  "valid pairs:",
  format(sum(keep_analysis), big.mark = ","),
  "of",
  format(length(keep_analysis), big.mark = ","),
  "\n"
)

if (!any(keep_analysis)) {
  stop(
    "No regions have enough valid matched pairs for analysis."
  )
}

delta_M_analysis <- delta_M[
  keep_analysis,
  ,
  drop = FALSE
]

delta_beta_analysis <- delta_beta[
  keep_analysis,
  ,
  drop = FALSE
]

tissue_beta_analysis <- tissue_beta[
  keep_analysis,
  ,
  drop = FALSE
]

plasma_beta_analysis <- plasma_beta[
  keep_analysis,
  ,
  drop = FALSE
]

analysis_annotation <- copy(
  region_annotation[
    keep_analysis
  ]
)

analysis_n_valid_pairs <- n_valid_pairs[
  keep_analysis
]

stopifnot(
  nrow(delta_M_analysis) > 0L,
  ncol(delta_M_analysis) ==
    expected_matched_patients,
  nrow(analysis_annotation) ==
    nrow(delta_M_analysis),
  identical(
    as.character(analysis_annotation$region_id),
    rownames(delta_M_analysis)
  )
)

##------------------------------------------------------------
## 7. RUN INTERCEPT-ONLY LIMMA MODEL
##------------------------------------------------------------

message_header(
  "RUNNING PAIRED LIMMA ANALYSIS ON DELTA M"
)

## Each column is one matched patient.
## The intercept tests whether the average within-patient delta M
## differs from zero.
design <- matrix(
  1,
  nrow = ncol(delta_M_analysis),
  ncol = 1,
  dimnames = list(
    colnames(delta_M_analysis),
    "Mean_Delta_M"
  )
)

fit <- limma::lmFit(
  delta_M_analysis,
  design = design
)

fit <- limma::eBayes(
  fit,
  robust = use_robust_ebayes,
  trend = use_trend_ebayes
)

limma_results <- limma::topTable(
  fit,
  coef = "Mean_Delta_M",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)

limma_results <- as.data.table(
  limma_results,
  keep.rownames = "region_id"
)

setnames(
  limma_results,
  old = c(
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B"
  ),
  new = c(
    "mean_delta_M",
    "average_delta_M",
    "moderated_t",
    "p_value",
    "FDR",
    "B_statistic"
  ),
  skip_absent = TRUE
)

stopifnot(
  nrow(limma_results) ==
    nrow(delta_M_analysis),
  !anyDuplicated(
    limma_results$region_id
  ),
  setequal(
    limma_results$region_id,
    rownames(delta_M_analysis)
  )
)

limma_results <- limma_results[
  match(
    rownames(delta_M_analysis),
    region_id
  )
]

stopifnot(
  identical(
    limma_results$region_id,
    rownames(delta_M_analysis)
  )
)

##------------------------------------------------------------
## 8. ADD BETA-SCALE EFFECT SIZES
##------------------------------------------------------------

message_header(
  "ADDING BETA-SCALE EFFECT SIZES"
)

effect_summary <- data.table(
  region_id =
    rownames(delta_beta_analysis),

  n_valid_pairs =
    analysis_n_valid_pairs,

  mean_tissue_beta =
    safe_row_means(
      tissue_beta_analysis
    ),

  mean_plasma_beta =
    safe_row_means(
      plasma_beta_analysis
    ),

  median_tissue_beta =
    safe_row_medians(
      tissue_beta_analysis
    ),

  median_plasma_beta =
    safe_row_medians(
      plasma_beta_analysis
    ),

  mean_delta_beta =
    safe_row_means(
      delta_beta_analysis
    ),

  median_delta_beta =
    safe_row_medians(
      delta_beta_analysis
    )
)

effect_summary[
  ,
  abs_median_delta_beta :=
    abs(median_delta_beta)
]

effect_summary[
  ,
  direction :=
    fifelse(
      median_delta_beta > 0,
      "Plasma_higher",
      fifelse(
        median_delta_beta < 0,
        "Tissue_higher",
        "No_change"
      )
    )
]

## Remove duplicated effect columns already carried over from script 04
effect_columns <- c(
  "n_valid_pairs",
  "mean_tissue_beta",
  "mean_plasma_beta",
  "median_tissue_beta",
  "median_plasma_beta",
  "mean_delta_beta",
  "median_delta_beta",
  "abs_median_delta_beta",
  "direction"
)

columns_to_remove <- intersect(
  effect_columns,
  names(analysis_annotation)
)

if (length(columns_to_remove) > 0L) {
  analysis_annotation[
    ,
    (columns_to_remove) := NULL
  ]
}

##------------------------------------------------------------
## 9. COMBINE ANNOTATION, EFFECTS, AND STATISTICS
##------------------------------------------------------------

message_header(
  "CREATING FINAL DMR RESULTS TABLE"
)

dmr_results <- merge(
  analysis_annotation,
  effect_summary,
  by = "region_id",
  all = FALSE,
  sort = FALSE
)

dmr_results <- merge(
  dmr_results,
  limma_results,
  by = "region_id",
  all = FALSE,
  sort = FALSE
)

## Restore original analysis-region order
dmr_results <- dmr_results[
  match(
    rownames(delta_M_analysis),
    region_id
  )
]

stopifnot(
  identical(
    dmr_results$region_id,
    rownames(delta_M_analysis)
  )
)

dmr_results[
  ,
  significance_class :=
    fifelse(
      FDR < fdr_threshold &
        median_delta_beta >=
          delta_beta_threshold,
      "Plasma_higher_DMR",
      fifelse(
        FDR < fdr_threshold &
          median_delta_beta <=
            -delta_beta_threshold,
        "Tissue_higher_DMR",
        "Not_significant"
      )
    )
]

dmr_results[
  ,
  minus_log10_FDR :=
    -log10(
      pmax(
        FDR,
        .Machine$double.xmin
      )
    )
]

## Rank strongest candidates first
setorder(
  dmr_results,
  FDR,
  -abs_median_delta_beta
)

##------------------------------------------------------------
## 10. CREATE CANDIDATE DMR TABLES
##------------------------------------------------------------

candidate_DMRs <- dmr_results[
  significance_class !=
    "Not_significant"
]

plasma_higher_DMRs <- candidate_DMRs[
  significance_class ==
    "Plasma_higher_DMR"
]

tissue_higher_DMRs <- candidate_DMRs[
  significance_class ==
    "Tissue_higher_DMR"
]

cat(
  "Candidate DMRs:",
  format(nrow(candidate_DMRs), big.mark = ","),
  "\n"
)

cat(
  "Plasma-higher DMRs:",
  format(nrow(plasma_higher_DMRs), big.mark = ","),
  "\n"
)

cat(
  "Tissue-higher DMRs:",
  format(nrow(tissue_higher_DMRs), big.mark = ","),
  "\n"
)

##------------------------------------------------------------
## 11. CREATE ANALYSIS SUMMARY
##------------------------------------------------------------

analysis_summary <- data.table(
  item = c(
    "regions_input_from_script_04",
    "regions_tested",
    "matched_patients",
    "minimum_valid_pairs",
    "FDR_threshold",
    "absolute_median_delta_beta_threshold",
    "candidate_DMRs",
    "plasma_higher_DMRs",
    "tissue_higher_DMRs"
  ),
  value = c(
    nrow(region_annotation),
    nrow(dmr_results),
    ncol(delta_M_analysis),
    minimum_valid_pairs,
    fdr_threshold,
    delta_beta_threshold,
    nrow(candidate_DMRs),
    nrow(plasma_higher_DMRs),
    nrow(tissue_higher_DMRs)
  )
)

print(
  analysis_summary
)

analysis_settings <- data.table(
  setting = c(
    "statistical_model",
    "response",
    "effect_definition",
    "multiple_testing",
    "FDR_threshold",
    "median_delta_beta_threshold",
    "minimum_valid_pairs",
    "robust_ebayes",
    "trend_ebayes"
  ),
  value = c(
    "intercept_only_limma_on_within_patient_difference",
    "delta_M",
    "plasma_minus_tissue",
    "Benjamini_Hochberg",
    fdr_threshold,
    delta_beta_threshold,
    minimum_valid_pairs,
    use_robust_ebayes,
    use_trend_ebayes
  )
)

##------------------------------------------------------------
## 11B. FINAL VALIDATION BEFORE SAVING
##------------------------------------------------------------

message_header(
  "FINAL VALIDATION BEFORE SAVING"
)

stopifnot(
  nrow(dmr_results) > 0L,
  nrow(dmr_results) ==
    nrow(delta_M_analysis),
  ncol(delta_M_analysis) ==
    expected_matched_patients,
  nrow(matched_metadata) ==
    expected_matched_patients,

  "region_id" %in%
    names(dmr_results),
  "mean_delta_M" %in%
    names(dmr_results),
  "median_delta_beta" %in%
    names(dmr_results),
  "p_value" %in%
    names(dmr_results),
  "FDR" %in%
    names(dmr_results),
  "significance_class" %in%
    names(dmr_results),

  !anyNA(dmr_results$region_id),
  !anyDuplicated(dmr_results$region_id),

  all(
    dmr_results$FDR >= 0 &
      dmr_results$FDR <= 1,
    na.rm = TRUE
  )
)

##------------------------------------------------------------
## 12. SAVE OUTPUTS
##------------------------------------------------------------

message_header(
  "SAVING PAIRED DMR RESULTS"
)

fwrite(
  dmr_results,
  file.path(
    dmr_dir,
    "SOLID_tissue_vs_plasma_all_DMR_results.tsv.gz"
  ),
  sep = "\t",
  compress = "gzip"
)

fwrite(
  candidate_DMRs,
  file.path(
    dmr_dir,
    "SOLID_tissue_vs_plasma_candidate_DMRs.tsv"
  ),
  sep = "\t"
)

fwrite(
  plasma_higher_DMRs,
  file.path(
    dmr_dir,
    "SOLID_tissue_vs_plasma_plasma_higher_DMRs.tsv"
  ),
  sep = "\t"
)

fwrite(
  tissue_higher_DMRs,
  file.path(
    dmr_dir,
    "SOLID_tissue_vs_plasma_tissue_higher_DMRs.tsv"
  ),
  sep = "\t"
)

fwrite(
  analysis_summary,
  file.path(
    dmr_dir,
    "SOLID_tissue_vs_plasma_DMR_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  analysis_settings,
  file.path(
    dmr_dir,
    "SOLID_tissue_vs_plasma_DMR_settings.tsv"
  ),
  sep = "\t"
)

saveRDS(
  fit,
  file.path(
    dmr_dir,
    "SOLID_tissue_vs_plasma_limma_fit.rds"
  )
)

dmr_object <- list(
  all_results =
    dmr_results,
  candidate_DMRs =
    candidate_DMRs,
  plasma_higher_DMRs =
    plasma_higher_DMRs,
  tissue_higher_DMRs =
    tissue_higher_DMRs,
  delta_M =
    delta_M_analysis,
  delta_beta =
    delta_beta_analysis,
  tissue_beta =
    tissue_beta_analysis,
  plasma_beta =
    plasma_beta_analysis,
  matched_metadata =
    matched_metadata,
  analysis_summary =
    analysis_summary,
  analysis_settings =
    analysis_settings
)

saveRDS(
  dmr_object,
  file.path(
    dmr_dir,
    "SOLID_tissue_vs_plasma_DMR_object.rds"
  )
)

## Reload saved object and verify
validation_object <- readRDS(
  file.path(
    dmr_dir,
    "SOLID_tissue_vs_plasma_DMR_object.rds"
  )
)

stopifnot(
  nrow(validation_object$all_results) > 0L,
  nrow(validation_object$all_results) ==
    nrow(delta_M_analysis),
  ncol(validation_object$delta_M) ==
    expected_matched_patients,
  nrow(validation_object$matched_metadata) ==
    expected_matched_patients
)

rm(validation_object)

cat(
  "\nScript 05 completed successfully.\n"
)

cat(
  "Regions tested:",
  format(nrow(dmr_results), big.mark = ","),
  "\n"
)

cat(
  "Candidate DMRs:",
  format(nrow(candidate_DMRs), big.mark = ","),
  "\n"
)

cat(
  "Plasma-higher DMRs:",
  format(nrow(plasma_higher_DMRs), big.mark = ","),
  "\n"
)

cat(
  "Tissue-higher DMRs:",
  format(nrow(tissue_higher_DMRs), big.mark = ","),
  "\n"
)

cat(
  "Outputs saved to:\n",
  dmr_dir,
  "\n"
)