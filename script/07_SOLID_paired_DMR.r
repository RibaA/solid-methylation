############################################################
## 07_SOLID_paired_DMR.r
##
## SOLID MATCHED TUMOR–PLASMA PAIRED DMR ANALYSIS
##
## Biological comparison:
##   Matched tumor EPIC vs plasma 5-base methylation
##
## Statistical strategy:
##   - Input regions were frozen in Script 02
##   - Downstream object and clinical metadata validated
##     in Script 03
##   - For each matched patient:
##
##       Delta M = Plasma M - Tumor M
##
##   - Intercept-only limma model tests whether mean
##     within-patient Delta M differs from zero
##
## Biological effect size:
##
##       Delta Beta = Plasma Beta - Tumor Beta
##
##   Positive Delta Beta:
##       Plasma higher methylation
##
##   Negative Delta Beta:
##       Tumor higher methylation
##
## Primary candidate DMR definition:
##
##   FDR < 0.05
##   AND
##   |median Delta Beta| >= 0.05
##
## Additional effect-size sensitivity:
##
##   |median Delta Beta| >= 0.05
##   |median Delta Beta| >= 0.10
##   |median Delta Beta| >= 0.20
##
## IMPORTANT:
##   - Existing valid_pair_mask is explicitly respected
##   - Minimum valid matched pairs = 10
##   - No additional sample or region QC is performed here
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
## 1. USER SETTINGS
############################################################

project_dir <- "C:/solid-methylation"

expected_matched_patients <- 13L

minimum_valid_pairs <- 10L

fdr_threshold <- 0.05

primary_delta_beta_threshold <- 0.05

delta_beta_sensitivity_thresholds <- c(
  0.05,
  0.10,
  0.20
)

use_robust_ebayes <- TRUE

use_trend_ebayes <- FALSE


############################################################
## 2. REQUIRED PACKAGES
############################################################

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
    "Missing package(s): ",
    paste(
      missing_packages,
      collapse = ", "
    ),
    "\nInstall CRAN packages with install.packages() ",
    "and Bioconductor packages with BiocManager::install()."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(matrixStats)
})


############################################################
## 3. DIRECTORIES
############################################################

input_file <- file.path(
  project_dir,
  "result",
  "03_matched_tissue_plasma",
  "SOLID_downstream_analysis_ready.rds"
)


dmr_dir <- file.path(
  project_dir,
  "result",
  "03_matched_tissue_plasma",
  "DMR"
)


dir.create(
  dmr_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


if (!file.exists(input_file)) {

  stop(
    "Analysis-ready input file not found:\n",
    input_file,
    "\nRun Script 03 first."
  )
}


############################################################
## 4. HELPER FUNCTIONS
############################################################

message_header <- function(text) {

  cat(
    "\n",
    paste(
      rep(
        "=",
        72
      ),
      collapse = ""
    ),
    "\n",
    text,
    "\n",
    paste(
      rep(
        "=",
        72
      ),
      collapse = ""
    ),
    "\n",
    sep = ""
  )
}


safe_row_means <- function(x) {

  out <- rowMeans(
    x,
    na.rm = TRUE
  )

  out[
    !is.finite(out)
  ] <- NA_real_

  out
}


safe_row_medians <- function(x) {

  out <- matrixStats::rowMedians(
    x,
    na.rm = TRUE
  )

  out[
    !is.finite(out)
  ] <- NA_real_

  out
}


############################################################
## 5. LOAD ANALYSIS-READY OBJECT
############################################################

message_header(
  "LOADING SOLID ANALYSIS-READY OBJECT"
)


obj <- readRDS(
  input_file
)


required_elements <- c(
  "tissue_beta",
  "plasma_beta",
  "delta_beta",
  "tissue_M",
  "plasma_M",
  "delta_M",
  "valid_pair_mask",
  "region_annotation",
  "patient_metadata",
  "matched_metadata",
  "settings"
)


missing_elements <- setdiff(
  required_elements,
  names(obj)
)


if (length(missing_elements) > 0L) {

  stop(
    "Analysis-ready object missing element(s): ",
    paste(
      missing_elements,
      collapse = ", "
    )
  )
}


tissue_beta <- obj$tissue_beta

plasma_beta <- obj$plasma_beta

delta_beta <- obj$delta_beta

tissue_M <- obj$tissue_M

plasma_M <- obj$plasma_M

delta_M <- obj$delta_M

valid_pair_mask <- obj$valid_pair_mask


region_annotation <- as.data.table(
  obj$region_annotation
)


patient_metadata <- as.data.table(
  obj$patient_metadata
)


matched_metadata <- as.data.table(
  obj$matched_metadata
)


############################################################
## 6. INPUT VALIDATION
############################################################

message_header(
  "VALIDATING INPUT MATRICES"
)


stopifnot(

  is.matrix(tissue_beta),

  is.matrix(plasma_beta),

  is.matrix(delta_beta),

  is.matrix(tissue_M),

  is.matrix(plasma_M),

  is.matrix(delta_M),

  is.matrix(valid_pair_mask),

  ncol(delta_M) ==
    expected_matched_patients,

  nrow(patient_metadata) ==
    expected_matched_patients,

  nrow(matched_metadata) ==
    expected_matched_patients,

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
    dim(tissue_beta),
    dim(valid_pair_mask)
  ),

  nrow(region_annotation) ==
    nrow(delta_M)
)


patient_ids <- colnames(
  delta_M
)


stopifnot(

  identical(
    patient_ids,
    colnames(tissue_beta)
  ),

  identical(
    patient_ids,
    colnames(plasma_beta)
  ),

  identical(
    patient_ids,
    colnames(delta_beta)
  ),

  identical(
    patient_ids,
    patient_metadata$patient_id
  )
)


############################################################
## Check region IDs
############################################################

if (is.null(rownames(delta_M))) {

  stop(
    "delta_M matrix has no row names."
  )
}


stopifnot(
  "region_id" %in%
    names(region_annotation)
)


stopifnot(
  identical(
    as.character(
      region_annotation$region_id
    ),
    rownames(
      delta_M
    )
  )
)


cat(
  "Regions loaded:",
  format(
    nrow(delta_M),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Matched patients:",
  ncol(delta_M),
  "\n"
)

cat(
  "Input validation: PASS\n"
)


############################################################
## 7. VERIFY DELTA DEFINITIONS
############################################################

message_header(
  "VERIFYING DELTA DEFINITIONS"
)


beta_error <- max(
  abs(
    delta_beta -
      (
        plasma_beta -
          tissue_beta
      )
  ),
  na.rm = TRUE
)


M_error <- max(
  abs(
    delta_M -
      (
        plasma_M -
          tissue_M
      )
  ),
  na.rm = TRUE
)


cat(
  "Maximum Delta Beta reconstruction error:",
  beta_error,
  "\n"
)

cat(
  "Maximum Delta M reconstruction error:",
  M_error,
  "\n"
)


stopifnot(
  beta_error < 1e-12,
  M_error < 1e-12
)


cat(
  "Delta definitions: PASS\n"
)


############################################################
## 8. APPLY FROZEN VALID-PAIR MASK
############################################################

message_header(
  "APPLYING FROZEN VALID-PAIR MASK"
)


delta_M_masked <- delta_M

delta_beta_masked <- delta_beta

tissue_beta_masked <- tissue_beta

plasma_beta_masked <- plasma_beta


delta_M_masked[
  !valid_pair_mask
] <- NA_real_


delta_beta_masked[
  !valid_pair_mask
] <- NA_real_


tissue_beta_masked[
  !valid_pair_mask
] <- NA_real_


plasma_beta_masked[
  !valid_pair_mask
] <- NA_real_


############################################################
## Count valid pairs per region
############################################################

n_valid_pairs <- rowSums(

  valid_pair_mask &

    is.finite(
      delta_M_masked
    ) &

    is.finite(
      delta_beta_masked
    )
)


cat(
  "\nValid-pair count summary:\n"
)

print(
  summary(
    n_valid_pairs
  )
)


############################################################
## Retain analysis-eligible regions
############################################################

keep_analysis <-
  n_valid_pairs >=
  minimum_valid_pairs


cat(
  "\nRegions with >=",
  minimum_valid_pairs,
  " valid matched pairs:",
  format(
    sum(keep_analysis),
    big.mark = ","
  ),
  "of",
  format(
    length(keep_analysis),
    big.mark = ","
  ),
  "\n"
)


if (!any(keep_analysis)) {

  stop(
    "No regions have enough valid matched pairs."
  )
}


############################################################
## 9. CREATE FINAL ANALYSIS MATRICES
############################################################

delta_M_analysis <- delta_M_masked[
  keep_analysis,
  ,
  drop = FALSE
]


delta_beta_analysis <- delta_beta_masked[
  keep_analysis,
  ,
  drop = FALSE
]


tissue_beta_analysis <- tissue_beta_masked[
  keep_analysis,
  ,
  drop = FALSE
]


plasma_beta_analysis <- plasma_beta_masked[
  keep_analysis,
  ,
  drop = FALSE
]


analysis_annotation <- copy(
  region_annotation[
    keep_analysis
  ]
)


analysis_n_valid_pairs <-
  n_valid_pairs[
    keep_analysis
  ]


stopifnot(

  nrow(
    delta_M_analysis
  ) > 0L,

  ncol(
    delta_M_analysis
  ) ==
    expected_matched_patients,

  nrow(
    analysis_annotation
  ) ==
    nrow(
      delta_M_analysis
    ),

  identical(
    as.character(
      analysis_annotation$region_id
    ),
    rownames(
      delta_M_analysis
    )
  )
)


############################################################
## 10. PAIRED LIMMA MODEL
############################################################

message_header(
  "RUNNING PAIRED LIMMA MODEL ON DELTA M"
)


############################################################
## Each column is one matched patient.
##
## Intercept tests:
##
##      H0: mean Delta M = 0
##
## Because Delta M is already:
##
##      Plasma M - Tumor M
##
## this is a paired within-patient analysis.
############################################################


design <- matrix(
  1,
  nrow = ncol(
    delta_M_analysis
  ),
  ncol = 1,
  dimnames = list(
    colnames(
      delta_M_analysis
    ),
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
    "mean_delta_M_limma",
    "average_delta_M",
    "moderated_t",
    "p_value",
    "FDR",
    "B_statistic"
  ),
  skip_absent = TRUE
)


############################################################
## Restore matrix row order
############################################################

limma_results <- limma_results[
  match(
    rownames(
      delta_M_analysis
    ),
    region_id
  )
]


stopifnot(
  identical(
    limma_results$region_id,
    rownames(
      delta_M_analysis
    )
  )
)


############################################################
## 11. BETA-SCALE EFFECT SUMMARIES
############################################################

message_header(
  "CALCULATING BETA-SCALE EFFECT SIZES"
)


effect_summary <- data.table(

  region_id =
    rownames(
      delta_beta_analysis
    ),

  n_valid_pairs =
    analysis_n_valid_pairs,

  mean_tumor_beta =
    safe_row_means(
      tissue_beta_analysis
    ),

  mean_plasma_beta =
    safe_row_means(
      plasma_beta_analysis
    ),

  median_tumor_beta =
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
  abs_mean_delta_beta :=
    abs(
      mean_delta_beta
    )
]


effect_summary[
  ,
  abs_median_delta_beta :=
    abs(
      median_delta_beta
    )
]


############################################################
## Direction
############################################################

effect_summary[
  ,
  direction :=
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


############################################################
## 12. REMOVE DUPLICATED EFFECT COLUMNS FROM ANNOTATION
############################################################

effect_columns <- c(
  "n_valid_pairs",
  "mean_tissue_beta",
  "mean_tumor_beta",
  "mean_plasma_beta",
  "median_tissue_beta",
  "median_tumor_beta",
  "median_plasma_beta",
  "mean_delta_beta",
  "median_delta_beta",
  "abs_mean_delta_beta",
  "abs_median_delta_beta",
  "direction"
)


columns_to_remove <- intersect(
  effect_columns,
  names(
    analysis_annotation
  )
)


if (length(columns_to_remove) > 0L) {

  analysis_annotation[
    ,
    (columns_to_remove) := NULL
  ]
}


############################################################
## 13. COMBINE ANNOTATION + EFFECT + STATISTICS
############################################################

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


############################################################
## Restore original region order
############################################################

dmr_results <- dmr_results[
  match(
    rownames(
      delta_M_analysis
    ),
    region_id
  )
]


stopifnot(
  identical(
    dmr_results$region_id,
    rownames(
      delta_M_analysis
    )
  )
)


############################################################
## 14. PRIMARY SIGNIFICANCE CLASS
##
## Primary candidate definition:
##
## FDR < 0.05
## AND
## |median Delta Beta| >= 0.05
############################################################

dmr_results[
  ,
  significance_class :=
    fifelse(
      FDR < fdr_threshold &
        median_delta_beta >=
        primary_delta_beta_threshold,

      "Plasma_higher_DMR",

      fifelse(
        FDR < fdr_threshold &
          median_delta_beta <=
          -primary_delta_beta_threshold,

        "Tumor_higher_DMR",

        "Not_significant"
      )
    )
]


############################################################
## 15. ADD EFFECT-SIZE FLAGS
############################################################

dmr_results[
  ,
  FDR_significant :=
    FDR <
    fdr_threshold
]


dmr_results[
  ,
  abs_delta_beta_ge_0.05 :=
    abs_median_delta_beta >=
    0.05
]


dmr_results[
  ,
  abs_delta_beta_ge_0.10 :=
    abs_median_delta_beta >=
    0.10
]


dmr_results[
  ,
  abs_delta_beta_ge_0.20 :=
    abs_median_delta_beta >=
    0.20
]


############################################################
## Volcano metric
############################################################

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


############################################################
## 16. CREATE PRIMARY CANDIDATE TABLES
############################################################

candidate_DMRs <- dmr_results[
  significance_class !=
    "Not_significant"
]


plasma_higher_DMRs <- candidate_DMRs[
  significance_class ==
    "Plasma_higher_DMR"
]


tumor_higher_DMRs <- candidate_DMRs[
  significance_class ==
    "Tumor_higher_DMR"
]


cat(
  "\nPrimary candidate DMRs:",
  format(
    nrow(candidate_DMRs),
    big.mark = ","
  ),
  "\n"
)


cat(
  "Plasma-higher DMRs:",
  format(
    nrow(plasma_higher_DMRs),
    big.mark = ","
  ),
  "\n"
)


cat(
  "Tumor-higher DMRs:",
  format(
    nrow(tumor_higher_DMRs),
    big.mark = ","
  ),
  "\n"
)


############################################################
## 17. EFFECT-SIZE SENSITIVITY SUMMARY
############################################################

message_header(
  "DMR EFFECT-SIZE SENSITIVITY"
)


sensitivity_summary <- rbindlist(
  lapply(
    delta_beta_sensitivity_thresholds,
    function(threshold) {

      significant_effect <- (
        dmr_results$FDR <
          fdr_threshold
      ) &
        (
          dmr_results$abs_median_delta_beta >=
            threshold
        )


      plasma_effect <- significant_effect &
        (
          dmr_results$median_delta_beta >
            0
        )


      tumor_effect <- significant_effect &
        (
          dmr_results$median_delta_beta <
            0
        )


      data.table(

        FDR_threshold =
          fdr_threshold,

        abs_median_delta_beta_threshold =
          threshold,

        total_DMRs =
          sum(
            significant_effect,
            na.rm = TRUE
          ),

        plasma_higher_DMRs =
          sum(
            plasma_effect,
            na.rm = TRUE
          ),

        tumor_higher_DMRs =
          sum(
            tumor_effect,
            na.rm = TRUE
          )
      )
    }
  )
)


print(
  sensitivity_summary
)


############################################################
## 18. FDR-ONLY SUMMARY
############################################################

n_fdr_significant <- sum(
  dmr_results$FDR <
    fdr_threshold,
  na.rm = TRUE
)


cat(
  "\nRegions with FDR <",
  fdr_threshold,
  ":",
  format(
    n_fdr_significant,
    big.mark = ","
  ),
  "\n"
)


############################################################
## 19. DIRECTION SUMMARY
############################################################

direction_summary <- data.table(

  direction = c(
    "Plasma_higher",
    "Tumor_higher"
  ),

  n_all_regions = c(

    sum(
      dmr_results$median_delta_beta > 0,
      na.rm = TRUE
    ),

    sum(
      dmr_results$median_delta_beta < 0,
      na.rm = TRUE
    )
  ),

  n_primary_DMRs = c(

    nrow(
      plasma_higher_DMRs
    ),

    nrow(
      tumor_higher_DMRs
    )
  )
)


############################################################
## 20. ANALYSIS SUMMARY
############################################################

analysis_summary <- data.table(

  item = c(
    "regions_input",
    "regions_tested",
    "matched_patients",
    "minimum_valid_pairs",
    "FDR_threshold",
    "primary_abs_median_delta_beta_threshold",
    "FDR_significant_regions",
    "primary_candidate_DMRs",
    "plasma_higher_DMRs",
    "tumor_higher_DMRs"
  ),

  value = c(
    nrow(
      region_annotation
    ),

    nrow(
      dmr_results
    ),

    ncol(
      delta_M_analysis
    ),

    minimum_valid_pairs,

    fdr_threshold,

    primary_delta_beta_threshold,

    n_fdr_significant,

    nrow(
      candidate_DMRs
    ),

    nrow(
      plasma_higher_DMRs
    ),

    nrow(
      tumor_higher_DMRs
    )
  )
)


cat(
  "\nAnalysis summary:\n"
)

print(
  analysis_summary
)


############################################################
## 21. ANALYSIS SETTINGS
############################################################

analysis_settings <- data.table(

  setting = c(
    "statistical_model",
    "response",
    "paired_effect_definition",
    "multiple_testing",
    "FDR_threshold",
    "primary_median_delta_beta_threshold",
    "minimum_valid_pairs",
    "robust_ebayes",
    "trend_ebayes",
    "genome_build"
  ),

  value = c(
    "intercept_only_limma_on_within_patient_delta_M",
    "delta_M",
    "plasma_minus_tumor",
    "Benjamini_Hochberg",
    fdr_threshold,
    primary_delta_beta_threshold,
    minimum_valid_pairs,
    use_robust_ebayes,
    use_trend_ebayes,
    obj$settings$genome_build
  )
)


############################################################
## 22. FINAL VALIDATION BEFORE SAVING
############################################################

message_header(
  "FINAL VALIDATION"
)


stopifnot(

  nrow(
    dmr_results
  ) > 0L,

  nrow(
    dmr_results
  ) ==
    nrow(
      delta_M_analysis
    ),

  ncol(
    delta_M_analysis
  ) ==
    expected_matched_patients,

  "region_id" %in%
    names(
      dmr_results
    ),

  "mean_delta_M_limma" %in%
    names(
      dmr_results
    ),

  "median_delta_beta" %in%
    names(
      dmr_results
    ),

  "p_value" %in%
    names(
      dmr_results
    ),

  "FDR" %in%
    names(
      dmr_results
    ),

  "significance_class" %in%
    names(
      dmr_results
    ),

  !anyNA(
    dmr_results$region_id
  ),

  !anyDuplicated(
    dmr_results$region_id
  ),

  all(
    dmr_results$FDR >= 0 &
      dmr_results$FDR <= 1,
    na.rm = TRUE
  )
)


cat(
  "Final validation: PASS\n"
)


############################################################
## 23. SORT RESULTS FOR OUTPUT
############################################################

setorder(
  dmr_results,
  FDR,
  -abs_median_delta_beta
)


setorder(
  candidate_DMRs,
  FDR,
  -abs_median_delta_beta
)


setorder(
  plasma_higher_DMRs,
  FDR,
  -median_delta_beta
)


setorder(
  tumor_higher_DMRs,
  FDR,
  median_delta_beta
)


############################################################
## 24. SAVE MAIN RESULTS
############################################################

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
  tumor_higher_DMRs,
  file.path(
    dmr_dir,
    "SOLID_tissue_vs_plasma_tumor_higher_DMRs.tsv"
  ),
  sep = "\t"
)


############################################################
## 25. SAVE SUMMARY TABLES
############################################################

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


fwrite(
  sensitivity_summary,
  file.path(
    dmr_dir,
    "SOLID_tissue_vs_plasma_DMR_effect_size_sensitivity.tsv"
  ),
  sep = "\t"
)


fwrite(
  direction_summary,
  file.path(
    dmr_dir,
    "SOLID_tissue_vs_plasma_DMR_direction_summary.tsv"
  ),
  sep = "\t"
)


############################################################
## 26. SAVE LIMMA FIT
############################################################

saveRDS(
  fit,
  file.path(
    dmr_dir,
    "SOLID_tissue_vs_plasma_limma_fit.rds"
  )
)


############################################################
## 27. CREATE REUSABLE DMR OBJECT
############################################################

dmr_object <- list(

  all_results =
    dmr_results,

  candidate_DMRs =
    candidate_DMRs,

  plasma_higher_DMRs =
    plasma_higher_DMRs,

  tumor_higher_DMRs =
    tumor_higher_DMRs,

  delta_M =
    delta_M_analysis,

  delta_beta =
    delta_beta_analysis,

  tissue_beta =
    tissue_beta_analysis,

  plasma_beta =
    plasma_beta_analysis,

  valid_pair_mask =
    valid_pair_mask[
      keep_analysis,
      ,
      drop = FALSE
    ],

  patient_metadata =
    patient_metadata,

  matched_metadata =
    matched_metadata,

  analysis_summary =
    analysis_summary,

  sensitivity_summary =
    sensitivity_summary,

  analysis_settings =
    analysis_settings,

  provenance = list(

    source_object =
      input_file,

    script =
      "07_SOLID_paired_DMR.r",

    delta_M_definition =
      "plasma_minus_tumor",

    delta_beta_definition =
      "plasma_minus_tumor",

    creation_date =
      as.character(
        Sys.Date()
      )
  )
)


dmr_object_file <- file.path(
  dmr_dir,
  "SOLID_tissue_vs_plasma_DMR_object.rds"
)


saveRDS(
  dmr_object,
  dmr_object_file
)


############################################################
## 28. RELOAD AND VALIDATE SAVED OBJECT
############################################################

validation_object <- readRDS(
  dmr_object_file
)


stopifnot(

  nrow(
    validation_object$all_results
  ) > 0L,

  nrow(
    validation_object$all_results
  ) ==
    nrow(
      delta_M_analysis
    ),

  ncol(
    validation_object$delta_M
  ) ==
    expected_matched_patients,

  nrow(
    validation_object$patient_metadata
  ) ==
    expected_matched_patients
)


rm(
  validation_object
)


cat(
  "\nSaved DMR object validation: PASS\n"
)


############################################################
## 29. FINAL CONSOLE SUMMARY
############################################################

cat(
  "\n============================================\n"
)

cat(
  "SOLID PAIRED DMR ANALYSIS COMPLETE\n"
)

cat(
  "============================================\n"
)


cat(
  "\nRegions tested:",
  format(
    nrow(
      dmr_results
    ),
    big.mark = ","
  ),
  "\n"
)


cat(
  "Matched patients:",
  ncol(
    delta_M_analysis
  ),
  "\n"
)


cat(
  "Minimum valid pairs:",
  minimum_valid_pairs,
  "\n"
)


cat(
  "\nFDR-significant regions:",
  format(
    n_fdr_significant,
    big.mark = ","
  ),
  "\n"
)


cat(
  "Primary candidate DMRs:",
  format(
    nrow(
      candidate_DMRs
    ),
    big.mark = ","
  ),
  "\n"
)


cat(
  "Plasma-higher DMRs:",
  format(
    nrow(
      plasma_higher_DMRs
    ),
    big.mark = ","
  ),
  "\n"
)


cat(
  "Tumor-higher DMRs:",
  format(
    nrow(
      tumor_higher_DMRs
    ),
    big.mark = ","
  ),
  "\n"
)


cat(
  "\nEffect-size sensitivity:\n"
)

print(
  sensitivity_summary
)


cat(
  "\nMedian Delta Beta among primary DMRs:\n"
)

if (nrow(candidate_DMRs) > 0L) {

  print(
    summary(
      candidate_DMRs$median_delta_beta
    )
  )

} else {

  cat(
    "No primary DMRs identified.\n"
  )
}


cat(
  "\nOutputs saved to:\n",
  dmr_dir,
  "\n",
  sep = ""
)


############################################################
## 30. SESSION INFORMATION
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    dmr_dir,
    "sessionInfo_07_SOLID_paired_DMR.txt"
  )
)