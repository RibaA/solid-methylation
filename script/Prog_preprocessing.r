##------------------------------------------------------------
# DNA METHYLATION PREPROCESSING
# Brain IDH-mutant FFPE tissue
#
# Platforms:
#   - Illumina MethylationEPIC v1: 3 samples
#   - Illumina MethylationEPIC v2: 26 samples
#
# Current objective:
#   1. Load raw RGChannelSet objects
#   2. Verify and repair sample metadata
#   3. Calculate detection P-values
#   4. Perform sample-level QC
#   5. Normalize EPIC v1 and EPIC v2 separately using Noob
#   6. Filter unreliable probes
#   7. Generate Beta and M-value matrices
#   8. Save processed objects
#
# Important:
#   EPIC v1 and EPIC v2 are NOT combined in this script.
##------------------------------------------------------------
##################################################
## library
#################################################
library(minfi)
library(SummarizedExperiment)
library(S4Vectors)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

#################################################
## functions
#################################################
extract_expected_metadata <- function(rgSet) {

  full_name <- sampleNames(rgSet)

  data.frame(
    full_sample_name = full_name,

    expected_Sample_Name = sub(
      "_.*$",
      "",
      full_name
    ),

    expected_Sentrix_ID = sub(
      "^[^_]+_([^_]+)_.*$",
      "\\1",
      full_name
    ),

    expected_Sentrix_Position = sub(
      "^.*_([^_]+)$",
      "\\1",
      full_name
    ),

    stringsAsFactors = FALSE
  )
}

compare_metadata <- function(rgSet) {

  existing <- as.data.frame(
    colData(rgSet)
  )

  expected <- extract_expected_metadata(
    rgSet
  )

  out <- data.frame(
    full_sample_name =
      sampleNames(rgSet),

    existing_Sample_Name =
      as.character(existing$Sample_Name),

    expected_Sample_Name =
      expected$expected_Sample_Name,

    existing_Sentrix_ID =
      as.character(existing$Sentrix_ID),

    expected_Sentrix_ID =
      expected$expected_Sentrix_ID,

    existing_Sentrix_Position =
      as.character(existing$Sentrix_Position),

    expected_Sentrix_Position =
      expected$expected_Sentrix_Position,

    stringsAsFactors = FALSE
  )

  out <- out %>%
    dplyr::mutate(
      Sample_Name_match =
        existing_Sample_Name ==
        expected_Sample_Name,

      Sentrix_ID_match =
        existing_Sentrix_ID ==
        expected_Sentrix_ID,

      Sentrix_Position_match =
        existing_Sentrix_Position ==
        expected_Sentrix_Position,

      any_mismatch =
        !Sample_Name_match |
        !Sentrix_ID_match |
        !Sentrix_Position_match
    )

  out
}

repair_array_metadata <- function(rgSet) {

  full_names <- sampleNames(rgSet)

  pd <- as.data.frame(
    colData(rgSet)
  )

  ## Preserve original metadata
  pd$Sample_Name_original <-
    as.character(pd$Sample_Name)

  pd$Sentrix_ID_original <-
    as.character(pd$Sentrix_ID)

  pd$Sentrix_Position_original <-
    as.character(pd$Sentrix_Position)

  ## Recover correct values from RGChannelSet column names
  pd$Sample_Name <-
    sub(
      "_.*$",
      "",
      full_names
    )

  pd$Sentrix_ID <-
    sub(
      "^[^_]+_([^_]+)_.*$",
      "\\1",
      full_names
    )

  pd$Sentrix_Position <-
    sub(
      "^.*_([^_]+)$",
      "\\1",
      full_names
    )

  rownames(pd) <- full_names

  colData(rgSet) <- S4Vectors::DataFrame(
    pd,
    row.names = full_names
  )

  rgSet
}

#################################################
## setup directory
#################################################
dir_input  <- "data/proc"
dir_output <- "result/data"
dir_figures <- "result/figures"

#################################################
## load data
#################################################
load(file.path(dir_input, 'RGSets_targets.RData'))
ls()

class(rgSet_EPICv1)
class(rgSet_EPICv2)

dim(rgSet_EPICv1)
dim(rgSet_EPICv2)

annotation(rgSet_EPICv1)
annotation(rgSet_EPICv2)

stopifnot(inherits(rgSet_EPICv1, "RGChannelSet"))
stopifnot(inherits(rgSet_EPICv2, "RGChannelSet"))
stopifnot(ncol(rgSet_EPICv1) == 3)
stopifnot(ncol(rgSet_EPICv2) == 26)

#################################################
## metadata issue
#################################################
expected_v1 <- extract_expected_metadata(
  rgSet_EPICv1
)

expected_v2 <- extract_expected_metadata(
  rgSet_EPICv2
)

expected_v1
head(expected_v2)

metadata_issue_v1 <- compare_metadata(
  rgSet_EPICv1
)

metadata_issue_v2 <- compare_metadata(
  rgSet_EPICv2
)

metadata_summary <- data.frame(
  Array = c(
    "EPIC",
    "EPICv2"
  ),

  total_samples = c(
    nrow(metadata_issue_v1),
    nrow(metadata_issue_v2)
  ),

  sample_name_mismatches = c(
    sum(!metadata_issue_v1$Sample_Name_match),
    sum(!metadata_issue_v2$Sample_Name_match)
  ),

  sentrix_id_mismatches = c(
    sum(!metadata_issue_v1$Sentrix_ID_match),
    sum(!metadata_issue_v2$Sentrix_ID_match)
  ),

  sentrix_position_mismatches = c(
    sum(!metadata_issue_v1$Sentrix_Position_match),
    sum(!metadata_issue_v2$Sentrix_Position_match)
  )
)

metadata_issue_all <- dplyr::bind_rows(
  metadata_issue_v1 %>%
    dplyr::mutate(Array_version = "EPIC"),

  metadata_issue_v2 %>%
    dplyr::mutate(Array_version = "EPICv2")
)

write.csv(metadata_issue_all, file = file.path(dir_output,"01_metadata_issue_before_repair.csv"),row.names = FALSE)

#################################################
## recover metadata
#################################################
rgSet_EPICv1_fixed <- repair_array_metadata(
  rgSet_EPICv1
)

rgSet_EPICv2_fixed <- repair_array_metadata(
  rgSet_EPICv2
)

metadata_check_v1_fixed <- compare_metadata(
  rgSet_EPICv1_fixed
)

metadata_check_v2_fixed <- compare_metadata(
  rgSet_EPICv2_fixed
)

stopifnot(
  all(metadata_check_v1_fixed$Sample_Name_match),
  all(metadata_check_v1_fixed$Sentrix_ID_match),
  all(metadata_check_v1_fixed$Sentrix_Position_match),

  all(metadata_check_v2_fixed$Sample_Name_match),
  all(metadata_check_v2_fixed$Sentrix_ID_match),
  all(metadata_check_v2_fixed$Sentrix_Position_match)
)

saveRDS(rgSet_EPICv1_fixed, file = file.path(dir_output, "rgSet_EPICv1_fixed.rds"))
saveRDS(rgSet_EPICv2_fixed, file = file.path(dir_output, "rgSet_EPICv2_fixed.rds"))

metadata_fixed <- rbind(
  cbind(
    Array_version = "EPIC",
    as.data.frame(
      colData(rgSet_EPICv1_fixed)
    )
  ),
  cbind(
    Array_version = "EPICv2",
    as.data.frame(
      colData(rgSet_EPICv2_fixed)
    )
  )
)

write.csv(metadata_fixed, file = file.path(dir_output,"metadata_after_repair.csv"), row.names = FALSE)

#################################################
## QC: calculate detection p-value
#################################################
detP_EPICv1 <- minfi::detectionP(rgSet_EPICv1_fixed)
detP_EPICv2 <- minfi::detectionP(rgSet_EPICv2_fixed)