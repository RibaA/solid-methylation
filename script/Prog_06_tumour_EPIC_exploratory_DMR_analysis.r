##------------------------------------------------------------
## TUMOUR EPIC METHYLATION
## Exploratory visualization, differential methylation,
## regional DMR analysis, and matched SOLID tumour evaluation
##
## Main objectives:
##   1. Load and validate harmonized EPIC v1/v2 tumour data
##   2. Identify tumour samples matched to SOLID cfDNA cases
##   3. Characterize methylation distributions and sample structure
##   4. Perform probe-level differential methylation analysis
##   5. Identify regional differentially methylated regions
##   6. Evaluate selected DMRs in matched tumour samples
##   7. Export genomic regions for projection into SOLID cfDNA
##------------------------------------------------------------
##############################################################
## load library
##############################################################
options(
  stringsAsFactors = FALSE,
  warn = 1
)

suppressPackageStartupMessages({
  library(minfi)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(matrixStats)
  library(limma)
  library(ComplexHeatmap)
  library(circlize)
  library(ggrepel)
  library(patchwork)
})

##############################################################
# set up directory
##############################################################
dir_input <- 'data/proc'
dir_output <- 'result/tumor-epic'

dir_figures <- file.path(
  dir_output,
  "figures"
)

dir_tables <- file.path(
  dir_output,
  "tables"
)

dir_objects <- file.path(
  dir_output,
  "objects"
)

dir.create(
  dir_figures,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dir_tables,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dir_objects,
  recursive = TRUE,
  showWarnings = FALSE
)

##------------------------------------------------------------
## 2. Inspect processed tumour EPIC objects
##------------------------------------------------------------
file_tumour_objects <- file.path(
  dir_input,
  "RGSets_targets.RData"
)

stopifnot(file.exists(file_tumour_objects))

objects_before_load <- ls()

load(file_tumour_objects)

objects_after_load <- setdiff(
  ls(),
  objects_before_load
)

cat("Objects loaded from RGSets_targets.RData:\n")
print(objects_after_load)