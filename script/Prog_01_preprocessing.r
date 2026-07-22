##------------------------------------------------------------
# DNA METHYLATION PREPROCESSING
# Brain IDH-mutant FFPE tissue
#
# Platforms:
# - Illumina MethylationEPIC v1: 3 samples
# - Illumina MethylationEPIC v2: 26 samples
#
# Input data:
# - Raw paired red and green IDAT files
# - IDAT-to-subject mapping file
# - Clinical and sample-level metadata
#
# Current objective:
# 1. Identify raw IDAT files and construct the targets table
# 2. Annotate array platform and annotation genome build
# 3. Merge IDAT identifiers with sample metadata
# 4. Load raw EPIC v1 and EPIC v2 RGChannelSet objects
# 5. Calculate detection P-values
# 6. Perform sample-level quality control
# 7. Normalize EPIC v1 and EPIC v2 separately using Noob
# 8. Filter unreliable probes
# 9. Generate Beta-value and M-value matrices
# 10. Save raw and processed objects
#
# Important:
# - EPIC v1 and EPIC v2 are processed separately because they
# use different array manifests and annotation resources.
# - The annotation resources currently use hg19 for EPIC v1
# and hg38 for EPIC v2.
# - EPIC v1 and EPIC v2 are NOT combined in this script.
##------------------------------------------------------------
##################################################
## Set up library and R version
#################################################
user_lib <- "C:/Users/farno/R/win-library/4.5"

.libPaths(
  c(
    user_lib,
    .libPaths()
  )
)

.libPaths()

Sys.setenv(
  COMSPEC = "C:\\Windows\\System32\\cmd.exe"
)

Sys.getenv("COMSPEC")

##################################################
## library
#################################################
library(minfi)
library(SummarizedExperiment)
library(S4Vectors)
library(GenomicRanges)
library(dplyr)
library(stringr)
library(tibble)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

#################################################
## setup directory
#################################################
dir_input  <- "data/raw"
dir_output <- "data/proc"

############################################################
## Identify all IDAT files
############################################################
idat_files <- list.files(
  dir_input,
  pattern = "_(Grn|Red)\\.idat$",
  full.names = TRUE
)

length(idat_files)
head(basename(idat_files))

############################################################
## Make targets file
############################################################
targets <- tibble(
  filename = basename(idat_files),
  full_path = idat_files
) %>%
  mutate(
    # Remove _Grn.idat or _Red.idat
    file_stem = str_remove(
      filename,
      "_(Grn|Red)\\.idat$"
    ),

    # GEO accession
    Sample_Name = str_extract(
      file_stem,
      "^GSM[0-9]+"
    ),

    # Keep the complete physical filename stem
    Basename = file.path(
      dir_input,
      file_stem
    ),

    # Extract Illumina array identifiers
    Sentrix_ID = str_extract(
      file_stem,
      "(?<=_)[0-9]+(?=_R[0-9]{2}C[0-9]{2}$)"
    ),

    Sentrix_Position = str_extract(
      file_stem,
      "R[0-9]{2}C[0-9]{2}$"
    )
  ) %>%
  distinct(
    Sample_Name,
    Basename,
    Sentrix_ID,
    Sentrix_Position
  ) %>%
  arrange(Sample_Name)

dim(targets)
head(targets)

###########################################################
## Annotate array type
###########################################################
targets <- targets %>%
  mutate(
    Array = if_else(
      Sample_Name %in% c(
        "GSM9325970",
        "GSM9325971",
        "GSM9325972"
      ),
      "EPIC",
      "EPICv2"
    ),
    Genome_build = if_else(
      Array == "EPIC",
      "hg19",
      "hg38"
    )
  )

targets %>%
  count(Array, Genome_build)

stopifnot(
  all(file.exists(paste0(targets$Basename, "_Grn.idat")))
)

stopifnot(
  all(file.exists(paste0(targets$Basename, "_Red.idat")))
)

stopifnot(nrow(targets) == 29)
stopifnot(!anyDuplicated(targets$Sample_Name))
stopifnot(!anyDuplicated(targets$Basename))

#################################################
## load metadata and mapping
#################################################
meta_data <- read.csv('data/meta_data_research.csv')
idat_meta <- read.table('data/idat_meta.tsv', header = TRUE)

idat_meta <- idat_meta %>%
  mutate(
    SampleID = trimws(as.character(SampleID)),
    SolidID  = trimws(as.character(SolidID))
  )

meta_data <- meta_data %>%
  mutate(
    Subject = trimws(as.character(Subject))
  )

idat_meta <- idat_meta %>%
  mutate(
    SolidID = if_else(
      SolidID == "SOLID-006",
      "SOLID-009",
      SolidID
    )
  )

sample_mapping <- idat_meta %>%
  transmute(
    Sample_Name = SampleID,
    Subject = SolidID
  )

all(meta_data$Subject == meta_data$X, na.rm = TRUE)

###########################################################
## Add phenotype data
###########################################################
targets <- targets %>%
  left_join(
    sample_mapping,
    by = "Sample_Name"
  ) %>%
  left_join(
    meta_data,
    by = "Subject"
  ) %>%
  relocate(
    Sample_Name,
    Subject,
    Array,
    Genome_build,
    Sentrix_ID,
    Sentrix_Position,
    Basename
  )

targets %>%
  filter(is.na(Subject)) %>%
  select(
    Sample_Name,
    Sentrix_ID,
    Sentrix_Position
  )

###########################################################
## Get RGSets
###########################################################
targets_EPICv1 <- targets %>%
  filter(Array == "EPIC")

targets_EPICv2 <- targets %>%
  filter(Array == "EPICv2")

stopifnot(nrow(targets_EPICv1) == 3)
stopifnot(nrow(targets_EPICv2) == 26)

rgSet_EPICv1 <- read.metharray.exp(
  targets = targets_EPICv1,
  extended = TRUE,
  verbose = TRUE
)

rgSet_EPICv2 <- read.metharray.exp(
  targets = targets_EPICv2,
  extended = TRUE,
  verbose = TRUE
)

sampleNames(rgSet_EPICv1) <- targets_EPICv1$Sample_Name
sampleNames(rgSet_EPICv2) <- targets_EPICv2$Sample_Name

sampleNames(rgSet_EPICv1)
sampleNames(rgSet_EPICv2)

annotation(rgSet_EPICv1)
annotation(rgSet_EPICv2)

dim(rgSet_EPICv1)
dim(rgSet_EPICv2)

###########################################################
## Save checkpoint
###########################################################
save(
  targets,
  targets_EPICv1,
  targets_EPICv2,
  rgSet_EPICv1,
  rgSet_EPICv2,
  file = file.path(
    dir_output,
    "RGSets_targets.RData"
  )
)
