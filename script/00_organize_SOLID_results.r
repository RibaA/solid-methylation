############################################################
## 00_organize_SOLID_results.R
##
## Purpose:
##   Reorganize the current SOLID result folder into a
##   cleaner downstream-analysis structure.
##
## IMPORTANT:
##   - COPY files only
##   - DO NOT delete or overwrite existing results
##   - Keep original folders intact for verification
############################################################


############################################################
## 0. CLEAN ENVIRONMENT
############################################################

rm(list = ls())

options(
  stringsAsFactors = FALSE,
  width = 180,
  warn = 1
)


############################################################
## 1. PROJECT DIRECTORY
############################################################

project_dir <- "C:/solid-methylation"

result_dir <- file.path(
  project_dir,
  "result"
)

stopifnot(
  dir.exists(result_dir)
)


############################################################
## 2. NEW RESULT STRUCTURE
############################################################

new_dirs <- c(

  file.path(
    result_dir,
    "01_tissue_EPIC",
    "data"
  ),

  file.path(
    result_dir,
    "01_tissue_EPIC",
    "qc"
  ),

  file.path(
    result_dir,
    "02_plasma_5base",
    "figures"
  ),

  file.path(
    result_dir,
    "02_plasma_5base",
    "tables"
  ),

  file.path(
    result_dir,
    "03_matched_tissue_plasma",
    "preparation"
  ),

  file.path(
    result_dir,
    "03_matched_tissue_plasma",
    "paired_filtered"
  ),

  file.path(
    result_dir,
    "03_matched_tissue_plasma",
    "concordance"
  ),

  file.path(
    result_dir,
    "03_matched_tissue_plasma",
    "PCA_clustering"
  ),

  file.path(
    result_dir,
    "03_matched_tissue_plasma",
    "DMR"
  ),

  file.path(
    result_dir,
    "03_matched_tissue_plasma",
    "figures"
  ),

  file.path(
    result_dir,
    "03_matched_tissue_plasma",
    "tables"
  ),

  file.path(
    result_dir,
    "04_clinical",
    "metadata"
  ),

  file.path(
    result_dir,
    "04_clinical",
    "tables"
  ),

  file.path(
    result_dir,
    "05_integration",
    "OCTANE_comparison"
  ),

  file.path(
    result_dir,
    "05_integration",
    "candidate_regions"
  ),

  file.path(
    result_dir,
    "05_integration",
    "annotation"
  ),

  file.path(
    result_dir,
    "05_integration",
    "figures"
  ),

  file.path(
    result_dir,
    "Archive"
  )
)


for (d in new_dirs) {

  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


############################################################
## 3. HELPER: SAFE COPY
############################################################

safe_copy <- function(
    from,
    to_dir
) {

  if (!file.exists(from)) {

    warning(
      "Source does not exist: ",
      from
    )

    return(
      FALSE
    )
  }

  dir.create(
    to_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  destination <- file.path(
    to_dir,
    basename(from)
  )

  if (file.exists(destination)) {

    message(
      "SKIP existing: ",
      destination
    )

    return(
      TRUE
    )
  }

  ok <- file.copy(
    from = from,
    to = destination,
    overwrite = FALSE,
    copy.mode = TRUE,
    copy.date = TRUE
  )

  if (ok) {

    message(
      "COPIED: ",
      from,
      " -> ",
      destination
    )

  } else {

    warning(
      "COPY FAILED: ",
      from
    )
  }

  ok
}


############################################################
## 4. HELPER: COPY CONTENTS OF DIRECTORY
############################################################

copy_directory_contents <- function(
    source_dir,
    destination_dir
) {

  if (!dir.exists(source_dir)) {

    warning(
      "Directory does not exist: ",
      source_dir
    )

    return(
      NULL
    )
  }

  files <- list.files(
    source_dir,
    full.names = TRUE,
    recursive = TRUE,
    include.dirs = FALSE
  )

  if (length(files) == 0L) {

    message(
      "No files found in: ",
      source_dir
    )

    return(
      NULL
    )
  }

  for (f in files) {

    relative_path <- substring(
      f,
      nchar(source_dir) + 2L
    )

    destination_file <- file.path(
      destination_dir,
      relative_path
    )

    destination_parent <- dirname(
      destination_file
    )

    dir.create(
      destination_parent,
      recursive = TRUE,
      showWarnings = FALSE
    )

    if (file.exists(destination_file)) {

      message(
        "SKIP existing: ",
        destination_file
      )

      next
    }

    ok <- file.copy(
      from = f,
      to = destination_file,
      overwrite = FALSE,
      copy.mode = TRUE,
      copy.date = TRUE
    )

    if (ok) {

      message(
        "COPIED: ",
        f
      )

    } else {

      warning(
        "COPY FAILED: ",
        f
      )
    }
  }
}


############################################################
## 5. COPY CURRENT EPIC RESULTS
############################################################

copy_directory_contents(
  source_dir = file.path(
    result_dir,
    "data"
  ),
  destination_dir = file.path(
    result_dir,
    "01_tissue_EPIC",
    "data"
  )
)

copy_directory_contents(
  source_dir = file.path(
    result_dir,
    "qc"
  ),
  destination_dir = file.path(
    result_dir,
    "01_tissue_EPIC",
    "qc"
  )
)


############################################################
## 6. COPY CURRENT PLASMA RESULTS
############################################################

copy_directory_contents(
  source_dir = file.path(
    result_dir,
    "solid-5base",
    "figures"
  ),
  destination_dir = file.path(
    result_dir,
    "02_plasma_5base",
    "figures"
  )
)

copy_directory_contents(
  source_dir = file.path(
    result_dir,
    "solid-5base",
    "tables"
  ),
  destination_dir = file.path(
    result_dir,
    "02_plasma_5base",
    "tables"
  )
)


############################################################
## 7. COPY MATCHED PREPARATION RESULTS
############################################################

copy_directory_contents(
  source_dir = file.path(
    result_dir,
    "tissue-plasma",
    "matched_preparation"
  ),
  destination_dir = file.path(
    result_dir,
    "03_matched_tissue_plasma",
    "preparation"
  )
)


############################################################
## 8. COPY PAIRED FILTERED RESULTS
############################################################

copy_directory_contents(
  source_dir = file.path(
    result_dir,
    "tissue-plasma",
    "paired_filtered"
  ),
  destination_dir = file.path(
    result_dir,
    "03_matched_tissue_plasma",
    "paired_filtered"
  )
)


############################################################
## 9. COPY PAIRED DMR RESULTS
############################################################

copy_directory_contents(
  source_dir = file.path(
    result_dir,
    "tissue-plasma",
    "paired_DMR"
  ),
  destination_dir = file.path(
    result_dir,
    "03_matched_tissue_plasma",
    "DMR"
  )
)


############################################################
## 10. COPY DMR FIGURES
############################################################

copy_directory_contents(
  source_dir = file.path(
    result_dir,
    "tissue-plasma",
    "DMR_figures"
  ),
  destination_dir = file.path(
    result_dir,
    "03_matched_tissue_plasma",
    "figures"
  )
)


############################################################
## 11. COPY DMR TABLES / ANNOTATION
############################################################

copy_directory_contents(
  source_dir = file.path(
    result_dir,
    "tissue-plasma",
    "DMR_annotation"
  ),
  destination_dir = file.path(
    result_dir,
    "03_matched_tissue_plasma",
    "tables"
  )
)


############################################################
## 12. COPY CLINICAL METADATA
############################################################

clinical_source <- file.path(
  project_dir,
  "data",
  "meta_data_research.csv"
)

safe_copy(
  from = clinical_source,
  to_dir = file.path(
    result_dir,
    "04_clinical",
    "metadata"
  )
)


############################################################
## 13. COPY OLD tumor-epic OUTPUTS TO ARCHIVE
############################################################

copy_directory_contents(
  source_dir = file.path(
    result_dir,
    "tumor-epic"
  ),
  destination_dir = file.path(
    result_dir,
    "Archive",
    "tumor-epic_old"
  )
)


############################################################
## 14. WRITE ORGANIZATION LOG
############################################################

organization_log <- data.frame(

  item = c(
    "01_tissue_EPIC",
    "02_plasma_5base",
    "03_matched_tissue_plasma",
    "04_clinical",
    "05_integration",
    "Archive"
  ),

  purpose = c(
    "EPIC preprocessing, harmonization and QC",
    "SOLID plasma 5-base QC and exploratory results",
    "Matched tumor-plasma preparation and downstream analysis",
    "Clinical metadata and derived clinical tables",
    "Future SOLID-OCTANE integrated analyses",
    "Legacy/older exploratory outputs"
  ),

  stringsAsFactors = FALSE
)

write.table(
  organization_log,
  file.path(
    result_dir,
    "RESULT_FOLDER_STRUCTURE.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
## 15. FINAL SUMMARY
############################################################

cat(
  "\n============================================\n"
)

cat(
  "SOLID RESULT ORGANIZATION COMPLETE\n"
)

cat(
  "============================================\n"
)

cat(
  "\nIMPORTANT:\n"
)

cat(
  "- Files were COPIED, not moved.\n"
)

cat(
  "- Original result folders are unchanged.\n"
)

cat(
  "- No existing destination files were overwritten.\n"
)

cat(
  "\nNew top-level result folders:\n"
)

print(
  list.files(
    result_dir
  )
)

cat(
  "\nNew structure written to:\n"
)

cat(
  file.path(
    result_dir,
    "RESULT_FOLDER_STRUCTURE.tsv"
  ),
  "\n"
)
