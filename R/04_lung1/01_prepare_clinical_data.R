# 01_prepare_clinical_data.R
#
# Prepare and audit Lung1 clinical and survival data.
#
# Run from an analysis workspace configured through environment variables.

options(stringsAsFactors = FALSE)
options(timeout = 1800)
options(download.file.method = "libcurl")


# Analysis paths are supplied through environment variables.
required_env_path <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    stop("Environment variable ", name, " is not set.")
  }
  normalizePath(value, winslash = "/", mustWork = TRUE)
}

project_dir <- required_env_path("RADIOGENOMICS_PROJECT_DIR")
setwd(project_dir)

lung1_root <- required_env_path("LUNG1_DATA_DIR")

############################################################
# 1. Packages
############################################################

required_packages <- c("openxlsx")
missing_packages <- required_packages[!vapply(
  required_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)]
if (length(missing_packages) > 0L) {
  stop(
    "Missing R package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install dependencies with environment/install_r_dependencies.R."
  )
}
invisible(lapply(required_packages, library, character.only = TRUE))

library(openxlsx)

############################################################
# 2. Paths
############################################################

# The analysis workspace contains scripts, metadata, and model outputs.

# Large Lung1 imaging files are stored in the external data directory.

clinical_dir <- file.path(lung1_root, "clinical")
manifest_dir <- file.path(lung1_root, "TCIA_manifest")
dicom_dir <- file.path(lung1_root, "TCIA_DICOM")
nrrd_dir <- file.path(lung1_root, "NRRD_full")
metadata_dir <- file.path(lung1_root, "metadata")
results_dir <- file.path(lung1_root, "results")

dir.create(lung1_root, recursive = TRUE, showWarnings = FALSE)
dir.create(clinical_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dicom_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(nrrd_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n===== Lung1 analysis folders created =====\n")
cat("Lung1 root:\n")
print(lung1_root)

cat("\nFolder list:\n")
print(list.dirs(lung1_root, recursive = FALSE))

############################################################
# 3. Clinical file paths
############################################################

clinical_csv <- file.path(clinical_dir, "Lung1_clinical.csv")

# Copy an existing clinical file into the configured Lung1 data directory when available.
old_c_clinical_csv <- file.path(
  project_dir,
  "01_raw_data",
  "NSCLC_Radiomics_Lung1",
  "clinical",
  "Lung1_clinical.csv"
)

if (!file.exists(clinical_csv) && file.exists(old_c_clinical_csv)) {
  old_size <- file.info(old_c_clinical_csv)$size
  
  if (is.finite(old_size) && old_size > 1000) {
    cat("\nDetected an existing clinical CSV; copying it to the configured Lung1 workspace...\n")
    file.copy(old_c_clinical_csv, clinical_csv, overwrite = TRUE)
  }
}

############################################################
# 4. Download clinical CSV
############################################################

clinical_url <- "https://wiki.cancerimagingarchive.net/download/attachments/16056854/NSCLC%20Radiomics%20Lung1.clinical-version3-Oct%202019.csv?api=v2&modificationDate=1572013183040&version=1"

tcia_page_url <- "https://wiki.cancerimagingarchive.net/display/Public/NSCLC-Radiomics"

download_ok <- FALSE

cat("\n===== Download Lung1 clinical CSV =====\n")
cat("Output path:\n")
print(clinical_csv)

if (file.exists(clinical_csv)) {
  fsize <- file.info(clinical_csv)$size
  
  if (is.finite(fsize) && fsize > 1000) {
    cat("Clinical CSV already exists and looks valid. Skip download.\n")
    download_ok <- TRUE
  } else {
    cat("Existing clinical CSV is too small. Removing incomplete file.\n")
    unlink(clinical_csv, force = TRUE)
  }
}

if (!download_ok) {
  
  cat("\nTrying download method 1: utils::download.file with libcurl...\n")
  
  try({
    utils::download.file(
      url = clinical_url,
      destfile = clinical_csv,
      mode = "wb",
      quiet = FALSE,
      method = "libcurl"
    )
  }, silent = TRUE)
  
  if (file.exists(clinical_csv)) {
    fsize <- file.info(clinical_csv)$size
    if (is.finite(fsize) && fsize > 1000) {
      download_ok <- TRUE
    } else {
      unlink(clinical_csv, force = TRUE)
    }
  }
}

if (!download_ok) {
  
  cat("\nTrying download method 2: curl::curl_download...\n")
  
  if (!requireNamespace("curl", quietly = TRUE)) {
  stop("Missing R package 'curl'. Install dependencies with environment/install_r_dependencies.R.")
}
  
  if (requireNamespace("curl", quietly = TRUE)) {
    try({
      curl::curl_download(
        url = clinical_url,
        destfile = clinical_csv,
        quiet = FALSE
      )
    }, silent = TRUE)
  }
  
  if (file.exists(clinical_csv)) {
    fsize <- file.info(clinical_csv)$size
    if (is.finite(fsize) && fsize > 1000) {
      download_ok <- TRUE
    } else {
      unlink(clinical_csv, force = TRUE)
    }
  }
}

############################################################
# 5. If automatic download fails, open browser for manual download
############################################################

if (!download_ok) {
  
  cat("\n############################################################\n")
  cat("Automatic download failed, likely due to SSL connection error in R.\n")
  cat("This is not a data problem. Please manually download the clinical CSV.\n\n")
  
  cat("Step 1: The script will open the TCIA NSCLC-Radiomics page.\n")
  cat("Step 2: On the page, find Data Access -> Lung1 clinical (CSV) -> Download.\n")
  cat("Step 3: Save the downloaded file as exactly:\n")
  cat(clinical_csv, "\n")
  cat("Step 4: Re-run this script after saving the file.\n")
  cat("############################################################\n\n")
  
  try({
    browseURL(tcia_page_url)
  }, silent = TRUE)
  
  stop(
    "Clinical CSV was not downloaded automatically. ",
    "Please manually save it to: ",
    clinical_csv,
    " then re-run this script."
  )
}

if (!file.exists(clinical_csv)) {
  stop("Clinical CSV still not found: ", clinical_csv)
}

cat("\nClinical CSV ready:\n")
print(clinical_csv)

cat("\nClinical CSV file size:\n")
print(file.info(clinical_csv)$size)

############################################################
# 6. Read and audit clinical table
############################################################

clinical_df <- read.csv(
  clinical_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\n===== Lung1 clinical table =====\n")
cat("Clinical table dim:\n")
print(dim(clinical_df))

cat("\nClinical columns:\n")
print(colnames(clinical_df))

cat("\nFirst 6 rows:\n")
print(head(clinical_df))

############################################################
# 7. Detect patient ID and survival columns
############################################################

id_candidates <- c(
  "PatientID",
  "patient",
  "Patient",
  "Patient.Name",
  "PatientName",
  "patient_id",
  "id",
  "ID",
  "case",
  "Case"
)

id_hit <- intersect(id_candidates, colnames(clinical_df))

if (length(id_hit) > 0) {
  patient_id_col <- id_hit[1]
} else {
  
  id_content_hits <- sapply(colnames(clinical_df), function(cc) {
    any(grepl("^LUNG1-[0-9]+$", toupper(trimws(as.character(clinical_df[[cc]])))))
  })
  
  if (any(id_content_hits)) {
    patient_id_col <- names(id_content_hits)[which(id_content_hits)[1]]
  } else {
    patient_id_col <- NA_character_
  }
}

cat("\nDetected patient ID column:\n")
print(patient_id_col)

surv_candidates <- grep(
  "surv|time|death|event|status|days|os|dead",
  colnames(clinical_df),
  ignore.case = TRUE,
  value = TRUE
)

cat("\nSurvival-related candidate columns:\n")
print(surv_candidates)

############################################################
# 8. Standardize patient ID if possible
############################################################

clinical_clean <- clinical_df

if (!is.na(patient_id_col)) {
  clinical_clean$patient_id <- toupper(trimws(as.character(clinical_clean[[patient_id_col]])))
} else {
  warning("Patient ID column not detected automatically.")
}

############################################################
# 9. Save audit outputs
############################################################

write.csv(
  clinical_df,
  file = file.path(clinical_dir, "Lung1_clinical_raw_copy.csv"),
  row.names = FALSE
)

write.csv(
  clinical_clean,
  file = file.path(clinical_dir, "Lung1_clinical_clean_initial.csv"),
  row.names = FALSE
)

saveRDS(
  clinical_clean,
  file = file.path(clinical_dir, "Lung1_clinical_clean_initial.rds")
)

audit_summary <- data.frame(
  item = c(
    "lung1_root",
    "clinical_csv_path",
    "n_rows",
    "n_columns",
    "detected_patient_id_column",
    "survival_candidate_columns",
    "dicom_download_dir",
    "manifest_dir",
    "nrrd_dir"
  ),
  value = c(
    lung1_root,
    clinical_csv,
    nrow(clinical_df),
    ncol(clinical_df),
    patient_id_col,
    paste(surv_candidates, collapse = "; "),
    dicom_dir,
    manifest_dir,
    nrrd_dir
  ),
  stringsAsFactors = FALSE
)

write.csv(
  audit_summary,
  file = file.path(metadata_dir, "lung1_clinical_audit_summary.csv"),
  row.names = FALSE
)

wb <- createWorkbook()

addWorksheet(wb, "audit_summary")
writeData(wb, "audit_summary", audit_summary)

addWorksheet(wb, "clinical_raw")
writeData(wb, "clinical_raw", clinical_df)

addWorksheet(wb, "clinical_clean")
writeData(wb, "clinical_clean", clinical_clean)

saveWorkbook(
  wb,
  file = file.path(metadata_dir, "lung1_clinical_audit.xlsx"),
  overwrite = TRUE
)

############################################################
# 10. Done
############################################################

cat("\n===== Lung1 clinical preparation completed =====\n")
cat("Main output files:\n")

cat("\nFor NBIA Data Retriever, use this D-drive download folder:\n")

cat("1. Clinical table dim\n")
cat("2. Clinical columns\n")
cat("3. Detected patient ID column\n")
cat("4. Survival-related candidate columns\n")
