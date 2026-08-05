# 03_verify_dicom_download.R
#
# Audit the local NSCLC-Radiogenomics DICOM download.
#
# Run from an analysis workspace configured through environment variables.

options(stringsAsFactors = FALSE)


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

############################################################
# 1. Packages
############################################################

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Missing R package 'openxlsx'. Install dependencies with environment/install_r_dependencies.R.")
}

library(openxlsx)

############################################################
# 2. Paths
############################################################

base_dir <- "01_raw_data/NSCLC_Radiogenomics/TCIA_DICOM"
metadata_dir <- file.path(base_dir, "metadata")
collection_dir <- file.path(base_dir, "nsclc_radiogenomics")

out_dir <- "02_metadata"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n===== Current working directory =====\n")
print(getwd())

cat("\n===== DICOM base directory =====\n")
print(normalizePath(base_dir, winslash = "/", mustWork = FALSE))

if (!dir.exists(base_dir)) {
  stop("The DICOM base directory does not exist. Verify the configured download path.")
}

############################################################
# 3. Basic folder audit
############################################################

cat("\n===== Top-level folders in TCIA_DICOM =====\n")
top_folders <- list.files(base_dir, recursive = FALSE, full.names = FALSE)
print(top_folders)

cat("\nMetadata folder exists:\n")
print(dir.exists(metadata_dir))

cat("\nCollection folder exists:\n")
print(dir.exists(collection_dir))

############################################################
# 4. DICOM file listing
############################################################

dicom_files <- list.files(
  base_dir,
  pattern = "\\.dcm$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

cat("\n===== Number of DICOM files =====\n")
print(length(dicom_files))

if (length(dicom_files) == 0) {
  stop("No .dcm files were found. Verify that the NBIA Data Retriever download completed.")
}

base_norm <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
dicom_norm <- normalizePath(dicom_files, winslash = "/", mustWork = TRUE)

relative_paths <- ifelse(
  startsWith(dicom_norm, paste0(base_norm, "/")),
  substring(dicom_norm, nchar(base_norm) + 2),
  dicom_norm
)

cat("\n===== First 30 relative DICOM paths =====\n")
print(head(relative_paths, 30))

############################################################
# 5. Extract patient IDs from paths
############################################################

extract_r01_id <- function(x) {
  m <- regexpr("R01-[0-9]{3}", x, ignore.case = TRUE)
  out <- ifelse(m > 0, regmatches(x, m), NA)
  toupper(out)
}

patient_id_from_path <- extract_r01_id(relative_paths)

dicom_path_table <- data.frame(
  dicom_file = dicom_norm,
  relative_path = relative_paths,
  patient_id_downloaded = patient_id_from_path,
  series_dir = dirname(relative_paths),
  study_dir = dirname(dirname(relative_paths)),
  stringsAsFactors = FALSE
)

cat("\n===== DICOM path table preview =====\n")
print(head(dicom_path_table, 30))

cat("\nNumber of DICOM files with R01 patient ID in path:\n")
print(sum(!is.na(dicom_path_table$patient_id_downloaded)))

cat("\nNumber of DICOM files without R01 patient ID in path:\n")
print(sum(is.na(dicom_path_table$patient_id_downloaded)))

############################################################
# 6. Patient-level DICOM count
############################################################

patient_counts <- aggregate(
  dicom_file ~ patient_id_downloaded,
  data = dicom_path_table[!is.na(dicom_path_table$patient_id_downloaded), ],
  FUN = length
)

colnames(patient_counts) <- c("patient_id_downloaded", "n_dicom_files")

series_counts <- aggregate(
  series_dir ~ patient_id_downloaded,
  data = unique(dicom_path_table[!is.na(dicom_path_table$patient_id_downloaded), c("patient_id_downloaded", "series_dir")]),
  FUN = length
)

colnames(series_counts) <- c("patient_id_downloaded", "n_series_dirs")

study_counts <- aggregate(
  study_dir ~ patient_id_downloaded,
  data = unique(dicom_path_table[!is.na(dicom_path_table$patient_id_downloaded), c("patient_id_downloaded", "study_dir")]),
  FUN = length
)

colnames(study_counts) <- c("patient_id_downloaded", "n_study_dirs")

patient_counts <- merge(patient_counts, series_counts, by = "patient_id_downloaded", all.x = TRUE)
patient_counts <- merge(patient_counts, study_counts, by = "patient_id_downloaded", all.x = TRUE)

patient_counts <- patient_counts[order(patient_counts$patient_id_downloaded), ]

cat("\n===== Downloaded R01 patient count =====\n")
print(length(unique(patient_counts$patient_id_downloaded)))

cat("\n===== Patient-level DICOM count preview =====\n")
print(head(patient_counts, 50))

cat("\n===== Patient-level DICOM count summary =====\n")
print(summary(patient_counts$n_dicom_files))

############################################################
# 7. Compare with 130 molecular patients
############################################################

mol_id_path <- "02_metadata/nsclc_radiogenomics_molecular_patient_ids.rds"

if (!file.exists(mol_id_path)) {
  stop("Cannot find molecular ID table from the molecular patient-ID step.")
}

mol_id <- readRDS(mol_id_path)
mol_ids <- unique(toupper(mol_id$patient_id_molecular))
downloaded_ids <- unique(patient_counts$patient_id_downloaded)

missing_from_download <- setdiff(mol_ids, downloaded_ids)
downloaded_not_in_molecular <- setdiff(downloaded_ids, mol_ids)

match_summary <- data.frame(
  item = c(
    "molecular_patients_total",
    "downloaded_R01_patients",
    "molecular_patients_found_in_download",
    "molecular_patients_missing_from_download",
    "downloaded_R01_patients_not_in_molecular",
    "total_DICOM_files"
  ),
  value = c(
    length(mol_ids),
    length(downloaded_ids),
    length(intersect(mol_ids, downloaded_ids)),
    length(missing_from_download),
    length(downloaded_not_in_molecular),
    length(dicom_files)
  ),
  stringsAsFactors = FALSE
)

cat("\n===== Molecular vs downloaded image matching summary =====\n")
print(match_summary)

cat("\n===== Molecular patients missing from downloaded DICOM =====\n")
print(missing_from_download)

cat("\n===== Downloaded R01 patients not in molecular table =====\n")
print(downloaded_not_in_molecular)

############################################################
# 8. Metadata folder audit
############################################################

metadata_files <- character(0)

if (dir.exists(metadata_dir)) {
  metadata_files <- list.files(
    metadata_dir,
    recursive = TRUE,
    full.names = TRUE
  )
}

cat("\n===== Metadata files =====\n")
print(metadata_files)

metadata_file_table <- data.frame(
  file = metadata_files,
  file_name = basename(metadata_files),
  size_bytes = ifelse(length(metadata_files) > 0, file.info(metadata_files)$size, numeric(0)),
  stringsAsFactors = FALSE
)

cat("\n===== Metadata file table =====\n")
print(metadata_file_table)

############################################################
# 9. Try reading metadata CSV/TXT files if present
############################################################

metadata_text_files <- metadata_files[
  grepl("\\.(csv|tsv|txt)$", metadata_files, ignore.case = TRUE)
]

metadata_read_audit <- data.frame()

if (length(metadata_text_files) > 0) {
  for (ff in metadata_text_files) {
    cat("\n===== Inspecting metadata file =====\n")
    cat(ff, "\n")
    
    first_lines <- tryCatch(
      readLines(ff, n = 5, warn = FALSE),
      error = function(e) paste("READ ERROR:", e$message)
    )
    
    print(first_lines)
    
    sep_guess <- ifelse(grepl("\\.tsv$|\\.txt$", ff, ignore.case = TRUE), "\t", ",")
    
    tmp <- tryCatch(
      read.delim(
        ff,
        header = TRUE,
        sep = sep_guess,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      error = function(e) NULL
    )
    
    if (!is.null(tmp)) {
      one <- data.frame(
        file_name = basename(ff),
        nrow = nrow(tmp),
        ncol = ncol(tmp),
        columns = paste(colnames(tmp), collapse = ";"),
        stringsAsFactors = FALSE
      )
      
      metadata_read_audit <- rbind(metadata_read_audit, one)
    }
  }
}

cat("\n===== Metadata read audit =====\n")
print(metadata_read_audit)

############################################################
# 10. Save outputs
############################################################

write.csv(
  dicom_path_table,
  file = "02_metadata/dicom_path_table.csv",
  row.names = FALSE
)

write.csv(
  patient_counts,
  file = "02_metadata/downloaded_patient_dicom_counts.csv",
  row.names = FALSE
)

write.csv(
  match_summary,
  file = "02_metadata/molecular_download_matching_summary.csv",
  row.names = FALSE
)

write.csv(
  data.frame(patient_id = missing_from_download),
  file = "02_metadata/molecular_patients_missing_from_download.csv",
  row.names = FALSE
)

write.csv(
  metadata_file_table,
  file = "02_metadata/dicom_metadata_file_table.csv",
  row.names = FALSE
)

if (nrow(metadata_read_audit) > 0) {
  write.csv(
    metadata_read_audit,
    file = "02_metadata/dicom_metadata_read_audit.csv",
    row.names = FALSE
  )
}

saveRDS(
  dicom_path_table,
  file = "02_metadata/dicom_path_table.rds"
)

saveRDS(
  patient_counts,
  file = "02_metadata/downloaded_patient_dicom_counts.rds"
)

saveRDS(
  match_summary,
  file = "02_metadata/molecular_download_matching_summary.rds"
)

wb <- createWorkbook()

addWorksheet(wb, "match_summary")
writeData(wb, "match_summary", match_summary)

addWorksheet(wb, "patient_counts")
writeData(wb, "patient_counts", patient_counts)

addWorksheet(wb, "missing_from_download")
writeData(wb, "missing_from_download", data.frame(patient_id = missing_from_download))

addWorksheet(wb, "downloaded_not_molecular")
writeData(wb, "downloaded_not_molecular", data.frame(patient_id = downloaded_not_in_molecular))

addWorksheet(wb, "metadata_files")
writeData(wb, "metadata_files", metadata_file_table)

if (nrow(metadata_read_audit) > 0) {
  addWorksheet(wb, "metadata_read_audit")
  writeData(wb, "metadata_read_audit", metadata_read_audit)
}

addWorksheet(wb, "path_preview")
writeData(wb, "path_preview", head(dicom_path_table, 5000))

saveWorkbook(
  wb,
  file = "02_metadata/dicom_download_summary.xlsx",
  overwrite = TRUE
)

cat("\n===== DICOM download audit completed =====\n")
cat("Main output files:\n")
cat("1. 02_metadata/dicom_download_summary.xlsx\n")
cat("2. 02_metadata/dicom_path_table.csv\n")
cat("3. 02_metadata/downloaded_patient_dicom_counts.csv\n")
cat("4. 02_metadata/molecular_download_matching_summary.csv\n")

cat("1. Number of DICOM files\n")
cat("2. First 30 relative DICOM paths\n")
cat("3. Downloaded R01 patient count\n")
cat("4. Patient-level DICOM count preview\n")
cat("5. Molecular vs downloaded image matching summary\n")
cat("6. Molecular patients missing from downloaded DICOM\n")
cat("7. Metadata files\n")
cat("8. Metadata read audit\n")
