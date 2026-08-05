# 04_build_ct_seg_pairs.R
#
# Build candidate CT and segmentation series pairs.
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

if (!requireNamespace("oro.dicom", quietly = TRUE)) {
  stop("Missing R package 'oro.dicom'. Install dependencies with environment/install_r_dependencies.R.")
}

library(openxlsx)
library(oro.dicom)

############################################################
# 2. Paths
############################################################

base_dir <- "01_raw_data/NSCLC_Radiogenomics/TCIA_DICOM"
metadata_csv <- file.path(base_dir, "metadata", "metadata.csv")

path_table_rds <- "02_metadata/dicom_path_table.rds"
patient_counts_rds <- "02_metadata/downloaded_patient_dicom_counts.rds"

out_dir <- "02_metadata"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(path_table_rds)) {
  stop("dicom_path_table.rds was not found. Run the DICOM download audit first.")
}

if (!file.exists(metadata_csv)) {
  stop("Cannot find metadata.csv in TCIA_DICOM/metadata.")
}

############################################################
# 3. Load DICOM download-audit outputs and TCIA metadata.csv
############################################################

dicom_path_table <- readRDS(path_table_rds)
patient_counts <- readRDS(patient_counts_rds)

metadata_df <- read.csv(
  metadata_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\n===== Loaded inputs =====\n")
cat("DICOM path table dim:\n")
print(dim(dicom_path_table))

cat("Patient counts dim:\n")
print(dim(patient_counts))

cat("metadata.csv dim:\n")
print(dim(metadata_df))

cat("\nmetadata.csv columns:\n")
print(colnames(metadata_df))

############################################################
# 4. Build local series table
############################################################

series_table <- aggregate(
  dicom_file ~ patient_id_downloaded + study_dir + series_dir,
  data = dicom_path_table,
  FUN = length
)

colnames(series_table)[colnames(series_table) == "dicom_file"] <- "n_dicom_files"

first_file_table <- aggregate(
  dicom_file ~ patient_id_downloaded + study_dir + series_dir,
  data = dicom_path_table,
  FUN = function(x) x[1]
)

colnames(first_file_table)[colnames(first_file_table) == "dicom_file"] <- "first_dicom_file"

series_table <- merge(
  series_table,
  first_file_table,
  by = c("patient_id_downloaded", "study_dir", "series_dir"),
  all.x = TRUE
)

series_table <- series_table[order(series_table$patient_id_downloaded, series_table$study_dir, series_table$series_dir), ]

cat("\n===== Local series table =====\n")
print(dim(series_table))
print(head(series_table, 20))

############################################################
# 5. Normalize metadata path and merge with local series table
############################################################

metadata_df$S5cmdManifestPath_norm <- gsub("\\\\", "/", metadata_df$S5cmdManifestPath)

metadata_df$series_dir <- sub(
  "^.*TCIA_DICOM/",
  "",
  metadata_df$S5cmdManifestPath_norm
)

metadata_keep <- metadata_df[, intersect(
  c(
    "PatientID",
    "StudyInstanceUID",
    "SeriesInstanceUID",
    "Collection",
    "FileSize",
    "DownloadURL",
    "S5cmdManifestPath",
    "OriginalS5cmdURI",
    "completion_status",
    "series_dir"
  ),
  colnames(metadata_df)
)]

series_table <- merge(
  series_table,
  metadata_keep,
  by = "series_dir",
  all.x = TRUE
)

cat("\n===== Series table after merging metadata.csv =====\n")
print(dim(series_table))
print(head(series_table, 20))

cat("\nCompletion status table:\n")
print(table(series_table$completion_status, useNA = "ifany"))

############################################################
# 6. Helper functions for DICOM header extraction
############################################################

clean_text <- function(x) {
  if (length(x) == 0 || is.null(x)) return(NA_character_)
  x <- as.character(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  x <- paste(unique(x), collapse = " | ")
  x <- gsub("[\r\n\t]", " ", x)
  x <- gsub(" +", " ", x)
  trimws(x)
}

normalize_header_name <- function(x) {
  x <- tolower(as.character(x))
  gsub("[^a-z0-9]", "", x)
}

get_dicom_value <- function(hdr, possible_names) {
  if (is.null(hdr)) return(NA_character_)
  
  hdr <- as.data.frame(hdr, stringsAsFactors = FALSE)
  
  name_col <- intersect(
    c("name", "Name", "description", "Description", "field", "Field"),
    colnames(hdr)
  )
  
  value_col <- intersect(
    c("value", "Value", "val", "Val"),
    colnames(hdr)
  )
  
  if (length(name_col) == 0 || length(value_col) == 0) {
    return(NA_character_)
  }
  
  name_col <- name_col[1]
  value_col <- value_col[1]
  
  header_names_norm <- normalize_header_name(hdr[[name_col]])
  targets_norm <- normalize_header_name(possible_names)
  
  idx <- which(header_names_norm %in% targets_norm)
  
  if (length(idx) == 0) {
    idx <- which(sapply(targets_norm, function(p) any(grepl(p, header_names_norm, fixed = TRUE))))
    if (length(idx) > 0) {
      matched_target <- targets_norm[idx[1]]
      idx <- which(grepl(matched_target, header_names_norm, fixed = TRUE))
    }
  }
  
  if (length(idx) == 0) return(NA_character_)
  
  clean_text(hdr[[value_col]][idx])
}

read_one_dicom_header <- function(file) {
  cat("Reading header:", file, "\n")
  
  obj <- tryCatch(
    suppressWarnings(oro.dicom::readDICOMFile(file, pixelData = FALSE)),
    error = function(e1) {
      tryCatch(
        suppressWarnings(oro.dicom::readDICOMFile(file)),
        error = function(e2) {
          return(list(error_message = paste(e1$message, " | ", e2$message)))
        }
      )
    }
  )
  
  if (!is.null(obj$error_message)) {
    return(data.frame(
      header_read_success = FALSE,
      header_error = obj$error_message,
      Modality = NA,
      SeriesDescription = NA,
      StudyDescription = NA,
      BodyPartExamined = NA,
      ImageType = NA,
      SOPClassUID = NA,
      PatientID_header = NA,
      StudyInstanceUID_header = NA,
      SeriesInstanceUID_header = NA,
      StudyDate = NA,
      SeriesDate = NA,
      SeriesNumber = NA,
      InstanceNumber = NA,
      Rows = NA,
      Columns = NA,
      PixelSpacing = NA,
      SliceThickness = NA,
      Manufacturer = NA,
      ConvolutionKernel = NA,
      KVP = NA,
      XRayTubeCurrent = NA,
      ExposureTime = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  hdr <- obj$hdr
  
  data.frame(
    header_read_success = TRUE,
    header_error = NA,
    Modality = get_dicom_value(hdr, c("Modality")),
    SeriesDescription = get_dicom_value(hdr, c("SeriesDescription", "Series Description")),
    StudyDescription = get_dicom_value(hdr, c("StudyDescription", "Study Description")),
    BodyPartExamined = get_dicom_value(hdr, c("BodyPartExamined", "Body Part Examined")),
    ImageType = get_dicom_value(hdr, c("ImageType", "Image Type")),
    SOPClassUID = get_dicom_value(hdr, c("SOPClassUID", "SOP Class UID")),
    PatientID_header = get_dicom_value(hdr, c("PatientID", "Patient ID")),
    StudyInstanceUID_header = get_dicom_value(hdr, c("StudyInstanceUID", "Study Instance UID")),
    SeriesInstanceUID_header = get_dicom_value(hdr, c("SeriesInstanceUID", "Series Instance UID")),
    StudyDate = get_dicom_value(hdr, c("StudyDate", "Study Date")),
    SeriesDate = get_dicom_value(hdr, c("SeriesDate", "Series Date")),
    SeriesNumber = get_dicom_value(hdr, c("SeriesNumber", "Series Number")),
    InstanceNumber = get_dicom_value(hdr, c("InstanceNumber", "Instance Number")),
    Rows = get_dicom_value(hdr, c("Rows")),
    Columns = get_dicom_value(hdr, c("Columns")),
    PixelSpacing = get_dicom_value(hdr, c("PixelSpacing", "Pixel Spacing")),
    SliceThickness = get_dicom_value(hdr, c("SliceThickness", "Slice Thickness")),
    Manufacturer = get_dicom_value(hdr, c("Manufacturer")),
    ConvolutionKernel = get_dicom_value(hdr, c("ConvolutionKernel", "Convolution Kernel")),
    KVP = get_dicom_value(hdr, c("KVP")),
    XRayTubeCurrent = get_dicom_value(hdr, c("XRayTubeCurrent", "X Ray Tube Current", "X-Ray Tube Current")),
    ExposureTime = get_dicom_value(hdr, c("ExposureTime", "Exposure Time")),
    stringsAsFactors = FALSE
  )
}

############################################################
# 7. Read one DICOM header per series
############################################################

cat("\n===== Reading representative DICOM header for each series =====\n")
cat("Total series to read:\n")
print(nrow(series_table))

header_list <- vector("list", nrow(series_table))

for (i in seq_len(nrow(series_table))) {
  cat("\n--- Series", i, "of", nrow(series_table), "---\n")
  
  header_list[[i]] <- read_one_dicom_header(series_table$first_dicom_file[i])
  
  if (i %% 20 == 0) {
    cat("\nProgress:", i, "/", nrow(series_table), "\n")
  }
}

header_df <- do.call(rbind, header_list)

series_inventory <- cbind(series_table, header_df)

cat("\n===== Header read success table =====\n")
print(table(series_inventory$header_read_success, useNA = "ifany"))

cat("\n===== Modality table =====\n")
print(table(series_inventory$Modality, useNA = "ifany"))

cat("\n===== Series inventory preview =====\n")
print(head(series_inventory[, c(
  "patient_id_downloaded",
  "n_dicom_files",
  "Modality",
  "SeriesDescription",
  "StudyDescription",
  "BodyPartExamined",
  "ImageType",
  "series_dir"
)], 50))

############################################################
# 8. Classify CT and SEG series
############################################################

modality_upper <- toupper(trimws(series_inventory$Modality))
series_desc_upper <- toupper(paste(
  series_inventory$SeriesDescription,
  series_inventory$StudyDescription,
  series_inventory$ImageType
))

series_inventory$is_CT <- modality_upper == "CT"

series_inventory$is_SEG <- modality_upper == "SEG" |
  grepl("SEGMENT", series_desc_upper) |
  grepl("SEGMENTATION", series_desc_upper)

series_inventory$series_class <- "OTHER"
series_inventory$series_class[series_inventory$is_CT] <- "CT"
series_inventory$series_class[series_inventory$is_SEG] <- "SEG"

cat("\n===== Series class table =====\n")
print(table(series_inventory$series_class, useNA = "ifany"))

############################################################
# 9. Patient-level CT/SEG availability
############################################################

all_patients <- sort(unique(series_inventory$patient_id_downloaded))

patient_availability <- data.frame(
  patient_id = all_patients,
  stringsAsFactors = FALSE
)

patient_availability$n_total_series <- sapply(all_patients, function(pid) {
  sum(series_inventory$patient_id_downloaded == pid)
})

patient_availability$n_CT_series <- sapply(all_patients, function(pid) {
  sum(series_inventory$patient_id_downloaded == pid & series_inventory$is_CT, na.rm = TRUE)
})

patient_availability$n_SEG_series <- sapply(all_patients, function(pid) {
  sum(series_inventory$patient_id_downloaded == pid & series_inventory$is_SEG, na.rm = TRUE)
})

patient_availability$n_OTHER_series <- sapply(all_patients, function(pid) {
  sum(series_inventory$patient_id_downloaded == pid & series_inventory$series_class == "OTHER", na.rm = TRUE)
})

patient_availability$n_dicom_files_total <- sapply(all_patients, function(pid) {
  sum(series_inventory$n_dicom_files[series_inventory$patient_id_downloaded == pid], na.rm = TRUE)
})

patient_availability$has_CT <- patient_availability$n_CT_series > 0
patient_availability$has_SEG <- patient_availability$n_SEG_series > 0
patient_availability$has_CT_and_SEG <- patient_availability$has_CT & patient_availability$has_SEG

cat("\n===== Patient-level CT/SEG availability summary =====\n")
print(table(patient_availability$has_CT, patient_availability$has_SEG))

cat("\nNumber of patients with both CT and SEG:\n")
print(sum(patient_availability$has_CT_and_SEG))

cat("\nPatient availability preview:\n")
print(head(patient_availability, 80))

############################################################
# 10. Candidate CT-SEG pairing
############################################################

candidate_pairs <- data.frame()

for (pid in all_patients) {
  sub <- series_inventory[series_inventory$patient_id_downloaded == pid, ]
  
  ct_sub <- sub[sub$is_CT, ]
  seg_sub <- sub[sub$is_SEG, ]
  
  if (nrow(ct_sub) == 0 || nrow(seg_sub) == 0) {
    one <- data.frame(
      patient_id = pid,
      pair_status = "NO_CT_OR_NO_SEG",
      selected_CT_series_dir = NA,
      selected_CT_StudyInstanceUID = NA,
      selected_CT_SeriesInstanceUID = NA,
      selected_CT_n_dicom_files = NA,
      selected_CT_SeriesDescription = NA,
      selected_SEG_series_dir = NA,
      selected_SEG_StudyInstanceUID = NA,
      selected_SEG_SeriesInstanceUID = NA,
      selected_SEG_n_dicom_files = NA,
      selected_SEG_SeriesDescription = NA,
      pairing_rule = NA,
      stringsAsFactors = FALSE
    )
    
    candidate_pairs <- rbind(candidate_pairs, one)
    next
  }
  
  seg_sub <- seg_sub[order(-seg_sub$n_dicom_files), ]
  selected_seg <- seg_sub[1, ]
  
  ct_same_study <- ct_sub[ct_sub$StudyInstanceUID == selected_seg$StudyInstanceUID, ]
  
  if (nrow(ct_same_study) > 0) {
    ct_same_study <- ct_same_study[order(-ct_same_study$n_dicom_files), ]
    selected_ct <- ct_same_study[1, ]
    pairing_rule <- "largest_CT_in_same_study_as_selected_SEG"
  } else {
    ct_sub <- ct_sub[order(-ct_sub$n_dicom_files), ]
    selected_ct <- ct_sub[1, ]
    pairing_rule <- "largest_CT_in_patient_no_same_study_match"
  }
  
  one <- data.frame(
    patient_id = pid,
    pair_status = "CANDIDATE_PAIR_CREATED",
    selected_CT_series_dir = selected_ct$series_dir,
    selected_CT_StudyInstanceUID = selected_ct$StudyInstanceUID,
    selected_CT_SeriesInstanceUID = selected_ct$SeriesInstanceUID,
    selected_CT_n_dicom_files = selected_ct$n_dicom_files,
    selected_CT_SeriesDescription = selected_ct$SeriesDescription,
    selected_SEG_series_dir = selected_seg$series_dir,
    selected_SEG_StudyInstanceUID = selected_seg$StudyInstanceUID,
    selected_SEG_SeriesInstanceUID = selected_seg$SeriesInstanceUID,
    selected_SEG_n_dicom_files = selected_seg$n_dicom_files,
    selected_SEG_SeriesDescription = selected_seg$SeriesDescription,
    pairing_rule = pairing_rule,
    stringsAsFactors = FALSE
  )
  
  candidate_pairs <- rbind(candidate_pairs, one)
}

cat("\n===== Candidate CT-SEG pairing status =====\n")
print(table(candidate_pairs$pair_status, useNA = "ifany"))

cat("\n===== Pairing rule table =====\n")
print(table(candidate_pairs$pairing_rule, useNA = "ifany"))

cat("\n===== Candidate pairs preview =====\n")
print(head(candidate_pairs, 50))

############################################################
# 11. Merge with molecular score table if available
############################################################

score_path <- "05_molecular_scores/NSCLC_Radiogenomics_patient_molecular_scores.rds"

if (file.exists(score_path)) {
  mol_scores <- readRDS(score_path)
  
  if ("patient_id" %in% colnames(mol_scores)) {
    patient_availability <- merge(
      patient_availability,
      mol_scores,
      by.x = "patient_id",
      by.y = "patient_id",
      all.x = TRUE
    )
    
    candidate_pairs <- merge(
      candidate_pairs,
      mol_scores,
      by.x = "patient_id",
      by.y = "patient_id",
      all.x = TRUE
    )
    
    cat("\nMolecular scores merged into patient availability and candidate pairs.\n")
  } else {
    cat("\nMolecular score file found, but no patient_id column detected. Skipping merge.\n")
  }
} else {
  cat("\nMolecular score file not found. Skipping merge.\n")
}

############################################################
# 12. Save outputs
############################################################

write.csv(
  series_inventory,
  file = "02_metadata/ct_seg_series_inventory.csv",
  row.names = FALSE
)

write.csv(
  patient_availability,
  file = "02_metadata/patient_ct_seg_availability.csv",
  row.names = FALSE
)

write.csv(
  candidate_pairs,
  file = "02_metadata/ct_seg_candidate_pairs.csv",
  row.names = FALSE
)

saveRDS(
  series_inventory,
  file = "02_metadata/ct_seg_series_inventory.rds"
)

saveRDS(
  patient_availability,
  file = "02_metadata/patient_ct_seg_availability.rds"
)

saveRDS(
  candidate_pairs,
  file = "02_metadata/ct_seg_candidate_pairs.rds"
)

wb <- createWorkbook()

addWorksheet(wb, "series_inventory")
writeData(wb, "series_inventory", series_inventory)

addWorksheet(wb, "patient_availability")
writeData(wb, "patient_availability", patient_availability)

addWorksheet(wb, "candidate_pairs")
writeData(wb, "candidate_pairs", candidate_pairs)

addWorksheet(wb, "modality_table")
modality_table <- as.data.frame(table(series_inventory$Modality, useNA = "ifany"))
colnames(modality_table) <- c("Modality", "n_series")
writeData(wb, "modality_table", modality_table)

addWorksheet(wb, "series_class_table")
series_class_table <- as.data.frame(table(series_inventory$series_class, useNA = "ifany"))
colnames(series_class_table) <- c("series_class", "n_series")
writeData(wb, "series_class_table", series_class_table)

addWorksheet(wb, "availability_summary")
availability_summary <- data.frame(
  item = c(
    "patients_total",
    "patients_with_CT",
    "patients_with_SEG",
    "patients_with_CT_and_SEG",
    "total_series",
    "CT_series",
    "SEG_series",
    "OTHER_series"
  ),
  value = c(
    nrow(patient_availability),
    sum(patient_availability$has_CT),
    sum(patient_availability$has_SEG),
    sum(patient_availability$has_CT_and_SEG),
    nrow(series_inventory),
    sum(series_inventory$is_CT, na.rm = TRUE),
    sum(series_inventory$is_SEG, na.rm = TRUE),
    sum(series_inventory$series_class == "OTHER", na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
writeData(wb, "availability_summary", availability_summary)

saveWorkbook(
  wb,
  file = "02_metadata/ct_seg_series_inventory_summary.xlsx",
  overwrite = TRUE
)

cat("\n===== CT-SEG series metadata audit completed =====\n")
cat("Main output files:\n")
cat("1. 02_metadata/ct_seg_series_inventory_summary.xlsx\n")
cat("2. 02_metadata/ct_seg_series_inventory.csv\n")
cat("3. 02_metadata/patient_ct_seg_availability.csv\n")
cat("4. 02_metadata/ct_seg_candidate_pairs.csv\n")

cat("1. Header read success table\n")
cat("2. Modality table\n")
cat("3. Series class table\n")
cat("4. Patient-level CT/SEG availability summary\n")
cat("5. Number of patients with both CT and SEG\n")
cat("6. Candidate CT-SEG pairing status\n")
cat("7. Pairing rule table\n")
