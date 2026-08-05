# 02_query_tcia_metadata.R
#
# Query TCIA series metadata for the NSCLC-Radiogenomics cohort.
#
# Run from an analysis workspace configured through environment variables.

options(stringsAsFactors = FALSE)
options(timeout = 1800)


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

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Missing R package 'jsonlite'. Install dependencies with environment/install_r_dependencies.R.")
}

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Missing R package 'openxlsx'. Install dependencies with environment/install_r_dependencies.R.")
}

library(jsonlite)
library(openxlsx)

############################################################
# 2. Folders
############################################################

dirs <- c(
  "02_metadata",
  "07_results",
  "10_logs"
)

for (d in dirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("\n===== Current working directory =====\n")
print(getwd())

############################################################
# 3. Load molecular ID table from the molecular patient-ID step
############################################################

id_path <- "02_metadata/nsclc_radiogenomics_molecular_patient_ids.rds"

if (!file.exists(id_path)) {
  stop(
    "Cannot find molecular ID table:\n",
    id_path,
    "\nRun the patient-ID mapping step first."
  )
}

mol_id <- readRDS(id_path)

standardize_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- toupper(x)
  return(x)
}

mol_id$patient_id_molecular <- standardize_id(mol_id$patient_id_molecular)

cat("\n===== Molecular ID table loaded =====\n")
print(dim(mol_id))

cat("\nMolecular patient ID preview:\n")
print(head(mol_id$patient_id_molecular, 30))

############################################################
# 4. Query TCIA/NBIA series metadata
############################################################

collection_name <- "NSCLC-Radiogenomics"

series_urls <- c(
  paste0(
    "https://services.cancerimagingarchive.net/nbia-api/services/v1/getSeries?",
    "Collection=", URLencode(collection_name, reserved = TRUE)
  ),
  paste0(
    "https://services.cancerimagingarchive.net/services/v4/TCIA/query/getSeries?",
    "Collection=", URLencode(collection_name, reserved = TRUE),
    "&format=json"
  ),
  paste0(
    "https://services.cancerimagingarchive.net/services/TCIA/TCIA/query/getSeries?",
    "Collection=", URLencode(collection_name, reserved = TRUE),
    "&format=json"
  )
)

query_success <- FALSE
series_df <- NULL
used_url <- NA
error_log <- character(0)

cat("\n===== Trying TCIA/NBIA API URLs =====\n")

for (u in series_urls) {
  cat("\nTrying URL:\n")
  cat(u, "\n")
  
  tmp <- tryCatch(
    {
      jsonlite::fromJSON(u, flatten = TRUE)
    },
    error = function(e) {
      error_log <<- c(error_log, paste("FAILED:", u, "ERROR:", e$message))
      return(NULL)
    }
  )
  
  if (!is.null(tmp)) {
    tmp_df <- as.data.frame(tmp, stringsAsFactors = FALSE)
    
    if (nrow(tmp_df) > 0 && ncol(tmp_df) > 0) {
      series_df <- tmp_df
      used_url <- u
      query_success <- TRUE
      break
    }
  }
}

if (!query_success) {
  cat("\n===== TCIA API query failed =====\n")
  print(error_log)
  
  stop(
    "The TCIA/NBIA metadata query failed.\n",
    "Check network, proxy, and institutional firewall settings.\n",
    "Review the error log above for the underlying request failure."
  )
}

cat("\n===== TCIA/NBIA query succeeded =====\n")
cat("Used URL:\n")
cat(used_url, "\n")

cat("\n===== Raw TCIA series metadata dim =====\n")
print(dim(series_df))

cat("\n===== Raw TCIA series metadata columns =====\n")
print(colnames(series_df))

cat("\n===== Raw TCIA series metadata preview =====\n")
print(head(series_df, 20))

############################################################
# 5. Flexible column detection
############################################################

find_col <- function(df_names, candidates) {
  for (cc in candidates) {
    hit <- which(toupper(df_names) == toupper(cc))
    if (length(hit) > 0) {
      return(df_names[hit[1]])
    }
  }
  
  # loose matching
  for (cc in candidates) {
    hit <- grep(toupper(cc), toupper(df_names), fixed = TRUE)
    if (length(hit) > 0) {
      return(df_names[hit[1]])
    }
  }
  
  return(NA)
}

patient_col <- find_col(
  colnames(series_df),
  c("PatientID", "Patient ID", "Subject ID", "SubjectID")
)

modality_col <- find_col(
  colnames(series_df),
  c("Modality")
)

study_uid_col <- find_col(
  colnames(series_df),
  c("StudyInstanceUID", "Study UID", "StudyInstanceUid")
)

series_uid_col <- find_col(
  colnames(series_df),
  c("SeriesInstanceUID", "Series ID", "SeriesInstanceUid")
)

series_desc_col <- find_col(
  colnames(series_df),
  c("SeriesDescription", "Series Description")
)

body_part_col <- find_col(
  colnames(series_df),
  c("BodyPartExamined", "Body Part Examined")
)

manufacturer_col <- find_col(
  colnames(series_df),
  c("Manufacturer")
)

num_images_col <- find_col(
  colnames(series_df),
  c("NumberOfImages", "Number of images", "ImageCount", "Image Count")
)

cat("\n===== Detected important TCIA metadata columns =====\n")
detected_cols <- data.frame(
  role = c(
    "patient_col",
    "modality_col",
    "study_uid_col",
    "series_uid_col",
    "series_desc_col",
    "body_part_col",
    "manufacturer_col",
    "num_images_col"
  ),
  column_name = c(
    patient_col,
    modality_col,
    study_uid_col,
    series_uid_col,
    series_desc_col,
    body_part_col,
    manufacturer_col,
    num_images_col
  ),
  stringsAsFactors = FALSE
)

print(detected_cols)

if (is.na(patient_col) || is.na(modality_col) || is.na(series_uid_col)) {
  stop(
    "Required columns were not detected: PatientID, Modality, or SeriesInstanceUID.\n",
    "Inspect the raw TCIA series metadata column names."
  )
}

############################################################
# 6. Build cleaned image metadata table
############################################################

get_col_or_na <- function(df, colname) {
  if (is.na(colname)) {
    return(rep(NA, nrow(df)))
  } else {
    return(df[[colname]])
  }
}

image_series <- data.frame(
  patient_id_tcia = standardize_id(series_df[[patient_col]]),
  modality = standardize_id(series_df[[modality_col]]),
  study_uid = as.character(get_col_or_na(series_df, study_uid_col)),
  series_uid = as.character(get_col_or_na(series_df, series_uid_col)),
  series_description = as.character(get_col_or_na(series_df, series_desc_col)),
  body_part_examined = as.character(get_col_or_na(series_df, body_part_col)),
  manufacturer = as.character(get_col_or_na(series_df, manufacturer_col)),
  number_of_images = suppressWarnings(as.numeric(get_col_or_na(series_df, num_images_col))),
  stringsAsFactors = FALSE
)

image_series <- unique(image_series)

cat("\n===== Cleaned image series table dim =====\n")
print(dim(image_series))

cat("\n===== Cleaned image series preview =====\n")
print(head(image_series, 30))

############################################################
# 7. Modality and patient audit
############################################################

modality_summary <- as.data.frame(table(image_series$modality), stringsAsFactors = FALSE)
colnames(modality_summary) <- c("modality", "series_count")
modality_summary <- modality_summary[order(-modality_summary$series_count), ]

patient_modality_table <- data.frame(
  patient_id_tcia = sort(unique(image_series$patient_id_tcia)),
  stringsAsFactors = FALSE
)

patient_modality_table$has_CT <- patient_modality_table$patient_id_tcia %in%
  unique(image_series$patient_id_tcia[image_series$modality == "CT"])

patient_modality_table$has_SEG <- patient_modality_table$patient_id_tcia %in%
  unique(image_series$patient_id_tcia[image_series$modality == "SEG"])

patient_modality_table$has_PT <- patient_modality_table$patient_id_tcia %in%
  unique(image_series$patient_id_tcia[image_series$modality == "PT"])

patient_modality_table$n_CT_series <- sapply(
  patient_modality_table$patient_id_tcia,
  function(pid) sum(image_series$patient_id_tcia == pid & image_series$modality == "CT")
)

patient_modality_table$n_SEG_series <- sapply(
  patient_modality_table$patient_id_tcia,
  function(pid) sum(image_series$patient_id_tcia == pid & image_series$modality == "SEG")
)

cat("\n===== TCIA modality summary =====\n")
print(modality_summary)

cat("\n===== TCIA patient modality table preview =====\n")
print(head(patient_modality_table, 30))

############################################################
# 8. Match molecular patients with TCIA image patients
############################################################

match_table <- merge(
  mol_id,
  patient_modality_table,
  by.x = "patient_id_molecular",
  by.y = "patient_id_tcia",
  all.x = TRUE,
  sort = FALSE
)

match_table$in_TCIA_metadata <- !is.na(match_table$has_CT)
match_table$has_CT[is.na(match_table$has_CT)] <- FALSE
match_table$has_SEG[is.na(match_table$has_SEG)] <- FALSE
match_table$has_PT[is.na(match_table$has_PT)] <- FALSE
match_table$n_CT_series[is.na(match_table$n_CT_series)] <- 0
match_table$n_SEG_series[is.na(match_table$n_SEG_series)] <- 0

match_table$eligible_for_CT_radiomics_next <- match_table$in_TCIA_metadata & match_table$has_CT

match_summary <- data.frame(
  item = c(
    "molecular_patients_total",
    "unique_TCIA_patients_in_collection",
    "molecular_patients_found_in_TCIA_metadata",
    "molecular_patients_with_CT",
    "molecular_patients_with_SEG",
    "molecular_patients_with_CT_and_SEG",
    "molecular_patients_missing_from_TCIA_metadata",
    "molecular_patients_without_CT"
  ),
  value = c(
    nrow(mol_id),
    length(unique(image_series$patient_id_tcia)),
    sum(match_table$in_TCIA_metadata),
    sum(match_table$has_CT),
    sum(match_table$has_SEG),
    sum(match_table$has_CT & match_table$has_SEG),
    sum(!match_table$in_TCIA_metadata),
    sum(!match_table$has_CT)
  ),
  stringsAsFactors = FALSE
)

cat("\n===== Molecular-image matching summary =====\n")
print(match_summary)

cat("\n===== Matched molecular-image table preview =====\n")
print(head(match_table, 40))

cat("\n===== Molecular patients missing from TCIA metadata =====\n")
print(match_table$patient_id_molecular[!match_table$in_TCIA_metadata])

cat("\n===== Molecular patients without CT =====\n")
print(match_table$patient_id_molecular[!match_table$has_CT])

############################################################
# 9. Candidate CT and SEG series for matched patients
############################################################

matched_patient_ids <- match_table$patient_id_molecular[match_table$has_CT]

candidate_ct_series <- image_series[
  image_series$patient_id_tcia %in% matched_patient_ids &
    image_series$modality == "CT",
]

candidate_seg_series <- image_series[
  image_series$patient_id_tcia %in% matched_patient_ids &
    image_series$modality == "SEG",
]

candidate_ct_series <- candidate_ct_series[
  order(candidate_ct_series$patient_id_tcia, -candidate_ct_series$number_of_images),
]

candidate_seg_series <- candidate_seg_series[
  order(candidate_seg_series$patient_id_tcia),
]

cat("\n===== Candidate CT series preview =====\n")
print(head(candidate_ct_series, 50))

cat("\n===== Candidate SEG series preview =====\n")
print(head(candidate_seg_series, 50))

############################################################
# 10. Save outputs
############################################################

write.csv(
  image_series,
  file = "02_metadata/tcia_all_series_metadata.csv",
  row.names = FALSE
)

write.csv(
  patient_modality_table,
  file = "02_metadata/tcia_patient_modality_summary.csv",
  row.names = FALSE
)

write.csv(
  match_table,
  file = "02_metadata/molecular_imaging_match_table.csv",
  row.names = FALSE
)

write.csv(
  candidate_ct_series,
  file = "02_metadata/candidate_ct_series.csv",
  row.names = FALSE
)

write.csv(
  candidate_seg_series,
  file = "02_metadata/candidate_seg_series.csv",
  row.names = FALSE
)

saveRDS(
  image_series,
  file = "02_metadata/tcia_all_series_metadata.rds"
)

saveRDS(
  match_table,
  file = "02_metadata/molecular_imaging_match_table.rds"
)

wb <- createWorkbook()

addWorksheet(wb, "all_series_metadata")
writeData(wb, "all_series_metadata", image_series)

addWorksheet(wb, "modality_summary")
writeData(wb, "modality_summary", modality_summary)

addWorksheet(wb, "patient_modality_table")
writeData(wb, "patient_modality_table", patient_modality_table)

addWorksheet(wb, "match_summary")
writeData(wb, "match_summary", match_summary)

addWorksheet(wb, "molecular_image_match")
writeData(wb, "molecular_image_match", match_table)

addWorksheet(wb, "candidate_CT_series")
writeData(wb, "candidate_CT_series", candidate_ct_series)

addWorksheet(wb, "candidate_SEG_series")
writeData(wb, "candidate_SEG_series", candidate_seg_series)

addWorksheet(wb, "detected_columns")
writeData(wb, "detected_columns", detected_cols)

saveWorkbook(
  wb,
  file = "02_metadata/tcia_image_matching_summary.xlsx",
  overwrite = TRUE
)

cat("\n===== TCIA image metadata matching completed =====\n")
cat("Main output files:\n")
cat("1. 02_metadata/tcia_image_matching_summary.xlsx\n")
cat("2. 02_metadata/molecular_imaging_match_table.csv\n")
cat("3. 02_metadata/candidate_ct_series.csv\n")
cat("4. 02_metadata/candidate_seg_series.csv\n")

cat("1. Raw TCIA series metadata dim\n")
cat("2. Detected important TCIA metadata columns\n")
cat("3. TCIA modality summary\n")
cat("4. Molecular-image matching summary\n")
cat("5. Molecular patients missing from TCIA metadata\n")
cat("6. Candidate CT series preview\n")
cat("7. Candidate SEG series preview\n")
