# 01_export_supplementary_tables.R
#
# Export source tables required for the supplementary material.
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

lung1_root <- required_env_path("LUNG1_DATA_DIR")

out_dir <- file.path(
  project_dir,
  "07_results",
  "supplementary_exports"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
# 1. Exact CRS nonzero-feature mapping for Supplementary S1
############################################################

crs_map_rds <- file.path(
  project_dir,
  "02_metadata",
  "crs_nonzero_feature_mapping.rds"
)

if (!file.exists(crs_map_rds)) {
  stop(
    "Missing exact CRS mapping file:\n",
    crs_map_rds
  )
}

crs_map <- readRDS(crs_map_rds)

required_crs_cols <- c(
  "ensembl", "coefficient", "symbol", "hgnc_id",
  "name", "locus_group", "locus_type", "symbol_in_patient"
)

missing_crs_cols <- setdiff(required_crs_cols, names(crs_map))
if (length(missing_crs_cols) > 0) {
  stop(
    "CRS mapping file is missing columns: ",
    paste(missing_crs_cols, collapse = ", ")
  )
}

crs_map <- crs_map[, required_crs_cols, drop = FALSE]

if (nrow(crs_map) != 112L) {
  stop("Expected 112 nonzero CRS features, found ", nrow(crs_map), ".")
}

write.csv(
  crs_map,
  file.path(out_dir, "crs_nonzero_feature_mapping.csv"),
  row.names = FALSE,
  na = ""
)

############################################################
# 2. Exact signed-log domain-shift tables for Supplementary S5
############################################################

signed_log_dir <- file.path(
  lung1_root,
  "signed_log_rrs"
)

shift_file <- file.path(
  signed_log_dir,
  "signed_log_domain_shift_features.csv"
)

nonzero_file <- file.path(
  signed_log_dir,
  "signed_log_model_nonzero_coefficients.csv"
)

extreme_file <- file.path(
  signed_log_dir,
  "signed_log_extreme_patients.csv"
)

for (ff in c(shift_file, nonzero_file, extreme_file)) {
  if (!file.exists(ff)) {
    stop(
      "Missing required Lung1 signed-log audit output:\n",
      ff,
      "\nRun R/04_lung1/05_train_signed_log_sensitivity_model.R first if needed."
    )
  }
}

shift_df <- read.csv(
  shift_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

nonzero_df <- read.csv(
  nonzero_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

extreme_df <- read.csv(
  extreme_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

write.csv(
  shift_df,
  file.path(out_dir, "signed_log_domain_shift_features.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  nonzero_df,
  file.path(out_dir, "signed_log_model_nonzero_coefficients.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  extreme_df,
  file.path(out_dir, "signed_log_extreme_patients.csv"),
  row.names = FALSE,
  na = ""
)

############################################################
# 3. Summary and zip package
############################################################

summary_lines <- c(
  "===== SUPPLEMENTARY SOURCE EXPORT SUMMARY =====",
  paste0("CRS nonzero mapping rows: ", nrow(crs_map)),
  paste0(
    "Mapped to HGNC symbol: ",
    sum(!is.na(crs_map$symbol) & nzchar(trimws(as.character(crs_map$symbol))))
  ),
  paste0(
    "Available in NSCLC-Radiogenomics: ",
    sum(as.logical(crs_map$symbol_in_patient), na.rm = TRUE)
  ),
  paste0("Signed-log domain-shift feature rows: ", nrow(shift_df)),
  paste0("Signed-log nonzero coefficient rows: ", nrow(nonzero_df)),
  paste0("Signed-log extreme-patient rows: ", nrow(extreme_df)),
  "",
  paste0("Output directory: ", out_dir)
)

writeLines(
  summary_lines,
  file.path(out_dir, "supplementary_export_summary.txt"),
  useBytes = TRUE
)

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(out_dir)

zip_name <- "supplementary_source_tables.zip"
files_to_zip <- c(
  "crs_nonzero_feature_mapping.csv",
  "signed_log_domain_shift_features.csv",
  "signed_log_model_nonzero_coefficients.csv",
  "signed_log_extreme_patients.csv",
  "supplementary_export_summary.txt"
)

utils::zip(zipfile = zip_name, files = files_to_zip)

cat(paste(summary_lines, collapse = "\n"))
cat("\n\nZIP file:\n", file.path(out_dir, zip_name), "\n")
