# 05_merge_radiomics_and_molecular_data.R
#
# Apply radiomic quality control and merge radiomic features with molecular scores.
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

feature_csv <- "07_results/radiomics_features/nsclc_radiomics_features.csv"
feature_log_csv <- "07_results/radiomics_features/nsclc_radiomics_features_log.csv"

molecular_score_rds <- "05_molecular_scores/NSCLC_Radiogenomics_patient_molecular_scores.rds"

out_dir <- "07_results/radiomics_features"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

metadata_out_dir <- "02_metadata"
dir.create(metadata_out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(feature_csv)) {
  stop("Cannot find PyRadiomics feature CSV: ", feature_csv)
}

if (!file.exists(molecular_score_rds)) {
  stop("Cannot find molecular score RDS: ", molecular_score_rds)
}

############################################################
# 3. Load data
############################################################

cat("\n===== Loading PyRadiomics feature table =====\n")

feature_df <- read.csv(
  feature_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Feature table dim:\n")
print(dim(feature_df))

cat("\nFirst 20 columns:\n")
print(colnames(feature_df)[1:min(20, ncol(feature_df))])

cat("\n===== Loading PyRadiomics log table =====\n")

feature_log <- read.csv(
  feature_log_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Feature log status table:\n")
print(table(feature_log$status, useNA = "ifany"))

cat("\n===== Loading molecular score table =====\n")

mol_scores <- readRDS(molecular_score_rds)

cat("Molecular score table dim:\n")
print(dim(mol_scores))

cat("\nMolecular score columns:\n")
print(colnames(mol_scores))

############################################################
# 4. Basic ID audit
############################################################

if (!"patient_id" %in% colnames(feature_df)) {
  stop("patient_id column not found in PyRadiomics feature table.")
}

############################################################
# Standardize the molecular patient identifier column
############################################################

cat("\n===== Molecular ID column auto-detection =====\n")
cat("Molecular score columns:\n")
print(colnames(mol_scores))

# 1. Use an existing patient_id column when available.
if ("patient_id" %in% colnames(mol_scores)) {
  cat("Detected patient_id column directly.\n")
  mol_scores$patient_id <- toupper(trimws(as.character(mol_scores$patient_id)))
  
} else {
  
  # 2. Check common identifier column names.
  candidate_id_cols <- c(
    "patient_id_molecular",
    "sample_id",
    "SampleID",
    "sample",
    "Sample",
    "ID",
    "id",
    "case_id",
    "CaseID"
  )
  
  hit_cols <- intersect(candidate_id_cols, colnames(mol_scores))
  
  if (length(hit_cols) > 0) {
    id_col <- hit_cols[1]
    cat("Detected molecular patient ID column:\n")
    print(id_col)
    mol_scores$patient_id <- toupper(trimws(as.character(mol_scores[[id_col]])))
    
  } else {
    
    # 3. Search column values for R01-style identifiers.
    r01_col_hits <- sapply(colnames(mol_scores), function(cc) {
      any(grepl("^R01-[0-9]{3}$", toupper(trimws(as.character(mol_scores[[cc]])))))
    })
    
    if (any(r01_col_hits)) {
      id_col <- names(r01_col_hits)[which(r01_col_hits)[1]]
      cat("Detected R01-style patient ID column by content:\n")
      print(id_col)
      mol_scores$patient_id <- toupper(trimws(as.character(mol_scores[[id_col]])))
      
    } else {
      
      # 4. Use row names as a final R01-style identifier source.
      rn <- rownames(mol_scores)
      if (!is.null(rn) && any(grepl("^R01-[0-9]{3}$", toupper(trimws(rn))))) {
        cat("Detected patient ID from rownames.\n")
        mol_scores$patient_id <- toupper(trimws(rn))
      } else {
        stop(
          "Cannot automatically detect patient ID column in molecular score table. ",
        )
      }
    }
  }
}

cat("\nMolecular patient_id preview:\n")
print(head(mol_scores$patient_id, 20))

cat("\nNumber of unique molecular patient IDs:\n")
print(length(unique(mol_scores$patient_id)))

feature_df$patient_id <- toupper(trimws(feature_df$patient_id))
mol_scores$patient_id <- toupper(trimws(mol_scores$patient_id))

id_summary <- data.frame(
  item = c(
    "radiomics_patients",
    "molecular_patients",
    "overlap_patients",
    "radiomics_not_in_molecular",
    "molecular_not_in_radiomics"
  ),
  value = c(
    length(unique(feature_df$patient_id)),
    length(unique(mol_scores$patient_id)),
    length(intersect(unique(feature_df$patient_id), unique(mol_scores$patient_id))),
    length(setdiff(unique(feature_df$patient_id), unique(mol_scores$patient_id))),
    length(setdiff(unique(mol_scores$patient_id), unique(feature_df$patient_id)))
  ),
  stringsAsFactors = FALSE
)

cat("\n===== ID matching summary =====\n")
print(id_summary)

cat("\nRadiomics patients not in molecular table:\n")
print(setdiff(unique(feature_df$patient_id), unique(mol_scores$patient_id)))

cat("\nMolecular patients not in radiomics table:\n")
print(setdiff(unique(mol_scores$patient_id), unique(feature_df$patient_id)))

############################################################
# 5. Separate diagnostics and original radiomics features
############################################################

diagnostics_cols <- grep("^diagnostics_", colnames(feature_df), value = TRUE)
radiomics_cols <- grep("^original_", colnames(feature_df), value = TRUE)

cat("\n===== Feature column audit =====\n")
cat("Diagnostics columns:\n")
print(length(diagnostics_cols))

cat("Original radiomics feature columns:\n")
print(length(radiomics_cols))

if (length(radiomics_cols) == 0) {
  stop("No original_ radiomics feature columns detected.")
}

radiomics_raw <- feature_df[, c("patient_id", radiomics_cols), drop = FALSE]

############################################################
# 6. Convert radiomics features to numeric
############################################################

radiomics_num <- radiomics_raw

for (cc in radiomics_cols) {
  radiomics_num[[cc]] <- suppressWarnings(as.numeric(radiomics_num[[cc]]))
}

############################################################
# 7. Feature QC
############################################################

feature_qc <- data.frame(
  feature = radiomics_cols,
  missing_n = NA_integer_,
  missing_rate = NA_real_,
  infinite_n = NA_integer_,
  finite_n = NA_integer_,
  mean = NA_real_,
  sd = NA_real_,
  min = NA_real_,
  q25 = NA_real_,
  median = NA_real_,
  q75 = NA_real_,
  max = NA_real_,
  zero_variance = NA,
  stringsAsFactors = FALSE
)

for (i in seq_along(radiomics_cols)) {
  cc <- radiomics_cols[i]
  x <- radiomics_num[[cc]]
  
  miss <- is.na(x)
  infv <- is.infinite(x)
  finite_x <- x[is.finite(x)]
  
  feature_qc$missing_n[i] <- sum(miss)
  feature_qc$missing_rate[i] <- mean(miss)
  feature_qc$infinite_n[i] <- sum(infv)
  feature_qc$finite_n[i] <- length(finite_x)
  
  if (length(finite_x) > 0) {
    feature_qc$mean[i] <- mean(finite_x)
    feature_qc$sd[i] <- sd(finite_x)
    feature_qc$min[i] <- min(finite_x)
    feature_qc$q25[i] <- as.numeric(quantile(finite_x, 0.25, names = FALSE))
    feature_qc$median[i] <- median(finite_x)
    feature_qc$q75[i] <- as.numeric(quantile(finite_x, 0.75, names = FALSE))
    feature_qc$max[i] <- max(finite_x)
    feature_qc$zero_variance[i] <- is.na(sd(finite_x)) || sd(finite_x) == 0
  } else {
    feature_qc$zero_variance[i] <- TRUE
  }
}

cat("\n===== Feature QC summary =====\n")
cat("Total original radiomics features:\n")
print(length(radiomics_cols))

cat("Features with any missing values:\n")
print(sum(feature_qc$missing_n > 0))

cat("Features with any infinite values:\n")
print(sum(feature_qc$infinite_n > 0))

cat("Zero-variance features:\n")
print(sum(feature_qc$zero_variance))

############################################################
# 8. Clean features
############################################################

drop_features <- feature_qc$feature[
  feature_qc$missing_rate > 0.20 |
    feature_qc$infinite_n > 0 |
    feature_qc$zero_variance
]

keep_features <- setdiff(radiomics_cols, drop_features)

cat("\n===== Clean feature selection =====\n")
cat("Dropped features:\n")
print(length(drop_features))
print(drop_features)

cat("Kept features:\n")
print(length(keep_features))

radiomics_clean <- radiomics_num[, c("patient_id", keep_features), drop = FALSE]

############################################################
# 9. Merge with molecular scores
############################################################

merged_df <- merge(
  mol_scores,
  radiomics_clean,
  by = "patient_id",
  all = FALSE
)

cat("\n===== Merged molecular + radiomics table =====\n")
print(dim(merged_df))

cat("\nMerged table columns preview:\n")
print(colnames(merged_df)[1:min(30, ncol(merged_df))])

############################################################
# 10. Identify likely target columns
############################################################

target_candidates <- grep(
  "CRS|RSI|Hypoxia|hypoxia|HYP",
  colnames(merged_df),
  value = TRUE
)

cat("\n===== Target candidate columns =====\n")
print(target_candidates)

############################################################
# 11. Save outputs
############################################################

write.csv(
  radiomics_raw,
  file = file.path(out_dir, "radiomics_features_raw.csv"),
  row.names = FALSE
)

write.csv(
  radiomics_clean,
  file = file.path(out_dir, "radiomics_features_clean.csv"),
  row.names = FALSE
)

write.csv(
  feature_qc,
  file = file.path(out_dir, "radiomics_feature_qc.csv"),
  row.names = FALSE
)

write.csv(
  merged_df,
  file = file.path(out_dir, "radiomics_molecular_merged.csv"),
  row.names = FALSE
)

saveRDS(
  radiomics_raw,
  file = file.path(out_dir, "radiomics_features_raw.rds")
)

saveRDS(
  radiomics_clean,
  file = file.path(out_dir, "radiomics_features_clean.rds")
)

saveRDS(
  feature_qc,
  file = file.path(out_dir, "radiomics_feature_qc.rds")
)

saveRDS(
  merged_df,
  file = file.path(out_dir, "radiomics_molecular_merged.rds")
)

qc_summary <- data.frame(
  item = c(
    "radiomics_patients",
    "molecular_patients",
    "merged_patients",
    "diagnostics_columns",
    "original_radiomics_features_raw",
    "features_with_missing",
    "features_with_infinite",
    "zero_variance_features",
    "dropped_features",
    "kept_features",
    "target_candidate_columns"
  ),
  value = c(
    length(unique(feature_df$patient_id)),
    length(unique(mol_scores$patient_id)),
    nrow(merged_df),
    length(diagnostics_cols),
    length(radiomics_cols),
    sum(feature_qc$missing_n > 0),
    sum(feature_qc$infinite_n > 0),
    sum(feature_qc$zero_variance),
    length(drop_features),
    length(keep_features),
    paste(target_candidates, collapse = "; ")
  ),
  stringsAsFactors = FALSE
)

write.csv(
  qc_summary,
  file = file.path(out_dir, "radiomics_molecular_merge_summary.csv"),
  row.names = FALSE
)

wb <- createWorkbook()

addWorksheet(wb, "QC_summary")
writeData(wb, "QC_summary", qc_summary)

addWorksheet(wb, "ID_summary")
writeData(wb, "ID_summary", id_summary)

addWorksheet(wb, "feature_QC")
writeData(wb, "feature_QC", feature_qc)

addWorksheet(wb, "dropped_features")
writeData(wb, "dropped_features", data.frame(feature = drop_features))

addWorksheet(wb, "radiomics_clean_preview")
writeData(wb, "radiomics_clean_preview", head(radiomics_clean, 100))

addWorksheet(wb, "merged_preview")
writeData(wb, "merged_preview", head(merged_df, 100))

saveWorkbook(
  wb,
  file = file.path(out_dir, "radiomics_molecular_merge_summary.xlsx"),
  overwrite = TRUE
)

cat("\n===== Radiomics quality control and molecular-data merge completed =====\n")
cat("Main output files:\n")
cat("1. 07_results/radiomics_features/radiomics_molecular_merged.csv\n")
cat("2. 07_results/radiomics_features/radiomics_molecular_merged.rds\n")
cat("3. 07_results/radiomics_features/radiomics_molecular_merge_summary.xlsx\n")

cat("1. ID matching summary\n")
cat("2. Feature QC summary\n")
cat("3. Clean feature selection\n")
cat("4. Merged molecular + radiomics table\n")
cat("5. Target candidate columns\n")
