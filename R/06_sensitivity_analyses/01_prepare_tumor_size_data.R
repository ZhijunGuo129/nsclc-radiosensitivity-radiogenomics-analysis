# 01_prepare_tumor_size_data.R
#
# Prepare tumor-size and transportability analysis data.
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

############################################################
# 1. Packages
############################################################

required_packages <- c("openxlsx", "ggplot2")
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
library(ggplot2)

############################################################
# 2. Paths
############################################################


train_merged_rds <- file.path(
  project_dir,
  "07_results",
  "radiomics_features",
  "radiomics_molecular_merged.rds"
)

train_merged_csv <- file.path(
  project_dir,
  "07_results",
  "radiomics_features",
  "radiomics_molecular_merged.csv"
)

original_train_rrs_csv_1 <- file.path(
  project_dir,
  "07_results",
  "rrs_model",
  "rrs_training_cohort_analysis.csv"
)

original_train_rrs_csv_2 <- file.path(
  project_dir,
  "07_results",
  "rrs_model",
  "patient_mean_predictions.csv"
)

signed_log_train_rrs_csv <- file.path(
  lung1_root,
  "signed_log_rrs",
  "signed_log_cross_validated_predictions.csv"
)

lung1_feature_csv <- file.path(
  lung1_root,
  "radiomics_features",
  "lung1_radiomics_features.csv"
)

lung1_original_rrs_rds <- file.path(
  lung1_root,
  "external_RRS",
  "lung1_rrs_survival_dataset.rds"
)

lung1_original_rrs_csv <- file.path(
  lung1_root,
  "external_RRS",
  "lung1_rrs_survival_dataset.csv"
)

lung1_signed_log_rrs_csv <- file.path(
  lung1_root,
  "signed_log_rrs",
  "lung1_signed_log_rrs_scores.csv"
)

out_dir <- file.path(
  project_dir,
  "07_results",
  "tumor_size_domain_shift"
)

fig_dir <- file.path(
  project_dir,
  "08_figures"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
# 3. Helper functions
############################################################

read_table_auto <- function(rds_path = NULL, csv_path = NULL) {
  if (!is.null(rds_path) && file.exists(rds_path)) {
    return(readRDS(rds_path))
  }
  if (!is.null(csv_path) && file.exists(csv_path)) {
    return(read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  return(NULL)
}

standardize_id <- function(x) {
  toupper(trimws(as.character(x)))
}

znum <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  as.numeric(scale(x))
}

safe_log1p <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  log1p(pmax(x, 0))
}

detect_patient_id_col <- function(df) {
  candidates <- c("patient_id", "PatientID", "sampleid", "sample_id", "SampleID", "patient_id_molecular")
  id_col <- candidates[candidates %in% colnames(df)][1]
  if (is.na(id_col) || length(id_col) == 0) {
    stop("Cannot detect patient ID column. Columns are: ", paste(colnames(df), collapse = "; "))
  }
  id_col
}

detect_score_col <- function(df, exact_candidates, regex_candidates = NULL, exclude_regex = "group|cutoff|median_group") {
  
  for (cc in exact_candidates) {
    if (cc %in% colnames(df)) {
      if (is.numeric(suppressWarnings(as.numeric(df[[cc]])))) {
        return(cc)
      }
    }
  }
  
  if (!is.null(regex_candidates)) {
    for (pat in regex_candidates) {
      cand <- colnames(df)[grepl(pat, colnames(df), ignore.case = TRUE)]
      cand <- cand[!grepl(exclude_regex, cand, ignore.case = TRUE)]
      
      for (cc in cand) {
        x <- suppressWarnings(as.numeric(df[[cc]]))
        if (sum(is.finite(x)) >= 20 && sd(x, na.rm = TRUE) > 0) {
          return(cc)
        }
      }
    }
  }
  
  stop(
    "Cannot detect score column. Available columns are:\n",
    paste(colnames(df), collapse = "\n")
  )
}

add_size_variables <- function(df) {
  
  df2 <- df
  
  size_feature_candidates <- c(
    "original_shape_MeshVolume",
    "original_shape_VoxelVolume",
    "original_shape_Maximum3DDiameter",
    "original_shape_SurfaceArea",
    "original_shape_SurfaceVolumeRatio",
    "original_shape_MajorAxisLength",
    "original_shape_MinorAxisLength",
    "original_shape_LeastAxisLength"
  )
  
  available <- intersect(size_feature_candidates, colnames(df2))
  
  if (length(available) == 0) {
    stop("No standard PyRadiomics shape features found.")
  }
  
  for (ff in available) {
    df2[[ff]] <- suppressWarnings(as.numeric(df2[[ff]]))
  }
  
  if ("original_shape_MeshVolume" %in% colnames(df2)) {
    df2$log_MeshVolume <- safe_log1p(df2$original_shape_MeshVolume)
  }
  
  if ("original_shape_VoxelVolume" %in% colnames(df2)) {
    df2$log_VoxelVolume <- safe_log1p(df2$original_shape_VoxelVolume)
  }
  
  if ("original_shape_Maximum3DDiameter" %in% colnames(df2)) {
    df2$log_Maximum3DDiameter <- safe_log1p(df2$original_shape_Maximum3DDiameter)
  }
  
  if ("original_shape_SurfaceArea" %in% colnames(df2)) {
    df2$log_SurfaceArea <- safe_log1p(df2$original_shape_SurfaceArea)
  }
  
  df2
}

cor_one <- function(data, x_col, y_col, label = "") {
  
  x <- suppressWarnings(as.numeric(data[[x_col]]))
  y <- suppressWarnings(as.numeric(data[[y_col]]))
  
  ok <- is.finite(x) & is.finite(y)
  
  if (sum(ok) < 10) {
    return(data.frame(
      label = label,
      x = x_col,
      y = y_col,
      n = sum(ok),
      Pearson_r = NA,
      Pearson_p = NA,
      Spearman_r = NA,
      Spearman_p = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  pear <- suppressWarnings(cor.test(x[ok], y[ok], method = "pearson"))
  spear <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  
  data.frame(
    label = label,
    x = x_col,
    y = y_col,
    n = sum(ok),
    Pearson_r = as.numeric(pear$estimate),
    Pearson_p = pear$p.value,
    Spearman_r = as.numeric(spear$estimate),
    Spearman_p = spear$p.value,
    stringsAsFactors = FALSE
  )
}

run_adjusted_lm <- function(data, outcome_col, score_col, adjust_cols, model_label) {
  
  use_cols <- c(outcome_col, score_col, adjust_cols)
  use_cols <- use_cols[use_cols %in% colnames(data)]
  
  dd <- data[, use_cols, drop = FALSE]
  
  dd$outcome_y <- suppressWarnings(as.numeric(dd[[outcome_col]]))
  dd$score_z <- znum(dd[[score_col]])
  
  covar_names <- character(0)
  
  for (cc in adjust_cols) {
    if (cc %in% colnames(dd)) {
      new_name <- paste0(cc, "_z")
      dd[[new_name]] <- znum(dd[[cc]])
      covar_names <- c(covar_names, new_name)
    }
  }
  
  model_vars <- c("outcome_y", "score_z", covar_names)
  dd <- dd[complete.cases(dd[, model_vars, drop = FALSE]), , drop = FALSE]
  
  if (nrow(dd) < 20) {
    return(data.frame(
      model = model_label,
      outcome = outcome_col,
      score = score_col,
      adjust_covariates = paste(adjust_cols, collapse = "+"),
      n = nrow(dd),
      beta_score = NA,
      p_score = NA,
      model_R2 = NA,
      adj_R2 = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  formula_txt <- paste(
    "outcome_y ~ score_z",
    ifelse(length(covar_names) > 0, paste("+", paste(covar_names, collapse = "+")), "")
  )
  
  fit <- lm(as.formula(formula_txt), data = dd)
  sm <- summary(fit)
  
  beta_score <- sm$coefficients["score_z", "Estimate"]
  p_score <- sm$coefficients["score_z", "Pr(>|t|)"]
  
  data.frame(
    model = model_label,
    outcome = outcome_col,
    score = score_col,
    adjust_covariates = ifelse(length(adjust_cols) > 0, paste(adjust_cols, collapse = "+"), "none"),
    n = nrow(dd),
    beta_score = beta_score,
    p_score = p_score,
    model_R2 = sm$r.squared,
    adj_R2 = sm$adj.r.squared,
    stringsAsFactors = FALSE
  )
}

run_nested_size_compare <- function(data, outcome_col, score_col, size_covs, model_label) {
  
  size_covs <- size_covs[size_covs %in% colnames(data)]
  
  dd <- data[, c(outcome_col, score_col, size_covs), drop = FALSE]
  
  dd$outcome_y <- suppressWarnings(as.numeric(dd[[outcome_col]]))
  dd$score_z <- znum(dd[[score_col]])
  
  covar_names <- character(0)
  
  for (cc in size_covs) {
    new_name <- paste0(cc, "_z")
    dd[[new_name]] <- znum(dd[[cc]])
    covar_names <- c(covar_names, new_name)
  }
  
  dd <- dd[complete.cases(dd[, c("outcome_y", "score_z", covar_names), drop = FALSE]), , drop = FALSE]
  
  if (nrow(dd) < 20 || length(covar_names) == 0) {
    return(data.frame(
      model = model_label,
      score = score_col,
      size_covariates = paste(size_covs, collapse = "+"),
      n = nrow(dd),
      size_only_R2 = NA,
      size_plus_score_R2 = NA,
      delta_R2 = NA,
      anova_p_for_added_score = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  fit_size <- lm(
    as.formula(paste("outcome_y ~", paste(covar_names, collapse = "+"))),
    data = dd
  )
  
  fit_score <- lm(
    as.formula(paste("outcome_y ~", paste(c(covar_names, "score_z"), collapse = "+"))),
    data = dd
  )
  
  sm_size <- summary(fit_size)
  sm_score <- summary(fit_score)
  
  an <- anova(fit_size, fit_score)
  
  data.frame(
    model = model_label,
    score = score_col,
    size_covariates = paste(size_covs, collapse = "+"),
    n = nrow(dd),
    size_only_R2 = sm_size$r.squared,
    size_plus_score_R2 = sm_score$r.squared,
    delta_R2 = sm_score$r.squared - sm_size$r.squared,
    anova_p_for_added_score = an$`Pr(>F)`[2],
    stringsAsFactors = FALSE
  )
}

run_partial_correlation <- function(data, outcome_col, score_col, adjust_cols, label) {
  
  adjust_cols <- adjust_cols[adjust_cols %in% colnames(data)]
  
  dd <- data[, c(outcome_col, score_col, adjust_cols), drop = FALSE]
  dd$outcome_y <- suppressWarnings(as.numeric(dd[[outcome_col]]))
  dd$score_y <- suppressWarnings(as.numeric(dd[[score_col]]))
  
  covar_names <- character(0)
  
  for (cc in adjust_cols) {
    new_name <- paste0(cc, "_z")
    dd[[new_name]] <- znum(dd[[cc]])
    covar_names <- c(covar_names, new_name)
  }
  
  dd <- dd[complete.cases(dd[, c("outcome_y", "score_y", covar_names), drop = FALSE]), , drop = FALSE]
  
  if (nrow(dd) < 20 || length(covar_names) == 0) {
    return(data.frame(
      label = label,
      outcome = outcome_col,
      score = score_col,
      adjusted_for = paste(adjust_cols, collapse = "+"),
      n = nrow(dd),
      partial_Pearson_r = NA,
      partial_Pearson_p = NA,
      partial_Spearman_r = NA,
      partial_Spearman_p = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  covar_formula <- paste(covar_names, collapse = "+")
  
  y_resid <- resid(lm(as.formula(paste("outcome_y ~", covar_formula)), data = dd))
  s_resid <- resid(lm(as.formula(paste("score_y ~", covar_formula)), data = dd))
  
  pear <- suppressWarnings(cor.test(y_resid, s_resid, method = "pearson"))
  spear <- suppressWarnings(cor.test(y_resid, s_resid, method = "spearman", exact = FALSE))
  
  data.frame(
    label = label,
    outcome = outcome_col,
    score = score_col,
    adjusted_for = paste(adjust_cols, collapse = "+"),
    n = nrow(dd),
    partial_Pearson_r = as.numeric(pear$estimate),
    partial_Pearson_p = pear$p.value,
    partial_Spearman_r = as.numeric(spear$estimate),
    partial_Spearman_p = spear$p.value,
    stringsAsFactors = FALSE
  )
}

############################################################
# 4. Load training cohort data
############################################################

cat("\n===== Loading NSCLC-Radiogenomics training data =====\n")

train_df <- read_table_auto(
  rds_path = train_merged_rds,
  csv_path = train_merged_csv
)

if (is.null(train_df)) {
  stop("Cannot find the merged training table.")
}

train_df <- as.data.frame(train_df, check.names = FALSE)

id_col_train <- detect_patient_id_col(train_df)
train_df$patient_id <- standardize_id(train_df[[id_col_train]])

if (!"CRS_z" %in% colnames(train_df)) {
  stop("Cannot find CRS_z in training table.")
}

train_df$CRS_z <- suppressWarnings(as.numeric(train_df$CRS_z))

train_df <- add_size_variables(train_df)

cat("Training table dim:\n")
print(dim(train_df))

cat("\nAvailable size variables in training data:\n")
print(intersect(
  c("original_shape_MeshVolume", "original_shape_VoxelVolume",
    "original_shape_Maximum3DDiameter", "original_shape_SurfaceArea",
    "log_MeshVolume", "log_VoxelVolume",
    "log_Maximum3DDiameter", "log_SurfaceArea"),
  colnames(train_df)
))

############################################################
# 5. Load original training RRS
############################################################

cat("\n===== Loading original training RRS =====\n")

original_rrs_df <- NULL

if (file.exists(original_train_rrs_csv_1)) {
  original_rrs_df <- read.csv(
    original_train_rrs_csv_1,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
} else if (file.exists(original_train_rrs_csv_2)) {
  original_rrs_df <- read.csv(
    original_train_rrs_csv_2,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

has_original_rrs <- FALSE

if (!is.null(original_rrs_df)) {
  
  id_col <- detect_patient_id_col(original_rrs_df)
  original_rrs_df$patient_id <- standardize_id(original_rrs_df[[id_col]])
  
  cat("Original RRS table columns:\n")
  print(colnames(original_rrs_df))
  
  original_rrs_col <- detect_score_col(
    original_rrs_df,
    exact_candidates = c(
      "RRS_CV",
      "RRS_CV_mean",
      "RRS_CV_mean_prediction",
      "pred_mean",
      "Predicted_CRS_z",
      "predicted_CRS_z",
      "mean_prediction",
      "CV_prediction_mean"
    ),
    regex_candidates = c("RRS", "pred", "prediction")
  )
  
  cat("Detected original training RRS column:\n")
  print(original_rrs_col)
  
  original_rrs_use <- original_rrs_df[, c("patient_id", original_rrs_col), drop = FALSE]
  colnames(original_rrs_use)[2] <- "Original_RRS_CV"
  original_rrs_use$Original_RRS_CV <- suppressWarnings(as.numeric(original_rrs_use$Original_RRS_CV))
  
  train_df <- merge(
    train_df,
    original_rrs_use,
    by = "patient_id",
    all.x = TRUE
  )
  
  has_original_rrs <- TRUE
  
} else {
  cat("Original RRS table not found. Original RRS analyses will be skipped.\n")
}

############################################################
# 6. Load signed-log training RRS
############################################################

cat("\n===== Loading signed-log training RRS =====\n")

has_signed_log_rrs <- FALSE

if (file.exists(signed_log_train_rrs_csv)) {
  
  signed_log_rrs_df <- read.csv(
    signed_log_train_rrs_csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  id_col <- detect_patient_id_col(signed_log_rrs_df)
  signed_log_rrs_df$patient_id <- standardize_id(signed_log_rrs_df[[id_col]])
  
  signed_log_col <- detect_score_col(
    signed_log_rrs_df,
    exact_candidates = c("SignedLog_RRS_CV_mean", "SignedLog_RRS", "RRS_CV_mean"),
    regex_candidates = c("SignedLog", "RRS")
  )
  
  cat("Detected signed-log training RRS column:\n")
  print(signed_log_col)
  
  signed_log_rrs_use <- signed_log_rrs_df[, c("patient_id", signed_log_col), drop = FALSE]
  colnames(signed_log_rrs_use)[2] <- "SignedLog_RRS_CV"
  signed_log_rrs_use$SignedLog_RRS_CV <- suppressWarnings(as.numeric(signed_log_rrs_use$SignedLog_RRS_CV))
  
  train_df <- merge(
    train_df,
    signed_log_rrs_use,
    by = "patient_id",
    all.x = TRUE
  )
  
  has_signed_log_rrs <- TRUE
  
} else {
  cat("Signed-log training RRS table not found. Signed-log RRS analyses will be skipped.\n")
}

score_cols_train <- c()
if (has_original_rrs) score_cols_train <- c(score_cols_train, "Original_RRS_CV")
if (has_signed_log_rrs) score_cols_train <- c(score_cols_train, "SignedLog_RRS_CV")

if (length(score_cols_train) == 0) {
  stop("No training RRS score columns detected.")
}

cat("\nTraining analysis score columns:\n")
print(score_cols_train)

############################################################
# 7. Training cohort correlations
############################################################

primary_size_vars <- intersect(
  c(
    "log_MeshVolume",
    "log_VoxelVolume",
    "log_Maximum3DDiameter",
    "log_SurfaceArea",
    "original_shape_SurfaceVolumeRatio"
  ),
  colnames(train_df)
)

cat("\nPrimary size variables:\n")
print(primary_size_vars)

train_cor_list <- list()

for (sv in primary_size_vars) {
  train_cor_list[[length(train_cor_list) + 1]] <- cor_one(
    train_df,
    x_col = "CRS_z",
    y_col = sv,
    label = "Training_CRSz_vs_size"
  )
}

for (sc in score_cols_train) {
  train_cor_list[[length(train_cor_list) + 1]] <- cor_one(
    train_df,
    x_col = "CRS_z",
    y_col = sc,
    label = paste0("Training_CRSz_vs_", sc)
  )
  
  for (sv in primary_size_vars) {
    train_cor_list[[length(train_cor_list) + 1]] <- cor_one(
      train_df,
      x_col = sc,
      y_col = sv,
      label = paste0("Training_", sc, "_vs_size")
    )
  }
}

train_cor_df <- do.call(rbind, train_cor_list)
train_cor_df$Pearson_FDR <- p.adjust(train_cor_df$Pearson_p, method = "BH")
train_cor_df$Spearman_FDR <- p.adjust(train_cor_df$Spearman_p, method = "BH")

cat("\n===== Training correlation summary =====\n")
print(train_cor_df)

############################################################
# 8. Training adjusted models
############################################################

adjustment_sets <- list(
  none = character(0),
  log_MeshVolume = intersect(c("log_MeshVolume"), colnames(train_df)),
  log_Maximum3DDiameter = intersect(c("log_Maximum3DDiameter"), colnames(train_df)),
  log_MeshVolume_plus_log_Maximum3DDiameter = intersect(
    c("log_MeshVolume", "log_Maximum3DDiameter"),
    colnames(train_df)
  )
)

lm_list <- list()

for (sc in score_cols_train) {
  for (adj_name in names(adjustment_sets)) {
    lm_list[[length(lm_list) + 1]] <- run_adjusted_lm(
      data = train_df,
      outcome_col = "CRS_z",
      score_col = sc,
      adjust_cols = adjustment_sets[[adj_name]],
      model_label = paste0(sc, "_adjusted_", adj_name)
    )
  }
}

train_lm_df <- do.call(rbind, lm_list)

cat("\n===== Training adjusted linear models: CRS_z ~ RRS + tumor size =====\n")
print(train_lm_df)

############################################################
# 9. Nested size-only vs size+RRS comparison
############################################################

nested_list <- list()

primary_adjust <- intersect(
  c("log_MeshVolume", "log_Maximum3DDiameter"),
  colnames(train_df)
)

for (sc in score_cols_train) {
  nested_list[[length(nested_list) + 1]] <- run_nested_size_compare(
    data = train_df,
    outcome_col = "CRS_z",
    score_col = sc,
    size_covs = primary_adjust,
    model_label = paste0(sc, "_adds_beyond_size")
  )
}

nested_df <- do.call(rbind, nested_list)

cat("\n===== Nested model comparison: size-only vs size+RRS =====\n")
print(nested_df)

############################################################
# 10. Partial correlation adjusted for tumor size
############################################################

partial_list <- list()

for (sc in score_cols_train) {
  partial_list[[length(partial_list) + 1]] <- run_partial_correlation(
    data = train_df,
    outcome_col = "CRS_z",
    score_col = sc,
    adjust_cols = primary_adjust,
    label = paste0(sc, "_partial_adjusted_for_size")
  )
}

partial_df <- do.call(rbind, partial_list)

cat("\n===== Partial correlation adjusted for tumor size =====\n")
print(partial_df)

############################################################
# 11. Load Lung1 external data
############################################################

cat("\n===== Loading Lung1 external RRS and features =====\n")

lung_feat <- read.csv(
  lung1_feature_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

id_col_lung <- detect_patient_id_col(lung_feat)
lung_feat$patient_id <- standardize_id(lung_feat[[id_col_lung]])
lung_feat <- add_size_variables(lung_feat)

lung_df <- lung_feat

has_lung_original <- FALSE
has_lung_signed_log <- FALSE

lung_orig <- read_table_auto(
  rds_path = lung1_original_rrs_rds,
  csv_path = lung1_original_rrs_csv
)

if (!is.null(lung_orig)) {
  id_col <- detect_patient_id_col(lung_orig)
  lung_orig$patient_id <- standardize_id(lung_orig[[id_col]])
  
  orig_col <- detect_score_col(
    lung_orig,
    exact_candidates = c("Lung1_RRS", "RRS", "Original_Lung1_RRS"),
    regex_candidates = c("Lung1_RRS", "RRS")
  )
  
  lung_orig_use <- lung_orig[, c("patient_id", orig_col), drop = FALSE]
  colnames(lung_orig_use)[2] <- "Original_Lung1_RRS"
  lung_orig_use$Original_Lung1_RRS <- suppressWarnings(as.numeric(lung_orig_use$Original_Lung1_RRS))
  
  lung_df <- merge(
    lung_df,
    lung_orig_use,
    by = "patient_id",
    all.x = TRUE
  )
  
  has_lung_original <- TRUE
  
  cat("Detected Lung1 original RRS column:\n")
  print(orig_col)
}

if (file.exists(lung1_signed_log_rrs_csv)) {
  
  lung_rob <- read.csv(
    lung1_signed_log_rrs_csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  id_col <- detect_patient_id_col(lung_rob)
  lung_rob$patient_id <- standardize_id(lung_rob[[id_col]])
  
  rob_col <- detect_score_col(
    lung_rob,
    exact_candidates = c("SignedLog_Lung1_RRS", "SignedLog_RRS", "RRS"),
    regex_candidates = c("SignedLog", "RRS")
  )
  
  lung_signed_log_use <- lung_rob[, c("patient_id", rob_col), drop = FALSE]
  colnames(lung_signed_log_use)[2] <- "SignedLog_Lung1_RRS"
  lung_signed_log_use$SignedLog_Lung1_RRS <- suppressWarnings(as.numeric(lung_signed_log_use$SignedLog_Lung1_RRS))
  
  lung_df <- merge(
    lung_df,
    lung_signed_log_use,
    by = "patient_id",
    all.x = TRUE
  )
  
  has_lung_signed_log <- TRUE
  
  cat("Detected Lung1 signed-log RRS column:\n")
  print(rob_col)
}

lung_score_cols <- c()
if (has_lung_original) lung_score_cols <- c(lung_score_cols, "Original_Lung1_RRS")
if (has_lung_signed_log) lung_score_cols <- c(lung_score_cols, "SignedLog_Lung1_RRS")

cat("\nLung1 analysis score columns:\n")
print(lung_score_cols)

############################################################
# 12. Lung1 RRS vs tumor size correlations
############################################################

lung_size_vars <- intersect(
  c(
    "log_MeshVolume",
    "log_VoxelVolume",
    "log_Maximum3DDiameter",
    "log_SurfaceArea",
    "original_shape_SurfaceVolumeRatio"
  ),
  colnames(lung_df)
)

lung_cor_list <- list()

if (length(lung_score_cols) > 0) {
  for (sc in lung_score_cols) {
    for (sv in lung_size_vars) {
      lung_cor_list[[length(lung_cor_list) + 1]] <- cor_one(
        lung_df,
        x_col = sc,
        y_col = sv,
        label = paste0("Lung1_", sc, "_vs_size")
      )
    }
  }
  
  lung_cor_df <- do.call(rbind, lung_cor_list)
  lung_cor_df$Pearson_FDR <- p.adjust(lung_cor_df$Pearson_p, method = "BH")
  lung_cor_df$Spearman_FDR <- p.adjust(lung_cor_df$Spearman_p, method = "BH")
  
} else {
  lung_cor_df <- data.frame()
}

cat("\n===== Lung1 RRS vs tumor size correlation summary =====\n")
print(lung_cor_df)

############################################################
# 13. Analysis summary
############################################################

interpret_rows <- list()

for (sc in score_cols_train) {
  
  unadj <- train_lm_df[
    train_lm_df$score == sc &
      train_lm_df$adjust_covariates == "none",
  ]
  
  adj <- train_lm_df[
    train_lm_df$score == sc &
      train_lm_df$adjust_covariates == paste(primary_adjust, collapse = "+"),
  ]
  
  nested <- nested_df[nested_df$score == sc, ]
  
  part <- partial_df[partial_df$score == sc, ]
  
  conclusion <- "Unable to determine."
  
  if (nrow(adj) > 0 && is.finite(adj$p_score[1])) {
    if (adj$p_score[1] < 0.05) {
      conclusion <- "RRS remained significantly associated with CRS_z after adjustment for tumor size, supporting that RRS is not merely a tumor-size surrogate."
    } else {
      conclusion <- "RRS association with CRS_z was weakened or non-significant after tumor-size adjustment, suggesting partial dependence on tumor size or limited independent signal."
    }
  }
  
  interpret_rows[[length(interpret_rows) + 1]] <- data.frame(
    score = sc,
    unadjusted_beta = ifelse(nrow(unadj) > 0, unadj$beta_score[1], NA),
    unadjusted_p = ifelse(nrow(unadj) > 0, unadj$p_score[1], NA),
    adjusted_size_beta = ifelse(nrow(adj) > 0, adj$beta_score[1], NA),
    adjusted_size_p = ifelse(nrow(adj) > 0, adj$p_score[1], NA),
    partial_Pearson_r = ifelse(nrow(part) > 0, part$partial_Pearson_r[1], NA),
    partial_Pearson_p = ifelse(nrow(part) > 0, part$partial_Pearson_p[1], NA),
    delta_R2_beyond_size = ifelse(nrow(nested) > 0, nested$delta_R2[1], NA),
    added_score_p = ifelse(nrow(nested) > 0, nested$anova_p_for_added_score[1], NA),
    interpretation = conclusion,
    stringsAsFactors = FALSE
  )
}

interpretation_df <- do.call(rbind, interpret_rows)

cat("\n===== Tumor-size and domain-shift analysis summary =====\n")
print(interpretation_df)

############################################################
# 14. Save outputs
############################################################

write.csv(
  train_df,
  file = file.path(out_dir, "training_size_analysis_dataset.csv"),
  row.names = FALSE
)

write.csv(
  train_cor_df,
  file = file.path(out_dir, "training_correlations_CRS_RRS_size.csv"),
  row.names = FALSE
)

write.csv(
  train_lm_df,
  file = file.path(out_dir, "training_adjusted_linear_models.csv"),
  row.names = FALSE
)

write.csv(
  nested_df,
  file = file.path(out_dir, "training_nested_size_vs_size_plus_RRS.csv"),
  row.names = FALSE
)

write.csv(
  partial_df,
  file = file.path(out_dir, "training_partial_correlations_adjusted_size.csv"),
  row.names = FALSE
)

write.csv(
  lung_df,
  file = file.path(out_dir, "lung1_size_analysis_dataset.csv"),
  row.names = FALSE
)

write.csv(
  lung_cor_df,
  file = file.path(out_dir, "lung1_correlations_RRS_size.csv"),
  row.names = FALSE
)

write.csv(
  interpretation_df,
  file = file.path(out_dir, "tumor_size_analysis_summary.csv"),
  row.names = FALSE
)

############################################################
# 15. Excel workbook
############################################################

wb <- createWorkbook()

addWorksheet(wb, "analysis_summary")
writeData(wb, "analysis_summary", interpretation_df)

addWorksheet(wb, "training_correlations")
writeData(wb, "training_correlations", train_cor_df)

addWorksheet(wb, "training_adjusted_LM")
writeData(wb, "training_adjusted_LM", train_lm_df)

addWorksheet(wb, "nested_size_plus_RRS")
writeData(wb, "nested_size_plus_RRS", nested_df)

addWorksheet(wb, "partial_correlations")
writeData(wb, "partial_correlations", partial_df)

addWorksheet(wb, "Lung1_RRS_size_cor")
writeData(wb, "Lung1_RRS_size_cor", lung_cor_df)

addWorksheet(wb, "training_dataset")
writeData(wb, "training_dataset", train_df)

addWorksheet(wb, "Lung1_dataset")
writeData(wb, "Lung1_dataset", lung_df)

for (sheet in names(wb)) {
  setColWidths(wb, sheet, cols = 1:30, widths = "auto")
}

saveWorkbook(
  wb,
  file = file.path(out_dir, "tumor_size_domain_shift_summary.xlsx"),
  overwrite = TRUE
)

############################################################
# 16. Figures
############################################################

# Training: CRS_z vs RRS
for (sc in score_cols_train) {
  
  p <- ggplot(
    train_df,
    aes(x = .data[[sc]], y = CRS_z)
  ) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE) +
    theme_bw(base_size = 12) +
    labs(
      title = paste0("Training cohort: CRS_z vs ", sc),
      x = sc,
      y = "CRS_z"
    )
  
  ggsave(
    filename = file.path(fig_dir, paste0("training_CRSz_vs_", sc, ".png")),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
}

# Training: RRS vs log tumor volume
if ("log_MeshVolume" %in% colnames(train_df)) {
  
  for (sc in score_cols_train) {
    
    p <- ggplot(
      train_df,
      aes(x = log_MeshVolume, y = .data[[sc]])
    ) +
      geom_point(size = 2, alpha = 0.8) +
      geom_smooth(method = "lm", se = TRUE) +
      theme_bw(base_size = 12) +
      labs(
        title = paste0("Training cohort: ", sc, " vs log tumor volume"),
        x = "log1p(MeshVolume)",
        y = sc
      )
    
    ggsave(
      filename = file.path(fig_dir, paste0("training_", sc, "_vs_log_MeshVolume.png")),
      plot = p,
      width = 6,
      height = 5,
      dpi = 300
    )
  }
}

# Training residual plot adjusted for size
if (length(primary_adjust) > 0) {
  
  for (sc in score_cols_train) {
    
    dd <- train_df[, c("CRS_z", sc, primary_adjust), drop = FALSE]
    dd <- dd[complete.cases(dd), , drop = FALSE]
    
    if (nrow(dd) >= 20) {
      
      for (cc in primary_adjust) {
        dd[[paste0(cc, "_z")]] <- znum(dd[[cc]])
      }
      
      covar_formula <- paste(paste0(primary_adjust, "_z"), collapse = "+")
      
      dd$CRS_z_resid_size <- resid(
        lm(as.formula(paste("CRS_z ~", covar_formula)), data = dd)
      )
      
      dd$RRS_resid_size <- resid(
        lm(as.formula(paste(sc, "~", covar_formula)), data = dd)
      )
      
      p <- ggplot(
        dd,
        aes(x = RRS_resid_size, y = CRS_z_resid_size)
      ) +
        geom_point(size = 2, alpha = 0.8) +
        geom_smooth(method = "lm", se = TRUE) +
        theme_bw(base_size = 12) +
        labs(
          title = paste0("Training cohort: size-adjusted CRS_z vs ", sc),
          x = paste0(sc, " residual adjusted for tumor size"),
          y = "CRS_z residual adjusted for tumor size"
        )
      
      ggsave(
        filename = file.path(fig_dir, paste0("training_size_adjusted_CRSz_vs_", sc, ".png")),
        plot = p,
        width = 6,
        height = 5,
        dpi = 300
      )
    }
  }
}

# Lung1: external RRS vs log tumor volume
if ("log_MeshVolume" %in% colnames(lung_df) && length(lung_score_cols) > 0) {
  
  for (sc in lung_score_cols) {
    
    p <- ggplot(
      lung_df,
      aes(x = log_MeshVolume, y = .data[[sc]])
    ) +
      geom_point(size = 1.8, alpha = 0.75) +
      geom_smooth(method = "lm", se = TRUE) +
      theme_bw(base_size = 12) +
      labs(
        title = paste0("Lung1: ", sc, " vs log tumor volume"),
        x = "log1p(MeshVolume)",
        y = sc
      )
    
    ggsave(
      filename = file.path(fig_dir, paste0("lung1_", sc, "_vs_log_MeshVolume.png")),
      plot = p,
      width = 6,
      height = 5,
      dpi = 300
    )
  }
}

############################################################
# 17. Manuscript text draft
############################################################

txt_lines <- c(
  "Tumor-size and domain-shift analysis summary",
  "",
  "Purpose:",
  "This analysis evaluates whether RRS is merely a surrogate of tumor volume or maximum tumor diameter.",
  "",
  "Primary statistical tests:",
  "1. Correlation between CRS_z, RRS and tumor size features.",
  "2. Linear regression models of CRS_z using RRS with and without adjustment for log-transformed tumor volume and maximum 3D diameter.",
  "3. Nested model comparison between size-only and size-plus-RRS models.",
  "4. Partial correlation between CRS_z and RRS adjusted for tumor size.",
  "",
  "Analysis summary:",
  paste(capture.output(print(interpretation_df)), collapse = "\n"),
  "",
  "Recommended wording depends on adjusted results:",
  "If RRS remains significant after size adjustment:",
  "RRS remained associated with CRS_z after adjustment for tumor volume and maximum diameter, suggesting that the radiomic surrogate captured CRS-related information beyond gross tumor size.",
  "",
  "If RRS becomes non-significant after size adjustment:",
  "The association between RRS and CRS_z was attenuated after adjustment for tumor volume and maximum diameter, suggesting that tumor size or tumor-burden-related spatial heterogeneity partially contributed to the radiomic CRS signal."
)

writeLines(
  txt_lines,
  con = file.path(out_dir, "tumor_size_key_results.txt"),
  useBytes = TRUE
)

############################################################
# 18. Done
############################################################

cat("\n===== Tumor-size and domain-shift analysis completed =====\n")
cat("Main output files:\n")
cat(file.path(out_dir, "tumor_size_domain_shift_summary.xlsx"), "\n")
cat(file.path(out_dir, "tumor_size_analysis_summary.csv"), "\n")
cat(file.path(out_dir, "tumor_size_key_results.txt"), "\n")
cat(file.path(fig_dir, "training_size_adjusted_CRSz_vs_Original_RRS_CV.png"), "\n")

cat("1. Training correlation summary\n")
cat("2. Training adjusted linear models: CRS_z ~ RRS + tumor size\n")
cat("3. Nested model comparison: size-only vs size+RRS\n")
cat("4. Partial correlation adjusted for tumor size\n")
cat("5. Lung1 RRS vs tumor size correlation summary\n")
cat("6. Tumor-size and domain-shift analysis summary\n")
