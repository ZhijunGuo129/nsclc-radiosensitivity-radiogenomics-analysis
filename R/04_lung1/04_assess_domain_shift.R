# 04_assess_domain_shift.R
#
# Quantify radiomic domain shift between training data and Lung1.
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

required_packages <- c("glmnet", "openxlsx", "survival")
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

library(glmnet)
library(openxlsx)
library(survival)

############################################################
# 2. Paths
############################################################


train_rds <- file.path(
  project_dir,
  "07_results",
  "radiomics_features",
  "radiomics_molecular_merged.rds"
)

train_csv <- file.path(
  project_dir,
  "07_results",
  "radiomics_features",
  "radiomics_molecular_merged.csv"
)

lung1_feature_csv <- file.path(
  lung1_root,
  "radiomics_features",
  "lung1_radiomics_features.csv"
)

analysis_rds <- file.path(
  lung1_root,
  "external_RRS",
  "lung1_rrs_survival_dataset.rds"
)

frozen_model_rds <- file.path(
  project_dir,
  "07_results",
  "rrs_model",
  "rrs_elastic_net_model.rds"
)

out_dir <- file.path(
  lung1_root,
  "domain_shift_diagnostics"
)

fig_dir <- file.path(
  lung1_root,
  "figures"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
# 3. Load data
############################################################

cat("\n===== Loading data =====\n")

if (file.exists(train_rds)) {
  train_df <- readRDS(train_rds)
} else if (file.exists(train_csv)) {
  train_df <- read.csv(train_csv, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  stop("Cannot find training radiomics table.")
}

lung_feat <- read.csv(
  lung1_feature_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

analysis_df <- readRDS(analysis_rds)
frozen_model <- readRDS(frozen_model_rds)

cat("Training table dim:\n")
print(dim(train_df))

cat("\nLung1 feature table dim:\n")
print(dim(lung_feat))

cat("\nLung1 analysis table dim:\n")
print(dim(analysis_df))

cat("\nFrozen model names:\n")
print(names(frozen_model))

############################################################
# 4. Standardize IDs and features
############################################################

train_df$patient_id <- toupper(trimws(as.character(train_df$patient_id)))
lung_feat$patient_id <- toupper(trimws(as.character(lung_feat$patient_id)))
analysis_df$patient_id <- toupper(trimws(as.character(analysis_df$patient_id)))

model_features <- frozen_model$final_features

cat("\n===== Model feature audit =====\n")
cat("Frozen model required features:\n")
print(length(model_features))

missing_train <- setdiff(model_features, colnames(train_df))
missing_lung1 <- setdiff(model_features, colnames(lung_feat))

cat("\nMissing in training table:\n")
print(missing_train)

cat("\nMissing in Lung1 table:\n")
print(missing_lung1)

if (length(missing_train) > 0) {
  stop("Model features missing in training table.")
}

if (length(missing_lung1) > 0) {
  stop("Model features missing in Lung1 table.")
}

############################################################
# 5. Extract coefficients
############################################################

coef_mat <- as.matrix(
  coef(
    frozen_model$final_cvfit,
    s = frozen_model$final_lambda
  )
)

coef_df <- data.frame(
  feature = rownames(coef_mat),
  coefficient = as.numeric(coef_mat[, 1]),
  stringsAsFactors = FALSE
)

coef_df$abs_coefficient <- abs(coef_df$coefficient)

nonzero_coef_df <- coef_df[
  coef_df$abs_coefficient > 0,
]

nonzero_coef_df <- nonzero_coef_df[
  order(-nonzero_coef_df$abs_coefficient),
]

cat("\n===== Non-zero coefficients =====\n")
print(nonzero_coef_df)

nonzero_features <- setdiff(nonzero_coef_df$feature, "(Intercept)")

############################################################
# 6. Feature distribution shift diagnosis
############################################################

center <- frozen_model$center
scale_vec <- frozen_model$scale

center <- center[model_features]
scale_vec <- scale_vec[model_features]
scale_vec[is.na(scale_vec) | scale_vec == 0] <- 1

feature_shift_rows <- list()

for (ff in model_features) {
  
  tr <- suppressWarnings(as.numeric(train_df[[ff]]))
  lu <- suppressWarnings(as.numeric(lung_feat[[ff]]))
  
  tr <- tr[is.finite(tr)]
  lu <- lu[is.finite(lu)]
  
  tr_min <- min(tr, na.rm = TRUE)
  tr_max <- max(tr, na.rm = TRUE)
  tr_q01 <- as.numeric(quantile(tr, 0.01, na.rm = TRUE))
  tr_q99 <- as.numeric(quantile(tr, 0.99, na.rm = TRUE))
  
  lu_min <- min(lu, na.rm = TRUE)
  lu_max <- max(lu, na.rm = TRUE)
  
  lu_z_by_train <- (lu - center[ff]) / scale_vec[ff]
  
  outside_range_n <- sum(lu < tr_min | lu > tr_max, na.rm = TRUE)
  below_train_min_n <- sum(lu < tr_min, na.rm = TRUE)
  above_train_max_n <- sum(lu > tr_max, na.rm = TRUE)
  
  outside_q01q99_n <- sum(lu < tr_q01 | lu > tr_q99, na.rm = TRUE)
  
  feature_shift_rows[[length(feature_shift_rows) + 1]] <- data.frame(
    feature = ff,
    is_nonzero_model_feature = ff %in% nonzero_features,
    coefficient = ifelse(ff %in% coef_df$feature, coef_df$coefficient[match(ff, coef_df$feature)], 0),
    
    train_n = length(tr),
    lung1_n = length(lu),
    
    train_min = tr_min,
    train_q01 = tr_q01,
    train_q25 = as.numeric(quantile(tr, 0.25, na.rm = TRUE)),
    train_median = median(tr, na.rm = TRUE),
    train_mean = mean(tr, na.rm = TRUE),
    train_q75 = as.numeric(quantile(tr, 0.75, na.rm = TRUE)),
    train_q99 = tr_q99,
    train_max = tr_max,
    train_sd = sd(tr, na.rm = TRUE),
    
    lung1_min = lu_min,
    lung1_q01 = as.numeric(quantile(lu, 0.01, na.rm = TRUE)),
    lung1_q25 = as.numeric(quantile(lu, 0.25, na.rm = TRUE)),
    lung1_median = median(lu, na.rm = TRUE),
    lung1_mean = mean(lu, na.rm = TRUE),
    lung1_q75 = as.numeric(quantile(lu, 0.75, na.rm = TRUE)),
    lung1_q99 = as.numeric(quantile(lu, 0.99, na.rm = TRUE)),
    lung1_max = lu_max,
    lung1_sd = sd(lu, na.rm = TRUE),
    
    lung1_z_by_train_min = min(lu_z_by_train, na.rm = TRUE),
    lung1_z_by_train_q01 = as.numeric(quantile(lu_z_by_train, 0.01, na.rm = TRUE)),
    lung1_z_by_train_median = median(lu_z_by_train, na.rm = TRUE),
    lung1_z_by_train_q99 = as.numeric(quantile(lu_z_by_train, 0.99, na.rm = TRUE)),
    lung1_z_by_train_max = max(lu_z_by_train, na.rm = TRUE),
    lung1_max_abs_z_by_train = max(abs(lu_z_by_train), na.rm = TRUE),
    
    outside_train_minmax_n = outside_range_n,
    outside_train_minmax_pct = outside_range_n / length(lu),
    below_train_min_n = below_train_min_n,
    above_train_max_n = above_train_max_n,
    outside_train_q01q99_n = outside_q01q99_n,
    outside_train_q01q99_pct = outside_q01q99_n / length(lu),
    
    stringsAsFactors = FALSE
  )
}

feature_shift_df <- do.call(rbind, feature_shift_rows)

feature_shift_df <- feature_shift_df[
  order(
    -feature_shift_df$is_nonzero_model_feature,
    -feature_shift_df$lung1_max_abs_z_by_train,
    -feature_shift_df$outside_train_minmax_pct
  ),
]

cat("\n===== Feature distribution shift summary: top shifted model features =====\n")
print(
  head(
    feature_shift_df[
      ,
      c(
        "feature",
        "is_nonzero_model_feature",
        "coefficient",
        "lung1_max_abs_z_by_train",
        "outside_train_minmax_n",
        "outside_train_minmax_pct",
        "below_train_min_n",
        "above_train_max_n"
      )
    ],
    30
  )
)

############################################################
# 7. Recalculate original RRS and contribution matrix
############################################################

x_lung <- lung_feat[, model_features, drop = FALSE]

for (ff in model_features) {
  x_lung[[ff]] <- suppressWarnings(as.numeric(x_lung[[ff]]))
}

x_lung_mat <- as.matrix(x_lung)

x_lung_scaled <- sweep(x_lung_mat, 2, center, "-")
x_lung_scaled <- sweep(x_lung_scaled, 2, scale_vec, "/")

rrs_original_recalc <- as.numeric(
  predict(
    frozen_model$final_cvfit,
    newx = x_lung_scaled,
    s = frozen_model$final_lambda
  )
)

beta <- coef_df$coefficient[match(model_features, coef_df$feature)]
beta[is.na(beta)] <- 0
names(beta) <- model_features

contrib_mat <- sweep(x_lung_scaled, 2, beta, "*")
colnames(contrib_mat) <- model_features

contrib_df <- data.frame(
  patient_id = lung_feat$patient_id,
  Lung1_RRS_recalc = rrs_original_recalc,
  stringsAsFactors = FALSE
)

contrib_df <- cbind(
  contrib_df,
  as.data.frame(contrib_mat, check.names = FALSE)
)

############################################################
# 8. Identify extreme RRS patients
############################################################

rrs_q1 <- quantile(rrs_original_recalc, 0.25, na.rm = TRUE)
rrs_q3 <- quantile(rrs_original_recalc, 0.75, na.rm = TRUE)
rrs_iqr <- rrs_q3 - rrs_q1

rrs_lower_3iqr <- rrs_q1 - 3 * rrs_iqr
rrs_upper_3iqr <- rrs_q3 + 3 * rrs_iqr

rrs_z <- as.numeric(scale(rrs_original_recalc))

rrs_audit <- data.frame(
  patient_id = lung_feat$patient_id,
  Lung1_RRS_recalc = rrs_original_recalc,
  RRS_z = rrs_z,
  extreme_3IQR = rrs_original_recalc < rrs_lower_3iqr | rrs_original_recalc > rrs_upper_3iqr,
  extreme_absz5 = abs(rrs_z) > 5,
  stringsAsFactors = FALSE
)

extreme_ids <- rrs_audit$patient_id[
  rrs_audit$extreme_3IQR | rrs_audit$extreme_absz5
]

extreme_df <- rrs_audit[
  rrs_audit$patient_id %in% extreme_ids,
]

extreme_df <- extreme_df[
  order(extreme_df$Lung1_RRS_recalc),
]

cat("\n===== Extreme RRS audit from recalculated frozen model =====\n")
print(extreme_df)

############################################################
# 9. Contribution diagnosis for extreme patients
############################################################

top_contrib_rows <- list()

for (pid in extreme_ids) {
  
  idx <- which(contrib_df$patient_id == pid)
  
  if (length(idx) != 1) next
  
  vals <- as.numeric(contrib_df[idx, model_features])
  names(vals) <- model_features
  
  vals_nonzero <- vals[names(vals) %in% nonzero_features]
  
  ord <- order(abs(vals_nonzero), decreasing = TRUE)
  vals_sorted <- vals_nonzero[ord]
  
  top_n <- min(10, length(vals_sorted))
  
  for (j in seq_len(top_n)) {
    ff <- names(vals_sorted)[j]
    
    raw_value <- as.numeric(lung_feat[idx, ff])
    train_z <- (raw_value - center[ff]) / scale_vec[ff]
    
    top_contrib_rows[[length(top_contrib_rows) + 1]] <- data.frame(
      patient_id = pid,
      Lung1_RRS_recalc = contrib_df$Lung1_RRS_recalc[idx],
      rank_abs_contribution = j,
      feature = ff,
      coefficient = beta[ff],
      raw_feature_value = raw_value,
      train_center = center[ff],
      train_scale = scale_vec[ff],
      z_by_training = train_z,
      contribution = vals_sorted[j],
      abs_contribution = abs(vals_sorted[j]),
      train_min = feature_shift_df$train_min[match(ff, feature_shift_df$feature)],
      train_max = feature_shift_df$train_max[match(ff, feature_shift_df$feature)],
      lung1_min = feature_shift_df$lung1_min[match(ff, feature_shift_df$feature)],
      lung1_max = feature_shift_df$lung1_max[match(ff, feature_shift_df$feature)],
      outside_training_range = raw_value < feature_shift_df$train_min[match(ff, feature_shift_df$feature)] |
        raw_value > feature_shift_df$train_max[match(ff, feature_shift_df$feature)],
      stringsAsFactors = FALSE
    )
  }
}

top_contrib_df <- do.call(rbind, top_contrib_rows)

cat("\n===== Top contribution features in extreme RRS patients =====\n")
print(head(top_contrib_df, 80))

driver_summary <- aggregate(
  abs_contribution ~ feature,
  data = top_contrib_df,
  FUN = function(x) c(n = length(x), mean_abs = mean(x), max_abs = max(x))
)

driver_summary_expanded <- data.frame(
  feature = driver_summary$feature,
  n_times_in_top_contrib = driver_summary$abs_contribution[, "n"],
  mean_abs_contribution = driver_summary$abs_contribution[, "mean_abs"],
  max_abs_contribution = driver_summary$abs_contribution[, "max_abs"],
  stringsAsFactors = FALSE
)

driver_summary_expanded <- driver_summary_expanded[
  order(-driver_summary_expanded$max_abs_contribution),
]

cat("\n===== Extreme RRS driver summary =====\n")
print(driver_summary_expanded)

############################################################
# 10. Exploratory capped-RRS sensitivity
############################################################

x_cap_train_minmax <- x_lung_mat

for (ff in model_features) {
  tr_min <- feature_shift_df$train_min[match(ff, feature_shift_df$feature)]
  tr_max <- feature_shift_df$train_max[match(ff, feature_shift_df$feature)]
  
  x_cap_train_minmax[, ff] <- pmin(
    pmax(x_cap_train_minmax[, ff], tr_min),
    tr_max
  )
}

x_cap_train_minmax_scaled <- sweep(x_cap_train_minmax, 2, center, "-")
x_cap_train_minmax_scaled <- sweep(x_cap_train_minmax_scaled, 2, scale_vec, "/")

rrs_cap_train_minmax <- as.numeric(
  predict(
    frozen_model$final_cvfit,
    newx = x_cap_train_minmax_scaled,
    s = frozen_model$final_lambda
  )
)

x_scaled_cap_z3 <- x_lung_scaled
x_scaled_cap_z3[x_scaled_cap_z3 > 3] <- 3
x_scaled_cap_z3[x_scaled_cap_z3 < -3] <- -3

rrs_cap_z3 <- as.numeric(
  predict(
    frozen_model$final_cvfit,
    newx = x_scaled_cap_z3,
    s = frozen_model$final_lambda
  )
)

rrs_sensitivity_df <- data.frame(
  patient_id = lung_feat$patient_id,
  RRS_original = rrs_original_recalc,
  RRS_cap_train_minmax = rrs_cap_train_minmax,
  RRS_cap_z3 = rrs_cap_z3,
  stringsAsFactors = FALSE
)

rrs_sensitivity_df$RRS_original_z <- as.numeric(scale(rrs_sensitivity_df$RRS_original))
rrs_sensitivity_df$RRS_cap_train_minmax_z <- as.numeric(scale(rrs_sensitivity_df$RRS_cap_train_minmax))
rrs_sensitivity_df$RRS_cap_z3_z <- as.numeric(scale(rrs_sensitivity_df$RRS_cap_z3))

cat("\n===== Correlation among original and capped RRS =====\n")
print(
  cor(
    rrs_sensitivity_df[
      ,
      c("RRS_original", "RRS_cap_train_minmax", "RRS_cap_z3")
    ],
    use = "pairwise.complete.obs",
    method = "spearman"
  )
)

cat("\n===== RRS sensitivity score summaries =====\n")
print(summary(rrs_sensitivity_df[, c("RRS_original", "RRS_cap_train_minmax", "RRS_cap_z3")]))

############################################################
# 11. Exploratory survival test for capped RRS
############################################################

surv_df <- merge(
  analysis_df,
  rrs_sensitivity_df,
  by = "patient_id",
  all = FALSE
)

surv_df$OS_time_months <- suppressWarnings(as.numeric(surv_df$OS_time_months))
surv_df$OS_event <- suppressWarnings(as.numeric(surv_df$OS_event))

surv_df <- surv_df[
  is.finite(surv_df$OS_time_months) &
    is.finite(surv_df$OS_event) &
    surv_df$OS_time_months > 0 &
    surv_df$OS_event %in% c(0, 1),
]

run_survival_one_score <- function(data, score_col) {
  
  score <- data[[score_col]]
  med <- median(score, na.rm = TRUE)
  
  group <- ifelse(score >= med, "high", "low")
  group <- factor(group, levels = c("low", "high"))
  
  tmp <- data
  tmp$score_group <- group
  tmp$score_z <- as.numeric(scale(score))
  tmp$score_rank_z <- as.numeric(scale(rank(score, ties.method = "average")))
  
  lr <- survdiff(
    Surv(OS_time_months, OS_event) ~ score_group,
    data = tmp
  )
  
  lr_p <- pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)
  
  fit_bin <- coxph(
    Surv(OS_time_months, OS_event) ~ score_group,
    data = tmp
  )
  
  fit_rank <- coxph(
    Surv(OS_time_months, OS_event) ~ score_rank_z,
    data = tmp
  )
  
  s_bin <- summary(fit_bin)
  s_rank <- summary(fit_rank)
  
  data.frame(
    score = score_col,
    n = nrow(tmp),
    events = sum(tmp$OS_event == 1, na.rm = TRUE),
    median_cutoff = med,
    low_n = sum(tmp$score_group == "low"),
    high_n = sum(tmp$score_group == "high"),
    logrank_chisq = lr$chisq,
    logrank_p = lr_p,
    
    HR_high_vs_low = s_bin$conf.int[1, "exp(coef)"],
    CI_lower_high_vs_low = s_bin$conf.int[1, "lower .95"],
    CI_upper_high_vs_low = s_bin$conf.int[1, "upper .95"],
    p_high_vs_low = s_bin$coefficients[1, "Pr(>|z|)"],
    
    HR_rank_z = s_rank$conf.int[1, "exp(coef)"],
    CI_lower_rank_z = s_rank$conf.int[1, "lower .95"],
    CI_upper_rank_z = s_rank$conf.int[1, "upper .95"],
    p_rank_z = s_rank$coefficients[1, "Pr(>|z|)"],
    
    stringsAsFactors = FALSE
  )
}

survival_sensitivity <- do.call(
  rbind,
  lapply(
    c("RRS_original", "RRS_cap_train_minmax", "RRS_cap_z3"),
    function(cc) run_survival_one_score(surv_df, cc)
  )
)

cat("\n===== Exploratory survival sensitivity for capped RRS =====\n")
print(survival_sensitivity)

############################################################
# 12. Manuscript-oriented diagnostic summary
############################################################

diagnostic_summary <- data.frame(
  item = c(
    "training_patients",
    "lung1_patients_with_features",
    "model_features",
    "nonzero_model_features_excluding_intercept",
    "model_features_with_any_Lung1_outside_training_minmax",
    "model_features_with_more_than_5pct_Lung1_outside_training_minmax",
    "nonzero_model_features_with_any_Lung1_outside_training_minmax",
    "max_Lung1_abs_z_by_training_among_model_features",
    "n_RRS_extreme_3IQR_or_absz5",
    "original_RRS_min",
    "original_RRS_max",
    "cap_train_minmax_RRS_min",
    "cap_train_minmax_RRS_max",
    "cap_z3_RRS_min",
    "cap_z3_RRS_max",
    "original_RRS_logrank_p",
    "cap_train_minmax_logrank_p",
    "cap_z3_logrank_p"
  ),
  value = c(
    nrow(train_df),
    nrow(lung_feat),
    length(model_features),
    length(nonzero_features),
    sum(feature_shift_df$outside_train_minmax_n > 0),
    sum(feature_shift_df$outside_train_minmax_pct > 0.05),
    sum(feature_shift_df$outside_train_minmax_n[feature_shift_df$is_nonzero_model_feature] > 0),
    max(feature_shift_df$lung1_max_abs_z_by_train, na.rm = TRUE),
    length(extreme_ids),
    min(rrs_sensitivity_df$RRS_original, na.rm = TRUE),
    max(rrs_sensitivity_df$RRS_original, na.rm = TRUE),
    min(rrs_sensitivity_df$RRS_cap_train_minmax, na.rm = TRUE),
    max(rrs_sensitivity_df$RRS_cap_train_minmax, na.rm = TRUE),
    min(rrs_sensitivity_df$RRS_cap_z3, na.rm = TRUE),
    max(rrs_sensitivity_df$RRS_cap_z3, na.rm = TRUE),
    survival_sensitivity$logrank_p[survival_sensitivity$score == "RRS_original"],
    survival_sensitivity$logrank_p[survival_sensitivity$score == "RRS_cap_train_minmax"],
    survival_sensitivity$logrank_p[survival_sensitivity$score == "RRS_cap_z3"]
  ),
  stringsAsFactors = FALSE
)

cat("\n===== Diagnostic summary =====\n")
print(diagnostic_summary)

############################################################
# 13. Save outputs
############################################################

write.csv(
  feature_shift_df,
  file = file.path(out_dir, "feature_distribution_shift.csv"),
  row.names = FALSE
)

write.csv(
  nonzero_coef_df,
  file = file.path(out_dir, "primary_model_nonzero_coefficients.csv"),
  row.names = FALSE
)

write.csv(
  extreme_df,
  file = file.path(out_dir, "lung1_extreme_rrs_patients.csv"),
  row.names = FALSE
)

write.csv(
  top_contrib_df,
  file = file.path(out_dir, "extreme_rrs_feature_contributions.csv"),
  row.names = FALSE
)

write.csv(
  driver_summary_expanded,
  file = file.path(out_dir, "extreme_rrs_driver_summary.csv"),
  row.names = FALSE
)

write.csv(
  rrs_sensitivity_df,
  file = file.path(out_dir, "rrs_capping_sensitivity_scores.csv"),
  row.names = FALSE
)

write.csv(
  survival_sensitivity,
  file = file.path(out_dir, "rrs_capping_survival_sensitivity.csv"),
  row.names = FALSE
)

write.csv(
  diagnostic_summary,
  file = file.path(out_dir, "domain_shift_summary.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    feature_shift_df = feature_shift_df,
    nonzero_coef_df = nonzero_coef_df,
    extreme_df = extreme_df,
    top_contrib_df = top_contrib_df,
    driver_summary = driver_summary_expanded,
    rrs_sensitivity_df = rrs_sensitivity_df,
    survival_sensitivity = survival_sensitivity,
    diagnostic_summary = diagnostic_summary
  ),
  file = file.path(out_dir, "domain_shift_results.rds")
)

wb <- createWorkbook()

addWorksheet(wb, "diagnostic_summary")
writeData(wb, "diagnostic_summary", diagnostic_summary)

addWorksheet(wb, "feature_shift")
writeData(wb, "feature_shift", feature_shift_df)

addWorksheet(wb, "nonzero_coefficients")
writeData(wb, "nonzero_coefficients", nonzero_coef_df)

addWorksheet(wb, "extreme_RRS_patients")
writeData(wb, "extreme_RRS_patients", extreme_df)

addWorksheet(wb, "extreme_contributions")
writeData(wb, "extreme_contributions", top_contrib_df)

addWorksheet(wb, "driver_summary")
writeData(wb, "driver_summary", driver_summary_expanded)

addWorksheet(wb, "RRS_sensitivity_scores")
writeData(wb, "RRS_sensitivity_scores", rrs_sensitivity_df)

addWorksheet(wb, "survival_sensitivity")
writeData(wb, "survival_sensitivity", survival_sensitivity)

saveWorkbook(
  wb,
  file = file.path(out_dir, "lung1_domain_shift_summary.xlsx"),
  overwrite = TRUE
)

############################################################
# 14. Figures
############################################################

png(
  filename = file.path(fig_dir, "rrs_original_vs_capped.png"),
  width = 1600,
  height = 1300,
  res = 200
)

plot(
  rrs_sensitivity_df$RRS_original,
  rrs_sensitivity_df$RRS_cap_train_minmax,
  xlab = "Original Lung1 RRS",
  ylab = "Training-range-capped RRS",
  main = "Original vs capped external RRS",
  pch = 19
)

abline(lm(RRS_cap_train_minmax ~ RRS_original, data = rrs_sensitivity_df), lwd = 2)

dev.off()

top_shift_plot <- head(
  feature_shift_df[
    order(-feature_shift_df$lung1_max_abs_z_by_train),
  ],
  20
)

png(
  filename = file.path(fig_dir, "top_feature_distribution_shift.png"),
  width = 1800,
  height = 1400,
  res = 200
)

barplot(
  top_shift_plot$lung1_max_abs_z_by_train,
  names.arg = top_shift_plot$feature,
  las = 2,
  cex.names = 0.55,
  ylab = "Max absolute Lung1 z-score using training mean/SD",
  main = "Top 20 model features with largest external distribution shift"
)

dev.off()

############################################################
# 15. Done
############################################################

cat("\n===== Lung1 RRS domain-shift diagnostics completed =====\n")
cat("Main output files:\n")

cat("1. Non-zero coefficients\n")
cat("2. Feature distribution shift summary: top shifted model features\n")
cat("3. Extreme RRS audit from recalculated frozen model\n")
cat("4. Top contribution features in extreme RRS patients\n")
cat("5. Extreme RRS driver summary\n")
cat("6. Exploratory survival sensitivity for capped RRS\n")
cat("7. Diagnostic summary\n")
