# 02_analyze_tumor_size_adjustment.R
#
# Audit the incremental association of primary RRS beyond tumor size.
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
# 1. Paths
############################################################


training_size_file <- file.path(
  project_dir,
  "07_results",
  "tumor_size_domain_shift",
  "training_size_analysis_dataset.csv"
)

strict_audit_file <- file.path(
  project_dir,
  "07_results",
  "rrs_permutation",
  "rrs_full_workflow_permutation.xlsx"
)

out_dir <- file.path(
  project_dir,
  "07_results",
  "rrs_size_adjustment"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(training_size_file)) {
  stop("Cannot find training size-analysis dataset:\n", training_size_file)
}
if (!file.exists(strict_audit_file)) {
  stop("Cannot find RRS audit workbook:\n", strict_audit_file)
}

############################################################
# 2. Read and merge exact patient-level data
############################################################

size_df <- read.csv(
  training_size_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

strict_pred <- read.xlsx(
  strict_audit_file,
  sheet = "observed_patient_predictions"
)

required_size <- c(
  "patient_id",
  "CRS_z",
  "log_MeshVolume",
  "log_Maximum3DDiameter",
  "SignedLog_RRS_CV"
)
required_strict <- c(
  "patient_id",
  "observed_CRS_z",
  "RRS_CV_mean"
)

missing_size <- setdiff(required_size, colnames(size_df))
missing_strict <- setdiff(required_strict, colnames(strict_pred))

if (length(missing_size) > 0) {
  stop("Missing columns in size dataset: ", paste(missing_size, collapse = ", "))
}
if (length(missing_strict) > 0) {
  stop("Missing columns in repeated-CV predictions: ", paste(missing_strict, collapse = ", "))
}

size_df$patient_id <- toupper(trimws(as.character(size_df$patient_id)))
strict_pred$patient_id <- toupper(trimws(as.character(strict_pred$patient_id)))

dat <- merge(
  size_df,
  strict_pred[, required_strict],
  by = "patient_id",
  all = FALSE
)

numeric_cols <- c(
  "CRS_z",
  "observed_CRS_z",
  "log_MeshVolume",
  "log_Maximum3DDiameter",
  "RRS_CV_mean",
  "SignedLog_RRS_CV"
)

for (nm in numeric_cols) {
  dat[[nm]] <- suppressWarnings(as.numeric(dat[[nm]]))
}

dat <- dat[
  complete.cases(dat[, numeric_cols]),
  ,
  drop = FALSE
]

if (nrow(dat) != 117) {
  stop("Expected 117 complete matched patients, found ", nrow(dat), ".")
}

max_outcome_difference <- max(
  abs(dat$CRS_z - dat$observed_CRS_z),
  na.rm = TRUE
)

if (!is.finite(max_outcome_difference) || max_outcome_difference > 1e-8) {
  stop(
    "CRS_z mismatch between the size dataset and repeated-CV predictions. Max difference = ",
    max_outcome_difference
  )
}

dat$RRS_CV <- dat$RRS_CV_mean
dat$RRS_CV_z <- as.numeric(scale(dat$RRS_CV))
dat$SignedLog_RRS_CV_z <- as.numeric(scale(dat$SignedLog_RRS_CV))

############################################################
# 3. Helpers
############################################################

extract_score_row <- function(fit, score_term, model_name) {
  sm <- summary(fit)
  cf <- sm$coefficients
  
  if (!score_term %in% rownames(cf)) {
    stop("Score term not found in model: ", score_term)
  }
  
  ci <- confint(fit, parm = score_term)
  
  data.frame(
    model = model_name,
    n = nobs(fit),
    beta = unname(cf[score_term, "Estimate"]),
    se = unname(cf[score_term, "Std. Error"]),
    ci_low = unname(ci[1]),
    ci_high = unname(ci[2]),
    p_value = unname(cf[score_term, "Pr(>|t|)"]),
    model_R2 = summary(fit)$r.squared,
    adjusted_R2 = summary(fit)$adj.r.squared,
    stringsAsFactors = FALSE
  )
}

partial_cor <- function(data, outcome, score, covariates, label) {
  cc <- complete.cases(data[, c(outcome, score, covariates)])
  d <- data[cc, c(outcome, score, covariates), drop = FALSE]
  
  cov_formula <- as.formula(
    paste("~", paste(covariates, collapse = " + "))
  )
  
  y_resid <- residuals(
    lm(
      as.formula(
        paste(outcome, paste(covariates, collapse = " + "), sep = " ~ ")
      ),
      data = d
    )
  )
  
  x_resid <- residuals(
    lm(
      as.formula(
        paste(score, paste(covariates, collapse = " + "), sep = " ~ ")
      ),
      data = d
    )
  )
  
  pearson <- cor.test(y_resid, x_resid, method = "pearson")
  spearman <- suppressWarnings(
    cor.test(y_resid, x_resid, method = "spearman", exact = FALSE)
  )
  
  data.frame(
    label = label,
    outcome = outcome,
    score = score,
    adjusted_for = paste(covariates, collapse = "+"),
    n = nrow(d),
    partial_Pearson_r = unname(pearson$estimate),
    partial_Pearson_p = pearson$p.value,
    partial_Spearman_r = unname(spearman$estimate),
    partial_Spearman_p = spearman$p.value,
    stringsAsFactors = FALSE
  )
}

############################################################
# 4. Primary RRS size-adjusted models
############################################################

fit_size_only <- lm(
  CRS_z ~ log_MeshVolume + log_Maximum3DDiameter,
  data = dat
)

fit_size_strict <- lm(
  CRS_z ~ log_MeshVolume + log_Maximum3DDiameter + RRS_CV_z,
  data = dat
)

# Existing signed-log sensitivity model is retained as a sensitivity analysis.
fit_size_signedlog <- lm(
  CRS_z ~ log_MeshVolume + log_Maximum3DDiameter + SignedLog_RRS_CV_z,
  data = dat
)

strict_row <- extract_score_row(
  fit_size_strict,
  "RRS_CV_z",
  "Size_plus_strict_RRS"
)

signedlog_row <- extract_score_row(
  fit_size_signedlog,
  "SignedLog_RRS_CV_z",
  "Size_plus_signed_log_RRS_sensitivity"
)

nested_strict <- anova(fit_size_only, fit_size_strict)
nested_signedlog <- anova(fit_size_only, fit_size_signedlog)

size_summary <- data.frame(
  model = c(
    "Size_only",
    "Size_plus_strict_RRS",
    "Size_plus_signed_log_RRS_sensitivity"
  ),
  n = c(
    nobs(fit_size_only),
    nobs(fit_size_strict),
    nobs(fit_size_signedlog)
  ),
  model_R2 = c(
    summary(fit_size_only)$r.squared,
    summary(fit_size_strict)$r.squared,
    summary(fit_size_signedlog)$r.squared
  ),
  adjusted_R2 = c(
    summary(fit_size_only)$adj.r.squared,
    summary(fit_size_strict)$adj.r.squared,
    summary(fit_size_signedlog)$adj.r.squared
  ),
  delta_R2_vs_size_only = c(
    NA_real_,
    summary(fit_size_strict)$r.squared - summary(fit_size_only)$r.squared,
    summary(fit_size_signedlog)$r.squared - summary(fit_size_only)$r.squared
  ),
  added_score_p = c(
    NA_real_,
    nested_strict$`Pr(>F)`[2],
    nested_signedlog$`Pr(>F)`[2]
  ),
  stringsAsFactors = FALSE
)

adjusted_score_results <- rbind(strict_row, signedlog_row)

partial_results <- rbind(
  partial_cor(
    dat,
    outcome = "CRS_z",
    score = "RRS_CV",
    covariates = c("log_MeshVolume", "log_Maximum3DDiameter"),
    label = "RRS_CV_partial_adjusted_for_size"
  ),
  partial_cor(
    dat,
    outcome = "CRS_z",
    score = "SignedLog_RRS_CV",
    covariates = c("log_MeshVolume", "log_Maximum3DDiameter"),
    label = "SignedLog_RRS_CV_partial_adjusted_for_size"
  )
)

############################################################
# 5. Key analysis results
############################################################

final_values <- data.frame(
  metric = c(
    "n",
    "size_only_R2",
    "size_plus_strict_RRS_R2",
    "strict_RRS_delta_R2",
    "strict_RRS_standardized_beta",
    "strict_RRS_beta_CI_low",
    "strict_RRS_beta_CI_high",
    "strict_RRS_p",
    "strict_RRS_partial_Pearson_r",
    "strict_RRS_partial_Pearson_p",
    "strict_RRS_partial_Spearman_r",
    "strict_RRS_partial_Spearman_p",
    "size_plus_signed_log_RRS_R2",
    "signed_log_RRS_delta_R2",
    "signed_log_RRS_standardized_beta",
    "signed_log_RRS_p"
  ),
  value = c(
    nrow(dat),
    summary(fit_size_only)$r.squared,
    summary(fit_size_strict)$r.squared,
    summary(fit_size_strict)$r.squared - summary(fit_size_only)$r.squared,
    strict_row$beta,
    strict_row$ci_low,
    strict_row$ci_high,
    strict_row$p_value,
    partial_results$partial_Pearson_r[1],
    partial_results$partial_Pearson_p[1],
    partial_results$partial_Spearman_r[1],
    partial_results$partial_Spearman_p[1],
    summary(fit_size_signedlog)$r.squared,
    summary(fit_size_signedlog)$r.squared - summary(fit_size_only)$r.squared,
    signedlog_row$beta,
    signedlog_row$p_value
  ),
  stringsAsFactors = FALSE
)

############################################################
# 6. Save outputs
############################################################

write.csv(
  dat,
  file.path(out_dir, "rrs_size_adjustment_dataset.csv"),
  row.names = FALSE
)

write.csv(
  size_summary,
  file.path(out_dir, "rrs_size_adjustment_nested_models.csv"),
  row.names = FALSE
)

write.csv(
  adjusted_score_results,
  file.path(out_dir, "rrs_size_adjustment_coefficients.csv"),
  row.names = FALSE
)

write.csv(
  partial_results,
  file.path(out_dir, "rrs_size_adjustment_partial_correlations.csv"),
  row.names = FALSE
)

write.csv(
  final_values,
  file.path(out_dir, "rrs_size_adjustment_key_values.csv"),
  row.names = FALSE
)

wb <- createWorkbook()

addWorksheet(wb, "final_values")
writeData(wb, "final_values", final_values)

addWorksheet(wb, "nested_models")
writeData(wb, "nested_models", size_summary)

addWorksheet(wb, "score_coefficients")
writeData(wb, "score_coefficients", adjusted_score_results)

addWorksheet(wb, "partial_correlations")
writeData(wb, "partial_correlations", partial_results)

addWorksheet(wb, "analysis_dataset")
writeData(wb, "analysis_dataset", dat)

addWorksheet(wb, "method_note")
writeData(
  wb,
  "method_note",
  data.frame(
    note = c(
      "RRS_CV was read from the RRS full-workflow permutation analysis observed_patient_predictions sheet.",
      "All 117 patients were matched exactly by patient_id.",
      "CRS_z values in the size dataset and repeated-CV predictions were required to match within 1e-8.",
      "Primary RRS and signed-log RRS were standardized before adjusted linear modeling; beta therefore represents change in CRS_z per 1 SD increase in the RRS.",
      "Size covariates were log_MeshVolume and log_Maximum3DDiameter.",
      "The signed-log model is a sensitivity analysis and is not a second primary RRS model."
    ),
    stringsAsFactors = FALSE
  )
)

for (sh in names(wb)) {
  setColWidths(wb, sh, cols = 1:40, widths = "auto")
}

audit_xlsx <- file.path(
  out_dir,
  "rrs_size_adjustment_summary.xlsx"
)

saveWorkbook(wb, audit_xlsx, overwrite = TRUE)

txt <- c(
  "===== RRS TUMOR-SIZE ADJUSTMENT =====",
  paste0("Patients: ", nrow(dat)),
  paste0("Size-only R2: ", sprintf("%.10f", summary(fit_size_only)$r.squared)),
  paste0("Size + primary RRS R2: ", sprintf("%.10f", summary(fit_size_strict)$r.squared)),
  paste0(
    "Primary RRS delta R2: ",
    sprintf(
      "%.10f",
      summary(fit_size_strict)$r.squared - summary(fit_size_only)$r.squared
    )
  ),
  paste0("Primary RRS standardized beta: ", sprintf("%.10f", strict_row$beta)),
  paste0(
    "Primary RRS beta 95% CI: ",
    sprintf("%.10f", strict_row$ci_low),
    " to ",
    sprintf("%.10f", strict_row$ci_high)
  ),
  paste0("Primary RRS p: ", sprintf("%.10f", strict_row$p_value)),
  paste0(
    "Strict partial Pearson r: ",
    sprintf("%.10f", partial_results$partial_Pearson_r[1]),
    "; p = ",
    sprintf("%.10f", partial_results$partial_Pearson_p[1])
  ),
  paste0(
    "Strict partial Spearman rho: ",
    sprintf("%.10f", partial_results$partial_Spearman_r[1]),
    "; p = ",
    sprintf("%.10f", partial_results$partial_Spearman_p[1])
  ),
  paste0(
    "Size + signed-log sensitivity R2: ",
    sprintf("%.10f", summary(fit_size_signedlog)$r.squared)
  ),
  paste0(
    "Signed-log sensitivity delta R2: ",
    sprintf(
      "%.10f",
      summary(fit_size_signedlog)$r.squared - summary(fit_size_only)$r.squared
    )
  ),
  paste0("Signed-log standardized beta: ", sprintf("%.10f", signedlog_row$beta)),
  paste0("Signed-log p: ", sprintf("%.10f", signedlog_row$p_value)),
  "",
  paste0("Audit workbook: ", audit_xlsx)
)

txt_file <- file.path(
  out_dir,
  "rrs_size_adjustment_key_values.txt"
)

writeLines(txt, txt_file, useBytes = TRUE)
cat(paste(txt, collapse = "\n"))
cat("\n")
