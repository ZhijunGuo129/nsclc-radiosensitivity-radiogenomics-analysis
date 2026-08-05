# 05_train_signed_log_sensitivity_model.R
#
# Train the signed-log feature-transformed RRS sensitivity model and evaluate Lung1 OS.
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

lung1_feature_log_csv <- file.path(
  lung1_root,
  "radiomics_features",
  "lung1_radiomics_features_log.csv"
)

clinical_rds <- file.path(
  lung1_root,
  "clinical",
  "Lung1_clinical_clean_initial.rds"
)

out_dir <- file.path(
  lung1_root,
  "signed_log_rrs"
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

lung_log <- read.csv(
  lung1_feature_log_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

clinical <- readRDS(clinical_rds)

cat("Training table dim:\n")
print(dim(train_df))

cat("\nLung1 feature table dim:\n")
print(dim(lung_feat))

cat("\nLung1 feature extraction status:\n")
print(table(lung_log$status, useNA = "ifany"))

cat("\nClinical table dim:\n")
print(dim(clinical))

############################################################
# 4. ID and outcome preparation
############################################################

train_df$patient_id <- toupper(trimws(as.character(train_df$patient_id)))
lung_feat$patient_id <- toupper(trimws(as.character(lung_feat$patient_id)))

if (!"patient_id" %in% colnames(clinical)) {
  if ("PatientID" %in% colnames(clinical)) {
    clinical$patient_id <- toupper(trimws(as.character(clinical$PatientID)))
  } else {
    stop("Cannot find patient_id or PatientID in clinical table.")
  }
}

clinical$patient_id <- toupper(trimws(as.character(clinical$patient_id)))

outcome_col <- "CRS_z"

if (!outcome_col %in% colnames(train_df)) {
  stop("Cannot find CRS_z in training table.")
}

train_df[[outcome_col]] <- suppressWarnings(as.numeric(train_df[[outcome_col]]))

train_df <- train_df[
  is.finite(train_df[[outcome_col]]),
]

cat("\nTraining outcome CRS_z summary:\n")
print(summary(train_df[[outcome_col]]))

############################################################
# 5. Signed-log transform
############################################################

signed_log1p <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  sign(x) * log1p(abs(x))
}

train_feature_cols <- grep("^original_", colnames(train_df), value = TRUE)
lung_feature_cols <- grep("^original_", colnames(lung_feat), value = TRUE)

common_features <- intersect(train_feature_cols, lung_feature_cols)
common_features <- sort(common_features)

cat("\n===== Radiomics feature audit =====\n")
cat("Training original features:", length(train_feature_cols), "\n")
cat("Lung1 original features:", length(lung_feature_cols), "\n")
cat("Common original features:", length(common_features), "\n")

if (length(common_features) < 50) {
  stop("Too few common radiomics features.")
}

x_train_raw <- train_df[, common_features, drop = FALSE]
x_lung_raw <- lung_feat[, common_features, drop = FALSE]

x_train_log <- as.data.frame(
  lapply(x_train_raw, signed_log1p),
  check.names = FALSE
)

x_lung_log <- as.data.frame(
  lapply(x_lung_raw, signed_log1p),
  check.names = FALSE
)

# Remove non-finite and zero-variance features in the training data.
feature_qc <- data.frame(
  feature = common_features,
  train_missing_n = sapply(common_features, function(ff) sum(!is.finite(x_train_log[[ff]]))),
  lung_missing_n = sapply(common_features, function(ff) sum(!is.finite(x_lung_log[[ff]]))),
  train_sd = sapply(common_features, function(ff) sd(x_train_log[[ff]], na.rm = TRUE)),
  stringsAsFactors = FALSE
)

feature_qc$keep <- feature_qc$train_missing_n == 0 &
  feature_qc$lung_missing_n == 0 &
  is.finite(feature_qc$train_sd) &
  feature_qc$train_sd > 0

model_features <- feature_qc$feature[feature_qc$keep]

cat("\nFeature QC summary:\n")
print(table(feature_qc$keep, useNA = "ifany"))

cat("\nModel features after signed-log QC:\n")
print(length(model_features))

if (length(model_features) < 30) {
  stop("Too few model features after QC.")
}

x_train_log <- x_train_log[, model_features, drop = FALSE]
x_lung_log <- x_lung_log[, model_features, drop = FALSE]

y <- train_df[[outcome_col]]
n <- length(y)

############################################################
# 6. Helper functions
############################################################

calc_metrics <- function(y_true, y_pred) {
  ok <- is.finite(y_true) & is.finite(y_pred)
  y_true <- y_true[ok]
  y_pred <- y_pred[ok]
  
  if (length(y_true) < 5) {
    return(data.frame(
      n = length(y_true),
      Pearson = NA,
      Spearman = NA,
      MAE = NA,
      RMSE = NA,
      R2 = NA
    ))
  }
  
  pearson <- suppressWarnings(cor(y_true, y_pred, method = "pearson"))
  spearman <- suppressWarnings(cor(y_true, y_pred, method = "spearman"))
  mae <- mean(abs(y_true - y_pred))
  rmse <- sqrt(mean((y_true - y_pred)^2))
  r2 <- 1 - sum((y_true - y_pred)^2) / sum((y_true - mean(y_true))^2)
  
  data.frame(
    n = length(y_true),
    Pearson = pearson,
    Spearman = spearman,
    MAE = mae,
    RMSE = rmse,
    R2 = r2
  )
}

scale_by_train <- function(x_train, x_test) {
  center <- colMeans(x_train, na.rm = TRUE)
  scale_vec <- apply(x_train, 2, sd, na.rm = TRUE)
  scale_vec[is.na(scale_vec) | scale_vec == 0] <- 1
  
  x_train_s <- sweep(x_train, 2, center, "-")
  x_train_s <- sweep(x_train_s, 2, scale_vec, "/")
  
  x_test_s <- sweep(x_test, 2, center, "-")
  x_test_s <- sweep(x_test_s, 2, scale_vec, "/")
  
  list(
    x_train_s = x_train_s,
    x_test_s = x_test_s,
    center = center,
    scale = scale_vec
  )
}

extract_cox <- function(fit, model_name) {
  s <- summary(fit)
  coef_tab <- as.data.frame(s$coefficients)
  ci_tab <- as.data.frame(s$conf.int)
  
  out <- data.frame(
    model = model_name,
    variable = rownames(coef_tab),
    HR = ci_tab$`exp(coef)`,
    CI_lower = ci_tab$`lower .95`,
    CI_upper = ci_tab$`upper .95`,
    p_value = coef_tab$`Pr(>|z|)`,
    concordance = s$concordance[1],
    n = fit$n,
    nevent = fit$nevent,
    stringsAsFactors = FALSE
  )
  
  rownames(out) <- NULL
  out
}

get_logrank <- function(data, group_col = "SignedLog_RRS_group") {
  f <- as.formula(paste0("Surv(OS_time_months, OS_event) ~ ", group_col))
  lr <- survdiff(f, data = data)
  p <- pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)
  
  data.frame(
    test = "logrank",
    group_col = group_col,
    chisq = lr$chisq,
    df = length(lr$n) - 1,
    p_value = p,
    stringsAsFactors = FALSE
  )
}

get_survival_rates <- function(km_fit, times_months) {
  sm <- summary(km_fit, times = times_months)
  
  data.frame(
    time_months = sm$time,
    strata = sm$strata,
    n_risk = sm$n.risk,
    n_event = sm$n.event,
    survival = sm$surv,
    lower = sm$lower,
    upper = sm$upper,
    stringsAsFactors = FALSE
  )
}

valid_covariates <- function(data, candidate_covs) {
  keep <- character(0)
  
  for (cc in candidate_covs) {
    if (!cc %in% colnames(data)) next
    
    x <- data[[cc]]
    
    if (is.numeric(x)) {
      if (sum(is.finite(x)) > 5 && length(unique(x[is.finite(x)])) > 1) {
        keep <- c(keep, cc)
      }
    } else {
      x <- factor(x)
      if (length(unique(x[!is.na(x)])) > 1) {
        keep <- c(keep, cc)
      }
    }
  }
  
  keep
}

############################################################
# 7. Repeated CV signed-log model
############################################################

set.seed(20260726)

alpha_grid <- c(0.25, 0.5, 0.75, 1)
n_repeats <- 20
k_outer <- 5

x_train_mat <- as.matrix(x_train_log)

cv_pred_mat <- matrix(
  NA_real_,
  nrow = n,
  ncol = n_repeats
)

colnames(cv_pred_mat) <- paste0("repeat_", seq_len(n_repeats))
rownames(cv_pred_mat) <- train_df$patient_id

repeat_summary_list <- list()
alpha_summary_list <- list()

cat("\n===== Signed-log repeated cross-validation started =====\n")
cat("Repeats:", n_repeats, "\n")
cat("Outer folds:", k_outer, "\n")
cat("Alpha grid:", paste(alpha_grid, collapse = ", "), "\n")

for (rr in seq_len(n_repeats)) {
  
  cat("\nRepeat", rr, "/", n_repeats, "\n")
  
  folds <- sample(rep(seq_len(k_outer), length.out = n))
  
  pred_by_alpha <- matrix(
    NA_real_,
    nrow = n,
    ncol = length(alpha_grid)
  )
  
  colnames(pred_by_alpha) <- paste0("alpha_", alpha_grid)
  
  for (aa in seq_along(alpha_grid)) {
    
    alpha_now <- alpha_grid[aa]
    
    cat("  Alpha:", alpha_now, "\n")
    
    pred_now <- rep(NA_real_, n)
    
    for (ffold in seq_len(k_outer)) {
      
      train_idx <- which(folds != ffold)
      test_idx <- which(folds == ffold)
      
      x_tr <- x_train_mat[train_idx, , drop = FALSE]
      x_te <- x_train_mat[test_idx, , drop = FALSE]
      
      sc <- scale_by_train(x_tr, x_te)
      
      cvfit <- cv.glmnet(
        x = sc$x_train_s,
        y = y[train_idx],
        family = "gaussian",
        alpha = alpha_now,
        nfolds = 5,
        standardize = FALSE
      )
      
      pred_now[test_idx] <- as.numeric(
        predict(
          cvfit,
          newx = sc$x_test_s,
          s = "lambda.min"
        )
      )
    }
    
    pred_by_alpha[, aa] <- pred_now
    
    m <- calc_metrics(y, pred_now)
    
    alpha_summary_list[[length(alpha_summary_list) + 1]] <- cbind(
      data.frame(
        repeat_id = rr,
        alpha = alpha_now,
        stringsAsFactors = FALSE
      ),
      m
    )
  }
  
  alpha_metrics_this <- do.call(
    rbind,
    alpha_summary_list[
      sapply(alpha_summary_list, function(z) z$repeat_id[1]) == rr
    ]
  )
  
  best_row <- alpha_metrics_this[
    order(
      -alpha_metrics_this$Spearman,
      -alpha_metrics_this$Pearson,
      alpha_metrics_this$RMSE
    ),
  ][1, ]
  
  best_alpha <- best_row$alpha[1]
  best_alpha_col <- paste0("alpha_", best_alpha)
  
  cv_pred_mat[, rr] <- pred_by_alpha[, best_alpha_col]
  
  repeat_summary_list[[length(repeat_summary_list) + 1]] <- cbind(
    data.frame(
      repeat_id = rr,
      selected_alpha = best_alpha,
      stringsAsFactors = FALSE
    ),
    calc_metrics(y, cv_pred_mat[, rr])
  )
  
  cat("  Selected alpha:", best_alpha, "\n")
  cat("  Spearman:", round(best_row$Spearman, 4), "\n")
}

cv_repeat_summary <- do.call(rbind, repeat_summary_list)
cv_alpha_summary <- do.call(rbind, alpha_summary_list)

train_cv_pred <- data.frame(
  patient_id = train_df$patient_id,
  CRS_z = y,
  SignedLog_RRS_CV_mean = rowMeans(cv_pred_mat, na.rm = TRUE),
  stringsAsFactors = FALSE
)

train_cv_metrics_mean_pred <- calc_metrics(
  train_cv_pred$CRS_z,
  train_cv_pred$SignedLog_RRS_CV_mean
)

cv_metrics_summary <- data.frame(
  metric = c("Pearson", "Spearman", "MAE", "RMSE", "R2"),
  repeat_mean = c(
    mean(cv_repeat_summary$Pearson, na.rm = TRUE),
    mean(cv_repeat_summary$Spearman, na.rm = TRUE),
    mean(cv_repeat_summary$MAE, na.rm = TRUE),
    mean(cv_repeat_summary$RMSE, na.rm = TRUE),
    mean(cv_repeat_summary$R2, na.rm = TRUE)
  ),
  repeat_sd = c(
    sd(cv_repeat_summary$Pearson, na.rm = TRUE),
    sd(cv_repeat_summary$Spearman, na.rm = TRUE),
    sd(cv_repeat_summary$MAE, na.rm = TRUE),
    sd(cv_repeat_summary$RMSE, na.rm = TRUE),
    sd(cv_repeat_summary$R2, na.rm = TRUE)
  ),
  mean_prediction_metric = c(
    train_cv_metrics_mean_pred$Pearson,
    train_cv_metrics_mean_pred$Spearman,
    train_cv_metrics_mean_pred$MAE,
    train_cv_metrics_mean_pred$RMSE,
    train_cv_metrics_mean_pred$R2
  ),
  stringsAsFactors = FALSE
)

cat("\n===== Signed-log training CV summary =====\n")
print(cv_metrics_summary)

cat("\nSelected alpha table:\n")
print(table(cv_repeat_summary$selected_alpha, useNA = "ifany"))

############################################################
# 8. Train final signed-log model on all training patients
############################################################

center_full <- colMeans(x_train_mat, na.rm = TRUE)
scale_full <- apply(x_train_mat, 2, sd, na.rm = TRUE)
scale_full[is.na(scale_full) | scale_full == 0] <- 1

x_train_scaled_full <- sweep(x_train_mat, 2, center_full, "-")
x_train_scaled_full <- sweep(x_train_scaled_full, 2, scale_full, "/")

final_fit_list <- list()
final_fit_summary <- data.frame()

for (aa in alpha_grid) {
  
  set.seed(20260726 + round(aa * 100))
  
  fit <- cv.glmnet(
    x = x_train_scaled_full,
    y = y,
    family = "gaussian",
    alpha = aa,
    nfolds = 5,
    standardize = FALSE
  )
  
  final_fit_list[[as.character(aa)]] <- fit
  
  final_fit_summary <- rbind(
    final_fit_summary,
    data.frame(
      alpha = aa,
      lambda_min = fit$lambda.min,
      min_cvm = min(fit$cvm, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  )
}

final_fit_summary <- final_fit_summary[
  order(final_fit_summary$min_cvm),
]

final_alpha <- final_fit_summary$alpha[1]
final_cvfit <- final_fit_list[[as.character(final_alpha)]]
final_lambda <- final_cvfit$lambda.min

coef_mat <- as.matrix(
  coef(final_cvfit, s = final_lambda)
)

coef_df <- data.frame(
  feature = rownames(coef_mat),
  coefficient = as.numeric(coef_mat[, 1]),
  abs_coefficient = abs(as.numeric(coef_mat[, 1])),
  stringsAsFactors = FALSE
)

coef_df <- coef_df[
  order(-coef_df$abs_coefficient),
]

nonzero_coef_df <- coef_df[
  coef_df$abs_coefficient > 0,
]

cat("\n===== Final signed-log model =====\n")
cat("Final alpha:\n")
print(final_alpha)

cat("\nFinal lambda:\n")
print(final_lambda)

cat("\nFinal fit summary by alpha:\n")
print(final_fit_summary)

cat("\nNon-zero coefficients:\n")
print(nonzero_coef_df)

signed_log_model <- list(
  final_cvfit = final_cvfit,
  final_alpha = final_alpha,
  final_lambda = final_lambda,
  model_features = model_features,
  center = center_full,
  scale = scale_full,
  transform = "signed_log1p_abs: sign(x) * log1p(abs(x))",
  outcome_col = outcome_col,
  cv_metrics_summary = cv_metrics_summary,
  cv_repeat_summary = cv_repeat_summary,
  final_fit_summary = final_fit_summary,
  coefficients = coef_df,
  nonzero_coefficients = nonzero_coef_df,
  note = "Signed-log feature-transformed radiomics model trained to predict CRS_z"
)

saveRDS(
  signed_log_model,
  file = file.path(out_dir, "signed_log_rrs_model.rds")
)

############################################################
# 9. Apply signed-log model to Lung1
############################################################

x_lung_mat <- as.matrix(x_lung_log[, model_features, drop = FALSE])

x_lung_scaled <- sweep(x_lung_mat, 2, center_full, "-")
x_lung_scaled <- sweep(x_lung_scaled, 2, scale_full, "/")

lung_rrs <- as.numeric(
  predict(
    final_cvfit,
    newx = x_lung_scaled,
    s = final_lambda
  )
)

rrs_df <- data.frame(
  patient_id = lung_feat$patient_id,
  SignedLog_Lung1_RRS = lung_rrs,
  stringsAsFactors = FALSE
)

rrs_df$SignedLog_Lung1_RRS_z <- as.numeric(scale(rrs_df$SignedLog_Lung1_RRS))
rrs_median <- median(rrs_df$SignedLog_Lung1_RRS, na.rm = TRUE)

rrs_df$SignedLog_RRS_group <- ifelse(
  rrs_df$SignedLog_Lung1_RRS >= rrs_median,
  "RRS_high",
  "RRS_low"
)

rrs_df$SignedLog_RRS_group <- factor(
  rrs_df$SignedLog_RRS_group,
  levels = c("RRS_low", "RRS_high")
)

cat("\n===== Signed-log Lung1 RRS summary =====\n")
print(summary(rrs_df$SignedLog_Lung1_RRS))

cat("\nSigned-log RRS group table:\n")
print(table(rrs_df$SignedLog_RRS_group, useNA = "ifany"))

############################################################
# 10. Signed-log domain shift audit
############################################################

nonzero_features <- setdiff(nonzero_coef_df$feature, "(Intercept)")

shift_rows <- list()

for (ff in model_features) {
  
  tr <- x_train_log[[ff]]
  lu <- x_lung_log[[ff]]
  
  tr_min <- min(tr, na.rm = TRUE)
  tr_max <- max(tr, na.rm = TRUE)
  
  lu_z <- (lu - center_full[ff]) / scale_full[ff]
  
  outside_n <- sum(lu < tr_min | lu > tr_max, na.rm = TRUE)
  
  shift_rows[[length(shift_rows) + 1]] <- data.frame(
    feature = ff,
    is_nonzero_model_feature = ff %in% nonzero_features,
    coefficient = ifelse(
      ff %in% coef_df$feature,
      coef_df$coefficient[match(ff, coef_df$feature)],
      0
    ),
    train_min_log = tr_min,
    train_median_log = median(tr, na.rm = TRUE),
    train_max_log = tr_max,
    lung1_min_log = min(lu, na.rm = TRUE),
    lung1_median_log = median(lu, na.rm = TRUE),
    lung1_max_log = max(lu, na.rm = TRUE),
    lung1_max_abs_z_by_train_log = max(abs(lu_z), na.rm = TRUE),
    outside_train_minmax_n_log = outside_n,
    outside_train_minmax_pct_log = outside_n / length(lu),
    stringsAsFactors = FALSE
  )
}

shift_df <- do.call(rbind, shift_rows)

shift_df <- shift_df[
  order(
    -shift_df$is_nonzero_model_feature,
    -shift_df$lung1_max_abs_z_by_train_log
  ),
]

signed_log_outlier_q1 <- quantile(rrs_df$SignedLog_Lung1_RRS, 0.25, na.rm = TRUE)
signed_log_outlier_q3 <- quantile(rrs_df$SignedLog_Lung1_RRS, 0.75, na.rm = TRUE)
signed_log_outlier_iqr <- signed_log_outlier_q3 - signed_log_outlier_q1
signed_log_lower_3iqr <- signed_log_outlier_q1 - 3 * signed_log_outlier_iqr
signed_log_upper_3iqr <- signed_log_outlier_q3 + 3 * signed_log_outlier_iqr

rrs_df$SignedLog_RRS_extreme_3IQR <- rrs_df$SignedLog_Lung1_RRS < signed_log_lower_3iqr |
  rrs_df$SignedLog_Lung1_RRS > signed_log_upper_3iqr

rrs_df$SignedLog_RRS_extreme_absz5 <- abs(rrs_df$SignedLog_Lung1_RRS_z) > 5

signed_log_extreme_df <- rrs_df[
  rrs_df$SignedLog_RRS_extreme_3IQR | rrs_df$SignedLog_RRS_extreme_absz5,
]

cat("\n===== Signed-log domain shift summary =====\n")
cat("Model features:", length(model_features), "\n")
cat("Nonzero features:", length(nonzero_features), "\n")
cat("Features with any Lung1 outside training min-max after signed-log:\n")
print(sum(shift_df$outside_train_minmax_n_log > 0))

cat("\nFeatures with >5% Lung1 outside training min-max after signed-log:\n")
print(sum(shift_df$outside_train_minmax_pct_log > 0.05))

cat("\nNonzero features with any outside training min-max after signed-log:\n")
print(sum(shift_df$outside_train_minmax_n_log[shift_df$is_nonzero_model_feature] > 0))

cat("\nMax Lung1 abs z by training after signed-log:\n")
print(max(shift_df$lung1_max_abs_z_by_train_log, na.rm = TRUE))

cat("\nSigned-log RRS extreme patients:\n")
print(signed_log_extreme_df)

cat("\nTop shifted signed-log features:\n")
print(
  head(
    shift_df[
      ,
      c(
        "feature",
        "is_nonzero_model_feature",
        "coefficient",
        "lung1_max_abs_z_by_train_log",
        "outside_train_minmax_n_log",
        "outside_train_minmax_pct_log"
      )
    ],
    30
  )
)

############################################################
# 11. Merge with clinical and survival variables
############################################################

analysis_df <- merge(
  clinical,
  rrs_df,
  by = "patient_id",
  all = FALSE
)

analysis_df$OS_time_days <- suppressWarnings(as.numeric(analysis_df$Survival.time))
analysis_df$OS_event <- suppressWarnings(as.numeric(analysis_df$deadstatus.event))
analysis_df$OS_time_months <- analysis_df$OS_time_days / 30.4375

analysis_df <- analysis_df[
  is.finite(analysis_df$OS_time_months) &
    is.finite(analysis_df$OS_event) &
    analysis_df$OS_time_months > 0 &
    analysis_df$OS_event %in% c(0, 1) &
    is.finite(analysis_df$SignedLog_Lung1_RRS),
]

analysis_df$SignedLog_RRS_group <- factor(
  analysis_df$SignedLog_RRS_group,
  levels = c("RRS_low", "RRS_high")
)

analysis_df$SignedLog_RRS_z <- as.numeric(scale(analysis_df$SignedLog_Lung1_RRS))
analysis_df$SignedLog_RRS_rank_z <- as.numeric(
  scale(rank(analysis_df$SignedLog_Lung1_RRS, ties.method = "average"))
)

cat("\n===== Signed-log RRS + OS analysis dataset =====\n")
print(dim(analysis_df))

cat("\nOS event table:\n")
print(table(analysis_df$OS_event, useNA = "ifany"))

cat("\nSigned-log RRS group table:\n")
print(table(analysis_df$SignedLog_RRS_group, useNA = "ifany"))

############################################################
# 12. KM and Cox survival analysis
############################################################

km_fit <- survfit(
  Surv(OS_time_months, OS_event) ~ SignedLog_RRS_group,
  data = analysis_df
)

logrank_main <- get_logrank(
  analysis_df,
  "SignedLog_RRS_group"
)

survival_rates <- get_survival_rates(
  km_fit,
  times_months = c(12, 24, 36, 60)
)

cat("\n===== Signed-log RRS Kaplan-Meier log-rank test =====\n")
print(logrank_main)

cat("\n===== Signed-log RRS 1/2/3/5-year OS rates =====\n")
print(survival_rates)

cox_results <- data.frame()

fit_uni_bin <- coxph(
  Surv(OS_time_months, OS_event) ~ SignedLog_RRS_group,
  data = analysis_df
)

cox_results <- rbind(
  cox_results,
  extract_cox(fit_uni_bin, "SignedLog_univariate_binary_RRS_high_vs_low")
)

fit_uni_cont <- coxph(
  Surv(OS_time_months, OS_event) ~ SignedLog_RRS_z,
  data = analysis_df
)

cox_results <- rbind(
  cox_results,
  extract_cox(fit_uni_cont, "SignedLog_univariate_continuous_RRS_z")
)

fit_uni_rank <- coxph(
  Surv(OS_time_months, OS_event) ~ SignedLog_RRS_rank_z,
  data = analysis_df
)

cox_results <- rbind(
  cox_results,
  extract_cox(fit_uni_rank, "SignedLog_univariate_rank_RRS_z")
)

analysis_df$age_numeric <- suppressWarnings(as.numeric(analysis_df$age))
analysis_df$gender_factor <- factor(analysis_df$gender)
analysis_df$Overall_stage_factor <- factor(analysis_df$Overall.Stage)
analysis_df$Histology_factor <- factor(analysis_df$Histology)

candidate_covs <- c(
  "age_numeric",
  "gender_factor",
  "Overall_stage_factor",
  "Histology_factor"
)

valid_covs <- valid_covariates(analysis_df, candidate_covs)

cat("\nValid adjusted covariates:\n")
print(valid_covs)

if (length(valid_covs) > 0) {
  
  formula_adj <- as.formula(
    paste(
      "Surv(OS_time_months, OS_event) ~ SignedLog_RRS_group +",
      paste(valid_covs, collapse = " + ")
    )
  )
  
  fit_adj <- coxph(
    formula_adj,
    data = analysis_df
  )
  
  cox_results <- rbind(
    cox_results,
    extract_cox(fit_adj, "SignedLog_adjusted_OverallStage_binary_RRS")
  )
}

cat("\n===== Signed-log RRS Cox results =====\n")
print(cox_results)

############################################################
# 13. Analysis summary
############################################################

main_rrs_row <- cox_results[
  cox_results$model == "SignedLog_univariate_binary_RRS_high_vs_low" &
    grepl("SignedLog_RRS_group", cox_results$variable),
]

cont_rrs_row <- cox_results[
  cox_results$model == "SignedLog_univariate_continuous_RRS_z" &
    grepl("SignedLog_RRS_z", cox_results$variable),
]

rank_rrs_row <- cox_results[
  cox_results$model == "SignedLog_univariate_rank_RRS_z" &
    grepl("SignedLog_RRS_rank_z", cox_results$variable),
]

adj_rrs_row <- cox_results[
  cox_results$model == "SignedLog_adjusted_OverallStage_binary_RRS" &
    grepl("SignedLog_RRS_group", cox_results$variable),
]

analysis_summary <- data.frame(
  item = c(
    "training_patients",
    "lung1_patients_with_features",
    "lung1_OS_analysis_patients",
    "OS_events",
    "signed_log_RRS_median_cutoff",
    "RRS_high_n",
    "RRS_low_n",
    "training_CV_mean_Pearson",
    "training_CV_mean_Spearman",
    "training_CV_mean_R2",
    "training_mean_prediction_Pearson",
    "training_mean_prediction_Spearman",
    "training_mean_prediction_R2",
    "final_alpha",
    "final_lambda",
    "final_nonzero_features_excluding_intercept",
    "signed_log_RRS_min",
    "signed_log_RRS_max",
    "signed_log_RRS_extreme_3IQR_or_absz5_n",
    "signedlog_features_any_outside_train_minmax",
    "signedlog_features_more_than_5pct_outside_train_minmax",
    "signedlog_nonzero_features_any_outside_train_minmax",
    "signedlog_max_Lung1_abs_z_by_training",
    "logrank_p",
    "univariate_HR_RRS_high_vs_low",
    "univariate_CI_lower",
    "univariate_CI_upper",
    "univariate_p",
    "continuous_RRS_z_HR",
    "continuous_RRS_z_p",
    "rank_RRS_z_HR",
    "rank_RRS_z_p",
    "adjusted_OverallStage_HR_RRS_high_vs_low",
    "adjusted_OverallStage_CI_lower",
    "adjusted_OverallStage_CI_upper",
    "adjusted_OverallStage_p"
  ),
  value = c(
    nrow(train_df),
    nrow(lung_feat),
    nrow(analysis_df),
    sum(analysis_df$OS_event == 1, na.rm = TRUE),
    rrs_median,
    sum(analysis_df$SignedLog_RRS_group == "RRS_high"),
    sum(analysis_df$SignedLog_RRS_group == "RRS_low"),
    cv_metrics_summary$repeat_mean[cv_metrics_summary$metric == "Pearson"],
    cv_metrics_summary$repeat_mean[cv_metrics_summary$metric == "Spearman"],
    cv_metrics_summary$repeat_mean[cv_metrics_summary$metric == "R2"],
    cv_metrics_summary$mean_prediction_metric[cv_metrics_summary$metric == "Pearson"],
    cv_metrics_summary$mean_prediction_metric[cv_metrics_summary$metric == "Spearman"],
    cv_metrics_summary$mean_prediction_metric[cv_metrics_summary$metric == "R2"],
    final_alpha,
    final_lambda,
    length(nonzero_features),
    min(rrs_df$SignedLog_Lung1_RRS, na.rm = TRUE),
    max(rrs_df$SignedLog_Lung1_RRS, na.rm = TRUE),
    nrow(signed_log_extreme_df),
    sum(shift_df$outside_train_minmax_n_log > 0),
    sum(shift_df$outside_train_minmax_pct_log > 0.05),
    sum(shift_df$outside_train_minmax_n_log[shift_df$is_nonzero_model_feature] > 0),
    max(shift_df$lung1_max_abs_z_by_train_log, na.rm = TRUE),
    logrank_main$p_value,
    ifelse(nrow(main_rrs_row) > 0, main_rrs_row$HR[1], NA),
    ifelse(nrow(main_rrs_row) > 0, main_rrs_row$CI_lower[1], NA),
    ifelse(nrow(main_rrs_row) > 0, main_rrs_row$CI_upper[1], NA),
    ifelse(nrow(main_rrs_row) > 0, main_rrs_row$p_value[1], NA),
    ifelse(nrow(cont_rrs_row) > 0, cont_rrs_row$HR[1], NA),
    ifelse(nrow(cont_rrs_row) > 0, cont_rrs_row$p_value[1], NA),
    ifelse(nrow(rank_rrs_row) > 0, rank_rrs_row$HR[1], NA),
    ifelse(nrow(rank_rrs_row) > 0, rank_rrs_row$p_value[1], NA),
    ifelse(nrow(adj_rrs_row) > 0, adj_rrs_row$HR[1], NA),
    ifelse(nrow(adj_rrs_row) > 0, adj_rrs_row$CI_lower[1], NA),
    ifelse(nrow(adj_rrs_row) > 0, adj_rrs_row$CI_upper[1], NA),
    ifelse(nrow(adj_rrs_row) > 0, adj_rrs_row$p_value[1], NA)
  ),
  stringsAsFactors = FALSE
)

cat("\n===== Signed-log analysis summary =====\n")
print(analysis_summary)

############################################################
# 14. Save outputs
############################################################

write.csv(
  feature_qc,
  file = file.path(out_dir, "signed_log_feature_qc.csv"),
  row.names = FALSE
)

write.csv(
  cv_repeat_summary,
  file = file.path(out_dir, "signed_log_repeated_cv_by_repeat.csv"),
  row.names = FALSE
)

write.csv(
  cv_alpha_summary,
  file = file.path(out_dir, "signed_log_repeated_cv_alpha_summary.csv"),
  row.names = FALSE
)

write.csv(
  cv_metrics_summary,
  file = file.path(out_dir, "signed_log_repeated_cv_metrics.csv"),
  row.names = FALSE
)

write.csv(
  train_cv_pred,
  file = file.path(out_dir, "signed_log_cross_validated_predictions.csv"),
  row.names = FALSE
)

write.csv(
  final_fit_summary,
  file = file.path(out_dir, "signed_log_final_fit_summary_by_alpha.csv"),
  row.names = FALSE
)

write.csv(
  coef_df,
  file = file.path(out_dir, "signed_log_final_model_coefficients.csv"),
  row.names = FALSE
)

write.csv(
  nonzero_coef_df,
  file = file.path(out_dir, "signed_log_model_nonzero_coefficients.csv"),
  row.names = FALSE
)

write.csv(
  rrs_df,
  file = file.path(out_dir, "lung1_signed_log_rrs_scores.csv"),
  row.names = FALSE
)

write.csv(
  shift_df,
  file = file.path(out_dir, "signed_log_domain_shift_features.csv"),
  row.names = FALSE
)

write.csv(
  signed_log_extreme_df,
  file = file.path(out_dir, "signed_log_extreme_patients.csv"),
  row.names = FALSE
)

write.csv(
  analysis_df,
  file = file.path(out_dir, "lung1_signed_log_rrs_survival_dataset.csv"),
  row.names = FALSE
)

write.csv(
  logrank_main,
  file = file.path(out_dir, "lung1_signed_log_rrs_logrank.csv"),
  row.names = FALSE
)

write.csv(
  survival_rates,
  file = file.path(out_dir, "lung1_signed_log_rrs_survival_rates.csv"),
  row.names = FALSE
)

write.csv(
  cox_results,
  file = file.path(out_dir, "lung1_signed_log_rrs_cox_results.csv"),
  row.names = FALSE
)

write.csv(
  analysis_summary,
  file = file.path(out_dir, "signed_log_analysis_summary.csv"),
  row.names = FALSE
)

saveRDS(
  analysis_df,
  file = file.path(out_dir, "lung1_signed_log_rrs_survival_dataset.rds")
)

############################################################
# 15. Excel workbook
############################################################

wb <- createWorkbook()

addWorksheet(wb, "analysis_summary")
writeData(wb, "analysis_summary", analysis_summary)

addWorksheet(wb, "training_CV_metrics")
writeData(wb, "training_CV_metrics", cv_metrics_summary)

addWorksheet(wb, "repeat_summary")
writeData(wb, "repeat_summary", cv_repeat_summary)

addWorksheet(wb, "alpha_summary")
writeData(wb, "alpha_summary", cv_alpha_summary)

addWorksheet(wb, "coefficients")
writeData(wb, "coefficients", coef_df)

addWorksheet(wb, "nonzero_coefficients")
writeData(wb, "nonzero_coefficients", nonzero_coef_df)

addWorksheet(wb, "Lung1_RRS_scores")
writeData(wb, "Lung1_RRS_scores", rrs_df)

addWorksheet(wb, "domain_shift")
writeData(wb, "domain_shift", shift_df)

addWorksheet(wb, "extreme_RRS")
writeData(wb, "extreme_RRS", signed_log_extreme_df)

addWorksheet(wb, "logrank")
writeData(wb, "logrank", logrank_main)

addWorksheet(wb, "OS_rates")
writeData(wb, "OS_rates", survival_rates)

addWorksheet(wb, "cox_results")
writeData(wb, "cox_results", cox_results)

addWorksheet(wb, "analysis_dataset")
writeData(wb, "analysis_dataset", analysis_df)

saveWorkbook(
  wb,
  file = file.path(out_dir, "signed_log_rrs_analysis.xlsx"),
  overwrite = TRUE
)

############################################################
# 16. Figures
############################################################

png(
  filename = file.path(fig_dir, "signed_log_training_CRSz_vs_SignedLog_RRS_CV.png"),
  width = 1500,
  height = 1300,
  res = 200
)

plot(
  train_cv_pred$CRS_z,
  train_cv_pred$SignedLog_RRS_CV_mean,
  xlab = "Observed CRS_z",
  ylab = "Cross-validated signed-log RRS",
  main = "Training cohort: observed CRS_z vs signed-log RRS",
  pch = 19
)

abline(lm(SignedLog_RRS_CV_mean ~ CRS_z, data = train_cv_pred), lwd = 2)

dev.off()

png(
  filename = file.path(fig_dir, "signed_log_Lung1_SignedLog_RRS_distribution.png"),
  width = 1500,
  height = 1200,
  res = 200
)

hist(
  rrs_df$SignedLog_Lung1_RRS,
  breaks = 40,
  main = "Distribution of signed-log Lung1 RRS",
  xlab = "Signed-log Lung1 RRS"
)

abline(v = rrs_median, lwd = 2, lty = 2)

dev.off()

png(
  filename = file.path(fig_dir, "signed_log_Lung1_SignedLog_RRS_KM_OS.png"),
  width = 1600,
  height = 1300,
  res = 200
)

plot(
  km_fit,
  col = c("black", "red"),
  lwd = 2,
  xlab = "Overall survival time (months)",
  ylab = "Overall survival probability",
  main = "Lung1 external validation: OS by signed-log RRS group",
  mark.time = TRUE
)

legend(
  "topright",
  legend = c(
    paste0("RRS_low, n=", sum(analysis_df$SignedLog_RRS_group == "RRS_low")),
    paste0("RRS_high, n=", sum(analysis_df$SignedLog_RRS_group == "RRS_high")),
    paste0("Log-rank p=", signif(logrank_main$p_value, 3))
  ),
  col = c("black", "red", NA),
  lwd = c(2, 2, NA),
  bty = "n"
)

dev.off()

top_shift_plot <- head(
  shift_df[
    order(-shift_df$lung1_max_abs_z_by_train_log),
  ],
  20
)

png(
  filename = file.path(fig_dir, "signed_log_top20_signedlog_feature_shift.png"),
  width = 1800,
  height = 1400,
  res = 200
)

barplot(
  top_shift_plot$lung1_max_abs_z_by_train_log,
  names.arg = top_shift_plot$feature,
  las = 2,
  cex.names = 0.55,
  ylab = "Max absolute Lung1 z-score using training mean/SD after signed-log",
  main = "Top 20 signed-log model features with largest external shift"
)

dev.off()

############################################################
# 17. Done
############################################################

cat("\n===== Signed-log RRS sensitivity analysis completed =====\n")
cat("Main output files:\n")

cat("1. Signed-log training CV summary\n")
cat("2. Final signed-log model\n")
cat("3. Signed-log Lung1 RRS summary\n")
cat("4. Signed-log domain shift summary\n")
cat("5. Signed-log RRS Kaplan-Meier log-rank test\n")
cat("6. Signed-log RRS Cox results\n")
cat("7. Signed-log analysis summary\n")
