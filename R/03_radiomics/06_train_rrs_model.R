# 06_train_rrs_model.R
#
# Train the repeated cross-validated radiomic surrogate model.
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

required_packages <- c("glmnet", "openxlsx")
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

############################################################
# 2. Paths
############################################################

merged_path <- "07_results/radiomics_features/radiomics_molecular_merged.rds"

out_dir <- "07_results/rrs_model"
fig_dir <- "08_figures_final/rrs_model"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(merged_path)) {
  stop("Cannot find merged molecular-radiomics table: ", merged_path)
}

############################################################
# 3. Load data
############################################################

dat <- readRDS(merged_path)

cat("\n===== Loaded merged table =====\n")
print(dim(dat))

cat("\nColumns preview:\n")
print(colnames(dat)[1:min(30, ncol(dat))])

############################################################
# 4. Define outcome and features
############################################################

outcome_col <- "CRS_z"

if (!outcome_col %in% colnames(dat)) {
  stop("Outcome column not found: ", outcome_col)
}

feature_cols <- grep("^original_", colnames(dat), value = TRUE)

cat("\n===== Modeling target and feature audit =====\n")
cat("Outcome:\n")
print(outcome_col)

cat("Number of original radiomics features:\n")
print(length(feature_cols))

if (length(feature_cols) == 0) {
  stop("No original_ radiomics features found.")
}

model_df <- dat[, c("patient_id", outcome_col, feature_cols), drop = FALSE]

model_df[[outcome_col]] <- as.numeric(model_df[[outcome_col]])

for (cc in feature_cols) {
  model_df[[cc]] <- suppressWarnings(as.numeric(model_df[[cc]]))
}

model_df <- model_df[complete.cases(model_df[, c(outcome_col, feature_cols)]), ]

cat("\nModeling data dim after complete-case filtering:\n")
print(dim(model_df))

cat("\nOutcome summary:\n")
print(summary(model_df[[outcome_col]]))

############################################################
# 5. Utility functions
############################################################

calc_metrics <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  
  pearson <- suppressWarnings(cor(obs, pred, method = "pearson"))
  spearman <- suppressWarnings(cor(obs, pred, method = "spearman"))
  mae <- mean(abs(obs - pred))
  rmse <- sqrt(mean((obs - pred)^2))
  
  ss_res <- sum((obs - pred)^2)
  ss_tot <- sum((obs - mean(obs))^2)
  
  r2 <- 1 - ss_res / ss_tot
  
  data.frame(
    n = length(obs),
    Pearson = pearson,
    Spearman = spearman,
    MAE = mae,
    RMSE = rmse,
    R2 = r2,
    stringsAsFactors = FALSE
  )
}

make_folds <- function(n, k = 5) {
  sample(rep(seq_len(k), length.out = n))
}

remove_high_correlation_train <- function(x_train, cutoff = 0.90) {
  if (ncol(x_train) <= 1) {
    return(colnames(x_train))
  }
  
  cmat <- suppressWarnings(cor(x_train, use = "pairwise.complete.obs"))
  cmat[is.na(cmat)] <- 0
  diag(cmat) <- 0
  
  remove <- character(0)
  cols <- colnames(x_train)
  
  repeat {
    cmat_abs <- abs(cmat)
    max_cor <- max(cmat_abs, na.rm = TRUE)
    
    if (!is.finite(max_cor) || max_cor < cutoff) break
    
    idx <- which(cmat_abs == max_cor, arr.ind = TRUE)[1, ]
    
    col1 <- rownames(cmat_abs)[idx[1]]
    col2 <- colnames(cmat_abs)[idx[2]]
    
    mean_cor_1 <- mean(abs(cmat[col1, ]), na.rm = TRUE)
    mean_cor_2 <- mean(abs(cmat[col2, ]), na.rm = TRUE)
    
    drop_col <- ifelse(mean_cor_1 >= mean_cor_2, col1, col2)
    
    remove <- c(remove, drop_col)
    
    keep <- setdiff(colnames(cmat), drop_col)
    cmat <- cmat[keep, keep, drop = FALSE]
    
    if (ncol(cmat) <= 1) break
  }
  
  setdiff(cols, remove)
}

scale_train_test <- function(x_train, x_test) {
  mu <- apply(x_train, 2, mean, na.rm = TRUE)
  sdv <- apply(x_train, 2, sd, na.rm = TRUE)
  
  sdv[is.na(sdv) | sdv == 0] <- 1
  
  x_train_scaled <- sweep(x_train, 2, mu, "-")
  x_train_scaled <- sweep(x_train_scaled, 2, sdv, "/")
  
  x_test_scaled <- sweep(x_test, 2, mu, "-")
  x_test_scaled <- sweep(x_test_scaled, 2, sdv, "/")
  
  list(
    train = x_train_scaled,
    test = x_test_scaled,
    center = mu,
    scale = sdv
  )
}

############################################################
# 6. Repeated CV modeling
############################################################

set.seed(20260724)

x_all <- as.matrix(model_df[, feature_cols, drop = FALSE])
y_all <- model_df[[outcome_col]]
patient_ids <- model_df$patient_id

n <- nrow(model_df)

n_repeats <- 20
k_folds <- 5

alpha_grid <- c(0, 0.25, 0.5, 0.75, 1)

all_pred <- data.frame()
repeat_metrics <- data.frame()
fold_model_info <- data.frame()

cat("\n===== Start repeated cross-validation =====\n")
cat("n =", n, "\n")
cat("Repeats =", n_repeats, "\n")
cat("Folds =", k_folds, "\n")
cat("Alpha grid =", paste(alpha_grid, collapse = ", "), "\n")

for (rep_i in seq_len(n_repeats)) {
  
  cat("\n==============================\n")
  cat("Repeat", rep_i, "of", n_repeats, "\n")
  
  fold_id <- make_folds(n, k = k_folds)
  
  pred_rep <- rep(NA_real_, n)
  alpha_rep <- rep(NA_real_, n)
  lambda_rep <- rep(NA_real_, n)
  n_features_rep <- rep(NA_integer_, n)
  
  for (fold in seq_len(k_folds)) {
    
    cat("  Fold", fold, "of", k_folds, "\n")
    
    test_idx <- which(fold_id == fold)
    train_idx <- setdiff(seq_len(n), test_idx)
    
    x_train_raw <- x_all[train_idx, , drop = FALSE]
    x_test_raw <- x_all[test_idx, , drop = FALSE]
    
    y_train <- y_all[train_idx]
    
    # Step 1: zero-variance filtering based on training set only
    train_sd <- apply(x_train_raw, 2, sd, na.rm = TRUE)
    keep1 <- names(train_sd)[is.finite(train_sd) & train_sd > 0]
    
    x_train_1 <- x_train_raw[, keep1, drop = FALSE]
    x_test_1 <- x_test_raw[, keep1, drop = FALSE]
    
    # Step 2: high-correlation filtering based on training set only
    keep2 <- remove_high_correlation_train(x_train_1, cutoff = 0.90)
    
    x_train_2 <- x_train_1[, keep2, drop = FALSE]
    x_test_2 <- x_test_1[, keep2, drop = FALSE]
    
    # Step 3: standardization based on training set only
    scaled <- scale_train_test(x_train_2, x_test_2)
    
    x_train <- scaled$train
    x_test <- scaled$test
    
    # Step 4: inner CV for alpha and lambda
    best_alpha <- NA
    best_lambda <- NA
    best_cvm <- Inf
    best_fit <- NULL
    
    for (aa in alpha_grid) {
      
      set.seed(100000 + rep_i * 100 + fold * 10 + round(aa * 100))
      
      cvfit <- cv.glmnet(
        x = x_train,
        y = y_train,
        family = "gaussian",
        alpha = aa,
        nfolds = 5,
        standardize = FALSE,
        type.measure = "mse"
      )
      
      min_cvm <- min(cvfit$cvm, na.rm = TRUE)
      
      if (is.finite(min_cvm) && min_cvm < best_cvm) {
        best_cvm <- min_cvm
        best_alpha <- aa
        best_lambda <- cvfit$lambda.min
        best_fit <- cvfit
      }
    }
    
    if (is.null(best_fit)) {
      stop("No glmnet model fitted.")
    }
    
    pred <- as.numeric(predict(best_fit, newx = x_test, s = "lambda.min"))
    
    pred_rep[test_idx] <- pred
    alpha_rep[test_idx] <- best_alpha
    lambda_rep[test_idx] <- best_lambda
    n_features_rep[test_idx] <- ncol(x_train)
    
    fold_model_info <- rbind(
      fold_model_info,
      data.frame(
        repeat_id = rep_i,
        fold = fold,
        n_train = length(train_idx),
        n_test = length(test_idx),
        n_features_after_filter = ncol(x_train),
        best_alpha = best_alpha,
        best_lambda = best_lambda,
        best_inner_cvm = best_cvm,
        stringsAsFactors = FALSE
      )
    )
  }
  
  one_pred <- data.frame(
    repeat_id = rep_i,
    patient_id = patient_ids,
    observed_CRS_z = y_all,
    predicted_CRS_z = pred_rep,
    alpha = alpha_rep,
    lambda = lambda_rep,
    n_features_after_filter = n_features_rep,
    stringsAsFactors = FALSE
  )
  
  all_pred <- rbind(all_pred, one_pred)
  
  met <- calc_metrics(y_all, pred_rep)
  met$repeat_id <- rep_i
  
  repeat_metrics <- rbind(repeat_metrics, met)
  
  cat("Repeat", rep_i, "metrics:\n")
  print(met)
}

############################################################
# 7. Aggregate CV prediction per patient
############################################################

patient_pred <- aggregate(
  predicted_CRS_z ~ patient_id + observed_CRS_z,
  data = all_pred,
  FUN = mean
)

colnames(patient_pred)[colnames(patient_pred) == "predicted_CRS_z"] <- "predicted_CRS_z_mean"

patient_pred_sd <- aggregate(
  predicted_CRS_z ~ patient_id + observed_CRS_z,
  data = all_pred,
  FUN = sd
)

colnames(patient_pred_sd)[colnames(patient_pred_sd) == "predicted_CRS_z"] <- "predicted_CRS_z_sd"

patient_pred <- merge(
  patient_pred,
  patient_pred_sd,
  by = c("patient_id", "observed_CRS_z"),
  all.x = TRUE
)

overall_metrics_mean_pred <- calc_metrics(
  obs = patient_pred$observed_CRS_z,
  pred = patient_pred$predicted_CRS_z_mean
)

cat("\n===== Repeated CV metrics by repeat =====\n")
print(repeat_metrics)

cat("\n===== Summary of repeated CV metrics =====\n")

metrics_summary <- data.frame(
  metric = c("Pearson", "Spearman", "MAE", "RMSE", "R2"),
  mean = c(
    mean(repeat_metrics$Pearson, na.rm = TRUE),
    mean(repeat_metrics$Spearman, na.rm = TRUE),
    mean(repeat_metrics$MAE, na.rm = TRUE),
    mean(repeat_metrics$RMSE, na.rm = TRUE),
    mean(repeat_metrics$R2, na.rm = TRUE)
  ),
  sd = c(
    sd(repeat_metrics$Pearson, na.rm = TRUE),
    sd(repeat_metrics$Spearman, na.rm = TRUE),
    sd(repeat_metrics$MAE, na.rm = TRUE),
    sd(repeat_metrics$RMSE, na.rm = TRUE),
    sd(repeat_metrics$R2, na.rm = TRUE)
  ),
  median = c(
    median(repeat_metrics$Pearson, na.rm = TRUE),
    median(repeat_metrics$Spearman, na.rm = TRUE),
    median(repeat_metrics$MAE, na.rm = TRUE),
    median(repeat_metrics$RMSE, na.rm = TRUE),
    median(repeat_metrics$R2, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

print(metrics_summary)

cat("\n===== Metrics using mean prediction across repeats =====\n")
print(overall_metrics_mean_pred)

############################################################
# 8. Fit final model on all samples for interpretation only
############################################################

cat("\n===== Fit final model on all samples =====\n")

train_sd_all <- apply(x_all, 2, sd, na.rm = TRUE)
keep1_all <- names(train_sd_all)[is.finite(train_sd_all) & train_sd_all > 0]

x_all_1 <- x_all[, keep1_all, drop = FALSE]

keep2_all <- remove_high_correlation_train(x_all_1, cutoff = 0.90)

x_all_2 <- x_all_1[, keep2_all, drop = FALSE]

scaled_all <- scale_train_test(x_all_2, x_all_2)
x_final <- scaled_all$train

best_alpha_final <- NA
best_lambda_final <- NA
best_cvm_final <- Inf
best_cvfit_final <- NULL

for (aa in alpha_grid) {
  
  set.seed(20260724 + round(aa * 1000))
  
  cvfit <- cv.glmnet(
    x = x_final,
    y = y_all,
    family = "gaussian",
    alpha = aa,
    nfolds = 5,
    standardize = FALSE,
    type.measure = "mse"
  )
  
  min_cvm <- min(cvfit$cvm, na.rm = TRUE)
  
  if (is.finite(min_cvm) && min_cvm < best_cvm_final) {
    best_cvm_final <- min_cvm
    best_alpha_final <- aa
    best_lambda_final <- cvfit$lambda.min
    best_cvfit_final <- cvfit
  }
}

final_coef <- coef(best_cvfit_final, s = "lambda.min")

coef_df <- data.frame(
  feature = rownames(final_coef),
  coefficient = as.numeric(final_coef),
  stringsAsFactors = FALSE
)

coef_df <- coef_df[coef_df$coefficient != 0, ]
coef_df <- coef_df[order(abs(coef_df$coefficient), decreasing = TRUE), ]

cat("\nFinal model alpha:\n")
print(best_alpha_final)

cat("\nFinal model lambda:\n")
print(best_lambda_final)

cat("\nNumber of final features after global filter:\n")
print(ncol(x_final))

cat("\nNumber of non-zero coefficients including intercept:\n")
print(nrow(coef_df))

cat("\nTop non-zero coefficients:\n")
print(head(coef_df, 30))

############################################################
# 9. Save outputs
############################################################

write.csv(
  all_pred,
  file = file.path(out_dir, "repeated_cv_predictions.csv"),
  row.names = FALSE
)

write.csv(
  patient_pred,
  file = file.path(out_dir, "patient_mean_predictions.csv"),
  row.names = FALSE
)

write.csv(
  repeat_metrics,
  file = file.path(out_dir, "repeated_cv_metrics_by_repeat.csv"),
  row.names = FALSE
)

write.csv(
  metrics_summary,
  file = file.path(out_dir, "repeated_cv_summary.csv"),
  row.names = FALSE
)

write.csv(
  fold_model_info,
  file = file.path(out_dir, "outer_fold_model_info.csv"),
  row.names = FALSE
)

write.csv(
  coef_df,
  file = file.path(out_dir, "final_model_coefficients.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    final_cvfit = best_cvfit_final,
    final_alpha = best_alpha_final,
    final_lambda = best_lambda_final,
    final_features = colnames(x_final),
    center = scaled_all$center,
    scale = scaled_all$scale,
    outcome_col = outcome_col,
    feature_cols_original = feature_cols,
    model_note = "Final model fitted on all samples for interpretation only; performance should be reported from repeated CV."
  ),
  file = file.path(out_dir, "rrs_elastic_net_model.rds")
)

wb <- createWorkbook()

addWorksheet(wb, "metrics_summary")
writeData(wb, "metrics_summary", metrics_summary)

addWorksheet(wb, "mean_prediction_metrics")
writeData(wb, "mean_prediction_metrics", overall_metrics_mean_pred)

addWorksheet(wb, "metrics_by_repeat")
writeData(wb, "metrics_by_repeat", repeat_metrics)

addWorksheet(wb, "patient_mean_predictions")
writeData(wb, "patient_mean_predictions", patient_pred)

addWorksheet(wb, "fold_model_info")
writeData(wb, "fold_model_info", fold_model_info)

addWorksheet(wb, "final_coefficients")
writeData(wb, "final_coefficients", coef_df)

saveWorkbook(
  wb,
  file = file.path(out_dir, "rrs_model_summary.xlsx"),
  overwrite = TRUE
)

############################################################
# 10. Figures
############################################################

png(
  filename = file.path(fig_dir, "observed_vs_predicted.png"),
  width = 1600,
  height = 1400,
  res = 200
)

plot(
  patient_pred$observed_CRS_z,
  patient_pred$predicted_CRS_z_mean,
  pch = 19,
  xlab = "Observed CRS_z",
  ylab = "Predicted CRS_z",
  main = "CT radiomics prediction of CRS_z"
)

abline(lm(predicted_CRS_z_mean ~ observed_CRS_z, data = patient_pred), lwd = 2)
abline(0, 1, lty = 2)

legend(
  "topleft",
  legend = paste0(
    "Spearman = ",
    round(overall_metrics_mean_pred$Spearman, 3),
    "\nPearson = ",
    round(overall_metrics_mean_pred$Pearson, 3)
  ),
  bty = "n"
)

dev.off()

png(
  filename = file.path(fig_dir, "repeated_cv_spearman_distribution.png"),
  width = 1600,
  height = 1200,
  res = 200
)

hist(
  repeat_metrics$Spearman,
  breaks = 10,
  main = "Repeated CV Spearman distribution",
  xlab = "Spearman correlation"
)

abline(v = mean(repeat_metrics$Spearman, na.rm = TRUE), lwd = 2)

dev.off()

############################################################
# 11. Done
############################################################

cat("\nRRS model training completed\n")
cat("Main output files:\n")
cat("1. 07_results/rrs_model/rrs_model_summary.xlsx\n")
cat("2. 07_results/rrs_model/repeated_cv_summary.csv\n")
cat("3. 07_results/rrs_model/patient_mean_predictions.csv\n")
cat("4. 07_results/rrs_model/final_model_coefficients.csv\n")
cat("5. 08_figures_final/rrs_model/observed_vs_predicted.png\n")

cat("1. Summary of repeated CV metrics\n")
cat("2. Metrics using mean prediction across repeats\n")
cat("3. Final model alpha/lambda\n")
cat("4. Number of final features and non-zero coefficients\n")
cat("5. Top non-zero coefficients\n")
