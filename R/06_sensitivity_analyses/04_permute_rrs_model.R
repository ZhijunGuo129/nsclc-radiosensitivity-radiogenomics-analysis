# 04_permute_rrs_model.R
#
# Run the 2,000-permutation strict repeated-CV RRS audit.
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
# 0. Packages
############################################################

required_packages <- c("glmnet", "future", "future.apply", "openxlsx", "ggplot2")
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

############################################################
# Configuration
############################################################


input_file <- file.path(
  project_dir,
  "07_results",
  "radiomics_features",
  "radiomics_molecular_merged.csv"
)

# Archived RRS training outputs are used to recover the exact original
# 20 x 5 outer-fold assignments. They are not used as outcomes
# or predictions in the new nested-CV performance calculation.
archived_prediction_file <- file.path(
  project_dir,
  "07_results",
  "rrs_model",
  "repeated_cv_predictions.csv"
)

archived_fold_info_file <- file.path(
  project_dir,
  "07_results",
  "rrs_model",
  "outer_fold_model_info.csv"
)

out_dir <- file.path(
  project_dir,
  "07_results",
  "rrs_permutation"
)

fig_dir <- file.path(
  project_dir,
  "08_figures_final",
  "rrs_permutation"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Final target. The script safely resumes if interrupted.
n_perm_target <- 2000L

# Reproduce the main RRS modeling design.
n_repeats <- 20L
outer_k <- 5L
inner_k <- 5L
alpha_grid <- c(0, 0.25, 0.50, 0.75, 1.00)
corr_cutoff <- 0.90

# Reproducible analysis seeds.
# Exact outer folds are recovered from the archived RRS model audit.
# The original inner-fold IDs were not archived, so this strict
# reanalysis defines and freezes new reproducible inner folds.
inner_seed_base <- 2026072500L
permutation_seed_base <- 2026080400L

# Small batches reduce loss if R/RStudio is interrupted.
physical_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (!is.finite(physical_cores)) physical_cores <- 4L
n_workers <- max(1L, min(6L, as.integer(physical_cores) - 1L))
batch_size <- max(1L, n_workers)

############################################################
# 2. Output paths
############################################################

progress_rds <- file.path(out_dir, "rrs_permutation_progress.rds")
perm_csv <- file.path(out_dir, "rrs_permutation_values.csv")
audit_xlsx <- file.path(out_dir, "rrs_full_workflow_permutation.xlsx")
key_txt <- file.path(out_dir, "rrs_permutation_key_values.txt")
session_txt <- file.path(out_dir, "rrs_permutation_session_info.txt")
hist_png <- file.path(fig_dir, "rrs_permutation_histogram.png")
hist_pdf <- file.path(fig_dir, "rrs_permutation_histogram.pdf")

############################################################
# 3. Helper functions
############################################################

make_numeric <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  suppressWarnings(as.numeric(x))
}

# Balanced regression folds, generated once and then frozen.
make_regression_folds <- function(y, k = 5L, seed = 1L) {
  set.seed(seed)
  n <- length(y)
  
  # Quantile bins approximate caret-style regression stratification.
  n_bins <- min(5L, max(2L, floor(n / k)))
  probs <- seq(0, 1, length.out = n_bins + 1L)
  cuts <- unique(as.numeric(stats::quantile(y, probs = probs, na.rm = TRUE, type = 7)))
  
  if (length(cuts) < 3L) {
    return(sample(rep(seq_len(k), length.out = n)))
  }
  
  cuts[1] <- -Inf
  cuts[length(cuts)] <- Inf
  bins <- cut(y, breaks = cuts, include.lowest = TRUE, labels = FALSE)
  
  fold_id <- integer(n)
  for (bb in sort(unique(bins))) {
    idx <- which(bins == bb)
    fold_id[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  
  fold_id
}

make_balanced_folds <- function(n, k = 5L, seed = 1L) {
  set.seed(seed)
  sample(rep(seq_len(k), length.out = n))
}

# Exact iterative absolute-correlation filter.
# For each pair above cutoff, remove the feature with the larger
# mean absolute correlation to the remaining features.
correlation_filter <- function(x, cutoff = 0.90) {
  if (ncol(x) <= 1L) return(colnames(x))
  
  cm <- suppressWarnings(abs(stats::cor(x, use = "pairwise.complete.obs")))
  cm[!is.finite(cm)] <- 0
  diag(cm) <- 0
  keep_names <- colnames(x)
  
  repeat {
    max_cor <- max(cm, na.rm = TRUE)
    if (!is.finite(max_cor) || max_cor <= cutoff || ncol(cm) <= 1L) break
    
    pair <- which(cm == max_cor, arr.ind = TRUE)[1, ]
    i <- pair[1]
    j <- pair[2]
    
    mean_i <- mean(cm[i, -i, drop = TRUE], na.rm = TRUE)
    mean_j <- mean(cm[j, -j, drop = TRUE], na.rm = TRUE)
    
    remove_idx <- if (mean_i >= mean_j) i else j
    cm <- cm[-remove_idx, -remove_idx, drop = FALSE]
    keep_names <- keep_names[-remove_idx]
  }
  
  keep_names
}

prepare_outer_split <- function(
    X_df,
    train_idx,
    test_idx,
    repeat_id,
    fold_id,
    corr_cutoff,
    inner_k,
    inner_seed_base
) {
  x_train <- X_df[train_idx, , drop = FALSE]
  x_test <- X_df[test_idx, , drop = FALSE]
  
  # Training-only missingness assessment and median imputation.
  missing_rate <- vapply(x_train, function(z) mean(!is.finite(z)), numeric(1))
  keep_missing <- is.finite(missing_rate) & missing_rate <= 0.30
  x_train <- x_train[, keep_missing, drop = FALSE]
  x_test <- x_test[, names(keep_missing)[keep_missing], drop = FALSE]
  
  medians <- vapply(x_train, function(z) {
    med <- stats::median(z[is.finite(z)], na.rm = TRUE)
    if (!is.finite(med)) 0 else med
  }, numeric(1))
  
  for (nm in names(x_train)) {
    x_train[[nm]][!is.finite(x_train[[nm]])] <- medians[[nm]]
    x_test[[nm]][!is.finite(x_test[[nm]])] <- medians[[nm]]
  }
  
  # Training-only zero-variance removal.
  sds <- vapply(x_train, stats::sd, numeric(1), na.rm = TRUE)
  keep_sd <- is.finite(sds) & sds > 0
  x_train <- x_train[, keep_sd, drop = FALSE]
  x_test <- x_test[, names(keep_sd)[keep_sd], drop = FALSE]
  
  # Training-only redundancy filtering.
  selected_features <- correlation_filter(x_train, cutoff = corr_cutoff)
  x_train <- as.matrix(x_train[, selected_features, drop = FALSE])
  x_test <- as.matrix(x_test[, selected_features, drop = FALSE])
  storage.mode(x_train) <- "double"
  storage.mode(x_test) <- "double"
  
  inner_foldid <- make_balanced_folds(
    n = length(train_idx),
    k = inner_k,
    seed = inner_seed_base + repeat_id * 100L + fold_id
  )
  
  list(
    repeat_id = repeat_id,
    fold = fold_id,
    train_idx = train_idx,
    test_idx = test_idx,
    x_train = x_train,
    x_test = x_test,
    inner_foldid = inner_foldid,
    selected_features = selected_features,
    n_features_after_filter = length(selected_features)
  )
}

fit_one_outer_split <- function(split_obj, y_current, alpha_grid) {
  y_train <- y_current[split_obj$train_idx]
  
  best_cvm <- Inf
  best_fit <- NULL
  best_alpha <- NA_real_
  best_lambda <- NA_real_
  
  for (aa in alpha_grid) {
    fit <- tryCatch(
      glmnet::cv.glmnet(
        x = split_obj$x_train,
        y = y_train,
        family = "gaussian",
        alpha = aa,
        foldid = split_obj$inner_foldid,
        standardize = TRUE,
        type.measure = "mse",
        intercept = TRUE
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit)) next
    
    idx_min <- which.min(fit$cvm)
    cvm_min <- fit$cvm[idx_min]
    lambda_min <- fit$lambda[idx_min]
    
    if (is.finite(cvm_min) && cvm_min < best_cvm) {
      best_cvm <- cvm_min
      best_fit <- fit
      best_alpha <- aa
      best_lambda <- lambda_min
    }
  }
  
  if (is.null(best_fit)) {
    pred <- rep(mean(y_train), length(split_obj$test_idx))
    nzero <- 0L
  } else {
    pred <- as.numeric(stats::predict(
      best_fit,
      newx = split_obj$x_test,
      s = best_lambda
    ))
    
    cc <- as.matrix(stats::coef(best_fit, s = best_lambda))
    nzero <- sum(cc[-1, 1] != 0)
  }
  
  list(
    pred = pred,
    best_alpha = best_alpha,
    best_lambda = best_lambda,
    best_inner_cvm = best_cvm,
    nzero = nzero
  )
}

run_full_repeated_cv <- function(y_current, split_cache, n_repeats, n_samples, alpha_grid,
                                 keep_fold_details = FALSE) {
  pred_matrix <- matrix(NA_real_, nrow = n_samples, ncol = n_repeats)
  fold_details <- vector("list", length(split_cache))
  
  for (ss in seq_along(split_cache)) {
    split_obj <- split_cache[[ss]]
    fit <- fit_one_outer_split(split_obj, y_current, alpha_grid)
    
    pred_matrix[split_obj$test_idx, split_obj$repeat_id] <- fit$pred
    
    if (keep_fold_details) {
      fold_details[[ss]] <- data.frame(
        repeat_id = split_obj$repeat_id,
        fold = split_obj$fold,
        n_train = length(split_obj$train_idx),
        n_test = length(split_obj$test_idx),
        n_features_after_filter = split_obj$n_features_after_filter,
        best_alpha = fit$best_alpha,
        best_lambda = fit$best_lambda,
        best_inner_cvm = fit$best_inner_cvm,
        nzero = fit$nzero,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (any(!is.finite(pred_matrix))) {
    stop("At least one out-of-fold prediction is missing.")
  }
  
  patient_mean_prediction <- rowMeans(pred_matrix)
  
  repeat_metrics <- do.call(rbind, lapply(seq_len(n_repeats), function(rr) {
    pp <- pred_matrix[, rr]
    data.frame(
      repeat_id = rr,
      Pearson = suppressWarnings(stats::cor(y_current, pp, method = "pearson")),
      Spearman = suppressWarnings(stats::cor(y_current, pp, method = "spearman")),
      MAE = mean(abs(y_current - pp)),
      RMSE = sqrt(mean((y_current - pp)^2)),
      R2 = suppressWarnings(stats::cor(y_current, pp, method = "pearson")^2),
      stringsAsFactors = FALSE
    )
  }))
  
  patient_mean_metrics <- data.frame(
    n = n_samples,
    Pearson = suppressWarnings(stats::cor(y_current, patient_mean_prediction, method = "pearson")),
    Spearman = suppressWarnings(stats::cor(y_current, patient_mean_prediction, method = "spearman")),
    MAE = mean(abs(y_current - patient_mean_prediction)),
    RMSE = sqrt(mean((y_current - patient_mean_prediction)^2)),
    R2 = suppressWarnings(stats::cor(y_current, patient_mean_prediction, method = "pearson")^2),
    stringsAsFactors = FALSE
  )
  
  list(
    patient_mean_prediction = patient_mean_prediction,
    pred_matrix = if (keep_fold_details) pred_matrix else NULL,
    repeat_metrics = if (keep_fold_details) repeat_metrics else NULL,
    patient_mean_metrics = patient_mean_metrics,
    fold_details = if (keep_fold_details) do.call(rbind, fold_details) else NULL
  )
}

run_one_permutation <- function(perm_id, y, split_cache, n_repeats, n_samples,
                                alpha_grid, permutation_seed_base) {
  set.seed(permutation_seed_base + perm_id)
  y_perm <- sample(y, size = length(y), replace = FALSE)
  
  fit <- run_full_repeated_cv(
    y_current = y_perm,
    split_cache = split_cache,
    n_repeats = n_repeats,
    n_samples = n_samples,
    alpha_grid = alpha_grid,
    keep_fold_details = FALSE
  )
  
  mm <- fit$patient_mean_metrics
  
  data.frame(
    perm_id = perm_id,
    Pearson = mm$Pearson,
    Spearman = mm$Spearman,
    R2 = mm$R2,
    MAE = mm$MAE,
    RMSE = mm$RMSE,
    stringsAsFactors = FALSE
  )
}

############################################################
# 4. Read and validate the exact training table
############################################################

if (!file.exists(input_file)) {
  stop("Exact training file was not found:\n", input_file)
}

input_md5 <- unname(tools::md5sum(input_file))
dat <- read.csv(input_file, check.names = FALSE)

required_cols <- c("patient_id", "CRS_z")
missing_required <- setdiff(required_cols, names(dat))
if (length(missing_required) > 0L) {
  stop("Missing required columns: ", paste(missing_required, collapse = ", "))
}

feature_cols <- grep("^original_", names(dat), value = TRUE)

if (nrow(dat) != 117L) {
  stop("Expected 117 patients, detected ", nrow(dat), ".")
}
if (length(feature_cols) != 107L) {
  stop("Expected exactly 107 original radiomic features, detected ", length(feature_cols), ".")
}
if (!("original_glcm_Id" %in% feature_cols)) {
  stop("Critical feature original_glcm_Id was not retained.")
}

patient_id <- as.character(dat$patient_id)
y <- make_numeric(dat$CRS_z)
X_df <- as.data.frame(dat[, feature_cols, drop = FALSE])
X_df[] <- lapply(X_df, make_numeric)

if (any(!is.finite(y))) stop("CRS_z contains non-finite values.")
if (anyDuplicated(patient_id)) stop("patient_id contains duplicates.")

cat("\n===== Exact Figure 4 training input =====\n")
cat("Input:", input_file, "\n")
cat("MD5:", input_md5, "\n")
cat("Patients:", nrow(dat), "\n")
cat("Original radiomic features:", length(feature_cols), "\n")
cat("original_glcm_Id retained: TRUE\n")
cat("Workers:", n_workers, "\n")

############################################################
# 5. Recover the exact archived RRS outer folds and prepare
#    training-only preprocessing
############################################################

if (!file.exists(archived_prediction_file)) {
  stop("Archived RRS prediction file was not found:\n", archived_prediction_file)
}
if (!file.exists(archived_fold_info_file)) {
  stop("Archived RRS fold-information file was not found:\n", archived_fold_info_file)
}

archived_pred <- read.csv(archived_prediction_file, check.names = FALSE)
archived_fold_info <- read.csv(archived_fold_info_file, check.names = FALSE)

required_pred_cols <- c(
  "repeat_id", "patient_id", "alpha", "lambda",
  "n_features_after_filter"
)
required_fold_cols <- c(
  "repeat_id", "fold", "best_alpha", "best_lambda",
  "n_features_after_filter"
)

if (length(setdiff(required_pred_cols, names(archived_pred))) > 0L) {
  stop("Archived prediction file is missing required columns.")
}
if (length(setdiff(required_fold_cols, names(archived_fold_info))) > 0L) {
  stop("Archived fold-info file is missing required columns.")
}
if (nrow(archived_pred) != nrow(dat) * n_repeats) {
  stop("Archived prediction file should contain 117 x 20 = 2340 rows.")
}
if (nrow(archived_fold_info) != n_repeats * outer_k) {
  stop("Archived fold-info file should contain 20 x 5 = 100 rows.")
}

recover_fold <- function(repeat_id, alpha, lambda, n_features, fold_info) {
  cand <- fold_info[
    fold_info$repeat_id == repeat_id &
      abs(fold_info$best_alpha - alpha) < 1e-10 &
      abs(fold_info$best_lambda - lambda) < 1e-8 &
      fold_info$n_features_after_filter == n_features,
    ,
    drop = FALSE
  ]
  
  if (nrow(cand) != 1L) {
    stop(
      "Could not uniquely recover archived outer fold for repeat ", repeat_id,
      ", alpha ", alpha, ", lambda ", lambda,
      ", n_features ", n_features, "."
    )
  }
  
  as.integer(cand$fold[1])
}

archived_pred$recovered_fold <- mapply(
  recover_fold,
  repeat_id = archived_pred$repeat_id,
  alpha = archived_pred$alpha,
  lambda = archived_pred$lambda,
  n_features = archived_pred$n_features_after_filter,
  MoreArgs = list(fold_info = archived_fold_info)
)

outer_fold_matrix <- matrix(
  NA_integer_,
  nrow = nrow(dat),
  ncol = n_repeats,
  dimnames = list(patient_id, paste0("repeat_", seq_len(n_repeats)))
)

for (rr in seq_len(n_repeats)) {
  one_repeat <- archived_pred[archived_pred$repeat_id == rr, , drop = FALSE]
  if (nrow(one_repeat) != nrow(dat)) {
    stop("Archived repeat ", rr, " does not contain 117 patients.")
  }
  if (anyDuplicated(one_repeat$patient_id)) {
    stop("Archived repeat ", rr, " contains duplicate patient IDs.")
  }
  
  match_idx <- match(patient_id, one_repeat$patient_id)
  if (anyNA(match_idx)) {
    stop("Patient IDs do not match between the merged radiomics input and archived RRS predictions.")
  }
  
  outer_fold_matrix[, rr] <- one_repeat$recovered_fold[match_idx]
}

if (anyNA(outer_fold_matrix)) {
  stop("Failed to recover all archived outer-fold assignments.")
}

split_cache <- list()
feature_audit <- list()
cache_counter <- 0L

for (rr in seq_len(n_repeats)) {
  for (ff in seq_len(outer_k)) {
    test_idx <- which(outer_fold_matrix[, rr] == ff)
    train_idx <- setdiff(seq_len(nrow(dat)), test_idx)
    
    cache_counter <- cache_counter + 1L
    split_obj <- prepare_outer_split(
      X_df = X_df,
      train_idx = train_idx,
      test_idx = test_idx,
      repeat_id = rr,
      fold_id = ff,
      corr_cutoff = corr_cutoff,
      inner_k = inner_k,
      inner_seed_base = inner_seed_base
    )
    
    archived_n_features <- archived_fold_info$n_features_after_filter[
      archived_fold_info$repeat_id == rr &
        archived_fold_info$fold == ff
    ]
    
    if (length(archived_n_features) != 1L) {
      stop("Archived feature count is not unique for repeat ", rr, ", fold ", ff, ".")
    }
    
    split_cache[[cache_counter]] <- split_obj
    feature_audit[[cache_counter]] <- data.frame(
      repeat_id = rr,
      fold = ff,
      n_train = length(train_idx),
      n_test = length(test_idx),
      n_features_after_filter = split_obj$n_features_after_filter,
      archived_n_features_after_filter = archived_n_features,
      filter_count_matches_archived =
        split_obj$n_features_after_filter == archived_n_features,
      selected_features = paste(split_obj$selected_features, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
}

feature_audit_df <- do.call(rbind, feature_audit)

cat("\nExact archived outer folds recovered from the archived RRS outputs.\n")
cat("Training-only feature count after |r| >", corr_cutoff, "filter:\n")
print(summary(feature_audit_df$n_features_after_filter))

if (!all(feature_audit_df$filter_count_matches_archived)) {
  bad <- feature_audit_df[!feature_audit_df$filter_count_matches_archived, ]
  print(bad)
  stop(
    "The reconstructed training-only correlation filter does not reproduce ",
    "the archived RRS feature counts."
  )
}

cat(
  "Feature-filter audit: 100/100 outer folds exactly reproduce ",
  "the archived RRS feature counts.\n"
)

############################################################
# 6. Strict observed repeated-CV performance
############################################################

cat("\n===== Running observed nested-CV repeated CV =====\n")
observed_start <- Sys.time()

observed_fit <- run_full_repeated_cv(
  y_current = y,
  split_cache = split_cache,
  n_repeats = n_repeats,
  n_samples = nrow(dat),
  alpha_grid = alpha_grid,
  keep_fold_details = TRUE
)

observed_elapsed <- as.numeric(difftime(Sys.time(), observed_start, units = "mins"))
observed_metrics <- observed_fit$patient_mean_metrics

cat("\nStrict observed patient-mean performance:\n")
print(observed_metrics)
cat("Observed run time (minutes):", sprintf("%.2f", observed_elapsed), "\n")
patient_prediction_df <- data.frame(
  patient_id = patient_id,
  observed_CRS_z = y,
  RRS_CV_mean = observed_fit$patient_mean_prediction,
  RRS_CV_sd = apply(observed_fit$pred_matrix, 1, stats::sd),
  stringsAsFactors = FALSE
)

############################################################
# 7. Resume-safe strict full model-building permutations
############################################################

archived_prediction_md5 <- unname(tools::md5sum(archived_prediction_file))
archived_fold_info_md5 <- unname(tools::md5sum(archived_fold_info_file))

analysis_signature <- paste(
  "input_md5", input_md5,
  "archived_prediction_md5", archived_prediction_md5,
  "archived_fold_info_md5", archived_fold_info_md5,
  "n", nrow(dat),
  "p", length(feature_cols),
  "repeats", n_repeats,
  "outer_k", outer_k,
  "inner_k", inner_k,
  "alpha", paste(alpha_grid, collapse = ","),
  "corr", corr_cutoff,
  "inner_seed", inner_seed_base,
  "perm_seed", permutation_seed_base,
  sep = "|"
)

perm_results <- data.frame(
  perm_id = integer(0),
  Pearson = numeric(0),
  Spearman = numeric(0),
  R2 = numeric(0),
  MAE = numeric(0),
  RMSE = numeric(0),
  stringsAsFactors = FALSE
)

if (file.exists(progress_rds)) {
  old_progress <- readRDS(progress_rds)
  if (!identical(old_progress$analysis_signature, analysis_signature)) {
    stop(
      "Existing progress file was generated under a different analysis specification.\n",
      "Move or rename this file before restarting:\n", progress_rds
    )
  }
  perm_results <- old_progress$perm_results
  cat("\nLoaded existing strict progress:", nrow(perm_results), "permutations.\n")
}

completed_ids <- unique(perm_results$perm_id)
pending_ids <- setdiff(seq_len(n_perm_target), completed_ids)

if (length(pending_ids) > 0L) {
  future::plan(future::multisession, workers = n_workers)
  options(future.globals.maxSize = 4 * 1024^3)
  
  permutation_start <- Sys.time()
  
  while (length(pending_ids) > 0L) {
    current_ids <- head(pending_ids, batch_size)
    batch_start <- Sys.time()
    
    cat(
      "\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      "| strict permutations", min(current_ids), "to", max(current_ids),
      "| completed", nrow(perm_results), "/", n_perm_target, "\n"
    )
    
    batch_list <- future.apply::future_lapply(
      current_ids,
      FUN = run_one_permutation,
      y = y,
      split_cache = split_cache,
      n_repeats = n_repeats,
      n_samples = nrow(dat),
      alpha_grid = alpha_grid,
      permutation_seed_base = permutation_seed_base,
      future.seed = TRUE,
      future.scheduling = 1
    )
    
    batch_df <- do.call(rbind, batch_list)
    perm_results <- rbind(perm_results, batch_df)
    perm_results <- perm_results[order(perm_results$perm_id), , drop = FALSE]
    rownames(perm_results) <- NULL
    
    saveRDS(
      list(
        analysis_signature = analysis_signature,
        perm_results = perm_results,
        updated_at = as.character(Sys.time())
      ),
      progress_rds
    )
    write.csv(perm_results, perm_csv, row.names = FALSE)
    
    batch_minutes <- as.numeric(difftime(Sys.time(), batch_start, units = "mins"))
    completed_now <- nrow(perm_results)
    elapsed_hours <- as.numeric(difftime(Sys.time(), permutation_start, units = "hours"))
    rate_per_hour <- if (elapsed_hours > 0) length(setdiff(seq_len(completed_now), completed_ids)) / elapsed_hours else NA_real_
    remaining_hours <- if (is.finite(rate_per_hour) && rate_per_hour > 0) {
      (n_perm_target - completed_now) / rate_per_hour
    } else {
      NA_real_
    }
    
    cat(
      "Batch saved | batch minutes =", sprintf("%.2f", batch_minutes),
      "| total completed =", completed_now,
      "| estimated remaining hours =", sprintf("%.1f", remaining_hours), "\n"
    )
    
    completed_ids <- unique(perm_results$perm_id)
    pending_ids <- setdiff(seq_len(n_perm_target), completed_ids)
  }
  
  future::plan(future::sequential)
} else {
  cat("\nAll", n_perm_target, "strict permutations were already completed.\n")
}

############################################################
# 8. Empirical P values based only on the observed nested-CV result
############################################################

perm_results <- perm_results[perm_results$perm_id <= n_perm_target, , drop = FALSE]
perm_results <- perm_results[order(perm_results$perm_id), , drop = FALSE]

if (nrow(perm_results) != n_perm_target) {
  stop("Permutation result count is not equal to n_perm_target.")
}

observed_rho <- observed_metrics$Spearman
observed_r <- observed_metrics$Pearson

# Upper-tail test is the prespecified primary test because better model
# performance corresponds to a larger positive cross-validated correlation.
empirical_p_spearman_greater <- (
  sum(perm_results$Spearman >= observed_rho, na.rm = TRUE) + 1
) / (n_perm_target + 1)

# The two-sided sensitivity P value doubles the smaller empirical tail.
# This retains the permutation distribution's observed asymmetry and avoids
# the inappropriate absolute-correlation calculation.
empirical_p_spearman_less <- (
  sum(perm_results$Spearman <= observed_rho, na.rm = TRUE) + 1
) / (n_perm_target + 1)

empirical_p_spearman_two_sided <- min(
  1,
  2 * min(empirical_p_spearman_greater, empirical_p_spearman_less)
)

empirical_p_pearson_greater <- (
  sum(perm_results$Pearson >= observed_r, na.rm = TRUE) + 1
) / (n_perm_target + 1)

summary_tbl <- data.frame(
  metric = c(
    "observed_patient_mean_Spearman",
    "observed_patient_mean_Pearson",
    "observed_patient_mean_R2",
    "observed_patient_mean_MAE",
    "observed_patient_mean_RMSE",
    "empirical_p_Spearman_greater_primary",
    "empirical_p_Spearman_less",
    "empirical_p_Spearman_two_sided_doubled_tail",
    "empirical_p_Pearson_greater_secondary",
    "permutations_completed",
    "patients",
    "original_radiomic_features",
    "repeated_outer_CV",
    "inner_CV",
    "alpha_grid",
    "correlation_filter_cutoff"
  ),
  value = c(
    observed_metrics$Spearman,
    observed_metrics$Pearson,
    observed_metrics$R2,
    observed_metrics$MAE,
    observed_metrics$RMSE,
    empirical_p_spearman_greater,
    empirical_p_spearman_less,
    empirical_p_spearman_two_sided,
    empirical_p_pearson_greater,
    n_perm_target,
    nrow(dat),
    length(feature_cols),
    paste0(n_repeats, " repeats x ", outer_k, " folds"),
    paste0(inner_k, " folds"),
    paste(alpha_grid, collapse = ", "),
    corr_cutoff
  ),
  stringsAsFactors = FALSE
)

############################################################
# 9. Publication-style permutation histogram
############################################################

p_hist <- ggplot2::ggplot(perm_results, ggplot2::aes(x = Spearman)) +
  ggplot2::geom_histogram(
    bins = 36,
    fill = "#A7CBDE",
    color = "white",
    linewidth = 0.25
  ) +
  ggplot2::geom_vline(
    xintercept = observed_rho,
    color = "#FF7F0E",
    linewidth = 1.05
  ) +
  ggplot2::annotate(
    "label",
    x = -Inf,
    y = Inf,
    hjust = -0.05,
    vjust = 1.15,
    label = paste0(
      "Observed rho = ", sprintf("%.3f", observed_rho), "\n",
      "Empirical P = ", sprintf("%.4f", empirical_p_spearman_greater), "\n",
      n_perm_target, " full model-building permutations"
    ),
    size = 3.0,
    lineheight = 0.92,
    fill = "white",
    color = "black"
  ) +
  ggplot2::theme_classic(base_family = "Arial", base_size = 10) +
  ggplot2::labs(
    x = "Permuted cross-validated Spearman rho",
    y = "Count"
  ) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(face = "bold", color = "black"),
    axis.text = ggplot2::element_text(color = "black"),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA)
  )

ggplot2::ggsave(hist_png, p_hist, width = 5.4, height = 4.3, dpi = 600, bg = "white")
ggplot2::ggsave(hist_pdf, p_hist, width = 5.4, height = 4.3, bg = "white")

############################################################
# 10. Save comprehensive audit workbook
############################################################

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "summary")
openxlsx::writeData(wb, "summary", summary_tbl)

openxlsx::addWorksheet(wb, "observed_patient_predictions")
openxlsx::writeData(wb, "observed_patient_predictions", patient_prediction_df)

openxlsx::addWorksheet(wb, "observed_repeat_metrics")
openxlsx::writeData(wb, "observed_repeat_metrics", observed_fit$repeat_metrics)

openxlsx::addWorksheet(wb, "observed_fold_models")
openxlsx::writeData(wb, "observed_fold_models", observed_fit$fold_details)

openxlsx::addWorksheet(wb, "fold_feature_audit")
openxlsx::writeData(wb, "fold_feature_audit", feature_audit_df)

openxlsx::addWorksheet(wb, "permutation_values")
openxlsx::writeData(wb, "permutation_values", perm_results)

openxlsx::addWorksheet(wb, "input_features_107")
openxlsx::writeData(
  wb,
  "input_features_107",
  data.frame(feature = feature_cols, stringsAsFactors = FALSE)
)

openxlsx::addWorksheet(wb, "analysis_specification")
openxlsx::writeData(
  wb,
  "analysis_specification",
  data.frame(
    item = c(
      "Input file",
      "Input MD5",
      "Outcome",
      "Primary statistic",
      "Outer-fold source",
      "Outer CV",
      "Inner CV",
      "Alpha candidates",
      "Lambda rule",
      "Missing-value handling",
      "Zero-variance handling",
      "Redundancy filter",
      "Feature-filter scope",
      "Permutation design",
      "Empirical upper-tail P formula",
      "Empirical two-sided P formula",
      "Primary tail",
      "Parallel workers",
      "Analysis signature"
    ),
    value = c(
      input_file,
      input_md5,
      "CRS_z",
      "Spearman correlation of CRS_z with patient-level mean OOF prediction across 20 repeats",
      "Exact 20 x 5 outer-fold assignments recovered from archived RRS predictions and fold audit",
      paste0(n_repeats, " repeated ", outer_k, "-fold CV"),
      paste0(inner_k, "-fold CV within each outer training split; newly defined and frozen because original inner-fold IDs were not archived"),
      paste(alpha_grid, collapse = ", "),
      "Minimum inner-CV MSE (lambda.min equivalent)",
      "Training-fold median; applied unchanged to outer test fold",
      "Removed using outer training fold only",
      paste0("Iterative absolute-correlation filter, cutoff = ", corr_cutoff),
      "Outer training fold only",
      "Permute CRS_z; rerun all outcome-dependent alpha/lambda tuning and model fitting",
      "(number of permutation statistics at least as large as observed + 1) / (B + 1)",
      "min[1, 2 x min(upper-tail P, lower-tail P)]",
      "Upper-tail (larger positive cross-validated correlation is better); doubled-tail two-sided sensitivity also reported",
      n_workers,
      analysis_signature
    ),
    stringsAsFactors = FALSE
  )
)

openxlsx::addWorksheet(wb, "method_note")
openxlsx::writeData(
  wb,
  "method_note",
  data.frame(
    note = c(
      "This audit implements the complete repeated nested-cross-validation permutation workflow.",
      "The exact merged molecular-radiomics table is read directly; no automatic file guessing is used.",
      "All 107 original PyRadiomics features are retained at input, including original_glcm_Id.",
      "The exact original 20 x 5 outer-fold assignments are reconstructed from the archived RRS patient predictions and fold-model audit.",
      "The training-only correlation filter reproduces the archived RRS feature count in all 100 outer folds.",
      "The original inner-fold IDs were not archived; therefore this nested-CV reanalysis defines new reproducible inner folds and uses them identically for the observed and all permuted outcomes.",
      "Preprocessing parameters are estimated from each outer training fold only and then applied to its held-out test fold.",
      "Alpha and lambda are selected inside each outer training fold by inner cross-validation.",
      "The observed statistic and all permutation statistics are generated by the same repeated nested-CV function.",
      "Outcome-independent feature preprocessing is cached for computational efficiency because X and frozen folds do not change across outcome permutations; this is mathematically identical to recomputing it for every permutation.",
      "",
      "The prespecified primary test is the positive upper-tail empirical P value. The two-sided sensitivity P value doubles the smaller empirical tail and is capped at one."
    ),
    stringsAsFactors = FALSE
  )
)

for (sheet in names(wb)) {
  openxlsx::setColWidths(wb, sheet, cols = 1:80, widths = "auto")
}

openxlsx::saveWorkbook(wb, audit_xlsx, overwrite = TRUE)

############################################################
# 11. Save simple key-values file and session information
############################################################

key_lines <- c(
  "===== RRS FULL-WORKFLOW PERMUTATION RESULTS =====",
  paste0("Observed Spearman: ", format(observed_rho, digits = 10)),
  paste0("Observed Pearson: ", format(observed_r, digits = 10)),
  paste0("Empirical P Spearman upper tail (primary): ", format(empirical_p_spearman_greater, digits = 10)),
  paste0("Empirical P Spearman lower tail: ", format(empirical_p_spearman_less, digits = 10)),
  paste0("Empirical P Spearman doubled-tail two-sided: ", format(empirical_p_spearman_two_sided, digits = 10)),
  paste0("Empirical P Pearson greater: ", format(empirical_p_pearson_greater, digits = 10)),
  paste0("Permutations: ", n_perm_target),
  paste0("Patients: ", nrow(dat)),
  paste0("Input features: ", length(feature_cols)),
  paste0("Outer CV: ", n_repeats, " repeats x ", outer_k, " folds"),
  paste0("Inner CV: ", inner_k, " folds"),
  paste0("Alpha grid: ", paste(alpha_grid, collapse = ", ")),
  paste0("Feature count after training-only filtering: ",
         min(feature_audit_df$n_features_after_filter), " to ",
         max(feature_audit_df$n_features_after_filter)),
  paste0("Archived outer-fold feature counts matched: ",
         sum(feature_audit_df$filter_count_matches_archived), "/",
         nrow(feature_audit_df)),
  paste0("Audit workbook: ", audit_xlsx),
  paste0("Histogram: ", hist_png)
)

writeLines(key_lines, key_txt)
writeLines(capture.output(sessionInfo()), session_txt)

cat("\n", paste(key_lines, collapse = "\n"), "\n", sep = "")
cat("\n===== RRS full-workflow permutation analysis completed =====\n")
