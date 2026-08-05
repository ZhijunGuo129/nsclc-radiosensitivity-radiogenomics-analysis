# 03_permute_crs_model.R
#
# Estimate the empirical significance of the CRS using 2,000 complete
# model-building permutations under the same repeated nested-CV design used
# for the observed analysis.

options(stringsAsFactors = FALSE)
options(timeout = 1800)

required_env_path <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    stop("Environment variable ", name, " is not set.")
  }
  normalizePath(value, winslash = "/", mustWork = TRUE)
}

required_packages <- c(
  "glmnet",
  "future",
  "future.apply",
  "ggplot2",
  "openxlsx",
  "parallel"
)
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

project_dir <- required_env_path("RADIOGENOMICS_PROJECT_DIR")
setwd(project_dir)

expression_file <- file.path(
  project_dir,
  "01_raw_data",
  "RadioGx",
  "RadioGx_rnaseq_model_gene_by_sample.rds"
)
outcome_file <- file.path(
  project_dir,
  "05_molecular_scores",
  "RadioGx_sensitivity_model_outcomes.rds"
)
observed_results_file <- file.path(
  project_dir,
  "07_results",
  "crs_model",
  "crs_model_results.xlsx"
)
out_dir <- file.path(project_dir, "07_results", "crs_permutation")
figure_dir <- file.path(project_dir, "08_figures_final", "crs_permutation")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

for (path in c(expression_file, outcome_file, observed_results_file)) {
  if (!file.exists(path)) {
    stop("Required input was not found: ", path)
  }
}

n_permutations <- 2000L
outer_repeats <- 3L
outer_folds <- 5L
inner_folds <- 5L
top_variable_genes <- 2000L
alpha_grid <- c(0.5, 1.0)
outer_seed <- 20260722L
inner_seed <- 2026072200L
permutation_seed <- 2026080500L

physical_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (!is.finite(physical_cores)) physical_cores <- 4L
n_workers <- max(1L, min(6L, as.integer(physical_cores) - 1L))
batch_size <- max(1L, n_workers)

progress_file <- file.path(out_dir, "crs_permutation_progress.rds")
permutation_csv <- file.path(out_dir, "crs_permutation_values.csv")
results_workbook <- file.path(out_dir, "crs_full_workflow_permutation.xlsx")
key_values_file <- file.path(out_dir, "crs_permutation_key_values.txt")
session_file <- file.path(out_dir, "crs_permutation_session_info.txt")
histogram_png <- file.path(figure_dir, "crs_permutation_histogram.png")
histogram_pdf <- file.path(figure_dir, "crs_permutation_histogram.pdf")

clean_expression_matrix <- function(expression_matrix) {
  ensembl_ids <- sub("\\..*$", "", rownames(expression_matrix))
  matrix_numeric <- matrix(
    as.numeric(expression_matrix),
    nrow = nrow(expression_matrix),
    ncol = ncol(expression_matrix),
    dimnames = list(ensembl_ids, colnames(expression_matrix))
  )

  keep <- !is.na(rownames(matrix_numeric)) & nzchar(rownames(matrix_numeric))
  matrix_numeric <- matrix_numeric[keep, , drop = FALSE]

  if (anyDuplicated(rownames(matrix_numeric))) {
    groups <- rownames(matrix_numeric)
    summed <- rowsum(matrix_numeric, group = groups, reorder = FALSE)
    counts <- as.numeric(table(factor(groups, levels = rownames(summed))))
    matrix_numeric <- summed / counts
  }

  finite_values <- as.numeric(matrix_numeric[is.finite(matrix_numeric)])
  if (length(finite_values) == 0L) {
    stop("The expression matrix contains no finite values.")
  }
  if (unname(stats::quantile(finite_values, 0.99, na.rm = TRUE)) > 50) {
    matrix_numeric <- log2(matrix_numeric + 1)
  }
  matrix_numeric
}

make_balanced_folds <- function(n, k, seed) {
  set.seed(seed)
  sample(rep(seq_len(k), length.out = n))
}

select_top_variable_genes <- function(x, top_n) {
  variances <- apply(x, 2, stats::var, na.rm = TRUE)
  variances[!is.finite(variances)] <- 0
  variances <- variances[variances > 0]
  top_n <- min(as.integer(top_n), length(variances))
  names(sort(variances, decreasing = TRUE))[seq_len(top_n)]
}

prepare_analysis_data <- function() {
  expression_raw <- readRDS(expression_file)
  outcomes <- readRDS(outcome_file)
  expression_clean <- clean_expression_matrix(expression_raw)

  if (!"sampleid" %in% names(outcomes)) {
    stop("The radiation-response table does not contain sampleid.")
  }
  if (!"AUC_recomputed" %in% names(outcomes)) {
    stop("The radiation-response table does not contain AUC_recomputed.")
  }

  common_samples <- intersect(colnames(expression_clean), outcomes$sampleid)
  x <- t(expression_clean[, common_samples, drop = FALSE])
  outcome_rows <- outcomes[match(common_samples, outcomes$sampleid), , drop = FALSE]
  y <- suppressWarnings(as.numeric(outcome_rows$AUC_recomputed))

  keep <- is.finite(y) & complete.cases(x)
  x <- x[keep, , drop = FALSE]
  y <- y[keep]
  sample_ids <- rownames(x)

  if (nrow(x) != 504L) {
    warning("Expected 504 modeling samples; found ", nrow(x), ".")
  }
  if (ncol(x) < top_variable_genes) {
    stop("Fewer than 2,000 expression features are available for preselection.")
  }

  list(x = x, y = y, sample_ids = sample_ids)
}

build_split_cache <- function(x) {
  cache <- vector("list", outer_repeats * outer_folds)
  record_id <- 0L

  for (repeat_id in seq_len(outer_repeats)) {
    outer_fold_id <- make_balanced_folds(
      nrow(x),
      outer_folds,
      outer_seed + repeat_id
    )

    for (fold_id in seq_len(outer_folds)) {
      record_id <- record_id + 1L
      test_index <- which(outer_fold_id == fold_id)
      train_index <- setdiff(seq_len(nrow(x)), test_index)
      selected_features <- select_top_variable_genes(
        x[train_index, , drop = FALSE],
        top_n = top_variable_genes
      )
      inner_fold_id <- make_balanced_folds(
        length(train_index),
        inner_folds,
        inner_seed + repeat_id * 100L + fold_id
      )

      cache[[record_id]] <- list(
        repeat_id = repeat_id,
        fold_id = fold_id,
        train_index = train_index,
        test_index = test_index,
        selected_features = selected_features,
        inner_fold_id = inner_fold_id
      )
    }
  }
  cache
}

fit_one_split <- function(split, x, y_current) {
  x_train <- as.matrix(
    x[split$train_index, split$selected_features, drop = FALSE]
  )
  x_test <- as.matrix(
    x[split$test_index, split$selected_features, drop = FALSE]
  )
  y_train <- y_current[split$train_index]

  best_fit <- NULL
  best_alpha <- NA_real_
  best_cvm <- Inf

  for (alpha_value in alpha_grid) {
    fit <- glmnet::cv.glmnet(
      x = x_train,
      y = y_train,
      family = "gaussian",
      alpha = alpha_value,
      foldid = split$inner_fold_id,
      standardize = TRUE,
      type.measure = "mse"
    )
    current_cvm <- min(fit$cvm, na.rm = TRUE)
    if (is.finite(current_cvm) && current_cvm < best_cvm) {
      best_fit <- fit
      best_alpha <- alpha_value
      best_cvm <- current_cvm
    }
  }

  if (is.null(best_fit)) {
    stop(
      "Elastic-net tuning failed in repeat ", split$repeat_id,
      ", fold ", split$fold_id, "."
    )
  }

  coefficient_matrix <- as.matrix(
    stats::coef(best_fit, s = "lambda.min")
  )

  list(
    predictions = as.numeric(
      stats::predict(best_fit, newx = x_test, s = "lambda.min")
    ),
    alpha = best_alpha,
    lambda = best_fit$lambda.min,
    nonzero_coefficients = sum(coefficient_matrix[-1, 1] != 0)
  )
}

run_repeated_nested_cv <- function(x, y_current, split_cache, keep_details = FALSE) {
  prediction_matrix <- matrix(
    NA_real_,
    nrow = nrow(x),
    ncol = outer_repeats
  )
  fold_rows <- vector("list", length(split_cache))

  for (index in seq_along(split_cache)) {
    split <- split_cache[[index]]
    fit <- fit_one_split(split, x, y_current)
    prediction_matrix[split$test_index, split$repeat_id] <- fit$predictions

    if (keep_details) {
      fold_rows[[index]] <- data.frame(
        repeat_id = split$repeat_id,
        fold_id = split$fold_id,
        n_train = length(split$train_index),
        n_test = length(split$test_index),
        preselected_features = length(split$selected_features),
        alpha = fit$alpha,
        lambda = fit$lambda,
        nonzero_coefficients = fit$nonzero_coefficients,
        stringsAsFactors = FALSE
      )
    }
  }

  if (any(!is.finite(prediction_matrix))) {
    stop("At least one out-of-fold prediction is missing.")
  }

  mean_prediction <- rowMeans(prediction_matrix)
  metrics <- data.frame(
    n = length(y_current),
    Pearson = suppressWarnings(stats::cor(y_current, mean_prediction, method = "pearson")),
    Spearman = suppressWarnings(stats::cor(y_current, mean_prediction, method = "spearman")),
    R2 = suppressWarnings(stats::cor(y_current, mean_prediction, method = "pearson")^2),
    MAE = mean(abs(y_current - mean_prediction)),
    RMSE = sqrt(mean((y_current - mean_prediction)^2)),
    stringsAsFactors = FALSE
  )

  list(
    mean_prediction = mean_prediction,
    prediction_matrix = if (keep_details) prediction_matrix else NULL,
    fold_summary = if (keep_details) do.call(rbind, fold_rows) else NULL,
    metrics = metrics
  )
}

run_one_permutation <- function(permutation_id, x, y, split_cache) {
  set.seed(permutation_seed + permutation_id)
  y_permuted <- sample(y, size = length(y), replace = FALSE)
  fit <- run_repeated_nested_cv(
    x = x,
    y_current = y_permuted,
    split_cache = split_cache,
    keep_details = FALSE
  )
  data.frame(
    permutation_id = permutation_id,
    Pearson = fit$metrics$Pearson,
    Spearman = fit$metrics$Spearman,
    R2 = fit$metrics$R2,
    MAE = fit$metrics$MAE,
    RMSE = fit$metrics$RMSE,
    stringsAsFactors = FALSE
  )
}

analysis <- prepare_analysis_data()
x <- analysis$x
y <- analysis$y
sample_ids <- analysis$sample_ids
split_cache <- build_split_cache(x)

observed_fit <- run_repeated_nested_cv(
  x = x,
  y_current = y,
  split_cache = split_cache,
  keep_details = TRUE
)
observed_metrics <- observed_fit$metrics

archived_cv <- openxlsx::read.xlsx(
  observed_results_file,
  sheet = "observed_cv",
  check.names = FALSE
)
required_archived_columns <- c("sampleid", "y_true", "y_pred")
missing_archived_columns <- setdiff(required_archived_columns, names(archived_cv))
if (length(missing_archived_columns) > 0L) {
  stop(
    "The archived observed-CV table is missing: ",
    paste(missing_archived_columns, collapse = ", ")
  )
}
archived_patient <- aggregate(
  cbind(y_true, y_pred) ~ sampleid,
  data = archived_cv,
  FUN = mean
)
archived_metrics <- data.frame(
  n = nrow(archived_patient),
  Pearson = stats::cor(archived_patient$y_true, archived_patient$y_pred, method = "pearson"),
  Spearman = stats::cor(archived_patient$y_true, archived_patient$y_pred, method = "spearman"),
  stringsAsFactors = FALSE
)

if (abs(observed_metrics$Spearman - archived_metrics$Spearman) > 0.002 ||
    abs(observed_metrics$Pearson - archived_metrics$Pearson) > 0.002) {
  stop(
    "The recomputed observed performance does not reproduce the archived ",
    "CRS cross-validation result within a tolerance of 0.002."
  )
}

analysis_signature <- paste(
  tools::md5sum(expression_file),
  tools::md5sum(outcome_file),
  n_permutations,
  outer_repeats,
  outer_folds,
  inner_folds,
  top_variable_genes,
  paste(alpha_grid, collapse = ","),
  outer_seed,
  inner_seed,
  permutation_seed,
  sep = "|"
)

permutation_results <- data.frame()
if (file.exists(progress_file)) {
  progress <- tryCatch(readRDS(progress_file), error = function(e) NULL)
  if (!is.null(progress) && identical(progress$analysis_signature, analysis_signature)) {
    permutation_results <- progress$permutation_results
  }
}

completed_ids <- if (nrow(permutation_results) > 0L) {
  unique(permutation_results$permutation_id)
} else {
  integer(0)
}
pending_ids <- setdiff(seq_len(n_permutations), completed_ids)

if (length(pending_ids) > 0L) {
  future::plan(future::multisession, workers = n_workers)
  on.exit(future::plan(future::sequential), add = TRUE)

  while (length(pending_ids) > 0L) {
    current_ids <- head(pending_ids, batch_size)
    batch <- future.apply::future_lapply(
      current_ids,
      FUN = run_one_permutation,
      x = x,
      y = y,
      split_cache = split_cache,
      future.seed = TRUE
    )
    permutation_results <- rbind(permutation_results, do.call(rbind, batch))
    permutation_results <- permutation_results[
      order(permutation_results$permutation_id),
      ,
      drop = FALSE
    ]
    rownames(permutation_results) <- NULL

    saveRDS(
      list(
        analysis_signature = analysis_signature,
        permutation_results = permutation_results,
        updated_at = as.character(Sys.time())
      ),
      progress_file
    )
    write.csv(permutation_results, permutation_csv, row.names = FALSE)

    completed_ids <- unique(permutation_results$permutation_id)
    pending_ids <- setdiff(seq_len(n_permutations), completed_ids)
    message(
      "Completed ", nrow(permutation_results), " of ", n_permutations,
      " permutations."
    )
  }
}

permutation_results <- permutation_results[
  permutation_results$permutation_id <= n_permutations,
  ,
  drop = FALSE
]
if (nrow(permutation_results) != n_permutations) {
  stop("The permutation result count does not equal 2,000.")
}

observed_spearman <- observed_metrics$Spearman
observed_pearson <- observed_metrics$Pearson
p_spearman_greater <- (
  sum(permutation_results$Spearman >= observed_spearman) + 1
) / (n_permutations + 1)
p_spearman_less <- (
  sum(permutation_results$Spearman <= observed_spearman) + 1
) / (n_permutations + 1)
p_spearman_two_sided <- min(
  1,
  2 * min(p_spearman_greater, p_spearman_less)
)
p_pearson_greater <- (
  sum(permutation_results$Pearson >= observed_pearson) + 1
) / (n_permutations + 1)

observed_prediction_table <- data.frame(
  sampleid = sample_ids,
  y_true = y,
  observed_pred_recomputed = observed_fit$mean_prediction,
  stringsAsFactors = FALSE
)

summary_table <- data.frame(
  metric = c(
    "observed_Spearman",
    "observed_Pearson",
    "observed_R2",
    "observed_MAE",
    "observed_RMSE",
    "empirical_p_Spearman_greater_primary",
    "empirical_p_Spearman_less",
    "empirical_p_Spearman_two_sided_doubled_tail",
    "empirical_p_Pearson_greater_secondary",
    "permutations_completed",
    "samples",
    "outer_CV",
    "inner_CV",
    "top_variable_genes_per_outer_training_set",
    "alpha_grid"
  ),
  value = c(
    observed_metrics$Spearman,
    observed_metrics$Pearson,
    observed_metrics$R2,
    observed_metrics$MAE,
    observed_metrics$RMSE,
    p_spearman_greater,
    p_spearman_less,
    p_spearman_two_sided,
    p_pearson_greater,
    n_permutations,
    nrow(x),
    paste0(outer_repeats, " repeats x ", outer_folds, " folds"),
    paste0(inner_folds, " folds"),
    top_variable_genes,
    paste(alpha_grid, collapse = ", ")
  ),
  stringsAsFactors = FALSE
)

plot_data <- permutation_results
p_histogram <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Spearman)) +
  ggplot2::geom_histogram(
    bins = 36,
    fill = "#9FC2D3",
    colour = "white",
    linewidth = 0.25
  ) +
  ggplot2::geom_vline(
    xintercept = observed_spearman,
    colour = "#E68613",
    linewidth = 1.0
  ) +
  ggplot2::annotate(
    "label",
    x = -Inf,
    y = Inf,
    hjust = -0.05,
    vjust = 1.15,
    label = paste0(
      "Observed rho = ", sprintf("%.3f", observed_spearman), "\n",
      "Empirical P = ", sprintf("%.4f", p_spearman_greater), "\n",
      n_permutations, " complete model-building permutations"
    ),
    family = "Arial",
    size = 3.0,
    fill = "white"
  ) +
  ggplot2::theme_classic(base_family = "Arial", base_size = 10) +
  ggplot2::labs(
    x = "Permuted cross-validated Spearman rho",
    y = "Count"
  ) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(face = "bold"),
    axis.text = ggplot2::element_text(colour = "black")
  )

ggplot2::ggsave(
  histogram_png,
  p_histogram,
  width = 5.4,
  height = 4.3,
  dpi = 600,
  bg = "white"
)
ggplot2::ggsave(
  histogram_pdf,
  p_histogram,
  width = 5.4,
  height = 4.3,
  device = grDevices::cairo_pdf,
  bg = "white"
)

workbook <- openxlsx::createWorkbook()
openxlsx::addWorksheet(workbook, "summary")
openxlsx::writeData(workbook, "summary", summary_table)
openxlsx::addWorksheet(workbook, "observed_recomputed")
openxlsx::writeData(workbook, "observed_recomputed", observed_metrics)
openxlsx::addWorksheet(workbook, "observed_archived")
openxlsx::writeData(workbook, "observed_archived", archived_metrics)
openxlsx::addWorksheet(workbook, "observed_recomputed_pred")
openxlsx::writeData(workbook, "observed_recomputed_pred", observed_prediction_table)
openxlsx::addWorksheet(workbook, "observed_fold_summary")
openxlsx::writeData(workbook, "observed_fold_summary", observed_fit$fold_summary)
openxlsx::addWorksheet(workbook, "permutation_values")
openxlsx::writeData(workbook, "permutation_values", permutation_results)
openxlsx::addWorksheet(workbook, "analysis_specification")
openxlsx::writeData(
  workbook,
  "analysis_specification",
  data.frame(
    item = c(
      "Expression input",
      "Outcome input",
      "Outcome",
      "Primary statistic",
      "Outer CV",
      "Inner CV",
      "Feature preselection",
      "Alpha candidates",
      "Permutation design",
      "Primary empirical P formula",
      "Two-sided sensitivity formula",
      "Analysis signature"
    ),
    value = c(
      expression_file,
      outcome_file,
      "AUC_recomputed",
      "Spearman correlation between observed AUC and patient-level mean out-of-fold prediction",
      paste0(outer_repeats, " repeated ", outer_folds, "-fold CV"),
      paste0(inner_folds, "-fold CV within each outer training set"),
      paste0("Top ", top_variable_genes, " variable genes selected within each outer training set"),
      paste(alpha_grid, collapse = ", "),
      "Permute AUC_recomputed and repeat feature selection, alpha/lambda tuning, model fitting, and out-of-fold prediction",
      "(b_greater + 1) / (B + 1)",
      "min(1, 2 * min(P_greater, P_less))",
      analysis_signature
    ),
    stringsAsFactors = FALSE
  )
)
openxlsx::saveWorkbook(workbook, results_workbook, overwrite = TRUE)

key_lines <- c(
  "CRS COMPLETE MODEL-BUILDING PERMUTATION RESULTS",
  paste0("Observed Spearman: ", format(observed_spearman, digits = 10)),
  paste0("Observed Pearson: ", format(observed_pearson, digits = 10)),
  paste0("Primary upper-tail empirical P: ", format(p_spearman_greater, digits = 10)),
  paste0("Two-sided doubled-tail sensitivity P: ", format(p_spearman_two_sided, digits = 10)),
  paste0("Permutations: ", n_permutations),
  paste0("Samples: ", nrow(x)),
  paste0("Results workbook: ", results_workbook)
)
writeLines(key_lines, key_values_file)
writeLines(capture.output(sessionInfo()), session_file)

message("CRS permutation analysis completed.")
message("Results workbook: ", results_workbook)
