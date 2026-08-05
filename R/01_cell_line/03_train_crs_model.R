# 03_train_crs_model.R

options(stringsAsFactors = FALSE)
options(timeout = 1800)

required_env_path <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    stop("Environment variable ", name, " is not set.")
  }
  normalizePath(value, winslash = "/", mustWork = TRUE)
}

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

project_dir <- required_env_path("RADIOGENOMICS_PROJECT_DIR")
setwd(project_dir)

set.seed(20260722)

input_expression <- file.path(
  project_dir,
  "01_raw_data",
  "RadioGx",
  "RadioGx_rnaseq_model_gene_by_sample.rds"
)
input_outcomes <- file.path(
  project_dir,
  "05_molecular_scores",
  "RadioGx_sensitivity_model_outcomes.rds"
)

result_dir <- file.path(project_dir, "07_results", "crs_model")
figure_dir <- file.path(project_dir, "08_figures_final", "crs_model")
model_dir <- file.path(project_dir, "06_models")
score_dir <- file.path(project_dir, "05_molecular_scores")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(score_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_expression)) {
  stop("Expression matrix not found: ", input_expression)
}
if (!file.exists(input_outcomes)) {
  stop("Radiation-response outcome table not found: ", input_outcomes)
}

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

select_top_variable_genes <- function(x, top_n = 2000L) {
  variances <- apply(x, 2, stats::var, na.rm = TRUE)
  variances[!is.finite(variances)] <- 0
  variances <- variances[variances > 0]
  top_n <- min(as.integer(top_n), length(variances))
  names(sort(variances, decreasing = TRUE))[seq_len(top_n)]
}

make_balanced_folds <- function(n, k, seed) {
  set.seed(seed)
  sample(rep(seq_len(k), length.out = n))
}

fit_elastic_net <- function(
    x,
    y,
    alpha_grid = c(0.5, 1),
    inner_folds = 5L,
    fold_seed = 20260722L
) {
  fold_id <- make_balanced_folds(length(y), inner_folds, fold_seed)
  best <- NULL
  best_cvm <- Inf

  for (alpha_value in alpha_grid) {
    fit <- glmnet::cv.glmnet(
      x = as.matrix(x),
      y = as.numeric(y),
      family = "gaussian",
      alpha = alpha_value,
      foldid = fold_id,
      standardize = TRUE,
      type.measure = "mse"
    )

    current_cvm <- min(fit$cvm, na.rm = TRUE)
    if (is.finite(current_cvm) && current_cvm < best_cvm) {
      best_cvm <- current_cvm
      best <- list(
        fit = fit,
        alpha = alpha_value,
        lambda = fit$lambda.min,
        cvm = current_cvm
      )
    }
  }

  if (is.null(best)) {
    stop("Elastic-net tuning failed for every alpha candidate.")
  }
  best
}

repeated_outer_cv <- function(
    x,
    y,
    repeats = 3L,
    outer_folds = 5L,
    inner_folds = 5L,
    top_n = 2000L,
    alpha_grid = c(0.5, 1),
    seed = 20260722L
) {
  rows <- vector("list", repeats * outer_folds)
  row_id <- 0L

  for (repeat_id in seq_len(repeats)) {
    outer_fold_id <- make_balanced_folds(
      nrow(x),
      outer_folds,
      seed + repeat_id
    )

    for (fold_id in seq_len(outer_folds)) {
      row_id <- row_id + 1L
      test_index <- which(outer_fold_id == fold_id)
      train_index <- setdiff(seq_len(nrow(x)), test_index)

      selected_features <- select_top_variable_genes(
        x[train_index, , drop = FALSE],
        top_n = top_n
      )

      tuned_model <- fit_elastic_net(
        x = x[train_index, selected_features, drop = FALSE],
        y = y[train_index],
        alpha_grid = alpha_grid,
        inner_folds = inner_folds,
        fold_seed = seed + repeat_id * 100L + fold_id
      )

      predictions <- as.numeric(
        stats::predict(
          tuned_model$fit,
          newx = as.matrix(x[test_index, selected_features, drop = FALSE]),
          s = "lambda.min"
        )
      )

      rows[[row_id]] <- data.frame(
        repeat_id = repeat_id,
        fold_id = fold_id,
        sample_index = test_index,
        sampleid = rownames(x)[test_index],
        y_true = y[test_index],
        y_pred = predictions,
        alpha = tuned_model$alpha,
        lambda = tuned_model$lambda,
        preselected_features = length(selected_features),
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, rows)
}

performance_summary <- function(cv_predictions) {
  sample_mean <- aggregate(
    cbind(y_true, y_pred) ~ sampleid,
    data = cv_predictions,
    FUN = mean
  )

  data.frame(
    n = nrow(sample_mean),
    Spearman = suppressWarnings(stats::cor(
      sample_mean$y_true,
      sample_mean$y_pred,
      method = "spearman"
    )),
    Pearson = suppressWarnings(stats::cor(
      sample_mean$y_true,
      sample_mean$y_pred,
      method = "pearson"
    )),
    MAE = mean(abs(sample_mean$y_true - sample_mean$y_pred)),
    RMSE = sqrt(mean((sample_mean$y_true - sample_mean$y_pred)^2)),
    stringsAsFactors = FALSE
  )
}

expression_raw <- readRDS(input_expression)
outcomes <- readRDS(input_outcomes)
expression_clean <- clean_expression_matrix(expression_raw)

common_samples <- intersect(colnames(expression_clean), outcomes$sampleid)
expression_samples <- t(expression_clean[, common_samples, drop = FALSE])
outcome_rows <- outcomes[match(common_samples, outcomes$sampleid), , drop = FALSE]

prepare_endpoint <- function(endpoint) {
  y <- as.numeric(outcome_rows[[endpoint]])
  keep <- is.finite(y) & complete.cases(expression_samples)
  list(
    x = expression_samples[keep, , drop = FALSE],
    y = y[keep],
    outcomes = outcome_rows[keep, , drop = FALSE]
  )
}

auc_data <- prepare_endpoint("AUC_recomputed")
sf2_data <- prepare_endpoint("SF2")

if (nrow(auc_data$x) != 504L) {
  warning("Expected 504 AUC modeling samples; found ", nrow(auc_data$x), ".")
}

auc_cv <- repeated_outer_cv(auc_data$x, auc_data$y)
sf2_cv <- repeated_outer_cv(sf2_data$x, sf2_data$y)
auc_performance <- performance_summary(auc_cv)
sf2_performance <- performance_summary(sf2_cv)

final_features <- select_top_variable_genes(auc_data$x, top_n = 2000L)
final_model <- fit_elastic_net(
  x = auc_data$x[, final_features, drop = FALSE],
  y = auc_data$y,
  alpha_grid = c(0.5, 1),
  inner_folds = 5L,
  fold_seed = 20260722L
)

final_predictions <- as.numeric(
  stats::predict(
    final_model$fit,
    newx = as.matrix(auc_data$x[, final_features, drop = FALSE]),
    s = "lambda.min"
  )
)

coefficient_matrix <- as.matrix(
  stats::coef(final_model$fit, s = "lambda.min")
)
coefficients <- data.frame(
  feature = rownames(coefficient_matrix),
  coefficient = as.numeric(coefficient_matrix[, 1]),
  stringsAsFactors = FALSE
)
nonzero_coefficients <- coefficients[coefficients$coefficient != 0, , drop = FALSE]

cell_line_scores <- data.frame(
  sampleid = rownames(auc_data$x),
  CRS_predicted_AUC = final_predictions,
  AUC_recomputed = auc_data$y,
  SF2 = auc_data$outcomes$SF2,
  CellLine = auc_data$outcomes$CellLine,
  Primarysite = auc_data$outcomes$Primarysite,
  Histology = auc_data$outcomes$Histology,
  Subhistology = auc_data$outcomes$Subhistology,
  is_lung_strict = auc_data$outcomes$is_lung_strict,
  stringsAsFactors = FALSE
)

lung_scores <- cell_line_scores[cell_line_scores$is_lung_strict, , drop = FALSE]
lung_summary <- data.frame(
  lung_n = nrow(lung_scores),
  Spearman = suppressWarnings(stats::cor(
    lung_scores$CRS_predicted_AUC,
    lung_scores$AUC_recomputed,
    method = "spearman",
    use = "complete.obs"
  )),
  Pearson = suppressWarnings(stats::cor(
    lung_scores$CRS_predicted_AUC,
    lung_scores$AUC_recomputed,
    method = "pearson",
    use = "complete.obs"
  )),
  stringsAsFactors = FALSE
)

model_object <- list(
  final_model = final_model,
  final_features = final_features,
  endpoint = "AUC_recomputed",
  direction = "Higher CRS indicates higher predicted radiation AUC.",
  gene_id_type = "Ensembl",
  seed = 20260722L
)

saveRDS(model_object, file.path(model_dir, "crs_elastic_net_model.rds"))
saveRDS(cell_line_scores, file.path(score_dir, "cell_line_crs_scores.rds"))
openxlsx::write.xlsx(
  cell_line_scores,
  file.path(score_dir, "cell_line_crs_scores.xlsx"),
  overwrite = TRUE
)

workbook <- openxlsx::createWorkbook()
for (sheet in c(
  "observed_cv",
  "observed_performance",
  "sf2_cv",
  "sf2_performance",
  "final_features",
  "final_coefficients",
  "final_nonzero_coefficients",
  "cell_line_scores",
  "lung_lineage_summary"
)) {
  openxlsx::addWorksheet(workbook, sheet)
}
openxlsx::writeData(workbook, "observed_cv", auc_cv)
openxlsx::writeData(workbook, "observed_performance", auc_performance)
openxlsx::writeData(workbook, "sf2_cv", sf2_cv)
openxlsx::writeData(workbook, "sf2_performance", sf2_performance)
openxlsx::writeData(
  workbook,
  "final_features",
  data.frame(feature = final_features, stringsAsFactors = FALSE)
)
openxlsx::writeData(workbook, "final_coefficients", coefficients)
openxlsx::writeData(workbook, "final_nonzero_coefficients", nonzero_coefficients)
openxlsx::writeData(workbook, "cell_line_scores", cell_line_scores)
openxlsx::writeData(workbook, "lung_lineage_summary", lung_summary)
openxlsx::saveWorkbook(
  workbook,
  file.path(result_dir, "crs_model_results.xlsx"),
  overwrite = TRUE
)

sample_mean_auc <- aggregate(
  cbind(y_true, y_pred) ~ sampleid,
  data = auc_cv,
  FUN = mean
)

png(
  file.path(figure_dir, "cross_validated_predictions.png"),
  width = 1800,
  height = 1600,
  res = 300
)
plot(
  sample_mean_auc$y_true,
  sample_mean_auc$y_pred,
  xlab = "Observed radiation AUC",
  ylab = "Mean out-of-fold prediction",
  pch = 16
)
abline(stats::lm(y_pred ~ y_true, data = sample_mean_auc), lwd = 2)
dev.off()

message("CRS model training completed.")
message("Model results: ", file.path(result_dir, "crs_model_results.xlsx"))
message("Model object: ", file.path(model_dir, "crs_elastic_net_model.rds"))
