# 03_create_rrs_figure.R
#
# Create the primary RRS figure from analysis outputs.

options(stringsAsFactors = FALSE)

required_env_path <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) stop("Environment variable ", name, " is not set.")
  normalizePath(value, winslash = "/", mustWork = TRUE)
}

project_dir <- required_env_path("RADIOGENOMICS_PROJECT_DIR")
setwd(project_dir)

required_packages <- c("ggplot2", "patchwork", "openxlsx", "pROC")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0L) {
  stop(
    "Missing R package(s): ", paste(missing_packages, collapse = ", "),
    ". Install dependencies with environment/install_r_dependencies.R."
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

figure_dir <- file.path(project_dir, "08_figures_final", "rrs")
source_dir <- file.path(project_dir, "07_results", "figure_source_data", "rrs")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

permutation_file <- file.path(
  project_dir, "07_results", "rrs_permutation", "rrs_full_workflow_permutation.xlsx"
)
coefficient_file <- file.path(
  project_dir, "07_results", "rrs_model", "final_model_coefficients.csv"
)
size_file <- file.path(
  project_dir, "07_results", "rrs_size_adjustment", "rrs_size_adjustment_key_values.csv"
)

for (path in c(permutation_file, coefficient_file, size_file)) {
  if (!file.exists(path)) stop("Required analysis output was not found: ", path)
}

metric_value <- function(table, metric) {
  index <- match(metric, table$metric)
  if (is.na(index)) stop("Metric not found: ", metric)
  value <- suppressWarnings(as.numeric(table$value[index]))
  if (!is.finite(value)) stop("Metric is not numeric: ", metric)
  value
}

summary_table <- openxlsx::read.xlsx(permutation_file, sheet = "summary")
predictions <- openxlsx::read.xlsx(
  permutation_file, sheet = "observed_patient_predictions"
)
coefficients <- read.csv(coefficient_file, stringsAsFactors = FALSE, check.names = FALSE)
size_values <- read.csv(size_file, stringsAsFactors = FALSE, check.names = FALSE)

required_prediction_columns <- c("patient_id", "observed_CRS_z", "RRS_CV_mean")
missing_prediction_columns <- setdiff(required_prediction_columns, names(predictions))
if (length(missing_prediction_columns) > 0L) {
  stop("Prediction table is missing: ", paste(missing_prediction_columns, collapse = ", "))
}
if (!all(c("feature", "coefficient") %in% names(coefficients))) {
  stop("Coefficient table must contain feature and coefficient columns.")
}

predictions$observed_CRS_z <- as.numeric(predictions$observed_CRS_z)
predictions$RRS_CV_mean <- as.numeric(predictions$RRS_CV_mean)
predictions <- predictions[
  is.finite(predictions$observed_CRS_z) & is.finite(predictions$RRS_CV_mean),
  , drop = FALSE
]

n_patients <- nrow(predictions)
pearson_r <- cor(predictions$observed_CRS_z, predictions$RRS_CV_mean, method = "pearson")
spearman_rho <- cor(predictions$observed_CRS_z, predictions$RRS_CV_mean, method = "spearman")
p_upper <- metric_value(summary_table, "empirical_p_Spearman_greater_primary")
p_two_sided <- metric_value(summary_table, "empirical_p_Spearman_two_sided_doubled_tail")

crs_high <- predictions$observed_CRS_z >= median(predictions$observed_CRS_z)
roc_fit <- pROC::roc(
  response = crs_high,
  predictor = predictions$RRS_CV_mean,
  levels = c(FALSE, TRUE),
  direction = "<",
  quiet = TRUE
)
auc_value <- as.numeric(pROC::auc(roc_fit))
auc_ci <- as.numeric(pROC::ci.auc(roc_fit, method = "delong"))

nonzero <- coefficients[
  coefficients$feature != "(Intercept)" & coefficients$coefficient != 0,
  , drop = FALSE
]
nonzero <- nonzero[order(abs(nonzero$coefficient), decreasing = TRUE), , drop = FALSE]
nonzero$feature_label <- gsub("_", " ", sub("^original_", "", nonzero$feature))
nonzero$feature_label <- factor(
  nonzero$feature_label,
  levels = rev(nonzero$feature_label)
)

size_only_r2 <- metric_value(size_values, "size_only_R2")
strict_r2 <- metric_value(size_values, "size_plus_strict_RRS_R2")
strict_delta <- metric_value(size_values, "strict_RRS_delta_R2")
signed_r2 <- metric_value(size_values, "size_plus_signed_log_RRS_R2")
signed_delta <- metric_value(size_values, "signed_log_RRS_delta_R2")

font_family <- if (.Platform$OS.type == "windows") "Arial" else "sans"
base_theme <- theme_classic(base_size = 9, base_family = font_family) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    plot.background = element_rect(fill = "white", colour = NA)
  )

save_panel <- function(plot, name, width = 4.8, height = 4.0) {
  ggsave(
    file.path(figure_dir, paste0(name, ".png")), plot,
    width = width, height = height, dpi = 600, bg = "white"
  )
  ggsave(
    file.path(figure_dir, paste0(name, ".pdf")), plot,
    width = width, height = height, device = grDevices::cairo_pdf, bg = "white"
  )
}

cards <- data.frame(
  label = c(
    "Patients", "Radiomic features", "Final nonzero features",
    "Pearson r", "Spearman rho", "Empirical P",
    "AUC", "AUC 95% CI"
  ),
  value = c(
    n_patients, 107, nrow(nonzero),
    sprintf("%.3f", pearson_r), sprintf("%.3f", spearman_rho),
    sprintf("%.4f", p_upper), sprintf("%.3f", auc_value),
    sprintf("%.3f-%.3f", auc_ci[1], auc_ci[3])
  ),
  x = rep(1:4, 2),
  y = rep(c(2, 1), each = 4),
  stringsAsFactors = FALSE
)

panel_summary <- ggplot(cards, aes(x = x, y = y)) +
  geom_tile(
    width = 0.88, height = 0.62,
    fill = "#D6E7EF", colour = "#6C9AB0", linewidth = 0.45
  ) +
  geom_text(
    aes(label = paste0(label, "\n", value)),
    lineheight = 0.90, size = 3.0, family = font_family
  ) +
  coord_cartesian(xlim = c(0.45, 4.55), ylim = c(0.55, 2.45), expand = FALSE) +
  theme_void(base_family = font_family)

panel_predictions <- ggplot(
  predictions,
  aes(x = observed_CRS_z, y = RRS_CV_mean)
) +
  geom_point(size = 2.0, alpha = 0.78, colour = "#6FA6C1") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "#E67E00") +
  annotate(
    "text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.15,
    label = sprintf(
      "Pearson r = %.3f\nSpearman rho = %.3f\nUpper-tail P = %.4f\nTwo-sided P = %.4f",
      pearson_r, spearman_rho, p_upper, p_two_sided
    ),
    size = 3.0, family = font_family
  ) +
  labs(x = "Observed CRS_z", y = "Mean out-of-fold RRS prediction") +
  base_theme

panel_coefficients <- ggplot(
  nonzero,
  aes(x = coefficient, y = feature_label, fill = coefficient > 0)
) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_col(width = 0.65, colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = c("FALSE" = "#6FA6C1", "TRUE" = "#D18472")) +
  labs(x = "Elastic-net coefficient", y = NULL) +
  base_theme +
  theme(legend.position = "none", axis.text.y = element_text(size = 7.4))

size_models <- data.frame(
  model = factor(
    c("Size only", "Size + primary RRS", "Size + signed-log RRS"),
    levels = c("Size only", "Size + primary RRS", "Size + signed-log RRS")
  ),
  r2 = c(size_only_r2, strict_r2, signed_r2),
  delta = c(NA_real_, strict_delta, signed_delta),
  stringsAsFactors = FALSE
)

panel_size <- ggplot(size_models, aes(x = model, y = r2, fill = model)) +
  geom_col(width = 0.65, colour = "black", linewidth = 0.3) +
  geom_text(
    aes(label = sprintf("R² = %.3f", r2)),
    vjust = -0.45, size = 3.0, family = font_family
  ) +
  geom_text(
    data = size_models[is.finite(size_models$delta), , drop = FALSE],
    aes(label = sprintf("ΔR² = %.3f", delta)),
    vjust = 1.55, size = 2.8, family = font_family
  ) +
  scale_fill_manual(values = c(
    "Size only" = "#D5D5D5",
    "Size + primary RRS" = "#6098B7",
    "Size + signed-log RRS" = "#AFC9D8"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Model R²") +
  base_theme +
  theme(legend.position = "none", axis.text.x = element_text(size = 7.7))

save_panel(panel_summary, "rrs_model_summary")
save_panel(panel_predictions, "rrs_cross_validated_predictions")
save_panel(panel_coefficients, "rrs_model_coefficients")
save_panel(panel_size, "rrs_tumor_size_adjustment")

combined <- (panel_summary | panel_predictions) / (panel_coefficients | panel_size) +
  patchwork::plot_annotation(tag_levels = "A")
save_panel(combined, "rrs_figure_combined", width = 10.2, height = 8.0)

write.csv(predictions, file.path(source_dir, "cross_validated_predictions.csv"), row.names = FALSE)
write.csv(nonzero, file.path(source_dir, "model_coefficients.csv"), row.names = FALSE)
write.csv(size_models, file.path(source_dir, "tumor_size_models.csv"), row.names = FALSE)
write.csv(
  data.frame(
    metric = c(
      "n", "pearson_r", "spearman_rho", "upper_tail_p", "two_sided_p",
      "auc", "auc_ci_lower", "auc_ci_upper"
    ),
    value = c(
      n_patients, pearson_r, spearman_rho, p_upper, p_two_sided,
      auc_value, auc_ci[1], auc_ci[3]
    )
  ),
  file.path(source_dir, "figure_metrics.csv"),
  row.names = FALSE
)
