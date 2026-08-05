# 04_create_lung1_figure.R
#
# Create the Lung1 survival and transportability figure from analysis outputs.

options(stringsAsFactors = FALSE)

required_env_path <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) stop("Environment variable ", name, " is not set.")
  normalizePath(value, winslash = "/", mustWork = TRUE)
}

project_dir <- required_env_path("RADIOGENOMICS_PROJECT_DIR")
lung1_root <- required_env_path("LUNG1_DATA_DIR")
setwd(project_dir)

required_packages <- c("ggplot2", "patchwork")
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

clinical_file <- file.path(lung1_root, "clinical", "Lung1_clinical_clean_initial.rds")
primary_cox_file <- file.path(lung1_root, "survival_validation", "primary_rrs_cox_results.csv")
primary_summary_file <- file.path(lung1_root, "survival_validation", "primary_rrs_os_summary.csv")
signed_cox_file <- file.path(lung1_root, "signed_log_rrs", "lung1_signed_log_rrs_cox_results.csv")
stage_file <- file.path(
  project_dir, "07_results", "lung1_stage_adjusted_survival", "stage_adjusted_rrs_cox_result.csv"
)
domain_summary_file <- file.path(
  lung1_root, "domain_shift_diagnostics", "domain_shift_summary.csv"
)
driver_file <- file.path(
  lung1_root, "domain_shift_diagnostics", "extreme_rrs_feature_contributions.csv"
)
figure_dir <- file.path(project_dir, "08_figures_final", "lung1")
source_dir <- file.path(project_dir, "07_results", "figure_source_data", "lung1")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

for (path in c(
  clinical_file, primary_cox_file, primary_summary_file, signed_cox_file,
  stage_file, domain_summary_file, driver_file
)) {
  if (!file.exists(path)) stop("Required analysis output was not found: ", path)
}

clinical_data <- readRDS(clinical_file)
primary_cox <- read.csv(primary_cox_file, stringsAsFactors = FALSE, check.names = FALSE)
primary_summary <- read.csv(primary_summary_file, stringsAsFactors = FALSE, check.names = FALSE)
signed_cox <- read.csv(signed_cox_file, stringsAsFactors = FALSE, check.names = FALSE)
stage_result <- read.csv(stage_file, stringsAsFactors = FALSE, check.names = FALSE)
domain_summary <- read.csv(domain_summary_file, stringsAsFactors = FALSE, check.names = FALSE)
driver <- read.csv(driver_file, stringsAsFactors = FALSE, check.names = FALSE)

summary_value <- function(item) {
  index <- match(item, primary_summary$item)
  if (is.na(index)) stop("Summary item not found: ", item)
  primary_summary$value[index]
}
domain_value <- function(item) {
  index <- match(item, domain_summary$item)
  if (is.na(index)) stop("Domain-shift item not found: ", item)
  as.numeric(domain_summary$value[index])
}

primary_row <- primary_cox[grepl("RRS_group", primary_cox$term), , drop = FALSE]
if (nrow(primary_row) != 1L) stop("Primary binary RRS Cox result was not uniquely identified.")

signed_term_column <- if ("term" %in% names(signed_cox)) {
  "term"
} else if ("variable" %in% names(signed_cox)) {
  "variable"
} else {
  stop("Signed-log Cox table lacks a term or variable column.")
}

signed_candidates <- signed_cox[
  signed_cox$model == "SignedLog_univariate_binary_RRS_high_vs_low" &
    grepl("SignedLog_RRS_group", signed_cox[[signed_term_column]]),
  ,
  drop = FALSE
]
if (nrow(signed_candidates) != 1L) {
  stop("Signed-log univariate binary Cox result was not uniquely identified.")
}
signed_row <- signed_candidates[1, , drop = FALSE]
if (nrow(stage_result) != 1L) stop("Stage-adjusted Cox result must contain one row.")

forest <- data.frame(
  model = factor(
    c("Primary RRS", "Stage-adjusted primary RRS", "Signed-log sensitivity RRS"),
    levels = rev(c("Primary RRS", "Stage-adjusted primary RRS", "Signed-log sensitivity RRS"))
  ),
  HR = c(primary_row$HR, stage_result$HR, signed_row$HR),
  lower = c(primary_row$CI_lower, stage_result$CI_lower, signed_row$CI_lower),
  upper = c(primary_row$CI_upper, stage_result$CI_upper, signed_row$CI_upper),
  p = c(primary_row$p_value, stage_result$p_value, signed_row$p_value),
  stringsAsFactors = FALSE
)

clinical_records <- nrow(clinical_data)
radiomics_success <- domain_value("lung1_patients_with_features")
os_events <- as.numeric(summary_value("OS events"))
median_os <- as.numeric(summary_value("Median OS, months"))
model_features <- domain_value("model_features")
nonzero_features <- domain_value("nonzero_model_features_excluding_intercept")
features_gt5 <- domain_value("model_features_with_more_than_5pct_Lung1_outside_training_minmax")
nonzero_affected <- domain_value("nonzero_model_features_with_any_Lung1_outside_training_minmax")

contribution_column <- if ("single_feature_contribution" %in% names(driver)) {
  "single_feature_contribution"
} else if ("contribution" %in% names(driver)) {
  "contribution"
} else {
  stop("Feature-contribution table lacks a contribution column.")
}
patient_column <- if ("patient_id" %in% names(driver)) "patient_id" else NULL
feature_column <- if ("feature" %in% names(driver)) "feature" else NULL
if (is.null(feature_column)) stop("Feature-contribution table lacks a feature column.")

driver$contribution_value <- as.numeric(driver[[contribution_column]])
driver <- driver[is.finite(driver$contribution_value), , drop = FALSE]
extreme_driver <- driver[which.max(abs(driver$contribution_value)), , drop = FALSE]
extreme_feature <- gsub("_", " ", sub("^original_", "", extreme_driver[[feature_column]]))
extreme_patient <- if (!is.null(patient_column)) extreme_driver[[patient_column]] else "Lung1 extreme case"

font_family <- if (.Platform$OS.type == "windows") "Arial" else "sans"
base_theme <- theme_classic(base_size = 9, base_family = font_family) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    plot.background = element_rect(fill = "white", colour = NA)
  )

save_panel <- function(plot, name, width = 4.8, height = 4.0) {
  ggsave(file.path(figure_dir, paste0(name, ".png")), plot, width = width, height = height, dpi = 600, bg = "white")
  ggsave(file.path(figure_dir, paste0(name, ".pdf")), plot, width = width, height = height, device = grDevices::cairo_pdf, bg = "white")
}

cohort_cards <- data.frame(
  label = c(
    "Clinical records", "Radiomics success", "OS events",
    "Median OS", "Model-related features", "Final nonzero features"
  ),
  value = c(
    clinical_records, radiomics_success, os_events,
    paste0(sprintf("%.2f", median_os), " months"), model_features, nonzero_features
  ),
  x = rep(1:3, 2),
  y = rep(c(2, 1), each = 3),
  stringsAsFactors = FALSE
)

panel_cohort <- ggplot(cohort_cards, aes(x = x, y = y)) +
  geom_tile(width = 0.86, height = 0.60, fill = "#D6E7EF", colour = "#6C9AB0", linewidth = 0.45) +
  geom_text(aes(label = paste0(label, "\n", value)), lineheight = 0.9, size = 3.0, family = font_family) +
  coord_cartesian(xlim = c(0.45, 3.55), ylim = c(0.55, 2.45), expand = FALSE) +
  theme_void(base_family = font_family)

panel_forest <- ggplot(forest, aes(y = model, x = HR)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey55") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.12, linewidth = 0.7) +
  geom_point(size = 3.1, shape = 21, fill = "#70A6C3") +
  geom_text(
    aes(
      x = max(upper, na.rm = TRUE) + 0.02,
      label = sprintf("HR %.3f (%.3f-%.3f); P %.3f", HR, lower, upper, p)
    ),
    hjust = 0, size = 2.6, family = font_family
  ) +
  coord_cartesian(
    xlim = c(min(lower, na.rm = TRUE) - 0.05, max(upper, na.rm = TRUE) + 0.42),
    clip = "off"
  ) +
  labs(x = "Hazard ratio for overall survival", y = NULL) +
  base_theme + theme(plot.margin = margin(5.5, 85, 5.5, 5.5))

shift <- data.frame(
  category = factor(
    c("Model-related features", "Final nonzero features"),
    levels = c("Final nonzero features", "Model-related features")
  ),
  affected = c(features_gt5, nonzero_affected),
  total = c(model_features, nonzero_features),
  stringsAsFactors = FALSE
)
shift$percent <- 100 * shift$affected / shift$total

panel_shift <- ggplot(shift, aes(y = category, x = percent)) +
  geom_col(width = 0.55, fill = "#E7C3A2") +
  geom_text(
    aes(label = paste0(affected, "/", total)),
    hjust = 1.15, size = 3.1, fontface = "bold"
  ) +
  scale_x_continuous(
    limits = c(0, 105), breaks = c(0, 25, 50, 75, 100),
    labels = function(x) paste0(x, "%")
  ) +
  labs(x = "Features outside the training range in >5% of Lung1 patients", y = NULL) +
  base_theme

extreme <- data.frame(feature = extreme_feature, contribution = extreme_driver$contribution_value)
panel_extreme <- ggplot(extreme, aes(y = feature, x = contribution)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_col(width = 0.52, fill = "#D79C68") +
  geom_text(
    aes(label = sprintf("%.3f", contribution)),
    hjust = ifelse(extreme$contribution < 0, 1.10, -0.10),
    size = 3.1, fontface = "bold"
  ) +
  labs(
    subtitle = extreme_patient,
    x = "Single-feature contribution to RRS",
    y = NULL
  ) +
  coord_cartesian(
    xlim = range(c(extreme$contribution * 1.12, 1)),
    clip = "off"
  ) +
  base_theme + theme(plot.subtitle = element_text(size = 9))

save_panel(panel_cohort, "lung1_cohort_summary")
save_panel(panel_forest, "lung1_os_forest")
save_panel(panel_shift, "lung1_domain_shift_burden")
save_panel(panel_extreme, "lung1_extreme_driver")

combined <- (panel_cohort | panel_forest) / (panel_shift | panel_extreme) +
  patchwork::plot_annotation(tag_levels = "A")
save_panel(combined, "lung1_figure_combined", width = 10.2, height = 8.0)

write.csv(forest, file.path(source_dir, "survival_models.csv"), row.names = FALSE)
write.csv(shift, file.path(source_dir, "domain_shift_burden.csv"), row.names = FALSE)
write.csv(extreme, file.path(source_dir, "extreme_feature_contribution.csv"), row.names = FALSE)
