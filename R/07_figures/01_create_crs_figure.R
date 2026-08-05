# 01_create_crs_figure.R
#
# Generate the final CRS model figure from strict audit outputs.
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

required_packages <- c("ggplot2", "dplyr", "openxlsx", "scales", "grid", "patchwork")
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

library(ggplot2)
library(dplyr)
library(openxlsx)
library(scales)
library(grid)
library(patchwork)

############################################################
# 1. Paths
############################################################

cell_score_file <- file.path(
  project_dir,
  "05_molecular_scores",
  "cell_line_crs_scores.rds"
)

permutation_workbook <- file.path(
  project_dir,
  "07_results",
  "crs_permutation",
  "crs_full_workflow_permutation.xlsx"
)

fig_dir <- file.path(
  project_dir,
  "08_figures_final",
  "figure2_crs"
)

audit_dir <- file.path(
  project_dir,
  "07_results",
  "figure2_crs"
)

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(cell_score_file)) {
  stop("Missing cell-line score file: ", cell_score_file)
}

if (!file.exists(permutation_workbook)) {
  stop(
    "Missing CRS full-workflow permutation workbook. Run the full-workflow CRS permutation script first: ",
    permutation_workbook
  )
}

# 2. Figure palette
############################################################

pal <- c(
  blue = "#1F77B4",
  orange = "#FF7F0E",
  grey_dark = "#333333",
  grey_mid = "#999999",
  grey_light = "#D9D9D9",
  grey_pale = "#F5F5F5",
  white = "#FFFFFF"
)

# Unified visual rules
col_blue_main <- pal["blue"]
col_orange_main <- pal["orange"]

col_blue_card_sample <- alpha(col_blue_main, 0.10)
col_blue_card_core <- alpha(col_blue_main, 0.17)
col_blue_card_p <- alpha(col_blue_main, 0.11)

col_blue_arrow <- alpha(col_blue_main, 0.55)
col_blue_p_border <- alpha(col_blue_main, 0.85)

col_point_B <- col_blue_main
col_point_D <- col_blue_main
col_hist_C <- col_blue_main
col_fit <- col_orange_main
col_ci <- col_orange_main
theme_fig2 <- function(base_size = 8.5, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title = element_text(
        face = "bold",
        size = base_size + 0.2,
        color = "black"
      ),
      axis.text = element_text(
        size = base_size,
        color = "black"
      ),
      axis.line = element_line(
        linewidth = 0.42,
        color = "black"
      ),
      axis.ticks = element_line(
        linewidth = 0.34,
        color = "black"
      ),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "none",
      plot.margin = margin(7, 7, 7, 7)
    )
}

save_panel <- function(plot_obj, filename_base, width = 4.9, height = 4.1) {
  
  pdf_path <- file.path(fig_dir, paste0(filename_base, ".pdf"))
  png_path <- file.path(fig_dir, paste0(filename_base, ".png"))
  
  tryCatch(
    {
      ggsave(
        filename = pdf_path,
        plot = plot_obj,
        width = width,
        height = height,
        device = cairo_pdf,
        bg = "white"
      )
    },
    error = function(e) {
      ggsave(
        filename = pdf_path,
        plot = plot_obj,
        width = width,
        height = height,
        bg = "white"
      )
    }
  )
  
  ggsave(
    filename = png_path,
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 600,
    bg = "white"
  )
  
  cat("Saved:\n")
  cat(pdf_path, "\n")
  cat(png_path, "\n\n")
}

############################################################
# 3. Load observed nested-CV and permutation results
############################################################

observed_df <- read.xlsx(
  permutation_workbook,
  sheet = "observed_recomputed_pred",
  check.names = FALSE
)

required_observed_cols <- c(
  "sampleid",
  "y_true",
  "observed_pred_recomputed"
)
missing_observed <- setdiff(required_observed_cols, colnames(observed_df))
if (length(missing_observed) > 0) {
  stop(
    "Strict observed prediction sheet is missing columns: ",
    paste(missing_observed, collapse = ", ")
  )
}

cv_plot <- data.frame(
  sampleid = as.character(observed_df$sampleid),
  observed = suppressWarnings(as.numeric(observed_df$y_true)),
  predicted = suppressWarnings(as.numeric(observed_df$observed_pred_recomputed)),
  stringsAsFactors = FALSE
)
cv_plot <- cv_plot[
  is.finite(cv_plot$observed) & is.finite(cv_plot$predicted),
  ,
  drop = FALSE
]

n_repeat <- 3L
n_fold <- 5L

cv_pearson <- as.numeric(
  suppressWarnings(cor.test(cv_plot$observed, cv_plot$predicted, method = "pearson"))$estimate
)
cv_spearman <- as.numeric(
  suppressWarnings(cor.test(
    cv_plot$observed,
    cv_plot$predicted,
    method = "spearman",
    exact = FALSE
  ))$estimate
)

perm_df <- read.xlsx(
  permutation_workbook,
  sheet = "permutation_values",
  check.names = FALSE
)
if (!"spearman" %in% colnames(perm_df)) {
  stop("Strict permutation sheet does not contain column: spearman")
}
perm_df$perm_spearman <- suppressWarnings(as.numeric(perm_df$spearman))
perm_df <- perm_df[is.finite(perm_df$perm_spearman), , drop = FALSE]
n_perm <- nrow(perm_df)

empirical_p_spearman_greater <- (
  sum(perm_df$perm_spearman >= cv_spearman, na.rm = TRUE) + 1
) / (n_perm + 1)

lower_tail_p <- (
  sum(perm_df$perm_spearman <= cv_spearman, na.rm = TRUE) + 1
) / (n_perm + 1)
empirical_p_spearman_two_sided <- min(
  1,
  2 * min(empirical_p_spearman_greater, lower_tail_p)
)

# 5. Load final CRS prediction for lung-lineage subset
############################################################

cell_df <- readRDS(cell_score_file)
cell_df <- as.data.frame(cell_df, check.names = FALSE)

required_cell_cols <- c("AUC_recomputed", "CRS_predicted_AUC", "Primarysite")
missing_cell <- setdiff(required_cell_cols, colnames(cell_df))

if (length(missing_cell) > 0) {
  stop("cell-line score table missing columns: ", paste(missing_cell, collapse = ", "))
}

lung_plot <- data.frame(
  observed = suppressWarnings(as.numeric(cell_df$AUC_recomputed)),
  predicted = suppressWarnings(as.numeric(cell_df$CRS_predicted_AUC)),
  Primarysite = as.character(cell_df$Primarysite),
  stringsAsFactors = FALSE
)

lung_plot <- lung_plot[
  is.finite(lung_plot$observed) &
    is.finite(lung_plot$predicted),
]

lung_plot <- lung_plot[
  grepl("lung|nsclc|non.small|pulmonary", lung_plot$Primarysite, ignore.case = TRUE),
]

lung_pearson <- as.numeric(
  suppressWarnings(cor.test(
    lung_plot$observed,
    lung_plot$predicted,
    method = "pearson"
  ))$estimate
)

lung_spearman <- as.numeric(
  suppressWarnings(cor.test(
    lung_plot$observed,
    lung_plot$predicted,
    method = "spearman",
    exact = FALSE
  ))$estimate
)

############################################################
# 6. Shared axes for Figure 2B and Figure 2D
############################################################

all_x <- c(cv_plot$observed, lung_plot$observed)
all_y <- c(cv_plot$predicted, lung_plot$predicted)

x_lim <- c(
  floor(min(all_x, na.rm = TRUE) * 10) / 10 - 0.1,
  ceiling(max(all_x, na.rm = TRUE) * 10) / 10 + 0.1
)

y_lim <- c(
  floor(min(all_y, na.rm = TRUE) * 10) / 10 - 0.1,
  ceiling(max(all_y, na.rm = TRUE) * 10) / 10 + 0.1
)

x_lim[1] <- max(0, x_lim[1])
y_lim[1] <- max(0, y_lim[1])

x_breaks <- seq(
  ceiling(x_lim[1]),
  floor(x_lim[2]),
  by = 1
)

y_breaks <- seq(
  ceiling(y_lim[1]),
  floor(y_lim[2]),
  by = 1
)

box_x <- x_lim[1] + 0.68 * diff(x_lim)
box_y <- y_lim[1] + 0.96 * diff(y_lim)

############################################################
# 7. Figure 2A: refined model summary cards, no panel tag
############################################################

cards <- data.frame(
  x = c(1, 2, 3, 1, 2, 3),
  y = c(2, 2, 2, 1, 1, 1),
  type = c(
    "sample", "sample", "sample",
    "core", "core", "pvalue"
  ),
  title = c(
    "Cleveland\nRadioSet",
    "Modeling\nsamples",
    "Lung-lineage\nsubset",
    "CV Spearman\nrho",
    "CV Pearson\nr",
    "Empirical\np"
  ),
  value = c(
    "540\ncell lines",
    "504\ncell lines",
    "95\ncell lines",
    sprintf("%.3f", cv_spearman),
    sprintf("%.3f", cv_pearson),
    sprintf("%.4f", empirical_p_spearman_greater)
  ),
  stringsAsFactors = FALSE
)

cards$fill <- ifelse(
  cards$type == "core",
  col_blue_card_core,
  ifelse(cards$type == "pvalue", col_blue_card_p, col_blue_card_sample)
)

cards$border <- ifelse(
  cards$type == "pvalue",
  col_blue_p_border,
  col_blue_main
)

cards$line_width <- ifelse(
  cards$type == "core",
  0.75,
  ifelse(cards$type == "pvalue", 0.58, 0.55)
)

cards$value_size <- ifelse(
  cards$type == "core",
  3.35,
  ifelse(cards$type == "pvalue", 2.98, 3.00)
)

cards$title_size <- ifelse(
  cards$type == "core",
  2.18,
  ifelse(cards$type == "pvalue", 2.12, 2.20)
)

pA <- ggplot() +
  annotate(
    "text",
    x = 2,
    y = 2.50,
    label = "Dataset composition",
    size = 2.85,
    fontface = "bold",
    color = pal["grey_dark"]
  ) +
  annotate(
    "text",
    x = 2,
    y = 1.50,
    label = "Model performance for predicting radiation AUC",
    size = 2.85,
    fontface = "bold",
    color = pal["grey_dark"]
  ) +
  geom_tile(
    data = cards,
    aes(
      x = x,
      y = y,
      fill = fill,
      color = border,
      linewidth = line_width
    ),
    width = 0.82,
    height = 0.62
  ) +
  geom_segment(
    aes(x = 1.44, xend = 1.56, y = 2, yend = 2),
    arrow = arrow(length = unit(0.045, "inches"), type = "closed"),
    linewidth = 0.24,
    color = col_blue_arrow,
    lineend = "round"
  ) +
  geom_segment(
    aes(x = 2.44, xend = 2.56, y = 2, yend = 2),
    arrow = arrow(length = unit(0.045, "inches"), type = "closed"),
    linewidth = 0.24,
    color = col_blue_arrow,
    lineend = "round"
  ) +
  geom_text(
    data = cards,
    aes(
      x = x,
      y = y + 0.105,
      label = title,
      size = title_size
    ),
    lineheight = 0.88,
    fontface = "bold",
    color = pal["grey_dark"]
  ) +
  geom_text(
    data = cards,
    aes(
      x = x,
      y = y - 0.145,
      label = value,
      size = value_size
    ),
    lineheight = 0.90,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_linewidth_identity() +
  scale_size_identity() +
  coord_cartesian(
    xlim = c(0.42, 3.58),
    ylim = c(0.55, 2.62),
    expand = FALSE,
    clip = "on"
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.margin = margin(7, 7, 7, 7)
  )

save_panel(
  pA,
  "Figure2_model_summary",
  width = 4.9,
  height = 3.55
)

############################################################
# 8. Figure 2B: mean CV prediction scatter, no panel tag
############################################################

label_B <- paste0(
  "Mean CV prediction\n",
  n_repeat, " x ", n_fold, "-fold CV\n",
  "Cell lines = ", nrow(cv_plot), "\n",
  "Pearson r = ", sprintf("%.3f", cv_pearson), "\n",
  "Spearman rho = ", sprintf("%.3f", cv_spearman), "\n",
  "Linear fit with 95% CI"
)

pB <- ggplot(cv_plot, aes(x = observed, y = predicted)) +
  geom_point(
    shape = 16,
    color = col_point_B,
    size = 1.25,
    alpha = 0.34
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = col_fit,
    fill = col_ci,
    alpha = 0.14,
    linewidth = 1.05
  ) +
  annotate(
    "label",
    x = box_x,
    y = box_y,
    label = label_B,
    hjust = 0,
    vjust = 1,
    size = 2.05,
    lineheight = 0.92,
    linewidth = 0.18,
    label.r = unit(0.08, "lines"),
    fill = alpha("white", 0.97),
    color = "black"
  ) +
  coord_cartesian(
    xlim = x_lim,
    ylim = y_lim,
    clip = "on"
  ) +
  scale_x_continuous(
    breaks = x_breaks,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = y_breaks,
    expand = c(0, 0)
  ) +
  theme_fig2(base_size = 8.5) +
  labs(
    x = "Observed radiation AUC",
    y = "Mean CV-predicted radiation AUC"
  )

save_panel(
  pB,
  "Figure2_cross_validated_performance",
  width = 4.9,
  height = 4.1
)

############################################################
# 9. Figure 2C: 2000-permutation histogram, no panel tag
############################################################

hist_info <- hist(
  perm_df$perm_spearman,
  breaks = 35,
  plot = FALSE
)

max_count <- max(hist_info$counts, na.rm = TRUE)

y_max_C <- ceiling(max_count * 1.16 / 10) * 10
y_max_C <- max(y_max_C, 330)

x_left_C <- min(
  -0.18,
  floor(min(perm_df$perm_spearman, na.rm = TRUE) * 10) / 10 - 0.02
)

x_right_C <- 0.45

label_C <- paste0(
  "Observed rho = ", sprintf("%.3f", cv_spearman), "\n",
  "Empirical p = ", sprintf("%.4f", empirical_p_spearman_greater), "\n",
  n_perm, " permutations"
)

main_label_x <- x_left_C + 0.02
main_label_y <- y_max_C * 0.94

observed_label_x <- cv_spearman - 0.012
observed_label_y <- y_max_C * 0.94

pC <- ggplot(perm_df, aes(x = perm_spearman)) +
  geom_histogram(
    bins = 35,
    fill = col_hist_C,
    color = "white",
    linewidth = 0.25,
    alpha = 0.38
  ) +
  geom_vline(
    xintercept = cv_spearman,
    color = col_orange_main,
    linewidth = 1.05
  ) +
  annotate(
    "label",
    x = main_label_x,
    y = main_label_y,
    label = label_C,
    hjust = 0,
    vjust = 1,
    size = 2.05,
    lineheight = 0.92,
    linewidth = 0.18,
    label.r = unit(0.08, "lines"),
    fill = alpha("white", 0.97),
    color = "black"
  ) +
  annotate(
    "label",
    x = observed_label_x,
    y = observed_label_y,
    label = "Observed",
    hjust = 1,
    vjust = 1,
    size = 2.05,
    linewidth = 0.18,
    label.r = unit(0.08, "lines"),
    fill = alpha("white", 0.97),
    color = col_orange_main,
    fontface = "bold"
  ) +
  coord_cartesian(
    xlim = c(x_left_C, x_right_C),
    ylim = c(0, y_max_C),
    clip = "on"
  ) +
  scale_x_continuous(
    breaks = c(-0.1, 0, 0.1, 0.2, 0.3, 0.4),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq(0, floor(y_max_C / 50) * 50, by = 50),
    expand = c(0, 0)
  ) +
  theme_fig2(base_size = 8.5) +
  theme(
    plot.margin = margin(10, 7, 7, 7)
  ) +
  labs(
    x = "Permuted Spearman rho",
    y = "Count"
  )

save_panel(
  pC,
  "Figure2_permutation_distribution",
  width = 4.9,
  height = 4.1
)

############################################################
# 10. Figure 2D: lung-lineage final CRS prediction, no panel tag
############################################################

label_D <- paste0(
  "Final CRS prediction\n",
  "Within-cohort subset\n",
  "Cell lines = ", nrow(lung_plot), "\n",
  "Pearson r = ", sprintf("%.3f", lung_pearson), "\n",
  "Spearman rho = ", sprintf("%.3f", lung_spearman), "\n",
  "Linear fit with 95% CI"
)

pD <- ggplot(lung_plot, aes(x = observed, y = predicted)) +
  geom_point(
    shape = 16,
    color = col_point_D,
    size = 1.60,
    alpha = 0.48
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = col_fit,
    fill = col_ci,
    alpha = 0.15,
    linewidth = 1.00
  ) +
  annotate(
    "label",
    x = box_x,
    y = box_y,
    label = label_D,
    hjust = 0,
    vjust = 1,
    size = 2.05,
    lineheight = 0.92,
    linewidth = 0.18,
    label.r = unit(0.08, "lines"),
    fill = alpha("white", 0.97),
    color = "black"
  ) +
  coord_cartesian(
    xlim = x_lim,
    ylim = y_lim,
    clip = "on"
  ) +
  scale_x_continuous(
    breaks = x_breaks,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = y_breaks,
    expand = c(0, 0)
  ) +
  theme_fig2(base_size = 8.5) +
  labs(
    x = "Observed radiation AUC",
    y = "Final CRS-predicted radiation AUC"
  )

save_panel(
  pD,
  "Figure2_lung_lineage_concordance",
  width = 4.9,
  height = 4.1
)

############################################################
# 11. Combined Figure 2
############################################################


pA_comb <- pA +
  theme(
    plot.margin = margin(
      t = 16,
      r = 10,
      b = 8,
      l = 12
    )
  )

pB_comb <- pB +
  theme(
    plot.margin = margin(
      t = 16,
      r = 8,
      b = 8,
      l = 30
    )
  )

pC_comb <- pC +
  theme(
    plot.margin = margin(
      t = 18,
      r = 10,
      b = 8,
      l = 22
    )
  )

pD_comb <- pD +
  theme(
    plot.margin = margin(
      t = 18,
      r = 8,
      b = 8,
      l = 30
    )
  )

fig2_combined <- (
  pA_comb + pB_comb
) / (
  pC_comb + pD_comb
) +
  plot_layout(
    widths = c(1.02, 1.00),
    heights = c(0.96, 1.04)
  ) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 15,
      color = "black",
      family = "Arial"
    ),
    plot.tag.position = c(0.012, 0.988)
  )

combined_pdf <- file.path(
  fig_dir,
  "figure2_crs_combined.pdf"
)

combined_png <- file.path(
  fig_dir,
  "figure2_crs_combined.png"
)

combined_tiff <- file.path(
  fig_dir,
  "figure2_crs_combined.tiff"
)

tryCatch(
  {
    ggsave(
      filename = combined_pdf,
      plot = fig2_combined,
      width = 10.8,
      height = 8.6,
      device = cairo_pdf,
      bg = "white"
    )
  },
  error = function(e) {
    ggsave(
      filename = combined_pdf,
      plot = fig2_combined,
      width = 10.8,
      height = 8.6,
      bg = "white"
    )
  }
)

ggsave(
  filename = combined_png,
  plot = fig2_combined,
  width = 10.8,
  height = 8.6,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = combined_tiff,
  plot = fig2_combined,
  width = 10.8,
  height = 8.6,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

cat("\nCombined Figure 2 saved:\n")
cat(combined_pdf, "\n")
cat(combined_png, "\n")
cat(combined_tiff, "\n")

############################################################
# 12. Save audit workbook
############################################################

summary_df <- data.frame(
  metric = c(
    "Figure2B_unique_cell_lines",
    "Figure2B_raw_prediction_rows",
    "Figure2B_repeats",
    "Figure2B_folds",
    "Figure2B_Pearson_from_plotted_504_points",
    "Figure2B_Spearman_from_plotted_504_points",
    "Figure2C_n_permutations",
    "Figure2C_empirical_p_spearman_greater",
    "Figure2C_empirical_p_spearman_two_sided",
    "Figure2C_max_histogram_count",
    "Figure2C_y_axis_upper_limit",
    "Figure2C_x_axis_left_limit",
    "Figure2C_x_axis_right_limit",
    "Figure2D_lung_lineage_n",
    "Figure2D_Pearson",
    "Figure2D_Spearman",
    "Scatter_x_axis_min",
    "Scatter_x_axis_max",
    "Scatter_y_axis_min",
    "Scatter_y_axis_max",
    "Main_blue",
    "Main_orange",
    "B_point_alpha",
    "D_point_alpha",
    "C_histogram_alpha",
    "B_CI_alpha",
    "D_CI_alpha",
    "clip_setting",
    "smooth_band_definition"
  ),
  value = c(
    nrow(cv_plot),
    nrow(cv_raw),
    n_repeat,
    n_fold,
    cv_pearson,
    cv_spearman,
    n_perm,
    empirical_p_spearman_greater,
    empirical_p_spearman_two_sided,
    max_count,
    y_max_C,
    x_left_C,
    x_right_C,
    nrow(lung_plot),
    lung_pearson,
    lung_spearman,
    x_lim[1],
    x_lim[2],
    y_lim[1],
    y_lim[2],
    col_blue_main,
    col_orange_main,
    0.34,
    0.48,
    0.38,
    0.14,
    0.15,
    "on",
    "Linear fit with 95% confidence interval"
  ),
  stringsAsFactors = FALSE
)

wb <- createWorkbook()

addWorksheet(wb, "Figure2B_raw_CV")
writeData(wb, "Figure2B_raw_CV", cv_raw)

addWorksheet(wb, "Figure2B_cellline_mean")
writeData(wb, "Figure2B_cellline_mean", cv_plot)

addWorksheet(wb, "Figure2C_permutation_values")
writeData(wb, "Figure2C_permutation_values", perm_df)

addWorksheet(wb, "Figure2D_lung_lineage")
writeData(wb, "Figure2D_lung_lineage", lung_plot)

addWorksheet(wb, "Figure2A_cards")
writeData(wb, "Figure2A_cards", cards)

addWorksheet(wb, "summary")
writeData(wb, "summary", summary_df)

addWorksheet(wb, "method_note")
writeData(
  wb,
  "method_note",
  data.frame(
    note = c(
      "Figure 2B uses one cell-line-level mean held-out CV prediction per sample.",
      "Figure 2C uses 2,000 full model-building permutations in which radiation AUC labels are permuted and the complete repeated cross-validation workflow is rerun.",
      "Each permutation reruns feature selection, tuning, and Elastic Net model fitting within the repeated cross-validation workflow.",
      "The smooth bands in Figure 2B and Figure 2D are 95% confidence intervals of the linear fit.",
      "Separated panels intentionally omit A/B/C/D panel labels.",
      "Combined Figure 2 uses patchwork-generated A/B/C/D labels."
    ),
    stringsAsFactors = FALSE
  )
)

for (sheet in names(wb)) {
  setColWidths(wb, sheet, cols = 1:60, widths = "auto")
}

audit_file <- file.path(
  audit_dir,
  "figure2_source_data.xlsx"
)

saveWorkbook(
  wb,
  file = audit_file,
  overwrite = TRUE
)

############################################################
# 13. Done
############################################################

cat("\n===== Figure 2 generation completed =====\n")

cat("\nFigure output folder:\n")
cat(fig_dir, "\n\n")

cat("Generated figure files:\n")
print(list.files(fig_dir, full.names = FALSE))

cat("\nAudit workbook:\n")
cat(audit_file, "\n\n")

cat("Key final values:\n")
cat("Figure2A/B/C Spearman rho:", cv_spearman, "\n")
cat("Figure2A/B Pearson r:", cv_pearson, "\n")
cat("Figure2C empirical p greater:", empirical_p_spearman_greater, "\n")
cat("Figure2C empirical p two-sided:", empirical_p_spearman_two_sided, "\n")
cat("Figure2C max histogram count:", max_count, "\n")
cat("Figure2C y-axis upper limit:", y_max_C, "\n")
cat("Figure2C x-axis:", x_left_C, "to", x_right_C, "\n")
cat("Figure2D Pearson r:", lung_pearson, "\n")
cat("Figure2D Spearman rho:", lung_spearman, "\n")
