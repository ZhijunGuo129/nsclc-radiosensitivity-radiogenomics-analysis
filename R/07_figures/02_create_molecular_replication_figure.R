# 02_create_molecular_replication_figure.R
#
# Create the molecular architecture and Lung3 replication figure from
# pathway-level analysis outputs.

options(stringsAsFactors = FALSE)

required_env_path <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) stop("Environment variable ", name, " is not set.")
  normalizePath(value, winslash = "/", mustWork = TRUE)
}

project_dir <- required_env_path("RADIOGENOMICS_PROJECT_DIR")
setwd(project_dir)

required_packages <- c("ggplot2", "patchwork", "tidyr")
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

custom_file <- file.path(
  project_dir, "07_results", "lung3_validation", "predefined_module_replication.csv"
)
hallmark_file <- file.path(
  project_dir, "07_results", "lung3_validation", "hallmark_replication.csv"
)
figure_dir <- file.path(project_dir, "08_figures_final", "molecular_replication")
source_dir <- file.path(project_dir, "07_results", "figure_source_data", "molecular_replication")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

for (path in c(custom_file, hallmark_file)) {
  if (!file.exists(path)) stop("Required analysis output was not found: ", path)
}

custom <- read.csv(custom_file, stringsAsFactors = FALSE, check.names = FALSE)
hallmark <- read.csv(hallmark_file, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c(
  "pathway", "Training_Spearman_r", "Training_Spearman_FDR",
  "Lung3_Spearman_r", "Lung3_Spearman_FDR", "direction_concordant"
)
for (object_name in c("custom", "hallmark")) {
  table <- get(object_name)
  missing <- setdiff(required_columns, names(table))
  if (length(missing) > 0L) {
    stop(object_name, " table is missing: ", paste(missing, collapse = ", "))
  }
}

custom$pathway_label <- gsub("_", " ", custom$pathway)
custom <- custom[order(custom$Training_Spearman_r), , drop = FALSE]
custom$pathway_label <- factor(custom$pathway_label, levels = custom$pathway_label)

hallmark$pathway_label <- gsub("_", " ", sub("^HALLMARK_", "", hallmark$pathway))
hallmark$discovery_selected <- hallmark$Training_Spearman_FDR < 0.10

all_spearman <- cor(
  hallmark$Training_Spearman_r,
  hallmark$Lung3_Spearman_r,
  method = "spearman",
  use = "complete.obs"
)
all_pearson <- cor(
  hallmark$Training_Spearman_r,
  hallmark$Lung3_Spearman_r,
  method = "pearson",
  use = "complete.obs"
)
same_direction <- sum(hallmark$direction_concordant, na.rm = TRUE)

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

panel_discovery <- ggplot(custom, aes(x = Training_Spearman_r, y = pathway_label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey65") +
  geom_col(width = 0.62, fill = "#6FA6C1") +
  geom_point(
    aes(fill = Training_Spearman_FDR < 0.05),
    shape = 21, size = 2.3, colour = "black"
  ) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#E67E00")) +
  labs(x = "NSCLC-Radiogenomics Spearman rho", y = NULL) +
  base_theme + theme(legend.position = "none", axis.text.y = element_text(size = 7.0))

panel_lung3 <- ggplot(custom, aes(x = Lung3_Spearman_r, y = pathway_label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey65") +
  geom_col(width = 0.62, fill = "#B7CFE0") +
  geom_point(
    aes(fill = Lung3_Spearman_FDR < 0.05),
    shape = 21, size = 2.3, colour = "black"
  ) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#E67E00")) +
  labs(x = "Lung3 Spearman rho", y = NULL) +
  base_theme + theme(legend.position = "none", axis.text.y = element_text(size = 7.0))

panel_concordance <- ggplot(
  hallmark,
  aes(x = Training_Spearman_r, y = Lung3_Spearman_r, fill = discovery_selected)
) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey70") +
  geom_point(shape = 21, size = 2.5, colour = "black", alpha = 0.85) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "#E67E00") +
  annotate(
    "text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.15,
    label = sprintf(
      "Spearman rho = %.3f\nPearson r = %.3f\nSame direction = %d/50",
      all_spearman, all_pearson, same_direction
    ),
    size = 3.0, family = font_family
  ) +
  scale_fill_manual(values = c("FALSE" = "grey80", "TRUE" = "#7EB08A")) +
  labs(x = "NSCLC-Radiogenomics Spearman rho", y = "Lung3 Spearman rho") +
  base_theme + theme(legend.position = "none")

representative <- hallmark[order(hallmark$Training_Spearman_FDR, -abs(hallmark$Training_Spearman_r)), , drop = FALSE]
representative <- head(representative, 12)
heat <- representative[, c("pathway_label", "Training_Spearman_r", "Lung3_Spearman_r")]
heat <- tidyr::pivot_longer(
  heat,
  cols = c("Training_Spearman_r", "Lung3_Spearman_r"),
  names_to = "cohort",
  values_to = "rho"
)
heat$cohort <- factor(
  heat$cohort,
  levels = c("Training_Spearman_r", "Lung3_Spearman_r"),
  labels = c("NSCLC-Radiogenomics", "Lung3")
)
heat$pathway_label <- factor(
  heat$pathway_label,
  levels = rev(unique(representative$pathway_label))
)

panel_heatmap <- ggplot(heat, aes(x = cohort, y = pathway_label, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.7, family = font_family) +
  scale_fill_gradient2(
    low = "#3B73B9", mid = "white", high = "#C95C54", midpoint = 0,
    name = "Spearman rho"
  ) +
  labs(x = NULL, y = NULL) +
  base_theme +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(size = 6.8),
    legend.position = "bottom"
  )

save_panel(panel_discovery, "discovery_module_correlations")
save_panel(panel_lung3, "lung3_module_correlations")
save_panel(panel_concordance, "hallmark_effect_concordance")
save_panel(panel_heatmap, "representative_hallmark_heatmap")

combined <- (panel_discovery | panel_lung3) / (panel_concordance | panel_heatmap) +
  patchwork::plot_annotation(tag_levels = "A")
save_panel(combined, "molecular_replication_figure_combined", width = 10.2, height = 8.2)

write.csv(custom, file.path(source_dir, "predefined_module_replication.csv"), row.names = FALSE)
write.csv(hallmark, file.path(source_dir, "hallmark_replication.csv"), row.names = FALSE)
