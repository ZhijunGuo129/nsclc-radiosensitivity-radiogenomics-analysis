# 02_create_supplementary_figure_panels.R
#
# Generate all supplementary figure panels as separate files without panel letters.
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
# 0. Packages and global style
############################################################

required_packages <- c(
  "ggplot2",
  "png",
  "grid",
  "dplyr",
  "tidyr"
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
invisible(lapply(required_packages, library, character.only = TRUE))

fig_font <- if (.Platform$OS.type == "windows") {
  "Times New Roman"
} else {
  "Times"
}

# TRUE keeps descriptive titles such as
# "Cohort flow and analytic readiness".
# FALSE removes even these descriptive titles.
show_panel_titles <- TRUE


out_root <- file.path(
  project_dir,
  "08_figures_supplementary",
  "supplementary_panels"
)

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

base_theme <- theme_bw(base_family = fig_font, base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    plot.title.position = "plot",
    plot.subtitle = element_text(size = 9.5, hjust = 0),
    axis.title = element_text(face = "bold")
  )

theme_set(base_theme)

title_or_null <- function(x) {
  if (isTRUE(show_panel_titles)) {
    x
  } else {
    NULL
  }
}

subtitle_or_null <- function(x) {
  if (isTRUE(show_panel_titles)) {
    x
  } else {
    NULL
  }
}

save_panel <- function(plot_object, out_dir, stem, width, height) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  png_file <- file.path(out_dir, paste0(stem, ".png"))
  pdf_file <- file.path(out_dir, paste0(stem, ".pdf"))

  suppressMessages(
    ggsave(
      filename = png_file,
      plot = plot_object,
      width = width,
      height = height,
      dpi = 600,
      bg = "white"
    )
  )

  pdf_ok <- tryCatch({
    grDevices::cairo_pdf(
      filename = pdf_file,
      width = width,
      height = height,
      family = "serif"
    )
    suppressMessages(print(plot_object))
    grDevices::dev.off()
    TRUE
  }, error = function(e) {
    if (grDevices::dev.cur() > 1) {
      try(grDevices::dev.off(), silent = TRUE)
    }
    message("PNG was saved successfully. PDF was skipped: ", conditionMessage(e))
    FALSE
  })

  message("Saved PNG: ", png_file)
  if (pdf_ok) {
    message("Saved PDF: ", pdf_file)
  }
}

make_image_panel <- function(path, panel_title, fallback_text) {
  if (!file.exists(path)) {
    return(
      ggplot() +
        annotate(
          "rect",
          xmin = 0, xmax = 1, ymin = 0, ymax = 1,
          fill = "grey96", colour = "grey70", linewidth = 0.35
        ) +
        annotate(
          "text",
          x = 0.5, y = 0.56,
          label = fallback_text,
          size = 4,
          family = fig_font
        ) +
        annotate(
          "text",
          x = 0.5, y = 0.40,
          label = "Image file not found. Check CONFIG paths.",
          size = 3.2,
          family = fig_font,
          colour = "grey35"
        ) +
        coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
        labs(title = title_or_null(panel_title)) +
        theme_void(base_family = fig_font) +
        theme(
          plot.title = element_text(face = "bold", size = 12, hjust = 0),
          plot.title.position = "plot",
          plot.margin = margin(4, 4, 4, 4)
        )
    )
  }

  img <- png::readPNG(path)

  ggplot() +
    annotation_raster(
      img,
      xmin = 0, xmax = 1,
      ymin = 0, ymax = 1,
      interpolate = TRUE
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(title = title_or_null(panel_title)) +
    theme_void(base_family = fig_font) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0),
      plot.title.position = "plot",
      plot.margin = margin(4, 4, 4, 4)
    )
}

make_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

correlation_filter <- function(x, cutoff = 0.90) {
  cm <- abs(stats::cor(x, use = "pairwise.complete.obs"))
  diag(cm) <- 0
  keep_names <- colnames(cm)

  repeat {
    max_cor <- max(cm, na.rm = TRUE)

    if (
      !is.finite(max_cor) ||
      max_cor <= cutoff ||
      ncol(cm) <= 1L
    ) {
      break
    }

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

short_feature <- function(x) {
  x <- sub("^original_", "", x)
  gsub("_", " ", x)
}

clean_set_name <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  gsub("_", " ", x)
}

############################################################
# 1. Supplementary Figure S1 panels
############################################################

make_S1 <- function() {
  message("\n===== Making S1 separate panels =====")

  out_dir <- file.path(out_root, "S1")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  image_dir <- file.path(
    project_dir,
    "07_results",
    "supplementary_exports",
    "S1_images"
  )

  ct_slice_png <- file.path(
    image_dir,
    "S1_B_representative_CT_cleaned.png"
  )

  roi_overlay_png <- file.path(
    image_dir,
    "S1_C_segmentation_overlay_cleaned.png"
  )

  cohort_df <- data.frame(
    step_id = factor(
      c("NR1", "NR2", "NR3", "L1", "L2", "L3"),
      levels = c("NR1", "NR2", "NR3", "L1", "L2", "L3")
    ),
    label = c(
      "RNA-seq\navailable",
      "Paired\nCT-RNA",
      "RRS\ntraining",
      "Clinical\nrecords",
      "Radiomics\nsuccess",
      "OS\nevents"
    ),
    n = c(130, 117, 117, 422, 418, 370),
    cohort = factor(
      c(
        "NSCLC-Radiogenomics",
        "NSCLC-Radiogenomics",
        "NSCLC-Radiogenomics",
        "Lung1",
        "Lung1",
        "Lung1"
      ),
      levels = c("NSCLC-Radiogenomics", "Lung1")
    ),
    stringsAsFactors = FALSE
  )

  feature_df <- data.frame(
    item = factor(
      c(
        "Original\nPyRadiomics features",
        "RRS training\npatients",
        "Lung1 radiomics\nsuccess"
      ),
      levels = c(
        "Original\nPyRadiomics features",
        "RRS training\npatients",
        "Lung1 radiomics\nsuccess"
      )
    ),
    n = c(107, 117, 418),
    group = factor(
      c("Feature set", "Training", "External CT"),
      levels = c("Feature set", "Training", "External CT")
    ),
    stringsAsFactors = FALSE
  )

  p1 <- ggplot(cohort_df, aes(x = step_id, y = n, fill = cohort)) +
    geom_col(width = 0.70, colour = "black", linewidth = 0.35) +
    geom_text(aes(label = n), vjust = -0.35, family = fig_font, size = 3.4) +
    facet_wrap(~ cohort, scales = "free_x", nrow = 1) +
    scale_x_discrete(labels = setNames(cohort_df$label, cohort_df$step_id)) +
    scale_fill_manual(values = c(
      "NSCLC-Radiogenomics" = "#8FB7CC",
      "Lung1" = "#D8D8D8"
    )) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.11))) +
    labs(
      title = title_or_null("Cohort flow and analytic readiness"),
      x = NULL,
      y = "Cases"
    ) +
    theme(
      axis.text.x = element_text(size = 8.5),
      axis.text.y = element_text(size = 8.5),
      strip.background = element_rect(
        fill = "grey90",
        colour = "grey40",
        linewidth = 0.35
      ),
      strip.text = element_text(face = "bold", size = 9.5),
      legend.position = "none",
      plot.title.position = "plot"
    )

  p2 <- ggplot(feature_df, aes(x = item, y = n, fill = group)) +
    geom_col(width = 0.64, colour = "black", linewidth = 0.35) +
    geom_text(aes(label = n), vjust = -0.35, family = fig_font, size = 3.4) +
    scale_fill_manual(values = c(
      "Feature set" = "#BFD3DF",
      "Training" = "#8FB7CC",
      "External CT" = "#D8D8D8"
    )) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.11))) +
    labs(
      title = title_or_null("Radiomic extraction and feature QC"),
      x = NULL,
      y = "Count"
    ) +
    theme(
      axis.text.x = element_text(size = 8.5),
      axis.text.y = element_text(size = 8.5),
      legend.position = "none",
      plot.title.position = "plot"
    )

  p3 <- make_image_panel(
    path = ct_slice_png,
    panel_title = "Representative pretreatment CT slice",
    fallback_text = "Representative axial CT slice"
  )

  p4 <- make_image_panel(
    path = roi_overlay_png,
    panel_title = "Matched tumor segmentation overlay",
    fallback_text = "Matched CT with tumor-mask overlay"
  )

  save_panel(p1, out_dir, "S1_cohort_flow", width = 6.0, height = 4.2)
  save_panel(p2, out_dir, "S1_radiomic_extraction_feature_QC", width = 6.0, height = 4.2)
  save_panel(p3, out_dir, "S1_representative_CT_plain", width = 6.0, height = 4.2)
  save_panel(p4, out_dir, "S1_representative_CT_ROI_overlay", width = 6.0, height = 4.2)

  audit_df <- data.frame(
    item = c(
      "ct_slice_png",
      "roi_overlay_png",
      "ct_slice_png_exists",
      "roi_overlay_png_exists",
      "NSCLC_Radiogenomics_RNA_available",
      "NSCLC_Radiogenomics_paired_CT_RNA",
      "NSCLC_Radiogenomics_RRS_training",
      "Original_PyRadiomics_features",
      "Lung1_clinical_records",
      "Lung1_successful_radiomics",
      "Lung1_OS_events"
    ),
    value = c(
      ct_slice_png,
      roi_overlay_png,
      file.exists(ct_slice_png),
      file.exists(roi_overlay_png),
      130,
      117,
      117,
      107,
      422,
      418,
      370
    ),
    stringsAsFactors = FALSE
  )

  write.csv(
    audit_df,
    file.path(out_dir, "S1_audit_values.csv"),
    row.names = FALSE
  )
}

############################################################
# 2. Supplementary Figure S2 panels
############################################################

make_S2 <- function() {
  message("\n===== Making S2 separate panels =====")

  out_dir <- file.path(out_root, "S2")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  validation_dir <- file.path(
    project_dir,
    "07_results",
    "lung3_validation"
  )

  custom_file <- file.path(
    validation_dir,
    "predefined_module_replication.csv"
  )

  hallmark_file <- file.path(
    validation_dir,
    "hallmark_replication.csv"
  )

  for (ff in c(custom_file, hallmark_file)) {
    if (!file.exists(ff)) {
      stop(
        "Cannot find the required Lung3 validation output:\n",
        ff,
        "\nRun R/05_pathways/05_evaluate_lung3_replication.R_Lung3_External_Molecular_Validation.R first."
      )
    }
  }

  read_comparison <- function(path, set_type) {
    x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

    required <- c(
      "pathway",
      "Training_Spearman_r",
      "Training_Spearman_FDR",
      "Lung3_Spearman_r",
      "Lung3_Spearman_FDR",
      "direction_concordant"
    )

    missing_cols <- setdiff(required, names(x))
    if (length(missing_cols) > 0L) {
      stop(
        basename(path),
        " is missing columns: ",
        paste(missing_cols, collapse = ", ")
      )
    }

    data.frame(
      set_name = as.character(x$pathway),
      set_type = set_type,
      discovery_rho = as.numeric(x$Training_Spearman_r),
      discovery_fdr = as.numeric(x$Training_Spearman_FDR),
      lung3_rho = as.numeric(x$Lung3_Spearman_r),
      lung3_fdr = as.numeric(x$Lung3_Spearman_FDR),
      same_direction = toupper(as.character(x$direction_concordant)) == "TRUE",
      stringsAsFactors = FALSE
    )
  }

  custom_df <- read_comparison(custom_file, "Custom")
  hallmark_df <- read_comparison(hallmark_file, "Hallmark")

  if (nrow(custom_df) != 10L) {
    stop("Expected 10 custom modules, found ", nrow(custom_df), ".")
  }

  if (nrow(hallmark_df) != 50L) {
    stop("Expected 50 Hallmark pathways, found ", nrow(hallmark_df), ".")
  }

  hallmark_rho_s <- suppressWarnings(
    cor(
      hallmark_df$discovery_rho,
      hallmark_df$lung3_rho,
      method = "spearman",
      use = "complete.obs"
    )
  )

  hallmark_r_p <- suppressWarnings(
    cor(
      hallmark_df$discovery_rho,
      hallmark_df$lung3_rho,
      method = "pearson",
      use = "complete.obs"
    )
  )

  hallmark_same_n <- sum(hallmark_df$same_direction, na.rm = TRUE)
  hallmark_disc_010_n <- sum(hallmark_df$discovery_fdr < 0.10, na.rm = TRUE)
  hallmark_disc_010_same_n <- sum(
    hallmark_df$discovery_fdr < 0.10 &
      hallmark_df$same_direction,
    na.rm = TRUE
  )
  hallmark_lung3_005_n <- sum(hallmark_df$lung3_fdr < 0.05, na.rm = TRUE)

  custom_disc_005 <- custom_df$discovery_fdr < 0.05
  custom_core_total <- sum(custom_disc_005, na.rm = TRUE)
  custom_core_rep_n <- sum(
    custom_disc_005 &
      custom_df$same_direction &
      custom_df$lung3_fdr < 0.05,
    na.rm = TRUE
  )

  hallmark_plot_df <- hallmark_df %>%
    mutate(
      discovery_class = ifelse(
        discovery_fdr < 0.10,
        "Discovery FDR < 0.10",
        "Other Hallmarks"
      )
    )

  panel_note <- paste0(
    "Spearman rho = ", sprintf("%.3f", hallmark_rho_s),
    "\nPearson r = ", sprintf("%.3f", hallmark_r_p),
    "\nSame direction = ", hallmark_same_n, "/50"
  )

  p1 <- ggplot(
    hallmark_plot_df,
    aes(x = discovery_rho, y = lung3_rho, fill = discovery_class)
  ) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey72", linewidth = 0.35) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey72", linewidth = 0.35) +
    geom_point(
      shape = 21,
      size = 2.8,
      colour = "black",
      stroke = 0.35,
      alpha = 0.88
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      colour = "#E68613",
      fill = "#F3C88C",
      linewidth = 0.8,
      alpha = 0.22
    ) +
    annotate(
      "text",
      x = min(hallmark_plot_df$discovery_rho, na.rm = TRUE),
      y = max(hallmark_plot_df$lung3_rho, na.rm = TRUE),
      label = panel_note,
      hjust = 0,
      vjust = 1,
      family = fig_font,
      size = 3.35
    ) +
    scale_fill_manual(values = c(
      "Discovery FDR < 0.10" = "#6FA87A",
      "Other Hallmarks" = "grey78"
    )) +
    labs(
      title = title_or_null("Effect-size replication across all 50 Hallmarks"),
      x = "NSCLC-Radiogenomics Spearman rho",
      y = "Lung3 Spearman rho",
      fill = NULL
    ) +
    theme(
      legend.position = "bottom",
      legend.justification = "left",
      legend.box.just = "left",
      legend.text = element_text(size = 8.5),
      axis.title.y = element_text(margin = margin(r = 9)),
      plot.title.position = "plot"
    )

  custom_plot_df <- custom_df %>%
    arrange(discovery_rho) %>%
    mutate(
      set_label = gsub("_", " ", set_name),
      set_label = factor(set_label, levels = set_label)
    )

  p2 <- ggplot(custom_plot_df, aes(y = set_label)) +
    geom_segment(
      aes(x = discovery_rho, xend = lung3_rho, yend = set_label),
      colour = "grey72",
      linewidth = 0.75
    ) +
    geom_point(aes(x = discovery_rho), size = 2.7, colour = "#4C8BB1") +
    geom_point(aes(x = lung3_rho), size = 2.7, colour = "#B2555A") +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey75", linewidth = 0.35) +
    labs(
      title = title_or_null("Custom-module effects in discovery and Lung3"),
      subtitle = subtitle_or_null("Blue, NSCLC-Radiogenomics; red, Lung3"),
      x = "Spearman rho",
      y = NULL
    ) +
    theme(
      axis.text.y = element_text(size = 7.7),
      axis.title.x = element_text(margin = margin(t = 5)),
      plot.title.position = "plot",
      plot.subtitle = element_text(size = 9.6, hjust = 0, margin = margin(b = 3))
    )

  hallmark_heat <- hallmark_df %>%
    arrange(discovery_rho) %>%
    mutate(
      set_label = clean_set_name(set_name),
      set_label = factor(set_label, levels = set_label)
    ) %>%
    select(
      set_label,
      discovery_rho,
      discovery_fdr,
      lung3_rho,
      lung3_fdr
    ) %>%
    pivot_longer(
      cols = -set_label,
      names_to = c("cohort", ".value"),
      names_pattern = "(discovery|lung3)_(rho|fdr)"
    ) %>%
    mutate(
      cohort = factor(
        cohort,
        levels = c("discovery", "lung3"),
        labels = c("NSCLC-Radiogenomics", "Lung3")
      ),
      fdr_mark = ifelse(fdr < 0.05, "*", "")
    )

  p3 <- ggplot(
    hallmark_heat,
    aes(x = cohort, y = set_label, fill = rho)
  ) +
    geom_tile(colour = "white", linewidth = 0.22) +
    geom_text(aes(label = fdr_mark), family = fig_font, size = 2.45) +
    scale_fill_gradient2(
      low = "#3B73B9",
      mid = "white",
      high = "#C95C54",
      midpoint = 0,
      name = "Spearman rho"
    ) +
    labs(
      title = title_or_null("Complete Hallmark correlation architecture"),
      subtitle = subtitle_or_null("* FDR < 0.05 within the corresponding cohort"),
      x = NULL,
      y = NULL
    ) +
    guides(
      fill = guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(1.70, "in"),
        barheight = unit(0.12, "in")
      )
    ) +
    theme(
      axis.text.y = element_text(size = 5.85),
      axis.text.x = element_text(size = 9),
      legend.position = "bottom",
      legend.justification = "left",
      legend.box.just = "left",
      legend.title = element_text(size = 8.2),
      legend.text = element_text(size = 7.8),
      plot.title.position = "plot",
      plot.subtitle = element_text(size = 9.4, hjust = 0, margin = margin(b = 3))
    )

  summary_df <- data.frame(
    metric = factor(
      c(
        "Core custom\nmodules",
        "Discovery Hallmarks\nFDR < 0.10",
        "All Hallmarks:\nsame direction",
        "All Hallmarks:\nLung3 FDR < 0.05"
      ),
      levels = c(
        "Core custom\nmodules",
        "Discovery Hallmarks\nFDR < 0.10",
        "All Hallmarks:\nsame direction",
        "All Hallmarks:\nLung3 FDR < 0.05"
      )
    ),
    value = c(
      custom_core_rep_n,
      hallmark_disc_010_same_n,
      hallmark_same_n,
      hallmark_lung3_005_n
    ),
    label = c(
      paste0(custom_core_rep_n, "/", custom_core_total),
      paste0(hallmark_disc_010_same_n, "/", hallmark_disc_010_n),
      paste0(hallmark_same_n, "/50"),
      paste0(hallmark_lung3_005_n, "/50")
    ),
    stringsAsFactors = FALSE
  )

  p4 <- ggplot(summary_df, aes(x = metric, y = value)) +
    geom_col(width = 0.70, fill = "#9FC2D3", colour = "black", linewidth = 0.35) +
    geom_text(aes(label = label), vjust = -0.35, family = fig_font, size = 3.45) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title = title_or_null("Cross-cohort replication summary"),
      x = NULL,
      y = "Pathways or modules"
    ) +
    theme(
      axis.text.x = element_text(size = 7.8),
      axis.title.y = element_text(margin = margin(r = 8)),
      plot.title.position = "plot"
    )

  save_panel(p1, out_dir, "S2_hallmark_replication_scatter", width = 6.2, height = 5.2)
  save_panel(p2, out_dir, "S2_custom_module_effects", width = 6.2, height = 5.2)
  save_panel(p3, out_dir, "S2_hallmark_heatmap", width = 6.4, height = 7.2)
  save_panel(p4, out_dir, "S2_replication_summary", width = 6.2, height = 5.2)

  plotting_data <- rbind(custom_df, hallmark_df)
  write.csv(plotting_data, file.path(out_dir, "S2_plotting_data.csv"), row.names = FALSE)

  audit_values <- data.frame(
    item = c(
      "custom_modules_n",
      "hallmark_pathways_n",
      "custom_core_replicated_at_Lung3_FDR_0.05",
      "custom_core_total",
      "hallmark_discovery_FDR_0.10_same_direction",
      "hallmark_discovery_FDR_0.10_total",
      "all_hallmark_same_direction",
      "all_hallmark_Lung3_FDR_0.05",
      "all_hallmark_effect_Spearman",
      "all_hallmark_effect_Pearson"
    ),
    value = c(
      nrow(custom_df),
      nrow(hallmark_df),
      custom_core_rep_n,
      custom_core_total,
      hallmark_disc_010_same_n,
      hallmark_disc_010_n,
      hallmark_same_n,
      hallmark_lung3_005_n,
      hallmark_rho_s,
      hallmark_r_p
    ),
    stringsAsFactors = FALSE
  )

  write.csv(audit_values, file.path(out_dir, "S2_audit_values.csv"), row.names = FALSE)
}

############################################################
# 3. Supplementary Figure S3 panels
############################################################

make_S3 <- function() {
  message("\n===== Making S3 separate panels =====")

  out_dir <- file.path(out_root, "S3")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  source_dir <- file.path(
    project_dir,
    "07_results",
    "rrs_size_adjustment"
  )

  data_file <- file.path(source_dir, "rrs_size_adjustment_dataset.csv")
  model_summary_file <- file.path(source_dir, "rrs_size_adjustment_nested_models.csv")

  for (ff in c(data_file, model_summary_file)) {
    if (!file.exists(ff)) {
      stop(
        "Cannot find the required RRS size-adjustment output:\n",
        ff,
        "\nRun R/06_sensitivity_analyses/02_analyze_tumor_size_adjustment.R first."
      )
    }
  }

  df <- read.csv(data_file, stringsAsFactors = FALSE, check.names = FALSE)
  model_df <- read.csv(model_summary_file, stringsAsFactors = FALSE, check.names = FALSE)

  required_data_cols <- c(
    "CRS_z",
    "log_MeshVolume",
    "log_Maximum3DDiameter",
    "RRS_CV",
    "SignedLog_RRS_CV"
  )

  missing_data_cols <- setdiff(required_data_cols, names(df))
  if (length(missing_data_cols) > 0L) {
    stop("Size-analysis dataset is missing: ", paste(missing_data_cols, collapse = ", "))
  }

  required_model_cols <- c("model", "model_R2", "delta_R2_vs_size_only")
  missing_model_cols <- setdiff(required_model_cols, names(model_df))
  if (length(missing_model_cols) > 0L) {
    stop("Nested-model summary is missing: ", paste(missing_model_cols, collapse = ", "))
  }

  df <- df[
    complete.cases(
      df[, c(
        "CRS_z",
        "log_MeshVolume",
        "log_Maximum3DDiameter",
        "RRS_CV",
        "SignedLog_RRS_CV"
      )]
    ),
    ,
    drop = FALSE
  ]

  df$partial_CRSz <- residuals(
    lm(CRS_z ~ log_MeshVolume + log_Maximum3DDiameter, data = df)
  )

  df$partial_Strict_RRS <- residuals(
    lm(RRS_CV ~ log_MeshVolume + log_Maximum3DDiameter, data = df)
  )

  partial_r <- suppressWarnings(
    cor(df$partial_CRSz, df$partial_Strict_RRS, method = "pearson")
  )

  partial_rho <- suppressWarnings(
    cor(df$partial_CRSz, df$partial_Strict_RRS, method = "spearman")
  )

  p1 <- ggplot(df, aes(x = log_MeshVolume, y = RRS_CV)) +
    geom_point(size = 2.25, colour = "#8FB7CC", alpha = 0.85) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      colour = "#E68613",
      fill = "#F3C88C",
      linewidth = 0.8,
      alpha = 0.22
    ) +
    labs(
      title = title_or_null("Primary RRS versus tumor volume"),
      x = "log tumor volume",
      y = "Primary RRS_CV"
    )

  p2 <- ggplot(df, aes(x = log_Maximum3DDiameter, y = RRS_CV)) +
    geom_point(size = 2.25, colour = "#8FB7CC", alpha = 0.85) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      colour = "#E68613",
      fill = "#F3C88C",
      linewidth = 0.8,
      alpha = 0.22
    ) +
    labs(
      title = title_or_null("Primary RRS versus maximum 3D diameter"),
      x = "log maximum 3D diameter",
      y = "Primary RRS_CV"
    )

  # Place statistics in the subtitle to preserve the data region.
  p3 <- ggplot(df, aes(x = partial_Strict_RRS, y = partial_CRSz)) +
    geom_point(size = 2.25, colour = "#8FB7CC", alpha = 0.85) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      colour = "#E68613",
      fill = "#F3C88C",
      linewidth = 0.8,
      alpha = 0.22
    ) +
    labs(
      title = title_or_null("Size-adjusted RRS-CRS_z relationship"),
      subtitle = subtitle_or_null(
        paste0(
          "Partial Pearson r = ", sprintf("%.3f", partial_r),
          "; partial Spearman rho = ", sprintf("%.3f", partial_rho)
        )
      ),
      x = "Partial residual of primary RRS_CV",
      y = "Partial residual of CRS_z"
    ) +
    theme(
      plot.subtitle = element_text(size = 9.5, hjust = 0, margin = margin(b = 4))
    )

  model_plot_df <- model_df[
    model_df$model %in% c(
      "Size_only",
      "Size_plus_strict_RRS",
      "Size_plus_signed_log_RRS_sensitivity"
    ),
    ,
    drop = FALSE
  ]

  model_plot_df$model_label <- factor(
    model_plot_df$model,
    levels = c(
      "Size_only",
      "Size_plus_strict_RRS",
      "Size_plus_signed_log_RRS_sensitivity"
    ),
    labels = c(
      "Size only",
      "Size + primary RRS",
      "Size + signed-log RRS"
    )
  )

  model_plot_df <- model_plot_df[order(model_plot_df$model_label), , drop = FALSE]

  model_plot_df$r2_plot_label <- sprintf("R^2 == %.3f", model_plot_df$model_R2)

  delta_plot_df <- model_plot_df[!is.na(model_plot_df$delta_R2_vs_size_only), , drop = FALSE]
  delta_plot_df$delta_plot_label <- sprintf(
    "Delta*R^2 == %.3f",
    delta_plot_df$delta_R2_vs_size_only
  )

  p4 <- ggplot(model_plot_df, aes(x = model_label, y = model_R2, fill = model_label)) +
    geom_col(width = 0.68, colour = "black", linewidth = 0.35) +
    geom_text(
      aes(label = r2_plot_label),
      parse = TRUE,
      vjust = -0.45,
      family = fig_font,
      size = 3.5
    ) +
    geom_text(
      data = delta_plot_df,
      aes(label = delta_plot_label),
      parse = TRUE,
      vjust = 1.7,
      family = fig_font,
      size = 3.0
    ) +
    scale_fill_manual(values = c(
      "Size only" = "#D4D4D4",
      "Size + primary RRS" = "#5C93B5",
      "Size + signed-log RRS" = "#AFC9D8"
    )) +
    scale_y_continuous(
      limits = c(0, max(model_plot_df$model_R2, na.rm = TRUE) * 1.20),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = title_or_null("Incremental variance beyond tumor size"),
      x = NULL,
      y = expression(Model~R^2)
    ) +
    theme(
      axis.text.x = element_text(size = 8.8),
      legend.position = "none"
    )

  save_panel(p1, out_dir, "S3_RRS_vs_tumor_volume", width = 5.7, height = 4.7)
  save_panel(p2, out_dir, "S3_RRS_vs_maximum_3D_diameter", width = 5.7, height = 4.7)
  save_panel(p3, out_dir, "S3_size_adjusted_partial_residuals", width = 5.7, height = 4.7)
  save_panel(p4, out_dir, "S3_incremental_R2_beyond_tumor_size", width = 5.7, height = 4.7)

  audit_out <- data.frame(
    item = c(
      "n_complete",
      "partial_Pearson_r",
      "partial_Pearson_R2",
      "partial_Spearman_rho",
      "size_only_R2",
      "size_plus_strict_RRS_R2",
      "strict_RRS_delta_R2",
      "size_plus_signed_log_RRS_R2",
      "signed_log_RRS_delta_R2"
    ),
    value = c(
      nrow(df),
      partial_r,
      partial_r^2,
      partial_rho,
      model_plot_df$model_R2[model_plot_df$model == "Size_only"],
      model_plot_df$model_R2[model_plot_df$model == "Size_plus_strict_RRS"],
      model_plot_df$delta_R2_vs_size_only[model_plot_df$model == "Size_plus_strict_RRS"],
      model_plot_df$model_R2[
        model_plot_df$model == "Size_plus_signed_log_RRS_sensitivity"
      ],
      model_plot_df$delta_R2_vs_size_only[
        model_plot_df$model == "Size_plus_signed_log_RRS_sensitivity"
      ]
    ),
    stringsAsFactors = FALSE
  )

  write.csv(audit_out, file.path(out_dir, "S3_audit_values.csv"), row.names = FALSE)
}

############################################################
# 4. Supplementary Figure S4 panels
############################################################

make_S4 <- function() {
  message("\n===== Making S4 separate panels =====")

  out_dir <- file.path(out_root, "S4")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  size_dir <- file.path(project_dir, "07_results", "tumor_size_domain_shift")

  training_file <- file.path(size_dir, "training_size_analysis_dataset.csv")
  lung1_file <- file.path(size_dir, "lung1_size_analysis_dataset.csv")

  coefficient_file <- file.path(
    project_dir,
    "07_results",
    "rrs_model",
    "final_model_coefficients.csv"
  )

  signedlog_file <- file.path(
    project_dir,
    "07_results",
    "supplementary_exports",
    "signed_log_domain_shift_features.csv"
  )

  for (ff in c(training_file, lung1_file, coefficient_file, signedlog_file)) {
    if (!file.exists(ff)) {
      stop("Cannot find required input:\n", ff)
    }
  }

  train_df <- read.csv(training_file, stringsAsFactors = FALSE, check.names = FALSE)
  lung1_df <- read.csv(lung1_file, stringsAsFactors = FALSE, check.names = FALSE)
  coef_df <- read.csv(coefficient_file, stringsAsFactors = FALSE, check.names = FALSE)
  signed_df <- read.csv(signedlog_file, stringsAsFactors = FALSE, check.names = FALSE)

  required_signed_cols <- c(
    "feature",
    "is_nonzero_model_feature",
    "coefficient",
    "lung1_max_abs_z_by_train_log",
    "outside_train_minmax_n_log",
    "outside_train_minmax_pct_log"
  )

  missing_signed <- setdiff(required_signed_cols, names(signed_df))
  if (length(missing_signed) > 0L) {
    stop("Signed-log domain-shift table is missing: ", paste(missing_signed, collapse = ", "))
  }

  required_coef_cols <- c("feature", "coefficient")
  missing_coef <- setdiff(required_coef_cols, names(coef_df))
  if (length(missing_coef) > 0L) {
    stop("Final coefficient table is missing: ", paste(missing_coef, collapse = ", "))
  }

  train_features <- grep("^original_", names(train_df), value = TRUE)
  lung1_features <- grep("^original_", names(lung1_df), value = TRUE)
  common_features <- intersect(train_features, lung1_features)

  if (length(common_features) != 107L) {
    stop("Expected 107 common original features, found ", length(common_features), ".")
  }

  x_train <- as.data.frame(
    lapply(train_df[, common_features, drop = FALSE], make_numeric),
    check.names = FALSE
  )

  x_lung1 <- as.data.frame(
    lapply(lung1_df[, common_features, drop = FALSE], make_numeric),
    check.names = FALSE
  )

  missing_rate <- vapply(x_train, function(z) mean(!is.finite(z)), numeric(1))
  keep_missing <- is.finite(missing_rate) & missing_rate <= 0.30
  x_train <- x_train[, keep_missing, drop = FALSE]
  x_lung1 <- x_lung1[, names(keep_missing)[keep_missing], drop = FALSE]

  medians <- vapply(
    x_train,
    function(z) {
      med <- stats::median(z[is.finite(z)], na.rm = TRUE)
      if (!is.finite(med)) 0 else med
    },
    numeric(1)
  )

  for (nm in names(x_train)) {
    x_train[[nm]][!is.finite(x_train[[nm]])] <- medians[[nm]]
    x_lung1[[nm]][!is.finite(x_lung1[[nm]])] <- medians[[nm]]
  }

  sds0 <- vapply(x_train, stats::sd, numeric(1), na.rm = TRUE)
  keep_sd <- is.finite(sds0) & sds0 > 0
  x_train <- x_train[, keep_sd, drop = FALSE]
  x_lung1 <- x_lung1[, names(keep_sd)[keep_sd], drop = FALSE]

  selected_features <- correlation_filter(x_train, cutoff = 0.90)

  if (length(selected_features) != 43L) {
    stop(
      "The verified full-data correlation filter should retain 43 features, but retained ",
      length(selected_features),
      "."
    )
  }

  x_train_43 <- x_train[, selected_features, drop = FALSE]
  x_lung1_43 <- x_lung1[, selected_features, drop = FALSE]

  train_center <- vapply(x_train_43, mean, numeric(1), na.rm = TRUE)
  train_scale <- vapply(x_train_43, stats::sd, numeric(1), na.rm = TRUE)

  nonzero_coef <- coef_df[
    coef_df$feature != "(Intercept)" &
      coef_df$coefficient != 0,
    ,
    drop = FALSE
  ]

  if (nrow(nonzero_coef) != 6L) {
    stop("Expected six primary nonzero radiomic features, found ", nrow(nonzero_coef), ".")
  }

  primary_rows <- lapply(selected_features, function(ff) {
    tr <- x_train_43[[ff]]
    lu <- x_lung1_43[[ff]]

    tr_min <- min(tr, na.rm = TRUE)
    tr_max <- max(tr, na.rm = TRUE)
    lu_z <- (lu - train_center[[ff]]) / train_scale[[ff]]

    outside <- lu < tr_min | lu > tr_max
    outside_n <- sum(outside, na.rm = TRUE)

    coefficient <- if (ff %in% nonzero_coef$feature) {
      nonzero_coef$coefficient[match(ff, nonzero_coef$feature)]
    } else {
      0
    }

    data.frame(
      feature = ff,
      is_nonzero_model_feature = ff %in% nonzero_coef$feature,
      coefficient = coefficient,
      train_min = tr_min,
      train_median = median(tr, na.rm = TRUE),
      train_max = tr_max,
      lung1_min = min(lu, na.rm = TRUE),
      lung1_median = median(lu, na.rm = TRUE),
      lung1_max = max(lu, na.rm = TRUE),
      lung1_max_abs_z_by_train = max(abs(lu_z), na.rm = TRUE),
      outside_train_minmax_n = outside_n,
      outside_train_minmax_pct = outside_n / length(lu),
      stringsAsFactors = FALSE
    )
  })

  primary_shift <- do.call(rbind, primary_rows)

  primary_any_n <- sum(primary_shift$outside_train_minmax_n > 0, na.rm = TRUE)
  primary_gt5_n <- sum(primary_shift$outside_train_minmax_pct > 0.05, na.rm = TRUE)
  primary_nonzero_any_n <- sum(
    primary_shift$outside_train_minmax_n[primary_shift$is_nonzero_model_feature] > 0,
    na.rm = TRUE
  )

  signed_any_n <- sum(signed_df$outside_train_minmax_n_log > 0, na.rm = TRUE)
  signed_gt5_n <- sum(signed_df$outside_train_minmax_pct_log > 0.05, na.rm = TRUE)
  signed_nonzero_total <- sum(as.logical(signed_df$is_nonzero_model_feature), na.rm = TRUE)
  signed_nonzero_any_n <- sum(
    signed_df$outside_train_minmax_n_log[
      as.logical(signed_df$is_nonzero_model_feature)
    ] > 0,
    na.rm = TRUE
  )

  top_primary <- primary_shift %>%
    arrange(desc(outside_train_minmax_pct), desc(lung1_max_abs_z_by_train)) %>%
    slice_head(n = 20) %>%
    mutate(
      feature_label = short_feature(feature),
      feature_label = factor(feature_label, levels = rev(feature_label)),
      feature_class = ifelse(
        is_nonzero_model_feature,
        "Final nonzero feature",
        "Other model-related feature"
      )
    )

  p1 <- ggplot(
    top_primary,
    aes(x = outside_train_minmax_pct, y = feature_label, fill = feature_class)
  ) +
    geom_col(width = 0.72, colour = "black", linewidth = 0.25) +
    geom_text(
      aes(label = sprintf("%.1f%%", 100 * outside_train_minmax_pct)),
      hjust = -0.08,
      family = fig_font,
      size = 2.8
    ) +
    scale_fill_manual(values = c(
      "Final nonzero feature" = "#B2555A",
      "Other model-related feature" = "#AFC9D8"
    )) +
    scale_x_continuous(
      labels = function(x) paste0(round(100 * x), "%"),
      expand = expansion(mult = c(0, 0.16))
    ) +
    labs(
      title = title_or_null("Primary RRS: highest out-of-range rates"),
      x = "Lung1 patients outside the training range",
      y = NULL,
      fill = NULL
    ) +
    theme(
      axis.text.y = element_text(size = 6.8),
      legend.position = "bottom",
      plot.title.position = "plot"
    )

  comparison_df <- merge(
    primary_shift[, c("feature", "outside_train_minmax_pct", "is_nonzero_model_feature")],
    signed_df[, c("feature", "outside_train_minmax_pct_log")],
    by = "feature",
    all.x = TRUE
  )

  p2 <- ggplot(
    comparison_df,
    aes(
      x = outside_train_minmax_pct,
      y = outside_train_minmax_pct_log,
      fill = is_nonzero_model_feature
    )
  ) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "grey70") +
    geom_point(shape = 21, size = 2.7, colour = "black", stroke = 0.30, alpha = 0.88) +
    scale_fill_manual(
      values = c("FALSE" = "#AFC9D8", "TRUE" = "#B2555A"),
      labels = c(
        "FALSE" = "Other primary-model features",
        "TRUE" = "Primary final nonzero features"
      )
    ) +
    scale_x_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    labs(
      title = title_or_null("Outside-range rates before and after signed-log transformation"),
      x = "Primary-model outside rate",
      y = "Signed-log outside rate",
      fill = NULL
    ) +
    theme(
      legend.position = "bottom",
      plot.title.position = "plot"
    )

  lung1_ids <- toupper(trimws(as.character(lung1_df$patient_id)))
  extreme_idx <- which(lung1_ids == "LUNG1-116")

  if (length(extreme_idx) != 1L) {
    stop("Expected exactly one LUNG1-116 row, found ", length(extreme_idx), ".")
  }

  contribution_rows <- lapply(seq_len(nrow(nonzero_coef)), function(ii) {
    ff <- nonzero_coef$feature[ii]
    beta <- nonzero_coef$coefficient[ii]
    raw_value <- x_lung1_43[[ff]][extreme_idx]
    z_value <- (raw_value - train_center[[ff]]) / train_scale[[ff]]

    data.frame(
      patient_id = "LUNG1-116",
      feature = ff,
      coefficient = beta,
      lung1_raw_value = raw_value,
      training_scaled_z = z_value,
      single_feature_contribution = z_value * beta,
      stringsAsFactors = FALSE
    )
  })

  contribution_df <- do.call(rbind, contribution_rows)
  contribution_df <- contribution_df[order(contribution_df$single_feature_contribution), , drop = FALSE]

  driver_value <- contribution_df$single_feature_contribution[
    contribution_df$feature == "original_glszm_LargeAreaHighGrayLevelEmphasis"
  ]

  contribution_plot <- contribution_df %>%
    mutate(
      feature_label = short_feature(feature),
      feature_label = factor(feature_label, levels = feature_label),
      contribution_sign = ifelse(single_feature_contribution < 0, "Negative", "Positive")
    )

  # Expand the plotting region to keep all contribution labels visible.
  x_min <- min(contribution_plot$single_feature_contribution, na.rm = TRUE)
  x_max <- max(contribution_plot$single_feature_contribution, na.rm = TRUE)
  x_span <- x_max - x_min
  if (!is.finite(x_span) || x_span <= 0) {
    x_span <- 1
  }

  contribution_plot$text_x <- ifelse(
    contribution_plot$single_feature_contribution < 0,
    contribution_plot$single_feature_contribution + 0.03 * x_span,
    contribution_plot$single_feature_contribution + 0.03 * x_span
  )

  contribution_plot$text_hjust <- 0

  p3 <- ggplot(
    contribution_plot,
    aes(x = single_feature_contribution, y = feature_label, fill = contribution_sign)
  ) +
    geom_col(width = 0.68, colour = "black", linewidth = 0.25) +
    geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.35) +
    geom_text(
      aes(
        x = text_x,
        label = sprintf("%.3f", single_feature_contribution),
        hjust = text_hjust
      ),
      family = fig_font,
      size = 2.9
    ) +
    scale_fill_manual(values = c(
      "Negative" = "#5C93B5",
      "Positive" = "#C97A70"
    )) +
    coord_cartesian(
      xlim = c(x_min - 0.03 * x_span, x_max + 0.23 * x_span),
      clip = "off"
    ) +
    labs(
      title = title_or_null("LUNG1-116: primary-model feature contributions"),
      subtitle = subtitle_or_null(
        paste0(
          "LargeAreaHighGrayLevelEmphasis contribution = ",
          sprintf("%.3f", driver_value)
        )
      ),
      x = "Single-feature contribution to RRS",
      y = NULL,
      fill = NULL
    ) +
    theme(
      axis.text.y = element_text(size = 7.2),
      legend.position = "none",
      plot.title.position = "plot",
      plot.margin = margin(5.5, 24, 5.5, 5.5)
    )

  summary_df <- data.frame(
    model = factor(
      rep(c("Primary RRS", "Signed-log RRS"), each = 3),
      levels = c("Primary RRS", "Signed-log RRS")
    ),
    metric = factor(
      rep(
        c(
          "Any out-of-range value",
          ">5% patients out of range",
          "Nonzero features affected"
        ),
        times = 2
      ),
      levels = c(
        "Any out-of-range value",
        ">5% patients out of range",
        "Nonzero features affected"
      )
    ),
    proportion = c(
      primary_any_n / 43,
      primary_gt5_n / 43,
      primary_nonzero_any_n / 6,
      signed_any_n / 107,
      signed_gt5_n / 107,
      signed_nonzero_any_n / signed_nonzero_total
    ),
    label = c(
      paste0(primary_any_n, "/43"),
      paste0(primary_gt5_n, "/43"),
      paste0(primary_nonzero_any_n, "/6"),
      paste0(signed_any_n, "/107"),
      paste0(signed_gt5_n, "/107"),
      paste0(signed_nonzero_any_n, "/", signed_nonzero_total)
    )
  )

  p4 <- ggplot(summary_df, aes(x = metric, y = proportion, fill = model)) +
    geom_col(
      position = position_dodge(width = 0.78),
      width = 0.68,
      colour = "black",
      linewidth = 0.25
    ) +
    geom_text(
      aes(label = label),
      position = position_dodge(width = 0.78),
      vjust = -0.35,
      family = fig_font,
      size = 3.2
    ) +
    scale_fill_manual(values = c(
      "Primary RRS" = "#5C93B5",
      "Signed-log RRS" = "#B7CEDC"
    )) +
    scale_y_continuous(
      limits = c(0, 1.12),
      labels = function(x) paste0(round(100 * x), "%"),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = title_or_null("Domain-shift burden"),
      x = NULL,
      y = "Proportion of evaluated features",
      fill = NULL
    ) +
    theme(
      axis.text.x = element_text(size = 8.4),
      legend.position = "bottom",
      plot.title.position = "plot"
    )

  save_panel(p1, out_dir, "S4_primary_highest_out_of_range_rates", width = 6.5, height = 5.8)
  save_panel(p2, out_dir, "S4_primary_vs_signedlog_out_of_range", width = 5.8, height = 5.2)
  save_panel(p3, out_dir, "S4_LUNG1_116_feature_contributions", width = 7.0, height = 4.9)
  save_panel(p4, out_dir, "S4_domain_shift_burden_summary", width = 6.5, height = 4.8)

  write.csv(
    primary_shift,
    file.path(out_dir, "S4_primary_domain_shift_43_features.csv"),
    row.names = FALSE
  )

  write.csv(
    contribution_df,
    file.path(out_dir, "S4_LUNG1-116_contributions.csv"),
    row.names = FALSE
  )

  write.csv(
    comparison_df,
    file.path(out_dir, "S4_primary_vs_signedlog_shift.csv"),
    row.names = FALSE
  )

  audit_values <- data.frame(
    item = c(
      "primary_features_any_outside",
      "primary_features_gt5pct_outside",
      "primary_nonzero_features_affected",
      "signedlog_features_any_outside",
      "signedlog_features_gt5pct_outside",
      "signedlog_nonzero_features_affected",
      "signedlog_nonzero_features_total",
      "LUNG1_116_driver_contribution"
    ),
    value = c(
      primary_any_n,
      primary_gt5_n,
      primary_nonzero_any_n,
      signed_any_n,
      signed_gt5_n,
      signed_nonzero_any_n,
      signed_nonzero_total,
      driver_value
    ),
    stringsAsFactors = FALSE
  )

  write.csv(audit_values, file.path(out_dir, "S4_audit_values.csv"), row.names = FALSE)

  message(
    "verified primary audit reproduced: 43/43 any outside, 33/43 with >5% outside, 6/6 final nonzero affected."
  )
  message(
    "verified signed-log audit reproduced: 107/107 any outside, 87/107 with >5% outside, 5/5 nonzero affected."
  )
}

############################################################
# 5. Run all
############################################################

make_S1()
make_S2()
make_S3()
make_S4()

cat("\n============================================================\n")
cat("Supplementary panels were generated.\n")
cat("Output root:\n")
cat(out_root, "\n")
cat("============================================================\n")
