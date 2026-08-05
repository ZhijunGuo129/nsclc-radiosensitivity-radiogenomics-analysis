# 07_summarize_training_cohort.R
#
# Summarize the RRS training cohort and model outputs.
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

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Missing R package 'openxlsx'. Install dependencies with environment/install_r_dependencies.R.")
}

library(openxlsx)

############################################################
# 2. Paths
############################################################

merged_path <- "07_results/radiomics_features/radiomics_molecular_merged.rds"
pred_path <- "07_results/rrs_model/patient_mean_predictions.csv"

out_dir <- "07_results/rrs_model"
fig_dir <- "08_figures_final/rrs_training_summary"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(merged_path)) {
  stop("Cannot find merged table: ", merged_path)
}

if (!file.exists(pred_path)) {
  stop("Cannot find the RRS patient prediction table: ", pred_path)
}

############################################################
# 3. Load data
############################################################

mol_rad <- readRDS(merged_path)

pred_df <- read.csv(
  pred_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\n===== Loaded data =====\n")
cat("Merged molecular-radiomics table dim:\n")
print(dim(mol_rad))

cat("\nPrediction table dim:\n")
print(dim(pred_df))

cat("\nPrediction columns:\n")
print(colnames(pred_df))

############################################################
# 4. Merge prediction with molecular scores
############################################################

mol_rad$patient_id <- toupper(trimws(mol_rad$patient_id))
pred_df$patient_id <- toupper(trimws(pred_df$patient_id))

analysis_df <- merge(
  mol_rad[, c(
    "patient_id",
    "CRS",
    "CRS_z",
    "RSI",
    "RSI_z",
    "Hypoxia_core10_score",
    "Hypoxia_core10_z",
    "CRS_group_median",
    "RSI_group_median",
    "Hypoxia_group_median"
  )],
  pred_df,
  by = "patient_id",
  all = FALSE
)

cat("\n===== Analysis table =====\n")
print(dim(analysis_df))

cat("\nColumns:\n")
print(colnames(analysis_df))

############################################################
# 5. Define RRS
############################################################

analysis_df$RRS_CV <- as.numeric(analysis_df$predicted_CRS_z_mean)
analysis_df$RRS_CV_sd <- as.numeric(analysis_df$predicted_CRS_z_sd)

analysis_df$observed_CRS_z <- as.numeric(analysis_df$observed_CRS_z)
analysis_df$CRS_z <- as.numeric(analysis_df$CRS_z)

# observed_CRS_z should equal CRS_z
analysis_df$CRS_z_diff_check <- analysis_df$observed_CRS_z - analysis_df$CRS_z

cat("\n===== CRS_z consistency check =====\n")
print(summary(analysis_df$CRS_z_diff_check))

############################################################
# 6. Utility functions
############################################################

calc_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  
  data.frame(
    n = sum(ok),
    Pearson = suppressWarnings(cor(x[ok], y[ok], method = "pearson")),
    Pearson_p = suppressWarnings(cor.test(x[ok], y[ok], method = "pearson")$p.value),
    Spearman = suppressWarnings(cor(x[ok], y[ok], method = "spearman")),
    Spearman_p = suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)$p.value),
    stringsAsFactors = FALSE
  )
}

manual_auc <- function(labels, scores) {
  ok <- is.finite(labels) & is.finite(scores)
  labels <- labels[ok]
  scores <- scores[ok]
  
  labels <- as.numeric(labels)
  
  if (!all(labels %in% c(0, 1))) {
    stop("Labels must be 0/1.")
  }
  
  n_pos <- sum(labels == 1)
  n_neg <- sum(labels == 0)
  
  if (n_pos == 0 || n_neg == 0) {
    return(NA_real_)
  }
  
  ranks <- rank(scores, ties.method = "average")
  auc <- (sum(ranks[labels == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
  
  auc
}

make_roc_curve <- function(labels, scores) {
  ok <- is.finite(labels) & is.finite(scores)
  labels <- as.numeric(labels[ok])
  scores <- scores[ok]
  
  ord <- order(scores, decreasing = TRUE)
  labels <- labels[ord]
  scores <- scores[ord]
  
  thresholds <- c(Inf, unique(scores), -Inf)
  
  out <- data.frame()
  
  for (thr in thresholds) {
    pred <- ifelse(scores >= thr, 1, 0)
    
    tp <- sum(pred == 1 & labels == 1)
    fp <- sum(pred == 1 & labels == 0)
    tn <- sum(pred == 0 & labels == 0)
    fn <- sum(pred == 0 & labels == 1)
    
    tpr <- ifelse((tp + fn) > 0, tp / (tp + fn), NA)
    fpr <- ifelse((fp + tn) > 0, fp / (fp + tn), NA)
    
    out <- rbind(
      out,
      data.frame(
        threshold = thr,
        TPR = tpr,
        FPR = fpr,
        TP = tp,
        FP = fp,
        TN = tn,
        FN = fn
      )
    )
  }
  
  out <- out[order(out$FPR, out$TPR), ]
  out
}

calc_binary_metrics <- function(labels, scores, threshold) {
  labels <- as.numeric(labels)
  pred <- ifelse(scores >= threshold, 1, 0)
  
  tp <- sum(pred == 1 & labels == 1)
  fp <- sum(pred == 1 & labels == 0)
  tn <- sum(pred == 0 & labels == 0)
  fn <- sum(pred == 0 & labels == 1)
  
  sensitivity <- ifelse((tp + fn) > 0, tp / (tp + fn), NA)
  specificity <- ifelse((tn + fp) > 0, tn / (tn + fp), NA)
  accuracy <- (tp + tn) / length(labels)
  balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
  
  data.frame(
    threshold = threshold,
    TP = tp,
    FP = fp,
    TN = tn,
    FN = fn,
    sensitivity = sensitivity,
    specificity = specificity,
    accuracy = accuracy,
    balanced_accuracy = balanced_accuracy,
    stringsAsFactors = FALSE
  )
}

bootstrap_auc_ci <- function(labels, scores, n_boot = 1000, seed = 20260724) {
  set.seed(seed)
  
  n <- length(labels)
  aucs <- rep(NA_real_, n_boot)
  
  for (i in seq_len(n_boot)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    
    labs_i <- labels[idx]
    scores_i <- scores[idx]
    
    if (length(unique(labs_i)) < 2) {
      next
    }
    
    aucs[i] <- manual_auc(labs_i, scores_i)
  }
  
  data.frame(
    AUC_boot_mean = mean(aucs, na.rm = TRUE),
    AUC_boot_sd = sd(aucs, na.rm = TRUE),
    AUC_CI_lower = as.numeric(quantile(aucs, 0.025, na.rm = TRUE)),
    AUC_CI_upper = as.numeric(quantile(aucs, 0.975, na.rm = TRUE)),
    n_boot_valid = sum(is.finite(aucs)),
    stringsAsFactors = FALSE
  )
}

############################################################
# 7. Correlation: RRS_CV vs molecular scores
############################################################

cor_CRS <- calc_cor(analysis_df$CRS_z, analysis_df$RRS_CV)
cor_RSI <- calc_cor(analysis_df$RSI_z, analysis_df$RRS_CV)
cor_Hypoxia <- calc_cor(analysis_df$Hypoxia_core10_z, analysis_df$RRS_CV)

cor_table <- rbind(
  data.frame(comparison = "RRS_CV_vs_CRS_z", cor_CRS),
  data.frame(comparison = "RRS_CV_vs_RSI_z", cor_RSI),
  data.frame(comparison = "RRS_CV_vs_Hypoxia_core10_z", cor_Hypoxia)
)

cat("\n===== Correlation table =====\n")
print(cor_table)

############################################################
# 8. CRS high/low discrimination
############################################################

crs_median <- median(analysis_df$CRS_z, na.rm = TRUE)

analysis_df$CRS_high_binary <- ifelse(analysis_df$CRS_z >= crs_median, 1, 0)
analysis_df$CRS_group_by_CRSz_median <- ifelse(
  analysis_df$CRS_high_binary == 1,
  "CRS_high",
  "CRS_low"
)

rrs_median <- median(analysis_df$RRS_CV, na.rm = TRUE)
analysis_df$RRS_group_by_RRS_median <- ifelse(
  analysis_df$RRS_CV >= rrs_median,
  "RRS_high",
  "RRS_low"
)

cat("\n===== Group counts =====\n")
print(table(analysis_df$CRS_group_by_CRSz_median))
print(table(analysis_df$RRS_group_by_RRS_median))

wilcox_rrs_by_crs_group <- wilcox.test(
  RRS_CV ~ CRS_group_by_CRSz_median,
  data = analysis_df
)

cat("\n===== Wilcoxon test: RRS_CV by observed CRS group =====\n")
print(wilcox_rrs_by_crs_group)

auc_CRS_high <- manual_auc(
  labels = analysis_df$CRS_high_binary,
  scores = analysis_df$RRS_CV
)

auc_ci <- bootstrap_auc_ci(
  labels = analysis_df$CRS_high_binary,
  scores = analysis_df$RRS_CV,
  n_boot = 1000,
  seed = 20260724
)

roc_df <- make_roc_curve(
  labels = analysis_df$CRS_high_binary,
  scores = analysis_df$RRS_CV
)

binary_metrics_median <- calc_binary_metrics(
  labels = analysis_df$CRS_high_binary,
  scores = analysis_df$RRS_CV,
  threshold = rrs_median
)

cat("\n===== AUC for identifying CRS-high patients =====\n")
print(auc_CRS_high)
print(auc_ci)

cat("\n===== Binary metrics using RRS median cutoff =====\n")
print(binary_metrics_median)

cat("\n===== Confusion table: observed CRS group vs RRS group =====\n")
print(table(
  Observed_CRS = analysis_df$CRS_group_by_CRSz_median,
  Predicted_RRS = analysis_df$RRS_group_by_RRS_median
))

############################################################
# 9. RRS summary table
############################################################

rrs_summary <- data.frame(
  item = c(
    "n_patients",
    "CRS_z_median_cutoff",
    "RRS_CV_median_cutoff",
    "RRS_vs_CRSz_Pearson",
    "RRS_vs_CRSz_Pearson_p",
    "RRS_vs_CRSz_Spearman",
    "RRS_vs_CRSz_Spearman_p",
    "AUC_for_CRS_high",
    "AUC_boot_CI_lower",
    "AUC_boot_CI_upper",
    "Wilcoxon_RRS_by_CRS_group_p",
    "median_RRS_in_CRS_low",
    "median_RRS_in_CRS_high",
    "sensitivity_RRS_median_cutoff",
    "specificity_RRS_median_cutoff",
    "accuracy_RRS_median_cutoff",
    "balanced_accuracy_RRS_median_cutoff"
  ),
  value = c(
    nrow(analysis_df),
    crs_median,
    rrs_median,
    cor_CRS$Pearson,
    cor_CRS$Pearson_p,
    cor_CRS$Spearman,
    cor_CRS$Spearman_p,
    auc_CRS_high,
    auc_ci$AUC_CI_lower,
    auc_ci$AUC_CI_upper,
    wilcox_rrs_by_crs_group$p.value,
    median(analysis_df$RRS_CV[analysis_df$CRS_high_binary == 0], na.rm = TRUE),
    median(analysis_df$RRS_CV[analysis_df$CRS_high_binary == 1], na.rm = TRUE),
    binary_metrics_median$sensitivity,
    binary_metrics_median$specificity,
    binary_metrics_median$accuracy,
    binary_metrics_median$balanced_accuracy
  ),
  stringsAsFactors = FALSE
)

cat("\n===== RRS summary =====\n")
print(rrs_summary)

############################################################
# 10. Save outputs
############################################################

write.csv(
  analysis_df,
  file = file.path(out_dir, "rrs_training_cohort_analysis.csv"),
  row.names = FALSE
)

write.csv(
  cor_table,
  file = file.path(out_dir, "rrs_correlation_summary.csv"),
  row.names = FALSE
)

write.csv(
  roc_df,
  file = file.path(out_dir, "rrs_roc_curve.csv"),
  row.names = FALSE
)

write.csv(
  binary_metrics_median,
  file = file.path(out_dir, "rrs_binary_metrics.csv"),
  row.names = FALSE
)

write.csv(
  rrs_summary,
  file = file.path(out_dir, "rrs_training_cohort_summary.csv"),
  row.names = FALSE
)

saveRDS(
  analysis_df,
  file = file.path(out_dir, "rrs_training_cohort_analysis.rds")
)

wb <- createWorkbook()

addWorksheet(wb, "RRS_summary")
writeData(wb, "RRS_summary", rrs_summary)

addWorksheet(wb, "correlation_table")
writeData(wb, "correlation_table", cor_table)

addWorksheet(wb, "analysis_table")
writeData(wb, "analysis_table", analysis_df)

addWorksheet(wb, "ROC_curve")
writeData(wb, "ROC_curve", roc_df)

addWorksheet(wb, "binary_metrics")
writeData(wb, "binary_metrics", binary_metrics_median)


saveWorkbook(
  wb,
  file = file.path(out_dir, "rrs_training_cohort_summary.xlsx"),
  overwrite = TRUE
)

############################################################
# 11. Figures
############################################################

png(
  filename = file.path(fig_dir, "rrs_vs_crs_scatter.png"),
  width = 1600,
  height = 1400,
  res = 200
)

plot(
  analysis_df$CRS_z,
  analysis_df$RRS_CV,
  pch = 19,
  xlab = "Observed CRS_z",
  ylab = "RRS_CV",
  main = "Cross-validated radiomic score vs CRS_z"
)

abline(lm(RRS_CV ~ CRS_z, data = analysis_df), lwd = 2)

legend(
  "topleft",
  legend = paste0(
    "Spearman = ",
    round(cor_CRS$Spearman, 3),
    "\nPearson = ",
    round(cor_CRS$Pearson, 3)
  ),
  bty = "n"
)

dev.off()

png(
  filename = file.path(fig_dir, "rrs_by_crs_group.png"),
  width = 1400,
  height = 1200,
  res = 200
)

boxplot(
  RRS_CV ~ CRS_group_by_CRSz_median,
  data = analysis_df,
  xlab = "Observed CRS group",
  ylab = "RRS_CV",
  main = "RRS_CV distribution by CRS high/low group"
)

stripchart(
  RRS_CV ~ CRS_group_by_CRSz_median,
  data = analysis_df,
  vertical = TRUE,
  method = "jitter",
  pch = 19,
  add = TRUE
)

legend(
  "topleft",
  legend = paste0(
    "Wilcoxon p = ",
    signif(wilcox_rrs_by_crs_group$p.value, 3)
  ),
  bty = "n"
)

dev.off()

png(
  filename = file.path(fig_dir, "rrs_roc_curve.png"),
  width = 1400,
  height = 1200,
  res = 200
)

plot(
  roc_df$FPR,
  roc_df$TPR,
  type = "l",
  lwd = 2,
  xlab = "False positive rate",
  ylab = "True positive rate",
  main = "RRS_CV for identifying CRS-high tumors",
  xlim = c(0, 1),
  ylim = c(0, 1)
)

abline(0, 1, lty = 2)

legend(
  "bottomright",
  legend = paste0(
    "AUC = ",
    round(auc_CRS_high, 3),
    "\n95% CI = ",
    round(auc_ci$AUC_CI_lower, 3),
    "-",
    round(auc_ci$AUC_CI_upper, 3)
  ),
  bty = "n"
)

dev.off()

cat("\n===== RRS training-cohort analysis completed =====\n")
cat("Main output files:\n")
cat("1. 07_results/rrs_model/rrs_training_cohort_summary.xlsx\n")
cat("2. 07_results/rrs_model/rrs_training_cohort_summary.csv\n")
cat("3. 07_results/rrs_model/rrs_training_cohort_analysis.csv\n")
cat("4. 08_figures_final/rrs_training_summary/rrs_vs_crs_scatter.png\n")
cat("5. 08_figures_final/rrs_training_summary/rrs_by_crs_group.png\n")
cat("6. 08_figures_final/rrs_training_summary/rrs_roc_curve.png\n")

cat("1. Correlation table\n")
cat("2. Wilcoxon test\n")
cat("3. AUC for identifying CRS-high patients\n")
cat("4. Binary metrics using RRS median cutoff\n")
cat("5. RRS summary\n")
