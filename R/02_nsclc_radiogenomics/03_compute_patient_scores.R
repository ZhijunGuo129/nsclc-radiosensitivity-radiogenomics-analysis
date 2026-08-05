# 03_compute_patient_scores.R
#
# Compute patient-level CRS and predefined molecular-module scores.
#
# Run from an analysis workspace configured through environment variables.

options(stringsAsFactors = FALSE)
options(timeout = 1800)


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
# 2. Folders
############################################################

dirs <- c(
  "05_molecular_scores",
  "07_results",
  file.path("08_figures_final", "patient_molecular_scores")
)

for (d in dirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("\n===== Current working directory =====\n")
print(getwd())

############################################################
# 3. Load patient expression matrix from the transcriptome audit
############################################################

patient_expr_path <- "01_raw_data/NSCLC_Radiogenomics/GSE103584_NSCLC_RNAseq_SYMBOL_gene_by_sample_RAW_AUDITED.rds"

if (!file.exists(patient_expr_path)) {
  stop(
    "Cannot find patient expression matrix:\n",
    patient_expr_path,
    "\nRun the transcriptome audit step first."
  )
}

expr_raw <- readRDS(patient_expr_path)

cat("\n===== Raw patient expression loaded =====\n")
cat("Raw expression dim, genes x samples:\n")
print(dim(expr_raw))

############################################################
# 4. Basic preprocessing
############################################################

finite_raw <- as.numeric(expr_raw[is.finite(expr_raw)])

raw_summary <- data.frame(
  min = min(finite_raw, na.rm = TRUE),
  q25 = as.numeric(quantile(finite_raw, 0.25, na.rm = TRUE)),
  median = median(finite_raw, na.rm = TRUE),
  mean = mean(finite_raw, na.rm = TRUE),
  q75 = as.numeric(quantile(finite_raw, 0.75, na.rm = TRUE)),
  q95 = as.numeric(quantile(finite_raw, 0.95, na.rm = TRUE)),
  q99 = as.numeric(quantile(finite_raw, 0.99, na.rm = TRUE)),
  max = max(finite_raw, na.rm = TRUE),
  stringsAsFactors = FALSE
)

cat("\n===== Raw expression value summary =====\n")
print(raw_summary)

# In this GEO file, NA represents undetected / missing expression values.
# For scoring, we set NA to 0, then apply log2(x + 1) because q99 is high.
expr_proc <- expr_raw
expr_proc[is.na(expr_proc)] <- 0

patient_q99 <- raw_summary$q99

if (patient_q99 > 50) {
  expr_proc <- log2(expr_proc + 1)
  preprocess_note <- "NA set to 0, then log2(x + 1) applied"
} else {
  preprocess_note <- "NA set to 0, no log transform applied"
}

finite_proc <- as.numeric(expr_proc[is.finite(expr_proc)])

proc_summary <- data.frame(
  min = min(finite_proc, na.rm = TRUE),
  q25 = as.numeric(quantile(finite_proc, 0.25, na.rm = TRUE)),
  median = median(finite_proc, na.rm = TRUE),
  mean = mean(finite_proc, na.rm = TRUE),
  q75 = as.numeric(quantile(finite_proc, 0.75, na.rm = TRUE)),
  q95 = as.numeric(quantile(finite_proc, 0.95, na.rm = TRUE)),
  q99 = as.numeric(quantile(finite_proc, 0.99, na.rm = TRUE)),
  max = max(finite_proc, na.rm = TRUE),
  preprocess_note = preprocess_note,
  stringsAsFactors = FALSE
)

cat("\n===== Processed expression value summary =====\n")
print(proc_summary)

saveRDS(
  expr_proc,
  file = "01_raw_data/NSCLC_Radiogenomics/GSE103584_NSCLC_RNAseq_SYMBOL_gene_by_sample_LOG2_PROCESSED.rds"
)

############################################################
# 5. Load CRS symbol coefficient model from the CRS feature-mapping step
############################################################

crs_symbol_model_path <- "06_models/crs_gene_symbol_coefficients.rds"

if (!file.exists(crs_symbol_model_path)) {
  stop(
    "Cannot find the CRS gene-symbol coefficient model:\n",
    crs_symbol_model_path,
    "\nRun the CRS feature-mapping step first."
  )
}

crs_symbol_model <- readRDS(crs_symbol_model_path)

crs_intercept <- as.numeric(crs_symbol_model$intercept)
crs_coef <- crs_symbol_model$symbol_coefficients

cat("\n===== CRS symbol coefficient model loaded =====\n")
cat("CRS intercept:\n")
print(crs_intercept)

cat("\nCRS available symbol coefficient number:\n")
print(nrow(crs_coef))

cat("\nCRS coefficient preview:\n")
print(head(crs_coef, 30))

############################################################
# 6. Calculate CRS
############################################################

crs_genes <- crs_coef$symbol
missing_crs_genes <- setdiff(crs_genes, rownames(expr_proc))

cat("\n===== CRS gene availability check =====\n")
cat("CRS genes required:\n")
print(length(crs_genes))

cat("CRS genes missing in processed expression:\n")
print(length(missing_crs_genes))

if (length(missing_crs_genes) > 0) {
  stop(
    "Some CRS genes from the CRS feature-mapping step coefficient table are missing in expr_proc.\n",
    "Inspect missing_crs_genes before continuing."
  )
}

X_crs <- expr_proc[crs_genes, , drop = FALSE]

# Ensure coefficient order matches expression matrix row order
coef_vec <- crs_coef$coefficient[match(rownames(X_crs), crs_coef$symbol)]

CRS <- as.numeric(crs_intercept + colSums(X_crs * coef_vec))

names(CRS) <- colnames(expr_proc)

CRS_z <- as.numeric(scale(CRS))
names(CRS_z) <- names(CRS)

cat("\n===== CRS score summary =====\n")
print(summary(CRS))

cat("\n===== CRS z-score summary =====\n")
print(summary(CRS_z))

############################################################
# 7. Calculate classic RSI
############################################################

rsi_coef <- c(
  AR    = -0.0098009,
  JUN   =  0.0128283,
  STAT1 =  0.0254552,
  PRKCB = -0.0017589,
  RELA  = -0.0038171,
  ABL1  =  0.1070213,
  SUMO1 = -0.0002509,
  PAK2  = -0.0092431,
  HDAC1 = -0.0204469,
  IRF1  = -0.0441683
)

rsi_genes <- names(rsi_coef)

missing_rsi_genes <- setdiff(rsi_genes, rownames(expr_proc))

cat("\n===== RSI gene availability check =====\n")
cat("RSI genes required:\n")
print(rsi_genes)

cat("Missing RSI genes:\n")
print(missing_rsi_genes)

if (length(missing_rsi_genes) > 0) {
  stop(
    "Missing RSI genes. Cannot calculate RSI.\n",
    "Inspect the missing RSI genes before continuing."
  )
}

rsi_expr <- expr_proc[rsi_genes, , drop = FALSE]

# Published RSI is a rank-based linear algorithm.
# For each patient, rank the 10 RSI genes within that patient:
# lowest expression = 1, highest expression = 10.
rsi_rank <- apply(rsi_expr, 2, function(x) {
  rank(x, ties.method = "average")
})

RSI <- as.numeric(colSums(rsi_rank * rsi_coef[rownames(rsi_rank)]))
names(RSI) <- colnames(expr_proc)

RSI_z <- as.numeric(scale(RSI))
names(RSI_z) <- names(RSI)

cat("\n===== RSI score summary =====\n")
print(summary(RSI))

cat("\n===== RSI z-score summary =====\n")
print(summary(RSI_z))

############################################################
# 8. Calculate exploratory hypoxia core-gene score
############################################################

hypoxia_genes_core <- c(
  "HIF1A", "VEGFA", "CA9", "SLC2A1", "LDHA",
  "PDK1", "BNIP3", "NDRG1", "EGLN1", "EGLN3"
)

missing_hypoxia_genes <- setdiff(hypoxia_genes_core, rownames(expr_proc))

cat("\n===== Hypoxia gene availability check =====\n")
cat("Hypoxia core genes required:\n")
print(hypoxia_genes_core)

cat("Missing hypoxia genes:\n")
print(missing_hypoxia_genes)

if (length(missing_hypoxia_genes) > 0) {
  stop(
    "Missing hypoxia core genes. Cannot calculate hypoxia score.\n",
    "Inspect the missing hypoxia genes before continuing."
  )
}

hypoxia_expr <- expr_proc[hypoxia_genes_core, , drop = FALSE]

# Gene-wise z-score across patients, then average across genes.
hypoxia_z_by_gene <- t(scale(t(hypoxia_expr)))

Hypoxia_core10_score <- colMeans(hypoxia_z_by_gene, na.rm = TRUE)
Hypoxia_core10_z <- as.numeric(scale(Hypoxia_core10_score))

names(Hypoxia_core10_score) <- colnames(expr_proc)
names(Hypoxia_core10_z) <- colnames(expr_proc)

cat("\n===== Hypoxia core-10 score summary =====\n")
print(summary(Hypoxia_core10_score))

cat("\n===== Hypoxia core-10 z-score summary =====\n")
print(summary(Hypoxia_core10_z))

############################################################
# 9. Build molecular score table
############################################################

score_df <- data.frame(
  sampleid = colnames(expr_proc),
  CRS = CRS[colnames(expr_proc)],
  CRS_z = CRS_z[colnames(expr_proc)],
  RSI = RSI[colnames(expr_proc)],
  RSI_z = RSI_z[colnames(expr_proc)],
  Hypoxia_core10_score = Hypoxia_core10_score[colnames(expr_proc)],
  Hypoxia_core10_z = Hypoxia_core10_z[colnames(expr_proc)],
  CRS_group_median = ifelse(
    CRS[colnames(expr_proc)] >= median(CRS, na.rm = TRUE),
    "CRS_high_radioresistant_like",
    "CRS_low_radiosensitive_like"
  ),
  RSI_group_median = ifelse(
    RSI[colnames(expr_proc)] >= median(RSI, na.rm = TRUE),
    "RSI_high_radioresistant_like",
    "RSI_low_radiosensitive_like"
  ),
  Hypoxia_group_median = ifelse(
    Hypoxia_core10_score[colnames(expr_proc)] >= median(Hypoxia_core10_score, na.rm = TRUE),
    "Hypoxia_high",
    "Hypoxia_low"
  ),
  stringsAsFactors = FALSE
)

cat("\n===== Molecular score table preview =====\n")
print(head(score_df, 20))

cat("\n===== Molecular score table dim =====\n")
print(dim(score_df))

############################################################
# 10. Correlation audit
############################################################

cor_pairs <- data.frame(
  pair = c(
    "CRS_vs_RSI",
    "CRS_vs_Hypoxia",
    "RSI_vs_Hypoxia"
  ),
  pearson = c(
    cor(score_df$CRS, score_df$RSI, method = "pearson", use = "complete.obs"),
    cor(score_df$CRS, score_df$Hypoxia_core10_score, method = "pearson", use = "complete.obs"),
    cor(score_df$RSI, score_df$Hypoxia_core10_score, method = "pearson", use = "complete.obs")
  ),
  spearman = c(
    cor(score_df$CRS, score_df$RSI, method = "spearman", use = "complete.obs"),
    cor(score_df$CRS, score_df$Hypoxia_core10_score, method = "spearman", use = "complete.obs"),
    cor(score_df$RSI, score_df$Hypoxia_core10_score, method = "spearman", use = "complete.obs")
  ),
  stringsAsFactors = FALSE
)

cat("\n===== Molecular score correlation audit =====\n")
print(cor_pairs)

############################################################
# 11. Save tables
############################################################

write.csv(
  score_df,
  file = "05_molecular_scores/NSCLC_Radiogenomics_patient_molecular_scores.csv",
  row.names = FALSE
)

saveRDS(
  score_df,
  file = "05_molecular_scores/NSCLC_Radiogenomics_patient_molecular_scores.rds"
)

wb <- createWorkbook()

addWorksheet(wb, "patient_scores")
writeData(wb, "patient_scores", score_df)

addWorksheet(wb, "raw_expression_summary")
writeData(wb, "raw_expression_summary", raw_summary)

addWorksheet(wb, "processed_expression_summary")
writeData(wb, "processed_expression_summary", proc_summary)

addWorksheet(wb, "CRS_coefficients_used")
writeData(wb, "CRS_coefficients_used", crs_coef)

addWorksheet(wb, "RSI_coefficients")
writeData(
  wb,
  "RSI_coefficients",
  data.frame(gene = names(rsi_coef), coefficient = as.numeric(rsi_coef))
)

addWorksheet(wb, "hypoxia_core10_genes")
writeData(wb, "hypoxia_core10_genes", data.frame(gene = hypoxia_genes_core))

addWorksheet(wb, "correlation_audit")
writeData(wb, "correlation_audit", cor_pairs)

addWorksheet(wb, "notes")
writeData(
  wb,
  "notes",
  data.frame(
    item = c(
      "CRS_direction",
      "CRS_note",
      "RSI_direction",
      "Hypoxia_direction",
      "preprocessing"
    ),
    value = c(
      "Higher CRS indicates higher predicted AUC and greater radioresistance-associated molecular phenotype.",
      "CRS was approximately transferred from Ensembl to Gene Symbol using HGNC mapping; 92 of 112 non-zero CRS features were available in GSE103584.",
      "Higher RSI indicates more radioresistant-like phenotype; lower RSI indicates more radiosensitive-like phenotype.",
      "Higher hypoxia core-10 score indicates more hypoxia-like transcriptional phenotype.",
      preprocess_note
    ),
    stringsAsFactors = FALSE
  )
)

saveWorkbook(
  wb,
  file = "05_molecular_scores/NSCLC_Radiogenomics_patient_molecular_scores.xlsx",
  overwrite = TRUE
)

############################################################
# 12. Save simple figures
############################################################

png("08_figures_final/patient_molecular_scores/patient_crs_distribution.png", width = 1600, height = 1200, res = 200)
hist(score_df$CRS, breaks = 25, main = "Patient CRS distribution", xlab = "CRS")
abline(v = median(score_df$CRS, na.rm = TRUE), lwd = 2, lty = 2)
dev.off()

png("08_figures_final/patient_molecular_scores/patient_rsi_distribution.png", width = 1600, height = 1200, res = 200)
hist(score_df$RSI, breaks = 25, main = "Patient RSI distribution", xlab = "RSI")
abline(v = median(score_df$RSI, na.rm = TRUE), lwd = 2, lty = 2)
dev.off()

png("08_figures_final/patient_molecular_scores/patient_hypoxia_module_distribution.png", width = 1600, height = 1200, res = 200)
hist(score_df$Hypoxia_core10_score, breaks = 25, main = "Hypoxia core-10 score distribution", xlab = "Hypoxia core-10 score")
abline(v = median(score_df$Hypoxia_core10_score, na.rm = TRUE), lwd = 2, lty = 2)
dev.off()

png("08_figures_final/patient_molecular_scores/patient_crs_vs_rsi.png", width = 1600, height = 1200, res = 200)
plot(score_df$CRS, score_df$RSI, pch = 16, xlab = "CRS", ylab = "RSI", main = "CRS vs RSI")
abline(lm(RSI ~ CRS, data = score_df), lwd = 2)
dev.off()

png("08_figures_final/patient_molecular_scores/patient_crs_vs_hypoxia_module.png", width = 1600, height = 1200, res = 200)
plot(score_df$CRS, score_df$Hypoxia_core10_score, pch = 16, xlab = "CRS", ylab = "Hypoxia core-10 score", main = "CRS vs Hypoxia core-10")
abline(lm(Hypoxia_core10_score ~ CRS, data = score_df), lwd = 2)
dev.off()

############################################################
# 13. Finish
############################################################

cat("\n===== Patient molecular-score calculation completed =====\n")
cat("Main output files:\n")
cat("1. 05_molecular_scores/NSCLC_Radiogenomics_patient_molecular_scores.xlsx\n")
cat("2. 05_molecular_scores/NSCLC_Radiogenomics_patient_molecular_scores.csv\n")
cat("3. 05_molecular_scores/NSCLC_Radiogenomics_patient_molecular_scores.rds\n")
cat("4. 01_raw_data/NSCLC_Radiogenomics/GSE103584_NSCLC_RNAseq_SYMBOL_gene_by_sample_LOG2_PROCESSED.rds\n")
cat("5. 08_figures_final/patient_molecular_scores/patient_crs_distribution.png\n")
cat("6. 08_figures_final/patient_molecular_scores/patient_rsi_distribution.png\n")
cat("7. 08_figures_final/patient_molecular_scores/patient_hypoxia_module_distribution.png\n")
cat("8. 08_figures_final/patient_molecular_scores/patient_crs_vs_rsi.png\n")
cat("9. 08_figures_final/patient_molecular_scores/patient_crs_vs_hypoxia_module.png\n")

cat("1. Processed expression value summary\n")
cat("2. CRS score summary\n")
cat("3. RSI score summary\n")
cat("4. Hypoxia core-10 score summary\n")
cat("5. Molecular score correlation audit\n")
cat("6. Molecular score table preview\n")
