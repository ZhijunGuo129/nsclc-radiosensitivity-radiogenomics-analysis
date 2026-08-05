# 01_prepare_transcriptome.R
#
# Audit and preprocess the NSCLC-Radiogenomics transcriptomic matrix.
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

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Missing R package 'openxlsx'. Install dependencies with environment/install_r_dependencies.R.")
}

library(openxlsx)

dirs <- c(
  "01_raw_data/NSCLC_Radiogenomics",
  "02_metadata",
  "05_molecular_scores",
  "07_results",
  "10_logs"
)

for (d in dirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("\n===== Current working directory =====\n")
print(getwd())

############################################################
# 1. Check local file
############################################################

gse_dest <- "01_raw_data/NSCLC_Radiogenomics/GSE103584_R01_NSCLC_RNAseq.txt.gz"

cat("\n===== Checking local RNA-seq file =====\n")
cat("Expected file path:\n")
cat(gse_dest, "\n")

if (!file.exists(gse_dest)) {
  stop(
    "GSE103584_R01_NSCLC_RNAseq.txt.gz was not found locally.\n",
    "Place the file in 01_raw_data/NSCLC_Radiogenomics/."
  )
}

cat("\nRNA-seq file found locally.\n")
print(file.info(gse_dest))

############################################################
# 2. Load CRS model
############################################################

crs_model_path <- "06_models/crs_elastic_net_model.rds"

if (!file.exists(crs_model_path)) {
  stop(
    "Cannot find CRS model file:\n",
    crs_model_path,
    "\nRun the CRS model-training step before this script."
  )
}

crs_model_obj <- readRDS(crs_model_path)
final_features <- crs_model_obj$final_features

cat("\n===== CRS model loaded =====\n")
cat("CRS model feature number:\n")
print(length(final_features))

cat("\nCRS model endpoint:\n")
print(crs_model_obj$endpoint)

cat("\nCRS direction:\n")
print(crs_model_obj$direction)

cat("\nCRS gene ID type:\n")
if (!is.null(crs_model_obj$gene_id_type)) {
  print(crs_model_obj$gene_id_type)
} else {
  print("Not recorded")
}

############################################################
# 3. Read GSE103584 table
############################################################

cat("\n===== Reading GSE103584 RNA-seq table =====\n")

expr_df <- read.delim(
  gzfile(gse_dest),
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

cat("\n===== Expression data frame dim =====\n")
print(dim(expr_df))

cat("\n===== Raw column names preview =====\n")
print(head(colnames(expr_df), 50))

cat("\n===== Raw data preview: first 8 rows and 8 columns =====\n")
print(expr_df[1:min(8, nrow(expr_df)), 1:min(8, ncol(expr_df))])

############################################################
# 4. Identify gene column and sample columns
############################################################

# Important:
# The first column has an empty column name "".
# But it contains gene symbols such as A1BG, A2M.
gene_col_index <- 1
gene_col_name <- colnames(expr_df)[gene_col_index]

raw_gene_ids <- as.character(expr_df[[gene_col_index]])

clean_gene_symbol <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- toupper(x)
  return(x)
}

gene_ids_clean <- clean_gene_symbol(raw_gene_ids)

cat("\n===== Detected gene column =====\n")
cat("Gene column index:\n")
print(gene_col_index)

cat("Gene column name:\n")
print(gene_col_name)

cat("\n===== Gene ID preview =====\n")
print(head(raw_gene_ids, 30))

cat("\n===== Cleaned gene ID preview =====\n")
print(head(gene_ids_clean, 30))

ensembl_fraction <- mean(grepl("^ENSG", gene_ids_clean), na.rm = TRUE)
symbol_like_fraction <- mean(grepl("^[A-Z0-9\\-\\.]+$", gene_ids_clean), na.rm = TRUE)

gene_id_type <- ifelse(
  ensembl_fraction > 0.5,
  "Ensembl",
  "GeneSymbol"
)

cat("\n===== Gene ID type audit =====\n")
cat("Ensembl ID fraction:\n")
print(ensembl_fraction)

cat("Symbol-like fraction:\n")
print(symbol_like_fraction)

cat("Inferred patient gene ID type:\n")
print(gene_id_type)

sample_cols <- grep("^R01-", colnames(expr_df), value = TRUE)

cat("\n===== Number of detected patient sample columns =====\n")
print(length(sample_cols))

cat("\n===== Sample column preview =====\n")
print(head(sample_cols, 50))

if (length(sample_cols) == 0) {
  stop("No patient sample columns beginning with R01- were detected.")
}

############################################################
# 5. Build expression matrix
############################################################

expr_mat <- as.matrix(expr_df[, sample_cols, drop = FALSE])
mode(expr_mat) <- "numeric"

rownames(expr_mat) <- gene_ids_clean

keep_gene <- !is.na(rownames(expr_mat)) & rownames(expr_mat) != ""
expr_mat <- expr_mat[keep_gene, , drop = FALSE]

cat("\n===== Expression matrix before duplicated-gene collapsing =====\n")
print(dim(expr_mat))

cat("\nNumber of duplicated gene IDs:\n")
print(sum(duplicated(rownames(expr_mat))))

# Collapse duplicated genes by mean while handling NA
collapse_duplicate_genes <- function(mat) {
  gene_group <- rownames(mat)
  
  mat_zero <- mat
  mat_zero[is.na(mat_zero)] <- 0
  
  sums <- rowsum(mat_zero, group = gene_group, reorder = FALSE)
  counts <- rowsum(!is.na(mat) * 1, group = gene_group, reorder = FALSE)
  
  out <- sums / counts
  out[counts == 0] <- NA
  
  return(out)
}

if (any(duplicated(rownames(expr_mat)))) {
  cat("\nDuplicated gene IDs detected. Collapsing duplicated genes by mean...\n")
  expr_mat <- collapse_duplicate_genes(expr_mat)
}

cat("\n===== Final patient expression matrix dim, genes x samples =====\n")
print(dim(expr_mat))

cat("\n===== Final gene ID preview =====\n")
print(head(rownames(expr_mat), 30))

cat("\n===== Final patient/sample ID preview =====\n")
print(head(colnames(expr_mat), 50))

############################################################
# 6. Missingness audit
############################################################

sample_missing_rate <- colMeans(is.na(expr_mat))
gene_missing_rate <- rowMeans(is.na(expr_mat))

missingness_summary <- data.frame(
  item = c(
    "sample_missing_rate_min",
    "sample_missing_rate_median",
    "sample_missing_rate_max",
    "gene_missing_rate_min",
    "gene_missing_rate_median",
    "gene_missing_rate_max"
  ),
  value = c(
    min(sample_missing_rate, na.rm = TRUE),
    median(sample_missing_rate, na.rm = TRUE),
    max(sample_missing_rate, na.rm = TRUE),
    min(gene_missing_rate, na.rm = TRUE),
    median(gene_missing_rate, na.rm = TRUE),
    max(gene_missing_rate, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

cat("\n===== Missingness summary =====\n")
print(missingness_summary)

############################################################
# 7. Expression value audit
############################################################

finite_values <- as.numeric(expr_mat[is.finite(expr_mat)])

expr_value_summary <- data.frame(
  min = min(finite_values, na.rm = TRUE),
  q25 = as.numeric(quantile(finite_values, 0.25, na.rm = TRUE)),
  median = median(finite_values, na.rm = TRUE),
  mean = mean(finite_values, na.rm = TRUE),
  q75 = as.numeric(quantile(finite_values, 0.75, na.rm = TRUE)),
  q95 = as.numeric(quantile(finite_values, 0.95, na.rm = TRUE)),
  q99 = as.numeric(quantile(finite_values, 0.99, na.rm = TRUE)),
  max = max(finite_values, na.rm = TRUE),
  stringsAsFactors = FALSE
)

expr_value_summary$patient_log_recommendation <- ifelse(
  expr_value_summary$q99 > 50,
  "Likely raw/unlogged scale; consider log2(x + 1)",
  "Likely already normalized/log-scale"
)

cat("\n===== Expression value summary =====\n")
print(expr_value_summary)

############################################################
# 8. CRS direct feature overlap
############################################################

overlap_features_direct <- intersect(final_features, rownames(expr_mat))
missing_features_direct <- setdiff(final_features, rownames(expr_mat))

feature_overlap_direct <- data.frame(
  item = c(
    "CRS_model_features",
    "Patient_expression_genes",
    "Direct_overlapping_CRS_features",
    "Direct_missing_CRS_features",
    "Direct_overlap_fraction",
    "CRS_model_gene_id_type",
    "Patient_expression_gene_id_type"
  ),
  value = c(
    length(final_features),
    nrow(expr_mat),
    length(overlap_features_direct),
    length(missing_features_direct),
    round(length(overlap_features_direct) / length(final_features), 4),
    ifelse(!is.null(crs_model_obj$gene_id_type), crs_model_obj$gene_id_type, "Unknown"),
    gene_id_type
  ),
  stringsAsFactors = FALSE
)

cat("\n===== CRS model feature direct overlap =====\n")
print(feature_overlap_direct)

cat("\n===== Direct missing CRS features preview =====\n")
print(head(missing_features_direct, 60))

############################################################
# 9. RSI and hypoxia gene availability
############################################################

# Note:
# RSI exact published implementation may require careful formula checking later.
# Here we only audit symbol availability.
rsi_genes_common <- c(
  "AR", "JUN", "STAT1", "PRKCB", "RELA",
  "ABL1", "SUMO1", "PAK2", "HDAC1", "IRF1"
)

hypoxia_genes_small <- c(
  "HIF1A", "VEGFA", "CA9", "SLC2A1", "LDHA",
  "PDK1", "BNIP3", "NDRG1", "EGLN1", "EGLN3"
)

rsi_hypoxia_availability <- data.frame(
  gene_set = c("RSI_common_symbol_check", "Hypoxia_small_symbol_check"),
  n_defined = c(length(rsi_genes_common), length(hypoxia_genes_small)),
  n_present = c(
    length(intersect(rsi_genes_common, rownames(expr_mat))),
    length(intersect(hypoxia_genes_small, rownames(expr_mat)))
  ),
  genes_present = c(
    paste(intersect(rsi_genes_common, rownames(expr_mat)), collapse = ";"),
    paste(intersect(hypoxia_genes_small, rownames(expr_mat)), collapse = ";")
  ),
  stringsAsFactors = FALSE
)

cat("\n===== RSI / hypoxia gene symbol availability =====\n")
print(rsi_hypoxia_availability)

############################################################
# 10. Molecular flow audit
############################################################

molecular_flow <- data.frame(
  step = c(
    "GSE103584 detected RNA-seq patient columns",
    "Patient RNA-seq samples after expression matrix construction",
    "Patient RNA-seq genes after duplicate collapsing",
    "Direct CRS feature overlap before ID mapping",
    "Patient samples eligible for RSI/hypoxia scoring",
    "Patient samples eligible for CRS scoring after Ensembl-Symbol mapping"
  ),
  n = c(
    length(sample_cols),
    ncol(expr_mat),
    nrow(expr_mat),
    length(overlap_features_direct),
    ncol(expr_mat),
    NA
  ),
  note = c(
    "R01-* columns detected from GEO processed matrix",
    "Expression matrix successfully constructed",
    "Gene symbols used as row names",
    "Expected to be low because CRS uses Ensembl and patient matrix uses Gene Symbol",
    "RSI and hypoxia use gene symbols, so these can be computed next",
    "Requires the CRS Ensembl-to-gene-symbol mapping"
  ),
  stringsAsFactors = FALSE
)

cat("\n===== Molecular sample-flow audit =====\n")
print(molecular_flow)

############################################################
# 11. Save outputs
############################################################

wb <- createWorkbook()

addWorksheet(wb, "raw_preview")
writeData(
  wb,
  "raw_preview",
  expr_df[1:min(100, nrow(expr_df)), 1:min(50, ncol(expr_df))]
)

addWorksheet(wb, "gene_id_type_audit")
writeData(
  wb,
  "gene_id_type_audit",
  data.frame(
    item = c("ensembl_fraction", "symbol_like_fraction", "gene_id_type"),
    value = c(ensembl_fraction, symbol_like_fraction, gene_id_type),
    stringsAsFactors = FALSE
  )
)

addWorksheet(wb, "missingness_summary")
writeData(wb, "missingness_summary", missingness_summary)

addWorksheet(wb, "expression_value_summary")
writeData(wb, "expression_value_summary", expr_value_summary)

addWorksheet(wb, "feature_overlap_direct")
writeData(wb, "feature_overlap_direct", feature_overlap_direct)

addWorksheet(wb, "missing_CRS_features_direct")
writeData(wb, "missing_CRS_features_direct", data.frame(feature = missing_features_direct))

addWorksheet(wb, "sample_ids")
writeData(wb, "sample_ids", data.frame(sampleid = colnames(expr_mat)))

addWorksheet(wb, "rsi_hypoxia_availability")
writeData(wb, "rsi_hypoxia_availability", rsi_hypoxia_availability)

addWorksheet(wb, "molecular_flow")
writeData(wb, "molecular_flow", molecular_flow)

saveWorkbook(
  wb,
  file = "02_metadata/nsclc_radiogenomics_transcriptome_summary.xlsx",
  overwrite = TRUE
)

saveRDS(
  expr_mat,
  file = "01_raw_data/NSCLC_Radiogenomics/GSE103584_NSCLC_RNAseq_SYMBOL_gene_by_sample_RAW_AUDITED.rds"
)

write.csv(
  data.frame(sampleid = colnames(expr_mat)),
  file = "02_metadata/GSE103584_detected_RNAseq_sample_ids.csv",
  row.names = FALSE
)

cat("\n===== NSCLC-Radiogenomics RNA-seq audit completed =====\n")
cat("Main output files:\n")
cat("1. 02_metadata/nsclc_radiogenomics_transcriptome_summary.xlsx\n")
cat("2. 01_raw_data/NSCLC_Radiogenomics/GSE103584_NSCLC_RNAseq_SYMBOL_gene_by_sample_RAW_AUDITED.rds\n")
cat("3. 02_metadata/GSE103584_detected_RNAseq_sample_ids.csv\n")

cat("1. Expression data frame dim\n")
cat("2. Gene ID type audit\n")
cat("3. Number of detected patient sample columns\n")
cat("4. Final patient expression matrix dim\n")
cat("5. Missingness summary\n")
cat("6. Expression value summary\n")
cat("7. CRS model feature direct overlap\n")
cat("8. RSI / hypoxia gene symbol availability\n")
cat("9. Molecular sample-flow audit\n")
