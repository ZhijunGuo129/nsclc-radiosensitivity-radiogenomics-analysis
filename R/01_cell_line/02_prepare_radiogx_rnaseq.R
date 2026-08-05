# 02_prepare_radiogx_rnaseq.R
#
# Match RadioGx RNA-seq profiles to cell lines and radiation-response outcomes.
#
# Run from an analysis workspace configured through environment variables.

options(stringsAsFactors = FALSE)
options(timeout = 1800)

required_packages <- c("RadioGx", "CoreGx", "openxlsx")
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

library(RadioGx)
library(CoreGx)
library(openxlsx)

dir.create("02_metadata", recursive = TRUE, showWarnings = FALSE)
dir.create("01_raw_data/RadioGx", recursive = TRUE, showWarnings = FALSE)
dir.create("05_molecular_scores", recursive = TRUE, showWarnings = FALSE)

############################################################
# 1. Load Cleveland
############################################################

rds_path <- "01_raw_data/RadioGx/Cleveland_or_clevelandSmall_RadioSet.rds"

if (file.exists(rds_path)) {
  Cleveland <- readRDS(rds_path)
} else {
  Cleveland <- downloadRSet(
    name = "Cleveland",
    saveDir = "01_raw_data/RadioGx"
  )
  saveRDS(Cleveland, rds_path)
}

cell_info <- as.data.frame(cellInfo(Cleveland))
sens_info <- as.data.frame(sensitivityInfo(Cleveland))
sens_prof <- as.data.frame(sensitivityProfiles(Cleveland))
sens_prof$profile_rowname <- rownames(sens_prof)

cat("Cleveland loaded.\n")
print(Cleveland)

############################################################
# 2. Helper functions
############################################################

norm_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]", "", x)
  return(x)
}

parse_rnaseq_cellline <- function(x) {
  x <- as.character(x)
  
  # Remove prefix like G20468.
  x <- sub("^G[0-9]+\\.", "", x)
  
  # Remove suffix like .2 or .1
  x <- sub("\\.[0-9]+$", "", x)
  
  return(x)
}

############################################################
# 3. Extract RNAseq matrix
############################################################

rnaseq_mat <- as.matrix(molecularProfiles(Cleveland, "rnaseq"))

cat("\nRNAseq matrix dim:\n")
print(dim(rnaseq_mat))

rnaseq_colnames <- colnames(rnaseq_mat)

rnaseq_parse_table <- data.frame(
  expr_colname = rnaseq_colnames,
  parsed_cellline = parse_rnaseq_cellline(rnaseq_colnames),
  parsed_norm = norm_id(parse_rnaseq_cellline(rnaseq_colnames)),
  expr_col_index = seq_along(rnaseq_colnames),
  stringsAsFactors = FALSE
)

cat("\nRNAseq parsed column preview:\n")
print(head(rnaseq_parse_table, 20))

############################################################
# 4. Match parsed RNAseq names to cell_info$CellLine
############################################################

cell_info$CellLine_norm <- norm_id(cell_info$CellLine)
cell_info$sampleid <- as.character(cell_info$sampleid)

rnaseq_mapping <- merge(
  rnaseq_parse_table,
  cell_info,
  by.x = "parsed_norm",
  by.y = "CellLine_norm",
  all.x = TRUE
)

rnaseq_mapping <- rnaseq_mapping[order(rnaseq_mapping$expr_col_index), ]

cat("\nRNAseq matching summary:\n")
cat("RNAseq columns:", nrow(rnaseq_parse_table), "\n")
cat("Matched to cell_info:", sum(!is.na(rnaseq_mapping$sampleid)), "\n")
cat("Unmatched:", sum(is.na(rnaseq_mapping$sampleid)), "\n")

cat("\nUnmatched RNAseq examples:\n")
print(head(rnaseq_mapping[is.na(rnaseq_mapping$sampleid), c("expr_colname", "parsed_cellline")], 30))

############################################################
# 5. Build matched RNAseq matrix
############################################################

matched_mapping <- rnaseq_mapping[!is.na(rnaseq_mapping$sampleid), ]

# Remove duplicated sampleid if any, keep first occurrence
matched_mapping <- matched_mapping[!duplicated(matched_mapping$sampleid), ]

matched_indices <- matched_mapping$expr_col_index
rnaseq_matched <- rnaseq_mat[, matched_indices, drop = FALSE]

colnames(rnaseq_matched) <- matched_mapping$sampleid

cat("\nMatched RNAseq matrix dim:\n")
print(dim(rnaseq_matched))

############################################################
# 6. Map sensitivity outcomes to sampleid
############################################################

sens_prof$outcome_cellline <- sub("_radiation_.*$", "", sens_prof$profile_rowname)
sens_prof$outcome_cellline_norm <- norm_id(sens_prof$outcome_cellline)

cell_info_for_sens <- cell_info
cell_info_for_sens$CellLine_norm <- norm_id(cell_info_for_sens$CellLine)

sens_outcome <- merge(
  sens_prof,
  cell_info_for_sens,
  by.x = "outcome_cellline_norm",
  by.y = "CellLine_norm",
  all.x = TRUE
)

cat("\nSensitivity outcome mapping:\n")
cat("Total sensitivity profiles:", nrow(sens_prof), "\n")
cat("Mapped to cell_info:", sum(!is.na(sens_outcome$sampleid)), "\n")
cat("Unmapped:", sum(is.na(sens_outcome$sampleid)), "\n")

############################################################
# 7. Create modeling dataset
############################################################

# Some cell lines may have multiple sensitivity experiments.
# For safety, average outcomes by sampleid.
sens_model <- aggregate(
  cbind(AUC_published, AUC_recomputed, SF2) ~ sampleid,
  data = sens_outcome,
  FUN = function(x) mean(as.numeric(x), na.rm = TRUE)
)

# Add alpha/beta separately because they may be stored oddly
alpha_numeric <- suppressWarnings(as.numeric(sens_outcome$alpha))
beta_numeric <- suppressWarnings(as.numeric(sens_outcome$beta))

sens_outcome$alpha_numeric <- alpha_numeric
sens_outcome$beta_numeric <- beta_numeric

alpha_beta_model <- aggregate(
  cbind(alpha_numeric, beta_numeric) ~ sampleid,
  data = sens_outcome,
  FUN = function(x) mean(as.numeric(x), na.rm = TRUE)
)

sens_model <- merge(
  sens_model,
  alpha_beta_model,
  by = "sampleid",
  all.x = TRUE
)

model_ids <- intersect(colnames(rnaseq_matched), sens_model$sampleid)

rnaseq_model_mat <- rnaseq_matched[, model_ids, drop = FALSE]
sens_model_final <- sens_model[match(model_ids, sens_model$sampleid), ]

# Add cell line annotation
cell_anno <- cell_info[, c("sampleid", "CellLine", "Primarysite", "Histology", "Subhistology")]
sens_model_final <- merge(
  sens_model_final,
  cell_anno,
  by = "sampleid",
  all.x = TRUE,
  sort = FALSE
)

# Restore order
sens_model_final <- sens_model_final[match(model_ids, sens_model_final$sampleid), ]

sens_model_final$Primarysite_lower <- tolower(as.character(sens_model_final$Primarysite))
sens_model_final$is_lung_strict <- grepl(
  "lung|bronch|pulmonary",
  sens_model_final$Primarysite_lower,
  ignore.case = TRUE
)

cat("\n===== FINAL MODELING SUMMARY =====\n")
cat("RNAseq matched expression samples:", ncol(rnaseq_matched), "\n")
cat("RNAseq + sensitivity modeling samples:", length(model_ids), "\n")
cat("Lung modeling samples:", sum(sens_model_final$is_lung_strict, na.rm = TRUE), "\n")

cat("\nOutcome columns:\n")
print(colnames(sens_model_final))

cat("\nOutcome summary:\n")
print(summary(sens_model_final[, c("AUC_published", "AUC_recomputed", "SF2")]))

cat("\nPrimarysite distribution among modeling samples:\n")
print(sort(table(sens_model_final$Primarysite), decreasing = TRUE))

############################################################
# 8. Export
############################################################

wb <- createWorkbook()

addWorksheet(wb, "rnaseq_parse_table")
writeData(wb, "rnaseq_parse_table", rnaseq_parse_table)

addWorksheet(wb, "rnaseq_mapping")
writeData(wb, "rnaseq_mapping", rnaseq_mapping)

addWorksheet(wb, "matched_mapping")
writeData(wb, "matched_mapping", matched_mapping)

addWorksheet(wb, "sens_outcome")
writeData(wb, "sens_outcome", sens_outcome)

addWorksheet(wb, "sens_model_final")
writeData(wb, "sens_model_final", sens_model_final)

saveWorkbook(
  wb,
  file = "02_metadata/RadioGx_RNAseq_Match_Fix.xlsx",
  overwrite = TRUE
)

saveRDS(
  rnaseq_matched,
  file = "01_raw_data/RadioGx/RadioGx_rnaseq_matched_all_samples.rds"
)

saveRDS(
  rnaseq_model_mat,
  file = "01_raw_data/RadioGx/RadioGx_rnaseq_model_gene_by_sample.rds"
)

saveRDS(
  sens_model_final,
  file = "05_molecular_scores/RadioGx_sensitivity_model_outcomes.rds"
)

write.xlsx(
  sens_model_final,
  file = "05_molecular_scores/RadioGx_sensitivity_model_outcomes.xlsx",
  overwrite = TRUE
)

cat("\n===== RadioGx RNA-seq matching completed =====\n")
cat("Main outputs:\n")
cat("1. 02_metadata/RadioGx_RNAseq_Match_Fix.xlsx\n")
cat("2. 01_raw_data/RadioGx/RadioGx_rnaseq_model_gene_by_sample.rds\n")
cat("3. 05_molecular_scores/RadioGx_sensitivity_model_outcomes.xlsx\n")
