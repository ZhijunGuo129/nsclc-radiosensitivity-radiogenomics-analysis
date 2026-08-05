# 04_compute_lung3_crs.R
#
# Map Lung3 probes to gene symbols and compute CRS.
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

required_packages <- c("openxlsx")
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

library(openxlsx)

############################################################
# 1. Paths
############################################################


lung3_dir <- file.path(
  project_dir,
  "01_raw_data",
  "Lung3"
)

out_dir <- file.path(
  project_dir,
  "07_results",
  "lung3_validation"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expr_rds <- file.path(
  project_dir,
  "07_results",
  "lung3_data_audit",
  "lung3_data_audit_GSE58661_expression_table_probe_by_sample.rds"
)

sample_meta_csv <- file.path(
  project_dir,
  "07_results",
  "lung3_data_audit",
  "lung3_data_audit_GSE58661_sample_metadata.csv"
)

clinical_csv <- file.path(
  project_dir,
  "07_results",
  "lung3_data_audit",
  "lung3_data_audit_Lung3_clinical_clean_initial.csv"
)

gpl_soft <- file.path(
  lung3_dir,
  "GPL15048_family.soft.gz"
)

crs_model_rds <- file.path(
  project_dir,
  "06_models",
  "crs_gene_symbol_coefficients.rds"
)

if (!file.exists(expr_rds)) {
  stop("Cannot find Lung3 expression table RDS from the Lung3 data audit: ", expr_rds)
}

if (!file.exists(gpl_soft)) {
  stop("Cannot find GPL15048_family.soft.gz. Please download it manually to: ", gpl_soft)
}

if (!file.exists(crs_model_rds)) {
  stop("Cannot find CRS symbol model: ", crs_model_rds)
}

############################################################
# 2. Helpers
############################################################

clean_symbol_one <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x == "" | is.na(x)] <- NA
  
  # GEO/Affy annotations may contain multiple symbols separated by ///, //, ;, or ,
  x <- sapply(x, function(z) {
    if (is.na(z)) return(NA_character_)
    parts <- unlist(strsplit(z, "///|//|;|,", perl = TRUE))
    parts <- toupper(trimws(parts))
    parts <- parts[parts != "" & !is.na(parts)]
    if (length(parts) == 0) return(NA_character_)
    parts[1]
  })
  
  x <- gsub("\\s+", "", x)
  x
}

clean_id <- function(x) {
  toupper(trimws(as.character(x)))
}

extract_lung_number <- function(x) {
  x <- clean_id(x)
  x <- gsub("LUNG3", "LUNG", x)
  x <- gsub("[^A-Z0-9]+", "_", x)
  
  num <- rep(NA_character_, length(x))
  
  hit1 <- grepl("LUNG[_]?0*[0-9]+", x)
  num[hit1] <- sub(".*LUNG[_]?0*([0-9]+).*", "\\1", x[hit1])
  
  hit2 <- is.na(num) & grepl("[0-9]+", x)
  num[hit2] <- sub(".*?([0-9]+).*", "\\1", x[hit2])
  
  num
}

############################################################
# 3. Load Lung3 expression table
############################################################

cat("\n===== Loading Lung3 probe-level expression =====\n")

expr_table <- readRDS(expr_rds)
expr_table <- as.data.frame(expr_table, check.names = FALSE)

cat("Probe-level expression table dim:\n")
print(dim(expr_table))

id_col <- colnames(expr_table)[1]

probe_ids <- as.character(expr_table[[id_col]])

expr_mat <- as.matrix(expr_table[, -1, drop = FALSE])
mode(expr_mat) <- "numeric"

rownames(expr_mat) <- probe_ids

cat("\nExpression numeric summary:\n")
print(summary(as.vector(expr_mat)))

q99 <- quantile(as.vector(expr_mat), 0.99, na.rm = TRUE)

log_note <- ifelse(q99 > 100, "Values look unlogged; log2(x+1) applied.", "Values look already log-scale; no log transform applied.")

if (q99 > 100) {
  expr_mat <- log2(expr_mat + 1)
}

cat("\nPreprocess note:\n")
cat(log_note, "\n")

cat("\nExpression summary after preprocessing:\n")
print(summary(as.vector(expr_mat)))

############################################################
############################################################
# 4. Parse GPL15048 annotation - FAST CHUNK VERSION
############################################################

cat("\n===== Parsing GPL15048 annotation: FAST CHUNK VERSION =====\n")

cat("GPL file path:\n")
cat(gpl_soft, "\n")

cat("GPL file size bytes:\n")
print(file.info(gpl_soft)$size)

extract_platform_table_from_soft <- function(soft_path, chunk_n = 5000) {
  
  con <- gzfile(soft_path, open = "rt")
  on.exit(close(con), add = TRUE)
  
  in_table <- FALSE
  platform_lines <- character()
  total_lines <- 0
  
  repeat {
    
    chunk <- readLines(con, n = chunk_n, warn = FALSE)
    
    if (length(chunk) == 0) break
    
    total_lines <- total_lines + length(chunk)
    
    if (total_lines %% 50000 < chunk_n) {
      cat("Read lines:", total_lines, "\n")
      flush.console()
    }
    
    if (!in_table) {
      
      begin_hit <- grep("^!platform_table_begin", chunk)
      
      if (length(begin_hit) > 0) {
        in_table <- TRUE
        
        start_pos <- begin_hit[1] + 1
        
        if (start_pos <= length(chunk)) {
          chunk2 <- chunk[start_pos:length(chunk)]
        } else {
          chunk2 <- character()
        }
        
      } else {
        next
      }
      
    } else {
      chunk2 <- chunk
    }
    
    if (in_table) {
      
      end_hit <- grep("^!platform_table_end", chunk2)
      
      if (length(end_hit) > 0) {
        
        if (end_hit[1] > 1) {
          platform_lines <- c(platform_lines, chunk2[1:(end_hit[1] - 1)])
        }
        
        cat("Found platform_table_end.\n")
        cat("Total lines read:", total_lines, "\n")
        cat("Platform table lines collected:", length(platform_lines), "\n")
        
        break
        
      } else {
        
        platform_lines <- c(platform_lines, chunk2)
      }
    }
  }
  
  if (length(platform_lines) == 0) {
    stop("No platform table lines extracted from SOFT file.")
  }
  
  platform_lines
}

platform_lines <- extract_platform_table_from_soft(gpl_soft)

cat("\nFirst 3 platform table lines:\n")
print(head(platform_lines, 3))

platform_table <- read.delim(
  text = paste(platform_lines, collapse = "\n"),
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  comment.char = "",
  sep = "\t"
)

cat("\nGPL15048 platform table dim:\n")
print(dim(platform_table))

cat("\nGPL15048 columns:\n")
print(colnames(platform_table))

if (!"ID" %in% colnames(platform_table)) {
  stop("GPL table has no ID column.")
}

symbol_col_candidates <- c(
  "GeneSymbol",
  "Gene Symbol",
  "GENE_SYMBOL",
  "Symbol",
  "gene_symbol"
)

symbol_col <- symbol_col_candidates[symbol_col_candidates %in% colnames(platform_table)][1]

if (is.na(symbol_col) || length(symbol_col) == 0) {
  
  cat("\nCannot directly detect GeneSymbol column. Searching columns containing symbol...\n")
  print(colnames(platform_table)[grepl("symbol", colnames(platform_table), ignore.case = TRUE)])
  
  stop("Cannot detect GeneSymbol column in GPL15048 table.")
}

cat("\nDetected symbol column:\n")
print(symbol_col)

annot <- data.frame(
  probe_id = as.character(platform_table$ID),
  gene_symbol_raw = as.character(platform_table[[symbol_col]]),
  stringsAsFactors = FALSE
)

annot$gene_symbol <- clean_symbol_one(annot$gene_symbol_raw)

annot <- annot[
  !is.na(annot$probe_id) &
    annot$probe_id != "" &
    !is.na(annot$gene_symbol) &
    annot$gene_symbol != "",
]

annot <- annot[!duplicated(annot$probe_id), ]

cat("\nClean annotation dim:\n")
print(dim(annot))

cat("\nFirst 10 clean annotation rows:\n")
print(head(annot, 10))
############################################################
# 5. Map probes to gene symbols
############################################################

cat("\n===== Mapping probes to gene symbols =====\n")

common_probes <- intersect(rownames(expr_mat), annot$probe_id)

cat("Expression probes:", nrow(expr_mat), "\n")
cat("Annotation probes with symbols:", nrow(annot), "\n")
cat("Common probes:", length(common_probes), "\n")

expr_use <- expr_mat[common_probes, , drop = FALSE]
annot_use <- annot[match(common_probes, annot$probe_id), , drop = FALSE]

stopifnot(all(rownames(expr_use) == annot_use$probe_id))

# remove controls and invalid gene names
valid_symbol <- grepl("^[A-Z0-9.-]+$", annot_use$gene_symbol) &
  !grepl("^AFFX", annot_use$gene_symbol)

expr_use <- expr_use[valid_symbol, , drop = FALSE]
annot_use <- annot_use[valid_symbol, , drop = FALSE]

cat("Valid symbol-mapped probes after filtering:", nrow(expr_use), "\n")
cat("Unique mapped gene symbols:", length(unique(annot_use$gene_symbol)), "\n")

# Collapse multiple probes per gene:
# choose the probe with highest mean expression for each gene
probe_mean <- rowMeans(expr_use, na.rm = TRUE)

probe_gene_df <- data.frame(
  probe_id = rownames(expr_use),
  gene_symbol = annot_use$gene_symbol,
  probe_mean = probe_mean,
  stringsAsFactors = FALSE
)

probe_gene_df <- probe_gene_df[
  order(probe_gene_df$gene_symbol, -probe_gene_df$probe_mean),
]

best_probe_df <- probe_gene_df[!duplicated(probe_gene_df$gene_symbol), ]

expr_gene <- expr_use[best_probe_df$probe_id, , drop = FALSE]
rownames(expr_gene) <- best_probe_df$gene_symbol

cat("\nGene-level expression dim, genes x samples:\n")
print(dim(expr_gene))

cat("\nFirst 10 gene symbols:\n")
print(head(rownames(expr_gene), 10))

############################################################
# ############################################################
# 6. Load CRS model and calculate Lung3 CRS - revised implementation
############################################################

cat("\n===== Calculating Lung3 CRS =====\n")

crs_model <- readRDS(crs_model_rds)

cat("\nCRS model object names:\n")
print(names(crs_model))

cat("\nStructure of crs_model$symbol_coefficients:\n")
print(str(crs_model$symbol_coefficients, max.level = 2))

extract_coef_vector <- function(crs_model) {
  
  coef_raw <- crs_model$symbol_coefficients
  
  # Case 1: already named numeric vector
  if (is.numeric(coef_raw) && !is.null(names(coef_raw))) {
    coef_vec <- coef_raw
    return(coef_vec)
  }
  
  # Case 2: matrix
  if (is.matrix(coef_raw)) {
    
    df <- as.data.frame(coef_raw, check.names = FALSE)
    df$gene_tmp <- rownames(coef_raw)
    
    numeric_cols <- colnames(df)[sapply(df, is.numeric)]
    numeric_cols <- setdiff(numeric_cols, "gene_tmp")
    
    if (length(numeric_cols) == 0) {
      stop("Cannot find numeric coefficient column in coefficient matrix.")
    }
    
    coef_col <- numeric_cols[which.max(sapply(numeric_cols, function(cc) {
      sum(is.finite(df[[cc]]) & df[[cc]] != 0, na.rm = TRUE)
    }))]
    
    coef_vec <- df[[coef_col]]
    names(coef_vec) <- df$gene_tmp
    
    return(coef_vec)
  }
  
  # Case 3: data.frame
  if (is.data.frame(coef_raw)) {
    
    df <- as.data.frame(coef_raw, check.names = FALSE)
    
    cat("\nCoefficient data.frame columns:\n")
    print(colnames(df))
    
    gene_col_candidates <- c(
      "gene_symbol",
      "symbol",
      "GeneSymbol",
      "gene",
      "feature",
      "features",
      "mapped_symbol",
      "patient_symbol"
    )
    
    gene_col <- gene_col_candidates[gene_col_candidates %in% colnames(df)][1]
    
    if (is.na(gene_col) || length(gene_col) == 0) {
      
      # fallback: use rownames if they look like gene symbols
      if (!is.null(rownames(df)) && !all(rownames(df) %in% as.character(seq_len(nrow(df))))) {
        df$gene_symbol_from_rownames <- rownames(df)
        gene_col <- "gene_symbol_from_rownames"
      } else {
        stop(
          "Cannot detect gene symbol column in coefficient data.frame. Columns are:\n",
          paste(colnames(df), collapse = "; ")
        )
      }
    }
    
    numeric_cols <- colnames(df)[sapply(df, function(x) {
      xx <- suppressWarnings(as.numeric(x))
      sum(is.finite(xx), na.rm = TRUE) > 0
    })]
    
    numeric_cols <- setdiff(numeric_cols, gene_col)
    
    coef_col_candidates <- c(
      "coefficient",
      "coef",
      "glmnet_coef",
      "beta",
      "estimate",
      "Coefficient"
    )
    
    coef_col <- coef_col_candidates[coef_col_candidates %in% numeric_cols][1]
    
    if (is.na(coef_col) || length(coef_col) == 0) {
      coef_col <- numeric_cols[which.max(sapply(numeric_cols, function(cc) {
        xx <- suppressWarnings(as.numeric(df[[cc]]))
        sum(is.finite(xx) & xx != 0, na.rm = TRUE)
      }))]
    }
    
    cat("\nDetected coefficient gene column:\n")
    print(gene_col)
    
    cat("\nDetected coefficient value column:\n")
    print(coef_col)
    
    coef_vec <- suppressWarnings(as.numeric(df[[coef_col]]))
    names(coef_vec) <- toupper(trimws(as.character(df[[gene_col]])))
    
    return(coef_vec)
  }
  
  # Case 4: list
  if (is.list(coef_raw)) {
    
    # Try unlist directly
    coef_try <- suppressWarnings(unlist(coef_raw))
    
    if (is.numeric(coef_try) && length(coef_try) > 0) {
      
      coef_vec <- coef_try
      
      # clean names from list-like names
      nm <- names(coef_vec)
      nm <- gsub(".*\\.", "", nm)
      nm <- gsub("^X", "", nm)
      names(coef_vec) <- nm
      
      return(coef_vec)
    }
  }
  
  stop("Unsupported CRS coefficient object format.")
}

coef_vec_raw <- extract_coef_vector(crs_model)

coef_names <- names(coef_vec_raw)
coef_values <- suppressWarnings(as.numeric(coef_vec_raw))

coef_vec <- coef_values
names(coef_vec) <- toupper(trimws(coef_names))

coef_vec <- coef_vec[
  !is.na(names(coef_vec)) &
    names(coef_vec) != "" &
    is.finite(coef_vec) &
    coef_vec != 0
]

# remove duplicated gene symbols in coefficient vector
coef_vec <- coef_vec[!duplicated(names(coef_vec))]

intercept <- suppressWarnings(as.numeric(crs_model$intercept)[1])

if (!is.finite(intercept)) {
  stop("CRS intercept is not numeric or not finite.")
}

cat("\nCRS intercept:\n")
print(intercept)

cat("\nCRS coefficient vector summary after extraction:\n")
cat("Number of non-zero coefficient genes:", length(coef_vec), "\n")
print(summary(coef_vec))

cat("\nFirst 20 coefficient genes:\n")
print(head(data.frame(gene = names(coef_vec), coefficient = coef_vec), 20))

available_genes <- intersect(names(coef_vec), rownames(expr_gene))
missing_genes <- setdiff(names(coef_vec), rownames(expr_gene))

cat("\nCRS coefficient genes:", length(coef_vec), "\n")
cat("Available in Lung3:", length(available_genes), "\n")
cat("Missing in Lung3:", length(missing_genes), "\n")
cat("Available fraction:", length(available_genes) / length(coef_vec), "\n")

if (length(available_genes) < 30) {
  warning("CRS gene overlap is low. Interpret Lung3 CRS cautiously.")
}

expr_crs <- expr_gene[available_genes, , drop = FALSE]
coef_use <- coef_vec[available_genes]

# calculate CRS
crs_raw <- as.numeric(intercept + crossprod(coef_use, expr_crs))
names(crs_raw) <- colnames(expr_crs)

lung3_scores <- data.frame(
  GSM = names(crs_raw),
  Lung3_CRS = crs_raw,
  Lung3_CRS_z = as.numeric(scale(crs_raw)),
  stringsAsFactors = FALSE
)

############################################################
# 7. Add sample metadata and clinical
############################################################

if (file.exists(sample_meta_csv)) {
  sample_meta <- read.csv(sample_meta_csv, stringsAsFactors = FALSE, check.names = FALSE)
  lung3_scores <- merge(
    sample_meta,
    lung3_scores,
    by = "GSM",
    all.y = TRUE
  )
}

if (file.exists(clinical_csv)) {
  
  clinical <- read.csv(clinical_csv, stringsAsFactors = FALSE, check.names = FALSE)
  
  if ("lung_number_from_clinical" %in% colnames(clinical) &&
      "lung_number_from_title" %in% colnames(lung3_scores)) {
    
    lung3_scores <- merge(
      lung3_scores,
      clinical,
      by.x = "lung_number_from_title",
      by.y = "lung_number_from_clinical",
      all.x = TRUE
    )
  }
}

cat("\nLung3 CRS score table dim:\n")
print(dim(lung3_scores))

cat("\nLung3 CRS summary:\n")
print(summary(lung3_scores$Lung3_CRS))

cat("\nLung3 CRS_z summary:\n")
print(summary(lung3_scores$Lung3_CRS_z))

############################################################
# 8. Basic clinical correlation with tumor maximum diameter
############################################################

diam_col <- "characteristics.tag.tumor.size.maximumdiameter"

clinical_cor <- data.frame()

if (diam_col %in% colnames(lung3_scores)) {
  
  lung3_scores$tumor_max_diameter <- suppressWarnings(as.numeric(lung3_scores[[diam_col]]))
  
  ok <- is.finite(lung3_scores$Lung3_CRS_z) & is.finite(lung3_scores$tumor_max_diameter)
  
  if (sum(ok) >= 20) {
    
    pear <- cor.test(
      lung3_scores$Lung3_CRS_z[ok],
      lung3_scores$tumor_max_diameter[ok],
      method = "pearson"
    )
    
    spear <- cor.test(
      lung3_scores$Lung3_CRS_z[ok],
      lung3_scores$tumor_max_diameter[ok],
      method = "spearman",
      exact = FALSE
    )
    
    clinical_cor <- data.frame(
      variable = "tumor_max_diameter",
      n = sum(ok),
      Pearson_r = as.numeric(pear$estimate),
      Pearson_p = pear$p.value,
      Spearman_r = as.numeric(spear$estimate),
      Spearman_p = spear$p.value,
      stringsAsFactors = FALSE
    )
  }
}

cat("\n===== Lung3 CRS clinical-size correlation =====\n")
print(clinical_cor)

############################################################
# 9. Audit summary
############################################################

audit_summary <- data.frame(
  item = c(
    "probe_level_expression_rows",
    "expression_samples",
    "GPL15048_annotation_rows_with_symbols",
    "common_probe_count",
    "valid_symbol_mapped_probe_count",
    "unique_gene_symbols_after_mapping",
    "gene_level_expression_rows",
    "CRS_coefficient_genes",
    "CRS_genes_available_in_Lung3",
    "CRS_genes_missing_in_Lung3",
    "CRS_gene_available_fraction",
    "preprocess_note"
  ),
  value = c(
    nrow(expr_mat),
    ncol(expr_mat),
    nrow(annot),
    length(common_probes),
    nrow(expr_use),
    length(unique(annot_use$gene_symbol)),
    nrow(expr_gene),
    length(coef_vec),
    length(available_genes),
    length(missing_genes),
    round(length(available_genes) / length(coef_vec), 4),
    log_note
  ),
  stringsAsFactors = FALSE
)

cat("\n===== Lung3 CRS calculation audit summary =====\n")
print(audit_summary)

interpretation <- data.frame(
  item = c(
    "Lung3 molecular validation feasibility",
    "CRS scoring feasibility",
    "Main limitation",
    "Next recommended step"
  ),
  value = c(
    "Feasible if CRS gene overlap is adequate after GPL15048 probe-to-symbol mapping.",
    paste0(
      "Available CRS genes in Lung3: ",
      length(available_genes),
      "/",
      length(coef_vec),
      " (",
      round(100 * length(available_genes) / length(coef_vec), 1),
      "%)."
    ),
    "Lung3 is microarray-based and cross-platform relative to the RNA-seq training molecular cohort; therefore Lung3 CRS should be treated as exploratory external molecular validation.",
    "Run R/05_pathways/05_evaluate_lung3_replication.R to test whether Lung3_CRS_z reproduces the CRS-associated Hallmark/custom pathway pattern observed in NSCLC-Radiogenomics."
  ),
  stringsAsFactors = FALSE
)

cat("\n===== Lung3 CRS calculation interpretation =====\n")
print(interpretation)

############################################################
# 10. Save outputs
############################################################

saveRDS(
  expr_gene,
  file = file.path(out_dir, "lung3_gene_expression.rds")
)

write.csv(
  best_probe_df,
  file = file.path(out_dir, "lung3_probe_to_gene_mapping.csv"),
  row.names = FALSE
)

write.csv(
  lung3_scores,
  file = file.path(out_dir, "lung3_crs_scores.csv"),
  row.names = FALSE
)

saveRDS(
  lung3_scores,
  file = file.path(out_dir, "lung3_crs_scores.rds")
)

write.csv(
  clinical_cor,
  file = file.path(out_dir, "lung3_crs_clinical_size_association.csv"),
  row.names = FALSE
)

write.csv(
  audit_summary,
  file = file.path(out_dir, "lung3_crs_mapping_summary.csv"),
  row.names = FALSE
)

write.csv(
  interpretation,
  file = file.path(out_dir, "lung3_crs_analysis_summary.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(missing_CRS_gene = missing_genes),
  file = file.path(out_dir, "lung3_missing_crs_genes.csv"),
  row.names = FALSE
)

wb <- createWorkbook()

addWorksheet(wb, "interpretation")
writeData(wb, "interpretation", interpretation)

addWorksheet(wb, "audit_summary")
writeData(wb, "audit_summary", audit_summary)

addWorksheet(wb, "CRS_scores")
writeData(wb, "CRS_scores", lung3_scores)

addWorksheet(wb, "clinical_size_cor")
writeData(wb, "clinical_size_cor", clinical_cor)

addWorksheet(wb, "best_probe_mapping")
writeData(wb, "best_probe_mapping", best_probe_df)

addWorksheet(wb, "missing_CRS_genes")
writeData(wb, "missing_CRS_genes", data.frame(missing_CRS_gene = missing_genes))

for (sheet in names(wb)) {
  setColWidths(wb, sheet, cols = 1:30, widths = "auto")
}

saveWorkbook(
  wb,
  file = file.path(out_dir, "lung3_crs_calculation.xlsx"),
  overwrite = TRUE
)

############################################################
# 11. Done
############################################################

cat("\n===== DONE: Lung3 CRS calculation Lung3 probe-to-symbol mapping and CRS calculation finished =====\n")
cat("Main output:\n")
cat(file.path(out_dir, "lung3_crs_calculation.xlsx"), "\n")

cat("1. GPL15048 platform table dim and columns\n")
cat("2. Clean annotation dim\n")
cat("3. Mapping probes to gene symbols\n")
cat("4. Calculating Lung3 CRS\n")
cat("5. Lung3 CRS summary\n")
cat("6. Lung3 CRS calculation audit summary\n")
cat("7. Lung3 CRS calculation interpretation\n")
