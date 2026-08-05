# 03_prepare_lung3_data.R
#
# Audit Lung3 expression data and sample annotations.
#
# Run from an analysis workspace configured through environment variables.

options(stringsAsFactors = FALSE)
options(timeout = 3600)


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

required_packages <- c("openxlsx", "readxl")
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
library(readxl)

############################################################
# 2. Paths
############################################################


lung3_dir <- file.path(
  project_dir,
  "01_raw_data",
  "Lung3"
)

out_dir <- file.path(
  project_dir,
  "07_results",
  "lung3_data_audit"
)

dir.create(lung3_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

geo_matrix_file <- file.path(
  lung3_dir,
  "GSE58661_series_matrix.txt.gz"
)

geo_matrix_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE58nnn/GSE58661/matrix/GSE58661_series_matrix.txt.gz"

clinical_candidates <- c(
  file.path(lung3_dir, "Lung3.metadata.xls"),
  file.path(lung3_dir, "Lung3.metadata.xlsx"),
  file.path(lung3_dir, "Lung3.csv"),
  file.path(lung3_dir, "Lung3.clinical.csv"),
  file.path(lung3_dir, "Lung3_metadata.csv")
)

crs_symbol_model_rds <- file.path(
  project_dir,
  "06_models",
  "crs_gene_symbol_coefficients.rds"
)

############################################################
# 3. Helpers
############################################################

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

strip_quotes <- function(x) {
  x <- gsub('^"', "", x)
  x <- gsub('"$', "", x)
  x
}

get_matrix_meta_line <- function(lines, prefix) {
  idx <- grep(paste0("^", prefix), lines)
  if (length(idx) == 0) return(NULL)
  vals <- strsplit(lines[idx[1]], "\t")[[1]]
  vals <- vals[-1]
  strip_quotes(vals)
}

guess_id_column <- function(df) {
  candidates <- c(
    "PatientID", "patient_id", "Patient.ID", "Patient",
    "id", "ID", "case", "Case", "name", "Name",
    "Sample", "sample", "sampleid", "SampleID"
  )
  
  hit <- candidates[candidates %in% colnames(df)][1]
  if (!is.na(hit) && length(hit) > 0) return(hit)
  
  # fallback: choose character column with most Lung-like IDs
  scores <- sapply(colnames(df), function(cc) {
    xx <- as.character(df[[cc]])
    sum(grepl("lung|LUNG|Lung", xx), na.rm = TRUE)
  })
  
  if (max(scores, na.rm = TRUE) > 0) {
    return(names(which.max(scores)))
  }
  
  return(colnames(df)[1])
}

read_clinical_auto <- function(path) {
  ext <- tolower(tools::file_ext(path))
  
  if (ext %in% c("xls", "xlsx")) {
    sheets <- readxl::excel_sheets(path)
    out <- readxl::read_excel(path, sheet = sheets[1])
    out <- as.data.frame(out, check.names = FALSE)
    attr(out, "sheet_used") <- sheets[1]
    return(out)
  }
  
  if (ext %in% c("csv", "txt")) {
    out <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    attr(out, "sheet_used") <- NA
    return(out)
  }
  
  stop("Unsupported clinical file type: ", path)
}

############################################################
# 4. Download GEO matrix if needed
############################################################

cat("\n===== Lung3 data audit Lung3 feasibility audit =====\n")

cat("\nLung3 raw data folder:\n")
cat(lung3_dir, "\n")

if (!file.exists(geo_matrix_file) || file.info(geo_matrix_file)$size < 1000000) {
  
  cat("\nGSE58661 series matrix not found or too small. Trying download from GEO FTP...\n")
  
  ok <- tryCatch(
    {
      download.file(
        url = geo_matrix_url,
        destfile = geo_matrix_file,
        mode = "wb",
        quiet = FALSE
      )
      TRUE
    },
    error = function(e) {
      cat("GEO matrix download failed:\n")
      cat(conditionMessage(e), "\n")
      FALSE
    }
  )
  
  if (!ok) {
    cat("\nPlease manually download GSE58661_series_matrix.txt.gz from GEO and place it here:\n")
    cat(geo_matrix_file, "\n")
  }
  
} else {
  cat("\nGSE58661 series matrix already exists. Skip download.\n")
}

############################################################
# 5. Parse GEO matrix
############################################################

geo_audit <- data.frame(
  item = character(),
  value = character(),
  stringsAsFactors = FALSE
)

geo_sample_meta <- NULL
expr_table <- NULL
expr_id_audit <- NULL
geo_ready <- FALSE

if (file.exists(geo_matrix_file) && file.info(geo_matrix_file)$size > 1000000) {
  
  cat("\n===== Reading GSE58661 series matrix =====\n")
  
  lines <- readLines(gzfile(geo_matrix_file), warn = FALSE)
  
  cat("Total lines in series matrix:\n")
  print(length(lines))
  
  sample_title <- get_matrix_meta_line(lines, "!Sample_title")
  sample_geo <- get_matrix_meta_line(lines, "!Sample_geo_accession")
  sample_source <- get_matrix_meta_line(lines, "!Sample_source_name_ch1")
  platform_id <- get_matrix_meta_line(lines, "!Sample_platform_id")
  
  geo_sample_meta <- data.frame(
    GSM = sample_geo,
    sample_title = sample_title,
    sample_source = ifelse(length(sample_source) == length(sample_geo), sample_source, NA),
    platform_id = ifelse(length(platform_id) == length(sample_geo), platform_id, NA),
    lung_number_from_title = extract_lung_number(sample_title),
    stringsAsFactors = FALSE
  )
  
  cat("\nGEO sample meta dim:\n")
  print(dim(geo_sample_meta))
  
  cat("\nFirst 10 GEO samples:\n")
  print(head(geo_sample_meta, 10))
  
  begin_idx <- grep("^!series_matrix_table_begin", lines)
  end_idx <- grep("^!series_matrix_table_end", lines)
  
  if (length(begin_idx) == 1 && length(end_idx) == 1 && end_idx > begin_idx) {
    
    matrix_lines <- lines[(begin_idx + 1):(end_idx - 1)]
    
    expr_table <- read.delim(
      text = paste(matrix_lines, collapse = "\n"),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    cat("\nExpression table dim, probes/genes x samples+ID:\n")
    print(dim(expr_table))
    
    cat("\nExpression table first 5 columns:\n")
    print(colnames(expr_table)[1:min(5, ncol(expr_table))])
    
    expr_ids <- as.character(expr_table[[1]])
    
    common_gene_symbols <- c(
      "TP53", "EGFR", "KRAS", "BRCA1", "BRCA2", "RAD51",
      "MKI67", "CDK1", "CCNB1", "MYC", "HIF1A", "VEGFA",
      "CD274", "STAT1", "IRF1"
    )
    
    expr_id_upper <- toupper(expr_ids)
    
    symbol_like_fraction <- mean(grepl("^[A-Z0-9.-]+$", expr_id_upper))
    common_symbol_hits <- sum(common_gene_symbols %in% expr_id_upper)
    
    crs_overlap_n <- NA
    crs_overlap_fraction <- NA
    
    if (file.exists(crs_symbol_model_rds)) {
      crs_model <- readRDS(crs_symbol_model_rds)
      crs_genes <- toupper(names(crs_model$symbol_coefficients))
      crs_overlap_n <- sum(crs_genes %in% expr_id_upper)
      crs_overlap_fraction <- crs_overlap_n / length(crs_genes)
    }
    
    expr_id_audit <- data.frame(
      item = c(
        "expression_rows",
        "expression_sample_columns",
        "first_id_column_name",
        "example_first_10_ids",
        "symbol_like_fraction",
        "common_known_symbol_hits_15",
        "CRS_symbol_overlap_n",
        "CRS_symbol_overlap_fraction",
        "likely_expression_id_type"
      ),
      value = c(
        nrow(expr_table),
        ncol(expr_table) - 1,
        colnames(expr_table)[1],
        paste(head(expr_ids, 10), collapse = "; "),
        round(symbol_like_fraction, 4),
        common_symbol_hits,
        crs_overlap_n,
        round(crs_overlap_fraction, 4),
        ifelse(
          is.finite(crs_overlap_fraction) && crs_overlap_fraction >= 0.5,
          "Likely gene symbols usable for direct CRS scoring",
          "Likely probe IDs; platform annotation/mapping needed before CRS scoring"
        )
      ),
      stringsAsFactors = FALSE
    )
    
    cat("\n===== Expression ID audit =====\n")
    print(expr_id_audit)
    
    geo_ready <- TRUE
    
  } else {
    cat("\nCannot find series_matrix_table_begin/end. Expression table not parsed.\n")
  }
  
  geo_audit <- rbind(
    geo_audit,
    data.frame(item = "GEO_matrix_file_exists", value = TRUE),
    data.frame(item = "GEO_matrix_file_size_bytes", value = file.info(geo_matrix_file)$size),
    data.frame(item = "GEO_sample_count", value = ifelse(is.null(geo_sample_meta), NA, nrow(geo_sample_meta))),
    data.frame(item = "GEO_platform_ids", value = ifelse(is.null(geo_sample_meta), NA, paste(unique(geo_sample_meta$platform_id), collapse = "; "))),
    data.frame(item = "GEO_expression_table_parsed", value = !is.null(expr_table))
  )
  
} else {
  
  geo_audit <- rbind(
    geo_audit,
    data.frame(item = "GEO_matrix_file_exists", value = FALSE),
    data.frame(item = "GEO_matrix_file_size_bytes", value = NA),
    data.frame(item = "GEO_sample_count", value = NA),
    data.frame(item = "GEO_platform_ids", value = NA),
    data.frame(item = "GEO_expression_table_parsed", value = FALSE)
  )
}

############################################################
# 6. Read Lung3 clinical file if available
############################################################

cat("\n===== Checking Lung3 clinical file =====\n")

clinical_file <- clinical_candidates[file.exists(clinical_candidates)][1]

clinical_df <- NULL
clinical_audit <- data.frame()
clinical_radiomics_candidates <- data.frame()
id_match_audit <- data.frame()

if (!is.na(clinical_file) && length(clinical_file) > 0) {
  
  cat("Clinical file found:\n")
  cat(clinical_file, "\n")
  
  clinical_df <- read_clinical_auto(clinical_file)
  
  cat("\nClinical table dim:\n")
  print(dim(clinical_df))
  
  cat("\nClinical columns:\n")
  print(colnames(clinical_df))
  
  clinical_id_col <- guess_id_column(clinical_df)
  
  cat("\nDetected clinical ID column:\n")
  print(clinical_id_col)
  
  clinical_df$clinical_patient_id_clean <- clean_id(clinical_df[[clinical_id_col]])
  clinical_df$lung_number_from_clinical <- extract_lung_number(clinical_df$clinical_patient_id_clean)
  
  survival_like_cols <- colnames(clinical_df)[
    grepl(
      "surv|dead|death|event|os|time|follow|recurr|progress|status",
      colnames(clinical_df),
      ignore.case = TRUE
    )
  ]
  
  radiomics_like_cols <- colnames(clinical_df)[
    grepl(
      "radiom|feature|texture|shape|firstorder|glcm|glrlm|glszm|gldm|ngtdm|volume|diameter|area|entropy|uniformity|compact|spheric",
      colnames(clinical_df),
      ignore.case = TRUE
    )
  ]
  
  clinical_audit <- data.frame(
    item = c(
      "clinical_file_exists",
      "clinical_file_path",
      "clinical_rows",
      "clinical_columns",
      "clinical_id_col",
      "unique_lung_numbers_in_clinical",
      "survival_like_columns",
      "radiomics_or_size_like_columns_n",
      "radiomics_or_size_like_columns"
    ),
    value = c(
      TRUE,
      clinical_file,
      nrow(clinical_df),
      ncol(clinical_df),
      clinical_id_col,
      length(unique(na.omit(clinical_df$lung_number_from_clinical))),
      paste(survival_like_cols, collapse = "; "),
      length(radiomics_like_cols),
      paste(radiomics_like_cols, collapse = "; ")
    ),
    stringsAsFactors = FALSE
  )
  
  if (length(radiomics_like_cols) > 0) {
    clinical_radiomics_candidates <- data.frame(
      candidate_column = radiomics_like_cols,
      example_values = sapply(radiomics_like_cols, function(cc) {
        paste(head(as.character(clinical_df[[cc]]), 5), collapse = "; ")
      }),
      stringsAsFactors = FALSE
    )
  }
  
  if (!is.null(geo_sample_meta)) {
    
    id_match <- merge(
      geo_sample_meta,
      clinical_df[, c("clinical_patient_id_clean", "lung_number_from_clinical"), drop = FALSE],
      by.x = "lung_number_from_title",
      by.y = "lung_number_from_clinical",
      all = TRUE
    )
    
    id_match_audit <- data.frame(
      item = c(
        "GEO_samples",
        "clinical_rows",
        "matched_by_lung_number",
        "GEO_not_in_clinical",
        "clinical_not_in_GEO"
      ),
      value = c(
        nrow(geo_sample_meta),
        nrow(clinical_df),
        sum(!is.na(id_match$GSM) & !is.na(id_match$clinical_patient_id_clean)),
        sum(!is.na(id_match$GSM) & is.na(id_match$clinical_patient_id_clean)),
        sum(is.na(id_match$GSM) & !is.na(id_match$clinical_patient_id_clean))
      ),
      stringsAsFactors = FALSE
    )
    
    cat("\n===== GEO-clinical ID matching audit =====\n")
    print(id_match_audit)
  }
  
} else {
  
  cat("\nClinical file not found.\n")
  cat("Please manually download Lung3 clinical file from TCIA and put it into:\n")
  cat(lung3_dir, "\n")
  cat("Accepted filenames include: Lung3.metadata.xls, Lung3.metadata.xlsx, Lung3.csv\n")
  
  clinical_audit <- data.frame(
    item = c(
      "clinical_file_exists",
      "clinical_file_path",
      "clinical_rows",
      "clinical_columns",
      "clinical_id_col",
      "unique_lung_numbers_in_clinical",
      "survival_like_columns",
      "radiomics_or_size_like_columns_n",
      "radiomics_or_size_like_columns"
    ),
    value = c(
      FALSE,
      NA,
      NA,
      NA,
      NA,
      NA,
      NA,
      NA,
      NA
    ),
    stringsAsFactors = FALSE
  )
}

############################################################
# 7. Feasibility conclusion
############################################################

feasibility <- data.frame(
  question = character(),
  answer = character(),
  interpretation = character(),
  stringsAsFactors = FALSE
)

feasibility <- rbind(
  feasibility,
  data.frame(
    question = "Can Lung3 support external molecular validation?",
    answer = ifelse(geo_ready, "Potentially yes", "Not yet"),
    interpretation = ifelse(
      geo_ready,
      "GSE58661 expression matrix was parsed. Next step is probe-to-gene-symbol mapping if expression IDs are not already gene symbols.",
      "GSE58661 matrix is missing or could not be parsed."
    )
  )
)

feasibility <- rbind(
  feasibility,
  data.frame(
    question = "Can Lung3 support direct frozen PyRadiomics RRS validation?",
    answer = "Uncertain / high risk",
    interpretation = "TCIA public table lists CT images but not clearly SEG/RTSTRUCT. Direct PyRadiomics replication requires tumor ROI or compatible precomputed radiomics features."
  )
)

feasibility <- rbind(
  feasibility,
  data.frame(
    question = "Can Lung3 support clinical OS validation?",
    answer = "Not primary for this project",
    interpretation = "Lung3 is a surgery-treated cohort, not a radiotherapy cohort. Survival analysis would be exploratory and not radiotherapy-specific."
  )
)

if (!is.null(clinical_df)) {
  
  feasibility <- rbind(
    feasibility,
    data.frame(
      question = "Does the clinical file contain radiomics-like fields?",
      answer = ifelse(nrow(clinical_radiomics_candidates) > 0, "Possibly", "No obvious fields detected"),
      interpretation = ifelse(
        nrow(clinical_radiomics_candidates) > 0,
        "Radiomics/size-like columns were detected in clinical file. Need inspect whether they are compatible with PyRadiomics features.",
        "No obvious radiomics feature columns were detected in clinical file."
      )
    )
  )
}

cat("\n===== Lung3 data audit feasibility conclusion =====\n")
print(feasibility)

############################################################
# 8. Save outputs
############################################################

if (!is.null(geo_sample_meta)) {
  write.csv(
    geo_sample_meta,
    file = file.path(out_dir, "lung3_data_audit_GSE58661_sample_metadata.csv"),
    row.names = FALSE
  )
}

if (!is.null(expr_table)) {
  # only save expression audit, not full matrix as CSV to avoid huge output
  saveRDS(
    expr_table,
    file = file.path(out_dir, "lung3_data_audit_GSE58661_expression_table_probe_by_sample.rds")
  )
}

if (!is.null(expr_id_audit)) {
  write.csv(
    expr_id_audit,
    file = file.path(out_dir, "lung3_data_audit_GSE58661_expression_ID_audit.csv"),
    row.names = FALSE
  )
}

if (!is.null(clinical_df)) {
  write.csv(
    clinical_df,
    file = file.path(out_dir, "lung3_data_audit_Lung3_clinical_clean_initial.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    clinical_df,
    file = file.path(out_dir, "lung3_data_audit_Lung3_clinical_clean_initial.rds")
  )
}

write.csv(
  geo_audit,
  file = file.path(out_dir, "lung3_data_audit_GEO_file_audit.csv"),
  row.names = FALSE
)

write.csv(
  clinical_audit,
  file = file.path(out_dir, "lung3_data_audit_Lung3_clinical_audit.csv"),
  row.names = FALSE
)

write.csv(
  clinical_radiomics_candidates,
  file = file.path(out_dir, "lung3_data_audit_Lung3_clinical_radiomics_like_columns.csv"),
  row.names = FALSE
)

write.csv(
  id_match_audit,
  file = file.path(out_dir, "lung3_data_audit_GEO_clinical_ID_matching_audit.csv"),
  row.names = FALSE
)

write.csv(
  feasibility,
  file = file.path(out_dir, "lung3_data_audit_Lung3_feasibility_conclusion.csv"),
  row.names = FALSE
)

wb <- createWorkbook()

addWorksheet(wb, "feasibility")
writeData(wb, "feasibility", feasibility)

addWorksheet(wb, "GEO_file_audit")
writeData(wb, "GEO_file_audit", geo_audit)

if (!is.null(expr_id_audit)) {
  addWorksheet(wb, "expression_ID_audit")
  writeData(wb, "expression_ID_audit", expr_id_audit)
}

if (!is.null(geo_sample_meta)) {
  addWorksheet(wb, "GEO_sample_meta")
  writeData(wb, "GEO_sample_meta", geo_sample_meta)
}

addWorksheet(wb, "clinical_audit")
writeData(wb, "clinical_audit", clinical_audit)

if (nrow(clinical_radiomics_candidates) > 0) {
  addWorksheet(wb, "radiomics_like_cols")
  writeData(wb, "radiomics_like_cols", clinical_radiomics_candidates)
}

if (nrow(id_match_audit) > 0) {
  addWorksheet(wb, "ID_matching")
  writeData(wb, "ID_matching", id_match_audit)
}

for (sheet in names(wb)) {
  setColWidths(wb, sheet, cols = 1:30, widths = "auto")
}

saveWorkbook(
  wb,
  file = file.path(out_dir, "lung3_data_audit_Lung3_Feasibility_Audit.xlsx"),
  overwrite = TRUE
)

############################################################
# 9. Done
############################################################

cat("\n===== DONE: Lung3 data audit Lung3 feasibility audit finished =====\n")
cat("Main output:\n")
cat(file.path(out_dir, "lung3_data_audit_Lung3_Feasibility_Audit.xlsx"), "\n")

cat("1. GEO sample meta dim and first 10 GEO samples\n")
cat("2. Expression table dim\n")
cat("3. Expression ID audit\n")
cat("4. Clinical table dim and clinical columns, if clinical file exists\n")
cat("5. GEO-clinical ID matching audit\n")
cat("6. Lung3 data audit feasibility conclusion\n")
