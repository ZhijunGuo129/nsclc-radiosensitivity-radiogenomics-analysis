# 05_evaluate_lung3_replication.R
#
# Evaluate cross-cohort replication of the CRS molecular architecture in Lung3.
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

required_packages <- c("openxlsx", "ggplot2")
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
library(ggplot2)

############################################################
# 2. Paths
############################################################


lung3_expr_rds <- file.path(
  project_dir,
  "07_results",
  "lung3_validation",
  "lung3_gene_expression.rds"
)

lung3_score_rds <- file.path(
  project_dir,
  "07_results",
  "lung3_validation",
  "lung3_crs_scores.rds"
)

lung3_score_csv <- file.path(
  project_dir,
  "07_results",
  "lung3_validation",
  "lung3_crs_scores.csv"
)

train_custom_assoc_csv <- file.path(
  project_dir,
  "07_results",
  "molecular_pathways",
  "predefined_module_crs_associations.csv"
)

train_hallmark_assoc_csv <- file.path(
  project_dir,
  "07_results",
  "hallmark_pathways",
  "hallmark_crs_associations.csv"
)

msigdb_dir <- file.path(
  project_dir,
  "01_raw_data",
  "MSigDB"
)

gmt_candidates <- c(
  file.path(msigdb_dir, "h.all.v2026.1.Hs.symbols.gmt"),
  file.path(msigdb_dir, "h.all.v2024.1.Hs.symbols.gmt")
)

gmt_file <- gmt_candidates[file.exists(gmt_candidates)][1]

out_dir <- file.path(
  project_dir,
  "07_results",
  "lung3_validation"
)

fig_dir <- file.path(
  project_dir,
  "08_figures"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(lung3_expr_rds)) {
  stop("Cannot find Lung3 gene-level expression RDS: ", lung3_expr_rds)
}

if (!file.exists(lung3_score_rds) && !file.exists(lung3_score_csv)) {
  stop("Cannot find Lung3 CRS score file.")
}

############################################################
# 3. Helper functions
############################################################

read_gmt <- function(gmt_path) {
  
  lines <- readLines(gmt_path, warn = FALSE)
  lines <- lines[nchar(lines) > 0]
  
  parts <- strsplit(lines, "\t")
  
  set_names <- sapply(parts, function(x) x[1])
  
  gene_sets <- lapply(parts, function(x) {
    genes <- x[-c(1, 2)]
    genes <- toupper(trimws(genes))
    genes <- genes[genes != "" & !is.na(genes)]
    unique(genes)
  })
  
  names(gene_sets) <- set_names
  
  gene_sets
}

score_gene_sets_mean_z <- function(expr_gene_by_sample, gene_sets, min_genes = 5) {
  
  gene_sd <- apply(expr_gene_by_sample, 1, sd, na.rm = TRUE)
  keep_gene <- is.finite(gene_sd) & gene_sd > 0
  
  expr_use <- expr_gene_by_sample[keep_gene, , drop = FALSE]
  
  z_mat <- t(scale(t(expr_use)))
  z_mat[!is.finite(z_mat)] <- NA
  
  pathway_scores <- list()
  pathway_gene_counts <- data.frame()
  
  all_genes <- rownames(z_mat)
  
  for (set_name in names(gene_sets)) {
    
    genes <- unique(toupper(trimws(gene_sets[[set_name]])))
    genes <- genes[genes != "" & !is.na(genes)]
    
    present <- intersect(genes, all_genes)
    
    pathway_gene_counts <- rbind(
      pathway_gene_counts,
      data.frame(
        pathway = set_name,
        input_genes = length(genes),
        present_genes = length(present),
        present_fraction = length(present) / max(length(genes), 1),
        stringsAsFactors = FALSE
      )
    )
    
    if (length(present) >= min_genes) {
      pathway_scores[[set_name]] <- colMeans(
        z_mat[present, , drop = FALSE],
        na.rm = TRUE
      )
    }
  }
  
  score_df <- as.data.frame(pathway_scores, check.names = FALSE)
  score_df$GSM <- colnames(expr_gene_by_sample)
  score_df <- score_df[, c("GSM", setdiff(colnames(score_df), "GSM")), drop = FALSE]
  
  list(
    scores = score_df,
    gene_counts = pathway_gene_counts
  )
}

associate_pathways_with_lung3_CRS <- function(pathway_score_df, score_df, pathway_cols) {
  
  dat <- merge(
    score_df[, c("GSM", "Lung3_CRS", "Lung3_CRS_z", "Lung3_CRS_group")],
    pathway_score_df,
    by = "GSM",
    all = FALSE
  )
  
  out_list <- list()
  
  for (pp in pathway_cols) {
    
    x <- suppressWarnings(as.numeric(dat[[pp]]))
    y <- suppressWarnings(as.numeric(dat$Lung3_CRS_z))
    
    ok <- is.finite(x) & is.finite(y)
    
    if (sum(ok) < 10) next
    
    pear <- suppressWarnings(cor.test(x[ok], y[ok], method = "pearson"))
    spear <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
    
    group_low <- x[dat$Lung3_CRS_group == "CRS_low" & is.finite(x)]
    group_high <- x[dat$Lung3_CRS_group == "CRS_high" & is.finite(x)]
    
    wil <- try(
      wilcox.test(group_high, group_low, exact = FALSE),
      silent = TRUE
    )
    
    wil_p <- ifelse(inherits(wil, "try-error"), NA, wil$p.value)
    
    out_list[[length(out_list) + 1]] <- data.frame(
      pathway = pp,
      n = sum(ok),
      Pearson_r = as.numeric(pear$estimate),
      Pearson_p = pear$p.value,
      Spearman_r = as.numeric(spear$estimate),
      Spearman_p = spear$p.value,
      median_CRS_low = median(group_low, na.rm = TRUE),
      median_CRS_high = median(group_high, na.rm = TRUE),
      median_diff_high_minus_low = median(group_high, na.rm = TRUE) - median(group_low, na.rm = TRUE),
      Wilcoxon_p = wil_p,
      stringsAsFactors = FALSE
    )
  }
  
  out <- do.call(rbind, out_list)
  
  out$Pearson_FDR <- p.adjust(out$Pearson_p, method = "BH")
  out$Spearman_FDR <- p.adjust(out$Spearman_p, method = "BH")
  out$Wilcoxon_FDR <- p.adjust(out$Wilcoxon_p, method = "BH")
  
  out <- out[order(out$Spearman_FDR, -abs(out$Spearman_r)), ]
  
  out
}

compare_with_training <- function(lung_assoc, training_csv, label) {
  
  if (!file.exists(training_csv)) {
    return(list(
      comparison = data.frame(),
      summary = data.frame(
        analysis = label,
        training_file_found = FALSE,
        n_compared = NA,
        n_training_FDR010 = NA,
        n_concordant_among_training_FDR010 = NA,
        concordance_rate_training_FDR010 = NA,
        stringsAsFactors = FALSE
      )
    ))
  }
  
  train_assoc <- read.csv(
    training_csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  train_use <- train_assoc[, c("pathway", "Spearman_r", "Spearman_FDR"), drop = FALSE]
  colnames(train_use) <- c("pathway", "Training_Spearman_r", "Training_Spearman_FDR")
  
  lung_use <- lung_assoc[, c("pathway", "Spearman_r", "Spearman_FDR"), drop = FALSE]
  colnames(lung_use) <- c("pathway", "Lung3_Spearman_r", "Lung3_Spearman_FDR")
  
  comp <- merge(
    train_use,
    lung_use,
    by = "pathway",
    all = FALSE
  )
  
  comp$Training_direction <- ifelse(comp$Training_Spearman_r > 0, "positive", "negative")
  comp$Lung3_direction <- ifelse(comp$Lung3_Spearman_r > 0, "positive", "negative")
  comp$direction_concordant <- comp$Training_direction == comp$Lung3_direction
  
  sig <- comp[is.finite(comp$Training_Spearman_FDR) & comp$Training_Spearman_FDR < 0.10, ]
  
  summary <- data.frame(
    analysis = label,
    training_file_found = TRUE,
    n_compared = nrow(comp),
    n_training_FDR010 = nrow(sig),
    n_concordant_among_training_FDR010 = ifelse(nrow(sig) > 0, sum(sig$direction_concordant), NA),
    concordance_rate_training_FDR010 = ifelse(nrow(sig) > 0, mean(sig$direction_concordant), NA),
    stringsAsFactors = FALSE
  )
  
  list(
    comparison = comp,
    summary = summary
  )
}

############################################################
# 4. Load Lung3 expression and CRS scores
############################################################

cat("\n===== Loading Lung3 gene-level expression and CRS scores =====\n")

expr_gene <- readRDS(lung3_expr_rds)
expr_gene <- as.matrix(expr_gene)

rownames(expr_gene) <- toupper(trimws(rownames(expr_gene)))

cat("Lung3 gene-level expression dim:\n")
print(dim(expr_gene))

cat("\nExpression value summary:\n")
print(summary(as.vector(expr_gene)))

if (file.exists(lung3_score_rds)) {
  lung3_scores <- readRDS(lung3_score_rds)
} else {
  lung3_scores <- read.csv(lung3_score_csv, stringsAsFactors = FALSE, check.names = FALSE)
}

lung3_scores <- as.data.frame(lung3_scores, check.names = FALSE)

if (!"GSM" %in% colnames(lung3_scores)) {
  stop("Lung3 score table does not contain GSM column.")
}

if (!"Lung3_CRS_z" %in% colnames(lung3_scores)) {
  stop("Lung3 score table does not contain Lung3_CRS_z column.")
}

lung3_scores$GSM <- trimws(as.character(lung3_scores$GSM))
lung3_scores$Lung3_CRS <- suppressWarnings(as.numeric(lung3_scores$Lung3_CRS))
lung3_scores$Lung3_CRS_z <- suppressWarnings(as.numeric(lung3_scores$Lung3_CRS_z))

if (all(is.na(lung3_scores$Lung3_CRS_z))) {
  stop("Lung3_CRS_z is all NA. Please rerun Lung3 CRS calculation after coefficient name fix.")
}

common_samples <- intersect(colnames(expr_gene), lung3_scores$GSM)

cat("\nSample matching:\n")
cat("Expression samples:", ncol(expr_gene), "\n")
cat("CRS score samples:", nrow(lung3_scores), "\n")
cat("Common samples:", length(common_samples), "\n")

expr_gene <- expr_gene[, common_samples, drop = FALSE]

lung3_scores <- lung3_scores[match(common_samples, lung3_scores$GSM), , drop = FALSE]

stopifnot(all(lung3_scores$GSM == colnames(expr_gene)))

lung3_scores$Lung3_CRS_group <- ifelse(
  lung3_scores$Lung3_CRS_z >= median(lung3_scores$Lung3_CRS_z, na.rm = TRUE),
  "CRS_high",
  "CRS_low"
)

lung3_scores$Lung3_CRS_group <- factor(
  lung3_scores$Lung3_CRS_group,
  levels = c("CRS_low", "CRS_high")
)

cat("\nLung3 CRS_z summary:\n")
print(summary(lung3_scores$Lung3_CRS_z))

cat("\nLung3 CRS group table:\n")
print(table(lung3_scores$Lung3_CRS_group, useNA = "ifany"))

############################################################
# 5. Custom radiosensitivity-related modules
############################################################

cat("\n===== Custom radiosensitivity-related modules in Lung3 =====\n")

custom_gene_sets <- list(
  
  DNA_DAMAGE_REPAIR = c(
    "ATM", "ATR", "CHEK1", "CHEK2", "TP53", "BRCA1", "BRCA2",
    "RAD51", "RAD50", "MRE11", "NBN", "XRCC5", "XRCC6",
    "PRKDC", "PARP1", "ERCC1", "FANCD2", "RPA1", "RPA2"
  ),
  
  HOMOLOGOUS_RECOMBINATION = c(
    "BRCA1", "BRCA2", "RAD51", "RAD51C", "RAD51D",
    "PALB2", "BARD1", "BRIP1", "XRCC2", "XRCC3",
    "MRE11", "RAD50", "NBN"
  ),
  
  NON_HOMOLOGOUS_END_JOINING = c(
    "XRCC5", "XRCC6", "PRKDC", "LIG4", "XRCC4", "NHEJ1", "DCLRE1C"
  ),
  
  CELL_CYCLE_G2M_CORE = c(
    "CDK1", "CCNB1", "CCNB2", "CDC25C", "CDC20",
    "AURKA", "AURKB", "BUB1", "BUB1B", "MAD2L1",
    "PLK1", "TOP2A", "MKI67"
  ),
  
  PROLIFERATION_CORE = c(
    "MKI67", "PCNA", "TOP2A", "MCM2", "MCM3",
    "MCM4", "MCM5", "MCM6", "MCM7", "CDK1", "CCNA2", "CCNB1"
  ),
  
  APOPTOSIS_CORE = c(
    "BAX", "BAK1", "BCL2", "BCL2L1", "CASP3", "CASP7",
    "CASP8", "CASP9", "FAS", "FASLG", "BBC3", "PMAIP1"
  ),
  
  EMT_CORE = c(
    "VIM", "CDH2", "CDH1", "SNAI1", "SNAI2", "TWIST1",
    "ZEB1", "ZEB2", "FN1", "ITGA5", "MMP2", "MMP9"
  ),
  
  HYPOXIA_CORE10 = c(
    "HIF1A", "VEGFA", "CA9", "SLC2A1", "LDHA",
    "PDK1", "BNIP3", "NDRG1", "EGLN1", "EGLN3"
  ),
  
  INTERFERON_GAMMA_CORE = c(
    "IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10",
    "IDO1", "GBP1", "HLA-A", "HLA-B", "HLA-C"
  ),
  
  IMMUNE_CHECKPOINT_CORE = c(
    "CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3",
    "TIGIT", "HAVCR2", "IDO1", "CD80", "CD86"
  )
)

custom_scored <- score_gene_sets_mean_z(
  expr_gene_by_sample = expr_gene,
  gene_sets = custom_gene_sets,
  min_genes = 5
)

custom_scores <- custom_scored$scores
custom_gene_counts <- custom_scored$gene_counts
custom_cols <- setdiff(colnames(custom_scores), "GSM")

custom_assoc <- associate_pathways_with_lung3_CRS(
  pathway_score_df = custom_scores,
  score_df = lung3_scores,
  pathway_cols = custom_cols
)

cat("\n===== Lung3 custom module gene counts =====\n")
print(custom_gene_counts)

cat("\n===== Lung3 custom module association with CRS_z =====\n")
print(custom_assoc)

############################################################
# 6. Hallmark pathways
############################################################

hallmark_available <- FALSE
hallmark_scores <- data.frame()
hallmark_gene_counts <- data.frame()
hallmark_assoc <- data.frame()

if (!is.na(gmt_file) && length(gmt_file) > 0 && file.exists(gmt_file)) {
  
  cat("\n===== Hallmark GMT in Lung3 =====\n")
  cat("GMT file used:\n")
  cat(gmt_file, "\n")
  
  hallmark_gene_sets <- read_gmt(gmt_file)
  
  cat("Number of Hallmark gene sets:\n")
  print(length(hallmark_gene_sets))
  
  hallmark_scored <- score_gene_sets_mean_z(
    expr_gene_by_sample = expr_gene,
    gene_sets = hallmark_gene_sets,
    min_genes = 10
  )
  
  hallmark_scores <- hallmark_scored$scores
  hallmark_gene_counts <- hallmark_scored$gene_counts
  hallmark_cols <- setdiff(colnames(hallmark_scores), "GSM")
  
  hallmark_assoc <- associate_pathways_with_lung3_CRS(
    pathway_score_df = hallmark_scores,
    score_df = lung3_scores,
    pathway_cols = hallmark_cols
  )
  
  hallmark_available <- TRUE
  
  cat("\n===== Lung3 Hallmark association with CRS_z: top positive Spearman =====\n")
  print(
    head(
      hallmark_assoc[order(-hallmark_assoc$Spearman_r), ],
      15
    )
  )
  
  cat("\n===== Lung3 Hallmark association with CRS_z: top negative Spearman =====\n")
  print(
    head(
      hallmark_assoc[order(hallmark_assoc$Spearman_r), ],
      15
    )
  )
  
  cat("\n===== Lung3 Hallmark significant pathways by Spearman_FDR < 0.10 =====\n")
  print(
    hallmark_assoc[hallmark_assoc$Spearman_FDR < 0.10, ]
  )
  
} else {
  cat("\nHallmark GMT file not found. Hallmark Lung3 analysis skipped.\n")
}

############################################################
# 7. Compare Lung3 with NSCLC-Radiogenomics directions
############################################################

cat("\n===== Comparing Lung3 with NSCLC-Radiogenomics directions =====\n")

custom_comp <- compare_with_training(
  lung_assoc = custom_assoc,
  training_csv = train_custom_assoc_csv,
  label = "Custom_modules"
)

cat("\nCustom module direction comparison summary:\n")
print(custom_comp$summary)

if (nrow(custom_comp$comparison) > 0) {
  cat("\nCustom module direction comparison table:\n")
  print(custom_comp$comparison)
}

if (hallmark_available) {
  
  hallmark_comp <- compare_with_training(
    lung_assoc = hallmark_assoc,
    training_csv = train_hallmark_assoc_csv,
    label = "Hallmark_pathways"
  )
  
  cat("\nHallmark direction comparison summary:\n")
  print(hallmark_comp$summary)
  
  cat("\nHallmark direction comparison table for training FDR < 0.10:\n")
  print(
    hallmark_comp$comparison[
      hallmark_comp$comparison$Training_Spearman_FDR < 0.10,
    ]
  )
  
} else {
  
  hallmark_comp <- list(
    comparison = data.frame(),
    summary = data.frame()
  )
}

############################################################
# 8. Key pathway table
############################################################

key_custom <- c(
  "CELL_CYCLE_G2M_CORE",
  "PROLIFERATION_CORE",
  "DNA_DAMAGE_REPAIR",
  "HOMOLOGOUS_RECOMBINATION",
  "HYPOXIA_CORE10",
  "INTERFERON_GAMMA_CORE",
  "EMT_CORE",
  "APOPTOSIS_CORE",
  "IMMUNE_CHECKPOINT_CORE"
)

key_hallmark <- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_UV_RESPONSE_DN",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_NOTCH_SIGNALING",
  "HALLMARK_HEDGEHOG_SIGNALING"
)

key_custom_table <- custom_assoc[custom_assoc$pathway %in% key_custom, ]
key_custom_table$set_type <- "Custom"

key_hallmark_table <- data.frame()

if (hallmark_available) {
  key_hallmark_table <- hallmark_assoc[hallmark_assoc$pathway %in% key_hallmark, ]
  key_hallmark_table$set_type <- "Hallmark"
}

key_pathway_table <- rbind(
  key_custom_table,
  key_hallmark_table
)

key_pathway_table <- key_pathway_table[
  order(key_pathway_table$set_type, key_pathway_table$Spearman_r),
]

cat("\n===== Key Lung3 CRS-associated pathways =====\n")
print(key_pathway_table)

############################################################
# 9. Analysis summary
############################################################

custom_sig <- custom_assoc[custom_assoc$Spearman_FDR < 0.10, ]

hallmark_sig <- data.frame()
if (hallmark_available) {
  hallmark_sig <- hallmark_assoc[hallmark_assoc$Spearman_FDR < 0.10, ]
}

custom_comp_sum <- custom_comp$summary
hallmark_comp_sum <- hallmark_comp$summary

interpretation_rows <- list()

interpretation_rows[[1]] <- data.frame(
  item = "Lung3 CRS scoring feasibility",
  value = paste0(
    "Lung3_CRS_z was calculated in ",
    nrow(lung3_scores),
    " patients after 91/92 CRS coefficient genes were available in the Lung3 gene-symbol expression matrix."
  ),
  stringsAsFactors = FALSE
)

interpretation_rows[[2]] <- data.frame(
  item = "Custom modules significant in Lung3",
  value = ifelse(
    nrow(custom_sig) > 0,
    paste(
      paste0(
        custom_sig$pathway,
        " (rho=",
        round(custom_sig$Spearman_r, 3),
        ", FDR=",
        signif(custom_sig$Spearman_FDR, 3),
        ")"
      ),
      collapse = "; "
    ),
    "None at Spearman_FDR < 0.10"
  ),
  stringsAsFactors = FALSE
)

interpretation_rows[[3]] <- data.frame(
  item = "Hallmark pathways significant in Lung3",
  value = ifelse(
    hallmark_available && nrow(hallmark_sig) > 0,
    paste(
      paste0(
        hallmark_sig$pathway,
        " (rho=",
        round(hallmark_sig$Spearman_r, 3),
        ", FDR=",
        signif(hallmark_sig$Spearman_FDR, 3),
        ")"
      ),
      collapse = "; "
    ),
    "None at Spearman_FDR < 0.10 or Hallmark not available"
  ),
  stringsAsFactors = FALSE
)

interpretation_rows[[4]] <- data.frame(
  item = "Direction concordance with NSCLC-Radiogenomics custom modules",
  value = ifelse(
    nrow(custom_comp_sum) > 0 && isTRUE(custom_comp_sum$training_file_found[1]),
    paste0(
      custom_comp_sum$n_concordant_among_training_FDR010[1],
      "/",
      custom_comp_sum$n_training_FDR010[1],
      " training-significant custom modules showed concordant direction in Lung3."
    ),
    "Training custom association file not found."
  ),
  stringsAsFactors = FALSE
)

interpretation_rows[[5]] <- data.frame(
  item = "Direction concordance with NSCLC-Radiogenomics Hallmark pathways",
  value = ifelse(
    nrow(hallmark_comp_sum) > 0 && isTRUE(hallmark_comp_sum$training_file_found[1]),
    paste0(
      hallmark_comp_sum$n_concordant_among_training_FDR010[1],
      "/",
      hallmark_comp_sum$n_training_FDR010[1],
      " training-significant Hallmark pathways showed concordant direction in Lung3."
    ),
    "Training Hallmark association file not found or Hallmark not available."
  ),
  stringsAsFactors = FALSE
)

interpretation_rows[[6]] <- data.frame(
  item = "Important limitation",
  value = "Lung3 is a microarray-based surgical NSCLC cohort and is cross-platform relative to NSCLC-Radiogenomics RNA-seq; therefore this should be interpreted as exploratory external molecular validation, not clinical radiotherapy validation.",
  stringsAsFactors = FALSE
)

interpretation_df <- do.call(rbind, interpretation_rows)

cat("\n===== Lung3 validation analysis summary =====\n")
print(interpretation_df)

############################################################
# 10. Save outputs
############################################################

write.csv(
  lung3_scores,
  file = file.path(out_dir, "lung3_crs_scores_used.csv"),
  row.names = FALSE
)

write.csv(
  custom_scores,
  file = file.path(out_dir, "lung3_custom_module_scores.csv"),
  row.names = FALSE
)

write.csv(
  custom_gene_counts,
  file = file.path(out_dir, "lung3_custom_module_gene_counts.csv"),
  row.names = FALSE
)

write.csv(
  custom_assoc,
  file = file.path(out_dir, "lung3_predefined_module_crs_associations.csv"),
  row.names = FALSE
)

write.csv(
  custom_comp$comparison,
  file = file.path(out_dir, "predefined_module_replication.csv"),
  row.names = FALSE
)

if (hallmark_available) {
  write.csv(
    hallmark_scores,
    file = file.path(out_dir, "lung3_hallmark_scores.csv"),
    row.names = FALSE
  )
  
  write.csv(
    hallmark_gene_counts,
    file = file.path(out_dir, "lung3_hallmark_gene_counts.csv"),
    row.names = FALSE
  )
  
  write.csv(
    hallmark_assoc,
    file = file.path(out_dir, "lung3_hallmark_crs_association.csv"),
    row.names = FALSE
  )
  
  write.csv(
    hallmark_comp$comparison,
    file = file.path(out_dir, "hallmark_replication.csv"),
    row.names = FALSE
  )
}

write.csv(
  key_pathway_table,
  file = file.path(out_dir, "lung3_key_pathway_results.csv"),
  row.names = FALSE
)

write.csv(
  interpretation_df,
  file = file.path(out_dir, "lung3_analysis_summary.csv"),
  row.names = FALSE
)

comparison_summary <- rbind(
  custom_comp$summary,
  hallmark_comp$summary
)

write.csv(
  comparison_summary,
  file = file.path(out_dir, "lung3_direction_concordance_summary.csv"),
  row.names = FALSE
)

############################################################
# 11. Excel workbook
############################################################

wb <- createWorkbook()

addWorksheet(wb, "interpretation")
writeData(wb, "interpretation", interpretation_df)

addWorksheet(wb, "direction_summary")
writeData(wb, "direction_summary", comparison_summary)

addWorksheet(wb, "key_pathways")
writeData(wb, "key_pathways", key_pathway_table)

addWorksheet(wb, "custom_assoc")
writeData(wb, "custom_assoc", custom_assoc)

addWorksheet(wb, "custom_gene_counts")
writeData(wb, "custom_gene_counts", custom_gene_counts)

addWorksheet(wb, "custom_direction_compare")
writeData(wb, "custom_direction_compare", custom_comp$comparison)

if (hallmark_available) {
  addWorksheet(wb, "hallmark_assoc")
  writeData(wb, "hallmark_assoc", hallmark_assoc)
  
  addWorksheet(wb, "hallmark_gene_counts")
  writeData(wb, "hallmark_gene_counts", hallmark_gene_counts)
  
  addWorksheet(wb, "hallmark_direction_compare")
  writeData(wb, "hallmark_direction_compare", hallmark_comp$comparison)
}

addWorksheet(wb, "Lung3_CRS_scores")
writeData(wb, "Lung3_CRS_scores", lung3_scores)

for (sheet in names(wb)) {
  setColWidths(wb, sheet, cols = 1:30, widths = "auto")
}

saveWorkbook(
  wb,
  file = file.path(out_dir, "lung3_molecular_replication.xlsx"),
  overwrite = TRUE
)

############################################################
# 12. Figures
############################################################

custom_plot_df <- custom_assoc
custom_plot_df$pathway <- factor(
  custom_plot_df$pathway,
  levels = custom_plot_df$pathway[order(custom_plot_df$Spearman_r)]
)

p_custom <- ggplot(
  custom_plot_df,
  aes(x = pathway, y = Spearman_r)
) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(
    title = "Lung3: custom modules associated with CRS_z",
    x = "",
    y = "Spearman correlation with Lung3_CRS_z"
  )

ggsave(
  filename = file.path(fig_dir, "lung3_custom_module_crs_association.png"),
  plot = p_custom,
  width = 8,
  height = 5,
  dpi = 300
)

if (hallmark_available) {
  
  hallmark_plot_df <- rbind(
    head(hallmark_assoc[order(-hallmark_assoc$Spearman_r), ], 10),
    head(hallmark_assoc[order(hallmark_assoc$Spearman_r), ], 10)
  )
  
  hallmark_plot_df <- hallmark_plot_df[!duplicated(hallmark_plot_df$pathway), ]
  
  hallmark_plot_df$pathway_clean <- gsub("^HALLMARK_", "", hallmark_plot_df$pathway)
  hallmark_plot_df$pathway_clean <- gsub("_", " ", hallmark_plot_df$pathway_clean)
  
  hallmark_plot_df$pathway_clean <- factor(
    hallmark_plot_df$pathway_clean,
    levels = hallmark_plot_df$pathway_clean[order(hallmark_plot_df$Spearman_r)]
  )
  
  p_hallmark <- ggplot(
    hallmark_plot_df,
    aes(x = pathway_clean, y = Spearman_r)
  ) +
    geom_col() +
    coord_flip() +
    theme_bw(base_size = 11) +
    labs(
      title = "Lung3: top Hallmark pathways associated with CRS_z",
      x = "",
      y = "Spearman correlation with Lung3_CRS_z"
    )
  
  ggsave(
    filename = file.path(fig_dir, "lung3_top_hallmark_crs_association.png"),
    plot = p_hallmark,
    width = 8,
    height = 6,
    dpi = 300
  )
}

if (nrow(custom_comp$comparison) > 0) {
  
  p_comp_custom <- ggplot(
    custom_comp$comparison,
    aes(x = Training_Spearman_r, y = Lung3_Spearman_r)
  ) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_bw(base_size = 11) +
    labs(
      title = "Direction comparison: custom modules",
      x = "NSCLC-Radiogenomics Spearman rho",
      y = "Lung3 Spearman rho"
    )
  
  ggsave(
    filename = file.path(fig_dir, "custom_module_effect_comparison.png"),
    plot = p_comp_custom,
    width = 6,
    height = 5,
    dpi = 300
  )
}

if (hallmark_available && nrow(hallmark_comp$comparison) > 0) {
  
  p_comp_hallmark <- ggplot(
    hallmark_comp$comparison,
    aes(x = Training_Spearman_r, y = Lung3_Spearman_r)
  ) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_bw(base_size = 11) +
    labs(
      title = "Direction comparison: Hallmark pathways",
      x = "NSCLC-Radiogenomics Spearman rho",
      y = "Lung3 Spearman rho"
    )
  
  ggsave(
    filename = file.path(fig_dir, "hallmark_effect_comparison.png"),
    plot = p_comp_hallmark,
    width = 6,
    height = 5,
    dpi = 300
  )
}

############################################################
# 13. Key results
############################################################

txt_lines <- c(
  "Lung3 molecular validation summary",
  "",
  "Purpose:",
  "This analysis evaluates whether the CRS-associated molecular pathway pattern discovered in NSCLC-Radiogenomics is reproduced in the independent Lung3/GSE58661 microarray cohort.",
  "",
  "Analysis summary:",
  paste(capture.output(print(interpretation_df)), collapse = "\n"),
  "",
  "Recommended cautious wording if direction concordance is good:",
  "In the independent Lung3 microarray cohort, CRS could be calculated using 91 of 92 CRS genes. Pathway-level analysis showed partial concordance with the NSCLC-Radiogenomics discovery cohort, supporting the cross-platform molecular relevance of CRS. Given the microarray platform and surgical nature of Lung3, this result should be interpreted as exploratory external molecular validation rather than radiotherapy outcome validation.",
  "",
  "Recommended cautious wording if direction concordance is weak:",
  "Although CRS could be calculated in Lung3 using 91 of 92 CRS genes, pathway-level associations were only partially reproduced, suggesting cross-platform and cohort-context dependence. Therefore, Lung3 was considered an exploratory molecular audit rather than definitive external validation."
)

writeLines(
  txt_lines,
  con = file.path(out_dir, "lung3_key_results.txt"),
  useBytes = TRUE
)

############################################################
# 14. Done
############################################################

cat("\n===== DONE: Lung3 validation Lung3 external molecular validation finished =====\n")
cat("Main output:\n")
cat(file.path(out_dir, "lung3_molecular_replication.xlsx"), "\n")
cat(file.path(out_dir, "lung3_analysis_summary.csv"), "\n")
cat(file.path(out_dir, "lung3_direction_concordance_summary.csv"), "\n")
cat(file.path(fig_dir, "lung3_custom_module_crs_association.png"), "\n")
cat(file.path(fig_dir, "lung3_top_hallmark_crs_association.png"), "\n")

cat("1. Lung3 custom module association with CRS_z\n")
cat("2. Lung3 Hallmark association with CRS_z: top positive Spearman\n")
cat("3. Lung3 Hallmark association with CRS_z: top negative Spearman\n")
cat("4. Lung3 Hallmark significant pathways by Spearman_FDR < 0.10\n")
cat("5. Custom module direction comparison summary\n")
cat("6. Hallmark direction comparison summary\n")
cat("7. Key Lung3 CRS-associated pathways\n")
cat("8. Lung3 validation analysis summary\n")
