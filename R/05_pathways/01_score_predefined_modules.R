# 01_score_predefined_modules.R
#
# Analyze predefined radiosensitivity-related molecular modules.
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


expr_rds <- file.path(
  project_dir,
  "01_raw_data",
  "NSCLC_Radiogenomics",
  "GSE103584_NSCLC_RNAseq_SYMBOL_gene_by_sample_LOG2_PROCESSED.rds"
)

score_rds <- file.path(
  project_dir,
  "05_molecular_scores",
  "NSCLC_Radiogenomics_patient_molecular_scores.rds"
)

score_csv <- file.path(
  project_dir,
  "05_molecular_scores",
  "NSCLC_Radiogenomics_patient_molecular_scores.csv"
)

out_dir <- file.path(
  project_dir,
  "07_results",
  "molecular_pathways"
)

fig_dir <- file.path(
  project_dir,
  "08_figures"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(expr_rds)) {
  stop("Cannot find processed expression RDS: ", expr_rds)
}

if (!file.exists(score_rds) && !file.exists(score_csv)) {
  stop("Cannot find molecular score table RDS or CSV.")
}

############################################################
# 3. Load data
############################################################

cat("\n===== Loading processed patient expression =====\n")

expr <- readRDS(expr_rds)

expr <- as.data.frame(expr, check.names = FALSE)

cat("Expression dim, genes x samples:\n")
print(dim(expr))

cat("\nFirst 10 gene IDs:\n")
print(head(rownames(expr), 10))

cat("\nFirst 10 sample IDs:\n")
print(head(colnames(expr), 10))

cat("\n===== Loading molecular scores =====\n")

if (file.exists(score_rds)) {
  scores <- readRDS(score_rds)
} else {
  scores <- read.csv(score_csv, stringsAsFactors = FALSE, check.names = FALSE)
}

scores <- as.data.frame(scores, check.names = FALSE)

cat("Score table dim:\n")
print(dim(scores))

cat("\nScore table columns:\n")
print(colnames(scores))

############################################################
# 4. Standardize expression gene symbols and sample IDs
############################################################

# gene symbols
rownames(expr) <- toupper(trimws(as.character(rownames(expr))))

# remove empty gene names
expr <- expr[rownames(expr) != "" & !is.na(rownames(expr)), , drop = FALSE]

# convert to numeric
for (cc in colnames(expr)) {
  expr[[cc]] <- suppressWarnings(as.numeric(expr[[cc]]))
}

# aggregate duplicated gene symbols by mean
if (any(duplicated(rownames(expr)))) {
  cat("\nDuplicated gene symbols detected. Aggregating by mean...\n")
  
  expr$gene_symbol_tmp <- rownames(expr)
  
  expr <- aggregate(
    . ~ gene_symbol_tmp,
    data = expr,
    FUN = function(x) mean(as.numeric(x), na.rm = TRUE)
  )
  
  rownames(expr) <- expr$gene_symbol_tmp
  expr$gene_symbol_tmp <- NULL
}

expr_mat <- as.matrix(expr)

cat("\nExpression dim after gene aggregation:\n")
print(dim(expr_mat))

############################################################
# 5. Standardize score sample ID
############################################################

candidate_id_cols <- c(
  "patient_id",
  "sampleid",
  "sample_id",
  "SampleID",
  "PatientID",
  "patient_id_molecular"
)

id_col <- candidate_id_cols[candidate_id_cols %in% colnames(scores)][1]

if (is.na(id_col) || length(id_col) == 0) {
  stop("Cannot detect sample ID column in molecular score table.")
}

scores$patient_id <- toupper(trimws(as.character(scores[[id_col]])))

if (!"CRS_z" %in% colnames(scores)) {
  stop("Cannot find CRS_z in molecular score table.")
}

if (!"CRS" %in% colnames(scores)) {
  stop("Cannot find CRS in molecular score table.")
}

scores$CRS <- suppressWarnings(as.numeric(scores$CRS))
scores$CRS_z <- suppressWarnings(as.numeric(scores$CRS_z))

############################################################
# 6. Match expression and scores
############################################################

common_samples <- intersect(colnames(expr_mat), scores$patient_id)

cat("\n===== Sample matching =====\n")
cat("Expression samples:", ncol(expr_mat), "\n")
cat("Score samples:", nrow(scores), "\n")
cat("Common samples:", length(common_samples), "\n")

if (length(common_samples) < 50) {
  stop("Too few matched samples.")
}

expr_mat <- expr_mat[, common_samples, drop = FALSE]

scores2 <- scores[match(common_samples, scores$patient_id), , drop = FALSE]

stopifnot(all(scores2$patient_id == colnames(expr_mat)))

# CRS median group
crs_median <- median(scores2$CRS_z, na.rm = TRUE)

scores2$CRS_group <- ifelse(
  scores2$CRS_z >= crs_median,
  "CRS_high",
  "CRS_low"
)

scores2$CRS_group <- factor(
  scores2$CRS_group,
  levels = c("CRS_low", "CRS_high")
)

cat("\nCRS_z summary:\n")
print(summary(scores2$CRS_z))

cat("\nCRS group table:\n")
print(table(scores2$CRS_group, useNA = "ifany"))

############################################################
# 7. Pathway scoring function
############################################################

score_gene_sets_mean_z <- function(expr_gene_by_sample, gene_sets, min_genes = 5) {
  
  # z-score each gene across samples
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
  score_df$patient_id <- colnames(expr_gene_by_sample)
  score_df <- score_df[, c("patient_id", setdiff(colnames(score_df), "patient_id")), drop = FALSE]
  
  list(
    scores = score_df,
    gene_counts = pathway_gene_counts
  )
}

associate_pathways_with_CRS <- function(pathway_score_df, scores_df, pathway_cols) {
  
  out_list <- list()
  
  dat <- merge(
    scores_df[, c("patient_id", "CRS", "CRS_z", "CRS_group")],
    pathway_score_df,
    by = "patient_id",
    all = FALSE
  )
  
  for (pp in pathway_cols) {
    
    x <- suppressWarnings(as.numeric(dat[[pp]]))
    y <- suppressWarnings(as.numeric(dat$CRS_z))
    
    ok <- is.finite(x) & is.finite(y)
    
    if (sum(ok) < 10) next
    
    pear <- suppressWarnings(cor.test(x[ok], y[ok], method = "pearson"))
    spear <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
    
    group_low <- x[dat$CRS_group == "CRS_low" & is.finite(x)]
    group_high <- x[dat$CRS_group == "CRS_high" & is.finite(x)]
    
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
  
  out <- out[
    order(out$Spearman_FDR, -abs(out$Spearman_r)),
  ]
  
  out
}

############################################################
# 8. Hallmark gene sets from msigdbr
############################################################

hallmark_available <- FALSE
hallmark_gene_sets <- list()
hallmark_scores <- NULL
hallmark_gene_counts <- NULL
hallmark_assoc <- NULL

cat("\n===== Loading Hallmark gene sets from msigdbr =====\n")

if (!requireNamespace("msigdbr", quietly = TRUE)) {
  stop("Missing R package 'msigdbr'. Install dependencies with environment/install_r_dependencies.R.")
}

if (requireNamespace("msigdbr", quietly = TRUE)) {
  
  hallmark_df <- tryCatch(
    {
      msigdbr::msigdbr(species = "Homo sapiens", category = "H")
    },
    error = function(e1) {
      tryCatch(
        {
          msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
        },
        error = function(e2) {
          NULL
        }
      )
    }
  )
  
  if (!is.null(hallmark_df)) {
    
    if (!"gene_symbol" %in% colnames(hallmark_df)) {
      stop("msigdbr result does not contain gene_symbol column.")
    }
    
    if (!"gs_name" %in% colnames(hallmark_df)) {
      stop("msigdbr result does not contain gs_name column.")
    }
    
    hallmark_df$gene_symbol <- toupper(trimws(hallmark_df$gene_symbol))
    
    hallmark_gene_sets <- split(
      hallmark_df$gene_symbol,
      hallmark_df$gs_name
    )
    
    cat("Hallmark gene sets loaded:\n")
    print(length(hallmark_gene_sets))
    
    hallmark_scored <- score_gene_sets_mean_z(
      expr_gene_by_sample = expr_mat,
      gene_sets = hallmark_gene_sets,
      min_genes = 10
    )
    
    hallmark_scores <- hallmark_scored$scores
    hallmark_gene_counts <- hallmark_scored$gene_counts
    
    hallmark_pathway_cols <- setdiff(colnames(hallmark_scores), "patient_id")
    
    hallmark_assoc <- associate_pathways_with_CRS(
      pathway_score_df = hallmark_scores,
      scores_df = scores2,
      pathway_cols = hallmark_pathway_cols
    )
    
    hallmark_available <- TRUE
    
    cat("\n===== Hallmark association with CRS_z: top positive Spearman =====\n")
    print(
      head(
        hallmark_assoc[order(-hallmark_assoc$Spearman_r), ],
        15
      )
    )
    
    cat("\n===== Hallmark association with CRS_z: top negative Spearman =====\n")
    print(
      head(
        hallmark_assoc[order(hallmark_assoc$Spearman_r), ],
        15
      )
    )
    
    cat("\n===== Hallmark significant pathways by Spearman_FDR < 0.10 =====\n")
    print(
      hallmark_assoc[hallmark_assoc$Spearman_FDR < 0.10, ]
    )
    
  } else {
    cat("msigdbr failed. Hallmark analysis skipped.\n")
  }
  
} else {
  cat("msigdbr not available. Hallmark analysis skipped.\n")
}

############################################################
# 9. Custom radiosensitivity-related gene modules
############################################################

cat("\n===== Custom radiosensitivity-related modules =====\n")

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
  expr_gene_by_sample = expr_mat,
  gene_sets = custom_gene_sets,
  min_genes = 5
)

custom_scores <- custom_scored$scores
custom_gene_counts <- custom_scored$gene_counts

custom_pathway_cols <- setdiff(colnames(custom_scores), "patient_id")

custom_assoc <- associate_pathways_with_CRS(
  pathway_score_df = custom_scores,
  scores_df = scores2,
  pathway_cols = custom_pathway_cols
)

cat("\n===== Custom module gene counts =====\n")
print(custom_gene_counts)

cat("\n===== Custom module association with CRS_z =====\n")
print(custom_assoc)

############################################################
# 10. Merge pathway scores with CRS table
############################################################

all_pathway_scores <- scores2[, c(
  "patient_id",
  "CRS",
  "CRS_z",
  "CRS_group"
)]

if (hallmark_available) {
  all_pathway_scores <- merge(
    all_pathway_scores,
    hallmark_scores,
    by = "patient_id",
    all.x = TRUE
  )
}

all_pathway_scores <- merge(
  all_pathway_scores,
  custom_scores,
  by = "patient_id",
  all.x = TRUE
)

############################################################
# 11. Manuscript-oriented analysis summary
############################################################

interpretation_rows <- list()

if (hallmark_available) {
  
  top_pos <- head(hallmark_assoc[order(-hallmark_assoc$Spearman_r), ], 10)
  top_neg <- head(hallmark_assoc[order(hallmark_assoc$Spearman_r), ], 10)
  
  interpretation_rows[[length(interpretation_rows) + 1]] <- data.frame(
    item = "Top positively CRS-associated Hallmark pathways",
    value = paste(
      paste0(top_pos$pathway, " (rho=", round(top_pos$Spearman_r, 3), ", FDR=", signif(top_pos$Spearman_FDR, 3), ")"),
      collapse = "; "
    ),
    stringsAsFactors = FALSE
  )
  
  interpretation_rows[[length(interpretation_rows) + 1]] <- data.frame(
    item = "Top negatively CRS-associated Hallmark pathways",
    value = paste(
      paste0(top_neg$pathway, " (rho=", round(top_neg$Spearman_r, 3), ", FDR=", signif(top_neg$Spearman_FDR, 3), ")"),
      collapse = "; "
    ),
    stringsAsFactors = FALSE
  )
  
  sig_hallmark <- hallmark_assoc[hallmark_assoc$Spearman_FDR < 0.10, ]
  
  interpretation_rows[[length(interpretation_rows) + 1]] <- data.frame(
    item = "Hallmark pathways with Spearman_FDR < 0.10",
    value = ifelse(
      nrow(sig_hallmark) > 0,
      paste(sig_hallmark$pathway, collapse = "; "),
      "None"
    ),
    stringsAsFactors = FALSE
  )
}

top_custom <- custom_assoc[order(-abs(custom_assoc$Spearman_r)), ]

interpretation_rows[[length(interpretation_rows) + 1]] <- data.frame(
  item = "Custom modules most associated with CRS_z",
  value = paste(
    paste0(
      top_custom$pathway,
      " (rho=",
      round(top_custom$Spearman_r, 3),
      ", FDR=",
      signif(top_custom$Spearman_FDR, 3),
      ")"
    ),
    collapse = "; "
  ),
  stringsAsFactors = FALSE
)

interpretation_rows[[length(interpretation_rows) + 1]] <- data.frame(
  item = "Important caution",
  value = "These pathway scores are exploratory mean-z signature scores and should be used for biological interpretation rather than causal inference.",
  stringsAsFactors = FALSE
)

interpretation_df <- do.call(rbind, interpretation_rows)

cat("\n===== Molecular-module analysis summary =====\n")
print(interpretation_df)

############################################################
# 12. Save outputs
############################################################

write.csv(
  all_pathway_scores,
  file = file.path(out_dir, "molecular_module_all_pathway_scores_by_patient.csv"),
  row.names = FALSE
)

write.csv(
  custom_gene_counts,
  file = file.path(out_dir, "molecular_module_custom_module_gene_counts.csv"),
  row.names = FALSE
)

write.csv(
  custom_assoc,
  file = file.path(out_dir, "predefined_module_crs_associations.csv"),
  row.names = FALSE
)

write.csv(
  interpretation_df,
  file = file.path(out_dir, "molecular_module_analysis_summary.csv"),
  row.names = FALSE
)

if (hallmark_available) {
  write.csv(
    hallmark_gene_counts,
    file = file.path(out_dir, "molecular_module_hallmark_gene_counts.csv"),
    row.names = FALSE
  )
  
  write.csv(
    hallmark_assoc,
    file = file.path(out_dir, "molecular_module_hallmark_CRS_association.csv"),
    row.names = FALSE
  )
}

############################################################
# 13. Excel workbook
############################################################

wb <- createWorkbook()

addWorksheet(wb, "analysis_summary")
writeData(wb, "analysis_summary", interpretation_df)

addWorksheet(wb, "custom_assoc")
writeData(wb, "custom_assoc", custom_assoc)

addWorksheet(wb, "custom_gene_counts")
writeData(wb, "custom_gene_counts", custom_gene_counts)

if (hallmark_available) {
  addWorksheet(wb, "hallmark_assoc")
  writeData(wb, "hallmark_assoc", hallmark_assoc)
  
  addWorksheet(wb, "hallmark_gene_counts")
  writeData(wb, "hallmark_gene_counts", hallmark_gene_counts)
}

addWorksheet(wb, "pathway_scores")
writeData(wb, "pathway_scores", all_pathway_scores)

for (sheet in names(wb)) {
  setColWidths(wb, sheet, cols = 1:20, widths = "auto")
}

saveWorkbook(
  wb,
  file = file.path(out_dir, "molecular_module_CRS_Molecular_Pathway_Analysis.xlsx"),
  overwrite = TRUE
)

############################################################
# 14. Figures
############################################################

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
  
  p1 <- ggplot(
    hallmark_plot_df,
    aes(x = pathway_clean, y = Spearman_r)
  ) +
    geom_col() +
    coord_flip() +
    theme_bw(base_size = 11) +
    labs(
      title = "Top CRS-associated Hallmark pathways",
      x = "",
      y = "Spearman correlation with CRS_z"
    )
  
  ggsave(
    filename = file.path(fig_dir, "molecular_module_Top_CRS_Associated_Hallmark_Pathways.png"),
    plot = p1,
    width = 8,
    height = 6,
    dpi = 300
  )
}

custom_plot_df <- custom_assoc
custom_plot_df$pathway <- factor(
  custom_plot_df$pathway,
  levels = custom_plot_df$pathway[order(custom_plot_df$Spearman_r)]
)

p2 <- ggplot(
  custom_plot_df,
  aes(x = pathway, y = Spearman_r)
) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(
    title = "Custom radiosensitivity-related modules associated with CRS",
    x = "",
    y = "Spearman correlation with CRS_z"
  )

ggsave(
  filename = file.path(fig_dir, "molecular_module_Custom_Module_CRS_Association.png"),
  plot = p2,
  width = 8,
  height = 5,
  dpi = 300
)

############################################################
# 15. Done
############################################################

cat("\n===== Molecular-module analysis completed =====\n")
cat("Main output files:\n")
cat(file.path(out_dir, "molecular_module_CRS_Molecular_Pathway_Analysis.xlsx"), "\n")
cat(file.path(out_dir, "molecular_module_analysis_summary.csv"), "\n")
cat(file.path(fig_dir, "molecular_module_Top_CRS_Associated_Hallmark_Pathways.png"), "\n")
cat(file.path(fig_dir, "molecular_module_Custom_Module_CRS_Association.png"), "\n")

cat("1. Hallmark association with CRS_z: top positive Spearman\n")
cat("2. Hallmark association with CRS_z: top negative Spearman\n")
cat("3. Hallmark significant pathways by Spearman_FDR < 0.10\n")
cat("4. Custom module association with CRS_z\n")
cat("5. Molecular-module analysis summary\n")
