# 02_score_hallmark_pathways.R
#
# Analyze Hallmark gene-set correlations with CRS.
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

msigdb_dir <- file.path(
  project_dir,
  "01_raw_data",
  "MSigDB"
)

out_dir <- file.path(
  project_dir,
  "07_results",
  "hallmark_pathways"
)

fig_dir <- file.path(
  project_dir,
  "08_figures"
)

dir.create(msigdb_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

gmt_file <- file.path(
  msigdb_dir,
  "h.all.v2026.1.Hs.symbols.gmt"
)

gmt_url <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/h.all.v2026.1.Hs.symbols.gmt"

############################################################
# 3. Download GMT if not exists
############################################################

download_gmt_safely <- function(url, destfile) {
  
  methods_to_try <- c("wininet", "libcurl", "auto")
  
  for (mm in methods_to_try) {
    
    cat("\nTrying download method:", mm, "\n")
    
    ok <- tryCatch(
      {
        if (mm == "auto") {
          download.file(
            url = url,
            destfile = destfile,
            mode = "wb",
            quiet = FALSE
          )
        } else {
          download.file(
            url = url,
            destfile = destfile,
            mode = "wb",
            method = mm,
            quiet = FALSE
          )
        }
        TRUE
      },
      error = function(e) {
        cat("Download failed with method", mm, ":\n")
        cat(conditionMessage(e), "\n")
        FALSE
      }
    )
    
    if (ok && file.exists(destfile) && file.info(destfile)$size > 100000) {
      cat("\nDownload succeeded with method:", mm, "\n")
      return(TRUE)
    }
  }
  
  FALSE
}

cat("\n===== Checking Hallmark GMT file =====\n")
cat("GMT path:\n")
cat(gmt_file, "\n")

if (!file.exists(gmt_file) || file.info(gmt_file)$size < 30000) {
  
  cat("\nGMT file not found or too small. Start downloading official MSigDB Hallmark GMT...\n")
  
  success <- download_gmt_safely(
    url = gmt_url,
    destfile = gmt_file
  )
  
  if (!success) {
    stop(
      "\nGMT download failed.\n",
      "Please manually download h.all.v2024.1.Hs.symbols.gmt from the official MSigDB release directory,\n",
      "then place it here:\n",
      gmt_file,
      "\nThen rerun this script.\n"
    )
  }
  
} else {
  cat("\nGMT file already exists. Skip download.\n")
}

cat("\nFinal GMT file size:\n")
print(file.info(gmt_file)$size)

############################################################
# 4. Read GMT file
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

cat("\n===== Reading Hallmark GMT =====\n")

hallmark_gene_sets <- read_gmt(gmt_file)

cat("Number of Hallmark gene sets:\n")
print(length(hallmark_gene_sets))

cat("\nFirst 10 Hallmark gene sets:\n")
print(head(names(hallmark_gene_sets), 10))

cat("\nGene count summary per Hallmark set:\n")
print(summary(sapply(hallmark_gene_sets, length)))

############################################################
# 5. Load expression and CRS scores
############################################################

cat("\n===== Loading processed expression =====\n")

expr <- readRDS(expr_rds)
expr <- as.data.frame(expr, check.names = FALSE)

rownames(expr) <- toupper(trimws(as.character(rownames(expr))))
expr <- expr[rownames(expr) != "" & !is.na(rownames(expr)), , drop = FALSE]

for (cc in colnames(expr)) {
  expr[[cc]] <- suppressWarnings(as.numeric(expr[[cc]]))
}

if (any(duplicated(rownames(expr)))) {
  
  cat("Duplicated genes detected. Aggregating by mean...\n")
  
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

cat("Expression dim after aggregation:\n")
print(dim(expr_mat))

cat("\n===== Loading CRS scores =====\n")

if (file.exists(score_rds)) {
  scores <- readRDS(score_rds)
} else {
  scores <- read.csv(score_csv, stringsAsFactors = FALSE, check.names = FALSE)
}

scores <- as.data.frame(scores, check.names = FALSE)

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
  stop("Cannot detect sample ID column.")
}

scores$patient_id <- toupper(trimws(as.character(scores[[id_col]])))
scores$CRS <- suppressWarnings(as.numeric(scores$CRS))
scores$CRS_z <- suppressWarnings(as.numeric(scores$CRS_z))

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
# 6. Pathway scoring
############################################################

score_gene_sets_mean_z <- function(expr_gene_by_sample, gene_sets, min_genes = 10) {
  
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

cat("\n===== Scoring Hallmark pathways =====\n")

hallmark_scored <- score_gene_sets_mean_z(
  expr_gene_by_sample = expr_mat,
  gene_sets = hallmark_gene_sets,
  min_genes = 10
)

hallmark_scores <- hallmark_scored$scores
hallmark_gene_counts <- hallmark_scored$gene_counts

hallmark_pathway_cols <- setdiff(colnames(hallmark_scores), "patient_id")

cat("\nHallmark pathway scores dim:\n")
print(dim(hallmark_scores))

cat("\nHallmark gene count summary:\n")
print(summary(hallmark_gene_counts$present_genes))

cat("\nNumber of scored Hallmark pathways:\n")
print(length(hallmark_pathway_cols))

############################################################
# 7. Association with CRS_z
############################################################

cat("\n===== Associating Hallmark pathways with CRS_z =====\n")

hallmark_assoc <- associate_pathways_with_CRS(
  pathway_score_df = hallmark_scores,
  scores_df = scores2,
  pathway_cols = hallmark_pathway_cols
)

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

############################################################
# 8. Analysis summary
############################################################

top_pos <- head(hallmark_assoc[order(-hallmark_assoc$Spearman_r), ], 10)
top_neg <- head(hallmark_assoc[order(hallmark_assoc$Spearman_r), ], 10)

sig_fdr_005 <- hallmark_assoc[hallmark_assoc$Spearman_FDR < 0.05, ]
sig_fdr_010 <- hallmark_assoc[hallmark_assoc$Spearman_FDR < 0.10, ]

interpretation_df <- data.frame(
  item = c(
    "Hallmark GMT version",
    "Number of Hallmark pathways scored",
    "Number of Hallmark pathways with Spearman_FDR < 0.05",
    "Number of Hallmark pathways with Spearman_FDR < 0.10",
    "Top positively CRS-associated Hallmark pathways",
    "Top negatively CRS-associated Hallmark pathways",
    "Important caution"
  ),
  value = c(
    "MSigDB Hallmark h.all.v2026.1.Hs.symbols.gmt",
    length(hallmark_pathway_cols),
    nrow(sig_fdr_005),
    nrow(sig_fdr_010),
    paste(
      paste0(
        top_pos$pathway,
        " (rho=",
        round(top_pos$Spearman_r, 3),
        ", FDR=",
        signif(top_pos$Spearman_FDR, 3),
        ")"
      ),
      collapse = "; "
    ),
    paste(
      paste0(
        top_neg$pathway,
        " (rho=",
        round(top_neg$Spearman_r, 3),
        ", FDR=",
        signif(top_neg$Spearman_FDR, 3),
        ")"
      ),
      collapse = "; "
    ),
    "Hallmark pathway scores were calculated as mean z-scores of expressed member genes and used for exploratory biological interpretation."
  ),
  stringsAsFactors = FALSE
)

cat("\n===== Hallmark pathway analysis Hallmark GMT analysis summary =====\n")
print(interpretation_df)

############################################################
# 9. Save outputs
############################################################

write.csv(
  hallmark_scores,
  file = file.path(out_dir, "hallmark_hallmark_pathway_scores_by_patient.csv"),
  row.names = FALSE
)

write.csv(
  hallmark_gene_counts,
  file = file.path(out_dir, "hallmark_hallmark_gene_counts.csv"),
  row.names = FALSE
)

write.csv(
  hallmark_assoc,
  file = file.path(out_dir, "hallmark_crs_associations.csv"),
  row.names = FALSE
)

write.csv(
  interpretation_df,
  file = file.path(out_dir, "hallmark_analysis_summary.csv"),
  row.names = FALSE
)

wb <- createWorkbook()

addWorksheet(wb, "analysis_summary")
writeData(wb, "analysis_summary", interpretation_df)

addWorksheet(wb, "hallmark_assoc")
writeData(wb, "hallmark_assoc", hallmark_assoc)

addWorksheet(wb, "hallmark_gene_counts")
writeData(wb, "hallmark_gene_counts", hallmark_gene_counts)

addWorksheet(wb, "hallmark_scores")
writeData(wb, "hallmark_scores", hallmark_scores)

for (sheet in names(wb)) {
  setColWidths(wb, sheet, cols = 1:30, widths = "auto")
}

saveWorkbook(
  wb,
  file = file.path(out_dir, "hallmark_Hallmark_GMT_CRS_Pathway_Analysis.xlsx"),
  overwrite = TRUE
)

############################################################
# 10. Figures
############################################################

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
  filename = file.path(fig_dir, "hallmark_Top_CRS_Associated_Hallmark_Pathways.png"),
  plot = p1,
  width = 8,
  height = 6,
  dpi = 300
)

sig_plot_df <- hallmark_assoc[hallmark_assoc$Spearman_FDR < 0.10, ]

if (nrow(sig_plot_df) > 0) {
  
  sig_plot_df$pathway_clean <- gsub("^HALLMARK_", "", sig_plot_df$pathway)
  sig_plot_df$pathway_clean <- gsub("_", " ", sig_plot_df$pathway_clean)
  
  sig_plot_df$pathway_clean <- factor(
    sig_plot_df$pathway_clean,
    levels = sig_plot_df$pathway_clean[order(sig_plot_df$Spearman_r)]
  )
  
  p2 <- ggplot(
    sig_plot_df,
    aes(x = pathway_clean, y = Spearman_r)
  ) +
    geom_col() +
    coord_flip() +
    theme_bw(base_size = 11) +
    labs(
      title = "Hallmark pathways associated with CRS_z, FDR < 0.10",
      x = "",
      y = "Spearman correlation with CRS_z"
    )
  
  ggsave(
    filename = file.path(fig_dir, "hallmark_Significant_Hallmark_Pathways_FDR010.png"),
    plot = p2,
    width = 8,
    height = 6,
    dpi = 300
  )
}

############################################################
# 11. Key results text
############################################################

txt_lines <- c(
  "Hallmark pathway analysis summary",
  "",
  "Database:",
  "MSigDB Hallmark h.all.v2026.1.Hs.symbols.gmt",
  "",
  "Top positive Hallmark pathways:",
  paste(
    paste0(
      top_pos$pathway,
      ": rho=",
      round(top_pos$Spearman_r, 3),
      ", FDR=",
      signif(top_pos$Spearman_FDR, 3)
    ),
    collapse = "\n"
  ),
  "",
  "Top negative Hallmark pathways:",
  paste(
    paste0(
      top_neg$pathway,
      ": rho=",
      round(top_neg$Spearman_r, 3),
      ", FDR=",
      signif(top_neg$Spearman_FDR, 3)
    ),
    collapse = "\n"
  ),
  "",
  "Significant pathways by Spearman_FDR < 0.10:",
  ifelse(
    nrow(sig_fdr_010) > 0,
    paste(
      paste0(
        sig_fdr_010$pathway,
        ": rho=",
        round(sig_fdr_010$Spearman_r, 3),
        ", FDR=",
        signif(sig_fdr_010$Spearman_FDR, 3)
      ),
      collapse = "\n"
    ),
    "None"
  ),
  "",
  "Recommended cautious wording:",
  "Hallmark pathway analysis suggested that CRS_z was associated with multiple canonical cancer-related biological programs. These exploratory pathway-level associations support the biological relevance of CRS, but should not be interpreted as causal mechanisms."
)

writeLines(
  txt_lines,
  con = file.path(out_dir, "hallmark_key_results.txt"),
  useBytes = TRUE
)

############################################################
# 12. Done
############################################################

cat("\n===== DONE: Hallmark pathway analysis Hallmark GMT CRS pathway analysis finished =====\n")
cat("Main output files:\n")
cat(file.path(out_dir, "hallmark_Hallmark_GMT_CRS_Pathway_Analysis.xlsx"), "\n")
cat(file.path(out_dir, "hallmark_analysis_summary.csv"), "\n")
cat(file.path(out_dir, "hallmark_key_results.txt"), "\n")
cat(file.path(fig_dir, "hallmark_Top_CRS_Associated_Hallmark_Pathways.png"), "\n")

cat("1. Number of Hallmark gene sets\n")
cat("2. Hallmark association with CRS_z: top positive Spearman\n")
cat("3. Hallmark association with CRS_z: top negative Spearman\n")
cat("4. Hallmark significant pathways by Spearman_FDR < 0.10\n")
cat("5. Hallmark pathway analysis Hallmark GMT analysis summary\n")
