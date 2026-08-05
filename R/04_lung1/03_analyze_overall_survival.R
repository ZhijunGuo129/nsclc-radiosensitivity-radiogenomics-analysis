# 03_analyze_overall_survival.R
#
# Evaluate exploratory overall-survival associations for the primary RRS in Lung1.

options(stringsAsFactors = FALSE)

required_env_path <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    stop("Environment variable ", name, " is not set.")
  }
  normalizePath(value, winslash = "/", mustWork = TRUE)
}

project_dir <- required_env_path("RADIOGENOMICS_PROJECT_DIR")
lung1_root <- required_env_path("LUNG1_DATA_DIR")
setwd(project_dir)

required_packages <- c("survival", "openxlsx")
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


input_file <- file.path(
  lung1_root,
  "external_RRS",
  "lung1_rrs_survival_dataset.rds"
)
out_dir <- file.path(lung1_root, "survival_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Lung1 RRS-survival analysis table was not found: ", input_file)
}

data <- readRDS(input_file)
required_columns <- c("patient_id", "OS_time_months", "OS_event", "Lung1_RRS")
missing_columns <- setdiff(required_columns, names(data))
if (length(missing_columns) > 0L) {
  stop("Missing column(s): ", paste(missing_columns, collapse = ", "))
}

data$patient_id <- toupper(trimws(as.character(data$patient_id)))
data$OS_time_months <- suppressWarnings(as.numeric(data$OS_time_months))
data$OS_event <- suppressWarnings(as.numeric(data$OS_event))
data$Lung1_RRS <- suppressWarnings(as.numeric(data$Lung1_RRS))

data <- data[
  is.finite(data$OS_time_months) &
    data$OS_time_months > 0 &
    data$OS_event %in% c(0, 1) &
    is.finite(data$Lung1_RRS),
  ,
  drop = FALSE
]

rrs_cutoff <- median(data$Lung1_RRS, na.rm = TRUE)
data$RRS_group <- factor(
  ifelse(data$Lung1_RRS >= rrs_cutoff, "RRS_high", "RRS_low"),
  levels = c("RRS_low", "RRS_high")
)
data$RRS_z <- as.numeric(scale(data$Lung1_RRS))
data$RRS_rank_z <- as.numeric(scale(rank(data$Lung1_RRS, ties.method = "average")))

extract_cox <- function(fit, model_name) {
  fit_summary <- summary(fit)
  coefficients <- as.data.frame(fit_summary$coefficients)
  intervals <- as.data.frame(fit_summary$conf.int)
  data.frame(
    model = model_name,
    term = rownames(coefficients),
    n = fit$n,
    events = fit$nevent,
    HR = intervals$`exp(coef)`,
    CI_lower = intervals$`lower .95`,
    CI_upper = intervals$`upper .95`,
    p_value = coefficients$`Pr(>|z|)`,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

binary_fit <- survival::coxph(
  survival::Surv(OS_time_months, OS_event) ~ RRS_group,
  data = data,
  ties = "efron"
)
continuous_fit <- survival::coxph(
  survival::Surv(OS_time_months, OS_event) ~ RRS_z,
  data = data,
  ties = "efron"
)
rank_fit <- survival::coxph(
  survival::Surv(OS_time_months, OS_event) ~ RRS_rank_z,
  data = data,
  ties = "efron"
)

cox_results <- rbind(
  extract_cox(binary_fit, "RRS high versus low"),
  extract_cox(continuous_fit, "RRS per standard deviation"),
  extract_cox(rank_fit, "Rank-standardized RRS")
)

km_fit <- survival::survfit(
  survival::Surv(OS_time_months, OS_event) ~ RRS_group,
  data = data
)
logrank <- survival::survdiff(
  survival::Surv(OS_time_months, OS_event) ~ RRS_group,
  data = data
)
logrank_p <- stats::pchisq(logrank$chisq, df = length(logrank$n) - 1L, lower.tail = FALSE)

q1 <- unname(stats::quantile(data$Lung1_RRS, 0.25, na.rm = TRUE))
q3 <- unname(stats::quantile(data$Lung1_RRS, 0.75, na.rm = TRUE))
iqr <- q3 - q1
data$extreme_RRS <- data$Lung1_RRS < q1 - 3 * iqr | data$Lung1_RRS > q3 + 3 * iqr | abs(data$RRS_z) > 5

sensitivity_data <- data[!data$extreme_RRS, , drop = FALSE]
sensitivity_fit <- survival::coxph(
  survival::Surv(OS_time_months, OS_event) ~ RRS_group,
  data = sensitivity_data,
  ties = "efron"
)
sensitivity_results <- extract_cox(
  sensitivity_fit,
  "RRS high versus low after excluding extreme scores"
)

summary_table <- data.frame(
  item = c(
    "Patients",
    "OS events",
    "Median OS, months",
    "Median RRS cutoff",
    "RRS-high patients",
    "RRS-low patients",
    "Log-rank P",
    "Extreme-score patients"
  ),
  value = c(
    nrow(data),
    sum(data$OS_event == 1),
    median(data$OS_time_months, na.rm = TRUE),
    rrs_cutoff,
    sum(data$RRS_group == "RRS_high"),
    sum(data$RRS_group == "RRS_low"),
    logrank_p,
    sum(data$extreme_RRS)
  ),
  stringsAsFactors = FALSE
)

write.csv(data, file.path(out_dir, "primary_rrs_os_analysis_data.csv"), row.names = FALSE)
write.csv(cox_results, file.path(out_dir, "primary_rrs_cox_results.csv"), row.names = FALSE)
write.csv(sensitivity_results, file.path(out_dir, "primary_rrs_outlier_sensitivity_results.csv"), row.names = FALSE)
write.csv(summary_table, file.path(out_dir, "primary_rrs_os_summary.csv"), row.names = FALSE)

workbook <- openxlsx::createWorkbook()
openxlsx::addWorksheet(workbook, "summary")
openxlsx::writeData(workbook, "summary", summary_table)
openxlsx::addWorksheet(workbook, "cox_results")
openxlsx::writeData(workbook, "cox_results", cox_results)
openxlsx::addWorksheet(workbook, "outlier_sensitivity")
openxlsx::writeData(workbook, "outlier_sensitivity", sensitivity_results)
openxlsx::saveWorkbook(
  workbook,
  file.path(out_dir, "primary_rrs_os_analysis.xlsx"),
  overwrite = TRUE
)

binary_row <- cox_results[grepl("RRS_group", cox_results$term), , drop = FALSE]
cat("Lung1 primary RRS OS analysis completed.\n")
cat("Patients:", nrow(data), "Events:", sum(data$OS_event == 1), "\n")
if (nrow(binary_row) == 1L) {
  cat(
    sprintf(
      "RRS high versus low: HR %.6f, 95%% CI %.6f-%.6f, P %.6f\n",
      binary_row$HR,
      binary_row$CI_lower,
      binary_row$CI_upper,
      binary_row$p_value
    )
  )
}
