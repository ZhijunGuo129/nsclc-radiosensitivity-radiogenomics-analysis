# 06_analyze_stage_adjusted_survival.R
#
# Evaluate the complete-case Overall.Stage-adjusted RRS association in Lung1.

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


clinical_file <- file.path(
  lung1_root,
  "clinical",
  "Lung1_clinical_clean_initial.rds"
)
rrs_candidates <- c(
  file.path(project_dir, "07_results", "tumor_size_domain_shift", "lung1_size_analysis_dataset.csv"),
  file.path(lung1_root, "tumor_size_domain_shift", "lung1_size_analysis_dataset.csv")
)
rrs_file <- rrs_candidates[file.exists(rrs_candidates)][1]
out_dir <- file.path(project_dir, "07_results", "lung1_stage_adjusted_survival")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(clinical_file)) {
  stop("Lung1 clinical data were not found: ", clinical_file)
}
if (length(rrs_file) == 0L || is.na(rrs_file)) {
  stop("Lung1 RRS data were not found. Checked: ", paste(rrs_candidates, collapse = "; "))
}

clean_id <- function(x) toupper(trimws(as.character(x)))
clean_stage <- function(x) {
  value <- trimws(as.character(x))
  value[value %in% c("", "NA", "N/A", "Unknown", "Not Available", "Not Reported", "[Not Available]", "[Not Reported]")] <- NA_character_
  value
}
find_column <- function(names_vector, candidates) {
  index <- match(tolower(candidates), tolower(names_vector))
  index <- index[!is.na(index)]
  if (length(index) == 0L) NA_character_ else names_vector[index[1]]
}

clinical <- readRDS(clinical_file)
rrs <- read.csv(rrs_file, stringsAsFactors = FALSE, check.names = FALSE)

if (!"patient_id" %in% names(clinical)) {
  id_column <- find_column(names(clinical), c("PatientID", "Patient.ID", "Patient_Id", "Case.ID"))
  if (is.na(id_column)) stop("Patient identifier column was not found in the clinical data.")
  clinical$patient_id <- clinical[[id_column]]
}
if (!all(c("patient_id", "Original_Lung1_RRS") %in% names(rrs))) {
  stop("The RRS data must contain patient_id and Original_Lung1_RRS.")
}

clinical$patient_id <- clean_id(clinical$patient_id)
rrs$patient_id <- clean_id(rrs$patient_id)
data <- merge(
  clinical,
  rrs[, c("patient_id", "Original_Lung1_RRS")],
  by = "patient_id",
  all = FALSE
)

time_column <- find_column(
  names(data),
  c("Survival.time", "Survival_time", "OS.time", "OS_time", "overall_survival_time")
)
event_column <- find_column(
  names(data),
  c("deadstatus.event", "deadstatus_event", "OS.event", "OS_event", "overall_survival_event")
)
stage_column <- find_column(
  names(data),
  c("Overall.Stage", "Overall_Stage", "Overall Stage", "Stage.Group", "Stage_group", "stage_group", "overall_stage")
)

if (any(is.na(c(time_column, event_column, stage_column)))) {
  stop("Required survival or Overall.Stage column was not found.")
}

data$OS_time_months <- suppressWarnings(as.numeric(data[[time_column]])) / 30.4375
data$OS_event <- suppressWarnings(as.numeric(data[[event_column]]))
data$Original_Lung1_RRS <- suppressWarnings(as.numeric(data$Original_Lung1_RRS))
data$Overall_Stage <- clean_stage(data[[stage_column]])

data <- data[
  is.finite(data$OS_time_months) & data$OS_time_months > 0 &
    data$OS_event %in% c(0, 1) &
    is.finite(data$Original_Lung1_RRS),
  ,
  drop = FALSE
]

rrs_cutoff <- median(data$Original_Lung1_RRS, na.rm = TRUE)
data$RRS_group <- factor(
  ifelse(data$Original_Lung1_RRS > rrs_cutoff, "RRS_high", "RRS_low"),
  levels = c("RRS_low", "RRS_high")
)
data$Overall_Stage <- factor(
  data$Overall_Stage,
  levels = c("I", "II", "IIIa", "IIIb")
)

analysis_data <- data[
  complete.cases(data[, c("OS_time_months", "OS_event", "RRS_group", "Overall_Stage")]),
  ,
  drop = FALSE
]

fit <- survival::coxph(
  survival::Surv(OS_time_months, OS_event) ~ RRS_group + Overall_Stage,
  data = analysis_data,
  ties = "efron",
  x = TRUE,
  model = TRUE
)

fit_summary <- summary(fit)
term <- grep("^RRS_group", rownames(fit_summary$coefficients), value = TRUE)
if (length(term) != 1L) stop("The RRS coefficient could not be identified in the Cox model.")
interval <- stats::confint(fit, parm = term)

result <- data.frame(
  model = "Overall.Stage-adjusted complete-case Cox model",
  stage_source_column = stage_column,
  n = nrow(analysis_data),
  events = sum(analysis_data$OS_event == 1),
  HR = unname(exp(fit_summary$coefficients[term, "coef"])),
  CI_lower = unname(exp(interval[1])),
  CI_upper = unname(exp(interval[2])),
  p_value = unname(fit_summary$coefficients[term, "Pr(>|z|)"]),
  stringsAsFactors = FALSE
)

stage_counts <- as.data.frame(table(analysis_data$Overall_Stage, useNA = "ifany"))
names(stage_counts) <- c("Overall_Stage", "n")

write.csv(result, file.path(out_dir, "stage_adjusted_rrs_cox_result.csv"), row.names = FALSE)
write.csv(stage_counts, file.path(out_dir, "overall_stage_counts.csv"), row.names = FALSE)
write.csv(analysis_data, file.path(out_dir, "stage_adjusted_analysis_data.csv"), row.names = FALSE)

workbook <- openxlsx::createWorkbook()
openxlsx::addWorksheet(workbook, "result")
openxlsx::writeData(workbook, "result", result)
openxlsx::addWorksheet(workbook, "stage_counts")
openxlsx::writeData(workbook, "stage_counts", stage_counts)
openxlsx::saveWorkbook(
  workbook,
  file.path(out_dir, "stage_adjusted_rrs_cox_summary.xlsx"),
  overwrite = TRUE
)

cat(
  sprintf(
    "Overall.Stage-adjusted RRS result: n=%d, events=%d, HR=%.6f, 95%% CI %.6f-%.6f, P=%.6f\n",
    result$n,
    result$events,
    result$HR,
    result$CI_lower,
    result$CI_upper,
    result$p_value
  )
)
