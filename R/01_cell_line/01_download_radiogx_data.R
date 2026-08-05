# 01_download_radiogx_data.R
#
# Download the Cleveland RadioSet and export source-data audit files.

options(stringsAsFactors = FALSE)
options(timeout = 1800)

required_env_path <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    stop("Environment variable ", name, " is not set.")
  }
  normalizePath(value, winslash = "/", mustWork = TRUE)
}

project_dir <- required_env_path("RADIOGENOMICS_PROJECT_DIR")
setwd(project_dir)

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
invisible(lapply(required_packages, library, character.only = TRUE))


raw_dir <- file.path(project_dir, "01_raw_data", "RadioGx")
metadata_dir <- file.path(project_dir, "02_metadata")
log_dir <- file.path(project_dir, "10_logs")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

rset_file <- file.path(raw_dir, "Cleveland_or_clevelandSmall_RadioSet.rds")

if (file.exists(rset_file)) {
  cleveland <- readRDS(rset_file)
} else {
  cleveland <- tryCatch(
    RadioGx::downloadRSet(name = "Cleveland", saveDir = raw_dir),
    error = function(e) NULL
  )

  if (is.null(cleveland)) {
    data("clevelandSmall", package = "RadioGx", envir = environment())
    cleveland <- get("clevelandSmall", envir = environment())
  }

  saveRDS(cleveland, rset_file)
}

cell_info <- as.data.frame(CoreGx::cellInfo(cleveland))
sensitivity_info <- as.data.frame(CoreGx::sensitivityInfo(cleveland))
sensitivity_profiles <- as.data.frame(CoreGx::sensitivityProfiles(cleveland))
molecular_data_names <- CoreGx::mDataNames(cleveland)

openxlsx::write.xlsx(
  cell_info,
  file.path(metadata_dir, "RadioGx_cell_info_raw.xlsx"),
  overwrite = TRUE
)
openxlsx::write.xlsx(
  sensitivity_info,
  file.path(metadata_dir, "RadioGx_sensitivity_info_raw.xlsx"),
  overwrite = TRUE
)
openxlsx::write.xlsx(
  sensitivity_profiles,
  file.path(metadata_dir, "RadioGx_sensitivity_profiles_raw.xlsx"),
  overwrite = TRUE
)
write.csv(
  data.frame(molecular_profile = molecular_data_names),
  file.path(metadata_dir, "RadioGx_molecular_profile_names.csv"),
  row.names = FALSE
)

profile_dimensions <- lapply(molecular_data_names, function(profile_name) {
  profile <- tryCatch(
    CoreGx::molecularProfiles(cleveland, profile_name),
    error = function(e) NULL
  )
  dims <- if (is.null(profile)) c(NA_integer_, NA_integer_) else dim(profile)
  data.frame(
    molecular_profile = profile_name,
    n_features = dims[1],
    n_samples = dims[2],
    stringsAsFactors = FALSE
  )
})
profile_dimensions <- do.call(rbind, profile_dimensions)

summary_table <- data.frame(
  item = c(
    "RadioSet class",
    "Cell lines",
    "Sensitivity records",
    "Sensitivity profiles",
    "Molecular profiles"
  ),
  value = c(
    paste(class(cleveland), collapse = "; "),
    nrow(cell_info),
    nrow(sensitivity_info),
    nrow(sensitivity_profiles),
    length(molecular_data_names)
  ),
  stringsAsFactors = FALSE
)

workbook <- openxlsx::createWorkbook()
openxlsx::addWorksheet(workbook, "summary")
openxlsx::writeData(workbook, "summary", summary_table)
openxlsx::addWorksheet(workbook, "profile_dimensions")
openxlsx::writeData(workbook, "profile_dimensions", profile_dimensions)
openxlsx::saveWorkbook(
  workbook,
  file.path(metadata_dir, "radiogx_source_data_summary.xlsx"),
  overwrite = TRUE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(log_dir, "RadioGx_source_data_session_info.txt")
)

cat("RadioGx source data audit completed.
")
cat("RadioSet file:", rset_file, "
")
