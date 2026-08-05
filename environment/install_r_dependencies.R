cran_packages <- c(
  "curl",
  "dplyr",
  "future",
  "future.apply",
  "ggplot2",
  "glmnet",
  "jsonlite",
  "msigdbr",
  "openxlsx",
  "oro.dicom",
  "patchwork",
  "png",
  "pROC",
  "readxl",
  "scales",
  "survival",
  "tidyr"
)

bioconductor_packages <- c(
  "CoreGx",
  "RadioGx"
)

missing_cran <- cran_packages[!vapply(
  cran_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)]
if (length(missing_cran) > 0L) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

missing_bioconductor <- bioconductor_packages[!vapply(
  bioconductor_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)]
if (length(missing_bioconductor) > 0L) {
  BiocManager::install(
    missing_bioconductor,
    ask = FALSE,
    update = FALSE
  )
}
