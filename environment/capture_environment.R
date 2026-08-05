output_dir <- file.path("environment", "captured")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "r_session_info.txt"))
