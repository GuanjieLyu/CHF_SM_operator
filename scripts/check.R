#!/usr/bin/env Rscript
# Run mathematical and design checks without archived result files.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE)))
source(file.path(root, "R", "bootstrap.R"), encoding = "UTF-8")
load_project(root)
args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--help")) {
  cat("Usage: Rscript scripts/check.R\nChecks mathematical identities, input validation, censoring calibration, and summaries.\n")
} else {
  if (length(args)) stop("Unknown option. Use --help.", call. = FALSE)
  setwd(root)
  tests <- sort(list.files(file.path(root, "tests"), pattern = "^test_.*\\.R$", full.names = TRUE))
  for (test in tests) {
    message("Checking ", basename(test))
    sys.source(test, envir = new.env(parent = globalenv()))
  }
  message("All checks passed (", length(tests), " test files).")
}
