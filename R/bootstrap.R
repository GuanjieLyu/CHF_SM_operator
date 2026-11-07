# Portable loading for command-line scripts and interactive examples.
# Sourcing this file only defines functions and does not start an analysis.

project_library <- function(root) {
  version <- paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = ".")
  file.path(root, ".library", paste0("R-", version), R.version$platform)
}

load_project <- function(root, envir = parent.frame(), use_local_library = TRUE) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  library <- project_library(root)
  if (use_local_library && dir.exists(library)) {
    .libPaths(c(library, .libPaths()))
  }
  modules <- c(
    "validation.R", "estimators.R", "design.R", "simulation.R",
    "pipeline.R", "plot_simulations.R", "plot_applications.R"
  )
  for (module in modules) {
    source(file.path(root, "R", module), local = envir, encoding = "UTF-8")
  }
  invisible(root)
}
