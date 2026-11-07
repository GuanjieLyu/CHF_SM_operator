#!/usr/bin/env Rscript
# Install dependencies only when this setup command is explicitly invoked.
# The dependency record is embedded so no separate lockfile is needed.

main <- function(args) {
  if (length(args) && !identical(args, "--help")) {
    stop("Unknown option. Use --help.", call. = FALSE)
  }
  if (identical(args, "--help")) {
    cat(paste(c(
      "Usage: Rscript scripts/setup.R",
      "Requires R 4.5.1 and Internet access.",
      "Restores muhaz 1.2.6.4, survival 3.8-3, Matrix 1.7-3,",
      "lattice 0.22-7, and KMsurv 0.1-6 to the project-local .library/.",
      "KMsurv and survival provide the public clinical example data.",
      "Source installation may require an R-compatible compiler toolchain.",
      "The temporary dependency lock is removed after setup.",
      "R and compiler toolchains are not installed by this script."
    ), collapse = "\n"), "\n")
    return(invisible(NULL))
  }
  if (getRversion() != package_version("4.5.1")) {
    stop("Use R 4.5.1 to restore this study's pinned environment.", call. = FALSE)
  }
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(script_arg) != 1L) {
    stop("Invoke this script with Rscript scripts/setup.R.", call. = FALSE)
  }
  script <- sub("^--file=", "", script_arg)
  root <- dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE)))
  source(file.path(root, "R", "bootstrap.R"), local = TRUE, encoding = "UTF-8")
  library <- project_library(root)
  dir.create(library, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(library)) {
    stop("Cannot create the project-local dependency library.", call. = FALSE)
  }
  .libPaths(c(library, .libPaths()))
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  if (!requireNamespace("renv", quietly = TRUE)) {
    utils::install.packages("renv", lib = library)
  }
  if (!requireNamespace("renv", quietly = TRUE)) {
    stop("Could not install renv. Check the setup output and Internet connection.",
      call. = FALSE
    )
  }

  versions <- c(
    muhaz = "1.2.6.4", survival = "3.8-3", Matrix = "1.7-3",
    lattice = "0.22-7", KMsurv = "0.1-6"
  )
  records <- sprintf(
    '    "%s": {"Package": "%s", "Version": "%s", "Source": "Repository", "Repository": "CRAN"}',
    names(versions), names(versions), unname(versions)
  )
  temporary_project <- tempfile("szasz-setup-", tmpdir = tempdir())
  if (!dir.create(temporary_project)) {
    stop("Cannot create temporary setup directory.", call. = FALSE)
  }
  temporary_project <- normalizePath(temporary_project, winslash = "/", mustWork = TRUE)
  temporary_root <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
  if (!startsWith(temporary_project, paste0(temporary_root, "/"))) {
    stop("Temporary setup path is outside the R session's temporary directory.",
      call. = FALSE
    )
  }
  on.exit(unlink(temporary_project, recursive = TRUE), add = TRUE)
  lockfile <- file.path(temporary_project, "renv.lock")
  writeLines(c(
    "{",
    '  "R": {"Version": "4.5.1", "Repositories": [{"Name": "CRAN", "URL": "https://cloud.r-project.org"}]},',
    '  "Packages": {',
    paste(records, collapse = ",\n"),
    "  }",
    "}"
  ), lockfile)

  # Keep renv project metadata in the temporary directory, not the source tree.
  renv::restore(
    project = temporary_project, library = library, lockfile = lockfile,
    prompt = FALSE, clean = FALSE, transactional = TRUE
  )
  installed <- vapply(names(versions), function(package) {
    # renv can reuse matching recommended packages from R's system library.
    as.character(utils::packageVersion(package))
  }, character(1))
  matches <- vapply(names(versions), function(package) {
    package_version(installed[[package]]) == package_version(versions[[package]])
  }, logical(1))
  if (!all(matches)) {
    stop("Setup did not restore the required versions of: ",
      paste(names(versions)[!matches], collapse = ", "),
      call. = FALSE
    )
  }
  cat("Pinned dependencies restored to:\n", library,
    "\nRun Rscript scripts/check.R next.\n",
    sep = ""
  )
  invisible(library)
}

tryCatch(main(commandArgs(trailingOnly = TRUE)), error = function(error) {
  message("ERROR: ", conditionMessage(error))
  quit(status = 1L)
})
