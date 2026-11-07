#!/usr/bin/env Rscript
# Run from any working directory. See --help for commands and profiles.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE)))
source(file.path(root, "R", "bootstrap.R"), encoding = "UTF-8")
load_project(root)

usage <- function() {
  cat(paste(c(
    "Usage: Rscript scripts/reproduce.R [options]",
    "",
    "  --stage all|simulate|summarize|figures|applications",
    "  --profile paper|smoke       Default: paper",
    "  --output PATH              Default: outputs/<profile>",
    "  --inputs PATH,PATH         Completed shard directories for summarize",
    "  --shard-index N            One-based index for simulate only",
    "  --shard-count N            Number of disjoint shards for simulate only",
    "  --help                     Print this message",
    "",
    "paper: 500 replications in each of 20 cells, 2001 integration points.",
    "smoke: 3 replications per cell and reduced resolution, for testing only.",
    "Application data are loaded from KMsurv and survival.",
    "Each shard must use its own output directory. Do not share outputs",
    "between profiles or code versions. No manuscript files are modified."
  ), collapse = "\n"), "\n")
}

main <- function(args) {
  if ("--help" %in% args) {
    usage()
    return(invisible(NULL))
  }
  options <- list(
    stage = "all", profile = "paper", output = NULL,
    inputs = NULL, shard_index = 1L, shard_count = 1L
  )
  allowed <- c("--stage", "--profile", "--output", "--inputs", "--shard-index", "--shard-count")
  if (length(args) %% 2L != 0L) stop("Options require a value. See --help.", call. = FALSE)
  if (length(args)) {
    keys <- args[seq(1L, length(args), by = 2L)]
    if (anyDuplicated(keys)) stop("Each option may be supplied only once.", call. = FALSE)
    for (i in seq(1L, length(args), by = 2L)) {
      if (!args[i] %in% allowed) stop("Unknown option: ", args[i], call. = FALSE)
      key <- gsub("-", "_", sub("^--", "", args[i]), fixed = TRUE)
      options[[key]] <- args[i + 1L]
    }
  }
  stages <- c("all", "simulate", "summarize", "figures", "applications")
  if (!options$stage %in% stages) stop("Unknown stage. See --help.", call. = FALSE)
  config <- study_config(options$profile)
  dependency_versions()
  shard_index <- suppressWarnings(as.numeric(options$shard_index))
  shard_count <- suppressWarnings(as.numeric(options$shard_count))
  if (options$stage != "simulate" && (!identical(shard_index, 1) || !identical(shard_count, 1))) {
    stop("Sharding is supported only with --stage simulate.", call. = FALSE)
  }
  output <- options$output
  if (is.null(output)) {
    output <- file.path(root, "outputs", config$profile)
  }
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  output <- normalizePath(output, winslash = "/", mustWork = TRUE)
  protected <- c(
    root, file.path(root, "R"), file.path(root, "scripts"),
    file.path(root, "reference"), file.path(root, "data"), file.path(root, "tests"),
    file.path(root, "docs"), file.path(root, ".github"), file.path(root, ".git"),
    file.path(root, ".library")
  )
  if (any(vapply(protected, function(path) {
    path <- normalizePath(path, winslash = "/", mustWork = FALSE)
    identical(tolower(output), tolower(path)) || startsWith(tolower(output), paste0(tolower(path), "/")) && path != root
  }, logical(1)))) {
    stop("Choose an output directory separate from source, input data, and reference results.", call. = FALSE)
  }
  simulation_dir <- file.path(output, "simulation")
  mode <- "live"
  prepare_output(root, output, config, mode, create = FALSE)
  if (config$profile == "smoke") message("SMOKE TEST: reduced calculations are not manuscript results.")
  with_output_lock(output, function() {
    prepare_output(root, output, config, mode)
    if (options$stage %in% c("all", "simulate")) {
      run_simulations(root, simulation_dir, config, shard_index, shard_count)
    }
    if (options$stage %in% c("all", "summarize")) {
      inputs <- if (is.null(options$inputs)) simulation_dir else strsplit(options$inputs, ",", fixed = TRUE)[[1]]
      aggregate_results(root, simulation_dir, config, inputs = inputs)
    }
    if (options$stage %in% c("all", "figures", "applications")) {
      with_staged_artifacts(output, function(staging) {
        if (options$stage %in% c("all", "figures")) {
          write_simulation_figures(simulation_dir, file.path(staging, "figures"),
            file.path(staging, "tables"),
            n_boot = config$n_boot
          )
        }
        if (options$stage %in% c("all", "applications")) {
          write_application_figures(NULL, file.path(staging, "figures"),
            file.path(staging, "applications"),
            illustration_reps = config$illustration_reps,
            kernel_grid = config$application_grid
          )
        }
      })
    }
    atomic_write(file.path(output, "sessionInfo.txt"), function(destination) {
      writeLines(capture.output(sessionInfo()), destination)
    })
    atomic_write(file.path(output, "run.txt"), function(destination) {
      writeLines(c(
        paste("profile:", config$profile), paste("stage:", options$stage),
        paste("mode:", mode),
        paste("completed UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE))
      ), destination)
    })
  }, lock_name = ".reproduce-lock")
  message("Completed. Outputs: ", output)
  invisible(output)
}

tryCatch(main(commandArgs(trailingOnly = TRUE)), error = function(error) {
  message("ERROR: ", conditionMessage(error))
  quit(status = 1L)
})
