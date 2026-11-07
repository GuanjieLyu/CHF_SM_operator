# Orchestration and provenance for reproducible experiments.
# The numerical functions themselves are defined in the other R modules.

study_config <- function(profile = c("paper", "smoke")) {
  profile <- match.arg(profile)
  list(
    schema_version = 1L,
    profile = profile,
    replications = if (profile == "paper") 500L else 3L,
    n_grid = if (profile == "paper") 2001L else 201L,
    n_boot = if (profile == "paper") 1000L else 50L,
    m_grid = c(1, seq(2, 100, 2)),
    tau = 3.7,
    illustration_reps = if (profile == "paper") 50L else 3L,
    application_grid = if (profile == "paper") 4001L else 401L
  )
}

expected_dependency_versions <- function() {
  c(muhaz = "1.2.6.4", survival = "3.8-3", Matrix = "1.7-3", lattice = "0.22-7", KMsurv = "0.1-6")
}

dependency_versions <- function() {
  expected <- expected_dependency_versions()
  packages <- names(expected)
  versions <- vapply(packages, function(package) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Missing dependency '", package, "'. Run Rscript scripts/setup.R.", call. = FALSE)
    }
    as.character(getNamespaceVersion(package))
  }, character(1))
  matches <- vapply(packages, function(package) {
    utils::packageVersion(package) == package_version(expected[[package]]) &&
      getNamespaceVersion(package) == package_version(expected[[package]])
  }, logical(1))
  if (!all(matches)) {
    details <- paste0(
      packages[!matches], " ", versions[!matches],
      " (required ", expected[!matches], ")"
    )
    stop("Dependency versions do not match the reproduction lock: ",
      paste(details, collapse = ", "),
      ". Run Rscript scripts/setup.R in a fresh R session.",
      call. = FALSE
    )
  }
  versions
}

source_hashes <- function(root, files) {
  paths <- file.path(root, files)
  if (!all(file.exists(paths))) {
    stop("Missing provenance input: ", paste(files[!file.exists(paths)], collapse = ", "),
      call. = FALSE
    )
  }
  hashes <- tools::md5sum(paths)
  if (anyNA(hashes)) stop("Cannot read a provenance input.", call. = FALSE)
  names(hashes) <- files
  hashes
}

analysis_signature <- function(root, config) {
  files <- c("validation.R", "estimators.R", "design.R", "simulation.R", "pipeline.R")
  hashes <- source_hashes(root, file.path("R", files))
  names(hashes) <- files
  list(
    schema_version = 1L,
    config = config,
    numerical_source_md5 = hashes,
    R_version = as.character(getRversion()),
    platform = R.version$platform,
    packages = dependency_versions()
  )
}

output_manifest_signature <- function(root, config, mode = "live") {
  mode <- match.arg(mode)
  list(
    schema_version = 1L,
    mode = mode,
    analysis = analysis_signature(root, config),
    artifact_source_md5 = source_hashes(root, c(
      "R/bootstrap.R", "R/plot_simulations.R", "R/plot_applications.R",
      "scripts/reproduce.R", "scripts/setup.R"
    ))
  )
}

# Validate before writing any artifacts. A new manifest is written only under
# the invocation lock, and existing outputs are never silently adopted.
prepare_output <- function(root, output_dir, config, mode = "live",
                           create = TRUE) {
  signature <- output_manifest_signature(root, config, match.arg(mode))
  destination <- file.path(output_dir, "run_manifest.rds")
  if (file.exists(destination)) {
    prior <- readRDS(destination)
    if (!is.list(prior) || !identical(prior$signature, signature)) {
      stop("This output directory belongs to a different profile, mode, source, data, or dependency version. Choose a new --output directory.",
        call. = FALSE
      )
    }
  } else {
    contents <- setdiff(
      list.files(output_dir, all.files = TRUE, no.. = TRUE),
      ".reproduce-lock"
    )
    if (length(contents)) {
      stop("Output directory is nonempty but has no run manifest. Choose a new --output directory.",
        call. = FALSE
      )
    }
    if (create) {
      save_rds_atomic(
        list(
          signature = signature,
          created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
        ),
        destination
      )
    }
  }
  invisible(signature)
}

# Replace an output only after its new contents have been written successfully.
# All writers are called while an output-directory lock is held.
atomic_write <- function(path, writer) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = ".write-", tmpdir = dirname(path))
  backup <- tempfile(pattern = ".previous-", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  writer(temporary)
  existed <- file.exists(path)
  if (existed && !file.rename(path, backup)) {
    stop("Cannot replace output: ", path, call. = FALSE)
  }
  if (!file.rename(temporary, path)) {
    if (existed && !file.rename(backup, path)) {
      stop("Output replacement failed. Original data are preserved at: ", backup, call. = FALSE)
    }
    stop("Cannot finalize output: ", path, call. = FALSE)
  }
  if (existed) unlink(backup)
  invisible(path)
}

save_rds_atomic <- function(object, path) {
  atomic_write(path, function(destination) saveRDS(object, destination, version = 3))
}

write_csv_atomic <- function(object, path) {
  atomic_write(path, function(destination) {
    utils::write.csv(object, destination, row.names = FALSE, na = "")
  })
}

# Complete every requested plot and application calculation before replacing
# any existing artifact. A plotting error leaves all prior artifacts intact.
with_staged_artifacts <- function(output_dir, action) {
  staging <- tempfile(pattern = ".artifacts-", tmpdir = output_dir)
  if (!dir.create(staging)) stop("Cannot create artifact staging directory.", call. = FALSE)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  action(staging)
  files <- list.files(staging, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  for (relative in files) {
    source <- file.path(staging, relative)
    atomic_write(file.path(output_dir, relative), function(destination) {
      if (!file.copy(source, destination, overwrite = FALSE)) {
        stop("Cannot stage output: ", relative, call. = FALSE)
      }
    })
  }
  invisible(file.path(output_dir, files))
}

with_output_lock <- function(output_dir, action, lock_name = ".run-lock") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) stop("Cannot create output directory.", call. = FALSE)
  lock <- file.path(output_dir, lock_name)
  if (!dir.create(lock, showWarnings = FALSE)) {
    stop(
      "Output is locked: ", lock,
      ". Do not run overlapping writers. If a previous run crashed, confirm it has stopped before removing this lock.",
      call. = FALSE
    )
  }
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  writeLines(
    c(paste("process", Sys.getpid()), format(Sys.time(), tz = "UTC", usetz = TRUE)),
    file.path(lock, "owner.txt")
  )
  action()
}

validate_result <- function(result, cell, m_grid, n_grid, signature = NULL) {
  required <- c("cell", "m_grid", "grid", "loss", "diagnostics", "protocol")
  if (!is.list(result) || !all(required %in% names(result))) {
    stop("Incomplete simulation result.", call. = FALSE)
  }
  expected_columns <- c("NA", "Kernel", paste0("SM_m", m_grid))
  valid <- identical(result$cell, cell) && identical(result$m_grid, m_grid) &&
    isTRUE(all.equal(result$grid, seq(0, cell$tau, length.out = n_grid), tolerance = 0)) &&
    is.matrix(result$loss) && identical(colnames(result$loss), expected_columns) &&
    nrow(result$loss) == cell$B && all(is.finite(result$loss)) &&
    all(result$loss >= 0) && nrow(result$diagnostics) == cell$B
  if (!is.null(signature)) valid <- valid && identical(result$signature, signature)
  if (!valid) stop("Invalid or incompatible result for ", cell$cell_id, ".", call. = FALSE)
  invisible(TRUE)
}

prepare_run <- function(root, output_dir, config) {
  signature <- analysis_signature(root, config)
  protocol_file <- file.path(output_dir, "protocol.rds")
  if (file.exists(protocol_file)) {
    prior <- readRDS(protocol_file)
    if (!identical(prior$signature, signature)) {
      stop("This output directory belongs to a different configuration or code version. Choose a new --output directory.", call. = FALSE)
    }
  } else {
    if (length(list.files(output_dir, pattern = "_result\\.rds$"))) {
      stop("Found cell results without their protocol. Choose a new --output directory.", call. = FALSE)
    }
    save_rds_atomic(list(signature = signature, created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)), protocol_file)
  }
  signature
}

run_simulations <- function(root, output_dir, config, shard_index = 1L, shard_count = 1L) {
  if (length(shard_index) != 1L || length(shard_count) != 1L ||
    !is.finite(shard_index) || !is.finite(shard_count) ||
    shard_index != as.integer(shard_index) || shard_count != as.integer(shard_count) ||
    shard_count < 1 || shard_index < 1 || shard_index > shard_count) {
    stop("Require integer 1 <= shard_index <= shard_count.", call. = FALSE)
  }
  # Shards are separate directories. Aggregation reads them only after completion.
  cells <- make_cells(B = config$replications, tau = config$tau)
  selected <- which((seq_len(nrow(cells)) - 1L) %% shard_count == shard_index - 1L)
  if (!length(selected)) stop("This shard has no cells.", call. = FALSE)
  with_output_lock(output_dir, function() {
    shard_file <- file.path(output_dir, "shard.rds")
    shard <- list(index = shard_index, count = shard_count)
    if (file.exists(shard_file) && !identical(readRDS(shard_file), shard)) {
      stop("This output directory belongs to a different shard. Choose a separate output directory.", call. = FALSE)
    }
    signature <- prepare_run(root, output_dir, config)
    save_rds_atomic(shard, shard_file)
    bases <- poisson_bases(seq(0, config$tau, length.out = config$n_grid), config$m_grid)
    for (i in selected) {
      cell <- cells[i, , drop = FALSE]
      destination <- file.path(output_dir, paste0(cell$cell_id, "_result.rds"))
      if (file.exists(destination)) {
        validate_result(readRDS(destination), cell, config$m_grid, config$n_grid, signature)
        message("Reusing ", cell$cell_id)
        next
      }
      message("Simulating ", cell$cell_id, " (", cell$B, " replications)")
      result <- simulate_cell(cell, config$m_grid, n_grid = config$n_grid, bases = bases)
      result$signature <- signature
      validate_result(result, cell, config$m_grid, config$n_grid, signature)
      save_rds_atomic(result, destination)
    }
    write_csv_atomic(cells[selected, , drop = FALSE], file.path(output_dir, "cells.csv"))
    invisible(signature)
  })
}

write_summaries <- function(results, output_dir, n_boot = 1000L) {
  summaries <- lapply(results, summarize_cell, n_boot = n_boot)
  collect <- function(field) do.call(rbind, lapply(summaries, `[[`, field))
  write_csv_atomic(collect("summary"), file.path(output_dir, "summary_table.csv"))
  write_csv_atomic(collect("curve"), file.path(output_dir, "curve_summary.csv"))
  write_csv_atomic(collect("bootstrap"), file.path(output_dir, "bootstrap_minimum.csv"))
  save_rds_atomic(results, file.path(output_dir, "cell_results.rds"))
  save_rds_atomic(summaries, file.path(output_dir, "summaries.rds"))
  invisible(summaries)
}

aggregate_results <- function(root, output_dir, config, inputs = output_dir) {
  cells <- make_cells(B = config$replications, tau = config$tau)
  signature <- analysis_signature(root, config)
  for (directory in inputs) {
    if (dir.exists(file.path(directory, ".run-lock"))) {
      stop("A simulation writer is still active: ", directory, call. = FALSE)
    }
    protocol <- file.path(directory, "protocol.rds")
    if (!file.exists(protocol) || !identical(readRDS(protocol)$signature, signature)) {
      stop("Incompatible or missing input protocol: ", directory, call. = FALSE)
    }
  }
  with_output_lock(output_dir, function() {
    prepare_run(root, output_dir, config)
    results <- lapply(seq_len(nrow(cells)), function(i) {
      cell <- cells[i, , drop = FALSE]
      candidates <- file.path(inputs, paste0(cell$cell_id, "_result.rds"))
      available <- candidates[file.exists(candidates)]
      if (length(available) != 1L) {
        stop("Require exactly one completed result for ", cell$cell_id, ".", call. = FALSE)
      }
      result <- readRDS(available)
      validate_result(result, cell, config$m_grid, config$n_grid, signature)
      result
    })
    write_csv_atomic(cells, file.path(output_dir, "cells.csv"))
    write_summaries(results, output_dir, n_boot = config$n_boot)
    invisible(results)
  })
}
