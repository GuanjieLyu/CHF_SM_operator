# Common-sample Monte Carlo losses and descriptive oracle-grid summaries.
# Source validation.R, estimators.R, and design.R before this file.

# Precompute Poisson weights. omitted_mass bounds the missing probability
# uniformly over grid, and multiplying it by the total NA increment bounds
# the absolute error in the smoothed cumulative hazard.
poisson_bases <- function(grid, m_grid, tail_tol = 1e-12) {
  .validate_grid(grid)
  .validate_m_grid(m_grid)
  .assert_scalar(tail_tol, "tail_tol", lower = 0, strictly = TRUE)
  if (tail_tol >= 1 || 1 - tail_tol == 1) {
    stop("tail_tol must be below 1 and large enough that 1 - tail_tol < 1.",
      call. = FALSE
    )
  }
  lapply(m_grid, function(m) {
    lambda_max <- m * max(grid)
    if (!is.finite(lambda_max)) {
      stop("Each m times max(grid) must be finite.", call. = FALSE)
    }
    endpoint <- qpois(1 - tail_tol, lambda_max)
    if (!is.finite(endpoint) || endpoint > .Machine$integer.max - 1) {
      stop("The requested Poisson basis exceeds the supported index range.",
        call. = FALSE
      )
    }
    jmax <- as.integer(endpoint)
    j <- 0:jmax
    list(
      m = m, j = j,
      P = outer(grid, j, function(t, k) dpois(k, m * t)),
      omitted_mass = ppois(jmax, lambda_max, lower.tail = FALSE)
    )
  })
}

.validate_bases <- function(bases, grid, m_grid) {
  if (!is.list(bases) || length(bases) != length(m_grid)) {
    stop("bases must contain one poisson_bases() entry per m_grid value.",
      call. = FALSE
    )
  }
  for (q in seq_along(m_grid)) {
    base <- bases[[q]]
    if (!is.list(base) || !all(c("m", "j", "P", "omitted_mass") %in% names(base))) {
      stop("Each basis must have m, j, P, and omitted_mass fields.", call. = FALSE)
    }
    .assert_scalar(base$m, "basis m", lower = 0, strictly = TRUE)
    .assert_numeric_vector(base$j, "basis j",
      min_length = 1L,
      nonnegative = TRUE
    )
    .assert_scalar(base$omitted_mass, "basis omitted_mass", lower = 0, upper = 1)
    if (base$m != m_grid[q] || !all(base$j == seq_along(base$j) - 1L) ||
      !is.matrix(base$P) || !is.numeric(base$P) ||
      !identical(dim(base$P), c(length(grid), length(base$j))) ||
      any(!is.finite(base$P)) || any(base$P < 0)) {
      stop("bases must match m_grid and the evaluation grid dimensions.",
        call. = FALSE
      )
    }
    if (!isTRUE(all.equal(base$P[, 1L], dpois(0, base$m * grid),
      check.attributes = FALSE
    ))) {
      stop("The supplied Poisson basis uses different evaluation coordinates.",
        call. = FALSE
      )
    }
  }
  invisible(bases)
}

.validate_output_dir <- function(output_dir) {
  if (!is.null(output_dir) &&
    (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !dir.exists(output_dir))) {
    stop("output_dir must be NULL or an existing directory.", call. = FALSE)
  }
  invisible(output_dir)
}

# Generate each dataset once, retain all datasets, and evaluate every method
# on the same fixed interval. Seeded calls restore the caller's RNG state,
# including when estimation fails. Optional snapshots use the reference schema.
simulate_cell <- function(cell, m_grid = c(1, seq(2, 100, 2)), n_grid = 2001L,
                          kernel_n_grid = n_grid, output_dir = NULL, bases = NULL) {
  .validate_cell(cell)
  .validate_m_grid(m_grid)
  .assert_scalar(n_grid, "n_grid",
    lower = 2, integer = TRUE,
    upper = .Machine$integer.max
  )
  .assert_scalar(kernel_n_grid, "kernel_n_grid",
    lower = 3, integer = TRUE,
    upper = .Machine$integer.max
  )
  .validate_output_dir(output_dir)
  grid <- seq(0, cell$tau, length.out = n_grid)
  if (!is.null(bases)) {
    .validate_bases(bases, grid, m_grid)
  }
  .with_reproducible_rng(cell$seed, {
    samples <- replicate(
      cell$B, generate_sample(cell$n, cell$scenario, cell$gamma),
      simplify = FALSE
    )
    if (!is.null(output_dir)) {
      saveRDS(
        list(cell = cell, samples = samples),
        file.path(output_dir, paste0(cell$cell_id, "_samples.rds"))
      )
    }
    fits <- lapply(samples, function(z) grouped_na(z$x, z$delta))
    weights <- trap_weights(grid)
    truth <- true_chf(grid, cell$scenario)
    na <- vapply(fits, function(f) na_curve(grid, f), numeric(length(grid)))
    kernel <- lapply(samples, function(z) {
      kernel_curve(z$x, z$delta, grid, kernel_n_grid)
    })
    kmat <- vapply(kernel, function(z) z$cumhaz, numeric(length(grid)))
    loss <- matrix(
      NA_real_, cell$B, 2L + length(m_grid),
      dimnames = list(NULL, c("NA", "Kernel", paste0("SM_m", m_grid)))
    )
    loss[, "NA"] <- colSums((na - truth)^2 * weights)
    loss[, "Kernel"] <- colSums((kmat - truth)^2 * weights)
    if (is.null(bases)) {
      bases <- poisson_bases(grid, m_grid)
    }
    bounds <- numeric(length(m_grid))
    for (q in seq_along(m_grid)) {
      base <- bases[[q]]
      h <- vapply(fits, function(f) {
        na_curve(base$j / base$m, f)
      }, numeric(length(base$j)))
      # A one-row basis is still a matrix (vapply otherwise simplifies it).
      if (is.null(dim(h))) {
        h <- matrix(h, nrow = length(base$j))
      }
      sm <- base$P %*% h
      loss[, q + 2L] <- colSums((sm - truth)^2 * weights)
      bounds[q] <- max(vapply(fits, function(f) sum(f$weight), 0)) *
        base$omitted_mass
    }
    if (any(!is.finite(loss))) {
      stop("A replicate loss is nonfinite. No sample is discarded or retried.",
        call. = FALSE
      )
    }
    diagnostics <- data.frame(
      replicate = seq_len(cell$B),
      observed_censor = vapply(samples, function(z) mean(z$delta == 0), 0),
      max_observed = vapply(samples, function(z) max(z$x), 0),
      kernel_fit_max = vapply(kernel, function(z) z$fit_max, 0),
      kernel_tail_extension = vapply(kernel, function(z) z$extension_used, FALSE),
      no_events = vapply(kernel, function(z) z$no_events, FALSE),
      kernel_warning = vapply(kernel, function(z) paste(z$warnings, collapse = "; "), "")
    )
    diagnostics$kernel_min_hazard <- vapply(kernel, function(z) min(z$fit_hazard), 0)
    kernel_version <- tryCatch(
      as.character(utils::packageVersion("muhaz")),
      error = function(e) NA_character_
    )
    ans <- list(
      cell = cell, m_grid = m_grid, grid = grid, loss = loss,
      diagnostics = diagnostics, max_sm_uniform_truncation_bound = max(bounds),
      protocol = list(
        quadrature = "Common fixed-grid composite trapezoid", n_grid = n_grid,
        kernel_n_grid = kernel_n_grid,
        sample_retention = "All generated samples; no rejection",
        kernel = list(
          package = "muhaz", version = kernel_version,
          bw.method = "local", b.cor = "both", kern = "epanechnikov",
          n.min.grid = 51L, min.time = 0,
          max.time = "min(tau,max observed time)",
          ties = "muhaz native sequential increments; failures before censorings",
          tail = "Constant integrated cumulative hazard beyond fitted endpoint"
        )
      )
    )
    if (!is.null(output_dir)) {
      saveRDS(ans, file.path(output_dir, paste0(cell$cell_id, "_result.rds")))
    }
    ans
  })
}

.validate_result <- function(result) {
  required <- c(
    "cell", "m_grid", "loss", "diagnostics",
    "max_sm_uniform_truncation_bound"
  )
  if (!is.list(result) || !all(required %in% names(result))) {
    stop("result must be a simulate_cell() result.", call. = FALSE)
  }
  .validate_cell(result$cell)
  .validate_m_grid(result$m_grid)
  z <- result$loss
  expected <- c("NA", "Kernel", paste0("SM_m", result$m_grid))
  if (!is.matrix(z) || !is.numeric(z) || nrow(z) != result$cell$B ||
    !identical(colnames(z), expected) || any(!is.finite(z)) || any(z < 0)) {
    stop("result$loss must be the finite, nonnegative paired-loss matrix.",
      call. = FALSE
    )
  }
  fields <- c(
    "observed_censor", "kernel_tail_extension", "no_events",
    "kernel_warning", "kernel_min_hazard"
  )
  diagnostics <- result$diagnostics
  if (!is.data.frame(diagnostics) || nrow(diagnostics) != nrow(z) ||
    !all(fields %in% names(diagnostics)) || anyNA(diagnostics[fields])) {
    stop("result$diagnostics must have a complete row for each replication.",
      call. = FALSE
    )
  }
  .assert_scalar(result$max_sm_uniform_truncation_bound,
    "max_sm_uniform_truncation_bound",
    lower = 0
  )
  invisible(result)
}

# Monte Carlo means and fixed-index/paired standard errors. The empirical
# minimum is selected using these same losses. Its fixed-m MCSE ignores this
# selection. The paired row bootstrap reselects m in every resample and reports
# descriptive spread, not confidence intervals for the population optimum.
summarize_cell <- function(result, n_boot = 1000L) {
  .validate_result(result)
  .assert_scalar(n_boot, "n_boot",
    lower = 2, integer = TRUE,
    upper = .Machine$integer.max
  )
  z <- result$loss
  B <- nrow(z)
  mm <- result$m_grid
  cell <- result$cell
  means <- colMeans(z)
  se <- apply(z, 2, sd) / sqrt(B)
  best <- which.min(means[-c(1, 2)]) + 2L
  diff_na <- z[, best] - z[, "NA"]
  diff_ker <- z[, best] - z[, "Kernel"]
  delta <- sweep(z[, -c(1, 2), drop = FALSE], 1, z[, best], "-")
  delta_mean <- colMeans(delta)
  delta_se <- apply(delta, 2, sd) / sqrt(B)
  near <- mm[delta_mean <= 1.96 * delta_se + 1e-14]
  .with_reproducible_rng(cell$seed + 500L, {
    boot_weights <- rmultinom(n_boot, B, rep(1 / B, B)) / B
    boot_means <- t(z[, -c(1, 2), drop = FALSE]) %*% boot_weights
    boot_index <- max.col(-t(boot_means), ties.method = "first")
    boot_min <- boot_means[cbind(boot_index, seq_len(n_boot))]
    summary <- data.frame(
      cell,
      NA_IMSE = means[1], NA_MCSE = se[1],
      Kernel_IMSE = means[2], Kernel_MCSE = se[2],
      m_grid_star = mm[best - 2L], SM_min_IMSE = means[best],
      SM_fixed_m_MCSE = se[best], SM_reselected_bootstrap_MCSE = sd(boot_min),
      SM_reselected_bootstrap_lo = unname(quantile(boot_min, 0.025)),
      SM_reselected_bootstrap_hi = unname(quantile(boot_min, 0.975)),
      SM_minus_NA = mean(diff_na), SM_minus_NA_MCSE = sd(diff_na) / sqrt(B),
      SM_minus_Kernel = mean(diff_ker),
      SM_minus_Kernel_MCSE = sd(diff_ker) / sqrt(B),
      reduction_vs_NA = 100 * (means[1] - means[best]) / means[1],
      reduction_vs_Kernel = 100 * (means[2] - means[best]) / means[2],
      near_optimal_descriptive_m = paste(near, collapse = ","),
      observed_censor = mean(result$diagnostics$observed_censor),
      tail_extension_count = sum(result$diagnostics$kernel_tail_extension),
      tail_extension_fraction = mean(result$diagnostics$kernel_tail_extension),
      zero_event_count = sum(result$diagnostics$no_events),
      kernel_warning_count = sum(nzchar(result$diagnostics$kernel_warning)),
      min_kernel_hazard = min(result$diagnostics$kernel_min_hazard),
      max_SM_truncation_bound = result$max_sm_uniform_truncation_bound,
      row.names = NULL, check.names = FALSE
    )
    curve <- data.frame(
      cell_id = cell$cell_id, scenario = cell$scenario, n = cell$n,
      target_censor = cell$target_censor,
      method = c("NA", "Kernel", rep("SM", length(mm))),
      m = c(NA, NA, mm), IMSE = unname(means), MCSE = unname(se),
      mean_difference_to_empirical_SM_min = c(NA, NA, unname(delta_mean)),
      paired_MCSE_to_empirical_SM_min = c(NA, NA, unname(delta_se))
    )
    boot <- data.frame(
      cell_id = cell$cell_id, m = mm,
      selection_frequency = tabulate(boot_index, nbins = length(mm)) / n_boot
    )
    list(
      summary = summary, curve = curve, bootstrap = boot,
      boot_minimum_values = boot_min
    )
  })
}
