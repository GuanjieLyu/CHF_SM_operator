# Recreate the introductory illustration and application figures.
#
# Load the estimator definitions before this module, then call:
# write_application_figures(figure_dir = "figures", results_dir = "results")
# For a fast check, use illustration_reps = 3L and kernel_grid = 401L.
# The default reads datasets from the installed KMsurv and survival packages.
# An explicit data_dir instead reads existing CSV files without modifying them.
# Nothing is downloaded or installed. Realdata_sensitivity.pdf is Figure B.2.

.app_plot_directory <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("Each output directory must be a nonempty path.", call. = FALSE)
  }
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE)) {
    stop("Cannot create output directory: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.app_plot_integer <- function(value, name, minimum) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
    value < minimum || value != floor(value) || value > .Machine$integer.max) {
    stop(name, " must be an integer of at least ", minimum, ".", call. = FALSE)
  }
  as.integer(value)
}

.app_plot_pdf <- function(path, height, draw) {
  previous_device <- grDevices::dev.cur()
  grDevices::pdf(path,
    width = 6.3, height = height, pointsize = 11,
    family = "Helvetica", useDingbats = FALSE
  )
  output_device <- grDevices::dev.cur()
  on.exit(
    {
      if (output_device %in% grDevices::dev.list()) {
        grDevices::dev.off(output_device)
      }
      if (previous_device > 1L && previous_device %in% grDevices::dev.list()) {
        grDevices::dev.set(previous_device)
      }
    },
    add = TRUE
  )
  draw()
  invisible(path)
}

.app_plot_seed <- function(seed, compute) {
  previous_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(
    {
      do.call(RNGkind, as.list(previous_kind))
      if (had_seed) {
        assign(".Random.seed", previous_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    },
    add = TRUE
  )
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  compute()
}

.app_plot_par <- function(mar = c(4.2, 4.5, 2, 0.8)) {
  graphics::par(
    mar = mar, mgp = c(2.65, 0.7, 0), tcl = -0.25, las = 1,
    family = "Helvetica", cex = 1, cex.axis = 1, cex.lab = 1,
    cex.main = 1, font.main = 2, xaxs = "i", yaxs = "i"
  )
}

.app_plot_frame <- function(xmax, ymax, xlab, title) {
  graphics::plot(
    NA_real_, NA_real_,
    xlim = c(0, xmax), ylim = c(0, ymax),
    xlab = xlab, ylab = "Cumulative hazard", main = title, bty = "l"
  )
  graphics::abline(h = pretty(c(0, ymax)), col = "#E6E6E6", lwd = 0.5)
}

.app_plot_na <- function(fit, xmax, col = "#252525", lwd = 1.4, lty = 1) {
  selected <- fit$time <= xmax
  times <- fit$time[selected]
  values <- fit$cumhaz[selected]
  endpoint <- if (length(values)) tail(values, 1L) else 0
  graphics::lines(c(0, times, xmax), c(0, values, endpoint),
    type = "s",
    col = col, lwd = lwd, lty = lty
  )
}

.app_plot_check_curve <- function(values, label) {
  if (any(!is.finite(values)) || any(values < -1e-9) ||
    any(diff(values) < -1e-8) || abs(values[1]) > 1e-12) {
    stop("Invalid cumulative hazard curve: ", label, call. = FALSE)
  }
}

.app_plot_read_data <- function(data_dir = NULL) {
  if (is.null(data_dir)) {
    required <- c("KMsurv", "survival")
    available <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
    if (!all(available)) {
      stop(
        "The default application datasets require the installed packages: ",
        paste(required[!available], collapse = ", "),
        ". Install them first or supply data_dir with existing CSV files.",
        call. = FALSE
      )
    }
    data_environment <- new.env(parent = baseenv())
    utils::data("larynx", package = "KMsurv", envir = data_environment)
    larynx <- data_environment$larynx
    pbc <- survival::pbc[, c("id", "time", "status")]
  } else {
    if (!is.character(data_dir) || length(data_dir) != 1L ||
      is.na(data_dir) || !nzchar(data_dir)) {
      stop("data_dir must be NULL or a nonempty directory path.", call. = FALSE)
    }
    paths <- file.path(data_dir, c("larynx.csv", "pbc_analysis.csv"))
    if (!all(file.exists(paths))) {
      stop("Missing input CSV files: ", paste(paths[!file.exists(paths)], collapse = ", "),
        call. = FALSE
      )
    }
    larynx <- utils::read.csv(paths[1], stringsAsFactors = FALSE)
    pbc <- utils::read.csv(paths[2], stringsAsFactors = FALSE)
  }
  if (!all(c("time", "delta") %in% names(larynx)) ||
    !all(c("time", "status") %in% names(pbc))) {
    stop("Expected time/delta in the laryngeal data and time/status in the PBC data.",
      call. = FALSE
    )
  }
  if (nrow(larynx) != 90L || anyNA(larynx$delta) ||
    !all(larynx$delta %in% c(0, 1)) || sum(larynx$delta == 0) != 40L ||
    nrow(pbc) != 418L || anyNA(pbc$status) ||
    !all(pbc$status %in% c(0, 1, 2)) || sum(pbc$status == 0) != 232L ||
    sum(pbc$status == 1) != 25L || sum(pbc$status == 2) != 161L) {
    stop("Application data do not match the documented cohorts.",
      call. = FALSE
    )
  }
  if (any(!is.finite(larynx$time)) || any(larynx$time <= 0) ||
    any(!is.finite(pbc$time)) || any(pbc$time <= 0)) {
    stop("Application follow-up times must be finite and positive.", call. = FALSE)
  }
  # Death and transplantation form the prespecified composite PBC endpoint.
  pbc$delta <- as.integer(pbc$status != 0)
  list(larynx = larynx, pbc = pbc)
}

.app_plot_fit <- function(data, name, unit, m, m_values, reference_time,
                          prefix, kernel_grid) {
  x <- data$time
  delta <- data$delta
  fit <- grouped_na(x, delta)
  xmax <- max(x)
  grid <- sort(unique(c(seq(0, xmax, length.out = kernel_grid), fit$time)))
  na <- na_curve(grid, fit)
  sm <- sm_curve(grid, fit, m)
  kernel <- kernel_curve(x, delta, grid, n_grid = kernel_grid)
  .app_plot_check_curve(na, paste(name, "Nelson-Aalen"))
  .app_plot_check_curve(sm, paste(name, "Szasz-Mirakyan"))
  .app_plot_check_curve(kernel$cumhaz, paste(name, "integrated kernel"))
  sensitivity <- vapply(
    m_values, function(value) sm_curve(grid, fit, value),
    numeric(length(grid))
  )
  for (j in seq_along(m_values)) {
    .app_plot_check_curve(sensitivity[, j], paste(name, "m =", m_values[j]))
  }
  curves <- data.frame(
    time = grid, Nelson_Aalen = na, integrated_kernel = kernel$cumhaz,
    Szasz_Mirakyan = sm
  )
  for (j in seq_along(m_values)) {
    curves[[paste0("SM_m", m_values[j])]] <- sensitivity[, j]
  }
  summary <- data.frame(
    dataset = name, time_unit = unit, n = length(x), events = sum(delta),
    censored = sum(delta == 0), censor_fraction = mean(delta == 0),
    distinct_observed_times = length(unique(x)),
    tied_observed_time_groups = sum(table(x) > 1),
    distinct_failure_times = length(fit$time),
    tied_failure_time_groups = sum(table(x[delta == 1]) > 1),
    mixed_event_censor_time_groups = length(intersect(
      unique(x[delta == 1]), unique(x[delta == 0])
    )),
    min_time = min(x), max_time = xmax, m = m, reference_time = reference_time,
    evaluation_sd = sqrt(reference_time / m), kernel_fit_max = kernel$fit_max,
    kernel_extension_used = kernel$extension_used,
    kernel_warning_count = length(kernel$warnings),
    kernel_warnings = paste(kernel$warnings, collapse = " | "),
    evaluation_points = length(grid),
    max_absolute_SM_NA_difference = max(abs(sm - na))
  )
  sensitivity_summary <- data.frame(
    dataset = name, m = m_values, reference_time = reference_time,
    evaluation_sd = sqrt(reference_time / m_values),
    max_absolute_SM_NA_difference = apply(abs(sensitivity - na), 2, max)
  )
  # This independent identity check is inexpensive. A second, doubled kernel
  # grid is deliberately not fitted by the plotting function.
  check_times <- sort(unique(c(seq(0, xmax, length.out = 101L), fit$time)))
  gamma_values <- sm_curve(check_times, fit, m)
  poisson_values <- as.vector(outer(
    check_times, ceiling(m * fit$time),
    function(time, shape) stats::ppois(shape - 1, m * time, lower.tail = FALSE)
  ) %*% fit$weight)
  identity_error <- max(abs(gamma_values - poisson_values))
  if (identity_error >= 1e-12) {
    stop("Gamma and Poisson representations disagree for ", name, call. = FALSE)
  }
  na_error <- NA_real_
  if (requireNamespace("survival", quietly = TRUE)) {
    reference <- survival::survfit(
      survival::Surv(time, delta) ~ 1,
      data = data, ctype = 1
    )
    na_error <- max(abs(
      fit$cumhaz - reference$cumhaz[match(fit$time, reference$time)]
    ))
    if (!is.finite(na_error) || na_error >= 1e-12) {
      stop("Nelson-Aalen reference check failed for ", name, call. = FALSE)
    }
  }
  list(
    name = name, unit = unit, m = m, m_values = m_values, prefix = prefix,
    xmax = xmax, fit = fit, time = grid, na = na, sm = sm,
    kernel = kernel, sensitivity = sensitivity, curves = curves,
    summary = summary, sensitivity_summary = sensitivity_summary,
    validation = data.frame(
      dataset = prefix, NA_vs_survfit_max_error = na_error,
      gamma_vs_Poisson_max_error = identity_error,
      kernel_grid_points = kernel_grid, kernel_refinement_checked = FALSE
    )
  )
}

.app_plot_main <- function(result, path) {
  .app_plot_pdf(path, 4.4, function() {
    graphics::layout(matrix(1:2, ncol = 1), heights = c(0.52, 3.88))
    graphics::par(mar = rep(0, 4), cex = 1)
    graphics::plot.new()
    graphics::legend(
      "center", c(
        "Nelson-Aalen", "Integrated kernel",
        paste0("Szasz-Mirakyan (m = ", result$m, ")")
      ),
      col = c("#252525", "#D55E00", "#0072B2"), lty = c(1, 2, 3),
      lwd = c(1.4, 1.7, 1.9), horiz = TRUE, bty = "n",
      x.intersp = 0.55, seg.len = 1.75, cex = 0.94
    )
    .app_plot_par(c(4.2, 4.5, 1.8, 0.8))
    .app_plot_frame(
      result$xmax, max(result$na, result$sm, result$kernel$cumhaz) * 1.07,
      paste0("Time (", result$unit, ")"), result$name
    )
    .app_plot_na(result$fit, result$xmax)
    graphics::lines(result$time, result$kernel$cumhaz,
      col = "#D55E00", lty = 2, lwd = 1.7
    )
    graphics::lines(result$time, result$sm, col = "#0072B2", lty = 3, lwd = 1.9)
  })
}

.app_plot_sensitivity <- function(results, path) {
  .app_plot_pdf(path, 6.8, function() {
    graphics::layout(matrix(1:4, ncol = 1), heights = c(0.43, 2.97, 0.43, 2.97))
    colors <- c("#009E73", "#0072B2", "#D55E00")
    patterns <- c(2, 3, 4)
    for (result in results) {
      graphics::par(mar = rep(0, 4), cex = 1)
      graphics::plot.new()
      graphics::legend(
        "center", c("Nelson-Aalen", paste0("m = ", result$m_values)),
        col = c("#252525", colors), lty = c(1, patterns), lwd = c(1.4, rep(1.8, 3)),
        horiz = TRUE, bty = "n", x.intersp = 0.65, seg.len = 1.8
      )
      .app_plot_par(c(4.1, 4.5, 1.75, 0.8))
      .app_plot_frame(
        result$xmax, max(result$na, result$sensitivity) * 1.07,
        paste0("Time (", result$unit, ")"), result$name
      )
      .app_plot_na(result$fit, result$xmax)
      for (j in seq_along(result$m_values)) {
        graphics::lines(result$time, result$sensitivity[, j],
          col = colors[j], lty = patterns[j], lwd = 1.8
        )
      }
    }
  })
}

.app_plot_illustration_data <- function(n_rep, kernel_grid) {
  seed <- 2026110701L
  .app_plot_seed(seed, function() {
    n <- 30L
    tau <- 3.7
    grid <- seq(0, tau, length.out = 801L)
    truth <- (grid / 2)^2
    target_censor <- 0.30
    probability_censored <- function(rate) {
      stats::integrate(function(time) {
        (1 - exp(-rate * time)) * stats::dweibull(time, shape = 2, scale = 2)
      }, lower = 0, upper = Inf, rel.tol = 1e-10)$value
    }
    censor_rate <- stats::uniroot(
      function(rate) probability_censored(rate) - target_censor,
      interval = c(1e-8, 10), tol = 1e-12
    )$root
    estimates <- vector("list", n_rep)
    samples <- vector("list", n_rep)
    curves <- vector("list", n_rep)
    summaries <- vector("list", n_rep)
    for (replicate in seq_len(n_rep)) {
      event_time <- stats::rweibull(n, shape = 2, scale = 2)
      censor_time <- stats::rexp(n, rate = censor_rate)
      observed_time <- pmin(event_time, censor_time)
      delta <- as.integer(event_time <= censor_time)
      fit <- grouped_na(observed_time, delta)
      na <- na_curve(grid, fit)
      sm <- sm_curve(grid, fit, 20)
      kernel <- kernel_curve(observed_time, delta, grid, n_grid = kernel_grid)
      .app_plot_check_curve(na, paste("illustration", replicate, "Nelson-Aalen"))
      .app_plot_check_curve(sm, paste("illustration", replicate, "Szasz-Mirakyan"))
      .app_plot_check_curve(kernel$cumhaz, paste("illustration", replicate, "kernel"))
      estimates[[replicate]] <- list(fit = fit, na = na, sm = sm, kernel = kernel$cumhaz)
      samples[[replicate]] <- data.frame(
        replicate = replicate, subject = seq_len(n), event_time = event_time,
        censor_time = censor_time, observed_time = observed_time, delta = delta
      )
      curves[[replicate]] <- data.frame(
        replicate = replicate, time = grid, truth = truth, Nelson_Aalen = na,
        integrated_kernel = kernel$cumhaz, Szasz_Mirakyan = sm
      )
      summaries[[replicate]] <- data.frame(
        replicate = replicate, n = n, events = sum(delta),
        censor_fraction = mean(delta == 0), kernel_fit_max = kernel$fit_max,
        kernel_extension_used = kernel$extension_used,
        kernel_warning_count = length(kernel$warnings),
        kernel_warnings = paste(kernel$warnings, collapse = " | ")
      )
    }
    summary <- do.call(rbind, summaries)
    list(
      estimates = estimates, time = grid, truth = truth, tau = tau,
      samples = do.call(rbind, samples), curves = do.call(rbind, curves),
      summary = summary, design = data.frame(
        seed = seed, replicates = n_rep, n = n, Weibull_shape = 2, Weibull_scale = 2,
        target_censor_fraction = target_censor, censor_distribution = "exponential",
        censor_rate = censor_rate, calibrated_censor_fraction = probability_censored(censor_rate),
        tau = tau, evaluation_points = length(grid), kernel_grid_points = kernel_grid,
        m = 20, mean_realized_censor_fraction = mean(summary$censor_fraction)
      )
    )
  })
}

.app_plot_illustration <- function(result, path) {
  ymax <- max(result$truth, unlist(lapply(result$estimates, function(estimate) {
    c(estimate$na, estimate$sm, estimate$kernel)
  }))) * 1.04
  .app_plot_pdf(path, 7.5, function() {
    graphics::layout(matrix(1:4, ncol = 1), heights = c(0.45, 2.35, 2.35, 2.35))
    graphics::par(mar = rep(0, 4), cex = 1)
    graphics::plot.new()
    graphics::legend(
      "center", c("True cumulative hazard", sprintf(
        "Individual estimates (%d samples)", length(result$estimates)
      )),
      col = c("#0072B2", "#777777"), lty = c(1, 2), lwd = c(2.4, 0.9),
      horiz = TRUE, bty = "n", x.intersp = 0.65, seg.len = 2.2
    )
    titles <- c(
      "(a) Nelson-Aalen", "(b) Integrated local kernel", "(c) Szasz-Mirakyan (m = 20)"
    )
    for (panel in seq_len(3L)) {
      .app_plot_par(c(3.65, 4.5, 1.8, 0.8))
      .app_plot_frame(result$tau, ymax, "Time", titles[panel])
      for (estimate in result$estimates) {
        if (panel == 1L) {
          .app_plot_na(estimate$fit, result$tau, col = "#77777788", lwd = 0.75, lty = 2)
        } else {
          values <- if (panel == 2L) estimate$kernel else estimate$sm
          graphics::lines(result$time, values, col = "#77777788", lwd = 0.75, lty = 2)
        }
      }
      graphics::lines(result$time, result$truth, col = "#0072B2", lwd = 2.4)
    }
  })
}

#' Write the application and introductory illustration figures.
#'
#' @param data_dir NULL by default, using larynx from the installed KMsurv
#'   package and pbc from survival. Alternatively, a directory containing
#'   larynx.csv and pbc_analysis.csv. Input files are never modified.
#' @param figure_dir Directory for Fig1.pdf, Fig5.pdf, Fig6.pdf and
#'   Realdata_sensitivity.pdf (the manuscript's Figure B.2).
#' @param results_dir Directory for plotted curves, cohort summaries,
#'   illustration samples and inexpensive identity-validation results.
#' @param illustration_reps Number of shared Figure 1 datasets. The paper uses
#'   50 datasets of 30 observations, with the fixed seed 2026110701.
#' @param kernel_grid Application fitting and evaluation grid size. The paper
#'   uses 4001. Figure 1 uses min(kernel_grid, 2001), preserving its published
#'   2001-point fitting grid while allowing faster smoke checks.
#' @return Invisibly, a list of generated figure and numerical-result paths.
#' @details The caller's RNG kind, RNG state and active graphics device are
#'   restored, including after an error. No finer-grid kernel check is run
#'   implicitly. Smoothing values are positive real parameters with no rounding.
write_application_figures <- function(data_dir = NULL, figure_dir, results_dir,
                                      illustration_reps = 50L, kernel_grid = 4001L) {
  illustration_reps <- .app_plot_integer(illustration_reps, "illustration_reps", 1L)
  kernel_grid <- .app_plot_integer(kernel_grid, "kernel_grid", 3L)
  for (name in c("grouped_na", "na_curve", "sm_curve", "kernel_curve")) {
    if (!is.function(get0(name, mode = "function", inherits = TRUE))) {
      stop("Load the estimator module before plotting. Missing function: ", name,
        call. = FALSE
      )
    }
  }
  data <- .app_plot_read_data(data_dir)
  figure_dir <- .app_plot_directory(figure_dir)
  results_dir <- .app_plot_directory(results_dir)
  results <- list(
    .app_plot_fit(
      data$larynx, "Laryngeal cancer", "years", 15, c(8, 15, 30),
      5, "laryngeal", kernel_grid
    ),
    .app_plot_fit(
      data$pbc, "Primary biliary cholangitis", "days", 50,
      c(25, 50, 100), 2500, "pbc", kernel_grid
    )
  )
  figure_paths <- file.path(
    figure_dir, c("Fig1.pdf", "Fig5.pdf", "Fig6.pdf", "Realdata_sensitivity.pdf")
  )
  .app_plot_main(results[[1]], figure_paths[2])
  .app_plot_main(results[[2]], figure_paths[3])
  .app_plot_sensitivity(results, figure_paths[4])
  illustration <- .app_plot_illustration_data(illustration_reps, min(kernel_grid, 2001L))
  .app_plot_illustration(illustration, figure_paths[1])
  numerical_paths <- character()
  write_result <- function(data, filename) {
    path <- file.path(results_dir, filename)
    utils::write.csv(data, path, row.names = FALSE)
    path
  }
  for (result in results) {
    numerical_paths <- c(
      numerical_paths,
      write_result(result$curves, paste0(result$prefix, "_curves.csv")),
      write_result(as.data.frame(result$fit), paste0(result$prefix, "_grouped_NA.csv")),
      write_result(data.frame(
        time = result$kernel$fit_grid, hazard = result$kernel$fit_hazard,
        cumhaz = result$kernel$fit_cumhaz
      ), paste0(result$prefix, "_kernel_native_grid.csv"))
    )
  }
  numerical_paths <- c(
    numerical_paths,
    write_result(do.call(rbind, lapply(results, function(x) x$summary)), "real_data_summary.csv"),
    write_result(
      do.call(rbind, lapply(results, function(x) x$sensitivity_summary)),
      "real_data_sensitivity_summary.csv"
    ),
    write_result(do.call(rbind, lapply(results, function(x) x$validation)), "real_data_validation.csv"),
    write_result(illustration$samples, "illustration_samples.csv"),
    write_result(illustration$curves, "illustration_curves.csv"),
    write_result(illustration$summary, "illustration_replicate_summary.csv"),
    write_result(illustration$design, "illustration_design.csv")
  )
  invisible(list(figures = figure_paths, results = numerical_paths))
}
