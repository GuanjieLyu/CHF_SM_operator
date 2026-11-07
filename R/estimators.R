# Nelson--Aalen, Szasz--Mirakyan, and integrated kernel estimators.
# Source validation.R before this file. Only kernel_curve() needs muhaz.

# Group failures at equal observed times. Subjects censored at the same time
# remain in its risk set. The result is also defined for a sample with no events.
grouped_na <- function(x, delta) {
  .validate_observations(x, delta)
  z <- sort(unique(x[delta == 1]))
  d <- vapply(z, function(t) sum(delta[x == t]), 0)
  risk <- vapply(z, function(t) sum(x >= t), 0L)
  a <- d / risk
  list(
    time = z, eventcount = d, risk = risk, weight = a,
    cumhaz = cumsum(a), n = length(x)
  )
}

# Evaluate the exact, right-continuous Nelson--Aalen step function.
na_curve <- function(t, fit) {
  .assert_numeric_vector(t, "t", nonnegative = TRUE)
  .validate_na_fit(fit)
  c(0, fit$cumhaz)[findInterval(t, fit$time) + 1L]
}

# Exact finite gamma-CDF representation, valid for every real m > 0.
# Positive observed failure times are required for the gamma shapes.
sm_curve <- function(t, fit, m) {
  .assert_numeric_vector(t, "t", nonnegative = TRUE)
  .validate_na_fit(fit)
  .assert_scalar(m, "m", lower = 0, strictly = TRUE)
  if (!length(fit$time)) {
    return(rep(0, length(t)))
  }
  shape <- ceiling(m * fit$time)
  if (any(!is.finite(shape)) || any(shape < 1)) {
    stop("m times each failure time must be finite and strictly positive.",
      call. = FALSE
    )
  }
  as.vector(outer(t, shape, function(tt, kk) {
    pgamma(tt, shape = kk, rate = m)
  }) %*% fit$weight)
}

# Cumulative trapezoidal integration on the supplied coordinates.
cumtrap <- function(x, y) {
  .validate_grid(x, "x")
  .assert_numeric_vector(y, "y", min_length = length(x))
  if (length(y) != length(x)) {
    stop("x and y must have the same length.", call. = FALSE)
  }
  c(0, cumsum(diff(x) * (head(y, -1) + tail(y, -1)) / 2))
}

# Weights whose scalar product with function values is their trapezoidal integral.
trap_weights <- function(t) {
  .validate_grid(t, "t")
  d <- diff(t)
  c(d[1] / 2, (head(d, -1) + tail(d, -1)) / 2, tail(d, 1) / 2)
}

# Local Epanechnikov hazard with both boundary corrections, integrated on the
# actual muhaz coordinates. The cumulative hazard is constant after fit_max.
# muhaz retains sequential risk weights within tied failures. Sorting failures
# before censorings makes this comparator deterministic but does not replace
# those native weights with the grouped increments used by NA and SM.
kernel_curve <- function(x, delta, t, n_grid = 2001L) {
  .validate_observations(x, delta)
  if (any(x[delta == 1] <= 0)) {
    stop("Kernel hazard estimation requires strictly positive failure times.",
      call. = FALSE
    )
  }
  .assert_numeric_vector(t, "t", min_length = 1L, nonnegative = TRUE)
  .assert_scalar(n_grid, "n_grid",
    lower = 3, integer = TRUE,
    upper = .Machine$integer.max
  )
  ord <- order(x, -delta)
  x <- x[ord]
  delta <- delta[ord]
  fit_max <- min(max(t), max(x))
  if (fit_max <= 0 || !any(delta == 1)) {
    return(list(
      cumhaz = rep(0, length(t)), fit_grid = c(0, fit_max),
      fit_hazard = c(0, 0), fit_cumhaz = c(0, 0), fit_max = fit_max,
      extension_used = any(t > fit_max), no_events = !any(delta == 1),
      warnings = character()
    ))
  }
  if (!requireNamespace("muhaz", quietly = TRUE)) {
    stop("kernel_curve() requires the R package 'muhaz'. Install it with ",
      "install.packages('muhaz').",
      call. = FALSE
    )
  }
  warnings <- character()
  f <- withCallingHandlers(
    muhaz::muhaz(
      x, delta,
      min.time = 0, max.time = fit_max,
      bw.method = "local", b.cor = "both", kern = "epanechnikov",
      n.min.grid = 51L, n.est.grid = as.integer(n_grid)
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (any(!is.finite(f$haz.est))) {
    stop("Nonfinite muhaz hazard. No sample is discarded or retried.",
      call. = FALSE
    )
  }
  grid <- f$est.grid
  h <- cumtrap(grid, f$haz.est)
  ans <- approx(grid, h, xout = t, rule = 2, ties = "ordered")$y
  list(
    cumhaz = ans, fit_grid = grid, fit_hazard = f$haz.est,
    fit_cumhaz = h, fit_max = tail(grid, 1),
    extension_used = any(t > tail(grid, 1)), no_events = FALSE,
    warnings = warnings
  )
}
