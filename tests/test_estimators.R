# Run from the repository root: Rscript --vanilla tests/test_estimators.R
# No optional package is needed for these mathematical checks.
local({
  root <- getwd()
  if (!file.exists(file.path(root, "R", "estimators.R"))) {
    script_arg <- grep("^--file=", commandArgs(), value = TRUE)
    if (length(script_arg) != 1L) {
      stop("Run this test from the repository root.")
    }
    root <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
  }
  test_env <- new.env(parent = globalenv())
  previous_libs <- .libPaths()
  previous_kind <- RNGkind()
  previous_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  previous_namespaces <- loadedNamespaces()
  for (file in c("validation.R", "estimators.R", "design.R", "simulation.R")) {
    sys.source(file.path(root, "R", file), envir = test_env)
  }
  stopifnot(identical(.libPaths(), previous_libs), identical(RNGkind(), previous_kind))
  current_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  stopifnot(identical(current_seed, previous_seed))
  stopifnot(identical(loadedNamespaces(), previous_namespaces))

  with(test_env, {
    assert_close <- function(actual, expected, tolerance = 1e-12) {
      stopifnot(isTRUE(all.equal(actual, expected,
        tolerance = tolerance,
        check.attributes = FALSE
      )))
    }
    expect_error <- function(expr, text) {
      error <- tryCatch(
        {
          force(expr)
          NULL
        },
        error = identity
      )
      stopifnot(inherits(error, "error"), grepl(text, conditionMessage(error), fixed = TRUE))
    }

    # Tied censoring stays in the risk set at the failure time.
    fit <- grouped_na(c(1, 1, 2, 3), c(1, 0, 1, 1))
    assert_close(fit$time, c(1, 2, 3))
    assert_close(fit$risk, c(4, 2, 1))
    assert_close(fit$weight, c(1 / 4, 1 / 2, 1))
    assert_close(
      na_curve(c(0, 0.9, 1, 1.5, 2, 3, 4), fit),
      c(0, 0, 0.25, 0.25, 0.75, 1.75, 1.75)
    )
    assert_close(grouped_na(rep(1, 5), rep(1, 5))$weight, 1)
    stopifnot(sum(fit$weight) <= sum(1 / seq_len(fit$n)))
    untied <- grouped_na(seq_len(8), rep(1, 8))
    assert_close(sum(untied$weight), sum(1 / seq_len(8)))

    # Integer and noninteger m*Z give the same finite Poisson/gamma identity.
    t <- c(0, 0.2, 1, 2, 4)
    for (m in c(0.37, 1, 2.5, 4)) {
      shape <- ceiling(m * fit$time)
      poisson <- as.vector(outer(t, shape, function(tt, r) {
        ppois(r - 1, lambda = m * tt, lower.tail = FALSE)
      }) %*% fit$weight)
      assert_close(sm_curve(t, fit, m), poisson)
      stopifnot(all(diff(sm_curve(t, fit, m)) >= 0))
      stopifnot(all(sm_curve(t, fit, m) <= sum(fit$weight)))
      positive_t <- c(0.3, 1.3, 2.7)
      step <- 1e-5
      derivative <- (sm_curve(positive_t + step, fit, m) -
        sm_curve(positive_t - step, fit, m)) / (2 * step)
      gamma_density <- as.vector(outer(positive_t, shape, function(tt, r) {
        dgamma(tt, shape = r, rate = m)
      }) %*% fit$weight)
      assert_close(derivative, gamma_density, tolerance = 1e-8)
    }

    # Changing units t -> t/a and m -> a*m preserves the entire curve.
    scale <- 7.3
    scaled <- grouped_na(c(1, 1, 2, 3) / scale, c(1, 0, 1, 1))
    assert_close(
      sm_curve(t / scale, scaled, scale * 0.37),
      sm_curve(t, fit, 0.37)
    )
    assert_close(na_curve(t / scale, scaled), na_curve(t, fit))

    # No-event estimates do not need muhaz and remain zero after follow-up ends.
    empty <- grouped_na(c(0, 1, 2), c(0, 0, 0))
    assert_close(na_curve(t, empty), rep(0, length(t)))
    assert_close(sm_curve(t, empty, 0.2), rep(0, length(t)))
    kernel <- kernel_curve(c(0, 1, 2), c(0, 0, 0), t)
    assert_close(kernel$cumhaz, rep(0, length(t)))
    stopifnot(kernel$no_events, kernel$extension_used, kernel$fit_max == 2)

    # Zero-time censoring is allowed, and an integral evaluated only at zero
    # is zero even when the sample includes strictly positive failure times.
    origin <- kernel_curve(c(0, 1, 2), c(0, 1, 0), 0)
    assert_close(origin$cumhaz, 0)
    stopifnot(!origin$no_events, !origin$extension_used)
    expect_error(kernel_curve(0, 1, c(0, 1)), "strictly positive failure times")
    expect_error(kernel_curve(c(0, 1), c(1, 0), c(0, 1)),
      "strictly positive failure times"
    )

    # Integration uses actual, potentially nonuniform coordinates.
    grid <- c(0, 0.1, 0.7, 1.5, 3)
    assert_close(cumtrap(grid, 2 * grid + 1), grid^2 + grid)
    assert_close(sum(trap_weights(grid) * (2 * grid + 1)), 12)
    assert_close(trap_weights(c(0, 2)), c(1, 1))

    # The truncated Poisson basis stays within its uniform tail bound.
    grid <- seq(0, 3.7, length.out = 51L)
    bases <- poisson_bases(grid, c(0.37, 2.5), tail_tol = 1e-8)
    for (base in bases) {
      approximation <- as.vector(base$P %*% na_curve(base$j / base$m, fit))
      exact <- sm_curve(grid, fit, base$m)
      bound <- sum(fit$weight) * base$omitted_mass
      stopifnot(max(abs(exact - approximation)) <= bound + 1e-14)
    }

    expect_error(grouped_na(c(1, NA), c(1, 0)), "finite numeric")
    expect_error(grouped_na(c(1, 2), c(1, NA)), "only 0 and 1")
    expect_error(grouped_na(c(1, 2), c(1, 2)), "only 0 and 1")
    expect_error(sm_curve(t, fit, c(1, 2)), "single finite number")
    expect_error(sm_curve(t, fit, Inf), "finite numeric")
    expect_error(sm_curve(t, fit, 0), "greater than 0")
    expect_error(sm_curve(t, grouped_na(0, 1), 1), "strictly positive")
    expect_error(na_curve(c(0, Inf), fit), "finite numeric")
    expect_error(trap_weights(c(0, 1, 1)), "strictly increasing")
    expect_error(cumtrap(c(0, 1), 1), "length at least 2")
    expect_error(kernel_curve(1, 1, t, n_grid = 3.5), "integer")
    expect_error(poisson_bases(grid, c(1, 1)), "distinct")
    expect_error(poisson_bases(grid, 1, tail_tol = 1e-30), "large enough")
  })
  cat("Estimator tests passed.\n")
})
