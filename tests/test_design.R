# Run from the repository root: Rscript --vanilla tests/test_design.R
# No optional package is needed for design, RNG, and summary checks.
local({
  root <- getwd()
  if (!file.exists(file.path(root, "R", "design.R"))) {
    script_arg <- grep("^--file=", commandArgs(), value = TRUE)
    if (length(script_arg) != 1L) {
      stop("Run this test from the repository root.")
    }
    root <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
  }
  test_env <- new.env(parent = globalenv())
  for (file in c("validation.R", "estimators.R", "design.R", "simulation.R")) {
    sys.source(file.path(root, "R", file), envir = test_env)
  }
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

    cells <- make_cells()
    stopifnot(
      nrow(cells) == 20L, sum(cells$supplement) == 2L,
      all(cells$B == 500L), all(cells$tau == 3.7),
      all(cells$seed == 251893938L + seq_len(20L) * 1000L)
    )
    stopifnot(identical(
      cells$cell_id[c(1, 18, 19, 20)],
      c("S1_n100_c10", "S3_n500_c50", "S4_n100_c30", "S4_n500_c30")
    ))
    stopifnot(max(abs(cells$population_censor - cells$target_censor)) < 1e-8)
    assert_close(population_censor("S1", 1), 0.5)
    assert_close(true_chf(c(0, 1, 2, 3), "S2"), c(0, 0.25, 1, 2.25))
    assert_close(true_chf(c(0, 1.5), "S3"), c(0, log(2)))
    assert_close(true_chf(c(0, 1, 2, 3), "S4"), c(0, 0.25, 0.5, 2))
    for (scenario in c("S1", "S2", "S3", "S4")) {
      stopifnot(population_censor(scenario, 0.1) < population_censor(scenario, 0.5))
    }

    # Preserve draw order: event draws precede independent censor draws.
    for (scenario in c("S1", "S2", "S3", "S4")) {
      sample <- .with_reproducible_rng(123L, generate_sample(12L, scenario, 0.4))
      reference <- .with_reproducible_rng(123L, {
        event <- switch(scenario,
          S1 = rexp(12L),
          S2 = rweibull(12L, shape = 2, scale = 2),
          S3 = {
            u <- runif(12L)
            1.5 * (u / (1 - u))^(1 / 3)
          },
          S4 = {
            e <- rexp(12L)
            ifelse(e <= 0.5, e / 0.25, 2 + (e - 0.5) / 1.5)
          }
        )
        censor <- rexp(12L, rate = 0.4)
        list(x = pmin(event, censor), delta = as.integer(event <= censor))
      })
      stopifnot(identical(sample, reference))
    }

    # Tiny hand-checkable losses verify paired MCSE and bootstrap selection.
    cell <- cells[1, , drop = FALSE]
    cell$B <- 4L
    loss <- cbind(
      "NA" = c(1, 2, 3, 4), Kernel = c(2, 2, 4, 4),
      SM_m1 = c(0.5, 1.5, 2.5, 3.5), SM_m2 = c(3, 1, 3, 1)
    )
    result <- list(
      cell = cell, m_grid = c(1, 2), loss = loss,
      diagnostics = data.frame(
        observed_censor = rep(0.1, 4), kernel_tail_extension = rep(FALSE, 4),
        no_events = rep(FALSE, 4), kernel_warning = rep("", 4),
        kernel_min_hazard = rep(0, 4)
      ),
      max_sm_uniform_truncation_bound = 1e-12
    )
    summary <- summarize_cell(result, n_boot = 100L)
    assert_close(summary$summary$NA_IMSE, 2.5)
    assert_close(summary$summary$NA_MCSE, sd(loss[, "NA"]) / 2)
    stopifnot(summary$summary$m_grid_star == 1)
    assert_close(summary$summary$SM_minus_NA, -0.5)
    assert_close(summary$summary$SM_minus_NA_MCSE, 0)
    assert_close(sum(summary$bootstrap$selection_frequency), 1)
    direct_minimum <- .with_reproducible_rng(cell$seed + 500L, {
      counts <- rmultinom(100L, 4L, rep(0.25, 4))
      apply(t(loss[, c("SM_m1", "SM_m2")]) %*% (counts / 4), 2, min)
    })
    assert_close(summary$boot_minimum_values, direct_minimum)

    # Seeded helpers restore state on success and errors, including absent seed.
    .with_reproducible_rng(77L, {
      RNGkind("L'Ecuyer-CMRG", "Box-Muller", "Rejection")
      set.seed(812L)
      saved_kind <- RNGkind()
      saved_seed <- get(".Random.seed", envir = .GlobalEnv)
      other <- summarize_cell(result, n_boot = 100L)
      stopifnot(
        identical(summary, other), identical(RNGkind(), saved_kind),
        identical(get(".Random.seed", envir = .GlobalEnv), saved_seed)
      )
      expect_error(
        .with_reproducible_rng(4L, stop("deliberate failure")),
        "deliberate failure"
      )
      stopifnot(
        identical(RNGkind(), saved_kind),
        identical(get(".Random.seed", envir = .GlobalEnv), saved_seed)
      )
      rm(".Random.seed", envir = .GlobalEnv)
      .with_reproducible_rng(5L, runif(3))
      stopifnot(
        !exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
        identical(RNGkind(), saved_kind)
      )
    })

    expect_error(make_cells(B = 1), "at least 2")
    expect_error(generate_sample(2.5, "S1", 1), "integer")
    expect_error(population_censor("unknown", 1), "one of S1")
    expect_error(population_censor("S1", 0), "greater than 0")
    expect_error(true_chf(c(0, Inf), "S1"), "finite numeric")
    expect_error(summarize_cell(result, n_boot = 1), "at least 2")
    broken <- result
    broken$loss[1, 1] <- NA_real_
    expect_error(summarize_cell(broken), "finite, nonnegative")
  })
  cat("Design and summary tests passed.\n")
})
