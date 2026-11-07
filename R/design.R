# Simulation truth, unconditional sampling, and censoring calibration.
# Source validation.R before this file. These functions need only base R.

true_chf <- function(t, scenario) {
  .assert_numeric_vector(t, "t", nonnegative = TRUE)
  .validate_scenario(scenario)
  switch(scenario,
    S1 = t,
    S2 = (t / 2)^2,
    S3 = log1p((t / 1.5)^3),
    S4 = 0.25 * pmin(t, 2) + 1.5 * pmax(t - 2, 0)
  )
}

# Draw event times, then independent exponential censoring times. This function
# uses the caller's RNG stream. simulate_cell() supplies its reproducible stream.
generate_sample <- function(n, scenario, gamma) {
  .assert_scalar(n, "n",
    lower = 1, integer = TRUE,
    upper = .Machine$integer.max
  )
  .validate_scenario(scenario)
  .assert_scalar(gamma, "gamma", lower = 0, strictly = TRUE)
  event <- switch(scenario,
    S1 = rexp(n),
    S2 = rweibull(n, shape = 2, scale = 2),
    S3 = {
      u <- runif(n)
      1.5 * (u / (1 - u))^(1 / 3)
    },
    S4 = {
      e <- rexp(n)
      ifelse(e <= 0.5, e / 0.25, 2 + (e - 0.5) / 1.5)
    }
  )
  censor <- rexp(n, rate = gamma)
  list(x = pmin(event, censor), delta = as.integer(event <= censor))
}

# P(C < T) = integral gamma exp(-gamma t) S_T(t) dt.
population_censor <- function(scenario, gamma) {
  .validate_scenario(scenario)
  .assert_scalar(gamma, "gamma", lower = 0, strictly = TRUE)
  if (scenario == "S1") {
    return(gamma / (1 + gamma))
  }
  integrate(
    function(t) gamma * exp(-gamma * t - true_chf(t, scenario)),
    lower = 0, upper = Inf, rel.tol = 1e-10
  )$value
}

# The 18 main design cells and two supplementary abrupt-hazard cells.
# Cell ordering and seeds are part of the reference simulation protocol.
make_cells <- function(B = 500L, tau = 3.7) {
  .assert_scalar(B, "B",
    lower = 2, integer = TRUE,
    upper = .Machine$integer.max
  )
  .assert_scalar(tau, "tau", lower = 0, strictly = TRUE)
  targets <- c(0.1, 0.3, 0.5)
  rates <- lapply(c("S1", "S2", "S3"), function(s) {
    vapply(targets, function(p) {
      if (s == "S1") {
        return(p / (1 - p))
      }
      uniroot(function(g) population_censor(s, g) - p,
        c(0.00001, 5),
        tol = 1e-11
      )$root
    }, 0)
  })
  names(rates) <- c("S1", "S2", "S3")
  out <- list()
  id <- 0L
  for (s in names(rates)) {
    for (n in c(100L, 500L)) {
      for (j in seq_len(3L)) {
        id <- id + 1L
        out[[id]] <- data.frame(
          cell_id = sprintf("%s_n%d_c%02d", s, n, c(10, 30, 50)[j]),
          scenario = s, n = n, target_censor = targets[j],
          gamma = rates[[s]][j], B = B, tau = tau,
          seed = 251893938L + id * 1000L, supplement = FALSE
        )
      }
    }
  }
  gamma4 <- uniroot(function(g) population_censor("S4", g) - 0.3,
    c(0.0001, 5),
    tol = 1e-11
  )$root
  for (n in c(100L, 500L)) {
    id <- id + 1L
    out[[id]] <- data.frame(
      cell_id = sprintf("S4_n%d_c30", n), scenario = "S4", n = n,
      target_censor = 0.3, gamma = gamma4, B = B, tau = tau,
      seed = 251893938L + id * 1000L, supplement = TRUE
    )
  }
  ans <- do.call(rbind, out)
  ans$population_censor <- mapply(population_censor, ans$scenario, ans$gamma)
  ans
}
