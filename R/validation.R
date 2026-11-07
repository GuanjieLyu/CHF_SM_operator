# Internal validation and scoped random-number utilities.
# Sourcing this file does not change options, library paths, or RNG state.

.assert_numeric_vector <- function(x, name, min_length = 0L,
                                   nonnegative = FALSE, positive = FALSE) {
  if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) ||
    length(x) < min_length || any(!is.finite(x))) {
    stop(name, " must be a finite numeric vector of length at least ",
      min_length, ".",
      call. = FALSE
    )
  }
  if (nonnegative && any(x < 0)) {
    stop(name, " must contain nonnegative values.", call. = FALSE)
  }
  if (positive && any(x <= 0)) {
    stop(name, " must contain positive values.", call. = FALSE)
  }
  invisible(x)
}

.assert_scalar <- function(x, name, lower = NULL, strictly = FALSE,
                           integer = FALSE, upper = NULL) {
  .assert_numeric_vector(x, name, min_length = 1L)
  if (length(x) != 1L) {
    stop(name, " must be a single finite number.", call. = FALSE)
  }
  if (!is.null(lower) && (if (strictly) x <= lower else x < lower)) {
    relation <- if (strictly) "greater than " else "at least "
    stop(name, " must be ", relation, lower, ".", call. = FALSE)
  }
  if (!is.null(upper) && x > upper) {
    stop(name, " must be at most ", upper, ".", call. = FALSE)
  }
  if (integer && x != floor(x)) {
    stop(name, " must be an integer.", call. = FALSE)
  }
  invisible(x)
}

.validate_observations <- function(x, delta) {
  .assert_numeric_vector(x, "x", min_length = 1L, nonnegative = TRUE)
  if (!(is.numeric(delta) || is.logical(delta)) || is.complex(delta) ||
    !is.null(dim(delta)) || length(delta) != length(x) ||
    anyNA(delta) || any(!delta %in% c(0, 1))) {
    stop("delta must have the same length as x and contain only 0 and 1.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

.validate_grid <- function(grid, name = "grid", min_length = 2L) {
  .assert_numeric_vector(grid, name, min_length, nonnegative = TRUE)
  if (any(diff(grid) <= 0)) {
    stop(name, " must be strictly increasing.", call. = FALSE)
  }
  invisible(grid)
}

.validate_m_grid <- function(m_grid) {
  .assert_numeric_vector(m_grid, "m_grid", min_length = 1L, positive = TRUE)
  if (anyDuplicated(m_grid)) {
    stop("m_grid must contain distinct positive values.", call. = FALSE)
  }
  invisible(m_grid)
}

.validate_na_fit <- function(fit) {
  fields <- c("time", "weight", "cumhaz")
  if (!is.list(fit) || !all(fields %in% names(fit))) {
    stop("fit must be a grouped_na result with time, weight, and cumhaz.",
      call. = FALSE
    )
  }
  for (field in fields) {
    .assert_numeric_vector(fit[[field]], paste0("fit$", field),
      nonnegative = TRUE
    )
  }
  if (length(fit$weight) != length(fit$time) ||
    length(fit$cumhaz) != length(fit$time) ||
    any(diff(fit$time) <= 0) || any(diff(fit$cumhaz) < 0)) {
    stop("fit must have matching lengths and increasing failure times.",
      call. = FALSE
    )
  }
  if (!isTRUE(all.equal(fit$cumhaz, cumsum(fit$weight),
    check.attributes = FALSE
  ))) {
    stop("fit$cumhaz must equal the cumulative sum of fit$weight.",
      call. = FALSE
    )
  }
  invisible(fit)
}

.validate_scenario <- function(scenario) {
  if (!is.character(scenario) || length(scenario) != 1L ||
    is.na(scenario) || !scenario %in% c("S1", "S2", "S3", "S4")) {
    stop("scenario must be one of S1, S2, S3, or S4.", call. = FALSE)
  }
  invisible(scenario)
}

.validate_cell <- function(cell) {
  required <- c(
    "cell_id", "scenario", "n", "target_censor", "gamma",
    "B", "tau", "seed", "supplement"
  )
  if (!is.data.frame(cell) || nrow(cell) != 1L ||
    !all(required %in% names(cell))) {
    stop("cell must be a single row returned by make_cells().", call. = FALSE)
  }
  if (!is.character(cell$cell_id) || is.na(cell$cell_id) ||
    !grepl("^[A-Za-z0-9_-]+$", cell$cell_id)) {
    stop("cell$cell_id must contain only letters, digits, underscores, or hyphens.",
      call. = FALSE
    )
  }
  .validate_scenario(cell$scenario)
  .assert_scalar(cell$n, "cell$n",
    lower = 1, integer = TRUE,
    upper = .Machine$integer.max
  )
  .assert_scalar(cell$B, "cell$B",
    lower = 2, integer = TRUE,
    upper = .Machine$integer.max
  )
  .assert_scalar(cell$tau, "cell$tau", lower = 0, strictly = TRUE)
  .assert_scalar(cell$gamma, "cell$gamma", lower = 0, strictly = TRUE)
  .assert_scalar(cell$target_censor, "cell$target_censor", lower = 0, upper = 1)
  .assert_scalar(cell$seed, "cell$seed",
    lower = 0, integer = TRUE,
    upper = .Machine$integer.max - 500L
  )
  if (!is.logical(cell$supplement) || is.na(cell$supplement)) {
    stop("cell$supplement must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(cell)
}

.with_reproducible_rng <- function(seed, code) {
  .assert_scalar(seed, "seed",
    lower = 0, integer = TRUE,
    upper = .Machine$integer.max
  )
  previous_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    previous_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
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
  force(code)
}
