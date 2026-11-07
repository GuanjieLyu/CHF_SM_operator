# Plot saved simulation summaries without rerunning the experiment.
#
# After loading this file, call:
# write_simulation_figures("results", "figures", "tables")
# The three input CSV files are read only. All figures and LaTeX tables are
# written to the supplied output directories. FigS1.pdf is Figure B.1 in the
# manuscript. The tabulated minimum is over the studied grid, while the
# estimator itself permits every positive real smoothing parameter.

.sim_plot_directory <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("Each output directory must be a nonempty path.", call. = FALSE)
  }
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE)) {
    stop("Cannot create output directory: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.sim_plot_require_columns <- function(data, columns, filename) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop(filename, " is missing columns: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

.sim_plot_pdf <- function(path, height, draw) {
  previous_device <- grDevices::dev.cur()
  grDevices::pdf(path,
    width = 6.3, height = height, family = "Times",
    pointsize = 11, useDingbats = FALSE
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

.sim_plot_legend <- function() {
  graphics::par(cex = 1, mar = rep(0, 4))
  graphics::plot.new()
  graphics::legend(
    "center", c("Nelson-Aalen", "Kernel", "Szasz-Mirakyan"),
    col = c("#333333", "#D55E00", "#0072B2"), lty = c(2, 3, 1),
    lwd = c(1.6, 1.8, 1.8), horiz = TRUE, bty = "n", cex = 1,
    x.intersp = 0.7, seg.len = 2
  )
}

.sim_plot_panel <- function(curves, cell, panel_letter) {
  rows <- curves[curves$cell_id == cell$cell_id, , drop = FALSE]
  sm <- rows[rows$method == "SM", , drop = FALSE]
  sm <- sm[order(sm$m), , drop = FALSE]
  na <- rows[rows$method == "NA", , drop = FALSE]
  kernel <- rows[rows$method == "Kernel", , drop = FALSE]
  if (!nrow(sm) || anyDuplicated(sm$m) || nrow(na) != 1L || nrow(kernel) != 1L) {
    stop("Expected one NA curve, one kernel curve and unique SM grid values for ",
      cell$cell_id,
      call. = FALSE
    )
  }
  if (any(!is.finite(sm$m)) || any(sm$m <= 0)) {
    stop("SM grid values must be finite positive real numbers.", call. = FALSE)
  }
  limits <- range(c(0, rows$IMSE + 1.96 * rows$MCSE))
  if (limits[2] == 0) limits[2] <- 1
  graphics::par(
    cex = 1, mar = c(3.7, 4, 2.3, 0.6), mgp = c(2.3, 0.7, 0),
    tcl = -0.25, cex.axis = 0.95, cex.lab = 1.05, cex.main = 0.98,
    las = 1, bty = "l"
  )
  graphics::plot(
    sm$m, sm$IMSE,
    type = "n", log = "x", xlim = range(sm$m),
    ylim = limits, xaxt = "n", xlab = expression(m ~ "(log scale)"),
    ylab = "IMSE", main = sprintf(
      "(%s) n = %d, %d%% censoring",
      panel_letter, cell$n, round(100 * cell$target_censor)
    )
  )
  ticks <- c(1, 2, 5, 10, 20, 50, 100)
  ticks <- ticks[ticks >= min(sm$m) & ticks <= max(sm$m)]
  if (length(ticks) < 2L) ticks <- range(sm$m)
  graphics::axis(1, at = ticks, labels = ticks, cex.axis = 0.95)
  graphics::polygon(
    c(sm$m, rev(sm$m)),
    c(pmax(0, sm$IMSE - 1.96 * sm$MCSE), rev(sm$IMSE + 1.96 * sm$MCSE)),
    col = grDevices::adjustcolor("#0072B2", alpha.f = 0.15), border = NA
  )
  graphics::abline(h = na$IMSE, col = "#333333", lty = 2, lwd = 1.6)
  graphics::abline(h = kernel$IMSE, col = "#D55E00", lty = 3, lwd = 1.8)
  graphics::lines(sm$m, sm$IMSE, col = "#0072B2", lwd = 1.8)
  best <- which.min(sm$IMSE)
  graphics::points(sm$m[best], sm$IMSE[best],
    pch = 21, bg = "white",
    col = "#0072B2", cex = 1, lwd = 1.4
  )
}

.sim_plot_table <- function(data, path, label, caption) {
  format_error <- function(value, error) sprintf("%.4f (%.4f)", value, error)
  lines <- c(
    "\\begin{table}[!htbp]", "\\centering", paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"), "\\footnotesize",
    "\\setlength{\\tabcolsep}{4pt}", "\\begin{tabular}{@{}lrccllr@{}}", "\\toprule",
    "Setting & $n$ & Censoring & Nelson--Aalen & Kernel & SM grid minimum & $m^\\star_{\\rm grid}$ \\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(data))) {
    row <- data[i, ]
    lines <- c(lines, sprintf(
      "%s & %d & $%d\\%%$ & %s & %s & %s & %g \\\\",
      sub("S", "S.", row$scenario), row$n, round(100 * row$target_censor),
      format_error(row$NA_IMSE, row$NA_MCSE),
      format_error(row$Kernel_IMSE, row$Kernel_MCSE),
      format_error(row$SM_min_IMSE, row$SM_fixed_m_MCSE), row$m_grid_star
    ))
    if (i < nrow(data) && data$scenario[i] != data$scenario[i + 1L]) {
      lines <- c(lines, "\\midrule")
    }
  }
  counts <- unique(data$B)
  replication_text <- if (length(counts) == 1L) {
    sprintf("based on %d replications per setting", counts)
  } else {
    "with replication counts given in the saved summaries"
  }
  lines <- c(
    lines, "\\bottomrule", "\\end{tabular}",
    "\\par\\smallskip\\begin{minipage}{\\textwidth}\\footnotesize",
    paste0(
      "Values are empirical IMSEs with Monte Carlo standard errors in parentheses, ",
      replication_text, ". The SM entry is the minimum over the studied grid. ",
      "Its displayed MCSE is the fixed-$m$ calculation evaluated at the selected ",
      "grid value and does not account for minimization. The kernel bandwidth ",
      "is data-selected. The grid's lower endpoint is a design choice, not a ",
      "lower bound on the parameter domain $m>0$."
    ), "\\end{minipage}", "\\end{table}"
  )
  writeLines(lines, path)
}

.sim_plot_selection_table <- function(summary, frequency, path, n_boot) {
  caption <- sprintf(
    "Monte Carlo stability of the empirical grid minimum under %d paired replicate-row bootstrap resamples.",
    n_boot
  )
  lines <- c(
    "\\begin{table}[!htbp]", "\\centering", paste0("\\caption{", caption, "}"),
    "\\label{tab:selection_stability}", "\\footnotesize",
    "\\setlength{\\tabcolsep}{5pt}", "\\begin{tabular}{@{}lrrrrrr@{}}", "\\toprule",
    "Setting & $n$ & Censoring & $m^\\star_{\\rm grid}$ & Reselection freq. & Bootstrap SD & Central range \\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(summary))) {
    row <- summary[i, ]
    selected <- frequency$selection_frequency[
      frequency$cell_id == row$cell_id & frequency$m == row$m_grid_star
    ]
    if (length(selected) != 1L || !is.finite(selected)) {
      stop("Missing or duplicate selection frequency for ", row$cell_id, call. = FALSE)
    }
    lines <- c(lines, sprintf(
      "%s & %d & $%d\\%%$ & %g & %.3f & %.4f & [%.4f, %.4f] \\\\",
      sub("S", "S.", row$scenario), row$n, round(100 * row$target_censor),
      row$m_grid_star, selected, row$SM_reselected_bootstrap_MCSE,
      row$SM_reselected_bootstrap_lo, row$SM_reselected_bootstrap_hi
    ))
  }
  lines <- c(
    lines, "\\bottomrule", "\\end{tabular}",
    "\\par\\smallskip\\begin{minipage}{\\textwidth}\\footnotesize",
    paste(
      "Reselection frequency is the fraction of bootstrap resamples selecting",
      "the grid value chosen in the original Monte Carlo sample. The bootstrap",
      "SD and central 95\\% range describe the reselected minimum IMSE. This",
      "range is not a confidence interval for the population oracle risk, and",
      "resampling does not remove minimization bias."
    ), "\\end{minipage}", "\\end{table}"
  )
  writeLines(lines, path)
}

#' Write simulation figures and tables from saved CSV summaries.
#'
#' @param results_dir Directory containing curve_summary.csv, summary_table.csv
#'   and bootstrap_minimum.csv from the common simulation engine.
#' @param figure_dir Output directory for Fig2.pdf, Fig3.pdf, Fig4.pdf and
#'   FigS1.pdf (the manuscript's Figure B.1).
#' @param table_dir Output directory for table1.tex, table_stress.tex and
#'   table_selection.tex. No manuscript directory is discovered or modified.
#' @param n_boot Number of bootstrap resamples used to generate the input
#'   summaries. It is used only in the selection-table caption.
#' @return Invisibly, a list of the generated figure and table paths.
write_simulation_figures <- function(results_dir, figure_dir, table_dir,
                                     n_boot = 1000L) {
  if (length(n_boot) != 1L || !is.numeric(n_boot) || !is.finite(n_boot) ||
    n_boot < 1L || n_boot != floor(n_boot)) {
    stop("n_boot must be a positive integer.", call. = FALSE)
  }
  input_files <- file.path(results_dir, c(
    "curve_summary.csv", "summary_table.csv", "bootstrap_minimum.csv"
  ))
  if (!all(file.exists(input_files))) {
    stop("Missing simulation summaries: ",
      paste(input_files[!file.exists(input_files)], collapse = ", "),
      call. = FALSE
    )
  }
  # Preserve the literal method label "NA". Only the numerical m column uses
  # the CSV token NA to denote an inapplicable parameter for the comparators.
  curves <- utils::read.csv(input_files[1], na.strings = "", stringsAsFactors = FALSE)
  summary <- utils::read.csv(input_files[2], stringsAsFactors = FALSE)
  frequency <- utils::read.csv(input_files[3], stringsAsFactors = FALSE)
  .sim_plot_require_columns(curves, c("cell_id", "method", "m", "IMSE", "MCSE"), input_files[1])
  .sim_plot_require_columns(summary, c(
    "cell_id", "scenario", "n", "target_censor", "B", "tau", "supplement",
    "NA_IMSE", "NA_MCSE", "Kernel_IMSE", "Kernel_MCSE", "m_grid_star",
    "SM_min_IMSE", "SM_fixed_m_MCSE", "SM_reselected_bootstrap_MCSE",
    "SM_reselected_bootstrap_lo", "SM_reselected_bootstrap_hi"
  ), input_files[2])
  .sim_plot_require_columns(frequency, c("cell_id", "m", "selection_frequency"), input_files[3])
  curves$m <- utils::type.convert(as.character(curves$m),
    na.strings = c("", "NA"), as.is = TRUE
  )
  if (!all(curves$method %in% c("NA", "Kernel", "SM")) ||
    any(!is.finite(curves$IMSE)) || any(!is.finite(curves$MCSE)) ||
    any(curves$IMSE < 0) || any(curves$MCSE < 0)) {
    stop("Curve summaries contain unknown methods or invalid errors.", call. = FALSE)
  }
  if (anyDuplicated(summary$cell_id) || anyNA(summary$supplement)) {
    stop("Summary rows must have unique cell identifiers and supplement flags.", call. = FALSE)
  }
  summary <- summary[order(summary$scenario, summary$n, summary$target_censor), , drop = FALSE]
  figure_dir <- .sim_plot_directory(figure_dir)
  table_dir <- .sim_plot_directory(table_dir)
  figure_paths <- character()
  for (scenario in c("S1", "S2", "S3", "S4")) {
    cells <- summary[summary$scenario == scenario, , drop = FALSE]
    if (!nrow(cells)) next
    supplement <- scenario == "S4"
    path <- file.path(figure_dir, if (supplement) {
      "FigS1.pdf"
    } else {
      paste0("Fig", as.integer(sub("S", "", scenario)) + 1L, ".pdf")
    })
    .sim_plot_pdf(path, if (supplement) 3.3 else 7.6, function() {
      if (supplement) {
        graphics::layout(rbind(c(1, 1), c(2, 3)), heights = c(0.45, 2.85))
      } else {
        graphics::layout(rbind(c(1, 1), c(2, 3), c(4, 5), c(6, 7)),
          heights = c(0.45, 2.38, 2.38, 2.38)
        )
      }
      .sim_plot_legend()
      panel <- 0L
      censoring <- if (supplement) 0.3 else c(0.1, 0.3, 0.5)
      for (target in censoring) {
        for (sample_size in c(100, 500)) {
          cell <- cells[cells$n == sample_size &
            abs(cells$target_censor - target) < 1e-8, , drop = FALSE]
          if (nrow(cell) != 1L) {
            stop("Incomplete figure design for ", scenario, ", n=", sample_size,
              ", censoring=", target,
              call. = FALSE
            )
          }
          panel <- panel + 1L
          .sim_plot_panel(curves, cell, letters[panel])
        }
      }
    })
    figure_paths <- c(figure_paths, path)
  }
  table_paths <- character()
  for (supplement in c(FALSE, TRUE)) {
    rows <- summary[summary$supplement == supplement, , drop = FALSE]
    if (!nrow(rows)) next
    path <- file.path(table_dir, if (supplement) "table_stress.tex" else "table1.tex")
    caption <- if (supplement) {
      paste(
        "Abrupt-hazard stress experiment S.4 with 30\\% target censoring.",
        "The hazard is bounded but discontinuous, so only the bounded-hazard",
        "case among the stated weak-convergence results applies at its",
        "specified smoothing rate. The local consistency assumptions also apply."
      )
    } else {
      sprintf(paste(
        "Recomputed empirical IMSEs on $[0,%g]$ for the three main scenarios.",
        "The Sz\\'asz--Mirakyan (SM) column is an oracle-grid benchmark,",
        "not an automatic data-driven procedure."
      ), unique(rows$tau))
    }
    .sim_plot_table(rows, path, if (supplement) "tab:stress" else "tab:imse_summary", caption)
    table_paths <- c(table_paths, path)
  }
  selection_path <- file.path(table_dir, "table_selection.tex")
  .sim_plot_selection_table(summary, frequency, selection_path, n_boot)
  invisible(list(figures = figure_paths, tables = c(table_paths, selection_path)))
}
