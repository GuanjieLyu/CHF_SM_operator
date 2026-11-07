#!/usr/bin/env Rscript
# A small deterministic example that does not require the kernel dependency.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE)))
source(file.path(root, "R", "bootstrap.R"), encoding = "UTF-8")
load_project(root)

fit <- grouped_na(x = c(0.5, 1, 1, 2, 3), delta = c(1, 1, 0, 1, 0))
time <- seq(0, 3, length.out = 7)
print(data.frame(
  time = time,
  Nelson_Aalen = na_curve(time, fit),
  Szasz_Mirakyan = sm_curve(time, fit, m = 2.5)
))
