#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output <- "perf-metrics.csv"

if (length(args) > 0L) {
  for (i in seq_along(args)) {
    if (identical(args[[i]], "--output") && i < length(args)) {
      output <- args[[i + 1L]]
    }
  }
}

if (!suppressWarnings(requireNamespace("eegwhiten", quietly = TRUE))) {
  if (!suppressWarnings(requireNamespace("pkgload", quietly = TRUE))) {
    stop("Either installed package 'eegwhiten' or package 'pkgload' is required.")
  }
  pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
} else {
  suppressPackageStartupMessages(library(eegwhiten))
}

bench_one <- function(n, d, k, method, eig_method, fast, reps = 3L) {
  set.seed(1000 + n + d + k + nchar(method))
  X <- matrix(rnorm(n * d), n, d)
  times <- numeric(reps)
  for (i in seq_len(reps)) {
    t <- system.time({
      model <- whiten_model(
        X,
        method = method,
        n_comp = k,
        lambda = 0.05,
        eig_method = eig_method,
        fast = fast
      )
      Z <- predict(model, X)
      invisible(check_whitening(Z))
    })
    times[[i]] <- unname(t[["elapsed"]])
  }
  data.frame(
    n = n,
    d = d,
    k = k,
    method = method,
    eig_method = eig_method,
    fast = fast,
    elapsed_median = median(times),
    elapsed_min = min(times),
    elapsed_max = max(times),
    stringsAsFactors = FALSE
  )
}

cases <- list(
  list(n = 1200, d = 128, k = 32, method = "PCA"),
  list(n = 900, d = 192, k = 48, method = "PCA-cor"),
  list(n = 800, d = 160, k = 40, method = "SVD")
)

rows <- list()
idx <- 1L
for (cfg in cases) {
  rows[[idx]] <- do.call(bench_one, c(cfg, list(eig_method = "base", fast = FALSE)))
  idx <- idx + 1L
  rows[[idx]] <- do.call(bench_one, c(cfg, list(eig_method = "auto", fast = TRUE)))
  idx <- idx + 1L
}

metrics <- do.call(rbind, rows)
row.names(metrics) <- NULL

if (!dir.exists(dirname(output))) {
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
}
write.csv(metrics, output, row.names = FALSE)
print(metrics)

# Soft regression guard: auto+fast should not be dramatically slower than base.
merged <- merge(
  subset(metrics, eig_method == "base"),
  subset(metrics, eig_method == "auto"),
  by = c("n", "d", "k", "method"),
  suffixes = c("_base", "_auto")
)
slow_ratio <- merged$elapsed_median_auto / pmax(merged$elapsed_median_base, .Machine$double.eps)
if (any(slow_ratio > 1.35)) {
  stop("Performance regression guard failed: auto+fast median time is >35% slower than base for at least one benchmark case.")
}
