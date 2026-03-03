#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output <- "tests/test-report.md"

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

suppressPackageStartupMessages(library(testthat))

results <- testthat::test_dir(
  "tests/testthat",
  reporter = testthat::SilentReporter$new()
)

expectations <- unlist(lapply(results, function(x) x$results), recursive = FALSE)

classify_expectation <- function(x) {
  if (inherits(x, "expectation_failure")) return("failure")
  if (inherits(x, "expectation_error")) return("error")
  if (inherits(x, "expectation_warning")) return("warning")
  if (inherits(x, "expectation_skip")) return("skip")
  if (inherits(x, "expectation_success")) return("success")
  "unknown"
}

status <- vapply(expectations, classify_expectation, character(1))
status_levels <- c("success", "failure", "error", "warning", "skip", "unknown")
status_counts <- setNames(integer(length(status_levels)), status_levels)
status_tab <- table(status)
status_counts[names(status_tab)] <- as.integer(status_tab)

total_tests <- length(results)
total_expectations <- length(expectations)

failed_rows <- list()
if (total_expectations > 0L) {
  failed_idx <- which(status %in% c("failure", "error", "warning"))
  if (length(failed_idx) > 0L) {
    k <- 1L
    for (res in results) {
      for (exp in res$results) {
        s <- classify_expectation(exp)
        if (s %in% c("failure", "error", "warning")) {
          failed_rows[[k]] <- data.frame(
            file = res$file,
            test = res$test,
            type = s,
            message = gsub("[\r\n]+", " ", exp$message),
            stringsAsFactors = FALSE
          )
          k <- k + 1L
        }
      }
    }
  }
}

report_lines <- c(
  "# Test Report",
  "",
  sprintf("- Generated at: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("- Total test blocks: %d", total_tests),
  sprintf("- Total expectations: %d", total_expectations),
  "",
  "## Summary",
  "",
  "| Status | Count |",
  "|---|---:|",
  sprintf("| success | %d |", status_counts[["success"]]),
  sprintf("| failure | %d |", status_counts[["failure"]]),
  sprintf("| error | %d |", status_counts[["error"]]),
  sprintf("| warning | %d |", status_counts[["warning"]]),
  sprintf("| skip | %d |", status_counts[["skip"]]),
  sprintf("| unknown | %d |", status_counts[["unknown"]]),
  ""
)

if (length(failed_rows) == 0L) {
  report_lines <- c(report_lines, "## Non-success Details", "", "No failures, errors, or warnings.")
} else {
  failed_df <- do.call(rbind, failed_rows)
  report_lines <- c(
    report_lines,
    "## Non-success Details",
    "",
    "| File | Test | Type | Message |",
    "|---|---|---|---|"
  )
  for (i in seq_len(nrow(failed_df))) {
    report_lines <- c(
      report_lines,
      sprintf(
        "| %s | %s | %s | %s |",
        failed_df$file[[i]],
        failed_df$test[[i]],
        failed_df$type[[i]],
        failed_df$message[[i]]
      )
    )
  }
}

out_dir <- dirname(output)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}
writeLines(report_lines, con = output, useBytes = TRUE)

cat(sprintf("Wrote test report to %s\n", output))

