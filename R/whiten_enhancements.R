############################################################
# Advanced whitening workflows: update, tuning, and reporting
############################################################

#' Incrementally Update a Whitening Model
#'
#' Updates a fitted \code{whiten_model} using only summary moments from
#' new observations, avoiding a full refit from all historical data.
#' Currently supports models fitted with
#' \code{cov_estimator = "empirical"}.
#'
#' @param model A \code{whiten_model} object.
#' @param X_new Numeric matrix of new samples with
#'   \code{ncol(X_new) = model$dim_in}.
#' @param sample_weight Optional non-negative sample weights for
#'   \code{X_new}.
#'
#' @return Updated \code{whiten_model} object.
#'
#' @seealso \code{\link{whiten_model}}
#'
#' @examples
#' set.seed(1)
#' X1 <- matrix(rnorm(200 * 8), 200, 8)
#' m <- whiten_model(X1, method = "ZCA", lambda = 0.1)
#' X2 <- matrix(rnorm(100 * 8), 100, 8)
#' m_updated <- whiten_model_update(m, X2)
#' m_updated$n_obs
#'
#' @export
whiten_model_update <- function(model, X_new, sample_weight = NULL) {
  if (!inherits(model, "whiten_model")) {
    stop("model must be a 'whiten_model' object.")
  }
  if (!is.matrix(X_new) || !is.numeric(X_new) || any(!is.finite(X_new))) {
    stop("X_new must be a finite numeric matrix.")
  }

  d <- if (!is.null(model$dim_in)) model$dim_in else ncol(model$W)
  if (ncol(X_new) != d) {
    stop(sprintf("Mismatch: model expects %d columns, but X_new has %d.", d, ncol(X_new)))
  }

  cov_estimator <- if (!is.null(model$cov_estimator)) model$cov_estimator else "empirical"
  if (!identical(cov_estimator, "empirical")) {
    stop("whiten_model_update currently supports only cov_estimator = 'empirical'.")
  }

  needed_stats <- c("stat_sum_w", "stat_sum_w2", "stat_sum_x", "stat_sum_xx")
  if (!all(vapply(needed_stats, function(nm) !is.null(model[[nm]]), logical(1)))) {
    stop("model is missing stored moments required for incremental updates; refit with current whiten_model().")
  }

  w_new <- .validate_sample_weight(sample_weight, nrow(X_new), allow_null = TRUE)
  m_new <- .weighted_moments(X_new, sample_weight = w_new)

  sum_w <- model$stat_sum_w + m_new$sum_w
  sum_w2 <- model$stat_sum_w2 + m_new$sum_w2
  sum_x <- model$stat_sum_x + m_new$sum_x
  sum_xx <- model$stat_sum_xx + m_new$sum_xx

  center_enabled <- if (!is.null(model$center_enabled)) {
    isTRUE(model$center_enabled)
  } else {
    TRUE
  }

  cfit <- .cov_from_moments(
    sum_w = sum_w,
    sum_w2 = sum_w2,
    sum_x = sum_x,
    sum_xx = sum_xx,
    center = center_enabled
  )

  mu <- cfit$mean
  names(mu) <- names(model$center)
  S <- cfit$cov

  lambda <- if (!is.null(model$lambda) && is.finite(model$lambda)) model$lambda else 0
  .check_unit_interval(lambda, "lambda")

  if (lambda > 0) {
    target_val <- mean(diag(S))
    target_mat <- diag(target_val, d)
    S <- (1 - lambda) * S + lambda * target_mat
  }

  method <- .normalize_whiten_method(model$method)
  is_reduction <- .is_reduction_method(method)

  eig_method <- if (!is.null(model$eig_method) && !is.na(model$eig_method) &&
    model$eig_method %in% c("auto", "base", "rspectra")) {
    model$eig_method
  } else {
    "auto"
  }
  fast <- if (!is.null(model$fast)) isTRUE(model$fast) else FALSE

  n_comp_final <- if (is_reduction) model$n_comp else NULL
  var_threshold <- if (!is.null(model$var_threshold)) model$var_threshold else NULL

  need_full_eigs <- !is_reduction || !is.null(var_threshold) || is.null(n_comp_final) || identical(eig_method, "base")
  eig_shrunk <- .check_symmetric_pd(S, require_pd = TRUE, return_eigen = need_full_eigs)
  attr(S, ".checked_spd") <- TRUE

  explained_var <- NA_real_
  spectrum_info <- NULL
  if (is_reduction) {
    if (is.null(n_comp_final) && !is.null(var_threshold)) {
      spectrum_info <- if (method %in% c("PCA", "SVD") && need_full_eigs) {
        list(values = pmax(eig_shrunk, 0), total = sum(diag(S)))
      } else {
        .component_spectrum_info(S, method, eig_method = eig_method, fast = fast)
      }
      n_comp_final <- .n_comp_from_threshold(spectrum_info$values, var_threshold)
    }

    if (is.null(n_comp_final)) {
      n_comp_final <- d
    }

    if (is.null(spectrum_info)) {
      k_spec <- min(n_comp_final, d)
      spectrum_info <- if (method %in% c("PCA", "SVD") && need_full_eigs) {
        list(values = pmax(eig_shrunk[seq_len(k_spec)], 0), total = sum(diag(S)))
      } else {
        .component_spectrum_info(
          S,
          method,
          k = k_spec,
          eig_method = eig_method,
          fast = fast
        )
      }
    }

    explained_var <- .explained_ratio(spectrum_info$values, n_comp_final, total = spectrum_info$total)
  }

  fit_args <- list(Sigma = S, returnW = TRUE, PhiPsi = FALSE, return_decomp = TRUE)
  if (is_reduction && !is.null(n_comp_final)) {
    fit_args$n_comp <- n_comp_final
  }
  if (is_reduction && !is.null(model$sign_basis) && is.matrix(model$sign_basis)) {
    fit_args$sign_ref <- model$sign_basis
  }
  if (method %in% c("PCA", "PCA-cor", "SVD", "ZCA", "ZCA-cor")) {
    fit_args$eig_method <- eig_method
    fit_args$fast <- fast
  }

  res <- do.call(.method_to_function(method), fit_args)
  decomp <- if (!is.null(res$decomp)) res$decomp else list()
  W <- res$W

  updated <- model
  updated$W <- W
  updated$Wt <- t(W)
  updated$center <- mu
  updated$n_comp <- n_comp_final
  updated$explained_var <- explained_var
  updated$U <- decomp$U
  updated$D <- decomp$D
  updated$inv_tW <- if (!is.null(decomp$inv_tW)) decomp$inv_tW else .pinv(t(W))
  updated$sign_basis <- if (!is.null(decomp$U)) decomp$U else updated$sign_basis
  updated$cov <- S

  updated$stat_sum_w <- sum_w
  updated$stat_sum_w2 <- sum_w2
  updated$stat_sum_x <- sum_x
  updated$stat_sum_xx <- sum_xx
  updated$n_obs <- if (!is.null(updated$n_obs)) updated$n_obs + nrow(X_new) else nrow(X_new)
  updated$update_count <- if (!is.null(updated$update_count)) updated$update_count + 1L else 1L
  updated$sample_weighted <- isTRUE(updated$sample_weighted) || !is.null(sample_weight)
  updated$dim_in <- d
  updated$dim_out <- nrow(W)

  class(updated) <- "whiten_model"
  updated
}

.normalize_lambda_grid <- function(lambda_grid) {
  vals <- if (is.list(lambda_grid)) lambda_grid else as.list(lambda_grid)
  if (length(vals) == 0L) {
    stop("lambda_grid must contain at least one value.")
  }

  lapply(vals, function(v) {
    if (is.character(v)) {
      if (length(v) != 1L || !(v %in% "auto")) {
        stop("Character values in lambda_grid must be exactly 'auto'.")
      }
      return("auto")
    }
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v)) {
      stop("Numeric values in lambda_grid must be finite scalars in [0, 1].")
    }
    .check_unit_interval(v, "lambda")
    as.numeric(v)
  })
}

.lambda_to_label <- function(x, digits = 4) {
  if (is.character(x)) return(x)
  format(signif(x, digits = digits), trim = TRUE, scientific = FALSE)
}

#' Automatically Tune Whitening Hyperparameters
#'
#' Performs cross-validation across candidate whitening configurations and
#' returns the best model fit on full data.
#'
#' @param X Numeric matrix \code{[n_samples x n_features]}.
#' @param y Optional labels for supervised tuning.
#' @param methods Candidate whitening methods.
#' @param n_comp_grid Optional integer grid for \code{n_comp}.
#' @param var_threshold_grid Optional numeric grid in \code{(0, 1]} for
#'   \code{var_threshold}.
#' @param lambda_grid Candidate shrinkage values. Can include numeric values
#'   in \code{[0, 1]} and \code{"auto"}.
#' @param lambda_method Method used when \code{lambda = "auto"}.
#' @param cov_estimator_grid Candidate covariance estimators.
#' @param sample_weight Optional sample weights (used in training folds).
#' @param cv_folds Number of folds.
#' @param seed Optional random seed.
#' @param scoring Scoring mode: \code{"auto"}, \code{"unsupervised"}, or
#'   \code{"accuracy"}.
#' @param score_fn Optional custom scoring function. For unsupervised tuning,
#'   called as \code{score_fn(Z_valid)}. For supervised tuning, called as
#'   \code{score_fn(Z_train, y_train, Z_valid, y_valid)}.
#' @param tyler_tol Convergence tolerance for Tyler covariance.
#' @param tyler_max_iter Maximum iterations for Tyler covariance.
#' @param tyler_eps Numerical stabilization floor for Tyler covariance.
#' @param eig_method Eigen solver backend.
#' @param fast Logical; fast mode for eigensolver.
#' @param sign_reference Optional reference vectors for sign stabilization.
#' @param top_n Number of top configurations retained in the summary table.
#'
#' @return A list of class \code{whiten_tune} containing the best model,
#'   best parameters, and cross-validation ranking.
#'
#' @seealso \code{\link{whiten_model}}, \code{\link{report_whitening}}
#'
#' @examples
#' set.seed(42)
#' X <- matrix(rnorm(150 * 6), 150, 6)
#' tuned <- auto_tune_whitening(
#'   X,
#'   methods = c("PCA", "ZCA"),
#'   n_comp_grid = c(3),
#'   lambda_grid = list(0, 0.1),
#'   cv_folds = 3, seed = 1, top_n = 4
#' )
#' head(tuned$ranking)
#'
#' @export
auto_tune_whitening <- function(X,
                                y = NULL,
                                methods = c("SVD", "ZCA", "ZCA-cor", "PCA", "PCA-cor", "Cholesky"),
                                n_comp_grid = NULL,
                                var_threshold_grid = NULL,
                                lambda_grid = list("auto", 0, 0.05, 0.1),
                                lambda_method = c("oas", "lw"),
                                cov_estimator_grid = "empirical",
                                sample_weight = NULL,
                                cv_folds = 5L,
                                seed = NULL,
                                scoring = c("auto", "unsupervised", "accuracy"),
                                score_fn = NULL,
                                tyler_tol = 1e-6,
                                tyler_max_iter = 200,
                                tyler_eps = 1e-6,
                                eig_method = c("auto", "base", "rspectra"),
                                fast = FALSE,
                                sign_reference = NULL,
                                top_n = 5L) {
  if (!is.matrix(X) || !is.numeric(X) || any(!is.finite(X))) {
    stop("X must be a finite numeric matrix.")
  }
  n <- nrow(X)
  d <- ncol(X)
  if (n < 4L) stop("X must have at least 4 rows for cross-validation.")
  if (d < 2L) stop("X must have at least 2 columns.")

  lambda_method <- match.arg(lambda_method)
  eig_method <- match.arg(eig_method)
  scoring <- match.arg(scoring)

  if (!is.null(score_fn) && !is.function(score_fn)) {
    stop("score_fn must be NULL or a function.")
  }

  methods <- unique(vapply(methods, .normalize_whiten_method, character(1)))

  if (!is.null(n_comp_grid)) {
    if (!is.numeric(n_comp_grid) || any(!is.finite(n_comp_grid))) {
      stop("n_comp_grid must be numeric.")
    }
    if (any(n_comp_grid != as.integer(n_comp_grid)) || any(n_comp_grid < 1L) || any(n_comp_grid > d)) {
      stop("n_comp_grid values must be integers between 1 and ncol(X).")
    }
    n_comp_grid <- sort(unique(as.integer(n_comp_grid)))
  }

  if (!is.null(var_threshold_grid)) {
    if (!is.numeric(var_threshold_grid) || any(!is.finite(var_threshold_grid))) {
      stop("var_threshold_grid must be numeric.")
    }
    invisible(vapply(var_threshold_grid, .check_var_threshold, numeric(0)))
    var_threshold_grid <- sort(unique(as.numeric(var_threshold_grid)))
  }

  lambda_values <- .normalize_lambda_grid(lambda_grid)

  cov_estimator_grid <- match.arg(
    cov_estimator_grid,
    choices = c("empirical", "mcd", "tyler"),
    several.ok = TRUE
  )
  cov_estimator_grid <- unique(cov_estimator_grid)

  sample_weight <- .validate_sample_weight(sample_weight, n, allow_null = TRUE)

  if (!is.null(y)) {
    if (length(y) != n) stop("y must have length nrow(X).")
  }

  if (identical(scoring, "auto")) {
    scoring <- if (is.null(y)) "unsupervised" else "accuracy"
  }
  if (identical(scoring, "accuracy") && is.null(y)) {
    stop("scoring='accuracy' requires non-NULL y.")
  }

  if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) || top_n < 1L) {
    stop("top_n must be a positive integer.")
  }
  top_n <- as.integer(top_n)

  fold_y <- if (identical(scoring, "accuracy")) y else NULL
  folds <- .make_cv_folds(n = n, k = cv_folds, y = fold_y, seed = seed)

  build_comp_grid <- function(method) {
    if (!.is_reduction_method(method)) {
      return(list(list(n_comp = NULL, var_threshold = NULL)))
    }

    out <- list(list(n_comp = NULL, var_threshold = NULL))

    if (!is.null(n_comp_grid)) {
      for (k in n_comp_grid) {
        out[[length(out) + 1L]] <- list(n_comp = k, var_threshold = NULL)
      }
    }

    if (!is.null(var_threshold_grid)) {
      for (vt in var_threshold_grid) {
        out[[length(out) + 1L]] <- list(n_comp = NULL, var_threshold = vt)
      }
    }

    tags <- vapply(out, function(x) {
      if (!is.null(x$n_comp)) {
        paste0("k:", x$n_comp)
      } else {
        paste0("vt:", ifelse(is.null(x$var_threshold), "none", format(x$var_threshold, digits = 8)))
      }
    }, character(1))

    out[!duplicated(tags)]
  }

  candidates <- list()
  idx <- 1L
  for (method in methods) {
    comp_grid <- build_comp_grid(method)
    for (comp_cfg in comp_grid) {
      for (lam in lambda_values) {
        for (cov_estimator in cov_estimator_grid) {
          if (is.character(lam) && identical(lam, "auto") && (!identical(cov_estimator, "empirical") || !is.null(sample_weight))) {
            next
          }

          candidates[[idx]] <- list(
            method = method,
            n_comp = comp_cfg$n_comp,
            var_threshold = comp_cfg$var_threshold,
            lambda = lam,
            cov_estimator = cov_estimator
          )
          idx <- idx + 1L
        }
      }
    }
  }

  if (length(candidates) == 0L) {
    stop("No valid tuning candidates were generated. Check lambda/cov_estimator/sample_weight compatibility.")
  }

  eval_candidate <- function(cfg) {
    fold_scores <- rep(NA_real_, length(folds))
    n_success <- 0L

    for (i in seq_along(folds)) {
      valid_idx <- folds[[i]]
      train_idx <- setdiff(seq_len(n), valid_idx)

      X_train <- X[train_idx, , drop = FALSE]
      X_valid <- X[valid_idx, , drop = FALSE]
      w_train <- if (is.null(sample_weight)) NULL else sample_weight[train_idx]

      fit <- tryCatch(
        whiten_model(
          X_train,
          center = TRUE,
          method = cfg$method,
          n_comp = cfg$n_comp,
          var_threshold = cfg$var_threshold,
          lambda = cfg$lambda,
          lambda_method = lambda_method,
          sample_weight = w_train,
          cov_estimator = cfg$cov_estimator,
          tyler_tol = tyler_tol,
          tyler_max_iter = tyler_max_iter,
          tyler_eps = tyler_eps,
          eig_method = eig_method,
          fast = fast,
          sign_reference = sign_reference
        ),
        error = function(e) e
      )

      if (inherits(fit, "error")) {
        next
      }

      Z_train <- tryCatch(stats::predict(fit, X_train), error = function(e) e)
      Z_valid <- tryCatch(stats::predict(fit, X_valid), error = function(e) e)
      if (inherits(Z_train, "error") || inherits(Z_valid, "error")) {
        next
      }

      score <- tryCatch({
        if (!is.null(score_fn)) {
          if (identical(scoring, "accuracy")) {
            score_fn(Z_train, y[train_idx], Z_valid, y[valid_idx])
          } else {
            score_fn(Z_valid)
          }
        } else if (identical(scoring, "accuracy")) {
          .default_supervised_score(Z_train, y[train_idx], Z_valid, y[valid_idx])
        } else {
          .default_whitening_score(Z_valid)
        }
      }, error = function(e) NA_real_)

      if (is.numeric(score) && length(score) == 1L && is.finite(score)) {
        fold_scores[[i]] <- as.numeric(score)
        n_success <- n_success + 1L
      }
    }

    valid <- is.finite(fold_scores)
    list(
      scores = fold_scores,
      mean_score = if (any(valid)) mean(fold_scores[valid]) else -Inf,
      sd_score = if (sum(valid) > 1L) stats::sd(fold_scores[valid]) else NA_real_,
      n_success = n_success
    )
  }

  evals <- lapply(candidates, eval_candidate)

  ranking <- do.call(
    rbind,
    lapply(seq_along(candidates), function(i) {
      cfg <- candidates[[i]]
      ev <- evals[[i]]
      data.frame(
        method = cfg$method,
        n_comp = if (is.null(cfg$n_comp)) NA_integer_ else as.integer(cfg$n_comp),
        var_threshold = if (is.null(cfg$var_threshold)) NA_real_ else as.numeric(cfg$var_threshold),
        lambda = .lambda_to_label(cfg$lambda),
        cov_estimator = cfg$cov_estimator,
        mean_score = ev$mean_score,
        sd_score = ev$sd_score,
        n_success = ev$n_success,
        stringsAsFactors = FALSE
      )
    })
  )

  ranking <- ranking[order(-ranking$mean_score, ranking$sd_score, na.last = TRUE), , drop = FALSE]
  rownames(ranking) <- NULL
  if (nrow(ranking) > top_n) {
    ranking <- ranking[seq_len(top_n), , drop = FALSE]
  }

  best_idx <- which.max(vapply(evals, function(e) e$mean_score, numeric(1)))
  if (!is.finite(evals[[best_idx]]$mean_score)) {
    stop("All tuning candidates failed during cross-validation.")
  }
  best_cfg <- candidates[[best_idx]]

  best_model <- whiten_model(
    X,
    center = TRUE,
    method = best_cfg$method,
    n_comp = best_cfg$n_comp,
    var_threshold = best_cfg$var_threshold,
    lambda = best_cfg$lambda,
    lambda_method = lambda_method,
    sample_weight = sample_weight,
    cov_estimator = best_cfg$cov_estimator,
    tyler_tol = tyler_tol,
    tyler_max_iter = tyler_max_iter,
    tyler_eps = tyler_eps,
    eig_method = eig_method,
    fast = fast,
    sign_reference = sign_reference
  )

  out <- list(
    best_model = best_model,
    best_params = best_cfg,
    ranking = ranking,
    scoring = scoring,
    cv_folds = as.integer(cv_folds),
    seed = seed,
    lambda_method = lambda_method
  )
  class(out) <- "whiten_tune"
  out
}

#' One-Click Whitening Report
#'
#' Generates a markdown report with model settings and whitening diagnostics.
#'
#' @param model A \code{whiten_model} or \code{whiten_tune} object.
#' @param data Optional original-space data matrix used to compute
#'   diagnostics.
#' @param Z Optional whitened matrix. If missing and \code{data} is provided,
#'   \code{predict(model, data)} is used.
#' @param tune_result Optional tuning result from
#'   \code{auto_tune_whitening()}.
#' @param file Output markdown path. If \code{NULL}, returns report text
#'   without writing to disk.
#' @param digits Numeric formatting precision.
#'
#' @return If \code{file} is non-NULL, invisibly returns the output path.
#'   Otherwise returns the markdown text.
#'
#' @seealso \code{\link{whiten_model}}, \code{\link{auto_tune_whitening}}
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(200 * 6), 200, 6)
#' m <- whiten_model(X, method = "ZCA", lambda = 0.1)
#' txt <- report_whitening(m, data = X, file = NULL)
#' cat(txt)
#'
#' @export
report_whitening <- function(model,
                             data = NULL,
                             Z = NULL,
                             tune_result = NULL,
                             file = "whitening-report.md",
                             digits = 4) {
  if (inherits(model, "whiten_tune")) {
    tune_result <- model
    model <- model$best_model
  }

  if (!inherits(model, "whiten_model")) {
    stop("model must be a 'whiten_model' or 'whiten_tune' object.")
  }

  if (!is.null(data)) {
    if (!is.matrix(data) || !is.numeric(data) || any(!is.finite(data))) {
      stop("data must be a finite numeric matrix.")
    }
  }

  if (is.null(Z) && !is.null(data)) {
    Z <- stats::predict(model, data)
  }

  diagnostics <- NULL
  if (!is.null(Z)) {
    diagnostics <- check_whitening(Z)
  }

  fmt <- function(x) {
    if (is.null(x)) return("NULL")
    if (length(x) == 1L && is.na(x)) return("NA")
    if (is.numeric(x)) return(format(signif(x, digits = digits), trim = TRUE, scientific = FALSE))
    as.character(x)
  }

  lines <- c(
    "# Whitening Report",
    "",
    sprintf("- Generated at: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("- Method: %s", model$method),
    sprintf("- Dimensions: %s -> %s", fmt(model$dim_in), fmt(model$dim_out)),
    sprintf("- n_comp: %s", fmt(model$n_comp)),
    sprintf("- var_threshold: %s", fmt(model$var_threshold)),
    sprintf("- explained_var: %s", fmt(model$explained_var)),
    sprintf("- lambda: %s", fmt(model$lambda)),
    sprintf("- lambda_input: %s", fmt(model$lambda_input)),
    sprintf("- lambda_method: %s", fmt(model$lambda_method)),
    sprintf("- cov_estimator: %s", fmt(model$cov_estimator)),
    sprintf("- sample_weighted: %s", fmt(model$sample_weighted)),
    sprintf("- eig_method: %s", fmt(model$eig_method)),
    sprintf("- fast: %s", fmt(model$fast)),
    sprintf("- n_obs: %s", fmt(model$n_obs)),
    sprintf("- update_count: %s", fmt(model$update_count))
  )

  lines <- c(lines, "")
  lines <- c(lines, "## Whitening Diagnostics", "")
  if (is.null(diagnostics)) {
    lines <- c(lines, "Diagnostics unavailable (provide `data` or `Z`).")
  } else {
    lines <- c(
      lines,
      sprintf("- diag_mean: %s", fmt(diagnostics$diag_mean)),
      sprintf("- diag_min: %s", fmt(min(diagnostics$diag_range))),
      sprintf("- diag_max: %s", fmt(max(diagnostics$diag_range))),
      sprintf("- offdiag_frob: %s", fmt(diagnostics$offdiag_frob))
    )
  }

  tune_obj <- if (!is.null(tune_result) && inherits(tune_result, "whiten_tune")) tune_result else NULL
  if (!is.null(tune_obj)) {
    lines <- c(lines, "", "## Tuning Summary", "")
    lines <- c(lines, sprintf("- Scoring: %s", tune_obj$scoring))
    lines <- c(lines, sprintf("- CV folds: %s", tune_obj$cv_folds))
    lines <- c(lines, "")

    tab <- tune_obj$ranking
    if (is.null(tab) || nrow(tab) == 0L) {
      lines <- c(lines, "No ranking rows available.")
    } else {
      lines <- c(lines, "| method | n_comp | var_threshold | lambda | cov_estimator | mean_score | sd_score | n_success |")
      lines <- c(lines, "|---|---:|---:|---|---|---:|---:|---:|")
      for (i in seq_len(nrow(tab))) {
        lines <- c(lines, sprintf(
          "| %s | %s | %s | %s | %s | %s | %s | %s |",
          fmt(tab$method[[i]]),
          fmt(tab$n_comp[[i]]),
          fmt(tab$var_threshold[[i]]),
          fmt(tab$lambda[[i]]),
          fmt(tab$cov_estimator[[i]]),
          fmt(tab$mean_score[[i]]),
          fmt(tab$sd_score[[i]]),
          fmt(tab$n_success[[i]])
        ))
      }
    }
  }

  text <- paste(lines, collapse = "\n")
  if (is.null(file)) {
    return(text)
  }

  out_dir <- dirname(file)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  writeLines(lines, con = file, useBytes = TRUE)
  invisible(file)
}
