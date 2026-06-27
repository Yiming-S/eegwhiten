############################################################
# S3 methods for whiten_model objects
############################################################

#' Print a whitening model
#'
#' @param x A \code{whiten_model} object.
#' @param ... Unused.
#'
#' @return The input object invisibly.
#'
#' @examples
#' m <- whiten_model(matrix(rnorm(200 * 6), 200, 6), method = "PCA",
#'                   lambda = 0.1)
#' print(m)
#'
#' @export
#' @method print whiten_model
print.whiten_model <- function(x, ...) {
  if (!inherits(x, "whiten_model")) stop("x must be a 'whiten_model' object.")

  show_val <- function(v, digits = NULL) {
    if (is.null(v)) return("NULL")
    if (length(v) == 1L && is.na(v)) return("NA")
    if (!is.null(digits) && is.numeric(v)) return(signif(v, digits))
    v
  }

  lambda_input <- if (!is.null(x$lambda_input)) x$lambda_input else as.character(x$lambda)
  lambda_method <- if (!is.null(x$lambda_method)) x$lambda_method else NA_character_
  eig_method <- if (!is.null(x$eig_method)) x$eig_method else NA_character_
  fast_mode <- if (!is.null(x$fast)) x$fast else NA
  cov_estimator <- if (!is.null(x$cov_estimator)) x$cov_estimator else NA_character_
  sample_weighted <- if (!is.null(x$sample_weighted)) x$sample_weighted else NA
  n_obs <- if (!is.null(x$n_obs)) x$n_obs else NA_integer_
  update_count <- if (!is.null(x$update_count)) x$update_count else 0L

  cat("<whiten_model>\n")
  cat("  method      :", x$method, "\n")
  cat("  dimensions  :", x$dim_in, "->", x$dim_out, "\n")
  cat("  n_comp      :", show_val(x$n_comp), "\n")
  cat("  var_threshold:", show_val(x$var_threshold), "\n")
  cat("  explained_var:", show_val(x$explained_var, digits = 4), "\n")
  cat("  lambda      :", show_val(x$lambda, digits = 4), "\n")
  cat("  shrink_target:", show_val(if (!is.null(x$shrink_target)) x$shrink_target else "identity"), "\n")
  cat("  eig_method  :", show_val(eig_method), "\n")
  cat("  fast_mode   :", show_val(fast_mode), "\n")
  cat("  cov_estimator:", show_val(cov_estimator), "\n")
  cat("  weighted_fit:", show_val(sample_weighted), "\n")
  cat("  n_obs       :", show_val(n_obs), "\n")
  cat("  updates     :", show_val(update_count), "\n")
  if (identical(lambda_input, "auto")) {
    cat("  lambda_mode :", "auto (", lambda_method, ")\n", sep = "")
  }

  invisible(x)
}

#' Summarize a whitening model
#'
#' @param object A \code{whiten_model} object.
#' @param data Optional numeric matrix in original feature space.
#' @param Z Optional whitened data matrix. If missing and \code{data} is
#'   provided, \code{predict(object, data)} will be used.
#' @param ... Unused.
#'
#' @return An object of class \code{"summary.whiten_model"}.
#'
#' @examples
#' X <- matrix(rnorm(200 * 6), 200, 6)
#' m <- whiten_model(X, method = "PCA", lambda = 0.1)
#' s <- summary(m, data = X)
#' print(s)
#'
#' @export
#' @method summary whiten_model
summary.whiten_model <- function(object, data = NULL, Z = NULL, ...) {
  if (!inherits(object, "whiten_model")) {
    stop("object must be a 'whiten_model' object.")
  }

  if (is.null(Z) && !is.null(data)) {
    Z <- stats::predict(object, data)
  }

  diagnostics <- NULL
  if (!is.null(Z)) {
    diagnostics <- check_whitening(Z)
  }

  out <- list(
    method = object$method,
    dim_in = object$dim_in,
    dim_out = object$dim_out,
    n_comp = object$n_comp,
    var_threshold = object$var_threshold,
    explained_var = object$explained_var,
    lambda = object$lambda,
    lambda_input = if (!is.null(object$lambda_input)) object$lambda_input else as.character(object$lambda),
    lambda_method = if (!is.null(object$lambda_method)) object$lambda_method else NA_character_,
    shrink_target = if (!is.null(object$shrink_target)) object$shrink_target else "identity",
    eig_method = if (!is.null(object$eig_method)) object$eig_method else NA_character_,
    fast = if (!is.null(object$fast)) object$fast else NA,
    cov_estimator = if (!is.null(object$cov_estimator)) object$cov_estimator else NA_character_,
    sample_weighted = if (!is.null(object$sample_weighted)) object$sample_weighted else NA,
    n_obs = if (!is.null(object$n_obs)) object$n_obs else NA_integer_,
    update_count = if (!is.null(object$update_count)) object$update_count else 0L,
    diagnostics = diagnostics
  )

  class(out) <- "summary.whiten_model"
  out
}

#' Print summary for a whitening model
#'
#' @param x A \code{summary.whiten_model} object.
#' @param ... Unused.
#'
#' @return The input object invisibly.
#' @export
#' @method print summary.whiten_model
print.summary.whiten_model <- function(x, ...) {
  show_val <- function(v, digits = NULL) {
    if (is.null(v)) return("NULL")
    if (length(v) == 1L && is.na(v)) return("NA")
    if (!is.null(digits) && is.numeric(v)) return(signif(v, digits))
    v
  }

  cat("<summary.whiten_model>\n")
  cat("  method      :", x$method, "\n")
  cat("  dimensions  :", x$dim_in, "->", x$dim_out, "\n")
  cat("  n_comp      :", show_val(x$n_comp), "\n")
  cat("  var_threshold:", show_val(x$var_threshold), "\n")
  cat("  explained_var:", show_val(x$explained_var, digits = 4), "\n")
  cat("  lambda      :", show_val(x$lambda, digits = 4), "\n")
  cat("  shrink_target:", show_val(if (!is.null(x$shrink_target)) x$shrink_target else "identity"), "\n")
  cat("  eig_method  :", show_val(x$eig_method), "\n")
  cat("  fast_mode   :", show_val(x$fast), "\n")
  cat("  cov_estimator:", show_val(x$cov_estimator), "\n")
  cat("  weighted_fit:", show_val(x$sample_weighted), "\n")
  cat("  n_obs       :", show_val(x$n_obs), "\n")
  cat("  updates     :", show_val(x$update_count), "\n")
  if (identical(x$lambda_input, "auto")) {
    cat("  lambda_mode :", "auto (", x$lambda_method, ")\n", sep = "")
  }

  if (!is.null(x$diagnostics)) {
    cat("  whitening   : diag_mean=", signif(x$diagnostics$diag_mean, 4),
        ", offdiag_frob=", signif(x$diagnostics$offdiag_frob, 4), "\n", sep = "")
  }

  invisible(x)
}

#' Plot covariance deviation from identity for whitened data
#'
#' @param x A \code{whiten_model} object.
#' @param data Optional numeric matrix in original feature space.
#' @param Z Optional whitened data matrix.
#' @param main Optional plot title.
#' @param ... Passed to \code{graphics::image()}.
#'
#' @return Invisibly returns \code{cov(Z) - I}.
#' @export
#' @method plot whiten_model
plot.whiten_model <- function(x, data = NULL, Z = NULL, main = NULL, ...) {
  if (!inherits(x, "whiten_model")) stop("x must be a 'whiten_model' object.")

  if (is.null(Z)) {
    if (is.null(data)) {
      stop("Provide either Z or data to plot whitening diagnostics.")
    }
    Z <- stats::predict(x, data)
  }
  if (!is.matrix(Z) || !is.numeric(Z) || any(!is.finite(Z))) {
    stop("Z must be a finite numeric matrix.")
  }

  S <- stats::cov(Z)
  d <- ncol(S)
  dev <- S - diag(d)

  if (is.null(main)) {
    main <- sprintf("Cov(Z) - I (%s)", x$method)
  }

  lim <- max(abs(dev))
  if (!is.finite(lim) || lim <= .Machine$double.eps) {
    lim <- 1e-12
  }

  cols <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(64)
  graphics::image(
    x = seq_len(d),
    y = seq_len(d),
    z = t(dev[d:1, , drop = FALSE]),
    zlim = c(-lim, lim),
    col = cols,
    xlab = "Component",
    ylab = "Component",
    main = main,
    ...
  )
  graphics::box()

  invisible(dev)
}
