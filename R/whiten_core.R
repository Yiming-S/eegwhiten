############################################################
# FILENAME: whiten_core.R
# DESCRIPTION: Core whitening API for EEG data
# UPDATES: Added Regularization (lambda) and Dim Reduction (n_comp)
############################################################

#' Fit a whitening matrix from a covariance matrix
#'
#' Compute a whitening matrix \code{W} from a symmetric
#' positive-definite covariance matrix using one of the
#' implemented methods.
#'
#' @param Sigma Symmetric positive-definite covariance matrix.
#' @param method Whitening method.
#'   Underscore aliases \code{"ZCA_cor"} and \code{"PCA_cor"} are also accepted.
#' @param n_comp Integer (optional); Number of components to keep
#'   (for PCA/PCA-cor/SVD).
#' @param var_threshold Numeric in (0, 1]; cumulative explained variance
#'   threshold used to choose \code{n_comp} automatically for
#'   \code{"PCA"}, \code{"PCA-cor"}, and \code{"SVD"}.
#' @param lambda Numeric in [0, 1]; optional shrinkage applied directly to
#'   \code{Sigma} toward a scaled identity target before whitening.
#' @param eig_method Eigen solver backend; one of \code{"auto"},
#'   \code{"base"}, or \code{"rspectra"}.
#' @param fast Logical; if \code{TRUE}, allow faster approximate settings
#'   for iterative eigensolvers.
#' @param sign_reference Optional reference vectors used for stable component
#'   sign orientation in \code{"PCA"}, \code{"PCA-cor"}, and \code{"SVD"}.
#'
#' @return Whitening matrix \code{W}.
#'
#' @export
whiten_fit <- function(Sigma,
                       method = c("SVD", "ZCA", "ZCA-cor",
                                  "PCA", "PCA-cor", "Cholesky"),
                       n_comp = NULL,
                       var_threshold = NULL,
                       lambda = 0,
                       eig_method = c("auto", "base", "rspectra"),
                       fast = FALSE,
                       sign_reference = NULL) {
  if (length(method) > 1L) method <- method[[1L]]
  method <- .normalize_whiten_method(method)
  eig_method <- match.arg(eig_method)
  .check_var_threshold(var_threshold)
  .check_unit_interval(lambda, "lambda")
  if (!is.logical(fast) || length(fast) != 1L || is.na(fast)) {
    stop("fast must be TRUE or FALSE.")
  }
  is_reduction <- .is_reduction_method(method)

  if (!is.null(n_comp)) {
    if (!is.numeric(n_comp) || length(n_comp) != 1L || !is.finite(n_comp)) {
      stop("n_comp must be a single finite numeric value.")
    }
    if (n_comp != as.integer(n_comp) || n_comp < 1L) {
      stop("n_comp must be a positive integer.")
    }
    n_comp <- as.integer(n_comp)
  }
  if (!is.null(n_comp) && !is.null(var_threshold)) {
    stop("Specify only one of n_comp or var_threshold.")
  }
  
  # All implemented methods expect symmetric PD covariance
  eig <- .check_symmetric_pd(Sigma, require_pd = TRUE)
  d <- ncol(Sigma)
  cond_num <- max(eig) / min(eig)
  if (lambda == 0 && is.finite(cond_num) && cond_num > 1e8) {
    warning(sprintf("Sigma is ill-conditioned (condition number ~ %.3e); consider lambda > 0 for stable whitening.", cond_num))
  }

  S <- Sigma
  if (lambda > 0) {
    target_val <- mean(diag(Sigma))
    target_mat <- diag(target_val, d)
    S <- (1 - lambda) * Sigma + lambda * target_mat
  }

  if (is_reduction) {
    if (is.null(n_comp) && !is.null(var_threshold)) {
      spectrum_info <- if (method %in% c("PCA", "SVD") && lambda == 0) {
        list(values = pmax(eig, 0), total = sum(diag(S)))
      } else {
        .component_spectrum_info(
          S,
          method,
          eig_method = eig_method,
          fast = fast
        )
      }
      n_comp <- .n_comp_from_threshold(spectrum_info$values, var_threshold)
    }
    if (!is.null(n_comp) && (n_comp < 1L || n_comp > d)) {
      stop("n_comp must be between 1 and ncol(Sigma).")
    }
  } else {
    if (!is.null(n_comp)) {
      warning("n_comp is ignored for ZCA/ZCA-cor/Cholesky.")
      n_comp <- NULL
    }
    if (!is.null(var_threshold)) {
      warning("var_threshold is ignored for ZCA/ZCA-cor/Cholesky.")
    }
  }
  attr(S, ".checked_spd") <- TRUE

  args <- list(Sigma = S, returnW = TRUE, PhiPsi = FALSE)
  
  # Only pass n_comp to methods that support it to avoid "unused argument" errors
  if (is_reduction && !is.null(n_comp)) {
    args$n_comp <- n_comp
  }
  if (is_reduction && !is.null(sign_reference)) {
    args$sign_ref <- sign_reference
  }
  if (method %in% c("PCA", "PCA-cor", "SVD", "ZCA", "ZCA-cor")) {
    args$eig_method <- eig_method
    args$fast <- fast
  }
  
  # Call the specific function (e.g., PCA(Sigma, ...))
  res <- do.call(.method_to_function(method), args)
  
  if (is.null(res$W)) stop("Whitening matrix W is NULL.")
  res$W
}

#' Fit a whitening model for EEG data
#'
#' The input matrix \code{X} is assumed to have trials/epochs in rows
#' and channels or features in columns. This function estimates
#' a centering vector and a whitening matrix, optionally using
#' regularization (shrinkage) and dimensionality reduction.
#'
#' @param X Numeric matrix [n_trials x n_channels].
#' @param center Logical; whether to subtract column means before
#'   computing the covariance.
#' @param method Whitening method; one of \code{"SVD"}, \code{"ZCA"},
#'   \code{"ZCA-cor"}, \code{"PCA"}, \code{"PCA-cor"}, \code{"Cholesky"}.
#'   Underscore aliases \code{"ZCA_cor"} and \code{"PCA_cor"} are also accepted.
#' @param n_comp Integer (optional); Number of components to keep.
#'   Only applicable for "PCA", "PCA-cor", and "SVD".
#' @param var_threshold Numeric in (0, 1]; cumulative explained variance
#'   threshold for automatic component selection. Ignored by non-reduction
#'   methods.
#' @param lambda Numeric (0 to 1) or \code{"auto"}; regularization parameter
#'   for covariance shrinkage.
#'   0 = No regularization (Empirical covariance).
#'   1 = Full shrinkage to target (scaled Identity).
#'   Values between 0 and 1 mix the two. Useful for high-dimensional/low-sample EEG data.
#' @param lambda_method Method used when \code{lambda = "auto"}.
#'   One of \code{"oas"} or \code{"lw"}.
#' @param sample_weight Optional non-negative sample weights with length
#'   \code{nrow(X)}.
#' @param cov_estimator Covariance estimator; one of \code{"empirical"},
#'   \code{"mcd"}, or \code{"tyler"}.
#' @param tyler_tol Convergence tolerance for \code{cov_estimator = "tyler"}.
#' @param tyler_max_iter Maximum iterations for Tyler covariance estimation.
#' @param tyler_eps Numerical stabilization floor for Tyler estimation.
#' @param eig_method Eigen solver backend; one of \code{"auto"},
#'   \code{"base"}, or \code{"rspectra"}.
#' @param fast Logical; if \code{TRUE}, allow faster approximate settings
#'   for iterative eigensolvers.
#' @param sign_reference Optional reference vectors used for stable component
#'   sign orientation in \code{"PCA"}, \code{"PCA-cor"}, and \code{"SVD"}.
#'
#' @return An object of class \code{"whiten_model"}.
#'
#' @export
whiten_model <- function(X, center = TRUE,
                         method = c("SVD", "ZCA", "ZCA-cor",
                                    "PCA", "PCA-cor", "Cholesky"),
                         n_comp = NULL,
                         var_threshold = NULL,
                         lambda = "auto",
                         lambda_method = c("oas", "lw"),
                         sample_weight = NULL,
                         cov_estimator = c("empirical", "mcd", "tyler"),
                         tyler_tol = 1e-6,
                         tyler_max_iter = 200,
                         tyler_eps = 1e-6,
                         eig_method = c("auto", "base", "rspectra"),
                         fast = FALSE,
                         sign_reference = NULL) {
  if (length(method) > 1L) method <- method[[1L]]
  method <- .normalize_whiten_method(method)
  lambda_method <- match.arg(lambda_method)
  cov_estimator <- match.arg(cov_estimator)
  eig_method <- match.arg(eig_method)
  if (!is.logical(fast) || length(fast) != 1L || is.na(fast)) {
    stop("fast must be TRUE or FALSE.")
  }
  is_reduction <- .is_reduction_method(method)
  .check_var_threshold(var_threshold)

  if (!is.matrix(X)) stop("X must be a matrix.")
  if (!is.numeric(X)) stop("X must be a numeric matrix.")
  if (any(!is.finite(X))) stop("X must contain only finite values.")
  if (!is.logical(center) || length(center) != 1L || is.na(center)) {
    stop("center must be TRUE or FALSE.")
  }

  n <- nrow(X)
  d <- ncol(X)

  if (n <= 1L) {
    stop("X must have at least 2 rows to estimate covariance.")
  }
  if (d <= 1L) {
    stop("X must have more than one column.")
  }

  .check_unit_interval(lambda, "lambda", allow_char_auto = TRUE)
  lambda_auto <- is.character(lambda)
  sample_weight <- .validate_sample_weight(sample_weight, n, allow_null = TRUE)
  if (lambda_auto && !is.null(sample_weight)) {
    stop("lambda='auto' is not supported with sample_weight. Set a numeric lambda in [0, 1].")
  }
  if (lambda_auto && cov_estimator != "empirical") {
    stop("lambda='auto' is only supported with cov_estimator='empirical'.")
  }

  # Validation for n_comp
  if (!is.null(n_comp)) {
    if (!is.numeric(n_comp) || length(n_comp) != 1L || !is.finite(n_comp)) {
      stop("n_comp must be a single finite numeric value.")
    }
    if (n_comp != as.integer(n_comp)) {
      stop("n_comp must be an integer.")
    }
    n_comp <- as.integer(n_comp)
    if (n_comp < 1L || n_comp > d) {
      stop("n_comp must be between 1 and ncol(X).")
    }
  }
  if (!is.null(n_comp) && !is.null(var_threshold)) {
    stop("Specify only one of n_comp or var_threshold.")
  }
  if (!is_reduction) {
    if (!is.null(n_comp)) {
      warning("n_comp is ignored for ZCA/ZCA-cor/Cholesky. Use PCA/PCA-cor/SVD for dimension reduction.")
      n_comp <- NULL
    }
    if (!is.null(var_threshold)) {
      warning("var_threshold is ignored for ZCA/ZCA-cor/Cholesky.")
      var_threshold <- NULL
    }
  }
  if (!is.null(sign_reference)) {
    if (!is_reduction) {
      warning("sign_reference is ignored for ZCA/ZCA-cor/Cholesky.")
      sign_reference <- NULL
    } else {
      if (!is.matrix(sign_reference)) stop("sign_reference must be a matrix.")
      if (nrow(sign_reference) != d) {
        stop("sign_reference must have nrow(sign_reference) == ncol(X).")
      }
    }
  }

  # Ensure column names exist
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("V", seq_len(d))
  }
  
  # 1. Covariance estimation
  cov_fit <- .estimate_covariance(
    X,
    center = center,
    sample_weight = sample_weight,
    cov_estimator = cov_estimator,
    tyler_tol = tyler_tol,
    tyler_max_iter = tyler_max_iter,
    tyler_eps = tyler_eps
  )
  mu <- cov_fit$mean
  names(mu) <- colnames(X)

  # 2. Optional shrinkage
  S <- cov_fit$cov
  if (lambda_auto) {
    Xc <- if (center) sweep(X, 2, mu, "-") else X
    lambda <- if (lambda_method == "oas") .lambda_oas(Xc) else .lambda_lw(Xc)
  }

  # Apply shrinkage toward scaled identity
  if (lambda > 0) {
    target_val <- mean(diag(S))
    target_mat <- diag(target_val, d)
    S <- (1 - lambda) * S + lambda * target_mat
  }

  need_full_eigs <- !is_reduction || !is.null(var_threshold) || is.null(n_comp) || identical(eig_method, "base")
  eig_shrunk <- .check_symmetric_pd(S, require_pd = TRUE, return_eigen = need_full_eigs)
  cond_raw <- if (need_full_eigs) {
    max(eig_shrunk) / min(pmax(eig_shrunk, .Machine$double.eps))
  } else {
    inv_rcond <- tryCatch(1 / rcond(S), error = function(e) NA_real_)
    inv_rcond
  }
  if (lambda == 0 && is.finite(cond_raw) && cond_raw > 1e8) {
    warning(sprintf("Sample covariance is ill-conditioned (condition number ~ %.3e); consider lambda='auto' or lambda > 0.", cond_raw))
  }
  attr(S, ".checked_spd") <- TRUE

  n_comp_final <- n_comp
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
  
  # 3. Fit Whitening Matrix
  fit_args <- list(Sigma = S, returnW = TRUE, PhiPsi = FALSE, return_decomp = TRUE)
  if (is_reduction && !is.null(n_comp_final)) {
    fit_args$n_comp <- n_comp_final
  }
  if (is_reduction && !is.null(sign_reference)) {
    fit_args$sign_ref <- sign_reference
  }
  if (method %in% c("PCA", "PCA-cor", "SVD", "ZCA", "ZCA-cor")) {
    fit_args$eig_method <- eig_method
    fit_args$fast <- fast
  }
  
  res <- do.call(.method_to_function(method), fit_args)
  W   <- res$W
  decomp <- if (!is.null(res$decomp)) res$decomp else list()
  inv_tW <- if (!is.null(decomp$inv_tW)) decomp$inv_tW else .pinv(t(W))
  sign_basis <- if (!is.null(decomp$U)) decomp$U else NULL
  
  structure(
    list(
      W      = W,
      Wt     = t(W),
      center = mu,
      center_enabled = center,
      method = method,
      n_comp = n_comp_final,
      var_threshold = var_threshold,
      explained_var = explained_var,
      lambda = lambda,
      lambda_input = if (lambda_auto) "auto" else as.character(lambda),
      lambda_method = if (lambda_auto) lambda_method else NA_character_,
      eig_method = eig_method,
      fast = fast,
      cov_estimator = cov_estimator,
      sample_weighted = !is.null(sample_weight),
      U = decomp$U,
      D = decomp$D,
      inv_tW = inv_tW,
      sign_basis = sign_basis,
      cov = S,
      n_obs = n,
      stat_sum_w = cov_fit$sum_w,
      stat_sum_w2 = cov_fit$sum_w2,
      stat_sum_x = cov_fit$sum_x,
      stat_sum_xx = cov_fit$sum_xx,
      update_count = 0L,
      dim_in = d,
      dim_out = nrow(W)
    ),
    class = "whiten_model"
  )
}

#' Apply a whitening model to new EEG data
#'
#' @param object A \code{whiten_model} object.
#' @param newdata Numeric matrix with the same number of columns
#'   as the data used to fit the model.
#' @param ... Unused.
#'
#' @return Whitened data matrix (possibly with reduced dimensions).
#'
#' @export
#' @method predict whiten_model
predict.whiten_model <- function(object, newdata, ...) {
  if (!inherits(object, "whiten_model")) stop("object must be a 'whiten_model'.")
  if (!is.matrix(newdata)) stop("newdata must be a matrix.")
  if (!is.numeric(newdata)) stop("newdata must be a numeric matrix.")
  if (any(!is.finite(newdata))) stop("newdata must contain only finite values.")

  d <- length(object$center)
  if (ncol(newdata) != d) {
    stop(sprintf("Mismatch: Model expects %d columns, but newdata has %d.", d, ncol(newdata)))
  }

  # Center using stored mean
  Xc <- sweep(newdata, 2, object$center, "-")

  # Z = centered X times whitening matrix
  Wt <- if (!is.null(object$Wt)) object$Wt else t(object$W)
  Z <- Xc %*% Wt

  # Update column names based on output dimension
  k <- ncol(Z)
  if (!is.null(rownames(object$W))) {
    colnames(Z) <- rownames(object$W)
  } else {
    colnames(Z) <- paste0("Comp", seq_len(k))
  }
  attr(Z, "method") <- object$method

  Z
}

#' Inverse transform from whitened space back to original EEG space
#'
#' @param Z Whitened data matrix.
#' @param model A \code{whiten_model} object.
#'
#' @return Approximate reconstruction of the original data matrix.
#'
#' @export
unwhiten <- function(Z, model) {
  if (!inherits(model, "whiten_model")) {
    stop("model must be a 'whiten_model' object.")
  }
  if (!is.matrix(Z)) stop("Z must be a matrix.")
  if (!is.numeric(Z)) stop("Z must be a numeric matrix.")
  if (any(!is.finite(Z))) stop("Z must contain only finite values.")

  k <- if (!is.null(model$dim_out)) model$dim_out else nrow(model$W)
  if (ncol(Z) != k) {
    stop(sprintf("Mismatch: model expects %d whitened components, but Z has %d.", k, ncol(Z)))
  }

  A_inv <- model$inv_tW
  if (is.null(A_inv)) {
    A <- t(model$W)
    if (nrow(A) == ncol(A)) {
      A_inv <- tryCatch(
        solve(A),
        error = function(e) {
          warning("Whitening matrix is ill-conditioned; using pseudo-inverse for reconstruction.")
          .pinv(A)
        }
      )
    } else {
      A_inv <- .pinv(A)
    }
  }

  Xc <- Z %*% A_inv
  X  <- sweep(Xc, 2, model$center, "+")

  colnames(X) <- names(model$center)
  X
}

#' Optimized inverse transform from whitened space
#'
#' Uses cached decomposition terms from \code{whiten_model} when available.
#'
#' @param Z Whitened data matrix.
#' @param model A \code{whiten_model} object.
#'
#' @return Approximate reconstruction of the original data matrix.
#' @export
unwhiten_fast <- function(Z, model) {
  if (!inherits(model, "whiten_model")) {
    stop("model must be a 'whiten_model' object.")
  }
  if (!is.matrix(Z)) stop("Z must be a matrix.")
  if (!is.numeric(Z)) stop("Z must be a numeric matrix.")
  if (any(!is.finite(Z))) stop("Z must contain only finite values.")

  k <- if (!is.null(model$dim_out)) model$dim_out else nrow(model$W)
  if (ncol(Z) != k) {
    stop(sprintf("Mismatch: model expects %d whitened components, but Z has %d.", k, ncol(Z)))
  }

  if (!is.null(model$inv_tW)) {
    Xc <- Z %*% model$inv_tW
    X <- sweep(Xc, 2, model$center, "+")
    colnames(X) <- names(model$center)
    return(X)
  }

  unwhiten(Z, model)
}

#' Diagnostic check for whitening quality
#'
#' @param Z Whitened data matrix.
#'
#' @return List of diagnostics.
#' @export
check_whitening <- function(Z) {
  if (!is.matrix(Z)) stop("Z must be a matrix.")
  if (!is.numeric(Z)) stop("Z must be a numeric matrix.")
  if (any(!is.finite(Z))) stop("Z must contain only finite values.")
  
  S <- stats::cov(Z)
  d <- ncol(Z)
  
  diag_vals  <- diag(S)
  diag_range <- range(diag_vals)
  
  off_diag <- S
  diag(off_diag) <- 0
  off_norm <- norm(off_diag, type = "F")
  
  list(
    diag_range   = diag_range,
    diag_mean    = mean(diag_vals),
    offdiag_frob = off_norm,
    dim          = d,
    cov_matrix   = S
  )
}
