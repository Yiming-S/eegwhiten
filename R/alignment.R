############################################################
# Euclidean / Riemannian alignment for batch EEG matrices
############################################################

#' Euclidean Alignment for EEG Matrices
#'
#' Aligns a list of raw EEG matrices or covariance matrices to a global
#' reference covariance. This avoids per-matrix independent whitening,
#' which can distort between-session/class covariance structure.
#'
#' @param X_list List of numeric matrices. Each element is either:
#'   raw data \code{[n_samples x n_channels]} or
#'   covariance \code{[n_channels x n_channels]}.
#' @param input Input interpretation:
#'   \code{"auto"}, \code{"raw"}, or \code{"cov"}.
#' @param mean_method Reference covariance mean:
#'   \code{"riemann"}, \code{"logeuclid"}, or \code{"euclid"}.
#' @param center Logical; whether to center each raw matrix before covariance
#'   computation when \code{input = "raw"}.
#' @param tol Tolerance for \code{"riemann"} mean convergence.
#' @param max_iter Maximum iterations for \code{"riemann"} mean.
#' @param eps Numeric stability floor for eigenvalues.
#'
#' @return A list with:
#'   \item{aligned}{Aligned matrices.}
#'   \item{W}{Global alignment matrix \code{C_ref^{-1/2}}.}
#'   \item{reference_cov}{Reference covariance matrix \code{C_ref}.}
#'   \item{covariances}{Covariances used to estimate \code{C_ref}.}
#'   \item{input_type}{Resolved input type (\code{"raw"} or \code{"cov"}).}
#'   \item{mean_method}{Mean method used.}
#'   \item{model}{Alignment model object (class \code{"ea_model"}).}
#'
#' @references
#' He, H., and Wu, D. (2019). Transfer learning for brain-computer interfaces:
#' A Euclidean space data alignment approach. IEEE Transactions on Biomedical
#' Engineering.
#'
#' @export
euclidean_alignment <- function(X_list,
                                input = c("auto", "raw", "cov"),
                                mean_method = c("logeuclid", "riemann", "euclid"),
                                center = TRUE,
                                tol = 1e-6,
                                max_iter = 50,
                                eps = 1e-12) {
  input <- match.arg(input)
  mean_method <- match.arg(mean_method)

  if (!is.list(X_list) || length(X_list) == 0L) {
    stop("X_list must be a non-empty list of matrices.")
  }
  if (!is.logical(center) || length(center) != 1L || is.na(center)) {
    stop("center must be TRUE or FALSE.")
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("tol must be a positive finite scalar.")
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1L || !is.finite(max_iter) || max_iter < 1) {
    stop("max_iter must be a positive integer.")
  }
  max_iter <- as.integer(max_iter)

  ok_mat <- vapply(X_list, function(x) is.matrix(x) && is.numeric(x) && all(is.finite(x)), logical(1))
  if (!all(ok_mat)) stop("All elements in X_list must be finite numeric matrices.")

  # Auto-detect whether inputs look like covariance matrices
  if (input == "auto") {
    all_square <- all(vapply(X_list, function(x) nrow(x) == ncol(x), logical(1)))
    all_sym <- all(vapply(X_list, function(x) isTRUE(all.equal(x, t(x), tolerance = 1e-8)), logical(1)))
    input_type <- if (all_square && all_sym) "cov" else "raw"
  } else {
    input_type <- input
  }

  if (input_type == "raw") {
    p <- ncol(X_list[[1]])
    same_p <- all(vapply(X_list, function(x) ncol(x) == p, logical(1)))
    if (!same_p) stop("All raw matrices must have the same number of columns (channels).")

    covariances <- lapply(X_list, function(x) {
      x_use <- if (center) sweep(x, 2, colMeans(x), "-") else x
      stats::cov(x_use)
    })
  } else {
    p <- ncol(X_list[[1]])
    same_dim <- all(vapply(X_list, function(x) nrow(x) == p && ncol(x) == p, logical(1)))
    if (!same_dim) stop("All covariance matrices must share the same dimensions.")

    covariances <- X_list
    invisible(lapply(covariances, .check_symmetric_pd, require_pd = TRUE))
  }

  # Reference covariance
  if (mean_method == "euclid") {
    C_ref <- Reduce("+", covariances) / length(covariances)
  } else if (mean_method == "logeuclid") {
    C_ref <- .expm_sym(Reduce("+", lapply(covariances, .logm_spd)) / length(covariances))
  } else {
    C_ref <- .riemann_mean_spd(covariances, tol = tol, max_iter = max_iter, eps = eps)
  }
  C_ref <- .symmetrize(C_ref)

  # Global alignment matrix
  W_ref <- .invsqrtm_spd(C_ref, eps = eps)

  aligned <- if (input_type == "raw") {
    lapply(X_list, function(x) {
      x_use <- if (center) sweep(x, 2, colMeans(x), "-") else x
      x_use %*% W_ref
    })
  } else {
    lapply(covariances, function(C) t(W_ref) %*% C %*% W_ref)
  }

  model <- structure(
    list(
      W = W_ref,
      reference_cov = C_ref,
      input_type = input_type,
      mean_method = mean_method,
      center = center,
      tol = tol,
      max_iter = max_iter,
      eps = eps
    ),
    class = "ea_model"
  )

  list(
    aligned = aligned,
    W = W_ref,
    reference_cov = C_ref,
    covariances = covariances,
    input_type = input_type,
    mean_method = mean_method,
    model = model
  )
}

#' Apply an EA Model to New Data
#'
#' @param object An \code{ea_model} object returned by
#'   \code{euclidean_alignment()}.
#' @param newdata A numeric matrix or a list of numeric matrices.
#' @param ... Unused.
#'
#' @return Aligned matrix or list of aligned matrices.
#' @export
#' @method predict ea_model
predict.ea_model <- function(object, newdata, ...) {
  if (!inherits(object, "ea_model")) {
    stop("object must be an 'ea_model' object.")
  }

  apply_one <- function(x) {
    if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x))) {
      stop("newdata must contain only finite numeric matrices.")
    }

    p <- ncol(object$W)
    if (identical(object$input_type, "raw")) {
      if (ncol(x) != p) {
        stop(sprintf("Raw matrix must have %d columns.", p))
      }
      x_use <- if (isTRUE(object$center)) sweep(x, 2, colMeans(x), "-") else x
      return(x_use %*% object$W)
    }

    if (nrow(x) != p || ncol(x) != p) {
      stop(sprintf("Covariance matrix must be %d x %d.", p, p))
    }
    if (!isTRUE(all.equal(x, t(x), tolerance = 1e-8))) {
      stop("Covariance matrix must be symmetric.")
    }
    t(object$W) %*% x %*% object$W
  }

  if (is.list(newdata)) {
    return(lapply(newdata, apply_one))
  }
  apply_one(newdata)
}

#' Inverse Transform for EA-Aligned Data
#'
#' @param Z An aligned matrix (or list of aligned matrices).
#' @param model An \code{ea_model} object.
#' @param center Optional center vector added back for raw-data inverse mode.
#'
#' @return Matrix (or list of matrices) mapped back from EA-aligned space.
#' @export
inverse_ea <- function(Z, model, center = NULL) {
  if (!inherits(model, "ea_model")) {
    stop("model must be an 'ea_model' object.")
  }

  A <- .sqrtm_spd(model$reference_cov, eps = model$eps)

  apply_one <- function(z) {
    if (!is.matrix(z) || !is.numeric(z) || any(!is.finite(z))) {
      stop("Z must contain only finite numeric matrices.")
    }

    p <- ncol(model$W)
    if (identical(model$input_type, "raw")) {
      if (ncol(z) != p) {
        stop(sprintf("Aligned raw matrix must have %d columns.", p))
      }
      x <- z %*% A
      if (!is.null(center)) {
        if (!is.numeric(center) || length(center) != p || any(!is.finite(center))) {
          stop("center must be a finite numeric vector with length ncol(model$W).")
        }
        x <- sweep(x, 2, center, "+")
      }
      return(x)
    }

    if (nrow(z) != p || ncol(z) != p) {
      stop(sprintf("Aligned covariance must be %d x %d.", p, p))
    }
    .symmetrize(t(A) %*% z %*% A)
  }

  if (is.list(Z)) {
    return(lapply(Z, apply_one))
  }
  apply_one(Z)
}
