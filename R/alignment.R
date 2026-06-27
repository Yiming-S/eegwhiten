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
#' @family alignment
#' @seealso \code{\link{predict.ea_model}}, \code{\link{inverse_ea}},
#'   \code{\link{whiten_batch}} with \code{mode = "ea"}.
#'
#' @references
#' He, H., and Wu, D. (2019). Transfer learning for brain-computer interfaces:
#' A Euclidean space data alignment approach. IEEE Transactions on Biomedical
#' Engineering.
#'
#' @examples
#' set.seed(1)
#' X1 <- matrix(rnorm(200 * 6), 200, 6)
#' X2 <- matrix(rnorm(180 * 6), 180, 6)
#' ea <- euclidean_alignment(list(X1, X2), input = "raw",
#'                           mean_method = "logeuclid")
#' str(ea$aligned)
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

#' Recenter a single EEG matrix to the identity
#'
#' Whitens one data matrix by the inverse square root of its own covariance
#' (or a supplied reference covariance), mapping that covariance to the
#' identity. Applied independently to each session/subject, this is the
#' per-domain recentering at the heart of Euclidean Alignment (He and Wu,
#' 2019) and Riemannian recentering (Zanini et al., 2018): every domain's mean
#' covariance is moved to a common point, which removes the bulk of
#' between-session / between-subject covariance shift before classification.
#'
#' Unlike \code{\link{euclidean_alignment}}, which applies a single shared
#' transform built from the grand-mean covariance, \code{recenter} aligns each
#' matrix by its own covariance, so each domain is individually whitened to the
#' identity.
#'
#' @param X Numeric matrix \code{[n_samples x n_channels]}.
#' @param reference Optional reference covariance \code{[p x p]}. If
#'   \code{NULL} (default), the covariance of \code{X} itself is used, mapping
#'   \code{X} to the identity.
#' @param center Logical; whether to subtract the column means before computing
#'   the covariance and before whitening.
#' @param lambda Numeric in \code{[0, 1]}; optional shrinkage of the covariance
#'   toward a scaled identity for stability with short epochs.
#' @param eps Numeric stability floor for eigenvalues.
#'
#' @return A list with:
#'   \item{Z}{Recentered data matrix.}
#'   \item{W}{Alignment matrix \code{reference^{-1/2}}.}
#'   \item{reference_cov}{Covariance that was mapped to the identity.}
#'   \item{center}{Column means subtracted (zeros if \code{center = FALSE}).}
#'
#' @references
#' He, H., and Wu, D. (2019). Transfer learning for brain-computer interfaces:
#' A Euclidean space data alignment approach. IEEE Transactions on Biomedical
#' Engineering.
#'
#' @family alignment
#' @seealso \code{\link{euclidean_alignment}}, \code{\link{whiten_batch}} with
#'   \code{mode = "recenter"}.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(200 * 6), 200, 6)
#' rc <- recenter(X)
#' diag(round(cov(rc$Z), 6))
#'
#' @export
recenter <- function(X, reference = NULL, center = TRUE, lambda = 0, eps = 1e-10) {
  if (!is.matrix(X) || !is.numeric(X) || any(!is.finite(X))) {
    stop("X must be a finite numeric matrix.")
  }
  if (!is.logical(center) || length(center) != 1L || is.na(center)) {
    stop("center must be TRUE or FALSE.")
  }
  .check_unit_interval(lambda, "lambda")

  p <- ncol(X)
  mu <- if (center) colMeans(X) else rep(0, p)
  Xc <- if (center) .center_cols(X, mu) else X

  if (is.null(reference)) {
    if (nrow(X) < 2L) {
      stop("X must have at least 2 rows to estimate a covariance, or pass an explicit reference.")
    }
    C <- .symmetrize(stats::cov(Xc))
    if (lambda > 0) C <- .shrink_cov(C, lambda, target = "identity")
  } else {
    if (!is.matrix(reference) || !is.numeric(reference) ||
        nrow(reference) != p || ncol(reference) != p) {
      stop(sprintf("reference must be a %d x %d numeric matrix.", p, p))
    }
    .check_symmetric_pd(reference, require_pd = TRUE, return_eigen = FALSE)
    C <- .symmetrize(reference)
    if (lambda > 0) C <- .shrink_cov(C, lambda, target = "identity")
  }

  # Guard against rank-deficient covariance (e.g. short epochs): a silent
  # eigenvalue floor would return non-whitened data, breaking the contract.
  ed <- eigen(C, symmetric = TRUE)
  if (min(ed$values) <= eps) {
    stop(sprintf(
      "recenter(): covariance is not positive-definite (min eigenvalue = %g). Increase lambda for short or rank-deficient epochs.",
      min(ed$values)))
  }
  W <- ed$vectors %*% (t(ed$vectors) * (1 / sqrt(ed$values)))
  Z <- Xc %*% W
  if (!is.null(colnames(X))) colnames(Z) <- colnames(X)

  structure(
    list(Z = Z, W = W, reference_cov = C, center = mu),
    class = "recenter"
  )
}

#' Print a recentering model
#'
#' @param x A \code{recenter} object.
#' @param ... Unused.
#' @return The input object invisibly.
#' @export
#' @method print recenter
print.recenter <- function(x, ...) {
  cat("<recenter>\n")
  cat("  channels   :", ncol(x$W), "\n")
  cat("  samples    :", nrow(x$Z), "\n")
  cat("  centered   :", any(x$center != 0), "\n")
  invisible(x)
}

#' Apply a recentering model to new data
#'
#' @param object A \code{recenter} object.
#' @param newdata Numeric matrix with the same number of channels.
#' @param ... Unused.
#' @return Recentered data matrix.
#' @export
#' @method predict recenter
predict.recenter <- function(object, newdata, ...) {
  if (!inherits(object, "recenter")) stop("object must be a 'recenter' object.")
  if (!is.matrix(newdata) || !is.numeric(newdata) || any(!is.finite(newdata))) {
    stop("newdata must be a finite numeric matrix.")
  }
  p <- ncol(object$W)
  if (ncol(newdata) != p) {
    stop(sprintf("Mismatch: model expects %d channels, but newdata has %d.", p, ncol(newdata)))
  }
  Xc <- .center_cols(newdata, object$center)
  Z <- Xc %*% object$W
  if (!is.null(colnames(newdata))) colnames(Z) <- colnames(newdata)
  Z
}

#' Print an EA model
#'
#' @param x An \code{ea_model} object.
#' @param ... Unused.
#' @return The input object invisibly.
#' @export
#' @method print ea_model
print.ea_model <- function(x, ...) {
  cat("<ea_model>\n")
  cat("  channels   :", ncol(x$W), "\n")
  cat("  input_type :", x$input_type, "\n")
  cat("  mean_method:", x$mean_method, "\n")
  cat("  centered   :", isTRUE(x$center), "\n")
  invisible(x)
}

#' Summarize an EA model
#'
#' @param object An \code{ea_model} object.
#' @param ... Unused.
#' @return The input object invisibly (prints reference-covariance conditioning).
#' @export
#' @method summary ea_model
summary.ea_model <- function(object, ...) {
  print(object)
  ev <- eigen(object$reference_cov, symmetric = TRUE, only.values = TRUE)$values
  ev <- ev[ev > 0]
  cond <- if (length(ev)) max(ev) / min(ev) else Inf
  cat("  ref cond   :", signif(cond, 4), "\n")
  invisible(object)
}

#' Apply an EA Model to New Data
#'
#' @param object An \code{ea_model} object returned by
#'   \code{euclidean_alignment()}.
#' @param newdata A numeric matrix or a list of numeric matrices.
#' @param ... Unused.
#'
#' @return Aligned matrix or list of aligned matrices.
#'
#' @seealso \code{\link{euclidean_alignment}}, \code{\link{inverse_ea}}
#'
#' @examples
#' set.seed(1)
#' X1 <- matrix(rnorm(200 * 5), 200, 5)
#' X2 <- matrix(rnorm(180 * 5), 180, 5)
#' ea <- euclidean_alignment(list(X1, X2), input = "raw")
#' Z_new <- predict(ea$model, matrix(rnorm(100 * 5), 100, 5))
#'
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
#' @param model An \code{ea_model} object.
#' @param Z An aligned matrix (or list of aligned matrices).
#' @param center Optional center vector added back for raw-data inverse mode.
#'
#' @return Matrix (or list of matrices) mapped back from EA-aligned space.
#'
#' @details The deprecated argument order \code{inverse_ea(Z, model)} is still
#'   accepted (with a warning) for backward compatibility.
#'
#' @family inverse transforms
#' @seealso \code{\link{euclidean_alignment}}, \code{\link{predict.ea_model}}
#'
#' @examples
#' set.seed(1)
#' X1 <- matrix(rnorm(200 * 5), 200, 5)
#' X2 <- matrix(rnorm(180 * 5), 180, 5)
#' ea <- euclidean_alignment(list(X1, X2), input = "raw", center = TRUE)
#' Z1 <- predict(ea$model, X1)
#' X1_rec <- inverse_ea(ea$model, Z1)
#'
#' @export
inverse_ea <- function(model, Z, center = NULL) {
  io <- .order_model_first(model, Z, "ea_model", "inverse_ea")
  model <- io$model
  Z <- io$data

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
