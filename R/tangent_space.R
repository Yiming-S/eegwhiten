############################################################
# Tangent space mapping for SPD covariance matrices
############################################################

# Coerce supported inputs into a list of square symmetric matrices.
# Accepts a list of [p x p] matrices or a 3D array [n, p, p].
.as_cov_list <- function(covs) {
  if (is.array(covs) && !is.matrix(covs) && length(dim(covs)) == 3L) {
    dims <- dim(covs)
    if (dims[2] != dims[3]) {
      stop("A 3D covariance array must have shape [n, p, p].")
    }
    # Force a matrix so a [n, 1, 1] array (p = 1) does not collapse to a vector.
    return(lapply(seq_len(dims[1]),
                  function(i) matrix(covs[i, , ], nrow = dims[2], ncol = dims[3])))
  }
  if (is.matrix(covs)) {
    return(list(covs))
  }
  if (is.list(covs)) {
    return(covs)
  }
  stop("covs must be a list of matrices or a 3D array [n, p, p].")
}

.check_cov_list <- function(cov_list) {
  if (length(cov_list) == 0L) stop("No covariance matrices provided.")
  p <- ncol(cov_list[[1]])
  for (C in cov_list) {
    if (!is.matrix(C) || !is.numeric(C) || any(!is.finite(C))) {
      stop("Each covariance must be a finite numeric matrix.")
    }
    if (nrow(C) != ncol(C)) stop("Each covariance must be square.")
    if (ncol(C) != p) stop("All covariances must share the same dimension.")
    if (!isTRUE(all.equal(C, t(C), tolerance = 1e-8))) {
      stop("Each covariance must be symmetric.")
    }
  }
  p
}

#' Per-epoch channel covariance matrices from an EEG tensor
#'
#' Computes one channel-by-channel covariance matrix per trial/epoch from a
#' 3D EEG tensor. These covariances are the natural input for Riemannian
#' pipelines such as \code{\link{tangent_space}} and Euclidean / Riemannian
#' recentering.
#'
#' @param X_tensor Numeric 3D array \code{[n_trials, n_channels, n_time]}.
#' @param lambda Numeric in \code{[0, 1]} or the string \code{"auto"}; optional
#'   shrinkage of each covariance toward a scaled identity
#'   (\code{(1 - lambda) C + lambda (tr(C)/p) I}). \code{"auto"} selects a
#'   per-epoch Ledoit-Wolf intensity from that epoch's time samples. Useful when
#'   the number of time samples is small relative to the number of channels.
#'   Defaults to 0 (no shrinkage). When \code{n_time <= n_channels}, covariances
#'   are rank-deficient (not positive-definite) unless \code{lambda > 0}; a
#'   warning is emitted in that case.
#' @param center Logical; whether to subtract each channel's mean over time
#'   before computing the covariance.
#' @param shrinkage Deprecated; use \code{lambda}.
#'
#' @return A list of \code{[n_channels x n_channels]} covariance matrices,
#'   one per trial. They are positive-definite (suitable as
#'   \code{\link{tangent_space}} input) only when \code{n_time > n_channels} or
#'   \code{lambda > 0}.
#'
#' @family riemannian
#' @seealso \code{\link{tangent_space}}, \code{\link{recenter}}
#'
#' @examples
#' set.seed(1)
#' X <- array(rnorm(20 * 8 * 64), dim = c(20, 8, 64))
#' covs <- epoch_covariances(X, lambda = 0.05)
#' length(covs)
#'
#' @export
epoch_covariances <- function(X_tensor, lambda = 0, center = TRUE, shrinkage = NULL) {
  if (!is.null(shrinkage)) {
    warning("'shrinkage' is deprecated; use 'lambda'.", call. = FALSE)
    lambda <- shrinkage
  }
  if (!is.array(X_tensor) || length(dim(X_tensor)) != 3L) {
    stop("X_tensor must be a 3D numeric array with shape [n_trials, n_channels, n_time].")
  }
  if (!is.numeric(X_tensor) || any(!is.finite(X_tensor))) {
    stop("X_tensor must contain only finite numeric values.")
  }
  .check_unit_interval(lambda, "lambda", allow_char_auto = TRUE)
  lambda_auto <- is.character(lambda)
  if (!is.logical(center) || length(center) != 1L || is.na(center)) {
    stop("center must be TRUE or FALSE.")
  }

  dims <- dim(X_tensor)
  n_trials <- dims[1]
  n_channels <- dims[2]
  n_time <- dims[3]
  if (n_time < 2L) stop("Each trial must have at least 2 time samples.")
  if (!lambda_auto && lambda == 0 && n_time <= n_channels) {
    warning("n_time <= n_channels with lambda = 0: per-epoch covariances are rank-deficient (not positive-definite). Set lambda > 0 (or \"auto\") for a stable tangent-space mapping.")
  }

  lapply(seq_len(n_trials), function(i) {
    # Trial slice is [n_channels x n_time].
    Xi <- matrix(X_tensor[i, , ], nrow = n_channels, ncol = n_time)
    if (center) Xi <- Xi - rowMeans(Xi)
    C <- .symmetrize(tcrossprod(Xi) / (n_time - 1L))
    # "auto" picks a per-epoch Ledoit-Wolf intensity from the time samples
    # (observations) of this epoch (channels as features).
    lam <- if (lambda_auto) .lambda_lw(t(Xi)) else lambda
    if (lam > 0) C <- .shrink_cov(C, lam, target = "identity")
    C
  })
}

#' Map SPD covariance matrices to the tangent space
#'
#' Projects a set of symmetric positive-definite (SPD) covariance matrices to
#' the tangent space at a reference point and vectorizes them. This linearizes
#' the curved SPD manifold so that ordinary Euclidean classifiers and
#' regressors can be applied to the resulting feature vectors -- the standard
#' Riemannian pipeline for EEG/BCI.
#'
#' For each covariance \code{C}, the map computes
#' \code{L = logm(ref^{-1/2} C ref^{-1/2})} and returns the upper triangle of
#' \code{L} with off-diagonal entries scaled by \code{sqrt(2)}, so that the
#' Euclidean norm of the feature vector equals the Riemannian distance from the
#' reference.
#'
#' @param covs A list of \code{[p x p]} SPD matrices, or a 3D array
#'   \code{[n, p, p]}.
#' @param ref Optional reference SPD matrix \code{[p x p]}. If \code{NULL}, the
#'   reference is the mean covariance computed with \code{mean_method}.
#' @param mean_method Method for the reference mean when \code{ref = NULL}; one
#'   of \code{"riemann"} (affine-invariant), \code{"logeuclid"}, or
#'   \code{"euclid"}.
#' @param tol Tolerance for \code{"riemann"} mean convergence.
#' @param max_iter Maximum iterations for \code{"riemann"} mean.
#' @param eps Numeric stability floor for eigenvalues.
#'
#' @return A numeric matrix \code{[n x p(p+1)/2]} of tangent-space features.
#'   The reference and metadata are attached as attributes \code{"reference"},
#'   \code{"p"}, and \code{"mean_method"}.
#'
#' @references
#' Barachant, A., Bonnet, S., Congedo, M., and Jutten, C. (2012). Multiclass
#' brain-computer interface classification by Riemannian geometry. IEEE
#' Transactions on Biomedical Engineering.
#'
#' @family riemannian
#' @seealso \code{\link{untangent_space}}, \code{\link{epoch_covariances}}
#'
#' @examples
#' set.seed(1)
#' X <- array(rnorm(30 * 6 * 64), dim = c(30, 6, 64))
#' covs <- epoch_covariances(X, lambda = 0.05)
#' V <- tangent_space(covs, mean_method = "logeuclid")
#' dim(V)
#'
#' @export
tangent_space <- function(covs,
                          ref = NULL,
                          mean_method = c("riemann", "logeuclid", "euclid"),
                          tol = 1e-6,
                          max_iter = 50,
                          eps = 1e-12) {
  mean_method <- match.arg(mean_method)
  cov_list <- .as_cov_list(covs)
  p <- .check_cov_list(cov_list)

  if (is.null(ref)) {
    ref <- switch(
      mean_method,
      euclid = Reduce("+", cov_list) / length(cov_list),
      logeuclid = .expm_sym(Reduce("+", lapply(cov_list, .logm_spd)) / length(cov_list)),
      riemann = .riemann_mean_spd(cov_list, tol = tol, max_iter = max_iter, eps = eps)
    )
    ref <- .symmetrize(ref)
  } else {
    if (!is.matrix(ref) || !is.numeric(ref) || nrow(ref) != p || ncol(ref) != p) {
      stop(sprintf("ref must be a %d x %d numeric matrix.", p, p))
    }
    .check_symmetric_pd(ref, require_pd = TRUE, return_eigen = FALSE)
    ref <- .symmetrize(ref)
  }

  ref_inv_half <- .invsqrtm_spd(ref, eps = eps)

  len_out <- p * (p + 1) / 2
  mapped <- vapply(cov_list, function(C) {
    M <- ref_inv_half %*% C %*% ref_inv_half
    .spd_tangent_vec(.logm_spd(.symmetrize(M), eps = eps))
  }, numeric(len_out))
  # vapply returns a plain vector when len_out == 1 (p == 1); force [n x len_out].
  V <- if (len_out == 1L) matrix(mapped, ncol = 1L) else t(mapped)

  attr(V, "reference") <- ref
  attr(V, "p") <- p
  attr(V, "mean_method") <- mean_method
  V
}

#' Map tangent-space vectors back to SPD covariance matrices
#'
#' Inverse of \code{\link{tangent_space}}: reconstructs SPD matrices from their
#' tangent-space feature vectors and a reference point.
#'
#' @param V A numeric matrix \code{[n x p(p+1)/2]} of tangent-space features, or
#'   a single numeric vector. If \code{ref} is \code{NULL}, the reference stored
#'   in \code{attr(V, "reference")} (as returned by \code{tangent_space}) is
#'   used.
#' @param ref Optional reference SPD matrix \code{[p x p]}.
#'
#' @return A list of \code{[p x p]} SPD matrices.
#'
#' @family riemannian
#' @seealso \code{\link{tangent_space}}
#'
#' @examples
#' set.seed(1)
#' X <- array(rnorm(10 * 5 * 64), dim = c(10, 5, 64))
#' covs <- epoch_covariances(X)
#' V <- tangent_space(covs, mean_method = "logeuclid")
#' covs_rec <- untangent_space(V)
#' max(abs(covs_rec[[1]] - covs[[1]]))
#'
#' @export
untangent_space <- function(V, ref = NULL) {
  if (is.null(ref)) {
    ref <- attr(V, "reference")
    if (is.null(ref)) {
      stop("ref is NULL and V has no 'reference' attribute. Provide ref.")
    }
  }
  if (!is.matrix(ref) || !is.numeric(ref) || nrow(ref) != ncol(ref)) {
    stop("ref must be a square numeric matrix.")
  }
  p <- ncol(ref)
  expected_len <- p * (p + 1) / 2

  if (is.vector(V) && !is.list(V)) {
    V <- matrix(V, nrow = 1L)
  }
  if (!is.matrix(V) || !is.numeric(V)) {
    stop("V must be a numeric matrix or vector.")
  }
  if (ncol(V) != expected_len) {
    stop(sprintf("V must have %d columns for a %d x %d reference.", expected_len, p, p))
  }

  ref_half <- .sqrtm_spd(.symmetrize(ref))

  lapply(seq_len(nrow(V)), function(i) {
    L <- .spd_tangent_unvec(V[i, ], p)
    .symmetrize(ref_half %*% .expm_sym(L) %*% ref_half)
  })
}
