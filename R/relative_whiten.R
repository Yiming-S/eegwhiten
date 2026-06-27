############################################################
# Relative (generalized) whitening of one covariance by another
############################################################

#' Relative whitening via a generalized eigenvalue decomposition
#'
#' Computes a linear transform that simultaneously whitens a reference
#' covariance to the identity and diagonalizes a target covariance. This is the
#' generalized eigenvalue problem \code{Sigma w = gamma * reference w} that
#' underlies spatial filtering methods such as CSP and xDAWN: the resulting
#' filters maximize the ratio of \code{Sigma} to \code{reference} variance.
#'
#' Given the reference inverse square root \code{R = reference^{-1/2}} and the
#' eigendecomposition \code{R Sigma R = V Gamma V^T}, the filters are
#' \code{W = V^T R}. They satisfy \code{W reference W^T = I} (reference is
#' whitened) and \code{W Sigma W^T = Gamma} (target is diagonal), with the
#' diagonal \code{Gamma} giving the per-component variance ratio.
#'
#' @param Sigma Target symmetric positive-definite covariance \code{[p x p]}
#'   (for example a signal- or class-conditional covariance).
#' @param reference Reference symmetric positive-definite covariance
#'   \code{[p x p]} (for example a noise or baseline covariance) that is mapped
#'   to the identity.
#' @param n_comp Optional integer; number of leading components (largest
#'   variance ratios) to keep. If \code{NULL}, all \code{p} are returned.
#' @param eps Numeric stability floor for eigenvalues.
#'
#' @return A list with:
#'   \item{W}{Filter matrix \code{[n_comp x p]}; apply as \code{Z = X \%*\% t(W)}.}
#'   \item{values}{Generalized eigenvalues (variance ratios), descending.}
#'   \item{patterns}{Forward patterns \code{[p x n_comp]} (columns), the spatial
#'     activity pattern of each component.}
#'   \item{reference_inv_sqrt}{The reference inverse square root \code{R}.}
#'
#' @references
#' Blankertz, B., Tomioka, R., Lemm, S., Kawanabe, M., and Mueller, K.-R.
#' (2008). Optimizing spatial filters for robust EEG single-trial analysis.
#' IEEE Signal Processing Magazine.
#'
#' @family alignment
#' @seealso \code{\link{whiten_model}}, \code{\link{euclidean_alignment}}
#'
#' @examples
#' set.seed(1)
#' p <- 6
#' Xa <- matrix(rnorm(200 * p), 200, p)
#' Xb <- matrix(rnorm(200 * p), 200, p) %*% diag(seq_len(p))
#' res <- whiten_relative(cov(Xb), cov(Xa))
#' round(res$values, 3)
#'
#' @export
whiten_relative <- function(Sigma, reference, n_comp = NULL, eps = 1e-10) {
  .check_symmetric_pd(Sigma, require_pd = TRUE, return_eigen = FALSE)
  .check_symmetric_pd(reference, require_pd = TRUE, return_eigen = FALSE)
  p <- ncol(Sigma)
  if (nrow(reference) != p || ncol(reference) != p) {
    stop("Sigma and reference must have the same dimensions.")
  }
  if (!is.null(n_comp)) {
    if (!is.numeric(n_comp) || length(n_comp) != 1L || !is.finite(n_comp) ||
        n_comp != as.integer(n_comp) || n_comp < 1L || n_comp > p) {
      stop("n_comp must be an integer between 1 and ncol(Sigma).")
    }
    n_comp <- as.integer(n_comp)
  }
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("eps must be a positive finite scalar.")
  }

  R <- .invsqrtm_spd(.symmetrize(reference), eps = eps)
  M <- .symmetrize(R %*% .symmetrize(Sigma) %*% R)
  ed <- eigen(M, symmetric = TRUE)
  V <- ed$vectors
  gamma <- ed$values

  k <- if (is.null(n_comp)) p else n_comp
  Vk <- V[, seq_len(k), drop = FALSE]
  gk <- gamma[seq_len(k)]

  # Filters: rows are components, columns are channels.
  W <- crossprod(Vk, R)
  ch_names <- colnames(Sigma)
  comp_names <- paste0("GED", seq_len(k))
  W <- .set_matrix_attr(W, row_names = comp_names, col_names = ch_names, method = "relative")

  # Forward patterns are columns of the (pseudo)inverse of the filter map.
  ref_half <- .sqrtm_spd(.symmetrize(reference), eps = eps)
  patterns <- ref_half %*% Vk
  if (!is.null(ch_names)) rownames(patterns) <- ch_names
  colnames(patterns) <- comp_names

  structure(
    list(
      W = W,
      values = gk,
      patterns = patterns,
      reference_inv_sqrt = R
    ),
    class = "relative_whiten"
  )
}

#' Print a relative-whitening model
#'
#' @param x A \code{relative_whiten} object.
#' @param ... Unused.
#' @return The input object invisibly.
#' @export
#' @method print relative_whiten
print.relative_whiten <- function(x, ...) {
  cat("<relative_whiten>\n")
  cat("  components :", nrow(x$W), "\n")
  cat("  channels   :", ncol(x$W), "\n")
  n_show <- min(6L, length(x$values))
  cat("  variance ratios:", paste(signif(x$values[seq_len(n_show)], 3), collapse = ", "),
      if (length(x$values) > 6L) "..." else "", "\n")
  invisible(x)
}

#' Apply a relative-whitening model to new data
#'
#' @param object A \code{relative_whiten} object.
#' @param newdata Numeric matrix with the same number of channels.
#' @param ... Unused.
#' @return Filtered data matrix \code{newdata \%*\% t(W)}.
#' @export
#' @method predict relative_whiten
predict.relative_whiten <- function(object, newdata, ...) {
  if (!inherits(object, "relative_whiten")) stop("object must be a 'relative_whiten' object.")
  if (!is.matrix(newdata) || !is.numeric(newdata) || any(!is.finite(newdata))) {
    stop("newdata must be a finite numeric matrix.")
  }
  if (ncol(newdata) != ncol(object$W)) {
    stop(sprintf("Mismatch: model expects %d channels, but newdata has %d.", ncol(object$W), ncol(newdata)))
  }
  newdata %*% t(object$W)
}
