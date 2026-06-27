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

  list(
    W = W,
    values = gk,
    patterns = patterns,
    reference_inv_sqrt = R
  )
}
