############################################################
# Riemannian geometry primitives for SPD covariance matrices
############################################################

#' Distance between two SPD covariance matrices
#'
#' Computes the distance between two symmetric positive-definite (SPD)
#' covariance matrices under one of three geometries. The affine-invariant
#' Riemannian metric is the natural distance on the SPD manifold and underlies
#' covariance-based alignment and classification.
#'
#' @param A,B Symmetric positive-definite matrices of the same dimension.
#' @param metric One of \code{"riemann"} (affine-invariant Riemannian / AIRM),
#'   \code{"logeuclid"} (log-Euclidean), or \code{"euclid"} (Frobenius).
#' @param eps Numeric stability floor for eigenvalues.
#'
#' @return A non-negative scalar distance.
#'
#' @details For the Riemannian metric the distance is
#'   \code{sqrt(sum(log(eigvals(A^{-1} B))^2))}; for log-Euclidean it is
#'   \code{||logm(A) - logm(B)||_F}; for Euclidean it is \code{||A - B||_F}.
#'
#' @references
#' Pennec, X., Fillard, P., and Ayache, N. (2006). A Riemannian framework for
#' tensor computing. International Journal of Computer Vision.
#'
#' @family riemannian
#' @seealso \code{\link{spd_mean}}, \code{\link{tangent_space}}
#'
#' @examples
#' set.seed(1)
#' A <- cov(matrix(rnorm(200 * 5), 200, 5))
#' B <- cov(matrix(rnorm(200 * 5), 200, 5))
#' riemann_distance(A, B)
#'
#' @export
riemann_distance <- function(A, B,
                             metric = c("riemann", "logeuclid", "euclid"),
                             eps = 1e-12) {
  metric <- match.arg(metric)
  if (!is.matrix(A) || !is.matrix(B) || !is.numeric(A) || !is.numeric(B)) {
    stop("A and B must be numeric matrices.")
  }
  if (nrow(A) != ncol(A) || nrow(B) != ncol(B) || nrow(A) != nrow(B)) {
    stop("A and B must be square and of the same dimension.")
  }
  if (any(!is.finite(A)) || any(!is.finite(B))) {
    stop("A and B must contain only finite values.")
  }
  A <- .symmetrize(A)
  B <- .symmetrize(B)

  if (metric == "euclid") {
    return(norm(A - B, type = "F"))
  }
  if (metric == "logeuclid") {
    return(norm(.logm_spd(A, eps = eps) - .logm_spd(B, eps = eps), type = "F"))
  }
  Ai2 <- .invsqrtm_spd(A, eps = eps)
  M <- .symmetrize(Ai2 %*% B %*% Ai2)
  ev <- pmax(eigen(M, symmetric = TRUE, only.values = TRUE)$values, eps)
  sqrt(sum(log(ev) ^ 2))
}

#' Mean of a set of SPD covariance matrices
#'
#' Computes the average of several symmetric positive-definite covariance
#' matrices under one of three geometries: the affine-invariant Riemannian
#' (geometric / Karcher) mean, the log-Euclidean mean, or the arithmetic
#' (Euclidean) mean.
#'
#' @param covs A list of \code{[p x p]} SPD matrices, or a 3D array
#'   \code{[n, p, p]}.
#' @param metric One of \code{"riemann"}, \code{"logeuclid"}, or
#'   \code{"euclid"}.
#' @param tol Convergence tolerance for the Riemannian mean.
#' @param max_iter Maximum iterations for the Riemannian mean.
#' @param eps Numeric stability floor for eigenvalues.
#'
#' @return A \code{[p x p]} SPD matrix.
#'
#' @family riemannian
#' @seealso \code{\link{riemann_distance}}, \code{\link{barycenter_whitener}},
#'   \code{\link{tangent_space}}
#'
#' @examples
#' set.seed(1)
#' covs <- lapply(1:5, function(i) cov(matrix(rnorm(200 * 4), 200, 4)))
#' M <- spd_mean(covs, metric = "riemann")
#' dim(M)
#'
#' @export
spd_mean <- function(covs,
                     metric = c("riemann", "logeuclid", "euclid"),
                     tol = 1e-6, max_iter = 50, eps = 1e-12) {
  metric <- match.arg(metric)
  cov_list <- .as_cov_list(covs)
  .check_cov_list(cov_list)
  n <- length(cov_list)

  M <- switch(
    metric,
    euclid = Reduce("+", cov_list) / n,
    logeuclid = .expm_sym(Reduce("+", lapply(cov_list, .logm_spd)) / n),
    riemann = .riemann_mean_spd(cov_list, tol = tol, max_iter = max_iter, eps = eps)
  )
  .symmetrize(M)
}

#' Fit a whitener that maps the barycenter of a set of covariances to identity
#'
#' Builds a single linear whitener \code{W = M^{-1/2}}, where \code{M} is the
#' mean covariance (\code{\link{spd_mean}}) of the supplied covariance matrices
#' under the chosen geometry. Applying \code{W} maps data whose covariance is
#' the barycenter \code{M} to the identity -- a covariance-space alternative to
#' \code{\link{euclidean_alignment}} that lets you align to the Riemannian
#' (geometric) mean of a set of session/epoch covariances.
#'
#' @param covs A list of \code{[p x p]} SPD matrices, or a 3D array
#'   \code{[n, p, p]} (e.g. from \code{\link{epoch_covariances}}).
#' @param metric Geometry of the barycenter; one of \code{"riemann"},
#'   \code{"logeuclid"}, or \code{"euclid"}.
#' @param tol,max_iter Convergence controls for the Riemannian mean.
#' @param eps Numeric stability floor for eigenvalues.
#'
#' @return An \code{ea_model} object (see \code{\link{euclidean_alignment}})
#'   with the whitening matrix \code{W} and reference covariance \code{M}; apply
#'   it with \code{predict()} and invert with \code{\link{inverse_ea}}.
#'
#' @family alignment
#' @seealso \code{\link{spd_mean}}, \code{\link{euclidean_alignment}},
#'   \code{\link{recenter}}
#'
#' @examples
#' set.seed(1)
#' X <- array(rnorm(20 * 6 * 128), dim = c(20, 6, 128))
#' covs <- epoch_covariances(X, lambda = 0.05)
#' bw <- barycenter_whitener(covs, metric = "riemann")
#' Z <- predict(bw, matrix(rnorm(100 * 6), 100, 6))
#'
#' @export
barycenter_whitener <- function(covs,
                                metric = c("riemann", "logeuclid", "euclid"),
                                tol = 1e-6, max_iter = 50, eps = 1e-12) {
  metric <- match.arg(metric)
  M <- spd_mean(covs, metric = metric, tol = tol, max_iter = max_iter, eps = eps)
  W <- .invsqrtm_spd(M, eps = eps)

  structure(
    list(
      W = W,
      reference_cov = M,
      input_type = "raw",
      mean_method = metric,
      center = FALSE,
      tol = tol,
      max_iter = max_iter,
      eps = eps
    ),
    class = "ea_model"
  )
}
