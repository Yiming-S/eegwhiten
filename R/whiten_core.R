############################################################
# Core whitening API for EEG data
# Rows = trials/epochs, columns = channels or features
############################################################

#' Fit a whitening matrix from a covariance matrix
#'
#' Compute a whitening matrix \code{W} from a symmetric
#' positive-definite covariance matrix using one of the
#' implemented methods (PCA, ZCA, SVD, Cholesky, etc.).
#'
#' @param Sigma Symmetric positive-definite covariance matrix.
#' @param method Whitening method; one of \code{"SVD"}, \code{"ZCA"},
#'   \code{"ZCA-cor"}, \code{"PCA"}, \code{"PCA-cor"}, \code{"Cholesky"}.
#'
#' @return Whitening matrix \code{W}.
#'
#' @export
whiten_fit <- function(Sigma,
                       method = c("SVD", "ZCA", "ZCA-cor",
                                  "PCA", "PCA-cor", "Cholesky")) {
  method <- match.arg(method)
  
  # All implemented methods expect symmetric PD covariance
  .check_symmetric_pd(Sigma, require_pd = TRUE)
  
  W <- switch(
    method,
    "SVD"      = SVD(Sigma)$W,
    "ZCA"      = ZCA(Sigma)$W,
    "ZCA-cor"  = ZCA_cor(Sigma)$W,
    "PCA"      = PCA(Sigma)$W,
    "PCA-cor"  = PCA_cor(Sigma)$W,
    "Cholesky" = Cholesky(Sigma)$W
  )
  
  if (is.null(W)) stop("Whitening matrix W is NULL.")
  W
}

#' Fit a whitening model for EEG data
#'
#' The input matrix \code{X} is assumed to have trials/epochs in rows
#' and channels or features in columns. The function estimates
#' a centering vector and a whitening matrix based on the empirical
#' covariance of the centered data.
#'
#' @param X Numeric matrix [n_trials x n_channels].
#' @param center Logical; whether to subtract column means before
#'   computing the covariance.
#' @param method Whitening method; one of \code{"SVD"}, \code{"ZCA"},
#'   \code{"ZCA-cor"}, \code{"PCA"}, \code{"PCA-cor"}, \code{"Cholesky"}.
#'
#' @return An object of class \code{"whiten_model"} containing:
#'   \item{W}{Whitening matrix.}
#'   \item{center}{Centering vector.}
#'   \item{method}{Method name.}
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(200 * 16), 200, 16)
#' colnames(X) <- paste0("Ch", seq_len(16))
#' wm <- whiten_model(X, center = TRUE, method = "ZCA")
#' Xw <- predict(wm, X)
#'
#' @export
whiten_model <- function(X, center = TRUE,
                         method = c("SVD", "ZCA", "ZCA-cor",
                                    "PCA", "PCA-cor", "Cholesky")) {
  method <- match.arg(method)
  
  if (!is.matrix(X)) stop("X must be a matrix.")
  n <- nrow(X)
  d <- ncol(X)
  
  if (d <= 1L) {
    stop("X must have more than one column.")
  }
  
  # Ensure column names exist
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("V", seq_len(d))
  }
  
  mu <- if (center) colMeans(X) else rep(0, d)
  names(mu) <- colnames(X)
  
  # Centered data
  Xc <- sweep(X, 2, mu, "-")
  
  # Empirical covariance of centered data
  S <- crossprod(Xc) / (n - 1L)
  
  W <- whiten_fit(S, method = method)
  
  structure(
    list(
      W      = W,
      center = mu,
      method = method
    ),
    class = "whiten_model"
  )
}

#' Apply a whitening model to new EEG data
#'
#' @param object A \code{whiten_model} object.
#' @param newdata Numeric matrix with the same number of columns
#'   as the data used to fit the model.
#' @param ... Unused, for S3 generic compatibility.
#'
#' @return Whitened data matrix.
#'
#' @export
#' @method predict whiten_model
predict.whiten_model <- function(object, newdata, ...) {
  if (!is.matrix(newdata)) stop("newdata must be a matrix.")
  
  d <- length(object$center)
  if (ncol(newdata) != d) {
    stop("Number of columns in newdata does not match the model.")
  }
  
  # Center using stored mean
  Xc <- sweep(newdata, 2, object$center, "-")
  
  # Z = centered X times whitening matrix
  Z <- tcrossprod(Xc, object$W)  # Xc %*% t(W)
  
  colnames(Z) <- paste0("L", seq_len(ncol(newdata)))
  attr(Z, "method") <- object$method
  
  Z
}

#' Inverse transform from whitened space back to original EEG space
#'
#' Given whitened data \code{Z} and a fitted \code{whiten_model},
#' approximately reconstruct the original-space features. This can
#' be used for interpretation or sanity checks.
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
  
  # For Z = Xc %*% t(W), with W square and invertible,
  # Xc = Z %*% solve(t(W))
  Xc <- Z %*% solve(t(model$W))
  X  <- sweep(Xc, 2, model$center, "+")
  
  colnames(X) <- names(model$center)
  X
}

#' Diagnostic check for whitening quality
#'
#' Compute basic diagnostics on the covariance of whitened data,
#' such as the range of diagonal elements and the Frobenius norm
#' of off-diagonal elements.
#'
#' @param Z Whitened data matrix.
#'
#' @return A list containing:
#'   \item{diag_range}{Range of diagonal elements of \code{cov(Z)}.}
#'   \item{diag_mean}{Mean of diagonal elements of \code{cov(Z)}.}
#'   \item{offdiag_frob}{Frobenius norm of off-diagonal entries of \code{cov(Z)}.}
#'   \item{dim}{Number of columns.}
#'   \item{cov_matrix}{The covariance matrix \code{cov(Z)}.}
#'
#' @export
check_whitening <- function(Z) {
  if (!is.matrix(Z)) stop("Z must be a matrix.")
  
  S <- cov(Z)
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