############################################################
# Cholesky-based whitening
############################################################

#' Cholesky-based whitening on a covariance matrix
#'
#' Compute a whitening matrix and factor loadings based on the
#' Cholesky factorization of a symmetric positive-definite
#' covariance matrix. This is useful for EEG data where the
#' covariance is estimated across channels or features.
#'
#' @param Sigma Symmetric positive-definite covariance matrix.
#' @param returnW Logical; if TRUE, return the whitening matrix \code{W}.
#' @param PhiPsi Logical; if TRUE, return factor loadings \code{Phi}
#'   and standardized loadings \code{Psi}.
#' @param return_decomp Logical; if TRUE, return decomposition terms used for
#'   fast inverse transformation.
#'
#' @return A list with some of the elements:
#'   \item{W}{Whitening matrix (if \code{returnW = TRUE}).}
#'   \item{Phi}{Factor loadings in the original space (if \code{PhiPsi = TRUE}).}
#'   \item{Psi}{Standardized loadings (if \code{PhiPsi = TRUE}).}
#'
#' @export
Cholesky <- function(Sigma, returnW = TRUE, PhiPsi = TRUE, return_decomp = FALSE) {
  # Sigma must be symmetric and positive-definite
  if (!isTRUE(attr(Sigma, ".checked_spd"))) {
    .check_symmetric_pd(Sigma, require_pd = TRUE)
  }
  
  v <- diag(Sigma)
  chol_Sigma <- chol(Sigma)  # upper triangular
  
  result <- list()
  
  if (returnW) {
    # W: components x original features
    W <- backsolve(chol_Sigma, diag(ncol(Sigma)), transpose = TRUE)
    result$W <- .set_matrix_attr(
      W,
      row_names = paste0("L", seq_len(ncol(Sigma))),
      col_names = colnames(Sigma),
      method    = "Cholesky"
    )
  }
  
  if (PhiPsi) {
    Phi <- t(chol_Sigma)                 # loadings in original space
    Psi <- sweep(Phi, 1L, 1 / sqrt(v), "*") # standardized loadings
    
    row_names <- colnames(Sigma)
    col_names <- paste0("L", seq_len(ncol(Sigma)))
    
    result$Phi <- .set_matrix_attr(Phi, row_names, col_names, "Cholesky")
    result$Psi <- .set_matrix_attr(Psi, row_names, col_names, "Cholesky")
  }

  if (return_decomp) {
    result$decomp <- list(
      U = chol_Sigma,
      D = diag(chol_Sigma),
      inv_tW = chol_Sigma
    )
  }
  
  result
}
