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
#'
#' @return A list with some of the elements:
#'   \item{W}{Whitening matrix (if \code{returnW = TRUE}).}
#'   \item{Phi}{Factor loadings in the original space (if \code{PhiPsi = TRUE}).}
#'   \item{Psi}{Standardized loadings (if \code{PhiPsi = TRUE}).}
#'
#' @export
Cholesky <- function(Sigma, returnW = TRUE, PhiPsi = TRUE) {
  # Sigma must be symmetric and positive-definite
  .check_symmetric_pd(Sigma, require_pd = TRUE)
  
  v <- diag(Sigma)
  chol_Sigma <- chol(Sigma)  # upper triangular
  
  result <- list()
  
  if (returnW) {
    # W: components x original features
    W <- solve(t(chol_Sigma))
    result$W <- .set_matrix_attr(
      W,
      row_names = paste0("L", seq_len(ncol(Sigma))),
      col_names = colnames(Sigma),
      method    = "Cholesky"
    )
  }
  
  if (PhiPsi) {
    Phi <- t(chol_Sigma)                 # loadings in original space
    Psi <- diag(1 / sqrt(v)) %*% Phi     # standardized loadings
    
    row_names <- colnames(Sigma)
    col_names <- paste0("L", seq_len(ncol(Sigma)))
    
    result$Phi <- .set_matrix_attr(Phi, row_names, col_names, "Cholesky")
    result$Psi <- .set_matrix_attr(Psi, row_names, col_names, "Cholesky")
  }
  
  result
}