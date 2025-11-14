############################################################
# ZCA whitening
############################################################

#' ZCA whitening on the covariance matrix
#'
#' Perform ZCA (Zero-phase Component Analysis) whitening on a
#' symmetric positive-definite covariance matrix. This transform
#' keeps the whitened data as close as possible (in L2 sense) to
#' the original data while decorrelating components.
#'
#' @param Sigma Symmetric positive-definite covariance matrix.
#' @param returnW Logical; if TRUE, return the whitening matrix \code{W}.
#' @param PhiPsi Logical; if TRUE, return factor loadings \code{Phi}
#'   and standardized loadings \code{Psi}.
#'
#' @return A list with some of the elements:
#'   \item{W}{ZCA whitening matrix.}
#'   \item{Phi}{Loadings in the original space.}
#'   \item{Psi}{Standardized loadings.}
#'
#' @export
ZCA <- function(Sigma, returnW = TRUE, PhiPsi = TRUE) {
  .check_symmetric_pd(Sigma, require_pd = TRUE)
  
  v <- diag(Sigma)
  
  eSigma <- eigen(Sigma, symmetric = TRUE)
  lambda <- eSigma$values
  U      <- eSigma$vectors
  
  if (any(lambda <= 0)) {
    stop("Sigma must be positive-definite (eigenvalues > 0) for ZCA.")
  }
  
  result <- list()
  
  # Precompute U * diag(1/sqrt(lambda)) * U^T
  U_diag <- U %*% diag(1 / sqrt(lambda)) %*% t(U)
  
  if (returnW) {
    W <- U_diag
    result$W <- .set_matrix_attr(
      W,
      row_names = paste0("L", seq_len(ncol(Sigma))),
      col_names = colnames(Sigma),
      method    = "ZCA"
    )
  }
  
  if (PhiPsi) {
    Phi <- U_diag                          # loadings in original space
    Psi <- diag(1 / sqrt(v)) %*% Phi       # standardized loadings
    
    row_names <- colnames(Sigma)
    col_names <- paste0("L", seq_len(ncol(Sigma)))
    
    result$Phi <- .set_matrix_attr(Phi, row_names, col_names, "ZCA")
    result$Psi <- .set_matrix_attr(Psi, row_names, col_names, "ZCA")
  }
  
  result
}