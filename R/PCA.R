############################################################
# PCA-based whitening on the covariance matrix
############################################################

#' PCA-based whitening on the covariance matrix
#'
#' Perform whitening based on the eigen-decomposition of a
#' symmetric positive-definite covariance matrix. This is a
#' standard PCA whitening transform where components correspond
#' to eigenvectors of \code{Sigma}.
#'
#' @param Sigma Symmetric positive-definite covariance matrix.
#' @param returnW Logical; if TRUE, return the whitening matrix \code{W}.
#' @param PhiPsi Logical; if TRUE, return factor loadings \code{Phi}
#'   and standardized loadings \code{Psi}.
#'
#' @return A list with some of the elements:
#'   \item{W}{Whitening matrix (if \code{returnW = TRUE}).}
#'   \item{Phi}{Factor loadings in the original EEG feature space.}
#'   \item{Psi}{Standardized loadings.}
#'
#' @export
PCA <- function(Sigma, returnW = TRUE, PhiPsi = TRUE) {
  .check_symmetric_pd(Sigma, require_pd = TRUE)
  
  v <- diag(Sigma)
  if (any(v <= 0)) {
    stop("Diagonal elements of Sigma must be positive.")
  }
  
  eSigma <- eigen(Sigma, symmetric = TRUE)
  U      <- eSigma$vectors
  lambda <- eSigma$values
  
  # Fix sign ambiguity in eigenvectors by making diagonal of U positive
  U <- sweep(U, 2, sign(diag(U)), "*")
  
  result <- list()
  
  if (returnW) {
    # W: components x original features
    W <- tcrossprod(diag(1 / sqrt(lambda)), U)
    result$W <- .set_matrix_attr(
      W,
      row_names = paste0("L", seq_len(ncol(Sigma))),
      col_names = colnames(Sigma),
      method    = "PCA"
    )
  }
  
  if (PhiPsi) {
    Phi <- U %*% diag(sqrt(lambda))         # loadings
    Psi <- diag(1 / sqrt(v)) %*% Phi        # standardized loadings
    
    row_names <- colnames(Sigma)
    col_names <- paste0("L", seq_len(ncol(Sigma)))
    
    result$Phi <- .set_matrix_attr(Phi, row_names, col_names, "PCA")
    result$Psi <- .set_matrix_attr(Psi, row_names, col_names, "PCA")
  }
  
  result
}