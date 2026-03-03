############################################################
# ZCA whitening on the correlation matrix
############################################################

#' ZCA whitening on the correlation matrix
#'
#' Perform ZCA whitening using the correlation matrix implied
#' by \code{Sigma}. This is a scale-invariant whitening transform
#' that equalizes variance and removes linear correlations.
#'
#' @param Sigma Symmetric positive-definite covariance matrix.
#' @param returnW Logical; if TRUE, return the whitening matrix \code{W}.
#' @param PhiPsi Logical; if TRUE, return factor loadings \code{Phi}
#'   and standardized loadings \code{Psi}.
#' @param return_decomp Logical; if TRUE, return decomposition terms used for
#'   fast inverse transformation.
#'
#' @return A list with some of the elements:
#'   \item{W}{Whitening matrix based on the correlation structure.}
#'   \item{Phi}{Factor loadings in the original space.}
#'   \item{Psi}{Standardized loadings.}
#'
#' @export
ZCA_cor <- function(Sigma, returnW = TRUE, PhiPsi = TRUE, return_decomp = FALSE) {
  if (!isTRUE(attr(Sigma, ".checked_spd"))) {
    .check_symmetric_pd(Sigma, require_pd = TRUE)
  }
  
  v <- diag(Sigma)
  if (any(v <= 0)) {
    stop("Diagonal elements of Sigma must be positive.")
  }
  
  R <- stats::cov2cor(Sigma)
  
  eR    <- eigen(R, symmetric = TRUE)
  theta <- eR$values
  G     <- eR$vectors
  
  if (any(theta <= 0)) {
    stop("Correlation matrix R must be positive-definite for ZCA-cor.")
  }
  
  result <- list()
  
  if (returnW) {
    # W: components x original features
    W <- sweep(G, 2L, 1 / sqrt(theta), "*") %*% t(G)
    W <- sweep(W, 2L, 1 / sqrt(v), "*")
    result$W <- .set_matrix_attr(
      W,
      row_names = paste0("L", seq_len(ncol(Sigma))),
      col_names = colnames(Sigma),
      method    = "ZCA-cor"
    )
  }
  
  if (PhiPsi) {
    Psi <- sweep(G, 2L, sqrt(theta), "*") %*% t(G)
    Phi <- sweep(Psi, 1L, sqrt(v), "*")
    
    row_names <- colnames(Sigma)
    col_names <- paste0("L", seq_len(ncol(Sigma)))
    
    result$Phi <- .set_matrix_attr(Phi, row_names, col_names, "ZCA-cor")
    result$Psi <- .set_matrix_attr(Psi, row_names, col_names, "ZCA-cor")
  }

  if (return_decomp) {
    inv_tW <- sweep(G, 2L, sqrt(theta), "*") %*% t(G)
    inv_tW <- sweep(inv_tW, 2L, sqrt(v), "*")
    result$decomp <- list(
      U = G,
      D = theta,
      scale_diag = v,
      inv_tW = inv_tW
    )
  }
  
  result
}
