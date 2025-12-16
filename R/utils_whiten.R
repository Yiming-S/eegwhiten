
############################################################
# FILENAME: utils_whiten.R
# DESCRIPTION: Internal utilities
# UPDATES: Added check_condition
############################################################

# Set dimnames and method attribute on a matrix
.set_matrix_attr <- function(x, row_names, col_names, method) {
  if (!is.null(row_names)) rownames(x) <- row_names
  if (!is.null(col_names)) colnames(x) <- col_names
  attr(x, "method") <- method
  x
}

# Basic check: Sigma must be symmetric; optionally positive-definite
.check_symmetric_pd <- function(Sigma, tol = 1e-8, require_pd = TRUE) {
  if (!is.matrix(Sigma)) stop("Sigma must be a matrix.")
  if (!isTRUE(all.equal(Sigma, t(Sigma)))) {
    stop("Sigma must be symmetric.")
  }
  eig <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  if (require_pd && any(eig <= tol)) {
    # More informative error for debugging
    min_eig <- min(eig)
    stop(sprintf("Sigma must be positive-definite. Min eigenvalue = %g. Try using lambda > 0 in whiten_model().", min_eig))
  }
  invisible(eig)
}

#' Check Condition Number of Data Matrix
#' 
#' @param X Numeric matrix
#' @return The ratio of max/min eigenvalue of cov(X). High values indicate ill-conditioning.
#' @export
check_condition <- function(X) {
  if(!is.matrix(X)) stop("X must be a matrix")
  ev <- eigen(cov(X), symmetric = TRUE, only.values = TRUE)$values
  # Filter small/negative values for ratio calculation
  ev <- ev[ev > 0]
  if(length(ev) == 0) return(Inf)
  max(ev) / min(ev)
}