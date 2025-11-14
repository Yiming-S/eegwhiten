############################################################
# Convenience wrappers for single-matrix and batch whitening
############################################################

#' One-shot whitening for a single EEG matrix
#'
#' Convenience wrapper around \code{whiten_model()} that fits
#' a whitening model on \code{X} and immediately applies it to
#' \code{X}, returning both the whitened data and the underlying
#' model parameters.
#'
#' @param X Numeric matrix [n_trials x n_channels].
#' @param center Logical; whether to subtract the column means.
#' @param method Whitening method; one of \code{"SVD"}, \code{"ZCA"},
#'   \code{"ZCA-cor"}, \code{"PCA"}, \code{"PCA-cor"}, \code{"Cholesky"}.
#'
#' @return A list with components:
#'   \item{Z}{Whitened data matrix.}
#'   \item{W}{Whitening matrix.}
#'   \item{center}{Centering vector used.}
#'   \item{method}{Method name.}
#'
#' @examples
#' set.seed(42)
#' X <- matrix(rnorm(200 * 16), 200, 16)
#' colnames(X) <- paste0("Ch", seq_len(16))
#' res <- whiten_matrix(X, center = TRUE, method = "PCA")
#' str(res$Z)
#'
#' @export
whiten_matrix <- function(X, center = TRUE,
                          method = c("SVD", "ZCA", "ZCA-cor",
                                     "PCA", "PCA-cor", "Cholesky")) {
  method <- match.arg(method)
  
  model <- whiten_model(X, center = center, method = method)
  Z     <- predict(model, X)
  
  list(
    Z      = Z,             # whitened data
    W      = model$W,       # whitening matrix
    center = model$center,  # centering vector
    method = model$method   # method name
  )
}

#' Batch whitening for a list of EEG matrices
#'
#' Each matrix in \code{X_list} is whitened independently using
#' \code{whiten_matrix}. For proper train/test usage, fit a model
#' on training data with \code{whiten_model()} and call
#' \code{predict()} on test data instead of using this function.
#'
#' @param X_list List of numeric matrices, each [n_trials x n_channels].
#' @param center Logical; whether to subtract the column means in each matrix.
#' @param method Whitening method; one of \code{"SVD"}, \code{"ZCA"},
#'   \code{"ZCA-cor"}, \code{"PCA"}, \code{"PCA-cor"}, \code{"Cholesky"}.
#'
#' @return A list of results as returned by \code{whiten_matrix()}.
#'
#' @export
whiten_batch <- function(X_list, center = TRUE,
                         method = c("SVD", "ZCA", "ZCA-cor",
                                    "PCA", "PCA-cor", "Cholesky")) {
  method <- match.arg(method)
  if (!is.list(X_list)) stop("X_list must be a list of matrices.")
  
  lapply(
    X_list,
    function(X) whiten_matrix(X, center = center, method = method)
  )
}