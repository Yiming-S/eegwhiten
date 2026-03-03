############################################################
# Tensor wrappers for EEG whitening
############################################################

# Convert [n_trials, n_channels, n_time] -> [n_trials * n_time, n_channels]
.tensor_to_matrix <- function(X_tensor) {
  if (!is.array(X_tensor) || length(dim(X_tensor)) != 3L) {
    stop("X_tensor must be a 3D numeric array with shape [n_trials, n_channels, n_time].")
  }
  if (!is.numeric(X_tensor) || any(!is.finite(X_tensor))) {
    stop("X_tensor must contain only finite numeric values.")
  }

  dims <- dim(X_tensor)
  n_trials <- dims[1]
  n_channels <- dims[2]
  n_time <- dims[3]

  X_perm <- aperm(X_tensor, c(1, 3, 2))
  X_mat <- array(X_perm, dim = c(n_trials * n_time, n_channels))

  list(
    X_mat = X_mat,
    n_trials = n_trials,
    n_channels = n_channels,
    n_time = n_time
  )
}

#' Fit a whitening model from a 3D EEG tensor
#'
#' Flattens \code{[n_trials, n_channels, n_time]} into
#' \code{[n_trials * n_time, n_channels]} and fits \code{whiten_model()}.
#'
#' @param X_tensor Numeric 3D array \code{[n_trials, n_channels, n_time]}.
#' @param ... Passed to \code{whiten_model()}.
#'
#' @return A \code{whiten_model} object.
#' @export
whiten_model_tensor <- function(X_tensor, ...) {
  tm <- .tensor_to_matrix(X_tensor)
  whiten_model(tm$X_mat, ...)
}

#' Apply a whitening model to a 3D EEG tensor
#'
#' @param model A \code{whiten_model} object.
#' @param X_tensor Numeric 3D array \code{[n_trials, n_channels, n_time]}.
#'
#' @return Whitened 3D array \code{[n_trials, n_components, n_time]}.
#' @export
predict_tensor <- function(model, X_tensor) {
  if (!inherits(model, "whiten_model")) {
    stop("model must be a 'whiten_model' object.")
  }

  tm <- .tensor_to_matrix(X_tensor)
  if (tm$n_channels != model$dim_in) {
    stop(sprintf("Channel mismatch: model expects %d channels, tensor has %d.", model$dim_in, tm$n_channels))
  }

  Z_mat <- stats::predict(model, tm$X_mat)
  Z_perm <- array(Z_mat, dim = c(tm$n_trials, tm$n_time, model$dim_out))
  aperm(Z_perm, c(1, 3, 2))
}

#' Inverse transform a whitened 3D EEG tensor
#'
#' @param Z_tensor Numeric 3D array \code{[n_trials, n_components, n_time]}.
#' @param model A \code{whiten_model} object.
#'
#' @return Reconstructed 3D array \code{[n_trials, n_channels, n_time]}.
#' @export
unwhiten_tensor <- function(Z_tensor, model) {
  if (!inherits(model, "whiten_model")) {
    stop("model must be a 'whiten_model' object.")
  }
  if (!is.array(Z_tensor) || length(dim(Z_tensor)) != 3L) {
    stop("Z_tensor must be a 3D numeric array with shape [n_trials, n_components, n_time].")
  }
  if (!is.numeric(Z_tensor) || any(!is.finite(Z_tensor))) {
    stop("Z_tensor must contain only finite numeric values.")
  }

  dims <- dim(Z_tensor)
  n_trials <- dims[1]
  n_comp <- dims[2]
  n_time <- dims[3]
  if (n_comp != model$dim_out) {
    stop(sprintf("Component mismatch: model expects %d components, tensor has %d.", model$dim_out, n_comp))
  }

  Z_perm <- aperm(Z_tensor, c(1, 3, 2))
  Z_mat <- array(Z_perm, dim = c(n_trials * n_time, n_comp))
  X_mat <- unwhiten_fast(Z_mat, model)

  X_perm <- array(X_mat, dim = c(n_trials, n_time, model$dim_in))
  aperm(X_perm, c(1, 3, 2))
}
