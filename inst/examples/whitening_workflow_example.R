
############################################################
# Example: Whitening simulated EEG epochs (trials x channels)
############################################################

# This script demonstrates how to:
#  1) Simulate correlated EEG-like data (trials x channels)
#  2) Fit a whitening model on training data
#  3) Apply whitening to train/test
#  4) Check whitening quality
#  5) Invert the transform (unwhiten) for interpretation

library(eegwhiten)

set.seed(123)

# ----------------------------------------------------------
# Simulate EEG-like data: 64 channels, AR(1)-like covariance
# ----------------------------------------------------------
n_train <- 200L   # number of training trials/epochs
n_test  <- 50L    # number of test trials
n_chan  <- 64L    # number of EEG channels

# Spatial covariance: channels with stronger correlation to neighbors
rho    <- 0.8
Sigma  <- toeplitz(rho ^ (0:(n_chan - 1L)))  # 64 x 64
chol_S <- chol(Sigma)

# Generate correlated Gaussian noise: X = Z %*% t(chol(Sigma))
Z_train <- matrix(rnorm(n_train * n_chan), n_train, n_chan)
Z_test  <- matrix(rnorm(n_test  * n_chan), n_test,  n_chan)

X_train <- Z_train %*% t(chol_S)  # training EEG features (trials x channels)
X_test  <- Z_test  %*% t(chol_S)  # test EEG features

colnames(X_train) <- paste0("Ch", seq_len(n_chan))
colnames(X_test)  <- colnames(X_train)

# ----------------------------------------------------------
# Fit whitening on training data and apply to train/test
# ----------------------------------------------------------
# Method can be "ZCA", "PCA", "PCA-cor", "ZCA-cor", "SVD", "Cholesky"
wmod <- whiten_model(X_train, center = TRUE, method = "ZCA")

Xw_train <- predict(wmod, X_train)  # whitened training data
Xw_test  <- predict(wmod, X_test)   # whitened test data

# ----------------------------------------------------------
# Check whitening quality on training data
# ----------------------------------------------------------
diag_res <- check_whitening(Xw_train)

str(diag_res)

# diag_res$diag_range  ~ c(1, 1) if perfect
# diag_res$offdiag_frob should be small if whitening is effective

# ----------------------------------------------------------
# Optional: invert the transform (unwhiten)
# ----------------------------------------------------------
X_train_recon <- unwhiten(Xw_train, wmod)

# Mean squared reconstruction error per channel
mse <- colMeans((X_train_recon - X_train) ^ 2)
summary(mse)

# If whitening model is numerically stable and invertible,
# MSE should be close to zero (up to floating-point precision).
