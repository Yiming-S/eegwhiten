# Package index

## Whitening models & transforms

Fit a whitening model, apply it, and invert it.

- [`whiten_model()`](https://yiming-s.github.io/eegwhiten/reference/whiten_model.md)
  : Fit a whitening model for EEG data
- [`whiten_fit()`](https://yiming-s.github.io/eegwhiten/reference/whiten_fit.md)
  : Fit a whitening matrix from a covariance matrix
- [`whiten_matrix()`](https://yiming-s.github.io/eegwhiten/reference/whiten_matrix.md)
  : One-shot whitening for a single EEG matrix
- [`predict(`*`<whiten_model>`*`)`](https://yiming-s.github.io/eegwhiten/reference/predict.whiten_model.md)
  : Apply a whitening model to new EEG data
- [`unwhiten()`](https://yiming-s.github.io/eegwhiten/reference/unwhiten.md)
  : Inverse transform from whitened space back to original EEG space
- [`unwhiten_fast()`](https://yiming-s.github.io/eegwhiten/reference/unwhiten_fast.md)
  : Optimized inverse transform from whitened space
- [`whitening_patterns()`](https://yiming-s.github.io/eegwhiten/reference/whitening_patterns.md)
  : Spatial filters and forward patterns of a whitening model

## Whitening methods

Low-level whitening transforms on a covariance matrix.

- [`PCA()`](https://yiming-s.github.io/eegwhiten/reference/PCA.md) :
  PCA-based whitening on the covariance matrix
- [`PCA_cor()`](https://yiming-s.github.io/eegwhiten/reference/PCA_cor.md)
  : PCA-based whitening on the correlation matrix
- [`ZCA()`](https://yiming-s.github.io/eegwhiten/reference/ZCA.md) : ZCA
  whitening on the covariance matrix
- [`ZCA_cor()`](https://yiming-s.github.io/eegwhiten/reference/ZCA_cor.md)
  : ZCA whitening on the correlation matrix
- [`SVD()`](https://yiming-s.github.io/eegwhiten/reference/SVD.md) :
  SVD-based whitening on the covariance matrix
- [`Cholesky()`](https://yiming-s.github.io/eegwhiten/reference/Cholesky.md)
  : Cholesky-based whitening on a covariance matrix

## Tuning & reporting

- [`auto_tune_whitening()`](https://yiming-s.github.io/eegwhiten/reference/auto_tune_whitening.md)
  : Automatically Tune Whitening Hyperparameters
- [`report_whitening()`](https://yiming-s.github.io/eegwhiten/reference/report_whitening.md)
  : One-Click Whitening Report

## Diagnostics

- [`check_whitening()`](https://yiming-s.github.io/eegwhiten/reference/check_whitening.md)
  : Diagnostic check for whitening quality
- [`check_condition()`](https://yiming-s.github.io/eegwhiten/reference/check_condition.md)
  : Check Condition Number of Data Matrix

## Online & batch

- [`whiten_model_update()`](https://yiming-s.github.io/eegwhiten/reference/whiten_model_update.md)
  : Incrementally Update a Whitening Model
- [`whiten_batch()`](https://yiming-s.github.io/eegwhiten/reference/whiten_batch.md)
  : Batch whitening for a list of EEG matrices

## Tensor data

Whitening for trial x channel x time arrays.

- [`whiten_model_tensor()`](https://yiming-s.github.io/eegwhiten/reference/whiten_model_tensor.md)
  : Fit a whitening model from a 3D EEG tensor
- [`predict_tensor()`](https://yiming-s.github.io/eegwhiten/reference/predict_tensor.md)
  : Apply a whitening model to a 3D EEG tensor
- [`unwhiten_tensor()`](https://yiming-s.github.io/eegwhiten/reference/unwhiten_tensor.md)
  : Inverse transform a whitened 3D EEG tensor

## Alignment & transfer learning

Cross-session / cross-subject covariance alignment.

- [`euclidean_alignment()`](https://yiming-s.github.io/eegwhiten/reference/euclidean_alignment.md)
  : Euclidean Alignment for EEG Matrices
- [`recenter()`](https://yiming-s.github.io/eegwhiten/reference/recenter.md)
  : Recenter a single EEG matrix to the identity
- [`predict(`*`<ea_model>`*`)`](https://yiming-s.github.io/eegwhiten/reference/predict.ea_model.md)
  : Apply an EA Model to New Data
- [`inverse_ea()`](https://yiming-s.github.io/eegwhiten/reference/inverse_ea.md)
  : Inverse Transform for EA-Aligned Data
- [`whiten_relative()`](https://yiming-s.github.io/eegwhiten/reference/whiten_relative.md)
  : Relative whitening via a generalized eigenvalue decomposition

## Riemannian geometry

Tangent-space mapping for SPD covariance matrices.

- [`epoch_covariances()`](https://yiming-s.github.io/eegwhiten/reference/epoch_covariances.md)
  : Per-epoch channel covariance matrices from an EEG tensor
- [`tangent_space()`](https://yiming-s.github.io/eegwhiten/reference/tangent_space.md)
  : Map SPD covariance matrices to the tangent space
- [`untangent_space()`](https://yiming-s.github.io/eegwhiten/reference/untangent_space.md)
  : Map tangent-space vectors back to SPD covariance matrices

## S3 methods

- [`print(`*`<whiten_model>`*`)`](https://yiming-s.github.io/eegwhiten/reference/print.whiten_model.md)
  : Print a whitening model
- [`summary(`*`<whiten_model>`*`)`](https://yiming-s.github.io/eegwhiten/reference/summary.whiten_model.md)
  : Summarize a whitening model
- [`plot(`*`<whiten_model>`*`)`](https://yiming-s.github.io/eegwhiten/reference/plot.whiten_model.md)
  : Plot covariance deviation from identity for whitened data
- [`print(`*`<summary.whiten_model>`*`)`](https://yiming-s.github.io/eegwhiten/reference/print.summary.whiten_model.md)
  : Print summary for a whitening model
- [`print(`*`<whiten_tune>`*`)`](https://yiming-s.github.io/eegwhiten/reference/print.whiten_tune.md)
  : Print a whitening tuning result
- [`predict(`*`<whiten_tune>`*`)`](https://yiming-s.github.io/eegwhiten/reference/predict.whiten_tune.md)
  : Apply the best tuned whitening model to new data
- [`print(`*`<recenter>`*`)`](https://yiming-s.github.io/eegwhiten/reference/print.recenter.md)
  : Print a recentering model
- [`predict(`*`<recenter>`*`)`](https://yiming-s.github.io/eegwhiten/reference/predict.recenter.md)
  : Apply a recentering model to new data
- [`print(`*`<relative_whiten>`*`)`](https://yiming-s.github.io/eegwhiten/reference/print.relative_whiten.md)
  : Print a relative-whitening model
- [`predict(`*`<relative_whiten>`*`)`](https://yiming-s.github.io/eegwhiten/reference/predict.relative_whiten.md)
  : Apply a relative-whitening model to new data
- [`print(`*`<ea_model>`*`)`](https://yiming-s.github.io/eegwhiten/reference/print.ea_model.md)
  : Print an EA model
- [`summary(`*`<ea_model>`*`)`](https://yiming-s.github.io/eegwhiten/reference/summary.ea_model.md)
  : Summarize an EA model

## Package

- [`eegwhiten`](https://yiming-s.github.io/eegwhiten/reference/eegwhiten-package.md)
  [`eegwhiten-package`](https://yiming-s.github.io/eegwhiten/reference/eegwhiten-package.md)
  : eegwhiten: Whitening Transforms for EEG Covariance Matrices
