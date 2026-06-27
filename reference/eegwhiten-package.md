# eegwhiten: Whitening Transforms for EEG Covariance Matrices

Whitening transforms (PCA, PCA on the correlation matrix, ZCA, ZCA on
the correlation matrix, SVD and Cholesky) with a simple fit/transform
API for EEG data, where rows are trials/epochs and columns are channels
or features. Provides reusable model objects, inverse transforms,
covariance shrinkage (Ledoit-Wolf and OAS) and robust covariance
estimators, cross-validated hyperparameter tuning, and diagnostics. Also
includes cross-session and cross-subject transfer-learning tools
(Euclidean alignment and per-domain recentering), Riemannian
tangent-space mapping of covariance matrices, and generalized-eigenvalue
("relative") whitening.

## See also

Useful links:

- <https://github.com/Yiming-S/eegwhiten>

- <https://yiming-s.github.io/eegwhiten/>

- Report bugs at <https://github.com/Yiming-S/eegwhiten/issues>

## Author

**Maintainer**: Yiming Shen <yiming.shen001@umb.edu>
