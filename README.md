# eegwhiten

Whitening transforms for EEG covariance matrices (trials × channels).

`eegwhiten` provides a small, focused set of whitening methods (PCA, ZCA, SVD, Cholesky, and their correlation-based variants) with a **model-based API** designed for EEG and other multichannel time series:

- Rows = trials / epochs  
- Columns = channels or features  
- Fit whitening on training data and reuse on validation / test data  
- Optionally invert the transform and check whitening quality

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("YimingShen/eegwhiten")
