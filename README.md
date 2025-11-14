# DA4BCI

**DA4BCI** is an R package that provides multiple domain adaptation methods tailored for EEG-based brain-computer interface (BCI) applications. It includes a unified interface for:

- **TCA (Transfer Component Analysis)**
- **SA (Subspace Alignment)**
- **MIDA (Maximum Independence Domain Adaptation)**
- **RD (Riemannian Distance)**
- **CORAL (Correlation Alignment)**
- **GFK (Geodesic Flow Kernel)**
- **ART (Aligned Riemannian Transport)**
- **PT (Parallel Transport)**
- **M3D (Manifold-based Multi-step Domain Adaptation)**
- **OT (Entropy-Regularized OT (Sinkhorn–Knopp) with Barycentric Mapping)**

These methods help align EEG data from different sessions or subjects, mitigating distributional shifts and enabling more robust learning.

## Distance Metrics and Evaluation

Additionally, DA4BCI implements various distance metrics and evaluation tools to quantitatively assess the effectiveness of domain adaptation:

- **Euclidean Distance Matrix:** Efficient computation of pairwise distances between datasets.
- **Wasserstein Distance:** Measures the minimal "cost" required to transform one distribution into another, emphasizing distribution alignment.
- **Maximum Mean Discrepancy (MMD):** Assesses differences between distributions using kernel methods, ideal for detecting subtle distributional shifts.
- **Energy Distance:** Captures differences between empirical distributions based on pairwise distances, useful for validating adaptation performance.

The `distanceSummary` function conveniently summarizes these metrics, providing a quick and comprehensive evaluation framework for domain adaptation results.

## Installation

1. Make sure you have R 3.5.0 or later.
2. In R, install the **remotes** (or **devtools**) package if you haven’t yet:
   ```r
   install.packages("remotes")
   remotes::install_github("Yiming-S/DA4BCI", force = TRUE)



   
