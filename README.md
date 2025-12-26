# GWCF: Geographically Weighted Causal Forest

An R package for estimating **spatially varying heterogeneous treatment effects** using causal forests with kernel-weighted spatial averaging.

## Overview

GWCF combines [Generalized Random Forests](https://github.com/grf-labs/grf) (Athey et al., 2019) with geographically weighted methods to address two critical challenges in spatial causal inference:

1. **Spatial validation leakage**: Standard cross-validation ignores spatial autocorrelation, leading to overoptimistic model selection
2. **Variance underestimation**: Ignoring spatial structure leads to confidence intervals that are too narrow

The package implements spatially-aware cross-validation for bandwidth selection and spatial block bootstrap for uncertainty quantification.

## Installation

```r
# Install from GitHub
devtools::install_github("your-username/GWCF", subdir = "gwcf")
```

### Dependencies

- `grf` (>= 2.3.0)
- `sf`
- `future.apply`
- `ggplot2`
- `terra` (optional, for raster data)

## Quick Start

```r
library(gwcf)
library(future)

# Enable parallel processing
plan(multisession, workers = 4)

# Prepare data from shapefile
data <- prep_from_sf(
  sf_obj = my_shapefile,
  y_col = "outcome",
  w_col = "treatment",
  x_cols = c("covariate1", "covariate2", "covariate3")
)

# Fit model (bandwidth auto-selected via spatial CV)
fit <- gw_causal_forest(
  Y = data$Y,
  W = data$W,
  X = data$X,
  coords = data$coords,
  treatment_type = "binary"
)

# View results
print(fit)
summary(fit)

# Predict treatment effects
predictions <- predict(fit, newX = data$X, newcoords = data$coords)

# Visualize
plot(fit, type = "effect_map", newX = data$X, newcoords = data$coords)
plot(fit, type = "ess_map")
plot(fit, type = "cv_curve")

# Diagnostics
diagnostics(fit)

# Clean up
plan(sequential)
```

## Key Features

### Data Preprocessing

```r
# From shapefile/sf object
data <- prep_from_sf(sf_obj, y_col, w_col, x_cols)

# From raster (with memory-efficient sampling for large rasters)
data <- prep_from_raster(
  y_raster, w_raster, x_rasters,
  sample_n = 10000,
  sample_type = "stratified"  # Spatially balanced sampling
)
```

### Model Fitting

```r
fit <- gw_causal_forest(
  Y, W, X, coords,
  treatment_type = "binary",    # or "continuous"
  bandwidth = NULL,             # Auto-select via spatial CV
  kernel = "bisquare",          # or "gaussian", "exponential"
  n_anchors = NULL,             # Auto: sqrt(n)
  num.trees = 2000
)
```

### Inference

```r
# Point estimates with analytical variance
pred <- predict(fit, newX, newcoords, estimate_variance = TRUE)

# Spatial block bootstrap for robust CIs
boot_ci <- spatial_bootstrap(
  fit, newX, newcoords,
  n_bootstrap = 200,
  n_blocks = 10,
  conf_level = 0.95
)
```

### Visualization

| Plot Type | Description |
|-----------|-------------|
| `effect_map` | Spatial distribution of estimated treatment effects |
| `ess_map` | Effective sample size at anchor points |
| `cv_curve` | Cross-validation error vs bandwidth |
| `overlap` | Treatment variable distribution |

## Methodology

### Geographically Weighted Estimation

The model fits local causal forests at anchor points, weighted by spatial proximity:

```
τ(s) = Σ w(s, aⱼ) · τⱼ(X)
```

where `w(s, aⱼ)` is the kernel weight based on distance to anchor `aⱼ`.

### Bandwidth Selection

Optimal bandwidth is selected via spatial cross-validation using R-loss:

```
R-loss = E[(Y - m̂(X)) - τ(X)(W - ê(X))]²
```

Spatial blocks ensure validation data is spatially separated from training data.

### Supported Kernels

- **Bisquare** (default): Compact support, smooth
- **Gaussian**: Infinite support, decreasing influence
- **Exponential**: Infinite support, faster decay

## Example Output

```
Geographically Weighted Causal Forest - Summary
================================================

Model Configuration:
  Treatment type: binary
  Number of anchors: 15
  Bandwidth: 30.4462 (bisquare kernel)
  Training observations: 500

Global Average Treatment Effect (ATE):
  ATE: 3.1991 (SE: 0.0133)
  95% CI: [3.1730, 3.2253]
  Effect range: [1.7852, 4.7768]
```

## Performance Notes

- Large rasters (>1M cells) use memory-efficient sampling
- Parallel processing via `future` framework
- Nuisance models cached during CV (10x speedup)
- Bootstrap can be slow for large datasets; use analytical variance for exploration

## Citation

If you use this package, please cite:

```bibtex
@software{gwcf2025,
  author = {Han, Bo},
  title = {GWCF: Geographically Weighted Causal Forest},
  year = {2025},
  url = {https://github.com/your-username/GWCF}
}
```

## References

- Athey, S., Tibshirani, J., & Wager, S. (2019). Generalized random forests. *The Annals of Statistics*, 47(2), 1148-1178.
- Brunsdon, C., Fotheringham, A. S., & Charlton, M. E. (1996). Geographically weighted regression. *Geographical Analysis*, 28(4), 281-298.

## License

MIT License. See [LICENSE](gwcf/LICENSE) for details.

## Status

⚠️ **Alpha release** - Core functionality tested, API may change.

Contributions and bug reports welcome!

