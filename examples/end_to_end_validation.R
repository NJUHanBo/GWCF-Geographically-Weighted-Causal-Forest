# =============================================================================
# GWCF End-to-End Validation Test
# 
# Purpose: Validate the complete workflow before real case analysis
# Method: Use simulated data with known true effects to verify estimation
# =============================================================================

cat("
================================================================================
                    GWCF End-to-End Validation Test
================================================================================
")

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------
cat("\n[0] Environment Setup\n")

library(devtools)
load_all("gwcf", reset = TRUE)

library(sf)
library(terra)
library(ggplot2)
library(future)

plan(multisession, workers = 2)
set.seed(42)

cat("  OK: Packages loaded\n")

# -----------------------------------------------------------------------------
# 1. Generate Simulated Data (with known true effects)
# -----------------------------------------------------------------------------
cat("\n[1] Generate Simulated Data\n")

n <- 500  # sample size
coords <- cbind(
  X = runif(n, 0, 100),
  Y = runif(n, 0, 100)
)

# Covariates
X1 <- rnorm(n, 50, 10)
X2 <- runif(n, 0, 1)
X3 <- rnorm(n, 0, 1)
X <- cbind(X1 = X1, X2 = X2, X3 = X3)

# Binary treatment: propensity depends on X1
propensity <- plogis(-2 + 0.04 * X1)
W <- rbinom(n, 1, propensity)

# True treatment effect: spatial heterogeneity
# Effect is larger in bottom-left, smaller in top-right
tau_true <- 2 + 3 * (1 - coords[,1]/100) * (1 - coords[,2]/100) + 0.5 * X2

# Potential outcomes
Y0 <- 10 + 0.1 * X1 + 2 * X2 + rnorm(n, 0, 1)
Y1 <- Y0 + tau_true
Y <- ifelse(W == 1, Y1, Y0)

cat(sprintf("  Sample size: %d\n", n))
cat(sprintf("  Treated: %d (%.1f%%)\n", sum(W), 100*mean(W)))
cat(sprintf("  True effect range: [%.2f, %.2f]\n", min(tau_true), max(tau_true)))

# -----------------------------------------------------------------------------
# 2. Test prep_from_sf
# -----------------------------------------------------------------------------
cat("\n[2] Test prep_from_sf\n")

# Create sf object
sim_sf <- st_as_sf(
  data.frame(
    outcome = Y,
    treatment = W,
    x1 = X1, x2 = X2, x3 = X3,
    lon = coords[,1],
    lat = coords[,2]
  ),
  coords = c("lon", "lat"),
  crs = 4326
)

data_sf <- prep_from_sf(
  sf_obj = sim_sf,
  y_col = "outcome",
  w_col = "treatment",
  x_cols = c("x1", "x2", "x3")
)

cat("  OK: prep_from_sf succeeded\n")
print(data_sf)

# -----------------------------------------------------------------------------
# 3. Test prep_from_raster (small raster)
# -----------------------------------------------------------------------------
cat("\n[3] Test prep_from_raster\n")

# Create small raster
r <- rast(nrows = 50, ncols = 50, xmin = 0, xmax = 100, ymin = 0, ymax = 100)
y_r <- r; values(y_r) <- rnorm(ncell(r), 10, 2)
w_r <- r; values(w_r) <- runif(ncell(r))
x_r <- c(r, r)
values(x_r) <- cbind(rnorm(ncell(r)), runif(ncell(r)))
names(x_r) <- c("var1", "var2")

data_raster <- prep_from_raster(
  y_raster = y_r,
  w_raster = w_r,
  x_rasters = x_r,
  sample_n = 300,
  sample_type = "stratified",
  seed = 123
)

cat("  OK: prep_from_raster succeeded\n")
print(data_raster)

# -----------------------------------------------------------------------------
# 4. Model Fitting (using sf data)
# -----------------------------------------------------------------------------
cat("\n[4] Model Fitting\n")

t_fit <- system.time({
  fit <- gw_causal_forest(
    Y = data_sf$Y,
    W = data_sf$W,
    X = data_sf$X,
    coords = data_sf$coords,
    treatment_type = "binary",
    bandwidth = NULL,  # auto-select
    n_anchors = 15,
    kernel = "bisquare",
    n_folds = 3,
    num.trees = 500,
    seed = 42
  )
})

cat(sprintf("  OK: Model fitted in %.1f seconds\n", t_fit["elapsed"]))
print(fit)

# -----------------------------------------------------------------------------
# 5. Prediction
# -----------------------------------------------------------------------------
cat("\n[5] Prediction\n")

preds <- predict(fit, newX = data_sf$X, newcoords = data_sf$coords)

cat(sprintf("  OK: Prediction succeeded\n"))
cat(sprintf("  Estimated effect range: [%.3f, %.3f]\n", min(preds$tau_hat), max(preds$tau_hat)))
cat(sprintf("  Mean SE: %.3f\n", mean(preds$se)))
cat(sprintf("  ESS range: [%.1f, %.1f]\n", min(preds$ess), max(preds$ess)))

# -----------------------------------------------------------------------------
# 6. Validation: Estimated vs True
# -----------------------------------------------------------------------------
cat("\n[6] Validation: Estimated vs True Effects\n")

correlation <- cor(preds$tau_hat, tau_true)
rmse <- sqrt(mean((preds$tau_hat - tau_true)^2))
bias <- mean(preds$tau_hat - tau_true)

cat(sprintf("  Correlation: %.3f\n", correlation))
cat(sprintf("  RMSE: %.3f\n", rmse))
cat(sprintf("  Bias: %.3f\n", bias))

# Threshold check
if (correlation > 0.5) {
  cat("  OK: Correlation test passed (r > 0.5)\n")
} else {
  cat("  WARNING: Low correlation, consider adjusting parameters\n")
}

# Visualization
p_validation <- ggplot(data.frame(true = tau_true, estimated = preds$tau_hat),
                       aes(x = true, y = estimated)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(title = sprintf("True vs Estimated Effects (r = %.3f)", correlation),
       x = "True Effect", y = "Estimated Effect")

ggsave("examples/plots/validation_true_vs_estimated.png", p_validation, 
       width = 6, height = 5, dpi = 150)
cat("  OK: Validation plot saved\n")

# -----------------------------------------------------------------------------
# 7. Diagnostics
# -----------------------------------------------------------------------------
cat("\n[7] Diagnostics\n")

diag <- diagnostics(fit)
print(diag)

# -----------------------------------------------------------------------------
# 8. Visualization
# -----------------------------------------------------------------------------
cat("\n[8] Visualization\n")

# Effect map
p_effect <- plot(fit, type = "effect_map", 
                 newX = data_sf$X, newcoords = data_sf$coords)
ggsave("examples/plots/validation_effect_map.png", p_effect, 
       width = 7, height = 6, dpi = 150)
cat("  OK: Effect map saved\n")

# ESS map
p_ess <- plot(fit, type = "ess_map")
ggsave("examples/plots/validation_ess_map.png", p_ess, 
       width = 7, height = 6, dpi = 150)
cat("  OK: ESS map saved\n")

# CV curve
if (!is.null(fit$cv_results)) {
  p_cv <- plot(fit, type = "cv_curve")
  ggsave("examples/plots/validation_cv_curve.png", p_cv, 
         width = 6, height = 4, dpi = 150)
  cat("  OK: CV curve saved\n")
}

# Overlap plot
p_overlap <- plot(fit, type = "overlap")
ggsave("examples/plots/validation_overlap.png", p_overlap, 
       width = 6, height = 4, dpi = 150)
cat("  OK: Overlap plot saved\n")

# -----------------------------------------------------------------------------
# 9. Summary
# -----------------------------------------------------------------------------
cat("\n[9] Summary\n")

summ <- summary(fit)

# -----------------------------------------------------------------------------
# 10. Edge Case Test
# -----------------------------------------------------------------------------
cat("\n[10] Edge Case Test\n")

# Test prediction for points outside study area
far_coords <- matrix(c(200, 200, 300, 300), ncol = 2, byrow = TRUE)
far_X <- matrix(rnorm(6), ncol = 3)
colnames(far_X) <- c("x1", "x2", "x3")

pred_far <- tryCatch({
  suppressWarnings(predict(fit, far_X, far_coords))
}, error = function(e) {
  cat("  FAIL: Prediction for external points failed:", e$message, "\n")
  NULL
})

if (!is.null(pred_far)) {
  cat("  OK: External point prediction succeeded (warning expected)\n")
}

# -----------------------------------------------------------------------------
# 11. Bootstrap (optional, time-consuming)
# -----------------------------------------------------------------------------
cat("\n[11] Bootstrap Inference (quick test, 5 iterations)\n")

t_boot <- system.time({
  boot_ci <- spatial_bootstrap(
    object = fit,
    newX = data_sf$X[1:10, ],
    newcoords = data_sf$coords[1:10, ],
    n_bootstrap = 5,  # quick test
    n_blocks = 5,
    conf_level = 0.95,
    seed = 123
  )
})

cat(sprintf("  OK: Bootstrap completed in %.1f seconds\n", t_boot["elapsed"]))
cat(sprintf("  CI width for first 3 points: %.3f, %.3f, %.3f\n", 
            boot_ci$ci_upper[1] - boot_ci$ci_lower[1],
            boot_ci$ci_upper[2] - boot_ci$ci_lower[2],
            boot_ci$ci_upper[3] - boot_ci$ci_lower[3]))

# -----------------------------------------------------------------------------
# Cleanup
# -----------------------------------------------------------------------------
plan(sequential)

# -----------------------------------------------------------------------------
# Final Summary
# -----------------------------------------------------------------------------
cat("\n")
cat("================================================================================\n")
cat("                         End-to-End Validation Complete\n")
cat("================================================================================\n")
cat("\n")
cat("Results Summary:\n")
cat(sprintf("  - Sample size: %d\n", n))
cat(sprintf("  - Treatment proportion: %.1f%%\n", 100*mean(W)))
cat(sprintf("  - Estimated vs true correlation: r = %.3f\n", correlation))
cat(sprintf("  - RMSE: %.3f\n", rmse))
cat(sprintf("  - Selected bandwidth: %.2f\n", fit$bandwidth))
cat(sprintf("  - Number of anchors: %d\n", length(fit$anchor_forests)))
cat("\n")

if (correlation > 0.5 && rmse < 2) {
  cat("PASSED: Validation successful. Ready for real case analysis.\n")
} else {
  cat("WARNING: Results need attention. Check parameter settings.\n")
}

cat("\nGenerated plots:\n")
cat("  - examples/plots/validation_true_vs_estimated.png\n")
cat("  - examples/plots/validation_effect_map.png\n")
cat("  - examples/plots/validation_ess_map.png\n")
cat("  - examples/plots/validation_cv_curve.png\n")
cat("  - examples/plots/validation_overlap.png\n")
