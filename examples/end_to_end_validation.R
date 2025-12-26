# =============================================================================
# GWCF 端到端验证测试
# 
# 目的：在进入真实案例分析前，验证整个工作流
# 方法：使用模拟数据（已知真实效应），检验估计是否合理
# =============================================================================

cat("
================================================================================
                    GWCF 端到端验证测试
================================================================================
")

# -----------------------------------------------------------------------------
# 0. 环境准备
# -----------------------------------------------------------------------------
cat("\n【0】环境准备\n")

library(devtools)
load_all("gwcf", reset = TRUE)

library(sf)
library(terra)
library(ggplot2)
library(future)

plan(multisession, workers = 2)
set.seed(42)

cat("  ✓ 包加载完成\n")

# -----------------------------------------------------------------------------
# 1. 生成模拟数据（已知真实效应）
# -----------------------------------------------------------------------------
cat("\n【1】生成模拟数据\n")

n <- 500  # 样本量
coords <- cbind(
  X = runif(n, 0, 100),
  Y = runif(n, 0, 100)
)

# 协变量
X1 <- rnorm(n, 50, 10)
X2 <- runif(n, 0, 1)
X3 <- rnorm(n, 0, 1)
X <- cbind(X1 = X1, X2 = X2, X3 = X3)

# 二元处理：倾向性得分依赖于X1
propensity <- plogis(-2 + 0.04 * X1)
W <- rbinom(n, 1, propensity)

# 真实处理效应：空间异质性
# 效应在左下角大，右上角小
tau_true <- 2 + 3 * (1 - coords[,1]/100) * (1 - coords[,2]/100) + 0.5 * X2

# 潜在结果
Y0 <- 10 + 0.1 * X1 + 2 * X2 + rnorm(n, 0, 1)
Y1 <- Y0 + tau_true
Y <- ifelse(W == 1, Y1, Y0)

cat(sprintf("  样本量: %d\n", n))
cat(sprintf("  处理组: %d (%.1f%%)\n", sum(W), 100*mean(W)))
cat(sprintf("  真实效应范围: [%.2f, %.2f]\n", min(tau_true), max(tau_true)))

# -----------------------------------------------------------------------------
# 2. 测试 prep_from_sf
# -----------------------------------------------------------------------------
cat("\n【2】测试 prep_from_sf\n")

# 创建 sf 对象
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

cat("  ✓ prep_from_sf 成功\n")
print(data_sf)

# -----------------------------------------------------------------------------
# 3. 测试 prep_from_raster（小栅格）
# -----------------------------------------------------------------------------
cat("\n【3】测试 prep_from_raster\n")

# 创建小栅格
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

cat("  ✓ prep_from_raster 成功\n")
print(data_raster)

# -----------------------------------------------------------------------------
# 4. 模型拟合（使用sf数据）
# -----------------------------------------------------------------------------
cat("\n【4】模型拟合\n")

t_fit <- system.time({
  fit <- gw_causal_forest(
    Y = data_sf$Y,
    W = data_sf$W,
    X = data_sf$X,
    coords = data_sf$coords,
    treatment_type = "binary",
    bandwidth = NULL,  # 自动选择
    n_anchors = 15,
    kernel = "bisquare",
    n_folds = 3,
    num.trees = 500,
    seed = 42
  )
})

cat(sprintf("  ✓ 模型拟合成功，耗时 %.1f 秒\n", t_fit["elapsed"]))
print(fit)

# -----------------------------------------------------------------------------
# 5. 预测
# -----------------------------------------------------------------------------
cat("\n【5】预测\n")

preds <- predict(fit, newX = data_sf$X, newcoords = data_sf$coords)

cat(sprintf("  ✓ 预测成功\n"))
cat(sprintf("  估计效应范围: [%.3f, %.3f]\n", min(preds$tau_hat), max(preds$tau_hat)))
cat(sprintf("  平均SE: %.3f\n", mean(preds$se)))
cat(sprintf("  ESS范围: [%.1f, %.1f]\n", min(preds$ess), max(preds$ess)))

# -----------------------------------------------------------------------------
# 6. 验证：估计 vs 真实
# -----------------------------------------------------------------------------
cat("\n【6】验证：估计 vs 真实效应\n")

correlation <- cor(preds$tau_hat, tau_true)
rmse <- sqrt(mean((preds$tau_hat - tau_true)^2))
bias <- mean(preds$tau_hat - tau_true)

cat(sprintf("  相关系数: %.3f\n", correlation))
cat(sprintf("  RMSE: %.3f\n", rmse))
cat(sprintf("  偏差: %.3f\n", bias))

# 判断标准
if (correlation > 0.5) {
  cat("  ✓ 相关性检验通过 (r > 0.5)\n")
} else {
  cat("  ⚠ 相关性较低，可能需要调整参数\n")
}

# 可视化
p_validation <- ggplot(data.frame(true = tau_true, estimated = preds$tau_hat),
                       aes(x = true, y = estimated)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(title = sprintf("真实 vs 估计效应 (r = %.3f)", correlation),
       x = "真实效应", y = "估计效应")

ggsave("examples/plots/validation_true_vs_estimated.png", p_validation, 
       width = 6, height = 5, dpi = 150)
cat("  ✓ 验证图已保存\n")

# -----------------------------------------------------------------------------
# 7. 诊断
# -----------------------------------------------------------------------------
cat("\n【7】诊断\n")

diag <- diagnostics(fit)
print(diag)

# -----------------------------------------------------------------------------
# 8. 可视化
# -----------------------------------------------------------------------------
cat("\n【8】可视化\n")

# 效应地图
p_effect <- plot(fit, type = "effect_map", 
                 newX = data_sf$X, newcoords = data_sf$coords)
ggsave("examples/plots/validation_effect_map.png", p_effect, 
       width = 7, height = 6, dpi = 150)
cat("  ✓ 效应地图已保存\n")

# ESS地图
p_ess <- plot(fit, type = "ess_map")
ggsave("examples/plots/validation_ess_map.png", p_ess, 
       width = 7, height = 6, dpi = 150)
cat("  ✓ ESS地图已保存\n")

# CV曲线
if (!is.null(fit$cv_results)) {
  p_cv <- plot(fit, type = "cv_curve")
  ggsave("examples/plots/validation_cv_curve.png", p_cv, 
         width = 6, height = 4, dpi = 150)
  cat("  ✓ CV曲线已保存\n")
}

# Overlap图
p_overlap <- plot(fit, type = "overlap")
ggsave("examples/plots/validation_overlap.png", p_overlap, 
       width = 6, height = 4, dpi = 150)
cat("  ✓ Overlap图已保存\n")

# -----------------------------------------------------------------------------
# 9. Summary
# -----------------------------------------------------------------------------
cat("\n【9】Summary\n")

summ <- summary(fit)

# -----------------------------------------------------------------------------
# 10. 边界情况测试
# -----------------------------------------------------------------------------
cat("\n【10】边界情况测试\n")

# 测试外部点预测
far_coords <- matrix(c(200, 200, 300, 300), ncol = 2, byrow = TRUE)
far_X <- matrix(rnorm(6), ncol = 3)
colnames(far_X) <- c("x1", "x2", "x3")

pred_far <- tryCatch({
  suppressWarnings(predict(fit, far_X, far_coords))
}, error = function(e) {
  cat("  ✗ 外部点预测失败:", e$message, "\n")
  NULL
})

if (!is.null(pred_far)) {
  cat("  ✓ 外部点预测成功（应有警告）\n")
}

# -----------------------------------------------------------------------------
# 11. Bootstrap（可选，耗时）
# -----------------------------------------------------------------------------
cat("\n【11】Bootstrap推断（快速测试，5次迭代）\n")

t_boot <- system.time({
  boot_ci <- spatial_bootstrap(
    object = fit,
    newX = data_sf$X[1:10, ],
    newcoords = data_sf$coords[1:10, ],
    n_bootstrap = 5,  # 快速测试
    n_blocks = 5,
    conf_level = 0.95,
    seed = 123
  )
})

cat(sprintf("  ✓ Bootstrap完成，耗时 %.1f 秒\n", t_boot["elapsed"]))
cat(sprintf("  前3个点的CI宽度: %.3f, %.3f, %.3f\n", 
            boot_ci$ci_upper[1] - boot_ci$ci_lower[1],
            boot_ci$ci_upper[2] - boot_ci$ci_lower[2],
            boot_ci$ci_upper[3] - boot_ci$ci_lower[3]))

# -----------------------------------------------------------------------------
# 清理
# -----------------------------------------------------------------------------
plan(sequential)

# -----------------------------------------------------------------------------
# 总结
# -----------------------------------------------------------------------------
cat("\n")
cat("================================================================================\n")
cat("                         端到端验证完成\n")
cat("================================================================================\n")
cat("\n")
cat("结果摘要:\n")
cat(sprintf("  - 样本量: %d\n", n))
cat(sprintf("  - 处理组比例: %.1f%%\n", 100*mean(W)))
cat(sprintf("  - 估计与真实相关: r = %.3f\n", correlation))
cat(sprintf("  - RMSE: %.3f\n", rmse))
cat(sprintf("  - 选定带宽: %.2f\n", fit$bandwidth))
cat(sprintf("  - Anchor数量: %d\n", length(fit$anchor_forests)))
cat("\n")

if (correlation > 0.5 && rmse < 2) {
  cat("✓ 验证通过！可以进行真实案例分析。\n")
} else {
  cat("⚠ 验证结果需要关注。建议检查参数设置。\n")
}

cat("\n生成的图表:\n")
cat("  - examples/plots/validation_true_vs_estimated.png\n")
cat("  - examples/plots/validation_effect_map.png\n")
cat("  - examples/plots/validation_ess_map.png\n")
cat("  - examples/plots/validation_cv_curve.png\n")
cat("  - examples/plots/validation_overlap.png\n")

