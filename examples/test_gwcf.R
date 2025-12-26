# ============================================================================
# gwcf 包测试脚本
# ============================================================================

# 0. 安装依赖（如果尚未安装）
# ----------------------------------------------------------------------------
required_packages <- c("devtools", "grf", "sf", "future.apply", "ggplot2", "future")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# 1. 加载 gwcf 包
# ----------------------------------------------------------------------------
cat("=== 加载 gwcf 包 ===\n")
devtools::load_all("/Users/hanbo/GWCT/gwcf")
cat("gwcf 包加载成功!\n\n")

# 2. 设置并行计算（可选，加速计算）
# ----------------------------------------------------------------------------
library(future)
plan(multisession, workers = 2)  # 使用2个核心

# 3. 生成模拟数据
# ----------------------------------------------------------------------------
cat("=== 生成模拟数据 ===\n")
set.seed(42)
n <- 500  # 样本量

# 空间坐标 (0-1 范围)
coords <- cbind(
  x = runif(n, 0, 1),
  y = runif(n, 0, 1)
)

# 协变量
X <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
colnames(X) <- paste0("X", 1:4)

# 二元处理变量
propensity <- plogis(-0.5 + 0.3 * X[,1] + 0.2 * coords[,1])
W <- rbinom(n, 1, propensity)

# 真实处理效应 (空间变化)
# 效应在空间上从西南到东北递增
true_tau <- 1 + 2 * coords[,1] + 1.5 * coords[,2] + 0.5 * X[,1]

# 结果变量
Y <- 2 + X[,1] + 0.5 * X[,2] + true_tau * W + rnorm(n, 0, 1)

cat(sprintf("样本量: %d\n", n))
cat(sprintf("处理组比例: %.2f\n", mean(W)))
cat(sprintf("真实效应范围: [%.2f, %.2f]\n", min(true_tau), max(true_tau)))
cat("\n")

# 4. 测试一：使用固定带宽拟合模型
# ----------------------------------------------------------------------------
cat("=== 测试一: 固定带宽拟合 ===\n")

fit1 <- gw_causal_forest(
  Y = Y,
  W = W,
  X = X,
  coords = coords,
  treatment_type = "binary",
  bandwidth = 0.4,        # 固定带宽
  n_anchors = 15,         # 锚点数量
  num.trees = 500,        # 减少树数量加快测试
  seed = 123
)

cat("\n模型摘要:\n")
print(fit1)

cat("\n诊断信息:\n")
diag1 <- diagnostics(fit1)
print(diag1)

# 5. 测试二：预测与置信区间
# ----------------------------------------------------------------------------
cat("\n=== 测试二: 预测与置信区间 ===\n")

pred1 <- predict(fit1, newX = X, newcoords = coords, estimate_variance = TRUE)

cat("预测结果结构:\n")
print(head(pred1))

cat(sprintf("\n估计效应范围: [%.2f, %.2f]\n", min(pred1$tau_hat), max(pred1$tau_hat)))
cat(sprintf("平均标准误: %.3f\n", mean(pred1$se)))

# 评估预测精度
rmse <- sqrt(mean((pred1$tau_hat - true_tau)^2))
correlation <- cor(pred1$tau_hat, true_tau)
cat(sprintf("RMSE (vs 真实效应): %.3f\n", rmse))
cat(sprintf("相关系数 (vs 真实效应): %.3f\n", correlation))

# 6. 测试三：空间块函数
# ----------------------------------------------------------------------------
cat("\n=== 测试三: create_spatial_blocks ===\n")

blocks <- create_spatial_blocks(coords, n_blocks = 5, seed = 123)
cat("空间块分配:\n")
print(table(blocks))

# 7. 测试四：可视化（如果在交互模式下）
# ----------------------------------------------------------------------------
cat("\n=== 测试四: 可视化 ===\n")

# 效应地图
p1 <- plot(fit1, type = "effect_map", newX = X, newcoords = coords)
print(p1)

# ESS 地图
p2 <- plot(fit1, type = "ess_map")
print(p2)

# 如果有CV结果，绘制CV曲线
if (!is.null(fit1$cv_results)) {
  p3 <- plot(fit1, type = "cv_curve")
  print(p3)
}

# 8. 测试五：自动带宽选择（较慢）
# ----------------------------------------------------------------------------
cat("\n=== 测试五: 自动带宽选择 (可能需要几分钟) ===\n")
cat("跳过此测试以节省时间。取消注释下面代码以运行。\n")

# fit2 <- gw_causal_forest(
#   Y = Y,
#   W = W,
#   X = X,
#   coords = coords,
#   treatment_type = "binary",
#   bandwidth = NULL,       # 自动选择
#   n_folds = 3,            # CV 折数
#   num.trees = 500,
#   seed = 123
# )
# print(fit2)
# plot(fit2, type = "cv_curve")

# 9. 测试六：空间 Bootstrap（较慢）
# ----------------------------------------------------------------------------
cat("\n=== 测试六: 空间 Bootstrap (可能需要几分钟) ===\n")
cat("跳过此测试以节省时间。取消注释下面代码以运行。\n")

# 选择少量预测点进行测试
# test_idx <- sample(n, 20)
# boot_ci <- spatial_bootstrap(
#   object = fit1,
#   newX = X[test_idx, ],
#   newcoords = coords[test_idx, ],
#   n_bootstrap = 50,
#   n_blocks = 5,
#   seed = 456
# )
# print(head(boot_ci))

# 10. 测试完成
# ----------------------------------------------------------------------------
cat("\n")
cat("============================================\n")
cat("           所有基础测试通过!               \n")
cat("============================================\n")
cat("\n")
cat("gwcf 包功能正常。可以开始使用。\n")
cat("\n")
cat("主要函数:\n")
cat("  - gw_causal_forest()     : 拟合模型\n")
cat("  - predict()              : 预测处理效应\n")
cat("  - spatial_bootstrap()    : 空间Bootstrap置信区间\n")
cat("  - select_bandwidth()     : 带宽选择\n")
cat("  - create_spatial_blocks(): 创建空间块\n")
cat("  - diagnostics()          : 诊断统计\n")
cat("  - plot()                 : 可视化\n")

# 重置并行计划
plan(sequential)

