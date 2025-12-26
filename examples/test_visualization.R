# ============================================================================
# gwcf 可视化功能测试
# ============================================================================

# 加载包
cat("=== 加载 gwcf 包 ===\n")
devtools::load_all("/Users/hanbo/GWCT/gwcf")
library(ggplot2)

# 生成模拟数据
cat("=== 生成模拟数据 ===\n")
set.seed(42)
n <- 400

coords <- cbind(x = runif(n, 0, 1), y = runif(n, 0, 1))
X <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
colnames(X) <- paste0("X", 1:4)

propensity <- plogis(-0.5 + 0.3 * X[,1] + 0.2 * coords[,1])
W <- rbinom(n, 1, propensity)
true_tau <- 1 + 2 * coords[,1] + 1.5 * coords[,2]
Y <- 2 + X[,1] + true_tau * W + rnorm(n, 0, 1)

cat(sprintf("样本量: %d\n\n", n))

# 拟合模型（使用自动带宽选择）
cat("=== 拟合模型 (自动带宽选择) ===\n")
fit <- gw_causal_forest(
  Y = Y,
  W = W,
  X = X,
  coords = coords,
  treatment_type = "binary",
  bandwidth = NULL,  # 自动选择
  n_folds = 3,
  n_anchors = 15,
  num.trees = 500,
  seed = 123
)

print(fit)

# 创建输出目录
output_dir <- "/Users/hanbo/GWCT/examples/plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# 测试一：效应地图
# ============================================================================
cat("\n=== 测试一: 效应地图 (effect_map) ===\n")

p1 <- plot(fit, type = "effect_map", newX = X, newcoords = coords)

# 增强样式
p1_enhanced <- p1 + 
  labs(
    title = "Geographically Weighted Treatment Effects",
    subtitle = sprintf("Bandwidth = %.3f, Anchors = %d", fit$bandwidth, length(fit$anchor_forests)),
    x = "X Coordinate",
    y = "Y Coordinate",
    color = "Treatment\nEffect (τ)"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray50")
  )

ggsave(file.path(output_dir, "effect_map.png"), p1_enhanced, width = 8, height = 6, dpi = 150)
cat("效应地图已保存: plots/effect_map.png\n")

# ============================================================================
# 测试二：ESS 地图
# ============================================================================
cat("\n=== 测试二: ESS 地图 (ess) ===\n")

p2 <- plot(fit, type = "ess_map")

# 增强样式
p2_enhanced <- p2 +
  labs(
    title = "Effective Sample Size at Anchor Points",
    subtitle = sprintf("Total anchors: %d", length(fit$anchor_forests)),
    x = "X Coordinate",
    y = "Y Coordinate",
    color = "ESS"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray50")
  )

ggsave(file.path(output_dir, "ess_map.png"), p2_enhanced, width = 8, height = 6, dpi = 150)
cat("ESS地图已保存: plots/ess_map.png\n")

# ============================================================================
# 测试三：CV 曲线
# ============================================================================
cat("\n=== 测试三: CV 曲线 (cv_curve) ===\n")

if (!is.null(fit$cv_results)) {
  p3 <- plot(fit, type = "cv_curve")
  
  # 增强样式
  best_bw <- fit$bandwidth
  best_score <- min(fit$cv_results$score)
  
  p3_enhanced <- p3 +
    geom_vline(xintercept = best_bw, linetype = "dashed", color = "red", alpha = 0.7) +
    annotate("text", x = best_bw, y = max(fit$cv_results$score), 
             label = sprintf("Best: %.3f", best_bw), 
             hjust = -0.1, color = "red", size = 3.5) +
    labs(
      title = "Spatial Cross-Validation for Bandwidth Selection",
      subtitle = sprintf("Best bandwidth: %.4f (CV Error: %.4f)", best_bw, best_score),
      x = "Bandwidth (log scale)",
      y = "Cross-Validation Error"
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray50")
    )
  
  ggsave(file.path(output_dir, "cv_curve.png"), p3_enhanced, width = 8, height = 5, dpi = 150)
  cat("CV曲线已保存: plots/cv_curve.png\n")
} else {
  cat("无CV结果 (模型使用固定带宽)\n")
}

# ============================================================================
# 测试四：真实效应 vs 估计效应对比
# ============================================================================
cat("\n=== 测试四: 真实效应 vs 估计效应对比 ===\n")

pred <- predict(fit, newX = X, newcoords = coords)

# 创建对比数据框
compare_df <- data.frame(
  x = coords[, 1],
  y = coords[, 2],
  true_tau = true_tau,
  estimated_tau = pred$tau_hat,
  error = pred$tau_hat - true_tau
)

# 4a: 散点对比图
p4a <- ggplot(compare_df, aes(x = true_tau, y = estimated_tau)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    title = "True vs Estimated Treatment Effects",
    subtitle = sprintf("Correlation: %.3f, RMSE: %.3f", 
                      cor(compare_df$true_tau, compare_df$estimated_tau),
                      sqrt(mean(compare_df$error^2))),
    x = "True Treatment Effect",
    y = "Estimated Treatment Effect"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray50")
  )

ggsave(file.path(output_dir, "true_vs_estimated.png"), p4a, width = 7, height = 6, dpi = 150)
cat("对比图已保存: plots/true_vs_estimated.png\n")

# 4b: 误差空间分布
p4b <- ggplot(compare_df, aes(x = x, y = y, color = error)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(
    title = "Spatial Distribution of Estimation Error",
    subtitle = "Blue = Underestimation, Red = Overestimation",
    x = "X Coordinate",
    y = "Y Coordinate",
    color = "Error\n(Est - True)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray50")
  )

ggsave(file.path(output_dir, "error_map.png"), p4b, width = 8, height = 6, dpi = 150)
cat("误差地图已保存: plots/error_map.png\n")

# ============================================================================
# 测试五：置信区间可视化
# ============================================================================
cat("\n=== 测试五: 置信区间可视化 ===\n")

# 选择一条横切线上的点
y_slice <- 0.5
tolerance <- 0.05
slice_idx <- which(abs(coords[,2] - y_slice) < tolerance)
slice_idx <- slice_idx[order(coords[slice_idx, 1])]

if (length(slice_idx) > 10) {
  slice_df <- data.frame(
    x = coords[slice_idx, 1],
    tau_hat = pred$tau_hat[slice_idx],
    ci_lower = pred$ci_lower[slice_idx],
    ci_upper = pred$ci_upper[slice_idx],
    true_tau = true_tau[slice_idx]
  )
  
  p5 <- ggplot(slice_df, aes(x = x)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "steelblue", alpha = 0.3) +
    geom_line(aes(y = tau_hat), color = "steelblue", size = 1) +
    geom_line(aes(y = true_tau), color = "red", linetype = "dashed", size = 1) +
    labs(
      title = sprintf("Treatment Effect Profile (Y ≈ %.1f)", y_slice),
      subtitle = "Blue: Estimated (with 95% CI), Red dashed: True",
      x = "X Coordinate",
      y = "Treatment Effect (τ)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray50")
    )
  
  ggsave(file.path(output_dir, "effect_profile.png"), p5, width = 8, height = 5, dpi = 150)
  cat("效应剖面图已保存: plots/effect_profile.png\n")
}

# ============================================================================
# 完成
# ============================================================================
cat("\n")
cat("============================================\n")
cat("         可视化测试完成!                    \n")
cat("============================================\n")
cat("\n")
cat("所有图片保存在: /Users/hanbo/GWCT/examples/plots/\n")
cat("  - effect_map.png        效应地图\n")
cat("  - ess_map.png           ESS地图\n")
cat("  - cv_curve.png          CV曲线\n")
cat("  - true_vs_estimated.png 真实vs估计对比\n")
cat("  - error_map.png         误差空间分布\n")
cat("  - effect_profile.png    效应剖面图\n")

