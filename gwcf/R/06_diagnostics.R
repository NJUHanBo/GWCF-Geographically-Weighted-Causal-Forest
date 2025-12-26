#' Print gw_causal_forest object
#' 
#' @param x A gw_causal_forest object.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns x.
#' @export
print.gw_causal_forest <- function(x, ...) {
  cat("Geographically Weighted Causal Forest\n")
  cat("-------------------------------------\n")
  cat(sprintf("Treatment type: %s\n", x$treatment_type))
  cat(sprintf("Number of anchors: %d\n", length(x$anchor_forests)))
  cat(sprintf("Bandwidth: %.4f (%s kernel)\n", x$bandwidth, x$kernel))
  cat(sprintf("Training observations: %d\n", x$n))
  cat(sprintf("ESS range at anchors: [%.1f, %.1f]\n", 
              min(x$ess_anchors), max(x$ess_anchors)))
  if (!is.null(x$cv_results)) {
    cat(sprintf("Spatial CV performed. Best score: %.4f\n", min(x$cv_results$score)))
  }
  invisible(x)
}

#' Summary method for gw_causal_forest
#' 
#' Provides comprehensive summary including global ATE estimation.
#' 
#' @param object A gw_causal_forest object.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns a list with summary statistics.
#' @export
summary.gw_causal_forest <- function(object, ...) {
  cat("Geographically Weighted Causal Forest - Summary\n")
  cat("================================================\n\n")
  
  # Basic info
  cat("Model Configuration:\n")
  cat(sprintf("  Treatment type: %s\n", object$treatment_type))
  cat(sprintf("  Number of anchors: %d\n", length(object$anchor_forests)))
  cat(sprintf("  Bandwidth: %.4f (%s kernel)\n", object$bandwidth, object$kernel))
  cat(sprintf("  Training observations: %d\n", object$n))
  cat("\n")
  
  # ESS distribution
  cat("Effective Sample Size (ESS) at Anchors:\n")
  ess_stats <- summary(object$ess_anchors)
  cat(sprintf("  Min: %.1f, Median: %.1f, Max: %.1f\n", 
              ess_stats["Min."], ess_stats["Median"], ess_stats["Max."]))
  low_ess <- sum(object$ess_anchors < 50)
  if (low_ess > 0) {
    cat(sprintf("  WARNING: %d anchors with ESS < 50\n", low_ess))
  }
  cat("\n")
  
  # Global ATE estimation
  cat("Global Average Treatment Effect (ATE):\n")
  if (!is.null(object$X_train) && !is.null(object$coords_train)) {
    # Predict at training locations
    preds <- predict(object, newX = object$X_train, newcoords = object$coords_train,
                     estimate_variance = TRUE)
    
    global_ate <- mean(preds$tau_hat)
    # SE of mean using delta method: SE(mean) = sqrt(mean(var) / n)
    global_se <- sqrt(mean(preds$se^2) / object$n)
    z <- qnorm(0.975)
    
    cat(sprintf("  ATE: %.4f (SE: %.4f)\n", global_ate, global_se))
    cat(sprintf("  95%% CI: [%.4f, %.4f]\n", 
                global_ate - z * global_se, global_ate + z * global_se))
    cat(sprintf("  Effect range: [%.4f, %.4f]\n", 
                min(preds$tau_hat), max(preds$tau_hat)))
  } else {
    cat("  (Training data not stored. Refit with store_data = TRUE for ATE.)\n")
  }
  cat("\n")
  
  # CV results
  if (!is.null(object$cv_results)) {
    cat("Spatial Cross-Validation:\n")
    cat(sprintf("  Best bandwidth: %.4f\n", object$bandwidth))
    cat(sprintf("  Best CV score: %.4f\n", min(object$cv_results$score)))
  }
  
  # Return summary list invisibly
  result <- list(
    n = object$n,
    n_anchors = length(object$anchor_forests),
    bandwidth = object$bandwidth,
    kernel = object$kernel,
    treatment_type = object$treatment_type,
    ess_summary = ess_stats
  )
  
  if (!is.null(object$X_train) && !is.null(object$coords_train)) {
    result$global_ate <- global_ate
    result$global_se <- global_se
    result$effect_range <- range(preds$tau_hat)
  }
  
  invisible(result)
}

#' Diagnostics for gw_causal_forest
#' 
#' Computes diagnostic statistics for model quality assessment.
#' 
#' @param object A gw_causal_forest object.
#' @param ess_threshold Numeric, ESS threshold for warnings (default 50).
#' @return List containing:
#'   \item{ess_summary}{Summary statistics of ESS across anchors}
#'   \item{ess_warnings}{Data frame of anchors with ESS below threshold}
#'   \item{overlap_summary}{Treatment variable distribution summary (if training data stored)}
#'   \item{cv_comparison}{CV results summary (if bandwidth was auto-selected)}
#'   \item{bandwidth}{Selected bandwidth value}
#' @export
diagnostics <- function(object, ess_threshold = 50) {
  # ESS summary
  ess_summary <- summary(object$ess_anchors)
  

  # ESS warnings: anchors below threshold
  low_ess_idx <- which(object$ess_anchors < ess_threshold)
  if (length(low_ess_idx) > 0) {
    ess_warnings <- data.frame(
      anchor_id = low_ess_idx,
      x = object$anchor_coords[low_ess_idx, 1],
      y = object$anchor_coords[low_ess_idx, 2],
      ess = object$ess_anchors[low_ess_idx]
    )
  } else {
    ess_warnings <- NULL
  }
  
  # Overlap summary (requires stored training data)
  overlap_summary <- NULL
  if (!is.null(object$W_train)) {
    W <- object$W_train
    if (object$treatment_type == "binary") {
      prop_treated <- mean(W)
      overlap_summary <- list(
        type = "binary",
        n_treated = sum(W == 1),
        n_control = sum(W == 0),
        prop_treated = prop_treated,
        balance_warning = prop_treated < 0.1 || prop_treated > 0.9
      )
    } else {
      overlap_summary <- list(
        type = "continuous",
        range = range(W),
        mean = mean(W),
        sd = sd(W),
        quantiles = quantile(W, c(0.05, 0.25, 0.5, 0.75, 0.95))
      )
    }
  }
  
  # CV comparison (if CV was performed)
  cv_comparison <- NULL
  if (!is.null(object$cv_results)) {
    cv_comparison <- list(
      best_bandwidth = object$bandwidth,
      best_score = min(object$cv_results$score),
      n_candidates = nrow(object$cv_results),
      bandwidth_range = range(object$cv_results$bandwidth)
    )
  }
  
  result <- list(
    ess_summary = ess_summary,
    ess_warnings = ess_warnings,
    overlap_summary = overlap_summary,
    cv_comparison = cv_comparison,
    bandwidth = object$bandwidth,
    n_anchors = length(object$anchor_forests),
    n_low_ess = length(low_ess_idx)
  )
  
  class(result) <- "gwcf_diagnostics"
  return(result)
}

#' Print method for gwcf_diagnostics
#' 
#' @param x A gwcf_diagnostics object.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns x.
#' @export
print.gwcf_diagnostics <- function(x, ...) {
  cat("GWCF Diagnostics\n")
  cat("================\n\n")
  
  cat("Bandwidth:", x$bandwidth, "\n")
  cat("Number of anchors:", x$n_anchors, "\n\n")
  
  cat("ESS Summary:\n")
  print(x$ess_summary)
  cat("\n")
  
  if (x$n_low_ess > 0) {
    cat(sprintf("WARNING: %d anchors with ESS < 50\n", x$n_low_ess))
    if (!is.null(x$ess_warnings)) {
      cat("Low ESS anchors:\n")
      print(x$ess_warnings)
    }
  } else {
    cat("All anchors have adequate ESS (>= 50)\n")
  }
  cat("\n")
  
  if (!is.null(x$overlap_summary)) {
    cat("Overlap Summary:\n")
    if (x$overlap_summary$type == "binary") {
      cat(sprintf("  Treated: %d, Control: %d (%.1f%% treated)\n",
                  x$overlap_summary$n_treated, 
                  x$overlap_summary$n_control,
                  x$overlap_summary$prop_treated * 100))
      if (x$overlap_summary$balance_warning) {
        cat("  WARNING: Severe imbalance in treatment assignment\n")
      }
    } else {
      cat(sprintf("  Treatment range: [%.3f, %.3f]\n", 
                  x$overlap_summary$range[1], x$overlap_summary$range[2]))
      cat(sprintf("  Mean: %.3f, SD: %.3f\n", 
                  x$overlap_summary$mean, x$overlap_summary$sd))
    }
  }
  cat("\n")
  
  if (!is.null(x$cv_comparison)) {
    cat("Cross-Validation:\n")
    cat(sprintf("  Best bandwidth: %.4f (score: %.4f)\n", 
                x$cv_comparison$best_bandwidth, x$cv_comparison$best_score))
    cat(sprintf("  Searched %d candidates in range [%.4f, %.4f]\n",
                x$cv_comparison$n_candidates,
                x$cv_comparison$bandwidth_range[1],
                x$cv_comparison$bandwidth_range[2]))
  }
  
  invisible(x)
}

#' @importFrom ggplot2 geom_histogram geom_density geom_vline scale_fill_manual
NULL

# Suppress R CMD check NOTEs for ggplot2 NSE variables
utils::globalVariables(c("bandwidth", "score", "x", "y", "ess", "tau", "W", "density"))

#' Plot method for gw_causal_forest
#' 
#' @param x A gw_causal_forest object.
#' @param type Character, plot type: "effect_map", "ess_map", "cv_curve", or "overlap".
#' @param newX Matrix of covariates (required for effect_map).
#' @param newcoords Matrix of coordinates (required for effect_map).
#' @param ... Additional arguments (unused).
#' @return A ggplot2 object.
#' @export
plot.gw_causal_forest <- function(x, type = c("effect_map", "ess_map", "cv_curve", "overlap"), 
                                   newX = NULL, newcoords = NULL, ...) {
  type <- match.arg(type)
  
  if (type == "cv_curve") {
    if (is.null(x$cv_results)) stop("No CV results available.")
    
    ggplot(x$cv_results, aes(x = bandwidth, y = score)) +
      geom_line() + geom_point() +
      scale_x_log10() +
      theme_minimal() +
      labs(title = "Spatial Cross-Validation", x = "Bandwidth (log scale)", y = "CV Error")
      
  } else if (type == "ess_map") {
    df <- data.frame(
      x = x$anchor_coords[, 1],
      y = x$anchor_coords[, 2],
      ess = x$ess_anchors
    )
    
    ggplot(df, aes(x = x, y = y, color = ess)) +
      geom_point(size = 3) +
      scale_color_viridis_c() +
      theme_minimal() +
      labs(title = "Effective Sample Size at Anchors")
      
  } else if (type == "effect_map") {
    if (is.null(newX) || is.null(newcoords)) stop("newX and newcoords required for effect map.")
    
    preds <- predict(x, newX, newcoords)
    df <- data.frame(
      x = newcoords[, 1],
      y = newcoords[, 2],
      tau = preds$tau_hat
    )
    
    ggplot(df, aes(x = x, y = y, color = tau)) +
      geom_point() +
      scale_color_viridis_c() +
      theme_minimal() +
      labs(title = "Estimated Treatment Effects")
      
  } else if (type == "overlap") {
    if (is.null(x$W_train)) {
      stop("Overlap plot requires training data. Refit model with store_data = TRUE.")
    }
    
    W <- x$W_train
    
    if (x$treatment_type == "binary") {
      df <- data.frame(
        W = factor(W, levels = c(0, 1), labels = c("Control", "Treated"))
      )
      
      ggplot(df, aes(x = W, fill = W)) +
        geom_bar(alpha = 0.8) +
        scale_fill_manual(values = c("Control" = "#2C3E50", "Treated" = "#E74C3C")) +
        theme_minimal() +
        labs(title = "Treatment Assignment Distribution",
             x = "Treatment Group", y = "Count") +
        theme(legend.position = "none")
        
    } else {
      df <- data.frame(W = W)
      
      ggplot(df, aes(x = W)) +
        geom_histogram(aes(y = after_stat(density)), bins = 30, 
                       fill = "#3498DB", alpha = 0.7) +
        geom_density(color = "#2C3E50", linewidth = 1) +
        geom_vline(xintercept = mean(W), color = "#E74C3C", 
                   linetype = "dashed", linewidth = 1) +
        theme_minimal() +
        labs(title = "Treatment Variable Distribution",
             x = "Treatment Value", y = "Density",
             caption = sprintf("Mean: %.3f, SD: %.3f", mean(W), sd(W)))
    }
  }
}

