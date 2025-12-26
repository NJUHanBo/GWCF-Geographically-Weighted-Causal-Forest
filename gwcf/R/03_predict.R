#' @importFrom stats qnorm
NULL

#' Predict method for gw_causal_forest
#'
#' @param object A gw_causal_forest object.
#' @param newX Matrix of new covariates.
#' @param newcoords Matrix of new coordinates.
#' @param estimate_variance Logical, whether to estimate variance for confidence intervals.
#' @param alpha Significance level for confidence intervals (default 0.05 for 95% CI).
#' @param ... Additional arguments.
#'
#' @return A data frame with columns: tau_hat (point estimates), ess (effective sample size),
#'   and if estimate_variance=TRUE: se (standard errors), ci_lower, ci_upper (confidence bounds).
#' @export
predict.gw_causal_forest <- function(object, newX, newcoords, 
                                     estimate_variance = TRUE, 
                                     alpha = 0.05, ...) {
  if (missing(newX) || missing(newcoords)) {
    stop("newX and newcoords must be provided.")
  }
  
  newX <- as.matrix(newX)
  newcoords <- as.matrix(newcoords)
  
  if (nrow(newX) != nrow(newcoords)) {
    stop("newX and newcoords must have the same number of rows.")
  }
  
  # 1. Compute interpolation weights
  dists <- compute_distance(newcoords, object$anchor_coords)
  weights_interp <- get_weights(dists, object$bandwidth, object$kernel)
  row_sums <- rowSums(weights_interp)
  
  # Handle points outside kernel support
  zero_sum_idx <- which(row_sums == 0)
  n_outside <- length(zero_sum_idx)
  
  if (n_outside > 0) {
    # Check how far outside these points are
    min_dists <- apply(dists[zero_sum_idx, , drop = FALSE], 1, min)
    max_reasonable_dist <- object$bandwidth * 3  # 3x bandwidth is very far
    
    # Points extremely far from any anchor (likely outside study area)
    extreme_idx <- zero_sum_idx[min_dists > max_reasonable_dist]
    
    if (length(extreme_idx) > 0) {
      warning(sprintf(
        "%d prediction point(s) are far outside the study area (>3x bandwidth from nearest anchor). ",
        length(extreme_idx)),
        "Predictions for these points are extrapolations and may be unreliable.",
        call. = FALSE)
    } else if (n_outside > nrow(newcoords) * 0.1) {
      # More than 10% of points are outside kernel support
      warning(sprintf(
        "%d prediction point(s) (%.1f%%) fall outside kernel support. ",
        n_outside, 100 * n_outside / nrow(newcoords)),
        "Consider using a larger bandwidth.",
        call. = FALSE)
    }
    
    # Fallback: assign to nearest anchor
    for (i in zero_sum_idx) {
      nearest <- which.min(dists[i, ])
      weights_interp[i, nearest] <- 1
      row_sums[i] <- 1
    }
  }
  weights_interp <- weights_interp / row_sums
  
  # 2. Get predictions from each anchor forest
  n_anchors <- length(object$anchor_forests)
  n_targets <- nrow(newX)
  pred_mat <- matrix(NA, nrow = n_targets, ncol = n_anchors)
  var_mat <- matrix(NA, nrow = n_targets, ncol = n_anchors)
  
  preds_list <- future.apply::future_lapply(1:n_anchors, function(j) {
    forest <- object$anchor_forests[[j]]
    raw_pred <- predict(forest, newX, estimate.variance = estimate_variance)
    
    preds <- as.vector(raw_pred$predictions)
    vars <- if (estimate_variance && !is.null(raw_pred$variance.estimates)) {
      as.vector(raw_pred$variance.estimates)
    } else {
      rep(NA, length(preds))
    }
    
    list(predictions = preds, variances = vars)
  }, future.seed = TRUE)
  
  for (j in 1:n_anchors) {
    pred_mat[, j] <- preds_list[[j]]$predictions
    var_mat[, j] <- preds_list[[j]]$variances
  }
  
  # 3. Compute weighted average of predictions
  tau_hat <- rowSums(weights_interp * pred_mat)
  ess <- apply(weights_interp, 1, compute_ess)
  
  # 4. Compute variance of weighted average
  # Var(sum(w_j * tau_j)) = sum(w_j^2 * Var(tau_j)) assuming independence
  if (estimate_variance && !all(is.na(var_mat))) {
    var_tau <- rowSums((weights_interp^2) * var_mat)
    se <- sqrt(var_tau)
    z <- qnorm(1 - alpha / 2)
    ci_lower <- tau_hat - z * se
    ci_upper <- tau_hat + z * se
    
    res <- data.frame(
      tau_hat = tau_hat,
      se = se,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      ess = ess
    )
  } else {
    res <- data.frame(
      tau_hat = tau_hat,
      ess = ess
    )
  }
  
  return(res)
}
