#' Select Optimal Bandwidth via Spatial CV
#'
#' Uses R-loss with pre-computed nuisance models and lightweight anchor fitting.
#' Does NOT recursively call gw_causal_forest() — fits anchor forests directly.
#'
#' @param Y Outcome vector.
#' @param W Treatment vector.
#' @param X Covariate matrix.
#' @param coords Coordinate matrix.
#' @param treatment_type "binary" or "continuous".
#' @param n_folds Number of spatial folds.
#' @param kernel Kernel type.
#' @param n_anchors Number of anchor points (passed from parent).
#' @param candidates Vector of candidate bandwidths.
#' @param seed Random seed.
#' @param ... Additional arguments (num.trees, honesty, min.node.size).
#'
#' @return List with bandwidth (optimal value) and cv_scores (data frame).
#' @export
select_bandwidth <- function(Y, W, X, coords, 
                             treatment_type = "binary", 
                             n_folds = 5, 
                             kernel = "bisquare",
                             n_anchors = NULL,
                             candidates = NULL,
                             seed = NULL,
                             ...) {
  
  if (!is.null(seed)) set.seed(seed)
  
  dots <- list(...)
  num.trees_cv <- min(dots$num.trees %||% 500, 500)
  honesty_cv <- dots$honesty %||% TRUE
  min.node.size_cv <- dots$min.node.size %||% 5
  
  # Spatial folds via k-means
  km <- kmeans(coords, centers = n_folds, nstart = 5)
  folds <- km$cluster
  
  # Default anchor count for CV: use user's value or heuristic on training size
  if (is.null(n_anchors)) {
    n_train_approx <- floor(nrow(coords) * (n_folds - 1) / n_folds)
    n_anchors <- max(10, floor(sqrt(n_train_approx)))
  }
  # CV uses fewer anchors for speed (at most user's value, at least 10)
  n_anchors_cv <- max(10, min(n_anchors, 30))
  
  # Candidate bandwidths from distance distribution
  if (is.null(candidates)) {
    n_sub <- min(nrow(coords), 1000)
    idx_sub <- sample(nrow(coords), n_sub)
    dists <- dist(coords[idx_sub, ])
    
    h_min <- as.numeric(quantile(dists, 0.05))
    h_max <- as.numeric(quantile(dists, 0.50))
    candidates <- exp(seq(log(h_min), log(h_max), length.out = 8))
  }
  
  message(sprintf("Running Spatial CV with %d folds for %d bandwidths...", 
                  n_folds, length(candidates)))
  
  # Pre-compute nuisance models (independent of bandwidth)
  message("  Pre-computing nuisance models for each fold...")
  
  nuisance_cache <- vector("list", n_folds)
  
  for (k in 1:n_folds) {
    train_idx <- which(folds != k)
    test_idx <- which(folds == k)
    
    X_train_nui <- cbind(X[train_idx, , drop = FALSE], coords[train_idx, , drop = FALSE])
    X_test_nui <- cbind(X[test_idx, , drop = FALSE], coords[test_idx, , drop = FALSE])
    
    rf_y <- grf::regression_forest(X_train_nui, Y[train_idx], num.trees = 500)
    rf_w <- grf::regression_forest(X_train_nui, W[train_idx], num.trees = 500)
    
    m_hat <- predict(rf_y, X_test_nui)$predictions
    e_hat <- predict(rf_w, X_test_nui)$predictions
    
    nuisance_cache[[k]] <- list(
      train_idx = train_idx,
      test_idx = test_idx,
      resid_y = Y[test_idx] - m_hat,
      resid_w = W[test_idx] - e_hat
    )
  }
  
  message("  Nuisance models computed. Evaluating bandwidths...")
  
  # Bandwidth evaluation: lightweight direct fitting (no recursive gw_causal_forest call)
  cv_scores <- numeric(length(candidates))
  
  for (i in seq_along(candidates)) {
    bw <- candidates[i]
    fold_losses <- numeric(n_folds)
    all_valid <- TRUE
    
    for (k in 1:n_folds) {
      cache <- nuisance_cache[[k]]
      train_idx <- cache$train_idx
      test_idx <- cache$test_idx
      
      # Select anchors from training data
      train_coords <- coords[train_idx, , drop = FALSE]
      anchors_cv <- select_anchors(train_coords, n_anchors_cv, seed = seed)
      
      # Compute weights for training data at anchors
      dists_train <- compute_distance(train_coords, anchors_cv)
      wts_train <- get_weights(dists_train, bw, kernel)
      
      # Fit local forests at each anchor (lightweight, sequential)
      anchor_forests <- vector("list", ncol(wts_train))
      for (j in seq_len(ncol(wts_train))) {
        w_j <- wts_train[, j]
        idx <- which(w_j > 0)
        if (length(idx) < min.node.size_cv * 4) next
        
        anchor_forests[[j]] <- tryCatch(
          grf::causal_forest(
            X[train_idx[idx], , drop = FALSE],
            Y[train_idx[idx]],
            W[train_idx[idx]],
            sample.weights = w_j[idx],
            num.trees = num.trees_cv,
            honesty = honesty_cv,
            min.node.size = min.node.size_cv,
            ci.group.size = 2
          ),
          error = function(e) NULL
        )
      }
      
      valid_j <- which(!sapply(anchor_forests, is.null))
      if (length(valid_j) == 0) {
        all_valid <- FALSE
        break
      }
      
      # Predict on test set: interpolate across valid anchors
      test_coords <- coords[test_idx, , drop = FALSE]
      dists_test <- compute_distance(test_coords, anchors_cv[valid_j, , drop = FALSE])
      wts_test <- get_weights(dists_test, bw, kernel)
      
      # Fallback for points outside kernel support
      row_sums <- rowSums(wts_test)
      zero_idx <- which(row_sums == 0)
      for (zi in zero_idx) {
        nearest <- which.min(dists_test[zi, ])
        wts_test[zi, nearest] <- 1
        row_sums[zi] <- 1
      }
      wts_test <- wts_test / row_sums
      
      # Collect predictions from each anchor
      X_test <- X[test_idx, , drop = FALSE]
      pred_mat <- matrix(0, nrow = length(test_idx), ncol = length(valid_j))
      for (jj in seq_along(valid_j)) {
        pred_mat[, jj] <- predict(anchor_forests[[valid_j[jj]]], X_test)$predictions
      }
      
      tau_hat <- rowSums(wts_test * pred_mat)
      
      # R-loss
      fold_losses[k] <- mean((cache$resid_y - tau_hat * cache$resid_w)^2)
    }
    
    if (!all_valid || any(is.na(fold_losses))) {
      cv_scores[i] <- Inf
    } else {
      cv_scores[i] <- mean(fold_losses)
    }
    message(sprintf("  Bandwidth %.1f: CV Score = %.4f", bw, cv_scores[i]))
  }
  
  # Select best (finite) bandwidth
  finite_idx <- which(is.finite(cv_scores))
  if (length(finite_idx) == 0) {
    warning("All bandwidths produced invalid CV scores. Using median candidate.")
    best_idx <- ceiling(length(candidates) / 2)
  } else {
    best_idx <- finite_idx[which.min(cv_scores[finite_idx])]
  }
  
  return(list(
    bandwidth = candidates[best_idx],
    cv_scores = data.frame(bandwidth = candidates, score = cv_scores)
  ))
}

# Null-coalescing operator (not exported)
`%||%` <- function(x, y) if (is.null(x)) y else x
