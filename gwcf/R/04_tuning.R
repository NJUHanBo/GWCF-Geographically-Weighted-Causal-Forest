#' Select Optimal Bandwidth via Spatial CV
#'
#' @param Y Outcome vector.
#' @param W Treatment vector.
#' @param X Covariate matrix.
#' @param coords Coordinate matrix.
#' @param treatment_type "binary" or "continuous".
#' @param n_folds Number of spatial folds.
#' @param kernel Kernel type.
#' @param candidates Vector of candidate bandwidths.
#' @param seed Random seed.
#' @param ... Additional arguments passed to gw_causal_forest (e.g. num.trees, honesty).
#'
#' @return List with bandwidth (optimal value) and cv_scores (data frame).
#' @export
select_bandwidth <- function(Y, W, X, coords, 
                             treatment_type = "binary", 
                             n_folds = 5, 
                             kernel = "bisquare",
                             candidates = NULL,
                             seed = NULL,
                             ...) {
  
  if (!is.null(seed)) set.seed(seed)
  
  km <- kmeans(coords, centers = n_folds, nstart = 5)
  folds <- km$cluster
  
  if (is.null(candidates)) {
    n_sub <- min(nrow(coords), 1000)
    idx_sub <- sample(nrow(coords), n_sub)
    dists <- dist(coords[idx_sub, ])
    
    h_min <- as.numeric(quantile(dists, 0.05))
    h_max <- as.numeric(quantile(dists, 0.50))
    candidates <- exp(seq(log(h_min), log(h_max), length.out = 10))
  }
  
  dots <- list(...)
  
  message(sprintf("Running Spatial CV with %d folds for %d bandwidths...", n_folds, length(candidates)))
  
  # ===========================================================================
  # OPTIMIZATION: Pre-compute nuisance models for each fold.
  # Nuisance models (m_hat, e_hat) are INDEPENDENT of bandwidth.
  # Previously: n_folds * n_candidates * 2 = 100 forest trainings
  # Now: n_folds * 2 = 10 forest trainings (10x speedup)
  # ===========================================================================
  
  message("  Pre-computing nuisance models for each fold...")
  
  nuisance_cache <- vector("list", n_folds)
  
  for (k in 1:n_folds) {
    train_idx <- which(folds != k)
    test_idx <- which(folds == k)
    
    # Features for nuisance: X + coords
    X_train_nui <- cbind(X[train_idx, , drop = FALSE], coords[train_idx, , drop = FALSE])
    X_test_nui <- cbind(X[test_idx, , drop = FALSE], coords[test_idx, , drop = FALSE])
    
    # Train nuisance forests (once per fold, not per bandwidth)
    rf_y <- grf::regression_forest(X_train_nui, Y[train_idx], num.trees = 500)
    rf_w <- grf::regression_forest(X_train_nui, W[train_idx], num.trees = 500)
    
    # Predict on test set
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
  
  # ===========================================================================
  # Bandwidth evaluation loop
  # ===========================================================================
  
  cv_scores <- numeric(length(candidates))
  
  for (i in seq_along(candidates)) {
    bw <- candidates[i]
    total_error <- 0
    
    for (k in 1:n_folds) {
      cache <- nuisance_cache[[k]]
      train_idx <- cache$train_idx
      test_idx <- cache$test_idx
      
      # Fit GWCF with this bandwidth
      args <- list(
        Y = Y[train_idx],
        W = W[train_idx],
        X = X[train_idx, , drop = FALSE],
        coords = coords[train_idx, , drop = FALSE],
        treatment_type = treatment_type,
        bandwidth = bw,
        kernel = kernel,
        n_anchors = max(5, floor(sqrt(length(train_idx)))),
        store_data = FALSE,  # Don't store for CV
        seed = seed
      )
      
      final_args <- c(args, dots)
      
      fit <- tryCatch({
        do.call(gw_causal_forest, final_args)
      }, error = function(e) return(NULL))
      
      if (is.null(fit)) {
        total_error <- Inf
        break
      }
      
      preds <- predict(fit, newX = X[test_idx, , drop = FALSE], 
                       newcoords = coords[test_idx, , drop = FALSE],
                       estimate_variance = FALSE)
      tau_hat <- preds$tau_hat
      
      # R-loss using cached nuisance residuals
      loss <- mean((cache$resid_y - tau_hat * cache$resid_w)^2)
      total_error <- total_error + loss
    }
    
    cv_scores[i] <- total_error / n_folds
    message(sprintf("  Bandwidth %.4f: CV Score = %.4f", bw, cv_scores[i]))
  }
  
  best_idx <- which.min(cv_scores)
  
  return(list(
    bandwidth = candidates[best_idx],
    cv_scores = data.frame(bandwidth = candidates, score = cv_scores)
  ))
}
