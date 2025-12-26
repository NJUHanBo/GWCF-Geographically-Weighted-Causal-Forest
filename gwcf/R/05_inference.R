#' Spatial Bootstrap for Confidence Intervals
#'
#' Computes confidence intervals via spatial block bootstrap. Requires the model
#' to be fitted with \code{store_data = TRUE} (default).
#'
#' @param object Fitted gw_causal_forest object.
#' @param newX Covariates for prediction.
#' @param newcoords Coordinates for prediction.
#' @param n_bootstrap Number of bootstrap iterations (default 200).
#' @param n_blocks Number of spatial blocks for resampling.
#' @param conf_level Confidence level (default 0.95 for 95% CI).
#' @param seed Random seed.
#'
#' @return Data frame with columns: ci_lower, ci_upper, se.
#' @export
spatial_bootstrap <- function(object, newX, newcoords, 
                              n_bootstrap = 200, 
                              n_blocks = 10, 
                              conf_level = 0.95,
                              seed = NULL) {
  
  # Check if training data is stored
  if (is.null(object$Y_train) || is.null(object$W_train) || 
      is.null(object$X_train) || is.null(object$coords_train)) {
    stop("Training data not stored in model object. ",
         "Refit with store_data = TRUE, or use spatial_bootstrap_with_data().")
  }
  
  # Delegate to spatial_bootstrap_with_data
  spatial_bootstrap_with_data(
    object = object,
    Y_train = object$Y_train,
    W_train = object$W_train,
    X_train = object$X_train,
    coords_train = object$coords_train,
    newX = newX,
    newcoords = newcoords,
    n_bootstrap = n_bootstrap,
    n_blocks = n_blocks,
    conf_level = conf_level,
    seed = seed
  )
}

#' Spatial Bootstrap with Data provided
#' 
#' @param object Fitted model.
#' @param Y_train Training outcome.
#' @param W_train Training treatment.
#' @param X_train Training covariates.
#' @param coords_train Training coordinates.
#' @param newX Prediction covariates.
#' @param newcoords Prediction coordinates.
#' @param n_bootstrap Number of bootstrap iterations.
#' @param n_blocks Number of spatial blocks.
#' @param conf_level Confidence level (default 0.95).
#' @param seed Random seed.
#' @param verbose Logical, whether to print progress.
#' 
#' @export
spatial_bootstrap_with_data <- function(object, Y_train, W_train, X_train, coords_train,
                                        newX, newcoords,
                                        n_bootstrap = 200,
                                        n_blocks = 10,
                                        conf_level = 0.95,
                                        seed = NULL,
                                        verbose = TRUE) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Convert conf_level to alpha for quantile calculation
  alpha <- 1 - conf_level
  
  # Create spatial blocks
  km <- kmeans(coords_train, centers = n_blocks, nstart = 5, iter.max = 50)
  blocks <- km$cluster
  block_ids <- sort(unique(blocks))
  
  # Extract parameters from original model
  treatment_type <- object$treatment_type
  bandwidth <- object$bandwidth
  n_anchors <- length(object$anchor_forests)
  kernel <- object$kernel
  num_trees <- if (!is.null(object$num.trees)) object$num.trees else 500
  
  # Bootstrap worker function
  bootstrap_one <- function(b) {
    # Resample blocks
    sampled_blocks <- sample(block_ids, replace = TRUE)
    
    # Construct indices
    boot_idx <- unlist(lapply(sampled_blocks, function(bid) which(blocks == bid)))
    
    # Refit model
    fit_b <- gw_causal_forest(
      Y = Y_train[boot_idx],
      W = W_train[boot_idx],
      X = X_train[boot_idx, , drop=FALSE],
      coords = coords_train[boot_idx, , drop=FALSE],
      treatment_type = treatment_type,
      bandwidth = bandwidth,
      n_anchors = n_anchors,
      kernel = kernel,
      num.trees = min(num_trees, 500),  # Cap for speed
      honesty = TRUE,
      store_data = FALSE,  # Don't store data in bootstrap samples
      seed = NULL
    )
    
    # Predict
    p <- predict(fit_b, newX, newcoords, estimate_variance = FALSE)$tau_hat
    return(p)
  }
  
  # Run bootstrap iterations sequentially (more reliable than parallel for dev)
  boot_preds <- vector("list", n_bootstrap)
  
  for (b in seq_len(n_bootstrap)) {
    if (verbose && b %% 10 == 0) {
      message(sprintf("  Bootstrap iteration %d/%d", b, n_bootstrap))
    }
    boot_preds[[b]] <- bootstrap_one(b)
  }
  
  # Aggregate results
  pred_mat <- do.call(cbind, boot_preds)
  
  ci_lower <- apply(pred_mat, 1, quantile, probs = alpha / 2)
  ci_upper <- apply(pred_mat, 1, quantile, probs = 1 - alpha / 2)
  se <- apply(pred_mat, 1, sd)
  
  return(data.frame(
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    se = se
  ))
}

