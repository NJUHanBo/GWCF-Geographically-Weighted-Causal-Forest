#' Geographically Weighted Causal Forest
#'
#' @param Y Numeric vector of outcomes.
#' @param W Numeric vector of treatment assignment (binary or continuous).
#' @param X Matrix of covariates.
#' @param coords Matrix of spatial coordinates.
#' @param treatment_type Character, either "binary" or "continuous".
#' @param bandwidth Numeric bandwidth. If NULL, selected via cross-validation.
#' @param n_anchors Integer, number of anchor points. If NULL, heuristic is used.
#' @param kernel Character, kernel function ("bisquare", "gaussian", "exponential").
#' @param n_folds Integer, number of folds for bandwidth selection.
#' @param num.trees Integer, number of trees per forest.
#' @param honesty Logical, whether to use honest splitting.
#' @param min.node.size Integer, minimum node size.
#' @param store_data Logical, whether to store training data for bootstrap (default TRUE).
#' @param seed Integer, random seed.
#'
#' @return An object of class \code{gw_causal_forest}.
#' @export
gw_causal_forest <- function(Y, W, X, coords,
                             treatment_type = c("binary", "continuous"),
                             bandwidth = NULL,
                             n_anchors = NULL,
                             kernel = "bisquare",
                             n_folds = 5,
                             num.trees = 2000,
                             honesty = TRUE,
                             min.node.size = 5,
                             store_data = TRUE,
                             seed = NULL) {
  
  treatment_type <- match.arg(treatment_type)
  
  if (nrow(X) != length(Y) || nrow(X) != length(W) || nrow(coords) != length(Y)) {
    stop("Dimensions of Y, W, X, and coords must match.")
  }
  if (any(is.na(Y)) || any(is.na(W)) || any(is.na(X)) || any(is.na(coords))) {
    stop("Missing values are not allowed.")
  }
  
  X <- as.matrix(X)
  coords <- as.matrix(coords)
  
  if (!is.null(seed)) set.seed(seed)
  
  # 1. Bandwidth Selection
  if (is.null(bandwidth)) {
    message("Selecting bandwidth via spatial cross-validation...")
    bw_res <- select_bandwidth(Y, W, X, coords, 
                               treatment_type = treatment_type,
                               n_folds = n_folds, kernel = kernel,
                               num.trees = num.trees, 
                               honesty = honesty,
                               min.node.size = min.node.size,
                               seed = seed)
    bandwidth <- bw_res$bandwidth
    cv_results <- bw_res$cv_scores
    message(sprintf("Selected bandwidth: %f", bandwidth))
  } else {
    cv_results <- NULL
  }
  
  # 2. Anchor Selection
  if (is.null(n_anchors)) {
    n_anchors <- max(10, floor(sqrt(nrow(X))))
  }
  anchors <- select_anchors(coords, n_anchors, seed = seed)
  
  # 3. Compute Weights
  message(sprintf("Computing weights for %d anchors...", nrow(anchors)))
  dists_to_anchors <- compute_distance(coords, anchors)
  weights_mat <- get_weights(dists_to_anchors, bandwidth, kernel)
  
  ess_anchors <- apply(weights_mat, 2, compute_ess)
  if (any(ess_anchors < 20)) {
    warning("Some anchors have very low ESS (< 20). Consider increasing bandwidth.")
  }
  
  # 4. Fit Local Forests (Unified Approach)
  # grf::causal_forest handles both binary and continuous W.
  # For continuous W, it implements the R-learner (orthogonalized residuals).
  # This replaces the manual orthogonalization step, reducing code complexity and errors.
  
  message("Fitting local forests...")
  weights_list <- as.list(as.data.frame(weights_mat))
  
  fit_anchor_worker <- function(w_i, X_local, Y_local, W_local, n_trees, h_val, mns) {
    subset_idx <- which(w_i > 0)
    if (length(subset_idx) < mns * 2) return(NULL) 
    
    X_sub <- X_local[subset_idx, , drop = FALSE]
    Y_sub <- Y_local[subset_idx]
    W_sub <- W_local[subset_idx]
    w_sub <- w_i[subset_idx]
    
    # causal_forest works for both binary and continuous W
    # ci.group.size >= 2 is required for variance estimation
    fit <- grf::causal_forest(X_sub, Y_sub, W_sub,
                              sample.weights = w_sub,
                              num.trees = n_trees,
                              honesty = h_val,
                              min.node.size = mns,
                              ci.group.size = 2)
    return(fit)
  }
  
  forest_list <- future.apply::future_lapply(weights_list, function(w) {
    fit_anchor_worker(w, X, Y, W, num.trees, honesty, min.node.size)
  }, future.seed = TRUE)
  
  valid_anchors_idx <- which(!sapply(forest_list, is.null))
  if (length(valid_anchors_idx) == 0) {
    stop("All anchor forests failed to fit.")
  }
  
  forest_list <- forest_list[valid_anchors_idx]
  anchors <- anchors[valid_anchors_idx, , drop = FALSE]
  ess_anchors <- ess_anchors[valid_anchors_idx]
  
  res <- list(
    anchor_forests = forest_list,
    anchor_coords = anchors,
    bandwidth = bandwidth,
    kernel = kernel,
    treatment_type = treatment_type,
    cv_results = cv_results,
    ess_anchors = ess_anchors,
    X_train_summary = apply(X, 2, range),
    n = length(Y),
    num.trees = num.trees,
    honesty = honesty,
    min.node.size = min.node.size,
    call = match.call()
  )
  
  # Store training data if requested (needed for spatial_bootstrap)
  if (store_data) {
    res$Y_train <- Y
    res$W_train <- W
    res$X_train <- X
    res$coords_train <- coords
  }
  
  class(res) <- "gw_causal_forest"
  return(res)
}
