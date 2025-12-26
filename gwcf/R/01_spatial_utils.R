#' @importFrom stats dist kmeans quantile predict sd
#' @importFrom methods is
#' @importFrom grf causal_forest regression_forest average_treatment_effect
#' @importFrom future.apply future_lapply
#' @importFrom sf st_distance st_as_sf st_coordinates
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_sf scale_color_viridis_c scale_x_log10 theme_minimal labs
NULL

#' Compute Euclidean distance matrix
#' 
#' @param coords1 Matrix of first set of coordinates (n x 2)
#' @param coords2 Matrix of second set of coordinates (m x 2)
#' @return Matrix of distances (n x m)
compute_distance <- function(coords1, coords2) {
  if (is.null(coords2)) {
    return(as.matrix(dist(coords1)))
  }
  
  c1 <- as.matrix(coords1)
  c2 <- as.matrix(coords2)
  
  # Vectorized calculation: d^2(x, y) = sum(x^2) + sum(y^2) - 2 * x %*% t(y)
  # This is much faster than looping in R
  
  sq_c1 <- rowSums(c1^2)
  sq_c2 <- rowSums(c2^2)
  
  # Outer sum: n x m matrix where (i,j) is sq_c1[i] + sq_c2[j]
  # We can do this efficiently
  cross_term <- tcrossprod(c1, c2)
  
  # Use outer to add squared sums
  # A_ij = sq_c1[i] + sq_c2[j]
  # Note: outer can be memory intensive for massive n*m.
  # For n=100k, m=50, it is 5M doubles = 40MB, totally fine.
  
  dists_sq <- outer(sq_c1, sq_c2, "+") - 2 * cross_term
  
  # Numerical stability: clamp negative small values to 0
  dists_sq[dists_sq < 0] <- 0
  
  return(sqrt(dists_sq))
}

#' Kernel functions
#' 
#' @param d Vector or matrix of distances
#' @param h Bandwidth
#' @return Vector or matrix of weights
kernel_bisquare <- function(d, h) {
  if (h <= 0) return(as.numeric(d == 0)) # degenerate case
  u <- d / h
  w <- (1 - u^2)^2
  w[abs(u) >= 1] <- 0
  return(w)
}

kernel_gaussian <- function(d, h) {
  if (h <= 0) return(as.numeric(d == 0))
  u <- d / h
  return(exp(-0.5 * u^2))
}

kernel_exponential <- function(d, h) {
  if (h <= 0) return(as.numeric(d == 0))
  u <- d / h
  return(exp(-abs(u)))
}

#' Get Kernel Weights
#' 
#' @param dists Distance matrix
#' @param bandwidth Bandwidth value
#' @param kernel_type Kernel name
#' @return Weight matrix
get_weights <- function(dists, bandwidth, kernel_type = "bisquare") {
  if (kernel_type == "bisquare") {
    return(kernel_bisquare(dists, bandwidth))
  } else if (kernel_type == "gaussian") {
    return(kernel_gaussian(dists, bandwidth))
  } else if (kernel_type == "exponential") {
    return(kernel_exponential(dists, bandwidth))
  } else {
    stop("Unknown kernel type")
  }
}

#' Select Anchors using K-means
#' 
#' @param coords Coordinate matrix
#' @param n_anchors Number of anchors
#' @param seed Random seed
#' @return Matrix of anchor coordinates
select_anchors <- function(coords, n_anchors, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # If n_anchors >= n_obs, use all points
  if (n_anchors >= nrow(coords)) {
    return(as.matrix(coords))
  }
  
  km <- kmeans(coords, centers = n_anchors, nstart = 5, iter.max = 50)
  return(km$centers)
}

#' Calculate Effective Sample Size
#' 
#' @param weights Numeric vector of weights
#' @return Numeric ESS value
compute_ess <- function(weights) {
  sum_w <- sum(weights)
  if (sum_w == 0) return(0)
  return((sum_w^2) / sum(weights^2))
}

#' Create Spatial Blocks
#' 
#' Partitions observations into spatial blocks using K-means clustering.
#' Used for spatial cross-validation and bootstrap.
#' 
#' @param coords Coordinate matrix (n x 2).
#' @param n_blocks Integer, number of spatial blocks.
#' @param seed Random seed for reproducibility.
#' @return Integer vector of length n, indicating block membership (1 to n_blocks).
#' @export
create_spatial_blocks <- function(coords, n_blocks = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  coords <- as.matrix(coords)
  
  if (n_blocks >= nrow(coords)) {
    return(seq_len(nrow(coords)))
  }
  
  km <- kmeans(coords, centers = n_blocks, nstart = 5, iter.max = 50)
  return(km$cluster)
}
