#' @importFrom sf st_coordinates st_centroid st_geometry_type st_transform st_crs st_drop_geometry
#' @importFrom stats complete.cases
NULL

# =============================================================================
# Internal helper: stratified sampling index
# Returns indices into the coordinate matrix, not the sampled data itself.
# Fixes O(n^2) loop by pre-allocating result vector.
# =============================================================================
.stratified_sample_idx <- function(xy, sample_n, sample_type, n_strata = NULL) {
  n_valid <- nrow(xy)
  
  if (sample_type == "random") {
    sample_idx <- sample(n_valid, sample_n)
    message(sprintf("Sampled %d cells using random sampling.", sample_n))
    return(sample_idx)
  }
  
  # Stratified sampling
  if (is.null(n_strata)) {
    n_strata <- min(100, max(9, floor(sqrt(sample_n))))
  }
  
  x_range <- range(xy[, 1])
  y_range <- range(xy[, 2])
  
  n_x <- ceiling(sqrt(n_strata))
  n_y <- ceiling(n_strata / n_x)
  
  x_breaks <- seq(x_range[1], x_range[2], length.out = n_x + 1)
  y_breaks <- seq(y_range[1], y_range[2], length.out = n_y + 1)
  
  x_bin <- findInterval(xy[, 1], x_breaks, all.inside = TRUE)
  y_bin <- findInterval(xy[, 2], y_breaks, all.inside = TRUE)
  stratum <- (y_bin - 1) * n_x + x_bin
  
  strata_counts <- table(stratum)
  strata_ids <- as.integer(names(strata_counts))
  strata_sizes <- as.integer(strata_counts)
  strata_props <- strata_sizes / sum(strata_sizes)
  
  # Allocate samples proportionally
  samples_per_stratum <- pmax(1L, round(strata_props * sample_n))
  
  # Adjust to match exact sample_n
  total_allocated <- sum(samples_per_stratum)
  if (total_allocated > sample_n) {
    excess <- total_allocated - sample_n
    order_idx <- order(samples_per_stratum, decreasing = TRUE)
    for (i in order_idx) {
      if (excess <= 0L) break
      reduce <- min(samples_per_stratum[i] - 1L, excess)
      samples_per_stratum[i] <- samples_per_stratum[i] - reduce
      excess <- excess - reduce
    }
  } else if (total_allocated < sample_n) {
    deficit <- sample_n - total_allocated
    order_idx <- order(strata_sizes, decreasing = TRUE)
    for (i in order_idx) {
      if (deficit <= 0L) break
      add <- min(strata_sizes[i] - samples_per_stratum[i], deficit)
      samples_per_stratum[i] <- samples_per_stratum[i] + add
      deficit <- deficit - add
    }
  }
  
  # PRE-ALLOCATE result vector (fixes O(n^2) append)
  total_samples <- sum(samples_per_stratum)
  sample_idx <- integer(total_samples)
  write_pos <- 1L
  
  for (i in seq_along(strata_ids)) {
    sid <- strata_ids[i]
    stratum_points <- which(stratum == sid)
    n_sample_stratum <- min(samples_per_stratum[i], length(stratum_points))
    
    if (n_sample_stratum > 0L) {
      sampled <- sample(stratum_points, n_sample_stratum)
      sample_idx[write_pos:(write_pos + n_sample_stratum - 1L)] <- sampled
      write_pos <- write_pos + n_sample_stratum
    }
  }
  
  # Trim if we allocated more than filled (shouldn't happen, but safe)
  sample_idx <- sample_idx[1:(write_pos - 1L)]
  
  n_strata_actual <- length(strata_ids)
  message(sprintf("Sampled %d cells using stratified sampling (%d strata, %dx%d grid).", 
                  length(sample_idx), n_strata_actual, n_x, n_y))
  
  return(sample_idx)
}

#' Prepare Data from sf Object
#'
#' Extracts Y, W, X, and coordinates from an sf object for use with gw_causal_forest.
#'
#' @param sf_obj An sf object (points or polygons).
#' @param y_col Character, name of the outcome variable column.
#' @param w_col Character, name of the treatment variable column.
#' @param x_cols Character vector, names of covariate columns.
#' @param coord_type Character, how to extract coordinates for polygons:
#'   "centroid" (default) uses polygon centroids,
#'   "point" assumes input is already points.
#' @param target_crs Optional CRS to transform coordinates to (e.g., EPSG code).
#'   If NULL, uses original CRS. Recommended: use projected CRS for local studies.
#' @param drop_na Logical, whether to remove rows with NA values (default TRUE).
#'
#' @return A list with components:
#'   \item{Y}{Numeric vector of outcomes}
#'   \item{W}{Numeric vector of treatment}
#'   \item{X}{Matrix of covariates}
#'   \item{coords}{Matrix of coordinates (n x 2)}
#'   \item{crs}{CRS of the coordinates}
#'   \item{n}{Number of observations}
#'   \item{removed_na}{Number of rows removed due to NA}
#'
#' @examples
#' \dontrun{
#' library(sf)
#' 
#' # From shapefile
#' shp <- st_read("districts.shp")
#' 
#' data <- prep_from_sf(
#'   sf_obj = shp,
#'   y_col = "gdp_growth",
#'   w_col = "policy_treated",
#'   x_cols = c("population", "education", "industry_share")
#' )
#' 
#' # Fit model
#' fit <- gw_causal_forest(
#'   Y = data$Y,
#'   W = data$W,
#'   X = data$X,
#'   coords = data$coords,
#'   treatment_type = "binary"
#' )
#' }
#'
#' @export
prep_from_sf <- function(sf_obj,
                         y_col,
                         w_col,
                         x_cols,
                         coord_type = c("centroid", "point"),
                         target_crs = NULL,
                         drop_na = TRUE) {
  
  coord_type <- match.arg(coord_type)
  
 
  # Validate inputs
  if (!inherits(sf_obj, "sf")) {
    stop("sf_obj must be an sf object. Use sf::st_read() to load shapefiles.")
  }
  
  all_cols <- c(y_col, w_col, x_cols)
  missing_cols <- setdiff(all_cols, names(sf_obj))
  if (length(missing_cols) > 0) {
    stop(sprintf("Columns not found in sf_obj: %s", paste(missing_cols, collapse = ", ")))
  }
  
  # Transform CRS if requested
  if (!is.null(target_crs)) {
    sf_obj <- sf::st_transform(sf_obj, target_crs)
  }
  
  # Extract coordinates based on geometry type
  geom_type <- unique(as.character(sf::st_geometry_type(sf_obj)))
  
  if (coord_type == "centroid" || any(geom_type %in% c("POLYGON", "MULTIPOLYGON"))) {
    # Use centroids for polygons
    coords <- sf::st_coordinates(sf::st_centroid(sf_obj, of_largest_polygon = TRUE))
  } else {
    # Direct coordinates for points
    coords <- sf::st_coordinates(sf_obj)
  }
  
  # Ensure 2D coordinates
  if (ncol(coords) > 2) {
    coords <- coords[, 1:2]
  }
  colnames(coords) <- c("X", "Y")
  
  # Extract attribute data
  df <- sf::st_drop_geometry(sf_obj)
  
  Y <- as.numeric(df[[y_col]])
  W <- as.numeric(df[[w_col]])
  X <- as.matrix(df[, x_cols, drop = FALSE])
  
  # Handle NA values
  n_original <- nrow(df)
  complete_idx <- complete.cases(Y, W, X, coords)
  n_complete <- sum(complete_idx)
  n_removed <- n_original - n_complete
  
 if (drop_na && n_removed > 0) {
    message(sprintf("Removed %d rows with NA values (%d remaining).", n_removed, n_complete))
    Y <- Y[complete_idx]
    W <- W[complete_idx]
    X <- X[complete_idx, , drop = FALSE]
    coords <- coords[complete_idx, , drop = FALSE]
  } else if (!drop_na && n_removed > 0) {
    warning(sprintf("%d rows contain NA values. Set drop_na = TRUE to remove them.", n_removed))
  }
  
  # Get CRS info
  crs_info <- sf::st_crs(sf_obj)
  
  result <- list(
    Y = Y,
    W = W,
    X = X,
    coords = coords,
    crs = crs_info,
    n = length(Y),
    removed_na = n_removed,
    geometry_type = geom_type[1]
  )
  
  class(result) <- "gwcf_data"
  return(result)
}


#' Prepare Data from Raster
#'
#' Extracts Y, W, X, and coordinates from raster data for use with gw_causal_forest.
#' Requires the terra package.
#'
#' @param y_raster A SpatRaster (terra) or file path for the outcome variable.
#' @param w_raster A SpatRaster or file path for the treatment variable.
#' @param x_rasters A SpatRaster with multiple layers, or a character vector of file paths
#'   for covariate rasters.
#' @param sample_n Optional integer. If specified, samples this many cells
#'   instead of using all valid cells. Useful for large rasters.
#' @param sample_type Character, sampling method:
#'   \itemize{
#'     \item "stratified" (default): Spatially stratified random sampling. Divides
#'       the study area into grid cells and samples proportionally from each.
#'       Ensures spatial representativeness.
#'     \item "random": Simple random sampling. May result in spatial clustering.
#'   }
#' @param n_strata Integer, number of spatial strata (grid cells) for stratified
#'   sampling. Default is sqrt(sample_n), capped at 100.
#' @param mask_raster Optional SpatRaster to mask the study area.
#' @param seed Random seed for reproducible sampling.
#'
#' @return A list with components:
#'   \item{Y}{Numeric vector of outcomes}
#'   \item{W}{Numeric vector of treatment}
#'   \item{X}{Matrix of covariates}
#'   \item{coords}{Matrix of coordinates (n x 2)}
#'   \item{crs}{CRS of the rasters}
#'   \item{n}{Number of observations}
#'   \item{resolution}{Raster resolution}
#'
#' @examples
#' \dontrun{
#' # Using terra
#' library(terra)
#' 
#' # Load rasters
#' y_rast <- rast("outcome.tif")
#' w_rast <- rast("treatment.tif")
#' x_rast <- rast(c("covar1.tif", "covar2.tif", "covar3.tif"))
#' 
#' # Prepare data (sample 10000 cells if raster is large)
#' data <- prep_from_raster(
#'   y_raster = y_rast,
#'   w_raster = w_rast,
#'   x_rasters = x_rast,
#'   sample_n = 10000,
#'   seed = 123
#' )
#' 
#' # Fit model
#' fit <- gw_causal_forest(
#'   Y = data$Y,
#'   W = data$W,
#'   X = data$X,
#'   coords = data$coords,
#'   treatment_type = "continuous"
#' )
#' }
#'
#' @export
prep_from_raster <- function(y_raster,
                             w_raster,
                             x_rasters,
                             sample_n = NULL,
                             sample_type = c("stratified", "random"),
                             n_strata = NULL,
                             mask_raster = NULL,
                             seed = NULL) {
  
  sample_type <- match.arg(sample_type)
  
  # Check if terra is available
 if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for raster processing. Install with: install.packages('terra')")
  }
  
  # Load rasters if file paths provided
  load_raster <- function(r) {
    if (is.character(r)) {
      terra::rast(r)
    } else if (inherits(r, "SpatRaster")) {
      r
    } else {
      stop("Raster inputs must be SpatRaster objects or file paths.")
    }
  }
  
  y_rast <- load_raster(y_raster)
  w_rast <- load_raster(w_raster)
  x_rast <- load_raster(x_rasters)
  
  # Stack all rasters
  all_rast <- c(y_rast, w_rast, x_rast)
  n_layers <- terra::nlyr(all_rast)
  
  # Apply mask if provided
  if (!is.null(mask_raster)) {
    mask_rast <- load_raster(mask_raster)
    all_rast <- terra::mask(all_rast, mask_rast)
  }
  
  total_cells <- terra::ncell(all_rast)
  
  # ==========================================================================
  # MEMORY-EFFICIENT SAMPLING
  # For large rasters (>1M cells), avoid loading all data before sampling.
  # Strategy: use terra::spatSample for random, or chunked processing for stratified.
  # ==========================================================================
  
  if (!is.null(sample_n) && total_cells > 1e6) {
    # Large raster path: memory-efficient
    if (!is.null(seed)) set.seed(seed)
    
    message(sprintf("Large raster (%d cells). Using memory-efficient sampling...", total_cells))
    
    if (sample_type == "random") {
      # Use terra::spatSample - samples directly without loading all data
      # Oversample by 50% to account for NA removal
      oversample <- ceiling(sample_n * 1.5)
      
      sampled <- terra::spatSample(all_rast, size = oversample, method = "random",
                                   na.rm = TRUE, xy = TRUE, as.df = TRUE)
      
      # Take exactly sample_n rows (or all if less available)
      if (nrow(sampled) > sample_n) {
        sampled <- sampled[sample(nrow(sampled), sample_n), ]
      }
      
      xy <- as.matrix(sampled[, c("x", "y")])
      vals <- as.matrix(sampled[, 3:ncol(sampled)])
      
      message(sprintf("Sampled %d cells using random sampling.", nrow(vals)))
      
    } else if (sample_type == "stratified") {
      # Stratified sampling for large rasters
      # Step 1: Get raster extent and create spatial grid
      ext <- terra::ext(all_rast)
      
      if (is.null(n_strata)) {
        n_strata <- min(100, max(9, floor(sqrt(sample_n))))
      }
      
      n_x <- ceiling(sqrt(n_strata))
      n_y <- ceiling(n_strata / n_x)
      
      # Create strata boundaries
      x_breaks <- seq(ext$xmin, ext$xmax, length.out = n_x + 1)
      y_breaks <- seq(ext$ymin, ext$ymax, length.out = n_y + 1)
      
      # Sample from each stratum using terra::spatSample with ext argument
      samples_per_stratum <- ceiling(sample_n / (n_x * n_y))
      
      sample_list <- vector("list", n_x * n_y)
      idx <- 1
      
      for (i in 1:n_x) {
        for (j in 1:n_y) {
          stratum_ext <- terra::ext(x_breaks[i], x_breaks[i+1], 
                                    y_breaks[j], y_breaks[j+1])
          
          stratum_rast <- terra::crop(all_rast, stratum_ext)
          
          # Sample from this stratum
          stratum_sample <- tryCatch({
            terra::spatSample(stratum_rast, size = samples_per_stratum,
                              method = "random", na.rm = TRUE, xy = TRUE, as.df = TRUE)
          }, error = function(e) NULL)
          
          if (!is.null(stratum_sample) && nrow(stratum_sample) > 0) {
            sample_list[[idx]] <- stratum_sample
          }
          idx <- idx + 1
        }
      }
      
      # Combine all stratum samples
      sampled <- do.call(rbind, sample_list[!sapply(sample_list, is.null)])
      
      # Trim to exact sample_n
      if (nrow(sampled) > sample_n) {
        sampled <- sampled[sample(nrow(sampled), sample_n), ]
      }
      
      xy <- as.matrix(sampled[, c("x", "y")])
      vals <- as.matrix(sampled[, 3:ncol(sampled)])
      
      n_strata_actual <- sum(!sapply(sample_list, is.null))
      message(sprintf("Sampled %d cells using stratified sampling (%d strata, %dx%d grid).", 
                      nrow(vals), n_strata_actual, n_x, n_y))
    }
    
  } else {
    # Small raster or no sampling: original path (load all, then optionally sample)
    vals <- terra::values(all_rast)
    xy <- terra::xyFromCell(all_rast, 1:total_cells)
    
    # Find complete cases
    complete_idx <- complete.cases(vals)
    n_valid <- sum(complete_idx)
    
    if (n_valid == 0) {
      stop("No valid (non-NA) cells found after combining all rasters.")
    }
    
    message(sprintf("Found %d valid cells out of %d total.", n_valid, nrow(vals)))
    
    # Subset to valid cells
    vals <- vals[complete_idx, , drop = FALSE]
    xy <- xy[complete_idx, , drop = FALSE]
    
    # Sample if requested
    if (!is.null(sample_n) && sample_n < n_valid) {
      if (!is.null(seed)) set.seed(seed)
      
      sample_idx <- .stratified_sample_idx(xy, sample_n, sample_type, n_strata)
      
      vals <- vals[sample_idx, , drop = FALSE]
      xy <- xy[sample_idx, , drop = FALSE]
    }
  }
  
  # Split into Y, W, X
  Y <- vals[, 1]
  W <- vals[, 2]
  X <- vals[, 3:n_layers, drop = FALSE]
  
  # Set column names for X
  x_names <- names(x_rast)
  if (is.null(x_names) || any(x_names == "")) {
    x_names <- paste0("X", seq_len(ncol(X)))
  }
  colnames(X) <- x_names
  colnames(xy) <- c("X", "Y")
  
  # Get metadata
  crs_info <- terra::crs(all_rast)
  res_info <- terra::res(all_rast)
  
  result <- list(
    Y = Y,
    W = W,
    X = X,
    coords = xy,
    crs = crs_info,
    n = length(Y),
    resolution = res_info,
    sampled = !is.null(sample_n)
  )
  
  class(result) <- "gwcf_data"
  return(result)
}


#' Print method for gwcf_data
#'
#' @param x A gwcf_data object.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns x.
#' @export
print.gwcf_data <- function(x, ...) {
  cat("GWCF Prepared Data\n")
  cat("==================\n\n")
  cat(sprintf("Observations: %d\n", x$n))
  cat(sprintf("Covariates: %d (%s)\n", ncol(x$X), paste(colnames(x$X), collapse = ", ")))
  cat(sprintf("Treatment range: [%.3f, %.3f]\n", min(x$W), max(x$W)))
  cat(sprintf("Outcome range: [%.3f, %.3f]\n", min(x$Y), max(x$Y)))
  cat(sprintf("Coordinate range X: [%.2f, %.2f]\n", min(x$coords[,1]), max(x$coords[,1])))
  cat(sprintf("Coordinate range Y: [%.2f, %.2f]\n", min(x$coords[,2]), max(x$coords[,2])))
  
  if (!is.null(x$geometry_type)) {
    cat(sprintf("Source geometry: %s\n", x$geometry_type))
  }
  if (!is.null(x$resolution)) {
    cat(sprintf("Raster resolution: %.2f x %.2f\n", x$resolution[1], x$resolution[2]))
  }
  if (!is.null(x$removed_na) && x$removed_na > 0) {
    cat(sprintf("Removed NA rows: %d\n", x$removed_na))
  }
  
  invisible(x)
}

