
#' Compute the aggregated empirical cumulative distribution function for all
#' data inside a circle with a given radius. Input variables are:
#' data: An (n x d)-dimensional matrix of observations.
#' coords: A (d x 2)-dimensional matrix with the coordinates of the data.
#' center: A 2-dimensional vector containing the coordinates of the center we use
#'   for computing the aggregated ECDF.
#' radius: The radius of the circle used for computing the aggregated ECDF.
#' @export
aggregated_ecdf = function(data, coords, center, radius, zero_threshold = 0) {
  stopifnot(ncol(data) == nrow(coords))
  window_index = extract_thinned_out_circles(
    coords = coords,
    center = center,
    n = 1,
    r = radius,
    index_only = TRUE)
  x = as.numeric(data[, window_index])
  marginal_distribution(x, zero_threshold)
}


#' Given a matrix of coordinates on a regular grid, extract a subset of the coordinates
#' that consists of circles centred around some center, with varying degrees of
#' densities. The input variables are:
#' coords: A matrix with two columns and multiple rows, containing (x, y) coordinates
#'   in a regular grid and with a distance-preserving projection.
#' center: One (x, y) tuple of coordinates, explaining the center of the circles we
#'   wish to extract.
#' n: A vector of densities for the different circles we wish to extract. If the ith
#'   value of n is k, then we only keep every kth of the original x-coordinates, and every
#'   kth of the original y-coordinates.
#' r: A vector of radii for the different circles we wish to extract. If e.g.
#'   n = (1, 3) and r = (2, 5), then we will extract every coordinate that is closer
#'   to the center than a radius of 2. Then, we will only keep every 3rd x- and y-coordinates,
#'   for all coordinates that have a distance between 2 and 5 from the center.
#'   Finally, all coords that are further away from the center than 5 will be dropped.
#' index_only: A bool describing if the function should return the extracted
#'   coordinates, or only their indices inside the coords object.
#' @export
extract_thinned_out_circles = function(coords,
                                       center,
                                       n,
                                       r = Inf,
                                       index_only = FALSE) {
  stopifnot(ncol(coords) == 2)
  stopifnot(length(n) == length(r))
  stopifnot(is.matrix(coords))

  # Locate all unique x-values and y-values
  x_vals = sort(unique(coords[, 1]))
  y_vals = sort(unique(coords[, 2]))

  # Compute distances and find the coord that is closest to the center
  dist_to_center = as.numeric(dist_euclid(center, coords))
  closest_to_center = coords[which.min(dist_to_center), ]
  x_center_index = which(x_vals == closest_to_center[1])
  y_center_index = which(y_vals == closest_to_center[2])

  # This will be the list of all indices we wish to extract from coords
  chosen_indices = NULL

  # Locate all the indices we wish to extract
  r = c(0, r)
  for (i in seq_along(n)) {
    # Indicies for all coords that have the correct distance from center
    dist_index = which(r[i] <= dist_to_center & dist_to_center <= r[i + 1])

    # Indices for all coords on a regular grid of resolution n[i] x n[i]
    thinned_x_index = c(seq(x_center_index, 1, by = -n[i]),
                        seq(x_center_index, length(x_vals), by = n[i]))
    thinned_y_index = c(seq(y_center_index, 1, by = -n[i]),
                        seq(y_center_index, length(y_vals), by = n[i]))
    thinned_x_vals = x_vals[thinned_x_index]
    thinned_y_vals = y_vals[thinned_y_index]
    thinned_index = which(coords[, 1] %in% thinned_x_vals & coords[, 2] %in% thinned_y_vals)

    # The indices of interest are the intersection of the two index sets above
    chosen_indices = c(chosen_indices, intersect(dist_index, thinned_index))
  }
  # Some times, duplicates may appear. Remove these
  chosen_indices = unique(chosen_indices)

  if (index_only) {
    res = chosen_indices
  } else {
    res = coords[chosen_indices, ]
  }
  res
}

#' Search through data after threshold exceedances at a set of conditioning sites.
#' Extract relevant information and observations for each of the threshold exceedances.
#'
#' The input variables are:
#' data: An (n x d)-dimensional matrix of observations.
#' coords: A (d x 2)-dimensional matrix with the coordinates of the data.
#' s0_index: a vector of indices describing which of the d coordinates are chosen
#'   as conditioning sites.
#' threshold: The threshold t for the conditional extremes model.
#' n, r: Used for only extracting a subset of the observations in a circle around
#'   the conditioning site. See the docs for extract_thinned_out_circles() for more info.
#'
#' The output of the function is a list where each element is a list or vector
#' of length equal, containing these variables:
#' s0_index: This is a subset of all the s0_indices given as input such that each
#'   conditioning site has at least one threshold exceedance.
#' s0: A list where element nr. i is a 2-dimensional vector with the coordinates of
#'   conditioning site nr. i.
#' obs_index: A list where element nr. i is a vector of indices for all the
#'   locations in coords we use for performing inference using conditioning site nr. i.
#' dist_to_s0: A list where element nr. i is a vector with all the distances between
#'   conditioning site nr. i and the locations used for inference with conditioning
#'   site nr. i.
#' y0: A list where element nr. i is a vector of all threshold exceedances for
#'   conditioning site nr. i.
#' y: A list where element nr. i is a matrix consisting of one column per threshold
#'   exceedance in y0[[i]]. Each column contains observations from the locations
#'   described in obs_index[[i]].
#' time_index: A list where element nr. i contains the times of all threshold
#'   exceedances from conditioning site nr. i.
#' n: A vector where element nr. i is the number of threshold exceedances at
#'   conditioning site nr. i.
#' @export
extract_extreme_fields = function(data,
                                  coords,
                                  s0_index,
                                  threshold,
                                  n,
                                  r,
                                  verbose = FALSE,
                                  n_cores = 1) {
  # s0_indices must be between 1 and nrow(coords)
  stopifnot(1 <= min(s0_index) && max(s0_index) <= nrow(coords))
  n_s0 = length(s0_index)

  # Look through the conditioning sites
  res = parallel::mclapply(
    X = seq_len(n_s0),
    mc.cores = n_cores,
    FUN = function(i) {
      res = list(s0_index = s0_index[i])

      # Extract all variables at s0 nr. i, and locate threshold exceedances
      y0 = data[, s0_index[i]]
      res$time_index = which(y0 > threshold)
      # If there are no threshold exceedances, skip to next conditioning site
      if (length(res$time_index) == 0) return(NULL)

      # Extract all threshold exceedances
      res$y0 = y0[res$time_index]

      # Locate the thinned out observations locations
      res$obs_index = extract_thinned_out_circles(
        coords = coords,
        center = coords[s0_index[i], , drop = FALSE],
        n = n,
        r = r,
        index_only = TRUE)
      res$obs_index = res$obs_index[res$obs_index != s0_index[i]]

      # Extract all observations of interest for the threshold exceedances
      res$y = t(data[res$time_index, res$obs_index, drop = FALSE])

      # Fill in information about s0 and the dist_to_s0
      res$s0 = coords[s0_index[i], , drop = FALSE]
      res$dist_to_s0 = as.numeric(dist_euclid(
        x = coords[s0_index[i], , drop = FALSE],
        y = coords[res$obs_index, , drop = FALSE]))

      if (verbose) message(i, " / ", n_s0)

      res
    })

  # Only keep information about the conditioning sites with at least one threshold exceedance
  good_index = which(sapply(res, length) > 0)
  res = purrr::transpose(res[good_index])

  res$s0_index = unlist(res$s0_index)

  # Find the number of threshold exceedances for each conditioning site
  res$n = sapply(res$y0, length)

  res
}


#' Compute the empirical estimators of χ_p(d) for a set of
#' threshold values and distances.
#' Function input:
#' - data is a list of data on the same format as that created by
#'   the function extract_extreme_fields().
#' - thresholds is a vector of threshold values used for estimating
#'   the extremal correlation coefficient.
#' - unique_dist is an integer vector of distances d for which we
#'   should estimate χ_p(d)
#' @export
empirical_chi = function(data, thresholds, unique_dist) {
  n_big = matrix(0, length(unique_dist), length(thresholds))
  n_all = matrix(0, length(unique_dist), length(thresholds))
  pb = progress_bar(length(data$y))
  for (j in seq_along(data$y)) {
    dist = round(data$dist_to_s0[[j]])
    for (d in seq_along(unique_dist)) {
      d_index = which(unique_dist[d] == dist)
      if (!any(d_index)) next
      for (i in seq_along(thresholds)) {
        y0_index = which(data$y0[[j]] > thresholds[i])
        if (!any(y0_index)) next
        n_big[d, i] = n_big[d, i] +
          sum(data$y[[j]][d_index, y0_index] > thresholds[i], na.rm = TRUE)
        n_all[d, i] = n_all[d, i] +
          sum(!is.na(data$y[[j]][d_index, y0_index]))
      }
    }
    pb$tick()
  }
  pb$terminate()
  chi = n_big / n_all
  attr(chi, "thresholds") = thresholds
  attr(chi, "unique_dist") = unique_dist
  chi
}

#' Plot the estimates of χ_p(d), created using the empirical_chi() function
#' @export
plot_chi = function(chi) {
  t = attr(chi, "thresholds")
  unique_dist = attr(chi, "unique_dist")
  data.frame(
    chi = as.numeric(rbind(1, chi)),
    dist = rep(c(0, unique_dist), length(t)),
    t = rep(t, each = length(unique_dist) + 1)) |>
    ggplot() +
    geom_line(aes(x = dist, y = chi, group = t, col = t))
}

#' Compute the empirical estimators of the first two conditional moments of a dataset
#' for a set of threshold exceedances and distances.
#' Function input:
#' - data is a list of data on the same format as that created by
#'   the function extract_extreme_fields().
#' - unique_dist is an integer vector of distances d for which we
#'   should estimate the conditional moments μ(d; y_0) and ζ(d; y_0)
#' - unique_y0 is a vector of threshold exceedances y0 for which we
#'   should estimate the conditional moments μ(d; y_0) and ζ(d; y_0)
#' @export
empirical_moments = function(data, unique_dist, unique_y0) {
  y_sum = matrix(0, length(unique_dist), length(unique_y0))
  y_sq_sum = matrix(0, length(unique_dist), length(unique_y0))
  y_n = matrix(0, length(unique_dist), length(unique_y0))
  pb = progress_bar(length(data$y))
  for (j in seq_along(data$y)) {
    dist = round(data$dist_to_s0[[j]])
    y0 = round(data$y0[[j]], 1)
    for (d in seq_along(unique_dist)) {
      d_index = which(unique_dist[d] == dist)
      if (!any(d_index)) next
      for (i in seq_along(unique_y0)) {
        y0_index = which(y0 == unique_y0[i])
        if (!any(y0_index)) next
        y_sum[d, i] = y_sum[d, i] + sum(data$y[[j]][d_index, y0_index], na.rm = TRUE)
        y_sq_sum[d, i] = y_sq_sum[d, i] + sum(data$y[[j]][d_index, y0_index]^2, na.rm = TRUE)
        y_n[d, i] = y_n[d, i] + sum(!is.na(data$y[[j]][d_index, y0_index]))
      }
    }
    pb$tick()
  }
  pb$terminate()
  y_mean = y_sum / y_n
  y_sq_mean = y_sq_sum / y_n
  y_sd = sqrt(y_sq_mean - y_mean^2)
  res = list(mean = y_mean, sd = y_sd)
  attr(res, "unique_dist") = unique_dist
  attr(res, "unique_y0") = unique_y0
  res
}

# Plot empirical moments created using the empirical_moments() function
#' @export
plot_moments = function(x) {
  unique_dist = attr(x, "unique_dist")
  unique_y0 = attr(x, "unique_y0")
  sd_plot = data.frame(
    sd = as.numeric(rbind(0, x$sd)),
    dist = rep(c(0, unique_dist), length(unique_y0)),
    y0 = rep(unique_y0, each = length(unique_dist) + 1)) |>
    ggplot() +
    geom_line(aes(x = dist, y = sd, group = y0, col = y0))
  mean_plot = data.frame(
    mean = as.numeric(rbind(unique_y0, x$mean)),
    dist = rep(c(0, unique_dist), length(unique_y0)),
    y0 = rep(unique_y0, each = length(unique_dist) + 1)) |>
    ggplot() +
    geom_line(aes(x = dist, y = mean, group = y0, col = y0))
  list(mean = mean_plot, sd = sd_plot)
}

#' Given a matrix of coordinates, return a sub-grid of the coordinates
#' with resolution (delta_s0 x delta_s0)
#' @export
get_s0_index = function(coords, delta_s0 = 1) {
  stopifnot(all(as.integer(delta_s0) == delta_s0))
  x_coords = coords[, 1] |> unique() |> sort()
  y_coords = coords[, 2] |> unique() |> sort()
  s0_locs = expand.grid(
    x = x_coords[seq(delta_s0, length(x_coords), by = delta_s0)],
    y = y_coords[seq(delta_s0, length(y_coords), by = delta_s0)]) |>
    as.matrix()
  s0_index = lapply(
    X = 1:nrow(s0_locs),
    FUN = function(i) which(coords[, 1] == s0_locs[i, 1] & coords[, 2] == s0_locs[i, 2]))
  s0_index = s0_index[sapply(s0_index, length) > 0] |>
    unlist() |>
    unname()
  s0_index
}
