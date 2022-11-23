#' Estimate the extremal correlation where we group together all observations
#' that have similar distances to s0
#' y0: n-dimensional vector of y(s0)
#' y: (nxd)-dimensional matrix of observations with Laplace marginals
#' dist_to_s0: d-dimensional vector of lenghts
#' dist_groups: vector of distances such that group nr. i contains all locations with
#'   dist_groups[i] <= dist_to_s0 < dist_groups[i + 1]
#' p: vector of probabilities for computing χ
#' @export
extremal_correlation_in_dist_groups = function(y, y0, dist_to_s0, dist_groups, p) {
  # Place all the distances into their respective dist_groups
  dist_groups = sort(dist_groups)
  if (is.finite(tail(dist_groups, 1))) dist_groups = c(dist_groups, Inf)
  dist_group_index = vector("list", length(dist_groups) - 1)
  for (j in seq_along(dist_groups[-1])) {
    dist_group_index[[j]] = which(dist_to_s0 >= dist_groups[j] & dist_to_s0 < dist_groups[j + 1])
  }

  # Compute the total number of large observations and the total number of observations
  # in each of the dist_groups, for each of the probabilities
  res = matrix(nrow = length(dist_groups) - 1, ncol = length(p))
  for (j in seq_along(p)) {
    big_index = which(y0 > qlaplace(p[j]))
    for (k in seq_along(dist_groups[-1])) {
      res[k, j] = mean(y[big_index, dist_group_index[[k]]] > qlaplace(p[j]), na.rm = TRUE)
    }
  }

  res
}

#' Calculate chi and chibar at the quantiles u for two random variables
#' This one is heavily based on the implementation in
#' https://github.com/harrysouthworth/texmex/blob/master/R/chi.R
#' @export
extremal_correlation = function(p, x, y, truncate = TRUE) {
  n = length(x)
  if (n != length(y)) stop("The data have different length")
  x_quantiles = rank(x, ties.method = "first") / (n + 1)
  y_quantiles = rank(y, ties.method = "first") / (n + 1)
  quantiles = cbind(x_quantiles, y_quantiles)

  rowmax = apply(quantiles, 1, max)
  rowmin = apply(quantiles, 1, min)
  if (min(p) < min(rowmax)) stop("lower quantile limit is too low")
  if (max(p) > max(rowmax)) stop("upper quantile limit is too high")

  # Find probability of both components being larger than p (C) and of
  # both components being less than p (C_bar)
  C = sapply(p, function(p) mean(rowmax < p))
  C_bar = sapply(p, function(p) mean(rowmin > p))

  # Get chi and chibar
  chi = 2 - log(C) / log(p) # 3.2 of Coles, Heffernan, Tawn
  chi_bar = 2 * log(1 - p) / log(C_bar) - 1 # Page 348 of Coles, Heffernan, Tawn

  # Estimate standard deviation using delta method
  σ_chi = sqrt((1 / log(p)^2) / C * (1 - C) / n)
  σ_chi_bar = sqrt((((4 * log(1 - p)^2) / (log(C_bar)^4 * C_bar^2)) *
                      C_bar * (1 - C_bar)) / n)

  # Truncate values of chi or chibar that are outside theoretical bounds
  if (truncate) {
    chivals = truncate_chi(chi, chi_bar, p)
    chi = chivals$chi
    chi_bar = chivals$chi_bar
  }

  data.frame(
    val = c(chi, chi_bar),
    sd = c(σ_chi, σ_chi_bar),
    par = rep(c("chi", "chibar"), each = length(p)),
    p = rep(p, 2))
}

truncate_chi = function(chi, chi_bar, p) {
  chi_lb = 2 - log(pmax(2 * p - 1, 0)) / log(p)
  chi_bar_lb = 2 * log(1 - p) / log(1 - 2 * p + pmax(2 * p - 1, 0)) - 1
  chi[chi > 1] = 1
  chi_bar[chi_bar > 1] = 1
  chi = sapply(seq_along(chi), function(i) pmax(chi[i], chi_lb[i]))
  chi_bar = sapply(seq_along(chi_bar), function(i) pmax(chi_bar[i], chi_bar_lb[i]))
  list(chi = chi, chi_bar = chi_bar)
}

