
#' Compute the empirical marginal distribution of the vector x and return
#' a function F(x) that is equal to the marginal distribution when x > zero_val
#' and returns NA when x <= zero_val.
#' @export
marginal_distribution = function(x, zero_threshold = 0) {
  positive_index = which(x > zero_threshold)
  stopifnot(any(positive_index))
  F = ecdf(x[positive_index])
  rm(x, positive_index)
  function(x) {
    res = F(x)
    x_not_positive = which(x <= zero_threshold)
    if (any(x_not_positive)) res[x_not_positive] = NA_real_
    res
  }
}

#' The quantile function of the Laplace distribution
#' @export
qlaplace = function(p) {
  res = p
  below_median = which(p <= .5 & p >= 0)
  above_median = which(p > .5 & p <= 1)
  if (any(below_median)) res[below_median] = log(2 * p[below_median])
  if (any(above_median)) res[above_median] = -log(2 * (1 - p[above_median]))
  res
}

#' @export
plaplace = function(x) {
  res = x
  below_zero = which(x <= 0)
  above_zero = which(x > 0)
  if (any(below_zero)) res[below_zero] = exp(x[below_zero]) / 2
  if (any(above_zero)) res[above_zero] = 1 - exp(-x[above_zero]) / 2
  res
}
