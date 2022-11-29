devtools::load_all()
library(microbenchmark)
library(sf)
library(dplyr)
library(tidyr)
library(INLA)
library(pbapply)

# In this script we evaluate the speed and correctness of the Rcpp functions used
# for computing log-likelihoods for the spatial conditional extremes model.
# The evaluation is performed using data from the case study, where we
# compute log-likelihoods both using the Rcpp functions, and using functions written
# in plain R. The values and computation speeds are then compared for the R and Rcpp functions.

# ==============================================================================
# Decide necessary model parameters:
# ==============================================================================
threshold = qlaplace(.95) # The threshold t for defining the conditional extremes model
n_cores = 15 # Run code in parallel. Check yourself that this is an ok number
r = 5 # Radius used for computing aggregated empirical distribution functions

# ==============================================================================
# Load the radar data
# ==============================================================================
radar = readRDS(file.path(downloads_dir(), "radar.rds"))
coords = st_coordinates(radar$coords)
n_loc = nrow(coords)

# ==============================================================================
# Transform the data to Laplace margins
# ==============================================================================
cl = parallel::makeForkCluster(n_cores)
transformed_data = pbapply::pblapply(
  X = 1:n_loc,
  cl = cl,
  FUN = function(i) {
    F = aggregated_ecdf(radar$data, coords, coords[i, ], radius = r)
    u = F(radar$data[, i])
    # Ensure that we don't get infinities
    if (any(u == 1, na.rm = TRUE)) {
      u[which(u == 1)] = (1 + max(u[which(u != 1)])) / 2
    }
    qlaplace(u)
  })
parallel::stopCluster(cl)
transformed_data = do.call(cbind, transformed_data)

# ==============================================================================
# Choose all conditioning sites for the global model
# ==============================================================================

# Create a grid of conditioning sites with resolution (delta_s0 x delta_s0)
delta_s0 = 4
x_coords = coords[, 1] |> unique() |> sort()
y_coords = coords[, 2] |> unique() |> sort()
s0_locs = expand.grid(
  x = x_coords[seq(delta_s0, length(x_coords), by = delta_s0)],
  y = y_coords[seq(delta_s0, length(y_coords), by = delta_s0)]) |>
  as.matrix()

# Find out which indices the chosen conditioning sites correspond to
s0_index = lapply(
  X = 1:nrow(s0_locs),
  FUN = function(i) which(coords[, 1] == s0_locs[i, 1] & coords[, 2] == s0_locs[i, 2]))
s0_index = s0_index[sapply(s0_index, length) > 0] |>
  unlist() |>
  unname()

# ==============================================================================
# Extract observations in a nonregular grid around conditioning sites that
# exceed the threshold
# ==============================================================================
thinning = c(1, 2, 4, 6, 8, 16, 32)
data = extract_extreme_fields(
  data = transformed_data,
  coords = coords,
  s0_index = s0_index,
  threshold = threshold,
  n = thinning,
  r = cumsum(4 * thinning))

# ==============================================================================
# Create the mesh and spde
# ==============================================================================
mesh = inla.mesh.2d(
  loc = coords,
  boundary = list(inla.nonconvex.hull(coords, convex = -.1),
                  inla.nonconvex.hull(coords, convex = -1)),
  max.edge = c(5, 30))
mesh$n
# plot(mesh)
# points(coords)

spde = inla.spde2.pcmatern(mesh, prior.range = c(60, .95), prior.sigma = c(5, .05))

# Compute the distances to the conditioning sites from the mesh nodes,
# and add the information to `data`
data$dist_to_s0_from_mesh = list()
for (i in seq_along(data$s0)) {
  data$dist_to_s0_from_mesh[[i]] = dist_euclid(mesh$loc[, 1:2], data$s0[[i]])
}

# ==============================================================================
# Compute log-likelihoods and compare computation times
# ==============================================================================

# Define the composite log-likelihood for the global conditional
# extremes model. The function below is a wrapper for the function loglik_conditional
#
# Input variables:
# theta: A vector of parameters of length 5 or 6, depending on the value of rho_b.
# y, y0, dist_to_s0, dist_to_s0_from_mesh, A, n_cores: See ?loglik_conditional for more info
# rho_b: Either NULL or a double. If is.null(rho_b), then theta has length 6, and
#   we try to estimate rho_b. Else, theta has length 5, and rho_b is fixed and equal
#   to the given value.
# sum_terms: A boolean describing wheter we should compute the sum of the log-likelihood,
#   or if we should remove one value for each threshold exceedance.
# no_beta: A boolean. Should we compute the composite log-likelihood using the
#   dconditional(_r)_no_beta() versions or using the dconditional(_r)() versions?
loglik = function(theta,
                  y,
                  y0,
                  dist_to_s0,
                  dist_to_s0_from_mesh,
                  A,
                  use_r,
                  rho_b = NULL,
                  sum_terms = TRUE,
                  no_beta = TRUE,
                  n_cores = 1) {
  if (is.null(rho_b)) {
    stopifnot(length(theta) == 6)
    rho_b = exp(theta[3])
    log_rho = theta[4]
    log_sigma = theta[5]
    tau = exp(theta[6])
  } else {
    stopifnot(length(theta) == 5)
    log_rho = theta[3]
    log_sigma = theta[4]
    tau = exp(theta[5])
  }
  lambda = exp(theta[1])
  kappa = exp(theta[2])
  Q = INLA::inla.spde2.precision(spde, c(log_rho, log_sigma))
  a_func = function(y, dist) {
    alpha = exp(- (dist / lambda)^kappa)
    matrix(rep(y, each = length(alpha)) * rep(alpha, length(y)),
           nrow = length(dist), ncol = length(y))
  }
  b_func = function(y, dist) {
    tmp = dist / rho_b
    tmp[tmp < 1e-9] = 1e-9
    b = sqrt(1 - exp(-2 * tmp))
    res = matrix(rep(b, length(y)), nrow = length(dist), ncol = length(y))
    if (!no_beta) {
      # Add some small random noise so we won't use the dconditional(_r)_no_beta() functions
      set.seed(1)
      res = res + rnorm(length(res), sd = 1e-8)
    }
    res
  }
  res = loglik_conditional2(
    y = y,
    y0 = y0,
    a_func = a_func,
    b_func = b_func,
    Q = Q,
    tau = tau,
    dist_to_s0 = dist_to_s0,
    dist_to_s0_from_mesh = dist_to_s0_from_mesh,
    A = A,
    n_cores = n_cores,
    use_r = use_r)
  if (sum_terms) res = sum(res)
  res
}

# Prerequisites for computing the log-likelihoods
theta = c(3.8, -.4, -.5, 2.3, .3, 4.2)
A = lapply(data$obs_index, function(x) inla.spde.make.A(mesh, coords[x, ]))
args = list(
  theta = theta,
  y = data$y,
  y0 = data$y0,
  dist_to_s0 = data$dist_to_s0,
  dist_to_s0_from_mesh = data$dist_to_s0_from_mesh,
  A = A,
  sum_terms = FALSE)

# Compute the log-likelihoods for all combinations of use_r and no_beta
results = list()
for (no_beta in c(TRUE, FALSE)) {
  for (use_r in c(TRUE, FALSE)) {
    args$use_r = use_r
    args$no_beta = no_beta

    time = proc.time()
    ll = do.call(loglik, args)
    time = proc.time() - time

    name = paste0("use_r = ", use_r, ", no_beta = ", no_beta)
    results[[name]] = list(ll = ll, time = time)
  }
}

# Compare log-scores. We see that the Rcpp functions and the R functions return the same values
ll = sapply(results, `[[`, "ll")
apply(ll[, 1:2], 1, function(x) max(x) - min(x)) |>
  summary()
apply(ll[, 3:4], 1, function(x) max(x) - min(x)) |>
  summary()

# Compare the computation times. The Rcpp functions are considerably faster than the R functions
times = sapply(results, `[[`, "time")
times
