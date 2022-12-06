devtools::load_all()
library(INLA)
library(ggplot2)
library(dplyr)
library(tidyr)
library(Matrix)
library(sf)
library(pbapply)
library(parallel)
library(purrr)

# Compile and link all the cgeneric scripts, if this has not already been done
make_cgeneric("all")

n_repl = 300 # Number of replications of the experiments
n_cores = 8 # Run code in parallel (this might be too many cores for you)
threshold = qlaplace(.999) # The threshold t for defining the conditional extremes model
rho = 40 # Range of the Matérn correlation
n_posterior_samples = 1e5 # Number of samples from posterior for estimating credible interval
spde_alpha = 1.5 # Value of alpha used for the SPDE model

# File used for saving all the results
filename = file.path(results_dir(), "simulation.rds")
if (!file.exists(filename)) saveRDS(list(), filename)

# Create a regular grid of coordinates for simulating the data
n_xy = 100
coords = expand.grid(x = 1:n_xy, y = 1:n_xy) |>
  as.matrix()
n_loc = nrow(coords)

# Create the mesh and the SPDE for simulating the data
mesh = inla.mesh.2d(
  loc = coords,
  boundary = list(inla.nonconvex.hull(coords, convex = -.1),
                  inla.nonconvex.hull(coords, convex = -1)),
  cutoff = 5,
  max.edge = c(10, 50))
spde = inla.spde2.pcmatern(
  mesh,
  prior.range = c(60, .95),
  prior.sigma = c(4, .05),
  alpha = spde_alpha)

# Create the projection matrix for the mesh and SPDE
A = inla.spde.make.A(mesh, coords)

# Compute the precision matrix of the SPDE approximation mesh nodes
Q = inla.spde2.precision(spde, c(log(rho), log(1)))

# Sample several realisations of the Gaussian random field, and
# transform them to Laplace marginals
set.seed(1)
n_samples = 1e4
transformed_data = local({
  samples = rnorm_spde(n_samples, Q)
  samples = as.matrix(A %*% samples)
  samples = samples + rnorm(length(samples), sd = .1)
  samples = qlaplace(pnorm(samples))
  any_big = apply(samples, 2, function(x) any(x > threshold))
  t(samples[, any_big])
})

# ==============================================================================
# Examine patterns in the data
# ==============================================================================

# Choose thresholds for computing χ_p(d)
p = 1 - 10^seq(-3, -5, by = -.1)
thresholds = qlaplace(p)

# Create a sub-grid of conditioning sites with resolution (2 x 2)
# for computing χ and conditional moments
delta_s0 = 2
s0_index = get_s0_index(coords, delta_s0)

# This describes how to remove some of the observations far away from the conditioning site
# to give more weight to close-by observations
thinning = c(1, 2, 4, 8, 16)

# Extract data from all time points where we observe a threshold exceedance at one of our
# conditioning sites
data = extract_extreme_fields(
  data = transformed_data,
  coords = coords,
  s0_index = s0_index,
  threshold = threshold,
  n_cores = n_cores,
  n = thinning,
  r = cumsum(8 * thinning))

# Choose values of d and y_0 for estimating conditional moments and χ_p(d)
dist_table = data$dist_to_s0 |> unlist() |> round() |> table()
unique_dist = dist_table[which(dist_table > 4000)] |> names() |> as.numeric()
y0_table = data$y0 |> unlist() |> round(1) |> table()
unique_y0 = y0_table[which(y0_table > 40)] |> names() |> as.numeric()

# Compute conditional moments and χ_p(d)
chi = empirical_chi(data, thresholds, unique_dist)
moments = empirical_moments(data, unique_dist, unique_y0)

#plot_chi(chi)
#plot_moments(moments)

# ==============================================================================
# Compute the MLE of our chosen model
# ==============================================================================

# Create a sub-grid of conditioning sites with resolution (6 x 6)
# for computing the MLE
delta_s0 = 6
s0_index = get_s0_index(coords, delta_s0)

# Extract data from all time points where we observe a threshold exceedance at one of our
# conditioning sites
data = extract_extreme_fields(
  data = transformed_data,
  coords = coords,
  s0_index = s0_index,
  threshold = threshold,
  n_cores = n_cores,
  n = thinning,
  r = cumsum(4 * thinning))

# This is the composite log-likelihood of our chosen model
loglik = function(theta,
                  y,
                  y0,
                  dist_to_s0,
                  dist_to_s0_from_mesh,
                  A,
                  spdes,
                  sum_terms = TRUE,
                  n_cores = 1) {
  stopifnot(all(sapply(dist_to_s0_from_mesh, min) == 0))

  # All model parameters
  tau = exp(theta[1])
  lambda = exp(theta[2])
  kappa = exp(theta[3])
  log_sigma = theta[4]
  log_rho = theta[5]
  beta0 = exp(theta[6]) / (1 + exp(theta[6]))
  lambda_b = exp(theta[7])
  kappa_b = exp(theta[8])

  # Compute precision matrices for all triangulated meshes
  Q = lapply(spdes, INLA::inla.spde2.precision, theta = c(log_rho, log_sigma))

  # Constrain all the precision matrices
  for (i in seq_along(Q)) {
    ii = which(dist_to_s0_from_mesh[[i]] == 0)
    Q[[i]] = Q[[i]][-ii, -ii]
    A[[i]] = A[[i]][, -ii]
    dist_to_s0_from_mesh[[i]] = dist_to_s0_from_mesh[[i]][-ii]
  }

  # Functions for computing a(d; y_0) and b(d; y_0)
  a_func = function(y, dist) {
    alpha = exp(- (dist / lambda)^kappa)
    a = rep(y, each = length(alpha)) * rep(alpha, length(y))
    matrix(a, nrow = length(dist), ncol = length(y))
  }
  b_func = function(y, dist) {
    beta = beta0 * exp(-(dist / lambda_b)^kappa_b)
    b = rep(y, each = length(beta)) ^ rep(beta, length(y))
    matrix(b, nrow = length(dist), ncol = length(y))
  }

  # Compute the composite log-likelihood
  res = loglik_conditional(
    y = y,
    y0 = y0,
    a_func = a_func,
    b_func = b_func,
    Q = Q,
    tau = tau,
    dist_to_s0 = dist_to_s0,
    dist_to_s0_from_mesh = dist_to_s0_from_mesh,
    A = A,
    n_cores = n_cores)

  if (sum_terms) res = sum(res)
  res
}

# Compute the gradient of the log-likelihood using numerical derivation
loglik_grad = function(theta, sum_terms = TRUE, ...) {
  res = numDeriv::jacobian(loglik, theta, ..., sum_terms = FALSE)
  if (sum_terms) res = apply(res, 2, sum)
  res
}

# Create a triangulated mesh for every of our chosen conditioning sites.
# The function INLA::inla.mesh.2d() creates triangulated meshes, but it will
# fail to create a mesh and run indefinitely for certain combinations of input arguments.
# When creating many different meshes, it is really hard to ensure that INLA::inla.mesh.2d()
# will work at a first try for all of them. Therefore, we need to create the meshes
# using a while loop that slightly changes the input arguments if the mesh creation fails
multimesh_data = parallel::mclapply(
  X = seq_along(data$obs_index),
  mc.cores = n_cores,
  mc.preschedule = FALSE,
  FUN = function(i) {
    convex = 40 # The original convexity used for creating the mesh boundary
    while (TRUE) {
      # The original mesh arguments
      args = list(
        loc = coords[c(data$obs_index[[i]], data$s0_index[[i]]), ],
        boundary = inla.nonconvex.hull(coords[data$obs_index[[i]], ], convex = convex),
        max.edge = 50)
      # Try to create a mesh using a function that automatically fails if the mesh
      # creation takes more than `timeout` seconds
      mesh = create_mesh(args, timeout = 5)
      if (!is.null(mesh)) break # Break out of the loop if this succeeded
      convex = convex + 5 # If not, increase the convexity range and try again
    }
    # Create an SPDE object, a projection matrix and the distance from all mesh nodes to
    # the conditioning site for the created mesh
    spde = inla.spde2.pcmatern(
      mesh,
      prior.range = c(60, .95),
      prior.sigma = c(4, .05),
      alpha = spde_alpha)
    dist_to_s0_from_mesh = dist_euclid(mesh$loc[, 1:2], data$s0[[i]])
    A = inla.spde.make.A(mesh, coords[data$obs_index[[i]], ])
    message("Completed mesh nr. ", i, " of ", length(data$obs_index))
    list(spde = spde, dist = dist_to_s0_from_mesh, A = A)
  })
multimesh_data = purrr::transpose(multimesh_data)

# Plot an example of the locations used for performing inference for a given
# conditioning site, together with a display of the triangular mesh
plot = local({
  i = 105
  df = as.data.frame(coords[data$obs_index[[i]], ]) |>
    dplyr::mutate(tag = "Used")
  df2 = as.data.frame(coords[-c(data$obs_index[[i]], data$s0_index[[i]]), ]) |>
    dplyr::mutate(tag = "Not used")
  s0_df = as.data.frame(coords[data$s0_index[[i]], , drop = FALSE]) |>
    dplyr::mutate(tag = "$\\bm s_0$")
  rbind(df, df2, s0_df) |>
    dplyr::mutate(tag = factor(tag, levels = c("$\\bm s_0$", "Used", "Not used"))) |>
    ggplot() +
    inlabru::gg(
      multimesh_data$spde[[i]]$mesh,
      int.color = "black",
      ext.color = "black",
      edge.color = "black") +
    geom_point(aes(x = x, y = y, shape = tag, size = tag, color = tag, alpha = tag)) +
    scale_color_manual(values = c("red", "black", "black")) +
    scale_shape_manual(values = c(17, 19, 19)) +
    scale_size_manual(values = c(4, 4, .5)) +
    scale_alpha_manual(values = c(1, 1, .7)) +
    theme_light() +
    guides(shape = "none", size = "none", color = "none", alpha = "none") +
    labs(x = "Easting", y = "Northing") +
    theme(
      line = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_text(size = 45),
      panel.grid = element_blank())
})

tikz_plot(file.path(image_dir(), "design-of-experiment.pdf"),
          plot, width = 15, height = 15)

# # Convert the heavy pdf file to a lighter compressed jpeg file
# library(magick)
# image = magick::image_read_pdf(file.path(image_dir(), "design-of-experiment.pdf"), density = 75)
# magick::image_write(
#   image = image,
#   path = file.path(image_dir(), "design-of-experiment.jpg"),
#   format = "jpg",
#   quality = 50)

# Compute the MLE
# ==============================================================================

# Initial values for the MLE computation
est = list(par = c(2, log(20), log(1), log(1.4), log(40), 0, log(20), log(1)), convergence = 1)
# Keep on optimising until we find a local or global maximum
while (est$convergence != 0) {
  est = optim(
    par = est$par,
    fn = loglik,
    y = data$y,
    y0 = data$y0,
    dist_to_s0 = data$dist_to_s0,
    dist_to_s0_from_mesh = multimesh_data$dist,
    A = multimesh_data$A,
    spdes = multimesh_data$spde,
    n_cores = 5,
    control = list(fnscale = -1, maxit = 500, trace = 6))
}

# Save the MLE
tmp = readRDS(filename)
tmp$mle = est
saveRDS(tmp, filename)
rm(tmp)

# ==============================================================================
# Fit the model to data using INLA
# ==============================================================================
est = readRDS(filename)$mle

# R-INLA requires the observations y to be on a vector format
y_inla = unlist(data$y)

# Our implementation of the cgeneric model for a(d; y_0) requires y0 and dist_to_s0 as input.
# However, it requires one value of y0 and dist_to_s0 for each of the observations y.
# We therefore replicate y0 and dist_to_s0 in the correct way so we get to vectors with
# equal length to y_inla, where y_inla[i] was observed when y(s0) was equal to
# y0_inla[i], and the distance between them was dist_to_s0_inla[i].
dist_to_s0_inla = data$dist_to_s0[rep(seq_along(data$n), data$n)]
y0_inla = rep(unlist(data$y0), sapply(dist_to_s0_inla, length))
dist_to_s0_inla = unlist(dist_to_s0_inla)

# Define the cgeneric model for a with priors.
# Use the MLE as initial values
a_priors = list(lambda = c(3, 4), kappa = c(0, 3))
a_model = a_model(
  y0 = y0_inla,
  dist_to_s0 = dist_to_s0_inla,
  init = est$par[2:3],
  priors = a_priors)

# Define the cgeneric model for Z_b with priors.
# Use the MLE as initial values
spde_priors = list(
  sigma = c(4, .05),
  rho = c(60, .95),
  beta0 = c(0, 2),
  lambda = c(3, 4),
  kappa = c(0, 3))
spde_model = spde_b_model_simulation(
  n = data$n,
  y0 = unlist(data$y0),
  spde = multimesh_data$spde,
  init = est$par[4:8],
  priors = spde_priors,
  dist_to_s0 = multimesh_data$dist)

# Create the data object and projection matrix necessary for performing inference.
# These are typically created using the INLA::inla.stack() function, but this function
# is not created for the type of cgeneric model we are implementing
inla_data = local({
  # Compute projection matrices A for each pair of a mesh and the observations used for inference
  A = list()
  pb = progress_bar(length(data$n))
  for (i in seq_along(data$n)) {
    location_index = rep(data$obs_index[[i]], data$n[i])
    A[[i]] = inla.spde.make.A(
      multimesh_data$spde[[i]]$mesh,
      loc = coords[data$obs_index[[i]], ],
      index = rep(seq_along(data$obs_index[[i]]), data$n[i]),
      repl = rep(seq_len(data$n[i]), each = length(data$obs_index[[i]])))
    # Account for the constraining method by removing columns corresponding to the conditioning site
    n_spde = multimesh_data$spde[[i]]$n.spde
    ii = which(multimesh_data$dist[[i]] == 0)
    ii = (0:(data$n[i] - 1)) * n_spde + ii
    A[[i]] = A[[i]][, -ii]
    pb$tick()
  }
  pb$terminate()
  A = Matrix::bdiag(A)

  # Remove all columns of A that don't contain any elements
  n_per_col = tail(A@p, -1) - head(A@p, -1)
  ii = which(n_per_col > 0)
  A = A[, ii]

  # Add a diagonal matrix to A for adding contributions to the linear predictor from a(d; y_0)
  A = cbind(A, Diagonal(length(y_inla), 1))

  # Create the data object used for performing inference with INLA
  data = list(
    y = y_inla,
    spatial = c(ii, rep(NA, length(y_inla))),
    idx = c(rep(NA, length(ii)), seq_along(y_inla)))

  list(A = A, data = data)
})

# Define the model formula for performing inference
formula = y ~ -1 +
  f(idx, model = a_model) +
  f(spatial, model = spde_model)

# Fit the model
fit = inla(
  formula = formula,
  data = inla_data$data,
  control.predictor = list(A = inla_data$A),
  control.family = list(hyper = list(prec = list(
    initial = est$par[1],
    prior = "pc.prec",
    param = c(1, .05)))),
  only.hyperparam = TRUE,
  control.inla = list(
    restart = 1,
    control.vb = list(enable = FALSE),
    int.strategy = "eb"),
  control.compute = list(config = TRUE),
  verbose = TRUE,
  num.threads = 1,
  inla.mode = "experimental")

# Rerun INLA to ensure that we have a good fit
fit = inla.rerun(fit)

# Keep only the relevant parts of the model fit, to reduce requirements on file storage
fit = list(misc = fit$misc,
           internal.marginals.hyperpar = fit$internal.marginals.hyperpar,
           mode = fit$mode,
           cpu = fit$cpu,
           summary.hyperpar = fit$summary.hyperpar)
fit$misc$reordering = NULL
fit$misc$family = NULL
fit$misc$linkfunctions = NULL
class(fit) = "inla"

# Save the results
tmp = readRDS(filename)
tmp$inla = fit
saveRDS(tmp, filename)
rm(tmp)

# ==============================================================================
# Compute the C matrix for performing post hoc adjustments
# ==============================================================================
fit = readRDS(filename)$inla

# Compute the Hessian matrix
H = solve(fit$misc$cov.intern)

# Compute numerical gradients for each of the composite log-likelihood terms
grads = loglik_grad(
  theta = fit$mode$theta,
  y = data$y,
  y0 = data$y0,
  dist_to_s0 = data$dist_to_s0,
  dist_to_s0_from_mesh = multimesh_data$dist,
  A = multimesh_data$A,
  n_cores = n_cores,
  spdes = multimesh_data$spde,
  sum_terms = FALSE)

# Estimate J by summing over all correlated pairs.
# Here, we know that observations from different time points are independent
times = unlist(data$time_index)
J = 0
for (k in seq_along(times)) {
  index = which(times == times[k])
  for (j in index) {
    J = J + grads[k, ] %*% grads[j, , drop = FALSE]
  }
}

# Compute the estimate for C
C = get_C(H, J)

# Save the results
tmp = readRDS(filename)
tmp$C = C
saveRDS(tmp, filename)
rm(tmp)

# ==============================================================================
# Simulate data from the unadjusted and adjusted model fits
# ==============================================================================

tmp = readRDS(filename)
C = tmp$C
fit = tmp$inla
rm(tmp)

# Simulate parameters from the unadjusted model fit and adjust them using
# the estimates of θ^* and C
set.seed(1)
n_theta = length(fit$mode$theta)
samples = list()
samples$unadjusted = inla.hyperpar.sample(
  n = n_posterior_samples,
  result = fit,
  intern = TRUE)
samples$adjusted = matrix(
  data = rep(fit$mode$theta, each = n_posterior_samples),
  ncol = n_theta)
samples$adjusted = samples$adjusted + (samples$unadjusted - samples$adjusted) %*% t(C)
for (i in seq_along(samples)) {
  samples[[i]] = exp(samples[[i]])
  samples[[i]][, 6] = samples[[i]][, 6] / (1 + samples[[i]][, 6])
}

# Plot the adjusted and unadjusted posteriors for all model parameters
theta_names = paste0(
  "$\\", c("tau", "lambda_a", "kappa_a", "sigma", "rho", "beta_0", "lambda_b", "kappa_b"), "$")
plot = local({
  df1 = as.data.frame(samples$unadjusted)
  names(df1) = theta_names
  df1$adjusted = "Unadjusted"
  df2 = as.data.frame(samples$adjusted)
  names(df2) = theta_names
  df2$adjusted = "Adjusted"
  rbind(df1, df2) |>
    tidyr::pivot_longer(-adjusted) |>
    dplyr::mutate(
      name = factor(name, levels = theta_names),
      adjusted = factor(adjusted, levels = c("Unadjusted", "Adjusted"))) |>
    ggplot() +
    geom_density(aes(x = value, linetype = adjusted)) +
    facet_wrap(~name, scales = "free", nrow = 2) +
    labs(linetype = "", col = "", x = "Parameter value", y = "Posterior density") +
    theme_light() +
    scale_linetype_manual(values = c("longdash", "solid")) +
    theme(
      text = element_text(size = 15),
      strip.text = element_text(size = 15, colour = "black"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      legend.position = "top")
})

# Export the plot to pdf format
tikz_plot(file.path(image_dir(), "simulation-posteriors.pdf"),
          plot, width = 10, height = 6)

# ==============================================================================
# Simulate data from the two model fits, and compute properties of the
# simulated data
# ==============================================================================

# Function for simulating data from a fitted model.
# Given a `samples` object of model parameter samples, as created further up,
# we use the first `n_samples` of these and simulate `n_samples_per_theta`
# realisations from the spatial conditional extremes model for each
# model parameter sample.
simulate_data = function(samples, n_samples, n_samples_per_theta) {
  n_samples = min(n_samples, nrow(samples)) # Ensure that n_samples <= nrow(samples)

  # Define the location we will use as conditioning site for all our simulated data
  my_s0_index = s0_index[105]

  # Simulate all threshold exceedances
  y0 = rexp(n_samples * n_samples_per_theta) + threshold

  # Define the coordinates where we will observe the simulated data
  my_coords = extract_thinned_out_circles(
    coords = coords,
    center = coords[my_s0_index, ],
    n = c(1, 2),
    r = c(6, Inf))

  # Create a mesh, an SPDE object and a projection matrix, for simulating data
  mesh = inla.mesh.2d(
    loc = my_coords,
    boundary = list(
      inla.nonconvex.hull(my_coords, convex = -.2),
      inla.nonconvex.hull(my_coords, convex = -1)),
    max.edge = c(10, 50))
  spde = inla.spde2.pcmatern(
    mesh = mesh,
    alpha = spde_alpha,
    prior.range = c(exp(fit$mode$theta[5]), .5),
    prior.sigma = c(exp(fit$mode$theta[4]), .5))
  A = inla.spde.make.A(mesh, my_coords)

  # Compute distances to the conditioning site
  dist_to_s0 = dist_euclid(coords[my_s0_index, ], my_coords)
  dist_to_s0_from_mesh = dist_euclid(coords[my_s0_index, ], mesh$loc[, -3])

  # Remove elements from the projection matrix that correspond to the conditioning site
  ii = which(dist_to_s0_from_mesh == 0)
  jj = which(dist_to_s0 == 0)
  dist_to_s0_from_mesh = dist_to_s0_from_mesh[-ii]
  dist_to_s0 = dist_to_s0[-jj]
  A = A[-jj, -ii]

  # Create a data-object of same type as that returned by the function extract_extreme_fields(),
  # for saving the simulated realisations
  data = list(
    y0 = list(y0),
    dist_to_s0 = list(dist_to_s0),
    n = n_samples * n_samples_per_theta)

  # Simulate from the fitted conditional extremes model
  data$y = pbapply::pblapply(
    X = 1:n_samples,
    cl = n_cores,
    FUN = function(i) {
      # These are the parameters from model-parameter sample nr. i
      params = samples[i, ]
      tau = params[1]
      lambda_a = params[2]
      kappa_a = params[3]
      sigma = params[4]
      rho = params[5]
      beta0 = params[6]
      lambda_b = params[7]
      kappa_b = params[8]

      # Create the functions for a(d; y_0) and b(d; y_0)
      a_func = function(y, dist) {
        alpha = exp(- (dist / lambda_a)^kappa_a)
        a = rep(y, each = length(alpha)) * rep(alpha, length(y))
        matrix(a, nrow = length(dist), ncol = length(y))
      }
      b_func = function(y, dist) {
        beta = beta0 * exp(-(dist / lambda_b)^kappa_b)
        b = rep(y, each = length(beta)) ^ rep(beta, length(y))
        matrix(b, nrow = length(dist), ncol = length(y))
      }

      # Create and constrain the precision matrix
      Q = inla.spde2.precision(spde, log(c(rho, sigma)))
      Q = Q[-ii, -ii]

      # Simulate from the spatial conditional distribution
      rconditional(
        y0 = list(y0[(i - 1) * n_samples_per_theta + 1:n_samples_per_theta]),
        a_func = a_func,
        b_func = b_func,
        Q = Q,
        tau = tau,
        dist_to_s0 = list(dist_to_s0),
        dist_to_s0_from_mesh = list(dist_to_s0_from_mesh),
        A = list(A))[[1]]
    })

  data$y = list(do.call(cbind, data$y))
  data
}

n_samples = 1e3 # The number of different model parameter samples we will use
n_samples_per_theta = 100 # The number of data realisations to simulate for each model parameter sample

# Simulate data from the unadjusted model
mydata = simulate_data(samples$unadjusted, n_samples, n_samples_per_theta)
# Estimate χ and conditional moments from the simulated data
unadjusted_chi = empirical_chi(mydata, attr(chi, "thresholds"), attr(chi, "unique_dist"))
unadjusted_moments = empirical_moments(
  data = mydata,
  unique_dist = attributes(moments)$unique_dist,
  unique_y0 = attributes(moments)$unique_y0)

# Simulate data from the adjusted model
mydata = simulate_data(samples$adjusted, n_samples, n_samples_per_theta)
# Estimate χ and conditional moments from the simulated data
adjusted_chi = empirical_chi(mydata, attr(chi, "thresholds"), attr(chi, "unique_dist"))
adjusted_moments = empirical_moments(
  data = mydata,
  unique_dist = attributes(moments)$unique_dist,
  unique_y0 = attributes(moments)$unique_y0)

rm(mydata)

# Plot the estimates of χ from the original data and the two model fits
chi_plot = local({
  nan_index = which(is.nan(adjusted_chi[, 1]))
  thresholds = attr(chi, "thresholds")
  p = plaplace(thresholds)
  unique_dist = attr(chi, "unique_dist")[-nan_index]
  tmp = list(Observations = chi, "Adjusted fit" = adjusted_chi, "Unadjusted fit" = unadjusted_chi)
  for (i in seq_along(tmp)) {
    tmp[[i]] = data.frame(
      chi = as.numeric(rbind(1, tmp[[i]][-nan_index, ])),
      dist = rep(c(0, unique_dist), length(p)),
      p = rep(-log10(1 - p), each = length(unique_dist) + 1),
      tag = names(tmp)[i])
  }
  do.call(rbind, tmp) |>
    dplyr::mutate(
      tag = factor(tag, levels = c("Observations", "Adjusted fit", "Unadjusted fit"))) |>
    ggplot() +
    geom_line(aes(x = dist, y = chi, group = p, col = p)) +
    labs(col = "$p$", y = "$\\widehat \\chi_p(d)$", x = "$d$") +
    scale_color_gradient(
      breaks = 3:5,
      limits = c(2.9, 5.1),
      labels = paste0("$1 - 10^{-", 3:5, "}$"),
      low = "black",
      high = "lightgray") +
    theme_light() +
    theme(axis.title.y = element_text(angle = 0, vjust = .5),
          strip.text = element_text(colour = "black"),
          strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")) +
    facet_wrap(~tag, nrow = 1) +
    lims(x = c(0, 40))
})

# Plot the estimates of conditional moments from the original data and the two model fits
moment_plots = local({
  nan_index = which(is.nan(adjusted_moments$mean[, 1]))
  unique_dist = attr(moments, "unique_dist")[-nan_index]
  unique_y0 = attr(moments, "unique_y0")
  tmp = list(
    Observations = moments,
    "Adjusted fit" = adjusted_moments,
    "Unadjusted fit" = unadjusted_moments)
  for (i in seq_along(tmp)) {
    tmp[[i]] = data.frame(
        sd = as.numeric(rbind(0, tmp[[i]]$sd[-nan_index, ])),
        mean = as.numeric(rbind(unique_y0, tmp[[i]]$mean[-nan_index, ])),
        dist = rep(c(0, unique_dist), length(unique_y0)),
        y0 = rep(unique_y0, each = length(unique_dist) + 1),
        tag = names(tmp)[i])
  }
  sd_plot = do.call(rbind, tmp) |>
    dplyr::mutate(
      tag = factor(tag, levels = c("Observations", "Adjusted fit", "Unadjusted fit"))) |>
    ggplot() +
    geom_line(aes(x = dist, y = sd, group = y0, col = y0)) +
    scale_color_gradient(low = "black", high = "lightgray") +
    labs(y = "$\\widehat \\zeta(d; y_0)$", col = "$y_0$", x = "$d$") +
    theme_light() +
    theme(axis.title.y = element_text(angle = 0, vjust = .5),
          strip.text = element_text(colour = "black"),
          strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")) +
    facet_wrap(~tag, nrow = 1) +
    lims(x = c(0, 60))
  mean_plot = do.call(rbind, tmp) |>
    dplyr::mutate(
      tag = factor(tag, levels = c("Observations", "Adjusted fit", "Unadjusted fit"))) |>
    ggplot() +
    geom_line(aes(x = dist, y = mean, group = y0, col = y0)) +
    labs(y = "$\\widehat \\mu(d; y_0)$", col = "$y_0$", x = "$d$") +
    scale_color_gradient(low = "black", high = "lightgray") +
    theme_light() +
    theme(axis.title.y = element_text(angle = 0, vjust = .5),
          strip.text = element_text(colour = "black"),
          strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")) +
    facet_wrap(~tag, nrow = 1) +
    lims(x = c(0, 60))
  list(mean = mean_plot, sd = sd_plot)
})

plot = patchwork::wrap_plots(chi_plot, moment_plots$mean, moment_plots$sd, ncol = 1)
plot = plot * theme(text = element_text(size = 16))
tikz_plot(file.path(image_dir(), "simulation-properties.pdf"), plot, width = 11, height = 8)

# ==============================================================================
# Compute out-of-sample log-scores
# ==============================================================================

# Simulate a new data set from the true spatial Gaussian random field
set.seed(1)
n_data_samples = 5e4
transformed_data = local({
  samples = rnorm_spde(n_data_samples, Q)
  samples = as.matrix(A %*% samples)
  samples = samples + rnorm(length(samples), sd = .1)
  samples = qlaplace(pnorm(samples))
  big_index = apply(samples, 2, function(x) any(x > threshold))
  res = t(samples[, big_index])
  attr(res, "big_index") = big_index
  res
})

# Extract extremes using the same design as we used for performing inference
data = extract_extreme_fields(
  data = transformed_data,
  coords = coords,
  s0_index = s0_index,
  threshold = threshold,
  n = thinning,
  r = cumsum(4 * thinning))

n_posterior_samples = 1e3 # The number n_s from the paper
# Compute composite log-likelihoods for n_s different model parameter samples using
# both the unadjusted and the adjusted model fit
ll = pbapply::pblapply(
  X = 1:n_posterior_samples,
  cl = n_cores,
  FUN = function(i) {
    ll = list()
    for (name in names(samples)) {
      param = samples[[name]][i, ]
      param[-6] = log(param[-6])
      param[6] = log(param[6]) - log(1 - param[6])
      ll[[name]] = loglik(
        theta = param,
        y = data$y,
        y0 = data$y0,
        dist_to_s0 = data$dist_to_s0,
        dist_to_s0_from_mesh = multimesh_data$dist[s0_index %in% data$s0_index],
        A = multimesh_data$A[s0_index %in% data$s0_index],
        sum_terms = FALSE,
        spdes = multimesh_data$spde[s0_index %in% data$s0_index])
    }
    ll
  })

# Save the results
tmp = readRDS(filename)
tmp$loglik = ll
saveRDS(tmp, filename)
rm(tmp)

# Compute the log-scores and compare them
tmp = purrr::transpose(ll)
tmp = sapply(tmp, function(x) {
  x = do.call(rbind, x)
  apply(x, 2, log_mean_exp)
  })

sum(tmp[, 2]) - sum(tmp[, 1]) # The adjusted model fit performs better

# Perform bootstrapping to see if the difference in log-score is significant
set.seed(1)
n_bootstraps = 5000
pb = progress_bar(n_bootstraps)
big_index = which(attr(transformed_data, "big_index"))
bootstrap_sums = sapply(
  X = seq_len(n_bootstraps),
  FUN = function(i) {
    index = sample.int(n_data_samples, n_data_samples, replace = TRUE)
    time_index = unlist(data$time_index)
    log_scores = lapply(
      X = unique(big_index[time_index]),
      FUN = function(j) {
        n = sum(index == j)
        ii = which(big_index[time_index] == j)
        if (length(ii) > 0 && n > 0) {
          do.call(rbind, rep(list(tmp[ii, ]), n))
        } else {
          NULL
        }
      })
    log_scores = do.call(rbind, log_scores)
    pb$tick()
    apply(log_scores, 2, sum)
  })
pb$terminate()

summary(bootstrap_sums[2, ] - bootstrap_sums[1, ]) # Yes, it is
