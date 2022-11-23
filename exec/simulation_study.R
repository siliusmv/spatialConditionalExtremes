devtools::load_all()
library(INLA)
library(ggplot2)
library(dplyr)
library(Matrix)
library(sf)
library(pbapply)
library(parallel)
library(mvnfast)

# Compile and link all the cgeneric scripts, if this has not already been done
make_cgeneric("all")

n_repl = 300 # Number of replications of the experiments
n_cores = 8 # Run code in parallel (this might be too many cores for you)
threshold = qlaplace(.999) # The threshold t for defining the conditional extremes model
rho = 40 # Range of the Matérn correlation
n_posterior_samples = 1e5 # Number of samples from posterior for estimating credible interval
spde_alpha = 1.5
thinning = c(1, 2, 4, 8, 16)
my_index = 105 # Used for fitting a single-site model to data, and performing some plotting

filename = file.path(results_dir(), "simulation.rds")
if (!file.exists(filename)) saveRDS(list(), filename)
#small_fits_filename = file.path(results_dir(), "simulation_small_fits.rds")
#if (!file.exists(small_fits_filename)) saveRDS(list(), small_fits_filename)

# Create a regular grid of coordinates, similar to the case study
n_xy = 100
coords = expand.grid(x = 1:n_xy, y = 1:n_xy) |>
  as.matrix()
n_loc = nrow(coords)

# Create the mesh and the SPDE
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

p = 1 - 10^seq(-3, -5, by = -.1)
thresholds = qlaplace(p)

delta_s0 = 2
s0_index = get_s0_index(coords, delta_s0)

data = extract_extreme_fields(
  data = transformed_data,
  coords = coords,
  s0_index = s0_index,
  threshold = threshold,
  n_cores = n_cores,
  n = thinning,
  r = cumsum(8 * thinning))

dist_table = data$dist_to_s0 |> unlist() |> round() |> table()
unique_dist = dist_table[which(dist_table > 4000)] |> names() |> as.numeric()
y0_table = data$y0 |> unlist() |> round(1) |> table()
unique_y0 = y0_table[which(y0_table > 40)] |> names() |> as.numeric()

chi = empirical_chi(data, thresholds, unique_dist)
moments = empirical_moments(data, unique_dist, unique_y0)

#plot_chi(chi)
#plot_moments(moments)

# ==============================================================================
# Estimate the shape of a() and b() using MLE
# ==============================================================================

delta_s0 = 6
s0_index = get_s0_index(coords, delta_s0)

data = extract_extreme_fields(
  data = transformed_data,
  coords = coords,
  s0_index = s0_index,
  threshold = threshold,
  n_cores = n_cores,
  n = thinning,
  r = cumsum(4 * thinning))

loglik = function(theta,
                  y,
                  y0,
                  dist_to_s0,
                  dist_to_s0_from_mesh,
                  A,
                  spdes,
                  use_beta = TRUE,
                  verbose = FALSE,
                  sum_terms = TRUE,
                  n_cores = 1) {
  stopifnot(all(sapply(dist_to_s0_from_mesh, min) == 0))
  tau = exp(theta[1])
  lambda = exp(theta[2])
  kappa = exp(theta[3])
  log_sigma = theta[4]
  log_rho = theta[5]
  beta0 = exp(theta[6]) / (1 + exp(theta[6]))
  lambda_b = exp(theta[7])
  kappa_b = exp(theta[8])
  Q = lapply(spdes, INLA::inla.spde2.precision, theta = c(log_rho, log_sigma))
  for (i in seq_along(Q)) {
    ii = which(dist_to_s0_from_mesh[[i]] == 0)
    Q[[i]] = Q[[i]][-ii, -ii]
    A[[i]] = A[[i]][, -ii]
    dist_to_s0_from_mesh[[i]] = dist_to_s0_from_mesh[[i]][-ii]
  }
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
  if (verbose) {
    message(
      "f = ", format(round(sum(res), 3), digits = 3), ", θ = ",
      paste(format(round(theta, 3), digits = 3), collapse = ", "))
  }
  res
}

# Compute the gradient of the log-likelihood using numerical derivation
loglik_grad = function(theta, sum_terms = TRUE, ...) {
  res = numDeriv::jacobian(loglik, theta, ..., sum_terms = FALSE)
  if (sum_terms) res = apply(res, 2, sum)
  res
}

# ==============================================================================
# Choose all conditioning sites for the global model
# ==============================================================================

multimesh_data = parallel::mclapply(
  X = seq_along(data$obs_index),
  mc.cores = n_cores,
  mc.preschedule = FALSE,
  FUN = function(i) {
    convex = 40
    while (TRUE) {
      args = list(
        loc = coords[c(data$obs_index[[i]], data$s0_index[[i]]), ],
        boundary = inla.nonconvex.hull(coords[data$obs_index[[i]], ], convex = convex),
        max.edge = 50)
      mesh = create_mesh(args)
      if (!is.null(mesh)) break
      convex = convex + 5
    }
    message(i)
    spde = inla.spde2.pcmatern(
      mesh,
      prior.range = c(60, .95),
      prior.sigma = c(4, .05),
      alpha = spde_alpha)
    dist_to_s0_from_mesh = dist_euclid(mesh$loc[, 1:2], data$s0[[i]])
    A = inla.spde.make.A(mesh, coords[data$obs_index[[i]], ])
    list(spde = spde, dist = dist_to_s0_from_mesh, A = A)
  })
multimesh_data = purrr::transpose(multimesh_data)

# Plot an example of the locations used for performing inference for a given
# conditioning site, together with a display of the triangular mesh
plot = local({
  i = my_index
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

# ==============================================================================
# Compute the MLE
# ==============================================================================
est = list(par = c(2, log(20), log(1), log(1.4), log(40), 0, log(20), log(1)), convergence = 1)
while (est$convergence != 0) {
  est = optim(
    par = est$par,
    fn = loglik,
    y = data$y,
    y0 = data$y0,
    dist_to_s0 = data$dist_to_s0,
    dist_to_s0_from_mesh = multimesh_data$dist,
    A = multimesh_data$A,
    verbose = TRUE,
    spdes = multimesh_data$spde,
    n_cores = 5,
    control = list(fnscale = -1, maxit = 500))
  message("done")
}

tmp = readRDS(filename)
tmp$mle = est
saveRDS(tmp, filename)

# ==============================================================================
# Fit the model to data using INLA
# ==============================================================================
est = readRDS(filename)$mle

# R-INLA requires the observations y to be on a vector format
y_inla = unlist(data$y)

# Our implementation of the cgeneric model for a requires y0 and dist_to_s0 as input.
# However, it requires one value of y0 and dist_to_s0 for each of the observations y.
# We therefore replicate y0 and dist_to_s0 in the correct way so we get to vectors with
# equal length to y_inla, where y_inla[i] was observed when y(s0) was equal to
# y0_inla[i], and the distance between them was dist_to_s0_inla[i].
dist_to_s0_inla = data$dist_to_s0[rep(seq_along(data$n), data$n)]
y0_inla = rep(unlist(data$y0), sapply(dist_to_s0_inla, length))
dist_to_s0_inla = unlist(dist_to_s0_inla)

# Define the cgeneric model for a
a_priors = list(lambda = c(3, 4), kappa = c(0, 3))
a_model = a_model(
  y0 = y0_inla,
  dist_to_s0 = dist_to_s0_inla,
  init = est$par[2:3],
  priors = a_priors)

# Define the cgeneric model for Z_b
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

inla_data = local({
  A = list()
  pb = progress_bar(length(data$n))
  for (i in seq_along(data$n)) {
    location_index = rep(data$obs_index[[i]], data$n[i])
    A[[i]] = inla.spde.make.A(
      multimesh_data$spde[[i]]$mesh,
      loc = coords[data$obs_index[[i]], ],
      index = rep(seq_along(data$obs_index[[i]]), data$n[i]),
      repl = rep(seq_len(data$n[i]), each = length(data$obs_index[[i]])))
    # Account for the constraining method
    n_spde = multimesh_data$spde[[i]]$n.spde
    ii = which(multimesh_data$dist[[i]] == 0)
    ii = (0:(data$n[i] - 1)) * n_spde + ii
    A[[i]] = A[[i]][, -ii]
    pb$tick()
  }
  pb$terminate()
  A = Matrix::bdiag(A)
  n_per_col = tail(A@p, -1) - head(A@p, -1)
  ii = which(n_per_col > 0)
  A = A[, ii]
  A = cbind(A, Diagonal(length(y_inla), 1))
  data = list(
    y = y_inla,
    spatial = c(ii, rep(NA, length(y_inla))),
    idx = c(rep(NA, length(ii)), seq_along(y_inla)))
  list(A = A, data = data)
})

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

tmp = readRDS(filename)
tmp$inla = fit
saveRDS(tmp, filename)

# # ==============================================================================
# # Fit a single-site model to the data
# # ==============================================================================
# 
# set.seed(1)
# single_site_indices = c(my_index, sample.int(length(s0_index), 5))
# 
# single_site_fits = list()
# for (i in seq_along(single_site_indices)) {
# 
#   # extract data from only one conditioning site
#   mydata = lapply(data, `[`, single_site_indices[i])
#   my_multimesh_data = lapply(multimesh_data, `[`, single_site_indices[i])
# 
#   y_inla = unlist(mydata$y)
#   dist_to_s0_inla = mydata$dist_to_s0[rep(seq_along(mydata$n), mydata$n)]
#   y0_inla = rep(unlist(mydata$y0), sapply(dist_to_s0_inla, length))
#   dist_to_s0_inla = unlist(dist_to_s0_inla)
# 
#   a_model = a_model(
#     y0 = y0_inla,
#     dist_to_s0 = dist_to_s0_inla,
#     init = est$par[2:3],
#     priors = a_priors)
# 
#   spde_model = spde_b_model_simulation(
#     n = mydata$n,
#     y0 = unlist(mydata$y0),
#     spde = my_multimesh_data$spde,
#     init = est$par[4:8],
#     priors = spde_priors,
#     dist_to_s0 = my_multimesh_data$dist)
# 
#   inla_data = local({
#     A = list()
#     pb = progress_bar(length(mydata$n))
#     for (i in seq_along(mydata$n)) {
#       location_index = rep(mydata$obs_index[[i]], mydata$n[i])
#       A[[i]] = inla.spde.make.A(
#         my_multimesh_data$spde[[i]]$mesh,
#         loc = coords[mydata$obs_index[[i]], ],
#         index = rep(seq_along(mydata$obs_index[[i]]), mydata$n[i]),
#         repl = rep(seq_len(mydata$n[i]), each = length(mydata$obs_index[[i]])))
#       # Account for the constraining method
#       n_spde = my_multimesh_data$spde[[i]]$n.spde
#       ii = which(my_multimesh_data$dist[[i]] == 0)
#       ii = (0:(mydata$n[i] - 1)) * n_spde + ii
#       A[[i]] = A[[i]][, -ii]
#       pb$tick()
#     }
#     pb$terminate()
#     A = Matrix::bdiag(A)
#     n_per_col = tail(A@p, -1) - head(A@p, -1)
#     ii = which(n_per_col > 0)
#     A = A[, ii]
#     A = cbind(A, Diagonal(length(y_inla), 1))
#     data = list(
#       y = y_inla,
#       spatial = c(ii, rep(NA, length(y_inla))),
#       idx = c(rep(NA, length(ii)), seq_along(y_inla)))
#     list(A = A, data = data)
#   })
# 
#   single_site_fits[[i]] = inla(
#     formula = formula,
#     data = inla_data$data,
#     control.predictor = list(A = inla_data$A),
#     control.family = list(hyper = list(prec = list(
#       initial = est$par[1],
#       prior = "pc.prec",
#       param = c(1, .05)))),
#     only.hyperparam = TRUE,
#     control.inla = list(
#       control.vb = list(enable = FALSE),
#       int.strategy = "eb"),
#     control.compute = list(config = TRUE),
#     num.threads = 1,
#     inla.mode = "experimental")
# 
#   single_site_fits[[i]] = inla.rerun(single_site_fits[[i]])
# 
#   # Keep only the relevant parts of the model fit, to reduce requirements on file storage
#   single_site_fits[[i]] = list(
#     misc = single_site_fits[[i]]$misc,
#     internal.marginals.hyperpar = single_site_fits[[i]]$internal.marginals.hyperpar,
#     mode = single_site_fits[[i]]$mode,
#     cpu = single_site_fits[[i]]$cpu,
#     index = single_site_indices[i],
#     s0_index = s0_index[single_site_indices[i]],
#     summary.hyperpar = single_site_fits[[i]]$summary.hyperpar)
#   single_site_fits[[i]]$misc$reordering = NULL
#   single_site_fits[[i]]$misc$family = NULL
#   single_site_fits[[i]]$misc$linkfunctions = NULL
#   class(single_site_fits[[i]]) = "inla"
# 
#   message("Progress: ", i, " / ", length(single_site_indices))
# }
# 
# tmp = readRDS(filename)
# tmp$single_site_fits = single_site_fits
# saveRDS(tmp, filename)
# 
# # # extract data from only one conditioning site
# # mydata = lapply(data, `[`, my_index)
# # my_multimesh_data = lapply(multimesh_data, `[`, my_index)
# # 
# # y_inla = unlist(mydata$y)
# # dist_to_s0_inla = mydata$dist_to_s0[rep(seq_along(mydata$n), mydata$n)]
# # y0_inla = rep(unlist(mydata$y0), sapply(dist_to_s0_inla, length))
# # dist_to_s0_inla = unlist(dist_to_s0_inla)
# # 
# # a_model = a_model(
# #   y0 = y0_inla,
# #   dist_to_s0 = dist_to_s0_inla,
# #   init = est$par[2:3],
# #   priors = a_priors)
# # 
# # spde_model = spde_b_model_simulation(
# #   n = mydata$n,
# #   y0 = unlist(mydata$y0),
# #   spde = my_multimesh_data$spde,
# #   init = est$par[4:8],
# #   priors = spde_priors,
# #   dist_to_s0 = my_multimesh_data$dist)
# # 
# # inla_data = local({
# #   A = list()
# #   pb = progress_bar(length(mydata$n))
# #   for (i in seq_along(mydata$n)) {
# #     location_index = rep(mydata$obs_index[[i]], mydata$n[i])
# #     A[[i]] = inla.spde.make.A(
# #       my_multimesh_data$spde[[i]]$mesh,
# #       loc = coords[mydata$obs_index[[i]], ],
# #       index = rep(seq_along(mydata$obs_index[[i]]), mydata$n[i]),
# #       repl = rep(seq_len(mydata$n[i]), each = length(mydata$obs_index[[i]])))
# #     # Account for the constraining method
# #     n_spde = my_multimesh_data$spde[[i]]$n.spde
# #     ii = which(my_multimesh_data$dist[[i]] == 0)
# #     ii = (0:(mydata$n[i] - 1)) * n_spde + ii
# #     A[[i]] = A[[i]][, -ii]
# #     pb$tick()
# #   }
# #   pb$terminate()
# #   A = Matrix::bdiag(A)
# #   n_per_col = tail(A@p, -1) - head(A@p, -1)
# #   ii = which(n_per_col > 0)
# #   A = A[, ii]
# #   A = cbind(A, Diagonal(length(y_inla), 1))
# #   data = list(
# #     y = y_inla,
# #     spatial = c(ii, rep(NA, length(y_inla))),
# #     idx = c(rep(NA, length(ii)), seq_along(y_inla)))
# #   list(A = A, data = data)
# # })
# # 
# # single_site_fit = inla(
# #   formula = formula,
# #   data = inla_data$data,
# #   control.predictor = list(A = inla_data$A),
# #   control.family = list(hyper = list(prec = list(
# #     initial = est$par[1],
# #     prior = "pc.prec",
# #     param = c(1, .05)))),
# #   only.hyperparam = TRUE,
# #   control.inla = list(
# #     control.vb = list(enable = FALSE),
# #     int.strategy = "eb"),
# #   control.compute = list(config = TRUE),
# #   verbose = TRUE,
# #   num.threads = 1,
# #   inla.mode = "experimental")
# # 
# # single_site_fit = inla.rerun(single_site_fit)
# # 
# # # Keep only the relevant parts of the model fit, to reduce requirements on file storage
# # single_site_fit = list(
# #   misc = single_site_fit$misc,
# #   internal.marginals.hyperpar = single_site_fit$internal.marginals.hyperpar,
# #   mode = single_site_fit$mode,
# #   cpu = single_site_fit$cpu,
# #   summary.hyperpar = single_site_fit$summary.hyperpar)
# # single_site_fit$misc$reordering = NULL
# # single_site_fit$misc$family = NULL
# # single_site_fit$misc$linkfunctions = NULL
# # class(single_site_fit) = "inla"
# # 
# # tmp = readRDS(filename)
# # tmp$single_site_inla = single_site_fit
# # saveRDS(tmp, filename)

# ==============================================================================
# Compute the C matrix
# ==============================================================================
fit = readRDS(filename)$inla
#single_site_fit = readRDS(filename)$single_site_inla
#single_site_fits = readRDS(filename)$single_site_fits

H = solve(fit$misc$cov.intern)
#H_single_site = solve(single_site_fit$misc$cov.intern)
#H_single_site = lapply(single_site_fits, function(x) solve(x$misc$cov.intern))

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
# Estimate J. We know that all time points are independent
times = unlist(data$time_index)
J = 0
for (k in seq_along(times)) {
  index = which(times == times[k])
  for (j in index) {
    J = J + grads[k, ] %*% grads[j, , drop = FALSE]
  }
}

# J_single_site = lapply(
#   X = seq_along(single_site_fits),
#   FUN = function(i) {
#     mydata = lapply(data, `[`, single_site_fits[[i]]$index)
#     my_multimesh_data = lapply(multimesh_data, `[`, single_site_fits[[i]]$index)
#     grads = loglik_grad(
#       theta = single_site_fits[[i]]$mode$theta,
#       y = mydata$y,
#       y0 = mydata$y0,
#       dist_to_s0 = mydata$dist_to_s0,
#       dist_to_s0_from_mesh = my_multimesh_data$dist,
#       A = my_multimesh_data$A,
#       n_cores = n_cores,
#       spdes = my_multimesh_data$spde,
#       sum_terms = FALSE)
#     mytimes = unlist(mydata$time_index)
#     J = 0
#     for (k in seq_along(mytimes)) {
#       index = which(mytimes == mytimes[k])
#       for (j in index) {
#         J = J + grads[k, ] %*% grads[j, , drop = FALSE]
#       }
#     }
#     J
#   })

#single_site_grads = loglik_grad(
#  theta = single_site_fit$mode$theta,
#  y = mydata$y,
#  y0 = mydata$y0,
#  dist_to_s0 = mydata$dist_to_s0,
#  dist_to_s0_from_mesh = my_multimesh_data$dist,
#  A = my_multimesh_data$A,
#  n_cores = n_cores,
#  spdes = my_multimesh_data$spde,
#  sum_terms = FALSE)
#
#mytimes = unlist(mydata$time_index)
#J_single_site = 0
#for (k in seq_along(mytimes)) {
#  index = which(mytimes == mytimes[k])
#  for (j in index) {
#    J_single_site = J_single_site + single_site_grads[k, ] %*% single_site_grads[j, , drop = FALSE]
#  }
#}

# Compute the estimate for C
C = get_C(H, J)
#C_single_site = get_C(H_single_site, J_single_site)
#C_single_site = lapply(
#  X = seq_along(single_site_fits),
#  FUN = function(i) {
#    get_C(H_single_site[[i]], J_single_site[[i]])
#  })

tmp = readRDS(filename)
tmp$C = C
#tmp$C_single_site = C_single_site
saveRDS(tmp, filename)

# ==============================================================================
# Simulate data from the unadjusted and adjusted model fits
# ==============================================================================

tmp = readRDS(filename)
C = tmp$C
#C_single_site = tmp$C_single_site
fit = tmp$inla
#single_site_fit = tmp$single_site_inla
#single_site_fits = tmp$single_site_fits
rm(tmp)

# Estimate credible intervals
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
#for (i in seq_along(single_site_fits)) {
#  name1 = paste0("unadjusted_single_", i)
#  name2 = paste0("adjusted_single_", i)
#  samples[[name1]] = inla.hyperpar.sample(
#    n = n_posterior_samples,
#    result = single_site_fits[[i]],
#    intern = TRUE)
#  samples[[name2]] = matrix(
#    data = rep(single_site_fits[[i]]$mode$theta, each = n_posterior_samples),
#    ncol = n_theta)
#  samples[[name2]] = samples[[name2]] +
#    (samples[[name1]] - samples[[name2]]) %*% t(C_single_site[[i]])
#}
#samples$unadjusted_single = inla.hyperpar.sample(
#  n = n_posterior_samples,
#  result = single_site_fit,
#  intern = TRUE)
#samples$adjusted_single = matrix(
#  data = rep(single_site_fit$mode$theta, each = n_posterior_samples),
#  ncol = n_theta)
#samples$adjusted_single = samples$adjusted_single +
#  (samples$unadjusted_single - samples$adjusted_single) %*% t(C_single_site)
for (i in seq_along(samples)) {
  samples[[i]] = exp(samples[[i]])
  samples[[i]][, 6] = samples[[i]][, 6] / (1 + samples[[i]][, 6])
}

theta_names = paste0(
  "$\\", c("tau", "lambda_a", "kappa_a", "sigma", "rho", "beta_0", "lambda_b", "kappa_b"), "$")
plot = local({
  df1 = as.data.frame(samples$unadjusted)
  names(df1) = theta_names
  df1$adjusted = "Unadjusted"
  #df1$composite = "Composite"
  df2 = as.data.frame(samples$adjusted)
  names(df2) = theta_names
  df2$adjusted = "Adjusted"
  #df2$composite = "Composite"
  #df3 = as.data.frame(samples$unadjusted_single)
  #names(df3) = theta_names
  #df3$adjusted = "Unadjusted"
  #df3$composite = "Single-site"
  #df4 = as.data.frame(samples$adjusted_single)
  #names(df4) = theta_names
  #df4$adjusted = "Adjusted"
  #df4$composite = "Single-site"
  rbind(df1, df2) |>
    #tidyr::pivot_longer(-c(adjusted, composite)) |>
    tidyr::pivot_longer(-adjusted) |>
    dplyr::mutate(
      name = factor(name, levels = theta_names),
      adjusted = factor(adjusted, levels = c("Unadjusted", "Adjusted"))) |>
    #dplyr::filter(
    #  !(name == theta_names[1] & (value < 15 | value > 28)),
    #  !(name == theta_names[2] & (value < 10 | value > 17.5)),
    #  !(name == theta_names[3] & (value < 1 | value > 1.5)),
    #  !(name == theta_names[4] & (value < 1.3 | value > 1.9)),
    #  !(name == theta_names[5] & (value < 58 | value > 82)),
    #  !(name == theta_names[6] & (value < .4 | value > 1)),
    #  !(name == theta_names[7] & (value < 5 | value > 14)),
    #  !(name == theta_names[8] & (value < .3 | value > 1.2))) |>
    ggplot() +
    #geom_density(aes(x = value, linetype = adjusted, col = composite), trim = TRUE) +
    geom_density(aes(x = value, linetype = adjusted)) +
    facet_wrap(~name, scales = "free", nrow = 2) +
    labs(linetype = "", col = "", x = "Parameter value", y = "Posterior density") +
    theme_light() +
    #scale_color_manual(values = c("black", "darkgray")) +
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

simulate_data = function(samples, n_samples, n_samples_per_theta) {
  n_samples = min(n_samples, nrow(samples))
  my_s0_index = s0_index[my_index]
  my_coords = extract_thinned_out_circles(
    coords = coords,
    center = coords[my_s0_index, ],
    n = c(1, 2),
    r = c(6, Inf))
  mesh = inla.mesh.2d(
    loc = my_coords,
    boundary = list(
      inla.nonconvex.hull(my_coords, convex = -.2),
      inla.nonconvex.hull(my_coords, convex = -1)),
    max.edge = c(10, 50))
  A = inla.spde.make.A(mesh, my_coords)
  dist_to_s0 = dist_euclid(coords[my_s0_index, ], my_coords)
  dist_to_s0_from_mesh = dist_euclid(coords[my_s0_index, ], mesh$loc[, -3])
  y0 = rexp(n_samples * n_samples_per_theta) + threshold
  ii = which(dist_to_s0_from_mesh == 0)
  jj = which(dist_to_s0 == 0)
  dist_to_s0_from_mesh = dist_to_s0_from_mesh[-ii]
  dist_to_s0 = dist_to_s0[-jj]
  A = A[-jj, -ii]
  spde = inla.spde2.pcmatern(
    mesh = mesh,
    alpha = spde_alpha,
    prior.range = c(exp(fit$mode$theta[5]), .5),
    prior.sigma = c(exp(fit$mode$theta[4]), .5))
  data = list(
    y0 = list(y0),
    dist_to_s0 = list(dist_to_s0),
    n = n_samples * n_samples_per_theta)
  data$y = pbapply::pblapply(
    X = 1:n_samples,
    cl = n_cores,
    FUN = function(i) {
      params = samples[i, ]
      tau = params[1]
      lambda_a = params[2]
      kappa_a = params[3]
      sigma = params[4]
      rho = params[5]
      beta0 = params[6]
      lambda_b = params[7]
      kappa_b = params[8]
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
      Q = inla.spde2.precision(spde, log(c(rho, sigma)))
      Q = Q[-ii, -ii]
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

n_samples = 1e3
n_samples_per_theta = 100

mydata = simulate_data(samples$unadjusted, n_samples, n_samples_per_theta)
unadjusted_chi = empirical_chi(mydata, attr(chi, "thresholds"), attr(chi, "unique_dist"))
unadjusted_moments = empirical_moments(
  data = mydata,
  unique_dist = attributes(moments)$unique_dist,
  unique_y0 = attributes(moments)$unique_y0)

mydata = simulate_data(samples$adjusted, n_samples, n_samples_per_theta)
adjusted_chi = empirical_chi(mydata, attr(chi, "thresholds"), attr(chi, "unique_dist"))
adjusted_moments = empirical_moments(
  data = mydata,
  unique_dist = attributes(moments)$unique_dist,
  unique_y0 = attributes(moments)$unique_y0)

rm(mydata)

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
# Compute out-of-sample log-likelihood to test if the adjustment method helps
# ==============================================================================

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

data = extract_extreme_fields(
  data = transformed_data,
  coords = coords,
  s0_index = s0_index,
  threshold = threshold,
  n = thinning,
  r = cumsum(4 * thinning))

n_posterior_samples = 1e3
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

tmp = readRDS(filename)
tmp$loglik = ll
saveRDS(tmp, filename)

tmp = purrr::transpose(ll)
tmp = sapply(tmp, function(x) {
  x = do.call(rbind, x)
  apply(x, 2, log_mean_exp)
  })

sum(tmp[, 2]) - sum(tmp[, 1])

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

summary(bootstrap_sums[2, ] - bootstrap_sums[1, ])
