devtools::load_all()
library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(INLA)
library(inlabru)

# Compile and link all the cgeneric scripts, if this has not already been done
make_cgeneric("all")

# ==============================================================================
# Decide necessary model parameters:
# ==============================================================================
threshold = qlaplace(.9997) # The threshold t for defining the conditional extremes model
n_cores = 15 # Run code in parallel. Check yourself that this is an ok number
r = 5 # Radius used for computing aggregated empirical distribution functions
n_posterior_samples = 1e5 # Number of samples from posterior for estimating credible interval
thinning = c(1, 2, 4, 8, 16)

filename = file.path(results_dir(), "case-study.rds")
if (!file.exists(filename)) saveRDS(list(), filename)

# ==============================================================================
# Load the radar data
# ==============================================================================
radar = readRDS(file.path(downloads_dir(), "radar.rds"))
coords = st_coordinates(radar$coords)
n_loc = nrow(coords)

test_index = which(lubridate::year(radar$times) >= 2020)

# ==============================================================================
# Transform the data to Laplace margins
# ==============================================================================
cl = parallel::makeForkCluster(n_cores)
transformed_data = pbapply::pblapply(
  X = 1:n_loc,
  cl = cl,
  FUN = function(i) {
    F = aggregated_ecdf(radar$data, coords, coords[i, ], radius = r, zero_threshold = .1)
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
# Examine patterns in the data
# ==============================================================================

first_threshold = qlaplace(.99)
delta_s0 = 1
s0_index = get_s0_index(coords, delta_s0)

first_data = extract_extreme_fields(
  data = transformed_data,
  coords = coords,
  s0_index = s0_index,
  threshold = first_threshold,
  n_cores = n_cores,
  n = 1,
  r = Inf)

dist_table = first_data$dist_to_s0 |> unlist() |> round() |> table()
unique_dist = dist_table[which(dist_table > 800)] |> names() |> as.numeric()
y0_table = first_data$y0 |> unlist() |> round(1) |> table()
unique_y0 = y0_table[which(y0_table > 40)] |> names() |> as.numeric()

moments = empirical_moments(first_data, unique_dist, unique_y0)

moment_plots = local({
  nan_index = is.nan(moments$mean[, 1])
  unique_dist = attr(moments, "unique_dist")[!nan_index]
  unique_y0 = attr(moments, "unique_y0")
  df = data.frame(
        sd = as.numeric(rbind(0, moments$sd[!nan_index, ])),
        mean = as.numeric(rbind(unique_y0, moments$mean[!nan_index, ])),
        dist = rep(c(0, unique_dist), length(unique_y0)),
        y0 = rep(unique_y0, each = length(unique_dist) + 1))
  sd_plot = ggplot(df) +
    geom_line(aes(x = dist, y = sd, group = y0, col = y0), size = 1) +
    scale_color_gradient(low = "black", high = "lightgray") +
    labs(y = "$\\widehat \\zeta(d; y_0)$", col = "$y_0$", x = "$d$") +
    theme_light() +
    theme(axis.title.y = element_text(angle = 0, vjust = .5))
  mean_plot = ggplot(df) +
    geom_line(aes(x = dist, y = mean, group = y0, col = y0), size = 1) +
    labs(y = "$\\widehat \\mu(d; y_0)$", col = "$y_0$", x = "$d$") +
    scale_color_gradient(low = "black", high = "lightgray") +
    theme_light() +
    theme(axis.title.y = element_text(angle = 0, vjust = .5),
          strip.text = element_text(colour = "black"),
          strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"))
  list(mean = mean_plot, sd = sd_plot)
})

plot = patchwork::wrap_plots(moment_plots$mean, moment_plots$sd, nrow = 1)
tikz_plot(file.path(image_dir(), "case-study_bad_properties.pdf"), plot, width = 10, height = 5)

# ==============================================================================
# Examine properties with a larger threshold
# ==============================================================================

delta_s0 = 1
s0_index = get_s0_index(coords, delta_s0)
data = extract_extreme_fields(
  data = transformed_data,
  coords = coords,
  s0_index = s0_index,
  threshold = threshold,
  n_cores = n_cores,
  n = 1,
  r = Inf)

p = 1 - 10^seq(-3, -5, by = -.1)
thresholds = qlaplace(p)
thresholds = thresholds[thresholds > threshold]

dist_table = data$dist_to_s0 |> unlist() |> round() |> table()
unique_dist = dist_table[which(dist_table > 800)] |> names() |> as.numeric()
y0_table = data$y0 |> unlist() |> round(1) |> table()
unique_y0 = y0_table[which(y0_table > 20)] |> names() |> as.numeric()
n_over_thresholds = sapply(thresholds, \(x) sum(unlist(data$y0) > x))
thresholds = thresholds[n_over_thresholds > 100]

chi = empirical_chi(data, thresholds, unique_dist)
moments = empirical_moments(data, unique_dist, unique_y0)

#plot_chi(chi)
#plot_moments(moments)

# ==============================================================================
# Compute the MLE
# ==============================================================================

delta_s0 = 2
s0_index = get_s0_index(coords, delta_s0)

data = extract_extreme_fields(
  data = transformed_data[-test_index, ],
  coords = coords,
  s0_index = s0_index,
  threshold = threshold,
  n_cores = n_cores,
  n = thinning,
  r = cumsum(6 * thinning))

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
  s0 = exp(theta[6])
  lambda_b = exp(theta[7])
  kappa_b = exp(theta[8])
  Q = lapply(spdes, INLA::inla.spde2.precision, theta = c(log_rho, log_sigma))
  if (length(Q) == 1) Q = rep(Q, length(dist_to_s0_from_mesh))
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
    sd = 1 + s0 * exp(-(dist / lambda_b)^kappa_b)
    matrix(rep(sd, length(y)), nrow = length(dist), ncol = length(y))
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

# Compute the Hessian of the log-likelihood using numerical derivation
loglik_hess = function(theta, ...) {
  res = numDeriv::hessian(loglik, theta, ...)
  res
}

# Choose all conditioning sites for computing the composite MLE
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
      prior.sigma = c(4, .05))
    dist_to_s0_from_mesh = dist_euclid(mesh$loc[, 1:2], data$s0[[i]])
    A = inla.spde.make.A(mesh, coords[data$obs_index[[i]], ])
    list(spde = spde, dist = dist_to_s0_from_mesh, A = A)
  })
multimesh_data = purrr::transpose(multimesh_data)

# Compute the MLE
# ------------------------------------------------------------------------------

est = list(par = c(0, log(30), log(.5), log(1.9), log(30), log(1), log(20), log(1)), convergence = 1)
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

delta_s0 = 1
s0_index = get_s0_index(coords, delta_s0)

data = extract_extreme_fields(
  data = transformed_data[-test_index, ],
  coords = coords,
  s0_index = s0_index,
  threshold = threshold,
  n_cores = n_cores,
  n = thinning,
  r = cumsum(6 * thinning))

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
      prior.sigma = c(4, .05))
    dist_to_s0_from_mesh = dist_euclid(mesh$loc[, 1:2], data$s0[[i]])
    A = inla.spde.make.A(mesh, coords[data$obs_index[[i]], ])
    list(spde = spde, dist = dist_to_s0_from_mesh, A = A)
  })
multimesh_data = purrr::transpose(multimesh_data)

est = readRDS(filename)$mle

# R-INLA requires the observations y to be on a vector format
y_inla = unlist(data$y)
# Keep track of all indices where y_inla is NA, so we can remove these
# before feeding the data into R-INLA, for reduced requirements on CPU and RAM.
na_index = which(is.na(y_inla))

# Our implementation of the cgeneric model for a requires y0 and dist_to_s0 as input.
# However, it requires one value of y0 and dist_to_s0 for each of the observations y.
# We therefore replicate y0 and dist_to_s0 in the correct way so we get to vectors with
# equal length to y_inla, where y_inla[i] was observed when y(s0) was equal to
# y0_inla[i], and the distance between them was dist_to_s0_inla[i].
dist_to_s0_inla = data$dist_to_s0[rep(seq_along(data$n), data$n)]
y0_inla = rep(unlist(data$y0), sapply(dist_to_s0_inla, length))
dist_to_s0_inla = unlist(dist_to_s0_inla)

# Remove data at all fo the NA indices
y_inla = y_inla[-na_index]
y0_inla = y0_inla[-na_index]
dist_to_s0_inla = dist_to_s0_inla[-na_index]

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
  b0 = c(0, 4),
  lambda = c(3, 4),
  kappa = c(0, 3))

spde_model = spde_b_model_case_study(
  n = data$n,
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
  A = A[-na_index, ii]
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

fit = inla(
  formula = formula,
  data = inla_data$data,
  control.predictor = list(A = inla_data$A),
  control.family = list(hyper = list(prec = list(
    initial = est$par[1],
    prior = "pc.prec",
    param = c(1, .05)))),
  only.hyperparam = TRUE,
  verbose = TRUE,
  control.inla = list(restart = 1), # INLA struggles a bit, so we restart the optimiser once
  num.threads = 1,
  inla.mode = "experimental")

fit = inla.rerun(fit) # We also rerun INLA once more, to improve the fit

# Keep only the relevant parts of the model fit, to reduce requirements on file storage
fit = list(
  misc = fit$misc,
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

# ==============================================================================
# Compute the C matrix
# ==============================================================================
fit = readRDS(filename)$inla

H = solve(fit$misc$cov.intern)

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

# Estimate J with a sliding window.
times = unlist(data$time_index)
J = 0
for (k in seq_along(times)) {
  index = which(abs(times - times[k]) <= 5)
  for (j in index) {
    J = J + grads[k, ] %*% grads[j, , drop = FALSE]
  }
}

# Compute the estimate for C
C = get_C(H, J)

tmp = readRDS(filename)
tmp$C = C
saveRDS(tmp, filename)

# ==============================================================================
# Simulate data from the unadjusted and adjusted model fits
# ==============================================================================

tmp = readRDS(filename)
C = tmp$C
fit = tmp$inla
rm(tmp)

# Estimate credible intervals
n_theta = length(fit$mode$theta)
samples = list()
set.seed(1)
samples$unadjusted = inla.hyperpar.sample(n_posterior_samples, fit, intern = TRUE)
samples$adjusted = matrix(
  data = rep(fit$mode$theta, each = n_posterior_samples),
  ncol = n_theta)
samples$adjusted = samples$adjusted + (samples$unadjusted - samples$adjusted) %*% t(C)
for (i in seq_along(samples)) {
  samples[[i]] = exp(samples[[i]])
}

theta_names = paste0(
  "$\\", c("tau", "lambda_a", "kappa_a", "sigma", "rho", "b_0", "lambda_b", "kappa_b"), "$")
theta_names[6] = "$b_0$"
plot = local({
  df1 = as.data.frame(samples$unadjusted)
  names(df1) = theta_names
  df1$tag = "Unadjusted"
  df2 = as.data.frame(samples$adjusted)
  names(df2) = theta_names
  df2$tag = "Adjusted"
  rbind(df1, df2) |>
    tidyr::pivot_longer(-tag) |>
    dplyr::mutate(name = factor(name, levels = theta_names)) |>
    dplyr::filter(
      !(name == theta_names[1] & (value < 0 | value > 150)),
      !(name == theta_names[2] & (value < 0 | value > 12)),
      !(name == theta_names[3] & (value < 0 | value > 2.5)),
      !(name == theta_names[4] & (value < 1.5 | value > 3)),
      !(name == theta_names[5] & (value < 8 | value > 21)),
      !(name == theta_names[6] & (value < 1 | value > 5)),
      !(name == theta_names[7] & (value < 2.5 | value > 7)),
      !(name == theta_names[8] & (value < .3 | value > 1.5))) |>
    ggplot() +
    geom_density(aes(x = value, linetype = tag), trim = TRUE) +
    facet_wrap(~name, scales = "free", nrow = 2) +
    labs(linetype = "", x = "Parameter value", y = "Posterior density") +
    theme_light() +
    theme(
      text = element_text(size = 15),
      strip.text = element_text(size = 15, colour = "black"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      legend.position = "top")
})

# Export the plot to pdf format
tikz_plot(file.path(image_dir(), "case-study_posteriors.pdf"),
          plot, width = 10, height = 6)

simulate_data = function(samples, n_samples, n_samples_per_theta) {
  n_samples = min(n_samples, nrow(samples))
  my_s0_index = s0_index[25]
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
      s0 = params[6]
      lambda_b = params[7]
      kappa_b = params[8]
      a_func = function(y, dist) {
        alpha = exp(- (dist / lambda_a)^kappa_a)
        a = rep(y, each = length(alpha)) * rep(alpha, length(y))
        matrix(a, nrow = length(dist), ncol = length(y))
      }
      b_func = function(y, dist) {
        sd = 1 + s0 * exp(-(dist / lambda_b)^kappa_b)
        matrix(rep(sd, length(y)), nrow = length(dist), ncol = length(y))
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
  nan_index = is.nan(adjusted_chi[, 1])
  thresholds = attr(chi, "thresholds")
  p = plaplace(thresholds)
  unique_dist = attr(chi, "unique_dist")[!nan_index]
  tmp = list(Empirical = chi, Adjusted = adjusted_chi, Unadjusted = unadjusted_chi)
  for (i in seq_along(tmp)) {
    tmp[[i]] = data.frame(
      chi = as.numeric(rbind(1, tmp[[i]][!nan_index, ])),
      dist = rep(c(0, unique_dist), length(p)),
      p = rep(-log10(1 - p), each = length(unique_dist) + 1),
      tag = names(tmp)[i])
  }
  do.call(rbind, tmp) |>
    dplyr::mutate(tag = factor(tag, levels = c("Empirical", "Adjusted", "Unadjusted"))) |>
    ggplot() +
    geom_line(aes(x = dist, y = chi, group = p, col = p)) +
    labs(col = "$p$", y = "$\\chi_p(d)$", x = "$d$") +
    scale_color_gradient(
      breaks = c(3.5, 4, 4.5),
      limits = c(3.4, 4.6),
      labels = paste0("$1 - 10^{-", c(3.5, 4, 4.5), "}$"),
      low = "black",
      high = "lightgray") +
    theme_light() +
    theme(axis.title.y = element_text(angle = 0, vjust = .5),
          strip.text = element_text(colour = "black"),
          strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")) +
    facet_wrap(~tag, nrow = 1) +
    lims(x = c(0, 20))
})

moment_plots = local({
  nan_index = is.nan(adjusted_moments$mean[, 1])
  unique_dist = attr(moments, "unique_dist")[!nan_index]
  unique_y0 = attr(moments, "unique_y0")
  tmp = list(Empirical = moments, Adjusted = adjusted_moments, Unadjusted = unadjusted_moments)
  for (i in seq_along(tmp)) {
    tmp[[i]] = data.frame(
        sd = as.numeric(rbind(0, tmp[[i]]$sd[!nan_index, ])),
        mean = as.numeric(rbind(unique_y0, tmp[[i]]$mean[!nan_index, ])),
        dist = rep(c(0, unique_dist), length(unique_y0)),
        y0 = rep(unique_y0, each = length(unique_dist) + 1),
        tag = names(tmp)[i])
  }
  sd_plot = do.call(rbind, tmp) |>
    dplyr::mutate(tag = factor(tag, levels = c("Empirical", "Adjusted", "Unadjusted"))) |>
    ggplot() +
    geom_line(aes(x = dist, y = sd, group = y0, col = y0)) +
    scale_color_gradient(low = "black", high = "lightgray") +
    labs(y = "$\\widehat \\zeta(d; y_0)$", col = "$y_0$", x = "$d$") +
    theme_light() +
    theme(axis.title.y = element_text(angle = 0, vjust = .5),
          strip.text = element_text(colour = "black"),
          strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")) +
    facet_wrap(~tag, nrow = 1)
  mean_plot = do.call(rbind, tmp) |>
    dplyr::mutate(tag = factor(tag, levels = c("Empirical", "Adjusted", "Unadjusted"))) |>
    ggplot() +
    geom_line(aes(x = dist, y = mean, group = y0, col = y0)) +
    labs(y = "$\\widehat \\mu(d; y_0)$", col = "$y_0$", x = "$d$") +
    scale_color_gradient(low = "black", high = "lightgray") +
    theme_light() +
    theme(axis.title.y = element_text(angle = 0, vjust = .5),
          strip.text = element_text(colour = "black"),
          strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")) +
    facet_wrap(~tag, nrow = 1)
  list(mean = mean_plot, sd = sd_plot)
})

plot = patchwork::wrap_plots(chi_plot, moment_plots$mean, moment_plots$sd, ncol = 1)
plot = plot * theme(text = element_text(size = 16))
tikz_plot(file.path(image_dir(), "case-study_properties.pdf"), plot, width = 11, height = 8)

# ==============================================================================
# Compute out-of-sample composite log-likelihoods
# ==============================================================================

# Create the test data, based on all conditioning sites
# that were not used for performing inference with the global models
eval_data = extract_extreme_fields(
  data = transformed_data[test_index, ],
  coords = coords,
  s0_index = s0_index,
  threshold = threshold,
  n = thinning,
  r = cumsum(10 * thinning))

multimesh_eval_data = parallel::mclapply(
  X = seq_along(eval_data$obs_index),
  mc.cores = n_cores,
  mc.preschedule = FALSE,
  FUN = function(i) {
    convex = 40
    while (TRUE) {
      args = list(
        loc = coords[c(eval_data$obs_index[[i]], eval_data$s0_index[[i]]), ],
        boundary = inla.nonconvex.hull(coords[eval_data$obs_index[[i]], ], convex = convex),
        max.edge = 50)
      mesh = create_mesh(args)
      if (!is.null(mesh)) break
      convex = convex + 5
    }
    message(i)
    spde = inla.spde2.pcmatern(
      mesh,
      prior.range = c(60, .95),
      prior.sigma = c(4, .05))
    dist_to_s0_from_mesh = dist_euclid(mesh$loc[, 1:2], eval_data$s0[[i]])
    A = inla.spde.make.A(mesh, coords[eval_data$obs_index[[i]], ])
    list(spde = spde, dist = dist_to_s0_from_mesh, A = A)
  })
multimesh_eval_data = purrr::transpose(multimesh_eval_data)

n_samples = 1e3
ll = pbapply::pblapply(
  X = 1:n_samples,
  cl = n_cores,
  FUN = function(i) {
    ll = list()
    for (name in names(samples)) {
      param = samples[[name]][i, ]
      param = log(param)
      ll[[name]] = loglik(
        theta = param,
        y = eval_data$y,
        y0 = eval_data$y0,
        dist_to_s0 = eval_data$dist_to_s0,
        dist_to_s0_from_mesh = multimesh_eval_data$dist,
        A = multimesh_eval_data$A,
        sum_terms = FALSE,
        spdes = multimesh_eval_data$spde)
    }
    ll
  })

tmp = readRDS(filename)
tmp$loglik = ll
saveRDS(tmp, filename)
rm(tmp)

ll = readRDS(filename)$loglik

tmp = purrr::transpose(ll)
tmp = sapply(tmp, function(x) {
  x = do.call(rbind, x)
  apply(x, 2, log_mean_exp)
  })

apply(tmp, 2, sum)
sum(tmp[, 2]) - sum(tmp[, 1])

set.seed(1)
n_bootstraps = 5000
pb = progress_bar(n_bootstraps)
bootstrap_diffs = sapply(
  X = seq_len(n_bootstraps),
  FUN = function(i) {
    index = sample.int(length(test_index), length(test_index), replace = TRUE)
    time_index = unlist(eval_data$time_index)
    log_scores = lapply(
      X = index,
      FUN = function(j) {
        ii = which(time_index == j)
        if (length(ii) > 0) {
          tmp[ii, ]
        } else {
          NULL
        }
      })
    log_scores = do.call(rbind, log_scores)
    sums = apply(log_scores, 2, sum)
    pb$tick()
    sums[2] - sums[1]
  })
pb$terminate()

summary(bootstrap_diffs)
mean(bootstrap_diffs < 0)
