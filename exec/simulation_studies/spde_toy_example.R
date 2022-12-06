devtools::load_all()
library(INLA)
library(mvtnorm)
library(ggplot2)

domain_size = 25 # Spatial domain is quadratic with lengths equal to domain_size
n_loc = 400 # Number of locations containing data
n_train = 200 # Number of replications of the Gaussian field in the training data
n_theta_star = 1e5 # Number of replications of the Gaussian field for estimating θ*
n_cores = 7 # Run code in parallel
n_repl = 300 # Number of replications of the experiment
n_posterior_samples = 1e5 # Number of samples from posterior for estimating credible interval
probs = c(.9, .95, .99) # Credible interval probabilities
filename = file.path(results_dir(), "low-rank.rds") # Filename for results
sigma = 1 # Standard deviance of the Gaussian random field
rho = 12 # Range of the Gaussian random field
tau = 100 # Precision of the nugget effect

# Draw all the locations randomly
set.seed(1)
loc = matrix(runif(n_loc * 2) * domain_size, ncol = 2)
dist = as.matrix(dist(loc))

# Matérn correlation function
matern_corr = function(dist, rho, nu = 1.5) {
  kappa = sqrt(8 * nu) / rho
  res = 2 ^ (1 - nu) / gamma(nu) * (kappa * dist) ^ nu * besselK(kappa * dist, nu)
  res[dist == 0] = 1
  res
}

# Compute covariance matrix
cov_mat = matern_corr(dist, rho) * sigma^2

# Create the coarse mesh
mesh = inla.mesh.2d(
  loc.domain = loc,
  boundary = list(inla.nonconvex.hull(loc, convex = -.3),
                  inla.nonconvex.hull(loc, convex = -1)),
  cutoff = 4,
  max.edge = c(5, 50))
#plot(mesh); points(loc)

# Define the SPDE and the projection matrix for the mesh
spde = inla.spde2.pcmatern(
  mesh,
  prior.range = c(rho, 0.5),
  prior.sigma = c(sigma, 0.5))
A = inla.spde.make.A(mesh, loc)

# Log-likelihood for the low-rank SPDE approximation
loglik = function(theta,
                  y,
                  spde,
                  A,
                  return_vec = FALSE) {
  tau = exp(theta[1])
  rho = exp(theta[2])
  sigma = exp(theta[3])
  Q = INLA::inla.spde2.precision(spde, theta = c(log(rho), log(sigma)))
  cov_mat = A %*% Matrix::solve(Q, Matrix::t(A))
  cov_mat = as.matrix(cov_mat) + diag(1 / tau, nrow(cov_mat))
  res = mvtnorm::dmvnorm(t(y), sigma = cov_mat, log = TRUE)
  if (!return_vec) res = sum(res)
  res
}

# Sample data for estimating θ*
set.seed(1)
y = mvtnorm::rmvnorm(n_theta_star, sigma = cov_mat) +
  rnorm(n_theta_star * n_loc, sd = tau^-.5)
y = t(y)

# Estimate θ*
est = list(par = c(log(tau), log(rho), log(sigma)), convergence = 1)
while (est$convergence != 0) {
  est = optim(
    par = est$par,
    fn = loglik,
    y = y,
    spde = spde,
    A = A,
    control = list(fnscale = -1, trace = 6))
}
theta_star = est$par

# Run all the n_repl experiments
set.seed(1)
model_fits = list()
attributes(model_fits)$theta_star = theta_star
pb = progress_bar(n_repl)
for (i in 1:n_repl) {

  # Simulate data from the true distribution
  y = mvtnorm::rmvnorm(n_train, sigma = cov_mat) + rnorm(n_train * n_loc, sd = tau^-.5)
  y = t(y)

  # Create the necessary objects for running R-INLA
  effects = list(spatial = inla.spde.make.index("spatial", n.spde = mesh$n, n.repl = n_train))
  y_inla = as.numeric(y)
  A_inla = inla.spde.make.A(
    mesh,
    loc,
    index = rep(1:n_loc, n_train),
    repl = rep(1:n_train, each = n_loc))
  stack = inla.stack(
    data = list(y = y_inla),
    A = list(spatial = A_inla),
    effects = effects)
  formula = y ~ -1 + f(spatial, model = spde, replicate = spatial.repl, nrep = n_train)

  # Compute the MLE of the data with the misspecified model,
  # and use this as initial values for R-INLA
  init = optim(
    par = c(log(tau), log(rho), log(sigma)),
    fn = loglik,
    y = y,
    spde = spde,
    A = A,
    control = list(fnscale = -1))

  # Fit the model with R-INLA
  fit = inla(
    formula = formula,
    data = inla.stack.data(stack),
    family = "gaussian",
    num.threads = n_cores,
    only.hyperparam = TRUE,
    control.mode = list(theta = init$par, restart = TRUE),
    control.inla = list(control.vb = list(enable = FALSE)),
    inla.mode = "experimental",
    control.predictor = list(A = inla.stack.A(stack)))

  # Estimate θ* using only the training data
  theta = fit$mode$theta

  # Estimate H, J and C
  H = solve(fit$misc$cov.intern)
  grads = numDeriv::jacobian(
    function(x) loglik(x, y, spde, A, return_vec = TRUE),
    x = theta)
  J = 0
  for (j in 1:n_train) {
    J = J + grads[j, ] %*% grads[j, , drop = FALSE]
  }
  C = get_C(H, J)

  # Keep only the relevant parts of the model fit, to reduce requirements on file storage
  fit = list(
    misc = fit$misc,
    internal.marginals.hyperpar = fit$internal.marginals.hyperpar,
    summary.hyperpar = fit$summary.hyperpar,
    model = list(theta = fit$mode$theta))
  fit$misc$reordering = NULL
  fit$misc$family = NULL
  fit$misc$linkfunctions = NULL
  class(fit) = "inla"

  # Save the temporary results
  model_fits[[i]] = list(fit = fit, C = C)
  saveRDS(model_fits, filename)
  pb$tick()
}
pb$terminate()

# Load the results
res = readRDS(filename)
theta_star = attributes(res)$theta_star
n_theta = length(theta_star)
theta_names = c("log_precision", "log_rho", "log_sigma")

# Estimate credible intervals
pb = progress_bar(length(res))
intervals = list()
for (i in seq_along(res)) {
  intervals[[i]] = list()
  theta_uncorrected = inla.hyperpar.sample(n_posterior_samples, res[[i]]$fit, intern = TRUE)
  for (k in 1:2) {
    if (k == 1) {
      theta_corrected = theta_uncorrected
      label = "Unadjusted"
    } else {
      theta_corrected = matrix(
        data = rep(res[[i]]$fit$mode$theta, each = n_posterior_samples),
        ncol = n_theta)
      theta_corrected = theta_corrected + (theta_uncorrected - theta_corrected) %*% t(res[[i]]$C)
      label = "Adjusted"
    }
    intervals[[i]][[length(intervals[[i]]) + 1]] = rbind(
      as.data.frame(apply(theta_corrected, 2, quantile, probs = (1 - probs) / 2)) |>
        dplyr::mutate(label = label, prob = probs, lim = "lower"),
      as.data.frame(apply(theta_corrected, 2, quantile, probs = 1 - (1 - probs) / 2)) |>
        dplyr::mutate(label = label, prob = probs, lim = "upper"))
    row.names(intervals[[i]][[length(intervals[[i]])]]) = NULL
    names(intervals[[i]][[length(intervals[[i]])]])[1:n_theta] = theta_names
  }
  intervals[[i]] = do.call(rbind, intervals[[i]]) |>
    dplyr::mutate(iter = i)
  pb$tick()
}
intervals = do.call(rbind, intervals)
pb$terminate()

# Find if θ* is inside the credible intervals or not
is_inside_interval = local({
  tmp1 = dplyr::filter(intervals, lim == "lower")
  tmp2 = dplyr::filter(intervals, lim == "upper")
  for (i in seq_along(theta_names)) {
    tmp1[[i]] = ifelse(tmp1[[i]] < theta_star[i], TRUE, FALSE)
    tmp2[[i]] = ifelse(tmp2[[i]] > theta_star[i], TRUE, FALSE)
  }
  tmp = (tmp1[, 1:n_theta] & tmp2[, 1:n_theta]) |>
    as.data.frame() |>
    dplyr::mutate(label = tmp1$label, prob = tmp1$prob, iter = tmp1$iter)
  tmp
})

# Compute coverage percentages
tmp = is_inside_interval |>
  tidyr::pivot_longer(all_of(theta_names)) |>
  dplyr::group_by(prob, label, name) |>
  dplyr::summarise(coverage = mean(value)) |>
  dplyr::mutate(
    name = factor(
      x = name,
      levels = theta_names,
      labels = paste0("$\\", c("tau", "rho", "sigma"))),
    label = factor(
      x = label,
      levels = c("Unadjusted", "Adjusted"),
      labels = c("$", "_\\text{adj}$"))) |>
  tidyr::pivot_wider(names_from = c(name, label), values_from = coverage)
tmp = tmp[, c(1, 5, 2, 6, 3, 7, 4)]
names(tmp) = sub("_", "", names(tmp))
print(tmp)

# ==============================================
# Reformat the results into latex tabular format
# ==============================================
table = paste(paste(c("Aim", names(tmp)[-1]), collapse = " & "), "\\\\")
table[2] = "\\hline"
j = 3
for (i in 1:nrow(tmp)) {
  table[j] = paste(
    c(paste0("$", round(100 * tmp$prob[i], digits = 0), "\\%$"),
      paste0("$", round(100 * tmp[i, -1], digits = 0), "\\%$")),
    collapse = " & ")
  table[j] = paste(table[j], "\\\\")
  j = j + 1
}
table = paste(paste(table, collapse = "\n"), "\n")

cat(table)
