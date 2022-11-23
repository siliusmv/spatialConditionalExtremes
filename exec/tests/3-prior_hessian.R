devtools::load_all()
library(numDeriv)

# Derivatives for the Gaussian prior part
# ------------------------------------------------------------------------------
log_dnorm_hess = function(sigma) -1 / sigma^2

mu = rnorm(1)
sigma = rexp(1)
x = rnorm(1, mu, sigma)

v1 = log_dnorm_hess(sigma)
v2 = numDeriv::hessian(function(x) log(dnorm(x, mu, sigma)), x)
v1 - v2

# Derivatives for the log-gamma prior part
# ------------------------------------------------------------------------------
dloggamma = function(x, a, b) {
  b^a / gamma(a) * exp(a * x - b * exp(x))
}

a = rexp(1)
b = rexp(1)
x = rexp(1)

integrate(dloggamma, -50, 50, a = a, b = b)

log_dloggamma_hess = function(x, b) -b * exp(x)

v1 = log_dloggamma_hess(x, b)
v2 = numDeriv::hessian(function(x) log(dloggamma(x, a, b)), x)
v1 - v2

# Derivatives for the PC prior part
# ------------------------------------------------------------------------------
dpc = function(theta, rho_0, sigma_0, p_rho, p_sigma) {
  if (is.matrix(theta)) {
    rho = exp(theta[, 1])
    sigma = exp(theta[, 2])
  } else {
    rho = exp(theta[1])
    sigma = exp(theta[2])
  }
  lambda_1 = -log(p_rho) * rho_0
  lambda_2 = -log(p_sigma) / sigma_0
  lambda_1 * lambda_2 * rho^-1 * sigma * exp(-lambda_1 / rho - lambda_2 * sigma)
}

rho_0 = rexp(1)
sigma_0 = rexp(1)
p_rho = runif(1)
p_sigma = runif(1)

integrate(
  f = function(x, ...) {
    sapply(
      X = x,
      FUN = function(x, ...) {
        integrate(
          f = function(y, ...) dpc(cbind(x, y), ...),
          lower = -50,
          upper = 50,
          ...)$value
      },
      ...)
  },
  lower = -50,
  upper = 50,
  rho_0 = rho_0,
  sigma_0 = sigma_0,
  p_rho = p_rho,
  p_sigma = p_sigma)

log_dpc_hess = function(theta, rho_0, sigma_0, p_rho, p_sigma) {
  lambda_1 = -log(p_rho) * rho_0
  lambda_2 = -log(p_sigma) / sigma_0
  - matrix(c(lambda_1 * exp(-theta[1]), 0, 0, lambda_2 * exp(theta[2])), 2, 2)
}

x = rexp(2)
v1 = log_dpc_hess(x, rho_0, sigma_0, p_rho, p_sigma)
v2 = numDeriv::hessian(function(x) log(dpc(x, rho_0, sigma_0, p_rho, p_sigma)), x)
v1 - v2
