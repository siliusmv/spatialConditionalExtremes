
#' Compute the matrix C that solves the equation
#' C H^-1 C^T = H^-1 J H^-1
#' @export
get_C = function(H, J) {
  if (length(H) == 1 && length(J) == 1) return(get_C_univ(H, J))
  tmp = svd(solve(H))
  M_1 = tmp$u %*% diag(sqrt(tmp$d)) %*% t(tmp$v)
  tmp = svd(solve(H, J) %*% solve(H))
  M_2 = tmp$u %*% diag(sqrt(tmp$d)) %*% t(tmp$v)
  C = t(solve(M_1, M_2))
  C
}

# Compute C when H and J are scalars
get_C_univ = function(H, J) {
  M_1 = sqrt(1 / H)
  M_2 = sqrt(J / H^2)
  C = M_2 / M_1
  C
}


#' The negative Hessian of the log of a Gaussian pdf
#' @export
h_gauss = function(sigma) sigma^-2

#' The negative Hessian of the log of a log-gamma pdf.
#' If tau ~ gamma(a, b), then x = log(tau) ~ loggamma(a, b)
#' @export
h_loggamma = function(x, rate) rate * exp(x)

#' The negative Hessian of the log of a Gaussian Matérn PC prior,
#' where theta = (log(ρ), log(σ))
#' @export
h_PC = function(theta, rho_0, sigma_0, p_rho, p_sigma) {
  lambda_1 = -log(p_rho) * rho_0
  lambda_2 = -log(p_sigma) / sigma_0
  c(lambda_1 * exp(-theta[1]), lambda_2 * exp(theta[2]))
}

#' @export
h_prior = function(priors, log_tau, log_rho, log_sigma) {
  res = NULL
  if (priors$lambda[1]) res = c(res, h_gauss(priors$lambda[3]))
  if (priors$kappa[1]) res = c(res, h_gauss(priors$kappa[3]))
  if (sum(as.logical(c(priors$rho[1], priors$sigma[1]))) == 1) {
    stop("I don't know how to deal with only one of ρ and σ being fixed")
  }
  if (priors$rho[1] && priors$sigma[1]) {
    res = c(
      res,
      h_PC(c(log_rho, log_sigma), priors$rho[2], priors$sigma[2], priors$rho[3], priors$sigma[3]))
  }
  if (priors$rho_b[1]) res = c(res, h_gauss(priors$rho_b[3]))
  if (priors$tau[1]) res = c(res, h_loggamma(log_tau, priors$tau[3]))
  diag(res)
}
