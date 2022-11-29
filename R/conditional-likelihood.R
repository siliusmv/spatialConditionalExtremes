
#' Comput the composite log-likelihood for the global conditional extremes model.
#'
#' The input variables are:
#' y: A list of matrices. Matrix nr. i has dimension (d_i x n_i) and contains n_i
#'   d_i-dimensional vectors of observations for conditioning site nr. i.
#' y0: A list of vectors. Vector nr. i has dimension n_i and contains the threshold
#'   exceedances from conditioning site nr. i.
#' a_func: A function a(y0, dist_to_s0) that takes in two vectors and returns a
#'   matrix where element (i, j) is the value of a(y0[j]), dist_to_s0[i]).
#' b_func: A function b(y0, dist_to_s0) that takes in two vectors and returns a
#'   matrix where element (i, j) is the value of b(y0[j], dist_to_s0[i]).
#'   b_func is different from a_func in that a_func requires the distance to s0
#'   for all locations with observations, while b_func requires the distance to s0
#'   for all locations on the triangular mesh used for defining the SPDE.
#' Q: Either a precision matrix, or a list of precision matrices, of the mesh for
#'   each conditioning site. The list must have length 1 or length(y0)
#' sigma: The covariance matrix of the weights used for building the SPDE approximation.
#' tau: The precision of the nugget effect.
#' dist_to_s0: List containing the distances to s0 for all observation locations,
#'   used for computing a. The list contains one vector of distances for each
#'   conditioning site. Vector nr. i has dimension d_i.
#' dist_to_s0_from_mesh: List containing the distances to s0 for all locations in
#'   the mesh, used for computing b. The list contains one vector of distances for each
#'   conditioning site.
#' A: List of projection matrices used for building the SPDE approximation. The list
#'   contains one projection matrix for each conditioning site.
#'   Each projection matrix must contain the same number of columns, equal to the
#'   lengths of dist_to_s0_from_mesh. Projection matrix nr. i has d_i rows.
#' n_cores: Number of cores to use for computing the log-likelihood in parallel. We
#'   cannot use more cores than the number of conditioning sites.
#' na.rm: Should NA's be removed or not? If na.rm = FALSE, then
#'   any column of x that returns an NA value will result in an NA value in the output.
#'   If na.rm = TRUE, then we remove the NA values before computing the likelihood.
#'   So if e.g. a column has 3 NA variables, then we remove these, and then we compute
#'   the likelihood for a (d-3)-dimensional Gaussian random variable.
#' use_r: A boolean, mostly used for testing the code. Should we perform the actual computations
#'   using R code or using Rcpp code?
#' @export
loglik_conditional = function(y,
                              y0,
                              a_func,
                              b_func,
                              Q,
                              tau,
                              dist_to_s0,
                              dist_to_s0_from_mesh,
                              A,
                              n_cores = 1,
                              na.rm = TRUE,
                              use_r = FALSE) {
  # Compute covariance matrices
  if (!is.list(Q)) Q = list(Q)
  Sigma = lapply(Q, function(q) as.matrix(Matrix::solve(q)))
  if (length(Sigma) == 1) Sigma = rep(Sigma, length(y0))
  # Ensure that the input has correct lengths/dimensions
  lengths = c(
    length(y0), length(y), length(A),
    length(dist_to_s0), length(dist_to_s0_from_mesh), length(Sigma))
  stopifnot(all(lengths == lengths[1]))
  stopifnot(all(sapply(A, nrow) == sapply(dist_to_s0, length)))
  stopifnot(all(sapply(A, ncol) == sapply(dist_to_s0_from_mesh, length)))

  # Compute the log-likelihood, possibly in parallel
  res = parallel::mclapply(
    X = seq_along(A),
    mc.cores = n_cores,
    FUN = function(i) {
      # Compute b
      b = b_func(y0[[i]], dist_to_s0_from_mesh[[i]])

      # Do the values of b depend on y or are they constant given the distance?
      const_b = all(apply(b, 1, function(x) length(unique(x)) == 1))

      # Create a list of arguments that are used for computing the log-likelihood
      args = list(
        x = y[[i]] - a_func(y0[[i]], dist_to_s0[[i]]),
        A = A[[i]],
        sigma0 = Sigma[[i]],
        nugget = 1 / tau,
        logd = TRUE,
        na_rm = na.rm)

      # If b does not depend on y, we can speed up computations considerably by using
      # a different function for computing the log-likelihood
      if (const_b) {
        args$b = b[, 1]
        func = ifelse(use_r, dconditional_r_no_beta, dconditional_no_beta)
      } else {
        args$B = b
        func = ifelse(use_r, dconditional_r, dconditional)
      }

      # Perform the actual computations
      as.numeric(do.call(func, args))
    })
  unlist(res)
}

#' Compute the (log-)likelihood of the conditional extremes distribution when
#' a and b are given. We assume that a has already been subtracted from x,
#' i.e. x = (y - a) for some observations y.
#' This function is mainly written to evaluate the speed and correctness of the Rcpp
#' function dconditional()
#'
#' The input variables are:
#' x: an (n x d)-dimensional matrix of observations where a has already been subtracted.
#' A: The projection matrix used for building the SPDE approximation
#' B: An (m x n)-dimensional matrix containing the values of b at the m triangular mesh nodes,
#'   for each of the n threshold exceedances at the conditioning sites of interest
#' sigma0: The covariance matrix of the m Gaussian random variables in the triangular mesh.
#' nugget: The variance of the nugget effect.
#' logd: Boolean explaining if we should return the log-likelihood or the likelihood.
#' na_rm: Boolean explaining how we should treat NA values. If na_rm = false, then
#'   any column of x that returns an NA value will result in an NA value in the output.
#'   If na_rm = true, then we remove the NA values before computing the likelihood.
#'   So if e.g. a column has 3 NA variables, then we remove these, and then we compute
#'   the likelihood for a (d-3)-dimensional Gaussian random variable.
dconditional_r = function(x,
                          A,
                          B,
                          sigma0,
                          nugget,
                          logd = TRUE,
                          na_rm = TRUE) {
  m = nrow(B)
  n = ncol(x)
  res = rep(NA_real_, n)

  # Loop over all columns in x
  for (i in 1:ncol(x)) {

    # Compute the value of sigma given B[, i]
    B_tmp = Matrix::Diagonal(m, B[, i])
    sigma = as.matrix(A %*% B_tmp %*% sigma0 %*% B_tmp %*% Matrix::t(A)) + diag(nugget, nrow(A))

    # Locate indices for all non-NA observations
    if (na_rm) {
      good_index = which(!is.na(x[, i]))
    } else {
      good_index = seq_along(x[, i])
    }

    # Compute the log-likelihood
    res[i] = mvtnorm::dmvnorm(x[good_index, i], sigma = sigma[good_index, good_index], log = TRUE)
  }

  res
}


#' Compute the (log-)likelihood of the conditional extremes distribution when
#' a and b are given. We assume that a has already been subtracted from x,
#' i.e. x = (y - a) for some observations y.
#' Further, we assume that b does not depend on y0, meaning that we only need one
#' value of b for each distance, instead of a matrix B that has one value for each
#' pair of distance and threshold exceedance y0.
#' This function is mainly written to evaluate the speed and correctness of the Rcpp
#' function dconditional_no_beta()
#'
#' The input variables are:
#' x: an (n x d)-dimensional matrix of observations where a has already been subtracted.
#' A: The projection matrix used for building the SPDE approximation
#' b: An m-dimensional vector containing the values of b at the m triangular mesh nodes.
#' sigma0: The covariance matrix of the m Gaussian random variables in the triangular mesh.
#' nugget: The variance of the nugget effect.
#' logd: Boolean explaining if we should return the log-likelihood or the likelihood.
#' na_rm: Boolean explaining how we should treat NA values. If na_rm = false, then
#'   any column of x that returns an NA value will result in an NA value in the output.
#'   If na_rm = true, then we remove the NA values before computing the likelihood.
#'   So if e.g. a column has 3 NA variables, then we remove these, and then we compute
#'   the likelihood for a (d-3)-dimensional Gaussian random variable.
dconditional_r_no_beta = function(x,
                                  A,
                                  b,
                                  sigma0,
                                  nugget,
                                  logd = TRUE,
                                  na_rm = TRUE) {
  m = length(b)
  n = ncol(x)
  res = rep(NA_real_, n)

  # Compute the value of sigma given b
  B = Matrix::Diagonal(m, b)
  sigma = as.matrix(A %*% B %*% sigma0 %*% B %*% Matrix::t(A)) + diag(nugget, nrow(A))

  # Loop over all columns in x
  for (i in 1:ncol(x)) {

    # Locate indices for all non-NA observations
    if (na_rm) {
      good_index = which(!is.na(x[, i]))
    } else {
      good_index = seq_along(x[, i])
    }

    # Compute the log-likelihood
    res[i] = mvtnorm::dmvnorm(x[good_index, i], sigma = sigma[good_index, good_index], log = TRUE)
  }

  res
}
