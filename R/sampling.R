
#' Sample from the global conditional extremes distribution of [x | max(x) > t],
#' using the simulation algorithm of Wadsworth et al. (2019).
#'
#' The arguments of the function are:
#' n: Number of samples.
#' n_thin: The algorithm is based on an importance sampling algorithm. First, we sample
#'   n_1 realisations, then all n_1 realisations get different weights, and
#'   we select n < n_1 of the original realisations using weighted sampling without
#'   replacement. We set n_1 = n_thin * n
#' a_func: Function a(y0, dist_to_s0) that takes in two vectors and returns a
#'   matrix where element (i, j) is the value of a(y0[j], dist_to_s0[i]).
#' b_func: Function b(y0, dist_to_s0) that takes in two vectors and returns a
#'   matrix where element (i, j) is the value of b(y0[j], dist_to_s0[i]).
#'   b_func is different from a_func in that a_func requires the distance to s0
#'   for all locations with observations, while b_func requires the distance to s0
#'   for all locations on the triangular mesh used for defining the SPDE.
#' Q: Precision matrix of the weights used for building the SPDE approximation.
#' tau: Precision of the nugget effect.
#' threshold: The threshold t used for defining the conditional extremes distribution.
#' dist_to_s0: List containing the distances to s0 for all observation locations,
#'   used for computing a. The list contains one vector of distances for each
#'   possible conditioning site in the domain of interest. Each of the distance
#'   vectors must have the same length.
#' dist_to_s0_from_mesh: List containing the distances to s0 for all locations in
#'   the mesh, used for computing b. The list contains one vector of distances for each
#'   possible conditioning site in the domain of interest. Each of the distance
#'   vectors must have the same length.
#' A: List of projection matrices used for building the SPDE approximation. The list
#'   contains one projection matrix for each possible conditioning site in the
#'   domain of interest. Each projection matrix must contain the same number of
#'   rows (equal to the lengths of dist_to_s0) and columns (equal to the lengths of
#'   dist_to_s0_from_mesh).
#' replace: bool that states if the n < n_1 chosen samples should be drawn with replacement or not.
#'   The original algorithm requires that replace = FALSE. However, weighted sampling without
#'   replacement can be extremely slow if n is large, and the algorithm is only "correct" when
#'   n / n_1 goes to infinity. However, when n / n_1 goes to infinity, sampling witout replacement
#'   is identical to sampling with replacement, which is manyfolds faster. Be aware that, for finite
#'   values of n, replace = TRUE runs the risk of sampling from distributions with point masses.
#'   However, this might not be a big problem, depending on use-case.
#' @export
wadsworth_sampling = function(n,
                              n_thin,
                              a_func,
                              b_func,
                              Q,
                              tau,
                              threshold,
                              dist_to_s0,
                              dist_to_s0_from_mesh,
                              A,
                              replace = FALSE) {

  # Compute n_1
  n_1 = as.integer(n * n_thin)
  stopifnot(n_1 > n)

  # Check that the input has the correct length
  n_loc = length(A)
  stopifnot(length(dist_to_s0) == n_loc)

  # Choose which locations the different simulations should use as conditioning sites
  s0_index = sort(sample.int(n_loc, n_1, replace = TRUE))

  # Find out how many "original" threshold exceedances we get at each conditioning site,
  # and sample them
  n_y0_per_s0 = sapply(1:n_loc, \(i) sum(s0_index == i))
  y0 = lapply(n_y0_per_s0, \(n) rexp(n) + threshold)

  # Sample from the conditional extremes model given the values of y(s0)
  res = rconditional(
    y0 = y0,
    a_func = a_func,
    b_func = b_func,
    Q = Q,
    tau = tau,
    dist_to_s0 = dist_to_s0,
    dist_to_s0_from_mesh = dist_to_s0_from_mesh,
    A = A)

  # Add the value of y0 to the correct location in the samples
  for (i in seq_along(res[-1])) {
    if (length(y0[[i]]) > 0) {
      res[[i]] = rbind(
        res[[i]][seq_len(i - 1), , drop = FALSE],
        y0[[i]],
        res[[i]][i:nrow(res[[i]]), , drop = FALSE])
    }
  }
  if (length(y0[[n_loc]]) > 0) res[[n_loc]] = rbind(res[[n_loc]], y0[[n_loc]])

  # Turn res into a matrix
  res = do.call(cbind, res)

  # Compute the importance weights
  importance_weights = apply(res, 2, function(x) 1 / sum(x > threshold))

  # Draw the n indices we are going to return
  sampling_index = sample.int(n_1, n, prob = importance_weights, replace = replace)

  # Return the n selected realisations
  res[, sampling_index]
}

#' Sample from the global conditional extremes distribution of [x | max(x) > t],
#' using the simulation algorithm of Keef et al. (2013).
#'
#' The arguments of the function are:
#' n: Number of samples.
#' a_func: Function a(y0, dist_to_s0) that takes in two vectors and returns a
#'   matrix where element (i, j) is the value of a(y0[j], dist_to_s0[i]).
#' b_func: Function b(y0, dist_to_s0) that takes in two vectors and returns a
#'   matrix where element (i, j) is the value of b(y0[j], dist_to_s0[i]).
#'   b_func is different from a_func in that a_func requires the distance to s0
#'   for all locations with observations, while b_func requires the distance to s0
#'   for all locations on the triangular mesh used for defining the SPDE.
#' Q: Precision matrix of the weights used for building the SPDE approximation.
#' tau: Precision of the nugget effect.
#' threshold: The threshold t used for defining the conditional extremes distribution.
#' dist_to_s0: List containing the distances to s0 for all observation locations,
#'   used for computing a. The list contains one vector of distances for each
#'   possible conditioning site in the domain of interest. Each of the distance
#'   vectors must have the same length.
#' dist_to_s0_from_mesh: List containing the distances to s0 for all locations in
#'   the mesh, used for computing b. The list contains one vector of distances for each
#'   possible conditioning site in the domain of interest. Each of the distance
#'   vectors must have the same length.
#' A: List of projection matrices used for building the SPDE approximation. The list
#'   contains one projection matrix for each possible conditioning site in the
#'   domain of interest. Each projection matrix must contain the same number of
#'   rows (equal to the lengths of dist_to_s0) and columns (equal to the lengths of
#'   dist_to_s0_from_mesh).
#' verbose: bool that states if the function should display progress or not.
#' @export
keef_sampling = function(n,
                         a_func,
                         b_func,
                         Q,
                         tau,
                         threshold,
                         dist_to_s0,
                         dist_to_s0_from_mesh,
                         A,
                         verbose = FALSE) {
  n_loc = length(A)
  # Check that the input has the correct length
  stopifnot(length(dist_to_s0) == n_loc && length(dist_to_s0_from_mesh) == n_loc)
  # Check that the A matrices have the same number of rows as the lengths of dist_to_s0
  stopifnot(all(sapply(dist_to_s0, length) == sapply(A, nrow)))
  stopifnot(length(unique(sapply(dist_to_s0, length))) == 1)
  # Check that all the dist_to_s0_from_mesh distances have length equal to the
  # number of nodes in the triangular mesh
  stopifnot(all(sapply(dist_to_s0_from_mesh, length) == nrow(Q)))

  # Preallocate the result matrix
  res = matrix(NA_real_, nrow = n_loc, ncol = n)

  # Choose which locations the different simulations should use as conditioning sites
  s0_index = sort(sample.int(n_loc, n, replace = TRUE))

  # Index over all samples where y(s0) is not the maximum observed value
  y0_not_max = seq_len(n)

  if (verbose) message("Start simulation")

  # Sample from the conditional model until y(s0) is the maximum value for all
  # the replications
  while (any(y0_not_max)) {

    if (verbose) {
      message("Samples remaining: ", length(y0_not_max), " / ", n)
    }

    # For each possible conditioning site, sample y(s0) for the replications
    # where y(s0) is not the maximum value
    n_y0_per_s0 = sapply(1:n_loc, function(i) sum(s0_index[y0_not_max] == i))
    y0 = lapply(n_y0_per_s0, function(n) rexp(n) + threshold)

    # Sample from the conditional extremes model given the values of y(s0)
    samples = rconditional(
      y0 = y0,
      a_func = a_func,
      b_func = b_func,
      Q = Q,
      tau = tau,
      dist_to_s0 = dist_to_s0,
      dist_to_s0_from_mesh = dist_to_s0_from_mesh,
      A = A)

    # Add the value of y0 to the correct location in the samples
    for (i in seq_along(samples[-1])) {
      if (length(y0[[i]]) > 0) {
        samples[[i]] = rbind(
          samples[[i]][seq_len(i - 1), , drop = FALSE],
          y0[[i]],
          samples[[i]][i:nrow(samples[[i]]), , drop = FALSE])
      }
    }
    if (length(y0[[n_loc]]) > 0) samples[[n_loc]] = rbind(samples[[n_loc]], y0[[n_loc]])

    # Update the rows of `res` where y(s0) is still not the maximum
    res[, y0_not_max] = do.call(cbind, samples)

    # Find the indices where y(s0) is still not the maximum value
    y0_not_max = which(sapply(
      X = seq_len(n),
      FUN = function(i) any(res[-s0_index[i], i] >= res[s0_index[i], i])))
  }

  res
}

#' Function for simulating from the conditional extremes model given the values of y(s0).
#' The arguments of the function are:
#' y0: A list of vectors containing the values of y(s0) for each of the possible
#'   conditioning sites s0. The vectors can have length 0.
#' a_func: Function a(y0, dist_to_s0) that takes in two vectors and returns a
#'   matrix where element (i, j) is the value of a(y0[j], dist_to_s0[i]).
#' b_func: Function b(y0, dist_to_s0) that takes in two vectors and returns a
#'   matrix where element (i, j) is the value of a(y0[j], dist_to_s0[i]).
#'   b_func is different from a_func in that a_func requires the distance to s0
#'   for all locations with observations, while b_func requires the distance to s0
#'   for all locations on the triangular mesh used for defining the SPDE.
#' Q: Precision matrix of the weights used for building the SPDE approximation.
#' tau: Precision of the nugget effect.
#' threshold: The threshold t used for defining the conditional extremes distribution.
#' dist_to_s0: List containing the distances to s0 for all observation locations,
#'   used for computing a. The list contains one vector of distances for each
#'   possible conditioning site in the domain of interest. Each of the distance
#'   vectors may be of different size.
#' dist_to_s0_from_mesh: List containing the distances to s0 for all locations in
#'   the mesh, used for computing b. The list contains one vector of distances for each
#'   possible conditioning site in the domain of interest. Each of the distance
#'   vectors must be of the same size.
#' A: List of projection matrices used for building the SPDE approximation. The list
#'   contains one projection matrix for each possible conditioning site in the
#'   domain of interest. Each projection matrix may contain a different number of
#'   rows, but all must have the same number of columns, which is equal to the number
#'   nodes in the triangular mesh.
#' @export
rconditional = function(y0,
                        a_func,
                        b_func,
                        Q,
                        tau,
                        dist_to_s0,
                        dist_to_s0_from_mesh,
                        A) {
  n = sapply(y0, length)
  n_loc = length(y0)
  # Check that the input has correct lengths
  lengths = c(length(A), length(dist_to_s0), length(dist_to_s0_from_mesh))
  stopifnot(all(lengths == n_loc))
  # Check that the A matrices have the same number of rows as the lengths of dist_to_s0
  stopifnot(all(sapply(dist_to_s0, length) == sapply(A, nrow)))
  # Check that all the dist_to_s0_from_mesh dist_to_s0ances have length equal to the
  # number of nodes in the triangular mesh
  stopifnot(all(sapply(dist_to_s0_from_mesh, length) == nrow(Q)))

  # Sample from the SPDE approximation
  if (is.list(Q)) {
    stopifnot(length(Q) == lengths[1])
    z = lapply(seq_along(Q), function(i) rnorm_spde(n = n[i], Q = Q[[i]]))
  } else {
    z = rnorm_spde(n = sum(n), Q = Q)
  }

  # Preallocate the result
  res = vector("list", n_loc)

  # Sample from the conditional extremes model, given the values of a, b and z
  for (i in 1:n_loc) {
    if (n[i] == 0) next
    if (!is.list(z)) {
      index = sum(n[seq_len(i - 1)]) + seq_len(n[i])
      res[[i]] = a_func(y0[[i]], dist_to_s0[[i]]) +
        as.matrix(A[[i]] %*% (b_func(y0[[i]], dist_to_s0_from_mesh[[i]]) * z[, index])) +
        rnorm(length(y0[[i]]) * length(dist_to_s0[[i]]), sd = tau^-.5)
    } else {
      res[[i]] = a_func(y0[[i]], dist_to_s0[[i]]) +
        as.matrix(A[[i]] %*% (b_func(y0[[i]], dist_to_s0_from_mesh[[i]]) * z[[i]])) +
        rnorm(length(y0[[i]]) * length(dist_to_s0[[i]]), sd = tau^-.5)
    }
  }

  res
}

#' This is a wrapper for the inla.qsample() function for sampling from a Gaussian
#' random field with sparse precision matrix Q. The wrapper allows us to get
#' reproduceable results if the seed outside the function scope is known.
#' It also suppresses warnings and removes some unnecessary column and row names
#' @export
rnorm_spde = function(n, Q, ...) {
  res = suppressWarnings(INLA::inla.qsample(n, Q, seed = round(runif(1, 0, 1e6)), ...))
  colnames(res) = NULL
  rownames(res) = NULL
  res
}
