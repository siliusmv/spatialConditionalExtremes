
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
