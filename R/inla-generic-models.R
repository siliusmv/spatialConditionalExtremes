#' Function for defining a cgeneric model for a.
#'
#' The input variables are:
#' y0: A vector of equal length to the data (y_inla) used for running R-INLA,
#'   such that y0[i] contains the value at the conditioning site for the same
#'   time as when y_inla[i] was observed.
#' dist_to_s0: A vector of equal length to y_inla, such that dist_to_s0[i] contains
#'   the distance from the location where y_inla[i] was observed, to the conditioning
#'   site where y0[i] was observed.
#' init: Initial values for the model parameters.
#' priors: Priors for the model parameters. This is a list with one element per model parameter,
#'   with names "lambda" and "kappa". Each element contains hyperparameters
#'   for the prior of that specific model parameter.
#' is_fixed: A vector of bools, stating whether each model parameter should be fixed to its
#'   initial value or if it should be estimated.
#' debug: A boolean stating if R-INLA should print debug information or not.
#' @export
a_model = function(y0,
                   dist_to_s0,
                   init,
                   priors,
                   is_fixed = rep(FALSE, 2),
                   debug = FALSE) {
  stopifnot(length(y0) == length(dist_to_s0))
  stopifnot(length(is_fixed) == length(init))
  stopifnot(is.logical(is_fixed))
  stopifnot(length(init) == 2)
  stopifnot(all(c("lambda", "kappa")[!is_fixed] %in% names(priors)))
  stopifnot(all(sapply(priors, length) == 2))

  args = list(debug = debug)

  # The name and location of the required c-function for defining the cgeneric model
  args$model = "a_model"
  args$shlib = file.path(cgeneric_dir(), "a.so")

  # Put all the arguments into args in the order defined in the c-function
  args$n = length(y0)
  args$is_fixed = as.integer(is_fixed)
  args$y0 = y0
  args$dist_to_s0 = dist_to_s0
  args$init = init
  args$lambda_prior = priors$lambda
  args$kappa_prior = priors$kappa

  # Define the model
  do.call(INLA::inla.cgeneric.define, args)
}

#' Function for defining a cgeneric model for Z_b.
#'
#' The input variables are:
#' n: A vector of length equal to the number of unique conditioning sites used
#'   for inference. Element nr. i of n contains the number of threshold exceedances
#'   found at conditioning site nr. i.
#' spde: A list with length equal to n, of inla.spde2 objects that each contains
#'   the matrices M0, M1, M2, B0, B1 and B2,
#'   which are necessary for computing the precision matrix of the Gaussian Matérn
#'   field used for creating the SPDE approximation.
#' init: Initial values for the model parameters
#' priors: Priors for the parameters of Z_b. This is a list with one element per model parameter,
#'   with names "sigma", "rho", "s0", "lambda" and "kappa". Each element contains hyperparameters
#'   for the prior of that specific model parameter.
#' dist_to_s0: A list of length equal to the n, where element nr. i contains the distances
#'   from all the mesh nodes in spde[i] to conditioning site nr. i.
#' is_fixed: A vector of bools, stating whether each model parameter should be fixed to its
#'   initial value or if it should be estimated.
#' debug: A boolean stating if R-INLA should print debug information or not.
#' @export
spde_b_model_case_study = function(n,
                                   spde,
                                   init,
                                   priors,
                                   dist_to_s0,
                                   is_fixed = rep(FALSE, 5),
                                   debug = FALSE) {
  if (!all(class(spde) == "list")) spde = list(spde)
  stopifnot(length(n) == length(spde))
  stopifnot(length(n) == length(dist_to_s0))
  stopifnot(all(sapply(spde, `[[`, "n.spde") == sapply(dist_to_s0, length)))
  stopifnot(is.logical(is_fixed))
  stopifnot(length(init) == 5)
  stopifnot(length(is_fixed) == length(init))
  stopifnot(all(c("sigma", "rho", "b0", "lambda", "kappa")[!is_fixed] %in% names(priors)))
  stopifnot(all(sapply(priors, length) == 2))
  stopifnot(all(sapply(dist_to_s0, min) == 0))

  # Locate and remove the s0_locs from dist_to_s0
  s0_location = sapply(dist_to_s0, function(x) which(x == 0))
  stopifnot(is.integer(s0_location))
  for (i in seq_along(dist_to_s0)) dist_to_s0[[i]] = dist_to_s0[[i]][-s0_location[i]]

  # Start building the argument list for the cgeneric model
  args = list(debug = debug)

  # The name and location of the required c-function for defining the cgeneric model
  args$model = "spde_b_model_case_study"
  args$shlib = file.path(cgeneric_dir(), "b.so")

  # Put all the arguments into args in the order defined in the c-function
  args$n = sum(sapply(dist_to_s0, length) * n)
  args$s0_index = rep(seq_along(n), n) - 1L # 0-indexation in C
  args$is_fixed = as.integer(is_fixed)
  args$s0_location = s0_location - 1L # 0-indexation in C
  args$init = init
  if (!is_fixed[1]) args$sigma_prior = priors$sigma
  if (!is_fixed[2]) args$rho_prior = priors$rho
  if (!is_fixed[3]) args$b0_prior = priors$b0
  if (!is_fixed[4]) args$lambda_prior = priors$lambda
  if (!is_fixed[5]) args$kappa_prior = priors$kappa

  for (i in seq_along(n)) {
    args[[paste0("dist_to_s0_", i)]] = dist_to_s0[[i]]
    args[[paste0("B0_", i)]] = spde[[i]]$param.inla$B0
    args[[paste0("B1_", i)]] = spde[[i]]$param.inla$B1
    args[[paste0("B2_", i)]] = spde[[i]]$param.inla$B2
    args[[paste0("M0_", i)]] = spde[[i]]$param.inla$M0
    args[[paste0("M1_", i)]] = spde[[i]]$param.inla$M1
    args[[paste0("M2_", i)]] = spde[[i]]$param.inla$M2
  }

  # Define the model
  do.call(INLA::inla.cgeneric.define, args)
}

#' Function for defining a cgeneric model for Z_b.
#'
#' The input variables are:
#' n: A vector of length equal to the number of unique conditioning sites used
#'   for inference. Element nr. i of n contains the number of threshold exceedances
#'   found at conditioning site nr. i.
#' spde: A list with length equal to n, of inla.spde2 objects that each contains
#'   the matrices M0, M1, M2, B0, B1 and B2,
#'   which are necessary for computing the precision matrix of the Gaussian Matérn
#'   field used for creating the SPDE approximation.
#' init: Initial values for the model parameters
#' priors: Priors for the model parameters. This is a list with one element per model parameter,
#'   with names "sigma", "rho", "beta0", "lambda" and "kappa". Each element contains hyperparameters
#'   for the prior of that specific model parameter.
#' dist_to_s0: A list of length equal to the n, where element nr. i contains the distances
#'   from all the mesh nodes in spde[i] to conditioning site nr. i.
#' is_fixed: A vector of bools, stating whether each model parameter should be fixed to its
#'   initial value or if it should be estimated.
#' debug: A boolean stating if R-INLA should print debug information or not.
#' @export
spde_b_model_simulation = function(n,
                                   y0,
                                   spde,
                                   init,
                                   priors,
                                   dist_to_s0,
                                   is_fixed = rep(FALSE, 5),
                                   debug = FALSE) {
  if (!all(class(spde) == "list")) spde = list(spde)
  stopifnot(length(n) == length(spde))
  stopifnot(length(n) == length(dist_to_s0))
  stopifnot(sum(n) == length(y0))
  stopifnot(all(sapply(spde, `[[`, "n.spde") == sapply(dist_to_s0, length)))
  stopifnot(is.logical(is_fixed))
  stopifnot(length(init) == 5)
  stopifnot(length(is_fixed) == length(init))
  stopifnot(all(c("sigma", "rho", "beta0", "lambda", "kappa")[!is_fixed] %in% names(priors)))
  stopifnot(all(sapply(priors, length) == 2))
  stopifnot(all(sapply(dist_to_s0, min) == 0))

  # Locate and remove the s0_locs from dist_to_s0
  s0_location = sapply(dist_to_s0, function(x) which(x == 0))
  stopifnot(is.integer(s0_location))
  for (i in seq_along(dist_to_s0)) dist_to_s0[[i]] = dist_to_s0[[i]][-s0_location[i]]

  # Start building the argument list for the cgeneric model
  args = list(debug = debug)

  # The name and location of the required c-function for defining the cgeneric model
  args$model = "spde_b_model_simulation"
  args$shlib = file.path(cgeneric_dir(), "b.so")

  # Put all the arguments into args in the order defined in the c-function
  args$n = sum(sapply(dist_to_s0, length) * n)
  args$s0_index = rep(seq_along(n), n) - 1L # 0-indexation in C
  args$is_fixed = as.integer(is_fixed)
  args$s0_location = s0_location - 1L # 0-indexation in C
  args$init = init
  args$y0 = y0
  if (!is_fixed[1]) args$sigma_prior = priors$sigma
  if (!is_fixed[2]) args$rho_prior = priors$rho
  if (!is_fixed[3]) args$beta0_prior = priors$beta0
  if (!is_fixed[4]) args$lambda_prior = priors$lambda
  if (!is_fixed[5]) args$kappa_prior = priors$kappa

  for (i in seq_along(n)) {
    args[[paste0("dist_to_s0_", i)]] = dist_to_s0[[i]]
    args[[paste0("B0_", i)]] = spde[[i]]$param.inla$B0
    args[[paste0("B1_", i)]] = spde[[i]]$param.inla$B1
    args[[paste0("B2_", i)]] = spde[[i]]$param.inla$B2
    args[[paste0("M0_", i)]] = spde[[i]]$param.inla$M0
    args[[paste0("M1_", i)]] = spde[[i]]$param.inla$M1
    args[[paste0("M2_", i)]] = spde[[i]]$param.inla$M2
  }

  # Define the model
  do.call(INLA::inla.cgeneric.define, args)
}

#' Function for defining a cgeneric model for the entire spatial conditional extremes model.
#'
#' The input variables are:
#' spde: An inla.spde2 object that contains the matrices M0, M1, M2, B0, B1 and B2,
#'   that are necessary for computing the precision matrix of the Gaussian Matérn
#'   field used for creating the SPDE approximation.
#' y0: A vector of equal length to the data (y_inla) used for running R-INLA,
#'   such that y0[i] contains the value at the conditioning site for the same
#'   time as when y_inla[i] was observed.
#' n: A vector of length equal to the number of unique conditioning sites used
#'   for inference. Element nr. i of n contains the number of threshold exceedances
#'   found at conditioning site nr. i.
#' init: Initial values for the subset of (log(λ), log(κ), log(ρ), log(σ), log(ρ_b), log(τ))
#'   that will be estimated. See the documentation for the prior variable for more information.
#' priors: Priors for the six parameters of the model. This is a list of length 6,
#'   with names "lambda", "kappa", "rho", "sigma", "rho_b" and "tau".
#'   The list elements are vectors of
#'   length 2 or 3. If the first element of a vector is 0, then we will not estimate
#'   that parameter, but we will fix it equal to the value given in the second element
#'   of the vector. If the first element is not 0, then the vector has length 3, and
#'   the two remaining values are used for placing a prior on that parameter.
#'   log(λ), log(κ) and log(ρ_b) are given Gaussian priors.
#'   ρ and σ are given a joint PC prior, see ?inla.spde2.pcmatern for more information.
#'   τ is given a gamma prior.
#' dist_to_s0: A vector of equal length to y_inla, such that dist_to_s0[i] contains
#'   the distance from the location where y_inla[i] was observed, to the conditioning
#'   site where y0[i] was observed.
#' dist_to_s0_from_mesh: A matrix of dimension (n x m), where m is the number of
#'   nodes in the triangular mesh, and n is the number of unique conditioning sites.
#'   Row nr. i of the matrix contains the distances from all mesh nodes to the
#'   ith conditioning site. dist_to_s0_from_mesh can also be a vector, if only
#'   one conditioning site is used.
#' C: NULL or a matrix that can be used to transform the prior.
#' theta_star: The mode of the likelihood. This is only used is C is non-null,
#'   and it must be non-null if C is non-null.
#' debug: A boolean stating if R-INLA should print debug information or not.
#' @export
spatial_conditional_extremes_cgeneric_model = function(spde,
                                                       y0,
                                                       n,
                                                       init,
                                                       priors,
                                                       dist_to_s0,
                                                       dist_to_s0_from_mesh,
                                                       C = NULL,
                                                       theta_star = NULL,
                                                       debug = FALSE) {
  stopifnot(all(c("lambda", "kappa", "rho", "sigma", "rho_b", "tau") %in% names(priors)))
  if (!is.matrix(dist_to_s0_from_mesh)) {
    dist_to_s0_from_mesh = matrix(dist_to_s0_from_mesh, nrow = 1)
  }
  stopifnot(length(n) == nrow(dist_to_s0_from_mesh))
  stopifnot(spde$n.spde == ncol(dist_to_s0_from_mesh))
  stopifnot(length(y0) == length(dist_to_s0))
  stopifnot(length(priors) == 6)
  args = list(debug = debug)

  # The name and location of the required c-function for defining the cgeneric model
  args$model = "inla_cgeneric_spatial_conditional_extremes"
  args$shlib = file.path(cgeneric_dir(), "spatial_conditional_extremes.so")

  # Check that the priors have the correct lengths
  for (i in seq_along(priors)) {
    if (priors[[i]][1] == 0) {
      stopifnot(length(priors[[i]]) == 2)
    } else {
      stopifnot(length(priors[[i]]) == 3)
    }
  }

  # Check that the init vector has the correct length
  n_fixed = sum(sapply(priors, `[`, 1) == 0)
  stopifnot(length(init) == 6 - n_fixed)

  # Put all the arguments into args in the order defined in the c-function
  args$n = spde$n.spde * sum(n) + length(y0) * 2
  args$s0_index = rep(seq_along(n), n) - 1L
  args$B0 = spde$param.inla$B0
  args$B1 = spde$param.inla$B1
  args$B2 = spde$param.inla$B2
  args$dist_to_s0_from_mesh = dist_to_s0_from_mesh
  args$M0 = spde$param.inla$M0
  args$M1 = spde$param.inla$M1
  args$M2 = spde$param.inla$M2
  args$y0 = y0
  args$dist_to_s0 = dist_to_s0
  args$init = init
  args$lambda_prior = priors$lambda
  args$kappa_prior = priors$kappa
  args$rho_prior = priors$rho
  args$sigma_prior = priors$sigma
  args$rho_b_prior = priors$rho_b
  args$tau_prior = priors$tau

  if (is.null(C)) {
    args$transform_prior = 0L
  } else {
    args$transform_prior = 1L
    stopifnot(!is.null(theta_star))
    stopifnot(all(dim(C) == length(init)))
    stopifnot(length(init) == length(theta_star))
    args$C = C
    args$theta_star = theta_star
    args$log_C_det = as.numeric(determinant(C)$modulus)
  }

  # Define the model
  do.call(INLA::inla.cgeneric.define, args)
}

# #' Function for defining a cgeneric model for Z_b.
# #'
# #' The input variables are:
# #' spde: An inla.spde2 object that contains the matrices M0, M1, M2, B0, B1 and B2,
# #'   that are necessary for computing the precision matrix of the Gaussian Matérn
# #'   field used for creating the SPDE approximation.
# #' n: A vector of length equal to the number of unique conditioning sites used
# #'   for inference. Element nr. i of n contains the number of threshold exceedances
# #'   found at conditioning site nr. i.
# #' init: Initial values for the subset of (log(ρ), log(σ), log(ρ_b)) that will be
# #'   estimated. See the documentation for the prior variable for more information.
# #' priors: Priors for the three parameters of Z_b. This is a list of length 3,
# #'   with names "rho", "sigma" and "rho_b". The three list elements are vectors of
# #'   length 2 or 3. If the first element of a vector is 0, then we will not estimate
# #'   that parameter, but we will fix it equal to the value given in the second element
# #'   of the vector. If the first element is not 0, then the vector has length 3, and
# #'   the two remaining values are used for placing a prior on that parameter.
# #'   ρ and σ are given a joint PC prior, see ?inla.spde2.pcmatern for more information.
# #'   log(ρ_b) is given a Gaussian prior. If e.g.
# #'   priors = list(rho = c(1, 60, .95), sigma = c(1, 5, .05), rho_b = c(1, 2, 3)), then
# #'   ρ and σ are given a PC prior such that P(ρ < 60) = .95 and P(σ > 5) = .05,
# #'   while log(ρ_b) is given a Gaussian prior with mean 2 and standard deviation 3.
# #'   If, on the other hand, priors$rho_b = c(0, 2), then the value of log(ρ_b) is
# #'   fixed equal to 2.
# #' dist_to_s0: A matrix of dimension (n x m), where m is the number of nodes in the
# #'   triangular mesh, and n is the number of unique conditioning sites. Row nr. i
# #'   of the matrix contains the distances from all mesh nodes to the ith conditioning
# #'   site. dist_to_s0 can also be a vector, if only one conditioning site is used.
# #' debug: A boolean stating if R-INLA should print debug information or not.
# #' @export
# spde_b_func = function(spde,
#                        n,
#                        init,
#                        priors,
#                        dist_to_s0,
#                        debug = FALSE) {
#   stopifnot(all(c("rho", "sigma", "rho_b") %in% names(priors)))
#   if (!is.matrix(dist_to_s0)) dist_to_s0 = matrix(dist_to_s0, nrow = 1)
#   stopifnot(length(n) == nrow(dist_to_s0))
#   stopifnot(spde$n.spde == ncol(dist_to_s0))
#   stopifnot(length(priors) == 3)
#   args = list(debug = debug)
# 
#   # The name and location of the required c-function for defining the cgeneric model
#   args$model = "inla_cgeneric_spde_model_with_b_func"
#   args$shlib = file.path(cgeneric_dir(), "b.so")
# 
#   # Check that the priors have the correct lengths
#   for (i in seq_along(priors)) {
#     if (priors[[i]][1] == 0) {
#       stopifnot(length(priors[[i]]) == 2)
#     } else {
#       stopifnot(length(priors[[i]]) == 3)
#     }
#   }
# 
#   # Check that the init vector has the correct length
#   n_fixed = sum(sapply(priors, `[`, 1) == 0)
#   stopifnot(length(init) == 3 - n_fixed)
# 
#   # Put all the arguments into args in the order defined in the c-function
#   args$n = spde$n.spde * sum(n)
#   args$B0 = spde$param.inla$B0
#   args$B1 = spde$param.inla$B1
#   args$B2 = spde$param.inla$B2
#   args$M0 = spde$param.inla$M0
#   args$M1 = spde$param.inla$M1
#   args$M2 = spde$param.inla$M2
#   args$init = init
#   args$rho_prior = priors$rho
#   args$sigma_prior = priors$sigma
#   args$rho_b_prior = priors$rho_b
#   args$dist_to_s0 = dist_to_s0
#   args$s0_index = rep(seq_along(n), n)
# 
#   # Define the model
#   do.call(INLA::inla.cgeneric.define, args)
# }
