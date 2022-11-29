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
