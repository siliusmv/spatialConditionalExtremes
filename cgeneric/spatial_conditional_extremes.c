
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "cgeneric.h"
#include "smat-operations.c"
#include "spde-precision.c"

/*
double * spde_ab_model(inla_cgeneric_cmd_tp cmd,
		       double * theta,
		       inla_cgeneric_data_tp * data) {

  double * ret = NULL;

  // ==========================================================
  // Assert that all the input variables look as they shold
  // ==========================================================

  // Number of observations
  assert(!strcmp(data->ints[0]->name, "n"));
  int n = data->ints[0]->ints[0];
  assert(n > 0);

  // Which s0 do the different replications come from?
  assert(!strcmp(data->ints[2]->name, "s0_index"));
  int n_repl = data->ints[2]->len;
  int *s0_index = data->ints[2]->ints;

  assert(!strcmp(data->ints[3]->name, "is_fixed"));
  int * is_fixed = data->ints[3]->ints;
  int n_theta = data->ints[3]->len;
  int n_free = 0;
  for (int i = 0; i < n_theta; ++i) {
    if (!is_fixed[i]) ++n_free;
  }

  // Look at the matrices required for building the meshes
  int n_mesh = data->n_smats / 3;
  assert(data->n_smats == n_mesh * 3); // Ensure that we did not get rounding when dividing
  assert(data->n_smats == data->n_mats); // We need the same number of mats and smats

  // Initial values and priors
  assert(!strcmp(data->doubles[0]->name, "init"));
  double * init = data->doubles[0]->doubles;
  assert(data->doubles[0]->len == n_theta);

  assert(!strcmp(data->doubles[1]->name, "y0"));
  double * y0 = data->doubles[1]->doubles;
  int n_y0 = data->doubles[1]->len;

  assert(!strcmp(data->doubles[2]->name, "unique_y0"));
  assert(data->doubles[2]->len == n_repl);
  double * unique_y0 = data->doubles[2]->doubles;

  assert(!strcmp(data->doubles[3]->name, "dist_to_s0"));
  assert(data->doubles[3]->len == n_y0);
  double * dist_to_s0 = data->doubles[3]->doubles;

  int prior_start = 4; // This is the first element of data->doubles that contains a prior
  int dist_start = prior_start + n_free; // This is the first element of data->doubles that contains dist_to_s0_from_mesh
  assert(data->n_doubles == n_mesh + dist_start);

  // ==========================================================
  // Initiate the correct parameter values, depending on whether
  // the parameters are fixed or estimated.
  // ==========================================================
  double lambda, kappa, log_sigma, rho_b, beta, log_rho;
  assert(n_theta == 6);
  if (theta) {
    double theta_full[6];
    int count = 0;
    for (int i = 0; i < n_theta; ++i) {
      if (is_fixed[i]) {
	theta_full[i] = init[i];
      } else {
	theta_full[i] = theta[count];
	++count;
      }
    }
    assert(count == n_free);
    lambda = exp(theta_full[0]);
    kappa = exp(theta_full[1]);
    log_sigma = theta_full[2];
    rho_b = exp(theta_full[3]);
    beta = exp(theta_full[4]);
    log_rho = theta_full[5];
  } else {
    lambda = kappa = log_sigma = rho_b = beta = log_rho = NAN;
  }

  // ==========================================================
  // This switch statement is the required method for implementing
  // cgeneric models in R-INLA
  // ==========================================================
  switch (cmd) {
  case INLA_CGENERIC_VOID:
    {
      assert(!(cmd == INLA_CGENERIC_VOID));
      break;
    }

  case INLA_CGENERIC_GRAPH:
    {
      // return a vector of indices with format
      // c(n, M, ii, jj)
      // where ii<=jj, ii is non-decreasing and jj is non-decreasing within each ii,
      // and M is the length of ii

      // Compute all the different precision matrices
      inla_cgeneric_smat_tp * all_precisions =
	spde_precision_multimesh(log_rho, log_sigma, data->mats, data->smats, n_mesh);

      inla_cgeneric_smat_tp precision = copy_smat(all_precisions + s0_index[0]);
      for (int i = 1; i < n_repl; ++i) {
	smat_smat_merge(&precision, all_precisions + s0_index[i]);
      }

      // Create the graph for a()
      diag_mat_tp diagonal_graph_part = diag(1, n_y0);

      // Merge the two graph matrices together
      smat_diag_merge(&precision, &diagonal_graph_part);

      // Extract all the necessary i and j values from the block matrix
      ret = Malloc(2 + 2 * precision.n, double);
      ret[0] = n;
      ret[1] = precision.n;
      for (int i = 0; i < precision.n; i++) {
	ret[2 + i] = precision.i[i];                                  // i
	ret[2 + precision.n + i] = precision.j[i];                    // j
      }

      // Avoid memory leaks
      for (int i = 0; i < n_mesh; ++i) {
	free_smat(all_precisions + i);
      }
      free_smat(&precision);
      free_diag(&diagonal_graph_part);

      break;
    }

  case INLA_CGENERIC_Q:
    {

      // Compute all the different precision matrices
      inla_cgeneric_smat_tp * all_precisions =
	spde_precision_multimesh(log_rho, log_sigma, data->mats, data->smats, n_mesh);

      inla_cgeneric_smat_tp precision = copy_smat(all_precisions + s0_index[0]);
      for (int i = 1; i < n_repl; ++i) {
	smat_smat_merge(&precision, all_precisions + s0_index[i]);
      }

      // Compute the values of b for each distance from a conditioning site
      double ** unique_sd_vals = Malloc(n_mesh, double *);
      double ** unique_alpha_vals = Malloc(n_mesh, double *);
      for (int i = 0; i < n_mesh; ++i) {
	unique_sd_vals[i] = Malloc(data->doubles[dist_start + i]->len, double);
	unique_alpha_vals[i] = Malloc(data->doubles[dist_start + i]->len, double);
	for (int j = 0; j < data->doubles[dist_start + i]->len; ++j) {
	  unique_sd_vals[i][j] = sqrt(1 - exp(-4 * data->doubles[dist_start + i]->doubles[j] / rho_b));
	  unique_alpha_vals[i][j] = exp(-beta * pow(data->doubles[dist_start + i]->doubles[j] / lambda, kappa));
	}
      }

      double * unique_transformed_y0 = Malloc(n_repl, double);
      for (int i = 0; i < n_repl; ++i) {
	unique_transformed_y0[i] = pow(unique_y0[i], beta);
      }

      // Create a diagonal matrix containing the values of 1 / b, that can be used
      // for rescaling the precision matrix with the correct b values
      diag_mat_tp b_inv;
      b_inv.x = Malloc(precision.nrow, double);
      b_inv.dim = precision.nrow;
      int count = 0;
      for (int i = 0; i < n_repl; ++i) {
	for (int j = 0; j < data->doubles[dist_start + s0_index[i]]->len; ++j) {
	  b_inv.x[count] = unique_sd_vals[s0_index[i]][j] *
	    (1 + unique_alpha_vals[s0_index[i]][j] * unique_transformed_y0[i]);
	  if (b_inv.x[count] < .0000000001) b_inv.x[count] = .0000000001; // Ensure that we don't divide by zero
	  b_inv.x[count] = 1 / b_inv.x[count];
	  ++count;
	}
      }

      // Multiply the diagonal matrix with the precision matrix: B_inv * Q * B_inv
      diag_smat_diag_mult(&precision, &b_inv);

      // Create the precision matrix for the a() part
      diag_mat_tp diag_precision = diag(exp(15), n_y0);

      // Add the a diagonal precision to the precision matrix
      smat_diag_merge(&precision, &diag_precision);

      // optimized format
      // return c(n, M, Qij) in the same order as defined in INLA_CGENERIC_GRAPH
      // where M is the length of Qij
      ret = Malloc(2 + precision.n, double);
      ret[0] = -1;                                   // code for optimized output
      ret[1] = precision.n;                          // number of (i <= j)
      for (int i = 0; i < precision.n; i++) {
	ret[2 + i] = precision.x[i];
      }

      // Avoid memory leaks
      for (int i = 0; i < n_mesh; ++i) {
	free(unique_sd_vals[i]);
	free(unique_alpha_vals[i]);
	free_smat(all_precisions + i);
      }
      free(unique_sd_vals);
      free(unique_alpha_vals);
      free(unique_transformed_y0);
      free_smat(&precision);
      free_diag(&b_inv);
      free_diag(&diag_precision);

      break;
    }

  case INLA_CGENERIC_MU:
    {
      // return (N, mu)
      // if N==0 then mu is not needed as its taken to be mu[]==0
      ret = Malloc(1 + n, double);
      ret[0] = n;

      int count = 1;
      for (int i = 0; i < n_repl; ++i) {
	for (int j = 0; j < data->doubles[dist_start + s0_index[i]]->len; ++j) {
	  ret[count] = 0;
	  ++count;
	}
      }
      for (int i = 0; i < n_y0; ++i) {
	ret[count] = y0[i] * exp(-pow(dist_to_s0[i] / lambda, kappa));
	++count;
      }

      break;
    }

  case INLA_CGENERIC_INITIAL:
    {
      // return c(M, initials)
      // where M is the number of hyperparameters

      // The initial values depend on how many parameters that are fixed or estimated.
      ret = Malloc(n_free + 1, double);
      ret[0] = n_free;
      int count = 0;
      for (int i = 0; i < n_theta; ++i) {
	if (!is_fixed[i]) {
	  ret[count + 1] = init[i];
	  ++count;
	}
      }
      assert(count == n_free);
      break;
    }

  case INLA_CGENERIC_LOG_NORM_CONST:
    {
      // return c(NORM_CONST) or a NULL-pointer if INLA should compute it by itself

      ret = NULL;
      break;
    }

  case INLA_CGENERIC_LOG_PRIOR:
    {
      // return c(LOG_PRIOR)

      ret = Calloc(1, double);

      // This is an approximation for -0.5 * log(2π)
      double gaussian_const = -0.91893853320467;

      int count = 0;
      double * prior;
      if (!is_fixed[0]) { // lambda
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(lambda) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[1]) { // kappa
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(kappa) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[2]) { // sigma
	prior = data->doubles[prior_start + count]->doubles;
	double lambda1 = -log(prior[1]) / prior[0];
	ret[0] += log(lambda1) - lambda1 * exp(log_sigma) + log_sigma;
	++count;
      }
      if (!is_fixed[3]) { // rho_b
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(rho_b) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[4]) { // beta
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(beta) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[5]) { // rho
	prior = data->doubles[prior_start + count]->doubles;
	double lambda0 = -log(prior[1]) * prior[0];
	ret[0] += log(lambda0) - lambda0 * exp(-log_rho) - log_rho;
	++count;
      }
      assert(count == n_free);

      break;
    }

  case INLA_CGENERIC_QUIT:
  default:
    break;
  }

  return (ret);
}

double *inla_cgeneric_spatial_conditional_extremes(inla_cgeneric_cmd_tp cmd,
						   double * theta,
						   inla_cgeneric_data_tp * data) {

  double *ret = NULL;

  // This is an approximation for -0.5 * log(2π)
  double gaussian_const = -0.91893853320467;
  // The fixed precision used in the precision matrix of a
  double high_prec = exp(15);

  // The total number of hyperparameters
  int n_theta_total = 6;

  // ==========================================================
  // Assert that all the input variables look as they shold
  // ==========================================================

  // Number of observations
  assert(!strcmp(data->ints[0]->name, "n"));
  int n = data->ints[0]->ints[0];
  assert(n > 0);

  // Which s0 do the different replications come from?
  assert(!strcmp(data->ints[2]->name, "s0_index"));
  int n_repl = data->ints[2]->len;
  int *s0_index = data->ints[2]->ints;

  // Should we transform the prior?
  assert(!strcmp(data->ints[3]->name, "transform_prior"));
  int transform_prior = data->ints[3]->ints[0];

  // B0, B1 and B2, required for building the SPDE precision
  assert(!strcmp(data->mats[0]->name, "B0"));
  int n_spde = data->mats[0]->nrow;
  int n_spde_total = n_spde * n_repl;
  assert(data->mats[0]->ncol == 3);
  assert(!strcmp(data->mats[1]->name, "B1"));
  assert(data->mats[1]->nrow == n_spde);
  assert(data->mats[1]->ncol == 3);
  assert(!strcmp(data->mats[2]->name, "B2"));
  assert(data->mats[2]->nrow == n_spde);
  assert(data->mats[2]->ncol == 3);

  // Distance to all conditioning sites from all mesh nodes
  assert(!strcmp(data->mats[3]->name, "dist_to_s0_from_mesh"));
  inla_cgeneric_mat_tp *dist_to_s0_from_mesh = data->mats[3];

  inla_cgeneric_mat_tp * C;
  if (transform_prior) {
    assert(!strcmp(data->mats[4]->name, "C"));
    C = data->mats[4];
  }

  // M0, M1 and M2, required for building the SPDE precision
  assert(!strcmp(data->smats[0]->name, "M0"));
  assert(data->smats[0]->nrow == n_spde);
  assert(data->smats[0]->ncol == n_spde);
  assert(!strcmp(data->smats[1]->name, "M1"));
  assert(data->smats[1]->nrow == n_spde);
  assert(data->smats[1]->ncol == n_spde);
  assert(!strcmp(data->smats[2]->name, "M2"));
  assert(data->smats[2]->nrow == n_spde);
  assert(data->smats[2]->ncol == n_spde);

  // y0
  assert(!strcmp(data->doubles[0]->name, "y0"));
  double *y0 = data->doubles[0]->doubles;
  int n_y0 = data->doubles[0]->len;
  assert(n == n_spde_total + 2 * n_y0);

  // dist_to_s0 from all observations locations
  assert(!strcmp(data->doubles[1]->name, "dist_to_s0"));
  assert(data->doubles[1]->len == n_y0);
  double *dist_to_s0 = data->doubles[1]->doubles;

  // Initial values and priors
  assert(!strcmp(data->doubles[2]->name, "init"));
  double *init = data->doubles[2]->doubles;
  int n_theta = data->doubles[2]->len;
  assert(!strcmp(data->doubles[3]->name, "lambda_prior"));
  double *lambda_prior = data->doubles[3]->doubles;
  assert(!strcmp(data->doubles[4]->name, "kappa_prior"));
  double *kappa_prior = data->doubles[4]->doubles;
  assert(!strcmp(data->doubles[5]->name, "rho_prior"));
  double *rho_prior = data->doubles[5]->doubles;
  assert(!strcmp(data->doubles[6]->name, "sigma_prior"));
  double *sigma_prior = data->doubles[6]->doubles;
  assert(!strcmp(data->doubles[7]->name, "rho_b_prior"));
  double *rho_b_prior = data->doubles[7]->doubles;
  assert(!strcmp(data->doubles[8]->name, "tau_prior"));
  double *tau_prior = data->doubles[8]->doubles;

  double log_C_det;
  double * theta_star;
  if (transform_prior) {
    assert(!strcmp(data->doubles[9]->name, "theta_star"));
    assert(n_theta == data->doubles[9]->len);
    theta_star = data->doubles[9]->doubles;
    assert(!strcmp(data->doubles[10]->name, "log_C_det"));
    log_C_det = data->doubles[10]->doubles[0];
  }

  // ==========================================================
  // Initiate the correct parameter values, depending on whether
  // the parameters are fixed or estimated.
  // ==========================================================
  double theta_full[n_theta_total];
  if (theta) {
    int count_theta = 0;
    if (lambda_prior[0] == 0) {
      theta_full[0] = lambda_prior[1];
    } else {
      theta_full[0] = theta[count_theta];
      ++count_theta;
    }
    if (kappa_prior[0] == 0) {
      theta_full[1] = kappa_prior[1];
    } else {
      theta_full[1] = theta[count_theta];
      ++count_theta;
    }
    if (rho_prior[0] == 0) {
      theta_full[2] = rho_prior[1];
    } else {
      theta_full[2] = theta[count_theta];
      ++count_theta;
    }
    if (sigma_prior[0] == 0) {
      theta_full[3] = sigma_prior[1];
    } else {
      theta_full[3] = theta[count_theta];
      ++count_theta;
    }
    if (rho_b_prior[0] == 0) {
      theta_full[4] = rho_b_prior[1];
    } else {
      theta_full[4] = theta[count_theta];
      ++count_theta;
    }
    if (tau_prior[0] == 0) {
      theta_full[5] = tau_prior[1];
    } else {
      theta_full[5] = theta[count_theta];
      ++count_theta;
    }
    assert(count_theta == n_theta);
  } else {
    for (int i = 0; i < n_theta_total; ++i) {
      theta_full[i] = NAN;
    }
  }
 
  // ==========================================================
  // This switch statement is the required method for implementing
  // cgeneric models in R-INLA
  // ==========================================================
  switch (cmd) {
  case INLA_CGENERIC_VOID:
    {
      assert(!(cmd == INLA_CGENERIC_VOID));
      break;
    }

  case INLA_CGENERIC_GRAPH:
    {
      // return a vector of indices with format
      // c(n, M, ii, jj)
      // where ii<=jj, ii is non-decreasing and jj is non-decreasing within each ii,
      // and M is the length of ii

      // Compute the precision matrix of the SPDE part, to locate non-zero elements
      inla_cgeneric_smat_tp precision = spde_precision(theta_full[2], theta_full[3], data->mats, data->smats);

      // Only keep the upper diagonal of the precision matrix
      upper_diag(&precision);

      // Sort the entries of the smat so they correspond to the R-INLA requirements
      sort_smat(&precision);

      // Replicate the precision matrix several times to create a block diagonal matrix
      block_diag_smat(&precision, n_repl);

      // Create the graph for the a part and the nugget part
      diag_mat_tp diagonal_graph_part = diag(1, 2 * n_y0);

      // Merge the two graph matrices together
      smat_diag_merge(&precision, &diagonal_graph_part);

      // Extract all the necessary i and j values from the block matrix
      ret = Malloc(2 + 2 * precision.n, double);
      ret[0] = n;
      ret[1] = precision.n;
      for (int i = 0; i < precision.n; i++) {
	ret[2 + i] = precision.i[i];                                  // i
	ret[2 + precision.n + i] = precision.j[i];                    // j
      }

      // Free unused matrices so we don't get memory leaks
      free_smat(&precision);
      free_diag(&diagonal_graph_part);

      break;
    }

  case INLA_CGENERIC_Q:
    {

      // Compute the precision matrix
      inla_cgeneric_smat_tp precision = spde_precision(theta_full[2], theta_full[3], data->mats, data->smats);

      // Only keep the upper diagonal of the precision matrix
      upper_diag(&precision);

      // Sort the entries of the smat so they correspond to the R-INLA requirements
      sort_smat(&precision);

      // Replicate the precision matrix several times to create a block diagonal matrix
      block_diag_smat(&precision, n_repl);

      // Compute the values of b for each distance from a conditioning site
      inla_cgeneric_mat_tp b_vals = copy_mat(dist_to_s0_from_mesh);
      double tmp;
      double rho_b = exp(theta_full[4]);
      for (int i = 0; i < b_vals.nrow * b_vals.ncol; i++) {
	// Ensure that we don't get numerical problems where b_vals.x[i] becomes INF or NAN
	tmp = dist_to_s0_from_mesh->x[i] / rho_b;
	if (tmp < 0.0000000001) {
	  tmp = 0.0000000001;
	}
	b_vals.x[i] = 1 / sqrt(1 - exp(-4 * tmp));
      }

      // Create a diagonal matrix containing the values of 1 / b, that can be used
      // for rescaling the precision matrix with the correct b values
      diag_mat_tp b_inv;
      b_inv.x = Malloc(precision.nrow, double);
      b_inv.dim = precision.nrow;
      for (int i = 0; i < n_repl; i++) {
	for (int j = 0; j < n_spde; j++) {
	  b_inv.x[j + i * n_spde] = b_vals.x[j + (s0_index[i]) * n_spde];
	}
      }
      // Multiply the diagonal matrix with the precision matrix: B_inv * Q * B_inv
      diag_smat_diag_mult(&precision, &b_inv);

      // Create the precision matrix for the a() part
      diag_mat_tp diag_precision = diag(high_prec, n_y0);

      // Add the a diagonal precision to the precision matrix
      smat_diag_merge(&precision, &diag_precision);

      // Change the diag_precision to that for the nugget effect
      double tau = exp(theta_full[5]);
      for (int i = 0; i < diag_precision.dim; ++i) {
	diag_precision.x[i] = tau;
      }

      // Add the nugget diagonal precision to the precision matrix
      smat_diag_merge(&precision, &diag_precision);

      // optimized format
      // return c(n, M, Qij) in the same order as defined in INLA_CGENERIC_GRAPH
      // where M is the length of Qij
      ret = Malloc(2 + precision.n, double);
      ret[0] = -1;                                   // code for optimized output
      ret[1] = precision.n;                          // number of (i <= j)
      for (int i = 0; i < precision.n; i++) {
	ret[2 + i] = precision.x[i];
      }

      // Avoid memory leaks
      free(b_inv.x);
      free_mat(&b_vals);
      free_smat(&precision);
      free_diag(&diag_precision);

      break;
    }

  case INLA_CGENERIC_MU:
    {
      // return (N, mu)
      // if N==0 then mu is not needed as its taken to be mu[]==0
      double lambda = exp(theta_full[0]);
      double kappa = exp(theta_full[1]);

      ret = Calloc(1 + n_spde_total + 2 * n_y0, double);
      ret[0] = n_spde_total + 2 * n_y0;
      //for (int i = 1; i <= n_spde_total; ++i) {
      //	ret[i] = 0;
      //}
      for (int i = 0; i < n_y0; ++i) {
	ret[n_spde_total + 1 + i] = y0[i] * exp(-pow(dist_to_s0[i] / lambda, kappa));
      }
      //for (int i = n_spde_total + n_y0 + 1; i <= ret[1]; ++i) {
      //	ret[i] = 0;
      //}
      break;
    }

  case INLA_CGENERIC_INITIAL:
    {
      // return c(M, initials)
      // where M is the number of hyperparameters

      // The initial values depend on how many parameters that are fixed or estimated.
      ret = Malloc(n_theta + 1, double);
      ret[0] = n_theta;
      for (int i = 0; i < n_theta; i++) {
	ret[i + 1] = init[i];
      }
      break;
    }

  case INLA_CGENERIC_LOG_NORM_CONST:
    {
      // return c(NORM_CONST) or a NULL-pointer if INLA should compute it by itself

      ret = NULL;
      break;
    }

  case INLA_CGENERIC_LOG_PRIOR:
    {
      // return c(LOG_PRIOR)
      double * theta_free = Malloc(n_theta, double);
      for (int i = 0; i < n_theta; ++i) {
	theta_free[i] = theta[i];
      }
      if (transform_prior) {
	for (int i = 0; i < n_theta; ++i) {
	  theta_free[i] -= theta_star[i];
	}
	mat_vec_mult_in_place(C, theta_free);
	for (int i = 0; i < n_theta; ++i) {
	  theta_free[i] += theta_star[i];
	}
      }

      // Only add terms for the parameters that are not fixed
      ret = Malloc(1, double);
      ret[0] = 0;
      int counter = 0;
      if (lambda_prior[0] != 0) {
      	ret[0] += gaussian_const - log(lambda_prior[2]) -
	  pow(theta_free[counter] - lambda_prior[1], 2) / (2 * pow(lambda_prior[2], 2));
	++counter;
      }
      if (kappa_prior[0] != 0) {
      	ret[0] += gaussian_const - log(kappa_prior[2]) -
	  pow(theta_free[counter] - kappa_prior[1], 2) / (2 * pow(kappa_prior[2], 2));
	++counter;
      }
      if (rho_prior[0] != 0) {
	double lambda0 = -log(rho_prior[2]) * rho_prior[1];
	ret[0] += log(lambda0) - lambda0 * exp(-theta_free[counter]) - theta_free[counter];
	++counter;
      }
      if (sigma_prior[0] != 0) {
	double lambda1 = -log(sigma_prior[2]) / sigma_prior[1];
	ret[0] += log(lambda1) - lambda1 * exp(theta_free[counter]) + theta_free[counter];
	++counter;
      }
      if (rho_b_prior[0] != 0) {
	ret[0] += gaussian_const - log(rho_b_prior[2]) -
	  pow(theta_free[counter] - rho_b_prior[1], 2) / (2 * pow(rho_b_prior[2], 2));
	++counter;
      }
      if (tau_prior[0] != 0) {
	ret[0] += tau_prior[1] * log(tau_prior[2]) - lgamma(tau_prior[1]) +
	  (tau_prior[1] - 1) * theta_free[counter] -
	  tau_prior[2] * exp(theta_free[counter]);
      }
      if (transform_prior) {
	ret[0] += log_C_det;
      }
      break;
    }

  case INLA_CGENERIC_QUIT:
  default:
    break;
  }

  return (ret);
}

*/

