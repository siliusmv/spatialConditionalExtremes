#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "cgeneric.h"

#define Calloc(n, type)  (type *) calloc(n, sizeof(type))
#define Malloc(n, type)  (type *) malloc(n * sizeof(type))

double * a_model(inla_cgeneric_cmd_tp cmd,
		 double * theta,
		 inla_cgeneric_data_tp * data) {
  double * ret = NULL;

  // This is an approximation for -0.5 * log(2 * pi)
  double gaussian_const = -0.91893853320467;
  // The fixed precision used in the precision matrix of a
  double high_prec = exp(15);

  // ==========================================================
  // Assert that all the input variables look as they shold
  // ==========================================================
  assert(!strcasecmp(data->doubles[0]->name, "y0"));
  int n = data->doubles[0]->len;
  assert(n == data->ints[0]->ints[0]);
  double *y0 = data->doubles[0]->doubles;

  assert(!strcmp(data->ints[2]->name, "is_fixed"));
  int * is_fixed = data->ints[2]->ints;
  int n_theta = data->ints[2]->len;
  int n_free = 0;
  for (int i = 0; i < n_theta; ++i) {
    if (!is_fixed[i]) ++n_free;
  }

  assert(!strcasecmp(data->doubles[1]->name, "dist_to_s0"));
  assert(data->doubles[1]->len == n);
  double *dist = data->doubles[1]->doubles;

  // Initial values for the hyperparameters
  assert(!strcasecmp(data->doubles[2]->name, "init"));
  double *init = data->doubles[2]->doubles;

  int prior_start = 3; // This is the first element of data->doubles that contains a prior

  // ==========================================================
  // Initiate the correct parameter values, depending on whether
  // the parameters are fixed or estimated.
  // ==========================================================
  double lambda, kappa;
  assert(n_theta == 2);
  if (theta) {
    double theta_full[2];
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
  } else {
    lambda = kappa = NAN;
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
      // The precision matrix is just a diagonal matrix with each element equal to high_prec
      ret = Malloc(2 + 2 * n, double);
      ret[0] = n;
      ret[1] = n;
      for (int i = 0; i < n; i++) {
      	ret[2 + i] = i;			       /* i */
      	ret[2 + n + i] = i;		       /* j */
      }
      break;
    }

  case INLA_CGENERIC_Q:
    {
      // optimized format
      // return c(n, M, Qij) in the same order as defined in INLA_CGENERIC_GRAPH
      // where M is the length of Qij

      // The precision matrix is just a diagonal matrix with each element equal to high_prec
      ret = Malloc(2 + n, double);
      ret[0] = -1;			       /* code for optimized output */
      ret[1] = n;			       /* number of (i <= j) */
      for (int i = 0; i < n; i++) {
	ret[2 + i] = high_prec;
      }
      break;
    }

  case INLA_CGENERIC_MU:
    {
      // return (n, mu)
      // if n==0 then mu is not needed as its taken to be mu[]==0
      ret = Malloc(1 + n, double);
      ret[0] = n;
      for (int i = 0; i < n; i++) {
	ret[1 + i] = y0[i] * exp(-pow(dist[i] / lambda, kappa));
      }
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
      ret = Malloc(1, double);
      ret[0] = n * (gaussian_const + 0.5 * log(high_prec));
      break;
    }

  case INLA_CGENERIC_LOG_PRIOR:
    {
      // return c(LOG_PRIOR)

      // Only add terms for the parameters that are not fixed
      ret = Calloc(1, double);
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

    }
  case INLA_CGENERIC_QUIT:
  default:
    break;
  }

  return ret;
}
