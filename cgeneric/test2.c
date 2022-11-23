#include <assert.h>
#if !defined(__FreeBSD__)
#endif
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>

#include "cgeneric.h"
#include "smat-operations.h"
#include "spde-precision.h"
#include "test-funcs.h"

#define Calloc(n, type)  (type *) calloc(n, sizeof(type))

int main(int argc, char * argv[]) {

  assert(argc == 3);
  double log_rho = atof(argv[1]);
  double log_sigma = atof(argv[2]);

  inla_cgeneric_data_tp data = read_cgeneric_data_from_dir("../cgeneric-data");

  inla_cgeneric_smat_tp Q = spde_precision(log_rho, log_sigma, data.mats, data.smats, -1);
  Q.name = Calloc(1, char);
  strcpy(Q.name, "Q");

  printf("Upper diagonal part of the precision matrix:\n");
  upper_diag(&Q);
  print_smat_as_mat(Q);

  printf("Upper diagonal part of the precision matrix, given in a sparse matrix format\n");
  printf("where all the elements are sorted such that i is nondecreasing and");
  printf("j is nondecreasing within each value of i:\n");
  sort_smat(&Q);
  print_smat(Q);

  printf("Block diagonal matrix that contains two copies of the upper diagonal of the precision matrix:\n");
  block_diag_smat(&Q, 2);
  print_smat_as_mat(Q);

  free_smat(&Q);

  return 0;
}
