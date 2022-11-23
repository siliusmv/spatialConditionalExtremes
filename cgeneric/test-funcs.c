#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <dirent.h>

#include "cgeneric.h"
#include "test-funcs.h"

inla_cgeneric_data_tp read_cgeneric_data_from_dir(const char *dir) {

  // Preallocate the result
  inla_cgeneric_data_tp res;
  res.n_mats = 3;
  res.n_smats = 3;
  res.n_ints = 0;
  res.n_doubles = 0;
  res.n_chars = 0;

  // Allocate the correct number of mats/smats in res
  res.mats = Calloc(res.n_mats, inla_cgeneric_mat_tp *);
  res.smats = Calloc(res.n_smats, inla_cgeneric_smat_tp *);

  char * mat_filenames[3] = {"B0.txt", "B1.txt", "B2.txt"};
  char * smat_filenames[3] = {"M0.txt", "M1.txt", "M2.txt"};
  char full_path[1024];

  // Read all the mats/smats from the files and into the res object
  for (int i = 0; i < 3; ++i) {
    strcpy(full_path, dir);
    strcat(full_path, "/");
    strcat(full_path, mat_filenames[i]);
    res.mats[i] = Calloc(1, inla_cgeneric_mat_tp);
    read_mat_from_file(full_path, res.mats[i]);

    strcpy(full_path, dir);
    strcat(full_path, "/");
    strcat(full_path, smat_filenames[i]);
    res.smats[i] = Calloc(1, inla_cgeneric_smat_tp);
    read_smat_from_file(full_path, res.smats[i]);
  }

  return res;
}

void read_mat_from_file(char const * filename, inla_cgeneric_mat_tp * res) {
  FILE *file = fopen(filename, "r");
  char line[1024], *tmp;
  assert(fgets(line, 1024, file) != NULL);
  assert(!strcmp(line, "mat\n"));
  assert(fgets(line, 1024, file) != NULL);
  res->nrow = atoi(line);
  assert(fgets(line, 1024, file) != NULL);
  res->ncol = atoi(line);
  res->x = Calloc(res->nrow * res->ncol, double);
  for (int i = 0; i < res->nrow * res->ncol; i++) {
    assert(fgets(line, 1024, file) != NULL);
    res->x[i] = strtod(line, &tmp);
  }
  fclose(file);
}

void read_smat_from_file(char const * filename, inla_cgeneric_smat_tp * res) {
  FILE *file = fopen(filename, "r");
  char line[1024];
  char *token, *tmp;
  assert(fgets(line, 1024, file) != NULL);
  assert(!strcmp(line, "smat\n"));
  assert(fgets(line, 1024, file) != NULL);
  res->n = atoi(line);
  assert(fgets(line, 1024, file) != NULL);
  res->nrow = atoi(line);
  assert(fgets(line, 1024, file) != NULL);
  res->ncol = atoi(line);
  res->x = Calloc(res->n, double);
  res->i = Calloc(res->n, int);
  res->j = Calloc(res->n, int);
  for (int i = 0; i < res->n; ++i) {
    assert(fgets(line, 1024, file) != NULL);
    token = strtok(line, ";");
    res->i[i] = atoi(token);
    token = strtok(NULL, ";");
    res->j[i] = atoi(token);
    token = strtok(NULL, ";");
    res->x[i] = strtod(token, &tmp);
  }
  fclose(file);
}

inla_cgeneric_mat_tp smat_to_mat(inla_cgeneric_smat_tp const A) {
  inla_cgeneric_mat_tp res;
  res.nrow = A.nrow;
  res.ncol = A.ncol;
  res.name = A.name;
  res.x = Calloc(res.nrow * res.ncol, double);
  for (int i = 0; i < A.n; i++) {
    res.x[A.j[i] + A.i[i] * res.ncol] = A.x[i];
  }
  return res;
}

void print_smat(inla_cgeneric_smat_tp const A) {
  printf("%s: (n: %d, nrow: %d, ncol: %d):\n", A.name, A.n, A.nrow, A.ncol);
  printf("i j x\n");
  for (int i = 0; i < A.n; i++) {
    printf("%d %d %f\n", A.i[i], A.j[i], A.x[i]);
  }
  printf("\n");
}

void print_mat(inla_cgeneric_mat_tp const A) {
  printf("%s: (nrow: %d, ncol: %d):\n", A.name, A.nrow, A.ncol);
  for (int i = 0; i < A.nrow; i++) {
    for (int j = 0; j < A.ncol; j++) {
      printf("%f ", A.x[j + i * A.ncol]);
    }
    printf("\n");
  }
  printf("\n");
}

void print_smat_as_mat(inla_cgeneric_smat_tp const A) {
  inla_cgeneric_mat_tp tmp = smat_to_mat(A);
  print_mat(tmp);
}
