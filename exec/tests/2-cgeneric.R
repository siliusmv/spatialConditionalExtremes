devtools::load_all()
library(INLA)
library(Matrix)

# In this script we evaluate the correctness of the cgeneric functions used for
# computing the precision matrix of the SPDE approximation.
# This is performed by creating a small mesh and an inla.spde2 object using that mesh.
# The precision matrix of the SPDE approximation is then computed an compared using both the
# cgeneric functions and using the inla.spde2.precision() function.

# Compile and link all the cgeneric test scripts, if this has not already been done
make_cgeneric("test")

# Create a small grid of locations
loc = as.matrix(expand.grid(x = 1:3, y = 1:2))

# Create a mesh on the grid
mesh = inla.mesh.2d(
  loc = loc,
  boundary = list(inla.nonconvex.hull(loc, convex = -.3)),
  cutoff = 1,
  max.edge = 4)
#plot(mesh)
#points(loc)
mesh$n

# Create the SPDE on the mesh
spde = inla.spde2.pcmatern(
  mesh,
  prior.range = c(1, .5),
  prior.sigma = c(1, .5))

# ====================================================================
# Save the B- and M-matrices from the SPDE object in a file format
# that is readable for our c code, so we can compare precision matrices
# made by our c functions with matrices made using inla.spde2.precision()
# ====================================================================

# Functions for saving the matrices
save_mat = function(m, path) {
  file.create(path)
  cat("mat", nrow(m), ncol(m), sep = "\n", file = path)
  cat(t(m), sep = "\n", file = path, append = TRUE)
}
save_smat = function(m, path) {
  file.create(path)
  n = length(m@x)
  cat("smat", n, nrow(m), ncol(m), sep = "\n", file = path)
  for (i in seq_len(n)) {
    cat(m@i[i], m@j[i], m@x[i], sep = ";", file = path, append = TRUE)
    cat("\n", file = path, append = TRUE)
  }
}

# Create the directory for saving the matrices
path = file.path(here::here(), "cgeneric-data")
dir.create(path, recursive = TRUE)

for (i in 0:2) {
  # Save all the B-mats
  filename = file.path(path, paste0("B", i, ".txt"))
  save_mat(spde$param.inla[[paste0("B", i)]], filename)
  # Save all the M-smats
  filename = file.path(path, paste0("M", i, ".txt"))
  save_smat(spde$param.inla[[paste0("M", i)]], filename)
}

# ============================================================
# Compute Q using c and using inla.spde2.precision(), and
# compare the results
# ============================================================

# Function for calling an executable from cgeneric_dir()
execute_c_script = function(name, ...) {
  current_path = getwd()
  on.exit(setwd(current_path))
  setwd(cgeneric_dir())
  execute_shell_script(name, ...)
}

# Compare R and C by computing and printing the precision matrix
# of the SPDE with parameters log_rho and log_sigma.
compare_r_and_c = function(log_rho, log_sigma) {
  Q = inla.spde2.precision(spde, c(log_rho, log_sigma))
  execute_c_script("./test1.o", c(log_rho, log_sigma))
  round(as.matrix(Q), digits = 6)
}

# Further test the implemented c functions for extracting the upper diagonal
# part of a matrix, sorting the elements of a sparse matrix,
# and creating a diagonal block matrix with the same matrix replicated multiple times
test_c_further = function(log_rho, log_sigma) {
  execute_c_script("./test2.o", c(log_rho, log_sigma))
}

# Draw some random values for log_rho and log_sigma
log_rho = runif(1, 0, 2)
log_sigma = runif(1, 0, 2)

# Compare R and C
compare_r_and_c(log_rho, log_sigma)
# Everything looks good!

# Test C further
test_c_further(log_rho, log_sigma)
# Everything looks good!

# Clean up after your self by removing the directory that contains
# textfiles of the matrices of the SPDE
unlink(path, recursive = TRUE)
