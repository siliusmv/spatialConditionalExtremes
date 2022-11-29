// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <cassert>

using namespace Rcpp;

static double const log2pi = std::log(2.0 * M_PI);

// Perform fast forward substitution for solving the equation
// Lx = b for x, where L is lower triangular
void forward_substitution_fast(arma::vec &res,
			       arma::mat const &L,
			       arma::vec const &b) {
  arma::uword const n = b.size();
  for (arma::uword i = 0; i < n; ++i) {
    res.at(i) = b.at(i);
    for (arma::uword j = 0; j < i; ++j) {
      res.at(i) -= L.at(i, j) * res.at(j);
    }
    res.at(i) /= L.at(i, i);
  }
}

// Compute dmvnorm quickly when x is not a matrix, but a vector
double dmvnorm_double_fast(arma::vec const &x,  
			   arma::mat const &sigma) {
  arma::uword const d = x.n_rows;
  arma::mat const chol = arma::chol(sigma, "lower");
  if (chol.n_cols != d || chol.n_rows != d) {
    stop("Dimensions of x and chol do not correspond");
  }
  double out =
    - arma::sum(log(chol.diag())) - (double)d * 0.5 * log2pi;
    
  arma::vec z(d);
  forward_substitution_fast(z, chol, x);
  out -= 0.5 * arma::dot(z, z);     

  return out;
}

// Compute Σ = A * B * Σ0 * B * A', where B is a diagonal
// matrix containing the values of b.
// This function is fast and unsafe, so if the length of b is different
// than the number of rows or columns of sigma, we get undefined behaviour.
arma::mat sigma_func(arma::sp_mat const &A,
		     arma::vec const &b,
		     arma::mat const &sigma,
		     double const nugget) {
  arma::mat res = sigma;
  arma::uword const d = b.size();
  for (arma::uword i = 0; i < d; ++i) {
    for (arma::uword j = 0; j < d; ++j) {
      res.at(j, i) *= b.at(i) * b.at(j);
    }
  }
  res = A * res * A.t();
  res.diag() += nugget;
  return res;
}

//' Compute the (log-)likelihood of the conditional extremes distribution when
//' a and b are given. We assume that a has already been subtracted from x,
//' i.e. x = (y - a) for some observations y.
//' Further, we assume that b does not depend on y0, meaning that we only need one
//' value of b for each distance, instead of a matrix B that has one value for each
//' pair of distance and threshold exceedance y0.
//' 
//' The input variables are:
//' x: an (n x d)-dimensional matrix of observations where a has already been subtracted.
//' A: The projection matrix used for building the SPDE approximation
//' b: An m-dimensional vector containing the values of b at the m triangular mesh nodes.
//' sigma0: The covariance matrix of the m Gaussian random variables in the triangular mesh.
//' nugget: The variance of the nugget effect.
//' logd: Boolean explaining if we should return the log-likelihood or the likelihood.
//' na_rm: Boolean explaining how we should treat NA values. If na_rm = false, then
//'   any column of x that returns an NA value will result in an NA value in the output.
//'   If na_rm = true, then we remove the NA values before computing the likelihood.
//'   So if e.g. a column has 3 NA variables, then we remove these, and then we compute
//'   the likelihood for a (d-3)-dimensional Gaussian random variable.
// [[Rcpp::export]]
arma::vec dconditional_no_beta(arma::mat const &x,  
			       arma::sp_mat const &A,
			       arma::vec const &b,
			       arma::mat const &sigma0,
			       double const nugget,
			       bool const logd = true,
			       bool const na_rm = true) { 
  // Check that the dimensions of the input are correct
  arma::uword const n = x.n_cols;
  arma::uword const d = x.n_rows;
  arma::uword const m = b.size();
  if (sigma0.n_cols != m || sigma0.n_rows != m) {
    stop("The dimensions of sigma0 and b do not agree with each other");
  }

  // Allocate the result
  arma::vec out(n);

  // Compute the covariance matrix of the Gaussian random field
  arma::mat sigma = sigma_func(A, b, sigma0, nugget);

  // Compute the cholesky factorisation of sigma, and the constant of the
  // Gaussian log-likelihood
  arma::mat const chol = arma::chol(sigma, "lower");
  double const lpdf_const =
    - arma::sum(log(chol.diag())) - (double)d * 0.5 * log2pi;

  // Allocate some variables used for computing the log-likelihood
  arma::uword count;
  arma::uvec good_index(d);
  arma::vec z(d);

  // Loop over all the columns of x and compute log-likelihood terms
  for (arma::uword i = 0; i < n; ++i) {
    z = x.col(i);

    if (na_rm) {
      // Count all the elements of z that are not NA
      count = 0;
      for (arma::uword j = 0; j < d; ++j) {
	if (!NumericVector::is_na(z(j))) {
	  good_index(count) = j;
	  ++count;
	}
      }

      if (count != d) {
	// If z contains NA variables, we need to remove these and compute the log-likelihood
	// of the non-NA variables. This is computationally demanding as it requires the
	// computation of another Cholesky factorisation
	out(i) = dmvnorm_double_fast(z(good_index.head(count)),
				     sigma(good_index.head(count), good_index.head(count)));
      }
    }

    if (!na_rm || count == d) {
      // If we don't care about removing NA variables, or z does not contain any NA variables,
      // we can use the already computed Cholesky factorisation to compute the log-likelihood of z
      forward_substitution_fast(z, chol, x.col(i));
      out(i) = lpdf_const - 0.5 * arma::dot(z, z);     
    }
  }
      
  if (logd) {
    return out;
  } else {
    return exp(out);
  }
}

//' Compute the (log-)likelihood of the conditional extremes distribution when
//' a and b are given. We assume that a has already been subtracted from x,
//' i.e. x = (y - a) for some observations y.
//' 
//' The input variables are:
//' x: an (n x d)-dimensional matrix of observations where a has already been subtracted.
//' A: The projection matrix used for building the SPDE approximation
//' B: An (m x n)-dimensional matrix containing the values of b at the m triangular mesh nodes,
//'   for each of the n threshold exceedances at the conditioning sites of interest
//' sigma0: The covariance matrix of the m Gaussian random variables in the triangular mesh.
//' nugget: The variance of the nugget effect.
//' logd: Boolean explaining if we should return the log-likelihood or the likelihood.
//' na_rm: Boolean explaining how we should treat NA values. If na_rm = false, then
//'   any column of x that returns an NA value will result in an NA value in the output.
//'   If na_rm = true, then we remove the NA values before computing the likelihood.
//'   So if e.g. a column has 3 NA variables, then we remove these, and then we compute
//'   the likelihood for a (d-3)-dimensional Gaussian random variable.
// [[Rcpp::export]]
arma::vec dconditional(arma::mat const &x,  
		       arma::sp_mat const &A,
		       arma::mat const &B,
		       arma::mat const &sigma0,
		       double const nugget,
		       bool const logd = true,
		       bool const na_rm = true) { 
  // Check that the dimensions of the input are correct
  arma::uword const n = x.n_cols;
  arma::uword const d = x.n_rows;
  arma::uword const m = B.n_rows;
  if (sigma0.n_cols != m || sigma0.n_rows != m) {
    stop("The dimensions of sigma0 and B do not agree with each other");
  }
  if (B.n_cols != n) stop("The dimensions of B and x do not agree with each other");

  // Allocate the result
  arma::vec out(n);

  // Allocate some variables used for computing the log-likelihood
  arma::uword count;
  arma::uvec good_index(d);
  arma::mat sigma(d, d);
  arma::vec z(d);

  // Loop over all the columns of x and compute log-likelihood terms
  for (arma::uword i = 0; i < n; ++i) {
    z = x.col(i);

    // Compute sigma using the ith column of B
    sigma = sigma_func(A, B.col(i), sigma0, nugget);

    if (na_rm) {
      // Count all the elements of z that are not NA
      count = 0;
      for (arma::uword j = 0; j < d; ++j) {
	if (!NumericVector::is_na(z(j))) {
	  good_index(count) = j;
	  ++count;
	}
      }

      if (count != d) {
	// If z contains NA variables, we need to remove these and compute the log-likelihood
	// of the non-NA variables. This is computationally demanding as it requires the
	// computation of another Cholesky factorisation
	out(i) = dmvnorm_double_fast(z(good_index.head(count)),
				     sigma(good_index.head(count), good_index.head(count)));
      }
    }

    if (!na_rm || count == d) {
      // If we don't care about removing NA variables, or z does not contain any NA variables,
      // we can use the already computed Cholesky factorisation to compute the log-likelihood of z
      out(i) = dmvnorm_double_fast(z, sigma);
    }
  }
      
  if (logd) {
    return out;
  } else {
    return exp(out);
  }
}
