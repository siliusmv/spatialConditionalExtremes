% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conditional-likelihood.R
\name{dconditional_r_no_beta}
\alias{dconditional_r_no_beta}
\title{Compute the (log-)likelihood of the conditional extremes distribution when
a and b are given. We assume that a has already been subtracted from x,
i.e. x = (y - a) for some observations y.
Further, we assume that b does not depend on y0, meaning that we only need one
value of b for each distance, instead of a matrix B that has one value for each
pair of distance and threshold exceedance y0.
This function is mainly written to evaluate the speed and correctness of the Rcpp
function dconditional_no_beta()}
\usage{
dconditional_r_no_beta(x, A, b, sigma0, nugget, logd = TRUE, na_rm = TRUE)
}
\description{
The input variables are:
x: an (n x d)-dimensional matrix of observations where a has already been subtracted.
A: The projection matrix used for building the SPDE approximation
b: An m-dimensional vector containing the values of b at the m triangular mesh nodes.
sigma0: The covariance matrix of the m Gaussian random variables in the triangular mesh.
nugget: The variance of the nugget effect.
logd: Boolean explaining if we should return the log-likelihood or the likelihood.
na_rm: Boolean explaining how we should treat NA values. If na_rm = false, then
any column of x that returns an NA value will result in an NA value in the output.
If na_rm = true, then we remove the NA values before computing the likelihood.
So if e.g. a column has 3 NA variables, then we remove these, and then we compute
the likelihood for a (d-3)-dimensional Gaussian random variable.
}
