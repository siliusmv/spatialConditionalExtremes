% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-processing.R
\name{empirical_chi}
\alias{empirical_chi}
\title{Compute the empirical estimators of χ_p(d) for a set of
threshold values and distances.
Function input:
\itemize{
\item data is a list of data on the same format as that created by
the function extract_extreme_fields().
\item thresholds is a vector of threshold values used for estimating
the extremal correlation coefficient.
\item unique_dist is an integer vector of distances d for which we
should estimate χ_p(d)
}}
\usage{
empirical_chi(data, thresholds, unique_dist)
}
\description{
Compute the empirical estimators of χ_p(d) for a set of
threshold values and distances.
Function input:
\itemize{
\item data is a list of data on the same format as that created by
the function extract_extreme_fields().
\item thresholds is a vector of threshold values used for estimating
the extremal correlation coefficient.
\item unique_dist is an integer vector of distances d for which we
should estimate χ_p(d)
}
}
