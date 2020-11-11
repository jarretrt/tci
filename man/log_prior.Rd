% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulations.R
\name{log_prior}
\alias{log_prior}
\title{Function to return the prior probability for a set of parameters assuming they are log-normally distributed.
It is assumed that the last value of of lpr is the prior mean for the variance parameter.}
\usage{
log_prior(lpr, mu, sig)
}
\arguments{
\item{lpr}{Logged parameter values to evaluate prior distribution at.}

\item{mu}{Mean of multivariate normal prior distribution.}

\item{sig}{Variance-covariance matrix of multivariate normal prior distribution.}
}
\description{
Function to return the prior probability for a set of parameters assuming they are log-normally distributed.
It is assumed that the last value of of lpr is the prior mean for the variance parameter.
}