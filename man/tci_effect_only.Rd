% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tci_algorithms.R
\name{tci_effect_only}
\alias{tci_effect_only}
\title{TCI algorithm for effect-site targeting}
\usage{
tci_effect_only(
  Ct,
  pkmod,
  dtm = 1/6,
  tmax_search = 10,
  maxrt = 1200,
  grid_len = 1200,
  ...
)
}
\arguments{
\item{Ct}{Numeric vector of target effect-site concentrations.}

\item{pkmod}{PK model}

\item{dtm}{Frequency of TCI updates. Default is 1/6 minutes = 10 seconds.}

\item{tmax_search}{Outer bound on times searched to find a maximum concentration
following an infusion of duration dtm. Defaults to 20 minutes. May need to be increased
if a drug has a slow elimination rate.}

\item{maxrt}{Maximum infusion rate of TCI pump. Defaults to 1200.}

\item{grid_len}{Number of time points used to identify time of maximum concentration.
Can be increased for more precision.}

\item{...}{List or vector of named values to be passed on to update.pkmod.}
}
\value{
Numeric value
}
\description{
Function for calculating a TCI infusion schedule corresponding to a set of target concentrations.
This function makes use of formulas described by Shafer and Gregg (1992) in "Algorithms to rapidly achieve
and maintain stable drug concentrations at the site of drug effect with a computer-controlled infusion pump"
}
\examples{
# three compartment model with effect site
my_mod <- pkmod(pars_pk = c(v1 = 8.995, v2 = 17.297, v3 = 120.963,
cl = 1.382, q2 = 0.919, q3 = 0.609, ke0 = 1.289))
tci_effect_only(Ct = 2, pkmod = my_mod)
# update parameters
tci_effect_only(Ct = 2, pkmod = my_mod, pars_pk = c(v1 = 12, cl = 2))
}
