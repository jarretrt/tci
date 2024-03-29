% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tci_algorithms.R
\name{inf_tci_pop}
\alias{inf_tci_pop}
\title{Calculate infusion rates to reach target values within a population}
\usage{
inf_tci_pop(
  mod,
  target_vals,
  target_tms,
  pop_fn = median,
  dtm = 1/6,
  inf_duration = NULL,
  cmpt = NULL,
  inittm = 0,
  maxrt = 1e+05
)
}
\arguments{
\item{mod}{Object with class `poppkmod`, created by `poppkmod()`.}

\item{target_vals}{A vector of numeric values indicating PK or PD targets for TCI algorithm.}

\item{target_tms}{A vector of numeric values indicating times at which the TCI algorithm should
begin targeting each value.}

\item{pop_fn}{Function to apply to distribution of response values. Defaults to setting
median value equal to targets.}

\item{dtm}{TCI update frequency. Defaults to 1/6, corresponding to 10-second
intervals if model parameters are in terms of minutes.}

\item{inf_duration}{Optional parameter to describe duration of infusions. Can be
less than or equal to `dtm`. Defaults to `dtm` if unspecified.}

\item{cmpt}{Compartment used for targeting. Defaults to PD response if applicable and
effect-site compartment if not.}

\item{inittm}{Initial time to start TCI algorithm. Cannot be greater than
the minimum value of `target_tms`.}

\item{maxrt}{Maximum infusion rate.}
}
\description{
This function will calculate infusion rates required to reach targets at a specified
position (e.g., median, 25th percentile) of the response distribution in a population.
The user supplies a `poppkmod` object, a set of values to be obtained, a set of
times at which to reach them, and a function to apply to the distribution of responses.
}
\examples{
nid = 100
data <- data.frame(ID = 1:nid, AGE = sample(10:70, nid, replace = TRUE),
TBW = sample(40:90, nid, TRUE), HGT = sample(130:210, nid, TRUE),
MALE = sample(c(TRUE,FALSE),nid,TRUE))
elvd_mod <- poppkmod(data, drug = "ppf", model = "eleveld")
# calculate infusions to keep 90\% of patients below 70
inf <- inf_tci_pop(elvd_mod, target_vals = c(70,70), target_tms = c(0,50),
pop_fn = function(...) quantile(..., 0.9), dtm = 10)
tms <- seq(0,50,0.1)
resp <- as.data.frame(predict(elvd_mod, inf, tms))
ggplot(resp, aes(x = time, y = pdresp, group = id)) +
geom_line(alpha = 0.1) +
geom_hline(yintercept = 70)
}
