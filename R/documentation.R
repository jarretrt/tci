#' tci package documentation
#' @title tci_documentation
#' @name tci_documentation
#' @description Functions to implement target-controlled infusion algorithms
#' @details This package contains functions to implement target-controlled
#' infusion (TCI) algorithms for compartmental PK models under intravenous administration.
#' TCI algorithms for plasma or effect-site targeting are included and can be
#' extended to pharmacodynamic responses. Custom PK-PD models and custom TCI
#' algorithms can be specified. Functions are provided to simulate responses from
#' PK/PK-PD models under open- or closed-loop control.

## usethis namespace: start
#' @useDynLib tci, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
