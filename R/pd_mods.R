# PD functions

#' Emax function. Assumes maximum value of Emx = 100 and maximum change in value of 100 (i.e. minimum value of 0).
#' @param pars Vector of parameters: (c50,gamma).
#' @param ce Effect-site concentration.
#' @param Emx Maximum value of function.
#' @param E0 Maximum difference between highest and lowest function value.
Emax <- function(pars, ce, Emx = 100, E0 = 100) E0 - Emx*(ce^pars[2] / (ce^pars[2] + pars[1]^pars[2]))


#' Emax function with only c50 to be estimated
#' @param pars c50
#' @param ce Effect site concentrations
#' @param gamma Slope at c50
#' @param E0 Effect at ce = 0
Emax1 <- function(pars, ce, gamma, E0) BISpred <- E0*(1 - ce^gamma/(ce^gamma + pars[1]^gamma))


# Corresponding inverse Emax function
#' @param pars c50
#' @param bis BIS values
#' @param E0 Effect at ce = 0
#' @param gamma Slope at c50
Hinv <- function(pars, bis, E0, gamma) unname(pars[1] * ((E0 - bis) /bis)^(1/gamma))
