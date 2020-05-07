# PD functions

#' Emax function. Assumes maximum value of Emx = 100 and maximum change in value of 100 (i.e. minimum value of 0).
#' @param pars Vector of parameters: (c50,gamma).
#' @param ce Effect-site concentration.
#' @param Emx Maximum value of function.
#' @param E0 Maximum difference between highest and lowest function value.
Emax <- function(pars, ce, Emx = 100, E0 = 100) E0 - Emx*(ce^pars[2] / (ce^pars[2] + pars[1]^pars[2]))

#' @examples
#' pk_pars <- eleveld_pk[eleveld_pk$ID == 403,c("V1","V2","V3","CL","Q2","Q3")]
#' pd_pars <- eleveld_pd[eleveld_pd$ID == 403,c("E50","KE0","EMAX","GAM","GAM1","RESD")]
#' pars <- c(k10 = pk_pars$CL / pk_pars$V1,
#'           k12 = pk_pars$Q2 / pk_pars$V1,
#'           k21 = pk_pars$Q2 / pk_pars$V2,
#'           k13 = pk_pars$Q3 / pk_pars$V1,
#'           k31 = pk_pars$Q3 / pk_pars$V3,
#'           v1 = pk_pars$V1,
#'           v2 = pk_pars$V2,
#'           v3 = pk_pars$V3,
#'           ke0 = pd_pars$KE0,
#'           c50 = pd_pars$E50,
#'           gamma = pd_pars$GAM,
#'           E0 = pd_pars$EMAX)
#' ivt <- list(list(begin=0.0, end=0.5, k_R=100),
#'             list(begin=0.5, end=32, k_R=3))
#' tms <- seq(0,32,0.1)
#' sol <- pk_solution_3cpt_metab(pars, ivt, init  = c(0,0,0,0))
#' plot(tms, Emax(pars = pars[c("c50","gamma")], ce = sol(tms)[4,]), type = "l", ylim = c(0,100), ylab = "BIS", xlab = "Minutes")


#' Decreasing Emax function with only c50 estimated and maximum effect assumed to be 0.
#' @param pars c50
#' @param ce Effect site concentrations
#' @param gamma Slope at c50
#' @param E0 Effect at ce = 0
Emax1 <- function(pars, ce, gamma, E0) BISpred <- E0*(1 - ce^gamma/(ce^gamma + pars[1]^gamma))
#' @examples
#' pk_pars <- eleveld_pk[eleveld_pk$ID == 403,c("V1","V2","V3","CL","Q2","Q3")]
#' pd_pars <- eleveld_pd[eleveld_pd$ID == 403,c("E50","KE0","EMAX","GAM","GAM1","RESD")]
#' pars <- c(k10 = pk_pars$CL / pk_pars$V1,
#'           k12 = pk_pars$Q2 / pk_pars$V1,
#'           k21 = pk_pars$Q2 / pk_pars$V2,
#'           k13 = pk_pars$Q3 / pk_pars$V1,
#'           k31 = pk_pars$Q3 / pk_pars$V3,
#'           v1 = pk_pars$V1,
#'           v2 = pk_pars$V2,
#'           v3 = pk_pars$V3,
#'           ke0 = pd_pars$KE0,
#'           c50 = pd_pars$E50,
#'           gamma = pd_pars$GAM,
#'           E0 = pd_pars$EMAX)
#' ivt <- list(list(begin=0.0, end=0.5, k_R=100),
#'             list(begin=0.5, end=32, k_R=3))
#' tms <- seq(0,32,0.1)
#' sol <- pk_solution_3cpt_metab(pars, ivt, init  = c(0,0,0,0))
#' plot(tms, Emax1(pars = pars["c50"], ce = sol(tms)[4,], gamma = pars["gamma"], E0 = pars["E0"]), type = "l", ylim = c(0,100), ylab = "BIS", xlab = "Minutes")


# Corresponding inverse Emax function
#' @param pars c50
#' @param bis BIS values
#' @param E0 Effect at ce = 0
#' @param gamma Slope at c50
Hinv <- function(pars, bis, E0, gamma) unname(pars[1] * ((E0 - bis) /bis)^(1/gamma))
#' @examples
#' Hinv(pars = 3, bis = 60, E0 = 100, gamma = 1.47)




