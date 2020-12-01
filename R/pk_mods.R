# --------------------------------------------------------------------------------------------------------------------------------
# - Library of PK functions for 1, 2, and 3 compartment models with effect-site compartment for 3 cmpt model ---------------------
# --------------------------------------------------------------------------------------------------------------------------------

#' One compartment IV infusion with first-order elimination.
#'
#' @param tm Vector of times to evaluate the PK function at
#' @param kR Infusion rate (e.g. ml/min).
#' @param pars Named vector of parameters with names ('ke','v') or ('cl').
#' @param init Initial concentration. Defaults to 0.
#' @param inittm Time of initiation of infusion. Defaults to 0.
#' @examples
#' pkmod1cpt(1,1,c(ke = 0.5, v = 1))
#' @export
pkmod1cpt <- function(tm, kR, pars, init = 0, inittm = 0){

  names(pars) <- tolower(names(pars))
  tm <- tm - inittm

  if(any(!(c("ke","v") %in% names(pars))) & any(!(c("cl","v") %in% names(pars))))
    stop("pars must have names ('ke','v') or ('cl','v')")

  v <- pars["v"]
  if(!("ke" %in% names(pars)))
    ke <- pars["cl"] / v
  else
    ke <- pars["ke"]

  return((kR/ke*(1-exp(-tm*ke)) + init*v * exp(-tm*ke)) / v)
}
class(pkmod1cpt) <- "pkmod"


#' Solution to three-compartment IV model
#'
#' 3 compartment IV infusion with first-order absorption between compartments and with an additional effect-site compartment.
#' The analytical solutions implemented in this function are provided in "ADVAN-style analytical solutions for common pharmacokinetic models" by
#' Abuhelwa et al. 2015.
#'
#' This function takes in arguments for each of the absorption and elimination rate constants of a three-compartment model
#' as well as initial concentrations, c0. ke0 gives the rate of elimination from the effect-site compartment into the
#' central compartment (i.e. k41). The rate of absorption into the effect-site compartment is set at 1/10,000 the value of ke0.
#' The function returns a set of functions that calculate the concentration in each of the four compartments as a function of
#' time.
#' @param tm Vector of times to evaluate the PK function at
#' @param kR Infusion rate (e.g. ml/min).
#' @param pars Named vector of parameters with names (k10,k12,k21,k13,k31,v1,v2,v3,ke0)
#' @param init Initial concentration
#' @param inittm Time of initiation of infusion
#' @param returncpt Optionally specify a single compartment to return concentrations for.
#' Defaults to returning all compartment concentrations.
#' @examples
#' pars_3cpt <- c(k10=1.5,k12=0.15,k21=0.09,k13=0.8,k31=0.8,v1=10,v2=15,v3=100,ke0=1)
#' pkmod3cptm(1,1,pars_3cpt)
#' @export
pkmod3cptm <- function(tm, kR, pars, init = c(0,0,0,0), inittm = 0,
                       returncpt = c("all","cpt1","cpt2","cpt3","cpt4")) {

  names(pars) <- tolower(names(pars))

  if(any(!(c("k10","k12","k21","k13","k31","v1","v2","v3","ke0") %in% names(pars))))
    stop("pars must include names ('k10','k12','k21','k13','k31','v1','v2','v3','ke0')")

  k10 <- pars["k10"]
  k12 <- pars["k12"]
  k21 <- pars["k21"]
  k13 <- pars["k13"]
  k31 <- pars["k31"]
  v1  <- pars["v1"]
  v2  <- pars["v1"]
  v3  <- pars["v1"]
  kme <- pars["ke0"]

  if(!("k20" %in% names(pars))){
    k20 <- 0
  } else{
    k20 <- pars["k20"]
  }
  if(!("k30" %in% names(pars))){
    k30 <- 0
  } else{
    k30 <- pars["k30"]
  }
  if(!("km" %in% names(pars))){
    km  <- kme / 1e5
  } else{
    km <- pars["km"]
  }
  if(!("v4" %in% names(pars))){
    v4  <- v1 / 1e5
  } else{
    v4 <- pars["v4"]
  }

  returncpt <- match.arg(returncpt)
  tm <- tm - inittm

  E1 <- k10+k12+k13+km
  E2 <- k21+k20
  E3 <- k31+k30

  a <- E1+E2+E3
  b <- E1*E2+E3*(E1+E2)-k12*k21-k13*k31
  c <- E1*E2*E3-E3*k12*k21-E2*k13*k31

  m <- (3*b - a^2)/3
  n <- (2*a^3 - 9*a*b + 27*c)/27
  Q <- (n^2)/4 + (m^3)/27

  alpha <- sqrt(-1*Q)
  beta <- -1*n/2
  gamma <- sqrt(beta^2+alpha^2)
  theta <- atan2(alpha,beta)

  lambda1 <- a/3 + gamma^(1/3)*(cos(theta/3) + sqrt(3)*sin(theta/3))
  lambda2 <- a/3 + gamma^(1/3)*(cos(theta/3) - sqrt(3)*sin(theta/3))
  lambda3 <- a/3 -(2*gamma^(1/3)*cos(theta/3))

  A1last <- init[1]*v1
  A2last <- init[2]*v2
  A3last <- init[3]*v3
  Amlast <- init[4]*v4

  B = A2last*k21+A3last*k31
  C = E3*A2last*k21+E2*A3last*k31
  I = A1last*k12*E3-A2last*k13*k31+A3last*k12*k31
  J = A1last*k13*E2+A2last*k13*k21-A3last*k12*k21

  if(returncpt %in% c("all", "cpt1")){
    A1term1 = A1last*(exp(-tm*lambda1)*(E2-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-tm*lambda2)*(E2-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-tm*lambda3)*(E2-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))
    A1term2 = exp(-tm*lambda1)*(C-B*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-tm*lambda2)*(B*lambda2-C)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-tm*lambda3)*(B*lambda3-C)/((lambda1-lambda3)*(lambda3-lambda2))
    A1term3 = kR*((E2*E3)/(lambda1*lambda2*lambda3)-exp(-tm*lambda1)*(E2-lambda1)*(E3-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-tm*lambda2)*(E2-lambda2)*(E3-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-tm*lambda3)*(E2-lambda3)*(E3-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)))
    A1term = A1term1+A1term2+A1term3
  } else A1term = NULL


  if(returncpt %in% c("all", "cpt2")){
    A2term1 = A2last*(exp(-tm*lambda1)*(E1-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-tm*lambda2)*(E1-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-tm*lambda3)*(E1-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))
    A2term2 = exp(-tm*lambda1)*(I-A1last*k12*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-tm*lambda2)*(A1last*k12*lambda2-I)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-tm*lambda3)*(A1last*k12*lambda3-I)/((lambda1-lambda3)*(lambda3-lambda2))
    A2term3 = kR*k12*(E3/(lambda1*lambda2*lambda3)-exp(-tm*lambda1)*(E3-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-tm*lambda2)*(E3-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-tm*lambda3)*(E3-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)))
    A2term = A2term1+A2term2+A2term3
  } else A2term = NULL


  if(returncpt %in% c("all", "cpt3")){
    A3term1 = A3last*(exp(-tm*lambda1)*(E1-lambda1)*(E2-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-tm*lambda2)*(E1-lambda2)*(E2-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-tm*lambda3)*(E1-lambda3)*(E2-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))
    A3term2 = exp(-tm*lambda1)*(J-A1last*k13*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-tm*lambda2)*(A1last*k13*lambda2-J)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-tm*lambda3)*(A1last*k13*lambda3-J)/((lambda1-lambda3)*(lambda3-lambda2))
    A3term3 = kR*k13*(E2/(lambda1*lambda2*lambda3)-exp(-tm*lambda1)*(E2-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-tm*lambda2)*(E2-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-tm*lambda3)*(E2-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)))
    A3term = A3term1+A3term2+A3term3
  } else A3term = NULL

  if(returncpt %in% c("all", "cpt4")){
    Amterm1 = Amlast*exp(-tm*kme) +km*A1last*(exp(-tm*lambda1)*(E2-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1)*(kme-lambda1))+exp(-tm*lambda2)*(E2-lambda2)*(E3-lambda2)/((kme-lambda2)*(lambda1-lambda2)*(lambda3-lambda2))+exp(-tm*lambda3)*(E2-lambda3)*(E3-lambda3)/((kme-lambda3)*(lambda1-lambda3)*(lambda2-lambda3))+exp(-tm*kme)*(E2-kme)*(E3-kme)/((lambda1-kme)*(lambda2-kme)*(lambda3-kme)))
    Amterm2 = km*(exp(-tm*lambda1)*(B*lambda1-C)/((lambda1-lambda2)*(lambda1-lambda3)*(lambda1-kme))+exp(-tm*lambda2)*(C-B*lambda2)/((lambda1-lambda2)*(lambda2-lambda3)*(lambda2-kme))+exp(-tm*lambda3)*(C-B*lambda3)/((lambda1-lambda3)*(lambda3-lambda2)*(lambda3-kme))-exp(-tm*kme)*(B*kme-C)/((lambda1-kme)*(kme-lambda2)*(kme-lambda3)))
    Amterm3 = km*kR*((E2*E3)/(lambda1*lambda2*lambda3*kme)-exp(-tm*lambda1)*(E2-lambda1)*(E3-lambda1)/(lambda1*(kme-lambda1)*(lambda2-lambda1)*(lambda3-lambda1))-exp(-tm*lambda2)*(E2-lambda2)*(E3-lambda2)/(lambda2*(kme-lambda2)*(lambda1-lambda2)*(lambda3-lambda2))-exp(-tm*lambda3)*(E2-lambda3)*(E3-lambda3)/(lambda3*(kme-lambda3)*(lambda1-lambda3)*(lambda2-lambda3))-exp(-tm*kme)*(E2-kme)*(E3-kme)/(kme*(lambda1-kme)*(lambda2-kme)*(lambda3-kme)))
    Amterm = Amterm1+Amterm2+Amterm3
  } else Amterm = NULL

  return(rbind(A1term/v1, A2term/v2, A3term/v3, Amterm/v4))

}
class(pkmod3cptm) <- "pkmod"



#' Solution to three-compartment IV model
#'
#' 3 compartment IV infusion with first-order absorption between compartments and with an additional effect-site compartment.
#' The analytical solutions implemented in this function are provided in "ADVAN-style analytical solutions for common pharmacokinetic models" by
#' Abuhelwa et al. 2015.
#'
#' This function takes in arguments for each of the absorption and elimination rate constants of a three-compartment model
#' as well as initial concentrations, c0. ke0 gives the rate of elimination from the effect-site compartment into the
#' central compartment (i.e. k41). The rate of absorption into the effect-site compartment is set at 1/10,000 the value of ke0.
#' The function returns a set of functions that calculate the concentration in each of the four compartments as a function of
#' time.
#' @param kR Infusion rate (e.g. ml/min).
#' @param k10 Rate of excretion from central compartment.
#' @param k12 Rate of transfer from compartment 1 to compartment 2.
#' @param k21 Rate of transfer from compartment 2 to compartment 1.
#' @param k13 Rate of transfer from compartment 1 to compartment 3.
#' @param k31 Rate of transfer from compartment 3 to compartment 1.
#' @param v1 Volume of compartment 1.
#' @param v2 Volume of compartment 2.
#' @param v3 Volume of compartment 3.
#' @param ke0 Rate of transfer from effect-site compartment to compartment 1.
#' @param c0 Initial concentrations. Defaults to 0 in each compartment.
#' @examples
#' data(eleveld_pk)
#' data(eleveld_pd)
#' pk_vars <- c("V1","V2","V3","CL","Q2","Q3")
#' pd_vars <- c("E50","KE0","EMAX","GAM","GAM1","RESD")
#' pk_pars <- subset(eleveld_pk, ID == 403, select = pk_vars)
#' pd_pars <- subset(eleveld_pd, ID == 403, select = pd_vars)
#'
#' sol <- pk_basic_solution_3cpt_metab(kR = 1,
#'                                     k10 = pk_pars$CL / pk_pars$V1,
#'                                     k12 = pk_pars$Q2 / pk_pars$V1,
#'                                     k21 = pk_pars$Q2 / pk_pars$V2,
#'                                     k13 = pk_pars$Q3 / pk_pars$V1,
#'                                     k31 = pk_pars$Q3 / pk_pars$V3,
#'                                     v1 = pk_pars$V1,
#'                                     v2 = pk_pars$V2,
#'                                     v3 = pk_pars$V3,
#'                                     ke0 = pd_pars$KE0,
#'                                     c0 = c(0,0,0,0))
#' # concentration in central and effect site compartments
#' tms <- seq(0,1,0.1)
#' cbind(sol$c_1(tms), sol$c_4(tms))
#' @export
pk_basic_solution_3cpt_metab <- function(kR,k10,k12,k21,k13,k31,v1,v2,v3,ke0,
                                         c0=c(0,0,0,0))
  {

  kme <- ke0 # k41
  # k14 Absorption into the effect site is set much slower than elimination. As
  # soon as any drug enters, it is eliminated
  km  <- ke0 / 1e5
  # volume of effect-site compartment is set small enough that it doesn't influence
  # PK of other compartments
  v4  = v1 / 1e5
  k20 <- 0
  k30 <- 0
  E1 <- k10+k12+k13+km
  E2 <- k21+k20
  E3 <- k31+k30

  a <- E1+E2+E3
  b <- E1*E2+E3*(E1+E2)-k12*k21-k13*k31
  c <- E1*E2*E3-E3*k12*k21-E2*k13*k31

  m <- (3*b - a^2)/3
  n <- (2*a^3 - 9*a*b + 27*c)/27
  Q <- (n^2)/4 + (m^3)/27

  alpha <- sqrt(-1*Q)
  beta <- -1*n/2
  gamma <- sqrt(beta^2+alpha^2)
  theta <- atan2(alpha,beta)

  lambda1 <- a/3 + gamma^(1/3)*(cos(theta/3) + sqrt(3)*sin(theta/3))
  lambda2 <- a/3 + gamma^(1/3)*(cos(theta/3) - sqrt(3)*sin(theta/3))
  lambda3 <- a/3 -(2*gamma^(1/3)*cos(theta/3))

  A1last <- c0[1]*v1
  A2last <- c0[2]*v2
  A3last <- c0[3]*v3
  Amlast <- c0[4]*v4
  Doserate <- kR

  B = A2last*k21+A3last*k31
  C = E3*A2last*k21+E2*A3last*k31
  I = A1last*k12*E3-A2last*k13*k31+A3last*k12*k31
  J = A1last*k13*E2+A2last*k13*k21-A3last*k12*k21

  # return concentration functions
  a_1 <- function(t) {
    A1term1 = A1last*(exp(-t*lambda1)*(E2-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E2-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E2-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))
    A1term2 = exp(-t*lambda1)*(C-B*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(B*lambda2-C)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(B*lambda3-C)/((lambda1-lambda3)*(lambda3-lambda2))
    A1term3 = Doserate*((E2*E3)/(lambda1*lambda2*lambda3)-exp(-t*lambda1)*(E2-lambda1)*(E3-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-t*lambda2)*(E2-lambda2)*(E3-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-t*lambda3)*(E2-lambda3)*(E3-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)))
    return(A1term1+A1term2+A1term3)
  }

  a_2 <- function(t) {
    A2term1 = A2last*(exp(-t*lambda1)*(E1-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E1-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E1-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))
    A2term2 = exp(-t*lambda1)*(I-A1last*k12*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(A1last*k12*lambda2-I)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(A1last*k12*lambda3-I)/((lambda1-lambda3)*(lambda3-lambda2))
    A2term3 = Doserate*k12*(E3/(lambda1*lambda2*lambda3)-exp(-t*lambda1)*(E3-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-t*lambda2)*(E3-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-t*lambda3)*(E3-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)))
    return(A2term1+A2term2+A2term3)
  }

  a_3 <- function(t) {
    A3term1 = A3last*(exp(-t*lambda1)*(E1-lambda1)*(E2-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E1-lambda2)*(E2-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E1-lambda3)*(E2-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))
    A3term2 = exp(-t*lambda1)*(J-A1last*k13*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(A1last*k13*lambda2-J)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(A1last*k13*lambda3-J)/((lambda1-lambda3)*(lambda3-lambda2))
    A3term3 = Doserate*k13*(E2/(lambda1*lambda2*lambda3)-exp(-t*lambda1)*(E2-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-t*lambda2)*(E2-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-t*lambda3)*(E2-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)))
    return(A3term1+A3term2+A3term3)
  }

  a_4 <- function(t) {
    Amterm1 = Amlast*exp(-t*kme) +km*A1last*(exp(-t*lambda1)*(E2-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1)*(kme-lambda1))+exp(-t*lambda2)*(E2-lambda2)*(E3-lambda2)/((kme-lambda2)*(lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E2-lambda3)*(E3-lambda3)/((kme-lambda3)*(lambda1-lambda3)*(lambda2-lambda3))+exp(-t*kme)*(E2-kme)*(E3-kme)/((lambda1-kme)*(lambda2-kme)*(lambda3-kme)))
    Amterm2 = km*(exp(-t*lambda1)*(B*lambda1-C)/((lambda1-lambda2)*(lambda1-lambda3)*(lambda1-kme))+exp(-t*lambda2)*(C-B*lambda2)/((lambda1-lambda2)*(lambda2-lambda3)*(lambda2-kme))+exp(-t*lambda3)*(C-B*lambda3)/((lambda1-lambda3)*(lambda3-lambda2)*(lambda3-kme))-exp(-t*kme)*(B*kme-C)/((lambda1-kme)*(kme-lambda2)*(kme-lambda3)))
    Amterm3 = km*Doserate*((E2*E3)/(lambda1*lambda2*lambda3*kme)-exp(-t*lambda1)*(E2-lambda1)*(E3-lambda1)/(lambda1*(kme-lambda1)*(lambda2-lambda1)*(lambda3-lambda1))-exp(-t*lambda2)*(E2-lambda2)*(E3-lambda2)/(lambda2*(kme-lambda2)*(lambda1-lambda2)*(lambda3-lambda2))-exp(-t*lambda3)*(E2-lambda3)*(E3-lambda3)/(lambda3*(kme-lambda3)*(lambda1-lambda3)*(lambda2-lambda3))-exp(-t*kme)*(E2-kme)*(E3-kme)/(kme*(lambda1-kme)*(lambda2-kme)*(lambda3-kme)))
    return(Amterm1+Amterm2+Amterm3)
  }

  c_1 <- function(t) a_1(t)/v1
  c_2 <- function(t) a_2(t)/v2
  c_3 <- function(t) a_3(t)/v3
  c_4 <- function(t) a_4(t)/v4

  return(list(c_1=c_1, c_2=c_2, c_3=c_3, c_4=c_4))
}



#' @name pk_solution_3cpt_metab
#' @title Iterate solution to three-compartment model
#' @description This function extends the function pk_basic_solution_3cpt_metab to a specified infusion schedule, rather than a single
#' infusion.
#' @param pars Named vector of parameters for a 3-compartment model with effect-site.
#' @param ivt Infusion schedule given in the form of a named list
#' (e.g. list(list(begin = 0, end = 2, k_R = 1), list(begin = 4, end = 6, k_R = 1)))
#' @param init inital concentrations for the 4 compartments.
#' @examples
#' data(eleveld_pk)
#' pk_pars <- subset(eleveld_pk, ID == 403, c("V1","V2","V3","CL","Q2","Q3"))
#' pd_pars <- subset(eleveld_pd, ID == 403, c("E50","KE0","EMAX","GAM","GAM1","RESD"))
#' pars <- c(k10 = pk_pars$CL / pk_pars$V1,
#'           k12 = pk_pars$Q2 / pk_pars$V1,
#'           k21 = pk_pars$Q2 / pk_pars$V2,
#'           k13 = pk_pars$Q3 / pk_pars$V1,
#'           k31 = pk_pars$Q3 / pk_pars$V3,
#'           v1 = pk_pars$V1,
#'           v2 = pk_pars$V2,
#'           v3 = pk_pars$V3,
#'           ke0 = pd_pars$KE0)
#' ivt <- list(list(begin=0.0, end=0.5, k_R=6),
#'             list(begin=8.0, end=8.5, k_R=6),
#'             list(begin=16.0, end=16.5, k_R=6),
#'             list(begin=24.0, end=24.5, k_R=6),
#'             list(begin=32.0, end=32.5, k_R=6))
#' init <- c(0,0,0,0)
#' sol <- pk_solution_3cpt_metab(pars, ivt, init)
#' sol(seq(0,32))
#' @export
pk_solution_3cpt_metab <- function(pars, ivt, init)
{
  k_10 = pars[1]
  k_12 = pars[2]
  k_21 = pars[3]
  k_13 = pars[4]
  k_31 = pars[5]
  v_1  = pars[6]
  v_2  = pars[7]
  v_3  = pars[8]
  k_e0 = pars[9]

  ## create a list of event times
  ibe <- sapply(ivt, `[[`, 'begin')
  ied <- sapply(ivt, `[[`, 'end')
  prd <- sort(unique(c(0, c(ibe,ied), Inf)))
  rits <- list()
  ## compute basic solution in each interval
  for(i in 1:(length(prd)-1)) {
    civt <- sapply(ivt, function(iv) {
      if(prd[i] >= iv$begin && prd[i] < iv$end) {
        iv$k_R
      } else { 0 }
    })

    rit <- list(begin=prd[i], end=prd[i+1],
                idose=sum(civt))
    if(i == 1) {
      rit$init <- init
    } else {
      rit$init <- c(rits[[i-1]]$c_1(rits[[i-1]]$end-rits[[i-1]]$begin),
                    rits[[i-1]]$c_2(rits[[i-1]]$end-rits[[i-1]]$begin),
                    rits[[i-1]]$c_3(rits[[i-1]]$end-rits[[i-1]]$begin),
                    rits[[i-1]]$c_4(rits[[i-1]]$end-rits[[i-1]]$begin))
    }

    sol <- pk_basic_solution_3cpt_metab(kR=rit$idose, k10 = k_10, k12 = k_12, k21 = k_21, k13 = k_13,
                                          k31 = k_31, v1 = v_1, v2 = v_2, v3 = v_3, ke0 = k_e0, c0=rit$init)

    rits[[i]] <- c(rit, sol)
  }

  ret <- function(tms) {
    sapply(tms, function(t) {
      val <- rep(NA,4)
      for(rit in rits) {
        if(t >= rit$begin && t <= rit$end) {
          val <- c(rit$c_1(t-rit$begin),
                   rit$c_2(t-rit$begin),
                   rit$c_3(t-rit$begin),
                   rit$c_4(t-rit$begin))
          break
        }
      }
      val
    })
  }
  return(ret)
}



#' @name pk_solution_3cpt_metab_singleinf
#' @title Single infusion 3-compartment PK solution
#' @description Piece-wise solution for a single infusion followed by a period with no infusion.
#' This function is similar to pk_solution_3cpt_metab, except that it accepts and
#' implements only the first infusion. This function exists primarily for
#' reducing computational speed when searching for time until maximum concentration.
#' @param pars Named vector of parameters for a 3-compartment model with effect-site.
#' @param ivt Infusion schedule given in the form of a named list
#' (e.g. list(list(begin = 0, end = 2, k_R = 1), list(begin = 4, end = 6, k_R = 1)))
#' @param init Inital concentrations for the 4 compartments.
#' @param ce_only Logical. Should only the effect-site concentration be returned.
#' Defaults to FALSE
#' @examples
#' data(eleveld_pk)
#' data(eleveld_pd)
#' pk_pars <- subset(eleveld_pk, ID == 403, select = c("V1","V2","V3","CL","Q2","Q3"))
#' pd_pars <- subset(eleveld_pd, ID == 403, select = c("E50","KE0","EMAX","GAM","GAM1","RESD"))
#'
#' pars <- c(k10 = pk_pars$CL / pk_pars$V1,
#'           k12 = pk_pars$Q2 / pk_pars$V1,
#'           k21 = pk_pars$Q2 / pk_pars$V2,
#'           k13 = pk_pars$Q3 / pk_pars$V1,
#'           k31 = pk_pars$Q3 / pk_pars$V3,
#'           v1 = pk_pars$V1,
#'           v2 = pk_pars$V2,
#'           v3 = pk_pars$V3,
#'           ke0 = pd_pars$KE0)
#' ivt <- list(begin = 0, end = 0.5, k_R = 1)
#' init <- c(0,0,0,0)
#' sol <- pk_solution_3cpt_metab_singleinf(pars, ivt, init)
#' sol(seq(0,5,0.1))
#' @export
pk_solution_3cpt_metab_singleinf <- function(pars, ivt, init, ce_only = FALSE){

  k_10 = pars[1]
  k_12 = pars[2]
  k_21 = pars[3]
  k_13 = pars[4]
  k_31 = pars[5]
  v_1  = pars[6]
  v_2  = pars[7]
  v_3  = pars[8]
  k_e0 = pars[9]

  begin <- ivt$begin
  end <- ivt$end
  delta <- end - begin

  sol_inf <- pk_basic_solution_3cpt_metab(kR=ivt$k_R, k10 = k_10, k12 = k_12, k21 = k_21, k13 = k_13,
                                            k31 = k_31, v1 = v_1, v2 = v_2, v3 = v_3, ke0 = k_e0, c0=init)
  c_end_inf <- c(sol_inf$c_1(delta), sol_inf$c_2(delta), sol_inf$c_3(delta), sol_inf$c_4(delta))
  sol_coast <- pk_basic_solution_3cpt_metab(kR=0, k10 = k_10, k12 = k_12, k21 = k_21, k13 = k_13,
                                              k31 = k_31, v1 = v_1, v2 = v_2, v3 = v_3, ke0 = k_e0, c0=c_end_inf)

  if(ce_only){
    calc_con <- function(tms){
      tms_inf <- tms[tms <= end]
      tms_coast <- tms[tms > end]
      return(c(sol_inf$c_4(tms_inf-begin), sol_coast$c_4(tms_coast-end)))
    }
  } else{
    calc_con <- function(tms){
      tms_inf <- tms[tms <= end]
      tms_coast <- tms[tms > end]
      val_inf <- cbind(sol_inf$c_1(tms_inf-begin), sol_inf$c_2(tms_inf-begin), sol_inf$c_3(tms_inf-begin), sol_inf$c_4(tms_inf-begin))
      val_coast <- cbind(sol_coast$c_1(tms_coast-end), sol_coast$c_2(tms_coast-end), sol_coast$c_3(tms_coast-end), sol_coast$c_4(tms_coast-end))
      return(t(rbind(val_inf, val_coast)))
    }
  }
  return(calc_con)
}





#' @name gen_eleveld_pd_pars
#' @title Eleveld model PD parameters
#' @description Function to generate PD parameters for Eleveld model.
#' @param theta Vector of fixed effects
#' @param eta Vector of random effects
#' @param patient_vars Named list of observed patient characteristics
#' @examples
#' data(eleveld_pd)
#' # PD fixed effect values and random effect variances from Eleveld et al. (2018)
#' eleveld_theta_pd_est <- c(3.08,0.146,93.0,1.47,8.03,0.0517,-0.00635,1.24,1.89)
#' eleveld_eta_pd_var <- c(0.242,0.702,0.230)
#' patient_covariates <- subset(eleveld_pd, ID == 403,
#' select = c("AGE","WGT","HGT","M1F2","PMA","TECH","A1V2"))
#' eta_obs <- c(mvtnorm::rmvnorm(1,sigma = diag(eleveld_eta_pd_var)))
#' gen_eleveld_pd_pars(theta = eleveld_theta_pd_est,
#'                     eta = eleveld_eta_pd_var,
#'                     patient_vars = patient_covariates)
#' @export
gen_eleveld_pd_pars <- function(theta, eta, patient_vars){
 AGE = patient_vars$AGE
 WGT = patient_vars$WGT
 ARTERIAL = ifelse(patient_vars$A1V2 == 1,1,0)
 AGEref=35
 faging <- function(x) exp(x*(AGE-AGEref))
 Ce50 = theta[1]*faging(theta[7])*exp(eta[1])
 ke0 = ifelse(ARTERIAL, theta[2]*(WGT/70)^(-0.25)*exp(eta[2]), theta[8]*(WGT/70)^(-0.25)*exp(eta[2]))
 BISbaseline = theta[3]
 gamma = theta[4]
 gamma2 = theta[9]
 sigma = theta[5]*exp(eta[3])
 bis_delay = 15 + exp(theta[6]*AGE)
 return(c(ke0 = ke0, c50 = Ce50, gamma = gamma, gamma2 = gamma2, E0 = BISbaseline, sigma = sigma, bis_delay = bis_delay))
}



#' @name gen_eleveld_pk_pars
#' @title Eleveld model PK parameters
#' @description Function to generate PK parameters for Eleveld model.
#' @param theta Vector of fixed effects
#' @param eta Vector of random effects
#' @param patient_vars Named list of observed patient characteristics
#' @param returnQ Logical. Should clearance be returned instead of rates
#'
#' @examples
#' data(eleveld_pk)
#' # PK fixed effect values and random effect variances from Eleveld et al. (2018)
#' eleveld_theta_pk_est <- c(6.28,25.5,273,1.79,1.75,1.11,0.191,42.3,9.06,-0.0156,
#' -0.00286,33.6,-0.0138,68.3,2.10,1.30,1.42,0.68)
#' eleveld_eta_pk_var <- c(0.610,0.565,0.597,0.265,0.346,0.209,0.463)
#'
#' # Example patient covariate values, fixed effects, and random effects
#' vars <- c("AGE","WGT","HGT","M1F2","PMA","TECH","BMI","FFM","A1V2")
#' patient_covariates <- subset(eleveld_pk, ID == 403, select = vars)
#' eta_obs <- c(mvtnorm::rmvnorm(1,sigma = diag(eleveld_eta_pk_var)))
#' gen_eleveld_pk_pars(theta = eleveld_theta_pk_est,
#'                     eta = eta_obs,
#'                     patient_vars = patient_covariates)
#' @export
gen_eleveld_pk_pars <- function(theta, eta, patient_vars, returnQ = FALSE){

  AGE  <- patient_vars$AGE
  PMA  <- patient_vars$PMA
  WGT  <- patient_vars$WGT
  HGT  <- patient_vars$HGT
  M1F2 <- patient_vars$M1F2
  TECH <- patient_vars$TECH
  A1V2 <- patient_vars$A1V2

  MALE <- ifelse(M1F2 == 1, 1, 0)
  OPIATE <- ifelse(TECH == 2, 1, 0)
  ARTERIAL <- ifelse(A1V2 == 1, 1, 0)
  AGEref=35
  WGTref=70
  HGTref=170

  if(is.na(PMA)) PMA = AGE + 40/52
  PMAref = AGEref + 40/52
  BMI = 10000 * WGT / HGT / HGT
  BMIref = 10000 * WGTref / HGTref / HGTref

  V1ref = theta[1]
  V2ref = theta[2]
  V3ref = theta[3]

  faging <- function(x) exp(x*(AGE-AGEref))
  fsigmoid <- function(x,E50,lambda) x^lambda / (x^lambda + E50^lambda)
  fcentral <- function(x) fsigmoid(x,theta[12],1)
  fCLmaturation <- fsigmoid(PMA*52,theta[8],theta[9])
  fCLmaturation_ref = fsigmoid(PMAref*52,theta[8],theta[9])
  fQ3maturation <- fsigmoid(AGE*52+40,theta[14],1)
  fQ3maturation_ref <- fsigmoid(AGEref*52+40,theta[14],1)
  fopiates <- function(x) ifelse(OPIATE,exp(x*AGE),1)

  fAlSallami <-     ifelse(MALE, (0.88+(1-0.88)/(1+(AGE/13.4)^(-12.7)))*(9270*WGT)/(6680+216*BMI),
                           (1.11+(1-1.11)/(1+(AGE/7.1)^(-1.1)))*(9270*WGT)/(8780+244*BMI))

  fAlSallami_ref <- (0.88+(1-0.88)/(1+(AGEref/13.4)^(-12.7))) *(9270*WGTref)/(6680+216*BMIref)
  V1arterial <- theta[1]*fcentral(WGT)/fcentral(WGTref)*exp(eta[1])
  V1venous <- V1arterial*(1+theta[17]*(1-fcentral(WGT)))
  V2 <- theta[2]*WGT/WGTref*faging(theta[10])*exp(eta[2])
  V3 <- theta[3]*(fAlSallami/fAlSallami_ref)*fopiates(theta[13])*exp(eta[3])

  CL <- c(t(c(MALE,1-MALE))%*%theta[c(4,15)]*(WGT/WGTref)^0.75*(fCLmaturation/fCLmaturation_ref)*fopiates(theta[11])*exp(eta[4]))
  Q2arterial <- theta[5]*(V2/V2ref)^(0.75)*(1+theta[16]*(1-fQ3maturation))*exp(eta[5])
  Q2venous <- Q2arterial*theta[18]
  Q3 <- theta[6]*(V3/V3ref)^(0.75)*(fQ3maturation/fQ3maturation_ref)*exp(eta[6])

  V1 = ifelse(ARTERIAL, V1arterial, V1venous)
  K10 = CL/V1
  K12 = ifelse(ARTERIAL, Q2arterial/V1, Q2venous/V1)
  K21 = ifelse(ARTERIAL, Q2arterial/V2, Q2venous/V2)
  K13 = Q3/V1
  K31 = Q3/V3

  if(!returnQ){
    return(c(k10 = K10, k12 = K12, k21 = K21, k13 = K13, k31 = K31, v1 = V1, v2 = V2, v3 = V3))
  } else{
    Q2 <- ifelse(ARTERIAL, Q2arterial, Q2venous)
    return(c(CL = CL, Q2 = Q2, Q3 = Q3, v1 = V1, v2 = V2, v3 = V3))
  }
}




#' @name gen_eleveld_pk_pars_nonmem
#' @title Generate Eleveld model PK parameters
#' @description R code adapted from NONMEM PK file provided in supplementary material of Eleveld et al. Function takes in fixed effect
#' parameter estimates and random effect variance estimates to return parameters for a 3 compartment pk model with an
#' effect site compartment.
#' @param THETA Vector of fixed effects
#' @param ETA Vector of random effect variances
#' @param PATIENT_VARS Named list of patient covariate values
#' @param returnQ Optional logical value to indicate if clearance values should be returned instead of elimination rate constants.
#' @examples
#' data(eleveld_pk)
#'
#' # PK fixed effect values and random effect variances from Eleveld et al. (2018)
#' eleveld_theta_pk_est <- c(6.28,25.5,273,1.79,1.75,1.11,0.191,42.3,9.06,
#' -0.0156,-0.00286,33.6,-0.0138,68.3,2.10,1.30,1.42,0.68)
#' eleveld_eta_pk_var <- c(0.610,0.565,0.597,0.265,0.346,0.209,0.463)
#'
#' # Example patient covariate values, fixed effects, and random effects
#' patient_covariates <- subset(eleveld_pk, ID == 403,
#' select = c("AGE","WGT","HGT","M1F2","PMA","TECH","A1V2"))
#' eta_obs <- c(mvtnorm::rmvnorm(1,sigma = diag(eleveld_eta_pk_var)))
#' gen_eleveld_pk_pars_nonmem(THETA = eleveld_theta_pk_est,
#'                            ETA = eta_obs,
#'                            PATIENT_VARS = patient_covariates)
#'
#' @export
gen_eleveld_pk_pars_nonmem <- function(THETA, ETA, PATIENT_VARS, returnQ = FALSE){

  # convert to logged values - add one when negative
  neg_ix <- which(THETA<0)
  THETA[-neg_ix] <- log(THETA[-neg_ix])
  THETA[neg_ix]  <- log(1+THETA[neg_ix])

  AGE  <- PATIENT_VARS$AGE
  PMA  <- PATIENT_VARS$PMA
  WGT  <- PATIENT_VARS$WGT
  HGT  <- PATIENT_VARS$HGT
  M1F2 <- PATIENT_VARS$M1F2
  TECH <- PATIENT_VARS$TECH
  A1V2 <- PATIENT_VARS$A1V2

  # Al-sallami FFM
  # 1.7 = reference height (m)
  # 70 = reference weight (kg)
  HT2=(HGT/100)*(HGT/100)
  MATM=0.88+((1-0.88)/(1+(AGE/13.4)^(-12.7)))
  MATF=1.11+((1-1.11)/(1+(AGE/7.1)^(-1.1)))
  MATR=0.88+((1-0.88)/(1+(35/13.4)^(-12.7)))
  FFMM=MATM*42.92*(HT2)*WGT/(30.93*(HT2)+WGT)
  FFMF=MATF*37.99*(HT2)*WGT/(35.98*(HT2)+WGT)
  FFMR=MATR*42.92*(1.7*1.7)*70/(30.93*(1.7*1.7)+70)
  FFM=FFMM*(2-M1F2) + FFMF*(M1F2-1)
  NFFM=FFM/FFMR

  # maturation
  DV1=1
  DV2=1
  DV3=1

  # sigmoidal maturation of CL
  PMW=PMA*52
  PMR=(35+40/52)*52
  ME50=exp(THETA[8])
  MGAM=exp(THETA[9])
  MCL=(PMW^MGAM)/(PMW^MGAM+ME50^MGAM)
  RCL=(PMR^MGAM)/(PMR^MGAM+ME50^MGAM)
  DCL=MCL/RCL
  DQ2=1

  # sigmoidal maturation of Q3 based on 40 weeks gestation
  PMEW=AGE*52+40
  PMER=35*52+40
  QE50=exp(THETA[14])
  MQ3=PMEW/(PMEW+QE50)
  RQ3=PMER/(PMER+QE50)
  DQ3=MQ3/RQ3

  # aging
  KV1=1
  KV2=exp(THETA[10]*(AGE-35))
  KV3=exp(THETA[13]*(AGE)*(TECH-1))
  KCL=exp(THETA[11]*(AGE)*(TECH-1))
  KQ2=1
  KQ3=1

  # covariate structure
  # V1 scales sigmoid with weight
  VV50=exp(THETA[12])
  CV1=WGT/(WGT+VV50)
  RV1=70/(70+VV50)
  M1 =(CV1/RV1) * KV1 * DV1
  VCV1=(A1V2-1)*(1-CV1)
  V1 =exp(THETA[1]+ETA[1]) * M1 * (1+VCV1*exp(THETA[17]))
  M2 =(WGT/70)^1 * KV2 * DV2
  V2 =exp(THETA[2]+ETA[2]) * M2
  M3 =(NFFM)^1 * KV3 * DV3 # this corresponds to the middle part of the V3 equation
  V3 =exp(THETA[3]+ETA[3]) * M3
  M4 =(WGT/70)^0.75 * KCL * DCL
  CL =exp((2-M1F2)*THETA[4]+(M1F2-1)*THETA[15]+ETA[4]) * M4
  RV2=exp(THETA[2])
  M5 =(V2/RV2)^0.75 * KQ2 * DQ2
  KM5=1+exp(THETA[16])*(1-MQ3)
  Q2 =exp(THETA[5]+ETA[5]+(A1V2-1)*THETA[18]) * M5 * KM5
  RV3=exp(THETA[3])
  M6 =(V3/RV3)^0.75 * KQ3 * DQ3
  Q3 =exp(THETA[6]+ETA[6]) * M6

  K10 = CL/V1
  K12 = Q2/V1
  K21 = Q2/V2
  K13 = Q3/V1
  K31 = Q3/V3

  if(returnQ) return(c(CL = CL, Q2 = Q2, Q3 = Q3, V1 = V1, V2 = V2, V3 = V3))
  else return(c(k10 = K10, k12 = K12, k21 = K21, k13 = K13, k31 = K31, V1 = V1, V2 = V2, V3 = V3))
}








