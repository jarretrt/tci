# --------------------------------------------------------------------------------------------------------------------------------
# - PK-PD model helper functions and methods -------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------

#' dosing schedule
# create_intvl <- function(dose, inittm = 0){
#   b <- cut2(dose$time +inittm, breaks = c(inittm,dose$time+inittm), include.lowest = TRUE, right = FALSE)
#   ss <- t(sapply(stringr::str_extract_all(levels(b),"-?[0-9.]+"), as.numeric))
#   setNames(data.frame(levels(b), dose$infrt, ss), c("intvl","infrt","begin","end"))
#   # infm <- cbind(dose$infrt, ss)
#   # dimnames(infm) = list(rep(NULL,nrow(ss)), c("infrt","begin","end"))
#   # infm
# }


#' Create dosing schedule
#'
#' Create a dosing schedule object with columns "infrt", "begin", "end" from
#' vectors of infusions and infusion end times. The argument "inittm" is used
#' to specify the starting time of the first infusion.
#'
#' @param dose Data frame with columns "time" and "infrt".
#' @param inittm Starting time of initial infusion
#' @export
create_intvl <- function(dose, inittm = 0){
  if(!all(c("time","infrt") %in% colnames(dose)))
    stop("dose must include columnnames 'time','infrt'")
  tms <- c(inittm,dose[,"time"])
  as.matrix(cbind(infrt = dose[,"infrt"],
                  begin = tms[-length(tms)],
                  end = tms[-1]))
}
#' @examples
#' dose <- data.frame(time = c(0.5,4,4.5,10), infrt = c(100,0,100,0))
#' create_intvl(dose)



#' @name restrict_sigmoid
#' @title Restrict target sigmoid values
#' @description Function to place restriction on gamma and E50 parameters of target sigmoid
#' such that it passes through point (tfinal, BISfinal+eps)
#'
#' @param t50 parameter of Emax model
#' @param tfinal end of the induction period
#' @param eps distance between BISfinal and the target function at tfinal
#' @param BIS0 starting BIS value
#' @param BISfinal asymptote of Emax model
#' @export
restrict_sigmoid <- function(t50, tfinal =10, eps = 1, BIS0 = 100, BISfinal = 50-eps){
  gamma <- log((BIS0-BISfinal)/eps - 1, base = tfinal/t50)
  c(c50 = t50, gamma = gamma, e0 = BIS0, emx = BIS0 - BISfinal)
}


#' Generate variance-covariance matrix for Eleveld PK-PD model
#'
#' Generate the variance-covariance matrix for Eleveld PK-PD model for an observation
#' via Monte Carlo sampling.
#'
#' @param dat Data frame of observed patient covariates
#' @param N Number of Monte Carlo samples
#' @param rates Logical. Should rate constants be calculated
#' @param varnames Column names of variables used to calculate variance-covariance matrix
#'
#' @export
eleveld_vcov <- function(dat,
                         N = 1000,
                         rates = TRUE,
                         varnames = c("K10","K12","K21","K13","K31","V1","V2","V3","KE0","CE50","SIGMA")){

  if(rates){
    varnames <- c("K10","K12","K21","K13","K31","V1","V2","V3","KE0","CE50","SIGMA")
  } else{
    varnames <- c("CL","Q2","Q3","V1","V2","V3","KE0","CE50","SIGMA")
  }
  vcv_list <- lapply(1:nrow(dat), function(i){
    mc_samples <- replicate(N, log(unlist(
      eleveld_poppk(dat[i,], rate = rates, PD = TRUE, rand = TRUE)[,varnames])
    ))
    vcv <- cov(t(mc_samples))
    if(rcond(vcv) < 1e-5)
      diag(vcv) <- diag(vcv) + 1e-3
    round(vcv,5)
  })

  vcv_list
}


# All parameters in the Eleveld model are log-normally distributed or constant within the population.
# For each patient, there is a fixed set of model parameters predicted with variability around it.
# The fixed set of model parameters includes an estimate of the residual error standard deviation.
# The estimate of the sd is given by omega5 = 8.03 in the PD model and omega7 = 0.191 in the pk model.
# The residual error also is log-normally distributed --> log(err) ~ N(log(omega7), )

#' Get logged parameters updated in Eleveld model
#'
#' Extract the logged parameter values to be updated within the Eleveld model
#' from a data frame of patient PK-PD values.
#'
#' @param x Vector or data frame with Eleveld PK-PD model parameters
#' @param pd Logical. Should PD parameters be returned in addition to PK parameters.
#' @export
elvdlpars <- function(x, pd = TRUE){

  if(pd){
    x <- x[c("K10","K12","K21","K13","K31","V1","V2","V3","KE0","CE50","SIGMA")]
    names(x) <-  c("k10","k12","k21","k13","k31","v1","v2","v3","ke0","c50","sigma")
  } else{
    x <- x[c("K10","K12","K21","K13","K31","V1","V2","V3","KE0","LN_SIGMA")]
    names(x) <-  c("k10","k12","k21","k13","k31","v1","v2","v3","ke0","ln_sigma")
  }

  if(nrow(x) == 1)
    x <- as.numeric(x)

  log(x)
}


#' Set default PK parameter values
#' Set default PK parameter values for a pkmod object.
#' @param pkmod pkmod object
#' @param pars PK parameters to assign as default values of pkmod
#' @export
assign_pars <- function(pkmod, pars){

  if(class(pkmod) != "pkmod")
    stop("Class of pkmod must be 'pkmod'")
  if(!("pars" %in% names(formals(pkmod))))
    stop("Object 'pkmod' must have argument 'pars'")
  formals(pkmod)$pars <- pars
  class(pkmod) <- "pkmod"

  return(pkmod)
}



