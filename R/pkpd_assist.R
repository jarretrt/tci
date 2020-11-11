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



#' Restrict target sigmoid values
#'
#' Function to place restriction on gamma and E50 parameters of target sigmoid
#' such that it passes through point (tfinal, BISfinal+eps)
#'
#' @param t50 parameter of Emax model
#' @param tfinal end of the induction period
#' @param eps distance between BISfinal and the target function at tfinal
#' @param BIS0 starting BIS value
#' @param BISfinal asymptote of Emax model
#'
#' @export
restrict_sigmoid <- function(t50, tfinal =10, eps = 1, BIS0 = 100, BISfinal = 50-eps){
  gamma <- log((BIS0-BISfinal)/eps - 1, base = tfinal/t50)
  c(c50 = t50, gamma = gamma, e0 = BIS0, emx = BIS0 - BISfinal)
}


#' Get variance-covariance matrix for population PK-PD models
#'
#' NOTE: This needs to be modified. There are correlations between parameters
#' induced by patient covariates.
#'
#' Function to extract covariance from population pk or pk-pd models
#' This gives the covariance for the random effects about the fixed effect estimates.
#' The variance estimates for the residual error terms give the variance of error terms
#' within the population about the logged
#'
#' @param poppk Selection of Schnnider PK model or Eleveld PK-PD model
#' @param pd Logical. Should PD parameters be returned for Eleveld model
#'
#' @export
poppk_cov <- function(poppk = c("Schnider","Eleveld"), pd = TRUE){
  warning("This function isn't working properly. Covariance terms should be non-zero
          for some parameters.")
  poppk <- match.arg(poppk)
  if(poppk == "Schnider") out <- diag(c(0.278,2.330,34.900,0.059,0.112,0.044,0.070,0.009,0.017,0.009,0.005))
  if(poppk == "Eleveld"){
    lvars <- c(0.610,0.565,0.597,0.265,0.346,0.209,0.463,0.242,0.702,0.230)
    names(lvars) <- c("v1","v2","v3","cl","q2","q3","resid_pk","ce50","ke0","resid_pd")
    lvars <- c(lvars, c(k10 = unname(lvars["cl"] + lvars["v1"]),
                        k12 = unname(lvars["q2"] + lvars["v1"]),
                        k21 = unname(lvars["q2"] + lvars["v2"]),
                        k13 = unname(lvars["q3"] + lvars["v1"]),
                        k31 = unname(lvars["q3"] + lvars["v3"])))
    if(pd){
      return(diag(lvars[c("k10","k12","k21","k13","k31","v1","v2","v3","ke0","ce50","resid_pd")]))
    } else{
      return(diag(lvars[c("k10","k12","k21","k13","k31","v1","v2","v3","ke0","resid_pk")]))
    }
  }
}
#' @examples poppk_cov("Eleveld", pd = TRUE)


# All parameters in the eleveld model are log-normally distributed or constant within the population.
# For each patient, there is a fixed set of model parameters predicted with variablilty around it.
# The fixed set of model parameters includes an estimate of the residual error standard deviation.
# The estimate of the sd is given by omega5 = 8.03 in the PD model and omega7 = 0.191 in the pk model.
# The residual error also is log-normally distributed --> log(err) ~ N(log(omega7), )

#' Get logged parameters updated in Eleveld model
#'
#' Extract the logged parameter values to be updated within the Eleveld model
#' from a data frame of patient PK-PD values.
#'
#' @param x vector or data frame with Eleveld PK-PD model parameters
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



