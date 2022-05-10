# --------------------------------------------------------------------------------------------------------------------------------
# - PK-PD model helper functions and methods -------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------



#' @name format_pars
#' @title Format parameters for use in Rcpp functions
#'
#' Order parameters for 1-4 compartment models to be used in Rcpp functions in
#' predict_pkmod method.
#'
#' @param pars Vector of named parameters. Names can be capitalized or lowercase
#' and can include variations of "V1" as "V" or clearance terms rather than
#' elimination rate constants.
#' @param ncmpt Number of compartments in the model. This should be a value
#' between 1 and 4. If ncmpt = 4, it assumes that the fourth compartment is an
#' effect-site without a corresponding volume parameter.
#' @return Numeric vector of transformed parameter values.
#' @examples
#' format_pars(c(V1 = 8.9, CL = 1.4, q2 = 0.9, v2 = 18), ncmpt = 2)
#' format_pars(c(V1 = 8.9, CL = 1.4, q2 = 0.9, v2 = 18, cl2 = 3), ncmpt = 2)
#' @export
format_pars <- function(pars, ncmpt = 3){

  names(pars) <- tolower(names(pars))

  if("v" %in% names(pars)){
    v1 <- pars["v"]
  } else{
    v1 <- pars["v1"]
  }
  if("cl" %in% names(pars)){
    k10 <- pars["cl"]/v1
  } else if("cl1" %in% names(pars)){
    k10 <- pars["cl1"]/v1
  } else if("ke" %in% names(pars))
    k10 <- pars["ke"]
  else{
    k10 <- pars["k10"]
  }
  if(ncmpt >= 2){
    v2 <- pars["v2"]
    if("q" %in% names(pars)){
      k12 = pars["q"]/v1
      k21 = pars["q"]/v2
    } else if("q2" %in% names(pars)){
      k12 = pars["q2"]/v1
      k21 = pars["q2"]/v2
    } else{
      k12 = pars["k12"]
      k21 = pars["k21"]
    }

    if("cl2" %in% names(pars)){
      k20 = pars["cl2"]/v2
    } else if("k20" %in% names(pars))
      k20 = pars["k20"]
    else
      k20 = 0

  }
  if(ncmpt >= 3){
    v3 <- pars["v3"]
    if("q3" %in% names(pars)){
      k13 = pars["q3"]/v1
      k31 = pars["q3"]/v3
    } else{
      k13 = pars["k13"]
      k31 = pars["k31"]
    }

    if("cl3" %in% names(pars)){
      k30 = pars["cl3"]/v3
    } else if("k30" %in% names(pars))
      k30 = pars["k30"]
    else
      k30 = 0
  }

  if(ncmpt >= 3){
    pars_out <- unname(c(k10,k20,k30,k12,k21,k13,k31,v1,v2,v3))
    names(pars_out) <- c("k10","k20","k30","k12","k21","k13","k31","v1","v2","v3")
  } else if(ncmpt == 2){
    pars_out <- unname(c(k10,k20,k12,k21,v1,v2))
    names(pars_out) <- c("k10","k20","k12","k21","v1","v2")
  } else{
    pars_out <- unname(c(k10,v1))
    names(pars_out) <- c("k10","v1")
  }

  if("ke0" %in% names(pars)){
    pars_out <- c(pars_out, pars["ke0"])
  }
  return(pars_out)
}


#' @name infer_pkfn
#' @title Identify pkfn from parameter names
#' @description Identify structural PK model function (i.e., `pkfn`) from parameter names.
#' Models available are 1-, 2-, and 3-compartment mammillary models, or 3-compartment with
#' an effect site, corresponding to functions `pkmod1cpt`, `pkmod2cpt`, `pkmod3cpt`, and
#' `pkmod3cptm`, respectively.
#' @param parnms Vector of parameter names.
#' @return Returns one of the following functions: `pkmod1cpt`, `pkmod2cpt`, `pkmod3cpt`,
#' or `pkmod3cptm` based on the parameter names entered.
#' @examples
#' # 1-compartment
#' infer_pkfn(c("CL","V"))
#' infer_pkfn(c("Cl","v1"))
#' # 2-compartment
#' infer_pkfn(c("CL","v","v2","q"))
#' # 3-compartment
#' infer_pkfn(c("CL","v","v2","q","Q2","V3"))
#' # 3-compartment with effect-site
#' infer_pkfn(c("CL","v","v2","q","Q2","V3","ke0"))
#' @export
infer_pkfn <- function(parnms){

  parnms <- tolower(parnms)

  if("ke0" %in% parnms)
    return(pkmod3cptm)

  if(any(c("v3","q3") %in% parnms))
    return(pkmod3cpt)

  if(any(c("v2","q2","q") %in% parnms))
    return(pkmod2cpt)

  if(any(c("cl","v","v1") %in% parnms))
    return(pkmod1cpt)

  warning("Could not infer pkfn from parameter names.")
}



#' @name list_parnms
#' @title Identify pkfn from parameter names
#' @description Identify structural PK model function (i.e., `pkfn`) from parameter names.
#' Models available are 1-, 2-, and 3-compartment mammillary models, or 3-compartment with
#' an effect site, corresponding to functions `pkmod1cpt`, `pkmod2cpt`, `pkmod3cpt`, and
#' `pkmod3cptm`, respectively.
#' @return Returns one of the following functions: `pkmod1cpt`, `pkmod2cpt`, `pkmod3cpt`,
#' or `pkmod3cptm` based on the parameter names entered.
#' @examples
#' list_parnms()
#' @export
list_parnms <- function(){
  cat("Acceptable names for 'pars_pk' vector (case-insensitive)", "\n")
  cat("--- First compartment options-------------", "\n")
  cat("Central volume: 'v','v1'","\n")
  cat("Elimination: 'cl','cl1','k10','ke'","\n")
  cat("--- Second compartment options------------", "\n")
  cat("Peripheral volume: 'v2'","\n")
  cat("Transfer: 'q','q2','k12','k21'","\n")
  cat("Elimination: 'cl2','k20'","\n")
  cat("--- Third compartment options-------------", "\n")
  cat("Second peripheral volume: 'v3'","\n")
  cat("Transfer: 'q3','k13','k31'","\n")
  cat("Elimination: 'cl3','k30'","\n")
  cat("--- Effect-site---------------------------", "\n")
  cat("Elimination: 'ke0'","\n")
}


#' @name list_pkmods
#' @title Print population PK models available in `tci`
#' @description Print population PK models available in `tci` for propofol (Marsh,
#' Schnider, Eleveld) and remifentanil (Minto, Kim, Eleveld).
#' @return Prints function names, model types, and required covariates for each
#' model.
#' @examples
#' list_pkmods()
#' @import knitr
#' @export
list_pkmods <- function(){
  tab <- data.frame(`Population model` = c("Marsh","Schnider","Eleveld (propofol)",
                              "Minto","Kim","Eleveld (remifentanil"),
                    Function = c("pkmod_marsh()","pkmod_schnider()","pkmod_eleveld_ppf()",
                                 "pkmod_minto()","pkmod_kim()","pkmod_eleveld_remi()"),
                    Drug = rep(c("Propofol","Remifentanil"), each = 3),
                    Type = c("PK","PK","PK/PKPD","PK/PKPD","PK","PK/PKPD"),
                    `Required covariates` = c("TBW","AGE, HGT, LBM or (TBW and MALE)",
                                   "AGE, HGT, MALE, TBW",
                                   "AGE, LBM or (MALE, TBW, and HGT)",
                                   "AGE, TBW, FFM or (MALE and BMI)",
                                   "AGE, MALE, BMI or (TBW and HGT)"))

  knitr::kable(tab, format = "pipe")
}


# #' @name restrict_sigmoid
# #' @title Restrict target sigmoid values
# #' @description Function to place restriction on gamma and E50 parameters of target sigmoid
# #' such that it passes through point (tfinal, BISfinal+eps)
# #'
# #' @param t50 parameter of Emax model
# #' @param tfinal end of the induction period
# #' @param eps distance between BISfinal and the target function at tfinal
# #' @param BIS0 starting BIS value
# #' @param BISfinal asymptote of Emax model
# #' @return Numeric vector of PD parameter values
# #' @export
# restrict_sigmoid <- function(t50, tfinal =10, eps = 1, BIS0 = 100, BISfinal = 50-eps){
#   gamma <- log((BIS0-BISfinal)/eps - 1, base = tfinal/t50)
#   c(c50 = t50, gamma = gamma, e0 = BIS0, emx = BIS0 - BISfinal)
# }


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
#' @return List of parameters used by Eleveld PK-PD model.
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
#' @return pkmod object
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



