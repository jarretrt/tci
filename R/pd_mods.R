# --------------------------------------------------------------------------------------------------------------------------------
# - Library of PD functions ------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------


#' @name emax
#' @title Emax function
#' @description Emax function. c50 is the concentration eliciting a 50% effect, gamma is the hill parameter
#' identifying the slope of the Emax curve at c50, E0 is the response value with no drug present,
#' Emx is the maximum effect size.
#' @param ce Vector of effect-site concentrations.
#' @param pars Named vector of parameter values with names (c50,gamma,e0,emx).
#' @return Numeric vector of same length as ce.
#' @examples
#' pars_emax <- c(c50 = 1.5, gamma = 1.47, e0 = 100, emx = 100)
#' ce_seq <- seq(0,4,0.1)
#' plot(ce_seq, emax(ce_seq, pars_emax), type = "l",
#' xlab = "Effect-site concentrtion (ug/mL)", ylab = "BIS")
#' @export
emax <- function(ce, pars){
  names(pars) <- tolower(names(pars))
  if("ce50" %in% names(pars)) names(pars)[names(pars) == "ce50"] <- "c50"
  pars["e0"] - pars["emx"]*(ce^pars["gamma"] / (ce^pars["gamma"] + pars["c50"]^pars["gamma"]))
}
class(emax) <- "pdmod"


#' @name emax_inv
#' @title Inverse Emax function
#' @description Inverse Emax function to return effect-site concentrations required to reach target effect.
#' @param pdresp PD response values
#' @param pars Named vector of parameter values with names (c50,gamma,E0,Emx).
#' @return Numeric vector of same length as pdresp.
#' @examples
#' pars_emax <- c(c50 = 1.5, gamma = 4, e0 = 100, emx = 100)
#' ce_seq <- seq(0,4,0.1)
#' all.equal(emax_inv(emax(ce_seq, pars_emax), pars_emax), ce_seq)
#' @export
emax_inv <- function(pdresp, pars){
  names(pars) <- tolower(names(pars))
  if("ce50" %in% names(pars)) names(pars)[names(pars) == "ce50"] <- "c50"
  eff <- abs(pdresp - pars["e0"])
  (eff*(pars["c50"]^pars["gamma"])/(pars["emx"]*(1-eff/pars["emx"])))^(1/pars["gamma"])
}

#' @name emax_eleveld
#' @title Emax function for Eleveld (2018) model.
#' @description The parameter gamma takes one of two values depending on whether ce <= c50.
#'
#' @param ce Vector of effect-site concentrations.
#' @param pars Vector of parameter values in order (c50,gamma,gamma2,e0,emx).
#' @return Numeric vector of same length as ce.
#' @examples
#' pars_emax_eleveld <- c(c50 = 1.5, e0 = 100, gamma = 1.47, gamma2 = 1.89)
#' ce_seq <- seq(0,4,0.1)
#' plot(ce_seq, emax_eleveld(ce_seq, pars_emax_eleveld), type = "l",
#' xlab = "Effect-site concentrtion (ug/mL)", ylab = "BIS")
#' @export
emax_eleveld <- function(ce, pars){
  names(pars) <- tolower(names(pars))
  if("ce50" %in% names(pars)) names(pars)[names(pars) == "ce50"] <- "c50"
  if("emax" %in% names(pars)) names(pars)[names(pars) == "emax"] <- "emx"
  c50    <- pars["c50"]
  e0     <- pars["e0"]
  gamma  <- pars["gamma"]
  gamma2 <- ifelse("gamma2" %in% names(pars), pars["gamma2"], pars["gamma"])
  emx    <- ifelse("emx" %in% names(pars), pars["emx"], pars["e0"])
  gam <- ifelse(ce <= c50, gamma, gamma2)

  e0 - emx*(ce^gam / (ce^gam + c50^gam))
}



#' @name emax_inv_eleveld
#' @title Inverse Emax function
#' @description Inverse of Emax function used by Eleveld population PK model.
#' @param pdresp PD response values
#' @param pars Named vector of parameter values with names (c50,gamma,E0,Emx).
#' @return Numeric vector of same length as pdresp.
#' @examples
#' pars_emax_eleveld <- c(c50 = 1.5, e0 = 100, gamma = 1.47, gamma2 = 1.89)
#' ce_seq <- seq(0,4,0.1)
#' all.equal(emax_inv_eleveld(emax_eleveld(ce_seq, pars_emax_eleveld), pars_emax_eleveld), ce_seq)
#' @export
emax_inv_eleveld <- function(pdresp, pars){

  names(pars) <- tolower(names(pars))
  if("ce50" %in% names(pars)) names(pars)[names(pars) == "ce50"] <- "c50"
  if("emax" %in% names(pars)) names(pars)[names(pars) == "emax"] <- "emx"
  if("emax" %in% names(pars)) names(pars)[names(pars) == "emax"] <- "emx"
  c50    <- pars["c50"]
  e0     <- pars["e0"]
  gamma  <- pars["gamma"]
  gamma2 <- ifelse("gamma2" %in% names(pars), pars["gamma2"], pars["gamma"])
  emx    <- ifelse("emx" %in% names(pars), pars["emx"], pars["e0"])

  mideff <- e0 - emx/2 # effect at c50
  eff <- abs(pdresp - e0)
  gam <- ifelse(pdresp >= mideff, gamma, gamma2)

  (eff*(c50^gam)/(emx*(1-eff/emx)))^(1/gam)
}



#' @name emax_remi
#' @title Emax function implemented by Eleveld remifentanil model
#' @description Emax function. c50 is the concentration eliciting a 50% effect, gamma is the hill parameter
#' identifying the slope of the Emax curve at c50, E0 is the response value with no drug present,
#' Emx is the maximum effect size.
#' @param ce Vector of effect-site concentrations.
#' @param pars Named vector of parameter values with names (c50,gamma,e0,emx).
#' @return Numeric vector of same length as ce.
#' @examples
#' pars_emax <- c(e0 = 19, emx = 5.6, c50 = 12.7, gamma = 2.87)
#' ce_seq <- seq(0,60,0.1)
#' plot(ce_seq, emax_remi(ce_seq, pars_emax), type = "l",
#' xlab = "Effect-site concentrtion (ug/mL)", ylab = "Spectral Edge Frequency (HZ)")
#' @export
emax_remi = function(ce,pars){
  names(pars) <- tolower(names(pars))
  if("ce50" %in% names(pars)) names(pars)[names(pars) == "ce50"] <- "c50"
  e0 = pars["e0"]
  emx = pars["emx"]
  c50 = pars["c50"]
  gam = pars["gamma"]
  e0 + (emx - e0)*(ce^gam / (ce^gam + c50^gam))
}


#' @name emax_inv_remi
#' @title Inverse Emax function implemented by Eleveld remifentanil model
#' @description Inverse Emax function to return effect-site concentrations required to reach target effect.
#' @param pdresp PD response values
#' @param pars Named vector of parameter values with names (c50,gamma,E0,Emx).
#' @return Numeric vector of same length as pdresp.
#' @examples
#' pars_emax <- c(e0 = 19, emx = 5.6, c50 = 12.7, gamma = 2.87)
#' emax_inv_remi(emax_remi(10, pars_emax), pars_emax)
#' @export
emax_inv_remi = function(pdresp, pars){
  names(pars) <- tolower(names(pars))
  if("ce50" %in% names(pars)) names(pars)[names(pars) == "ce50"] <- "c50"
  e0 = pars["e0"]
  emx = pars["emx"]
  c50 = pars["c50"]
  gam = pars["gamma"]
  c50*((emx-e0)/(pdresp-e0)-1)^(-1/gam)
}
