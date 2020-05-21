# TCI algorithms

#' TCI algorithm for plasma targeting, based on the algorithm described by Jacobs (1990)
#'
#' @param Cpt Target plasma concentration
#' @param pkmod PK model
#' @param dt Duration of the infusion
#' @param maxrt Maximum infusion rate. Defaults to 1,200 in reference to the maximum infusion rate in ml/h permitted by
#' existing TCI pumps (e.g. Anestfusor TCI program).
#' @param cmpt Compartment into which infusions are administered. Defaults to the first compartment.
tci_plasma <- function(Cpt, pkmod, dt, maxrt = 1200, cmpt = 1, ...){

  Cp1 <- pkmod(tm = dt, kR = 1, ...)
  Cp2 <- pkmod(tm = dt, kR = 2, ...)

  # for multi-compartment models, use only concentration in compartment 'cmpt'
  if(!is.null(dim(Cp1))){
    Cp1 <- Cp1[cmpt,]
    Cp2 <- Cp2[cmpt,]
  }

  m <- Cp2 - Cp1
  b <- Cp1 - m
  infrt <- (Cpt - b) / m
  if(infrt < 0)
    infrt <- 0
  if(infrt > maxrt)
    infrt <- maxrt
  return(c(kR = infrt, dt = dt))
}



#' Function for calculating a TCI infusion schedule corresponding to a set of target concentrations.
#' This function makes use of formulas described by Shafer and Gregg (1992) in "Algorithms to rapidly achieve
#' and maintain stable drug concentrations at the site of drug effect with a computer-controlled infusion pump"
#'
#' @param Cet Numeric vector of target effect-site concentrations.
#' @param pkmod PK model
#' @param dt Frequency of TCI updates. Default is 10 seconds expressed in terms of minutes = 1/6.
#' @param max_kR Maximum infusion rate of TCI pump. Defaults to 1200.
#' @param cmpt Effect
#'
#'
#' @param plasma_tol Maximum percent difference between predicted plasma concentration and target concentration permitted
#' in order to switch to plasma-targeting mode. Plasma-targeting mode is used primarily to increase computational efficiency
#' when the effect-site concentration is sufficiently stable and close to the target concentration.
#' @param effect_tol Maximum percent difference between predicted effect-site concentration and target concentration
#' permitted in order to switch to plasma-targeting mode.
tci_effect <- function(Cet, pkmod, dt = 1/6, max_kR = 1200, cmpt = NULL, tmax_search = 20, maxrt = 200, grid_len = 1200, ...){

  list2env(list(...), envir = environment())
  if(is.null(init)) init <- eval(formals(pkmod)$init)
  if(is.null(pars)) pars <- try(eval(formals(pkmod)$pars), silent = T)
  if(is.null(cmpt)) cmpt <- length(init)
  if(class(pars) == "try-error") stop("PK parameters must either be provided as arguments to the TCI algorithm or as defaults to the PK model.")

  cmpt_name <- paste0("c",cmpt)

  # infusions corresponding to unit infusion for duration and a null infusion
  unit_inf <- create_intvl(data.frame(time = c(dt, tmax_search), infrt = c(1,0)))
  null_inf <- create_intvl(data.frame(time = tmax_search, infrt = 0))

  # predict concentrations with no additional infusions
  B <- function(tm)
    predict(pkmod, inf = null_inf, pars = pars, init = init, tms = tm)[,cmpt_name]

  # predict concentrations with no additional infusions
  E <- function(tm)
    predict(pkmod, inf = unit_inf, pars = pars, init = rep(0,length(init)), tms = tm)[,cmpt_name]

  # predict to find the longest time of maximum concentration -- will always be shorter when any prior drug has been infused.
  grid_tmax <- seq(0,tmax_search,length.out = grid_len)
  con_proj <- E(grid_tmax)
  peak_ix <- which.max(con_proj)
  con_peak <- con_proj[peak_ix]
  tpeak <- grid_tmax[peak_ix]

  if(all(init == 0)){
    kR <- Cet / con_peak
  } else{
    tms <- seq(0, tpeak+0.5, length.out = grid_len)
    jpeak0 = tpeak - 0.1
    jpeak1 = jpeak0 + 0.1

    while(jpeak0 != jpeak1){
      jpeak0 = jpeak1
      I0 = (Cet - B(jpeak0)) / E(jpeak0)
      ceproj = B(tms) + E(tms)*I0
      jpeak1 = tms[which.max(ceproj)]
    }

    kR = unname((Cet-B(jpeak1)) / E(jpeak1))
  }

  if(kR < 0) kR = 0
  if(kR > max_kR) kR = max_kR

  return(c(kR = kR, dt = dt))
}



#' Modified effect-site TCI algorithm that follows Jacobs (1993) suggestion of
#' switching to plasma-targeting when the plasma concentration is within 10\%
#' of the target and the effect-site concentration is within 0.5\% of the target.
#' The modification decreases computation time and prevents oscilatory behavior
#' in the effect-site concentrations.
tci_comb <- function(Ct, pkmod, cptol = 0.1, cetol = 0.05, cp_cmpt = 1, ce_cmpt = 4, ...){

  list2env(list(...), envir = environment())

  if(is.null(init)) init <- eval(formals(pkmod)$init)

  if(Ct>0){
    if(abs((Ct-init[cp_cmpt]) / Ct) <= cptol & abs((Ct-init[ce_cmpt]) / Ct) <= cetol){
      tci_plasma(Cpt = Ct, pkmod = pkmod, ...)
    } else{
      tci_effect(Cet = Ct, pkmod = pkmod, ...)
    }
  } else 0

}




#' Function to iterate any arbitrary TCI algorithm to a series of points. By default,
#' the function will update infusion rates at fixed intervals (e.g. every 10 seconds);
#' however, users will have the option of waiting only calculating infusions after
#' the prior target has been obtained.
#'
#' The user passes the `iterate_tci` function a matrix of target concentrations and times
#' at which the target is set. This is translated into a step function that defines the
#' concentration target at all times.
#'
#' @param Ct Vector of target concentrations
#' @param tms Times at which the TCI algorithm should try to achieve the target concentrations
#' @param tci_alg TCI algorithm. Options are provided for effect-site (default) or plasma targeting.
#' Alternate algorithms can be specified through the 'tci_custom' argument.
#' @param
#'
tci <- function(Ct, tms, pkmod, pars, init = NULL,
                             tci_alg = c("effect","plasma"),
                             dt = 1/6, tci_custom = NULL, ...){


  tci_alg <- match.arg(tci_alg)

  if(!is.null(tci_custom)){
    tci_alg <- tci_custom
  } else{
    if(tci_alg == "effect") tci_alg <- tci_comb
    else tci_alg <- tci_plasma
  }

  # adjust times such that infusions start at tm = 0 and can be usd by stepfun
  inittm <- tms[1]
  tms <- tms[-1] - inittm

  # create step function to define targets at any point
  sf <- stepfun(tms, Ct)

  # define sequence of update times
  # updatetms <- seq(dt, max(tms)-dt, dt)
  updatetms <- seq(dt, max(tms), dt)

  ncpt <- length(eval(formals(pkmod)$init))
  if(is.null(init)) init <- rep(0,ncpt)
  inf <- matrix(NA, nrow = length(updatetms), ncol = 2)
  ini <- matrix(NA, nrow = ncpt, ncol = length(updatetms)+1)
  ini[,1] <- init

  # iterate through times
  for(i in 1:length(updatetms)){
    inf[i,] <- tci_alg(sf(updatetms[i]), pkmod = pkmod, pars = pars, dt = dt, init = ini[,i], ...)
    ini[,i+1] <- pkmod(tm = dt, kR = inf[i,1], pars = pars, init = ini[,i])
  }
  startcon <- data.frame(matrix(ini[,-ncol(ini)], ncol = nrow(ini), nrow = ncol(ini)-1, byrow = T))
  endcon <- data.frame(matrix(ini[,-1], ncol = nrow(ini), nrow = ncol(ini)-1, byrow = T))

  out <- create_intvl(dose = data.frame(time = seq(dt, max(tms), dt), infrt = inf[,1]), inittm = inittm)
  out <- cbind(out, inf[,2], sf(updatetms), startcon, endcon)
  names(out) <- c("intvl","infrt","begin","end","dt","Ct",paste0("c",1:ncpt, "_start"), paste0("c",1:ncpt, "_end"))

  class(out) <- c("tciinf",class(out))
  return(out)
}



#' Function to extend TCI grid to a set of PD targets
tci_pd <- function(pdresp, tms, pkmod, pdmod, pars_pk, pars_pd, pdinv, ecmpt = NULL, ...){
  Ct <- pdinv(pdresp, pars_pd)
  con <- tci(Ct = Ct, tms = tms, pkmod = pkmod, pars = pars_pk, ...)
  if(is.null(ecmpt))
    ecmpt <- length(eval(formals(pkmod)$init))
  con$pdt <- pdmod(con$Ct, pars_pd)
  con$pdresp <- pdmod(con[,paste0("c",ecmpt,"_start")], pars_pd)
  return(con)
}


