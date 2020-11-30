# --------------------------------------------------------------------------------------------------------------------------------
# - TCI algorithms ---------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------

#' TCI algorithm for plasma targeting
#'
#' TCI algorithm based on the algorithm described by Jacobs (1990).
#'
#' @param Cpt Target plasma concentration
#' @param pkmod PK model
#' @param dt Duration of the infusion
#' @param maxrt Maximum infusion rate. Defaults to 200 ml/min in reference to the
#' maximum infusion rate of 1200 ml/h permitted by
#' existing TCI pumps (e.g. Anestfusor TCI program).
#' @param cmpt Compartment into which infusions are administered. Defaults to the first compartment.
#' @param ... Arguments passed on to pkmod.
#'
#' @export
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
  return(c(kR = infrt))
}



#' TCI algorithm for effect-site targeting
#'
#' Function for calculating a TCI infusion schedule corresponding to a set of target concentrations.
#' This function makes use of formulas described by Shafer and Gregg (1992) in "Algorithms to rapidly achieve
#' and maintain stable drug concentrations at the site of drug effect with a computer-controlled infusion pump"
#'
#' @param Cet Numeric vector of target effect-site concentrations.
#' @param pkmod PK model
#' @param dt Frequency of TCI updates. Default is 1/6 minutes = 10 seconds.
#' @param ecmpt Effect site compartment number
#' @param tmax_search Outer bound on times searched to find a maximum concentration
#' following an infusion of duration dt. Defaults to 20 minutes. May need to be increased
#' if a drug has a slow elimination rate.
#' @param maxrt Maximum infusion rate of TCI pump. Defaults to 1200.
#' @param grid_len Number of time points used to identify time of maximum concentration.
#' Can be increased for more precision.
#' @param ... Arguments used by pkmod.
#' @export
tci_effect <- function(Cet, pkmod, dt = 1/6, ecmpt = NULL, tmax_search = 20,
                       maxrt = 1200, grid_len = 1200, ...){

  list2env(list(...), envir = environment())
  if(is.null(init)) init <- eval(formals(pkmod)$init)
  if(is.null(pars)) pars <- try(eval(formals(pkmod)$pars),
                                silent = TRUE)
  if(is.null(ecmpt)) ecmpt <- length(init)
  if(class(pars) == "try-error") stop("PK parameters must either be provided as arguments to the TCI algorithm or as defaults to the PK model.")

  ecmpt_name <- paste0("c",ecmpt)

  # infusions corresponding to unit infusion for duration and a null infusion
  unit_inf <- create_intvl(data.frame(time = c(dt, tmax_search), infrt = c(1,0)))
  null_inf <- create_intvl(data.frame(time = tmax_search, infrt = 0))

  # predict concentrations with no additional infusions
  B <- function(tm)
    predict(pkmod, inf = null_inf, pars = pars, init = init, tms = tm)[,ecmpt_name]

  # predict concentrations with no additional infusions
  E <- function(tm)
    predict(pkmod, inf = unit_inf, pars = pars, init = rep(0,length(init)), tms = tm)[,ecmpt_name]

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
    iter = 0

    while(jpeak0 != jpeak1){
      if(iter > 100) stop("Effect-site TCI algorithm did not converge.")
      jpeak0 = jpeak1
      I0 = (Cet - B(jpeak0)) / E(jpeak0)
      ceproj = B(tms) + E(tms)*I0
      jpeak1 = tms[which.max(ceproj)]
      iter = iter + 1
    }

    kR = unname((Cet-B(jpeak1)) / E(jpeak1))
  }

  if(kR < 0) kR = 0
  if(kR > maxrt) kR = maxrt

  return(c(kR = kR))
}



#' Effect-site TCI algorithm with plasma targeting within small range of target
#'
#' Modified effect-site TCI algorithm that follows Jacobs (1993) suggestion of
#' switching to plasma-targeting when the plasma concentration is within 10\%
#' of the target and the effect-site concentration is within 0.5\% of the target.
#' The modification decreases computation time and prevents oscillatory behavior
#' in the effect-site concentrations.
#'
#' @param Ct Numeric vector of target effect-site concentrations.
#' @param pkmod PK model
#' @param cptol Percentage of plasma concentration required to be within to switch
#' to plasma targeting.
#' @param cetol Percentage of effect-site concentration required to be within to switch
#' to plasma targeting.
#' @param cp_cmpt Position of central compartment. Defaults to first compartment.
#' @param ce_cmpt Position of effect-site compartment. Defaults to fourth compartment.
#' @param ... Arguments passed on to 'tci_plasma' and 'tci_effect' functions.
#'
#' @export
tci_comb <- function(Ct, pkmod, cptol = 0.1, cetol = 0.05, cp_cmpt = 1, ce_cmpt = 4, ...){

  list2env(list(...), envir = environment())

  if(is.null(init)) init <- eval(formals(pkmod)$init)

  if(Ct <= init[ce_cmpt])
    return(0)

  if(Ct>0){
    if(abs((Ct-init[cp_cmpt]) / Ct) <= cptol & abs((Ct-init[ce_cmpt]) / Ct) <= cetol){
      tci_plasma(Cpt = Ct, pkmod = pkmod, ...)
    } else{
      tci_effect(Cet = Ct, pkmod = pkmod, ...)
    }
  } else 0

}

#' Apply TCI algorithm
#'
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
#' @param tms Times at which the TCI algorithm should try to achieve the
#' target concentrations
#' @param pkmod PK model
#' @param pars PK model parameters
#' @param init Initial concentrations for PK model
#' @param tci_alg TCI algorithm. Options are provided for effect-site
#' (default) or plasma targeting. Alternate algorithms can be specified
#' through the 'tci_custom' argument.
#' @param tci_custom Custom TCI algorithm. Algorithm should have arguments
#' specifying target concentration, PK model, and duration of infusion to
#' reach the target.
#' @param dt Time difference between infusion rate updates.
#' @param ... Arguments passed on to TCI algorithm.
#' @export
tci <- function(Ct, tms, pkmod, pars, init = NULL,
                             tci_alg = c("effect","plasma"),
                             tci_custom = NULL, dt = 1/6, ...){


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
  updatetms <- seq(dt, max(tms), dt)
  # add epsilon so that sf evaluates the final step
  updatetms[length(updatetms)] <- updatetms[length(updatetms)] + 1e-5

  ncpt <- length(eval(formals(pkmod)$init))
  if(is.null(init)) init <- rep(0,ncpt)
  inf <- rep(NA, length(updatetms))
  ini <- matrix(NA, nrow = ncpt, ncol = length(updatetms)+1)
  ini[,1] <- init

  # iterate through times
  for(i in 1:length(updatetms)){
    inf[i] <- tci_alg(sf(updatetms[i]), pkmod = pkmod, pars = pars,
                      dt = dt, init = ini[,i], ...)
    ini[,i+1] <- pkmod(tm = dt, kR = inf[i], pars = pars, init = ini[,i])
  }

  startcon <- matrix(ini[,-ncol(ini)], ncol = nrow(ini), nrow = ncol(ini)-1, byrow = T)
  endcon <- matrix(ini[,-1], ncol = nrow(ini), nrow = ncol(ini)-1, byrow = T)
  dose <- create_intvl(cbind(time = seq(dt+inittm, max(tms)+inittm, dt), infrt = inf), inittm = inittm)
  out <- cbind(dose, dose[,"end"] - dose[,"begin"], sf(updatetms), startcon, endcon)
  colnames(out) <- c("infrt","begin","end","dt","Ct",paste0("c",1:ncpt, "_start"), paste0("c",1:ncpt, "_end"))
  class(out) <- c("tciinf",class(out))

  return(out)
}



#' Function to extend TCI grid to a set of PD targets
#'
#' @param pdresp PD targets to be passed on to the TCI algorithm.
#' @param tms Times corresponding to each PD target
#' @param pkmod PK model function
#' @param pdmod PD model function
#' @param pars_pk PK model parameters
#' @param pars_pd PD model parameters
#' @param pdinv PD inverse function
#' @param ecmpt Number corresponding to effect-site compartment. Defaults
#' to the last compartment.
#' @param ... Arguments to be passed on to 'tci'. These can include alternate
#' TCI algorithms if desired.
#'
#' @export
tci_pd <- function(pdresp, tms, pkmod, pdmod, pars_pk, pars_pd, pdinv,
                   ecmpt = NULL, ...){
  Ct <- pdinv(pdresp, pars_pd)
  con <- tci(Ct = Ct, tms = tms, pkmod = pkmod, pars = pars_pk, ...)
  con_class <- class(con)

  if(is.null(ecmpt))
    ecmpt <- length(eval(formals(pkmod)$init))

  pdt <- pdmod(con[,"Ct"], pars_pd)
  pdresp_start <- pdmod(con[,paste0("c",ecmpt,"_start")], pars_pd)
  pdresp_end <- pdmod(con[,paste0("c",ecmpt,"_end")], pars_pd)

  con <- cbind(con,
               pdt = pdt,
               pdresp_start = pdresp_start,
               pdresp_end = pdresp_end)
  class(con) <- con_class
  return(con)
}


