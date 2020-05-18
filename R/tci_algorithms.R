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
#' @param tci TCI algorithm. Options are provided for effect-site (default) or plasma targeting.
#' Alternate algorithms can be specified through the 'tci_custom' argument.
#' @param
#'
iterate_tci_grid <- function(Ct, tms, pkmod, pars, init = NULL,
                             tci = c("effect","plasma"),
                             dt = 1/6, tci_custom = NULL, ...){


  tci <- match.arg(tci)

  if(!is.null(tci_custom)){
    tci <- tci_custom
  } else{
    if(tci == "effect") tci <- tci_comb
    else tci <- tci_plasma
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
    inf[i,] <- tci(sf(updatetms[i]), pkmod = pkmod, pars = pars, dt = dt, init = ini[,i], ...)
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




#' #' Function for calculating a TCI infusion schedule corresponding to a set of target concentrations. Primary inputs are a vector
#' #' of target effect-site concentrations and a set of parameters for a 3-compartment PK model. Initial concentrations (init)
#' #' can be set, as can the frequency with which the TCI device updates the target (delta), the maximum infusion rate, and
#' #' the sensitivity associated with switching to plasma-targeting mode (plasma_tol, effect_tol).
#' #'
#' #' This function makes use of formulas described by Shafer and Gregg (1992) in "Algorithms to rapidly achieve
#' #' and maintain stable drug concentrations at the site of drug effect with a computer-controlled infusion pump", as well as
#' #' the plasma-targeting algorithm described by Jacobs (1990) in "Algorithm for Optimal Linear Model-Based Control with
#' #' Application to Pharmacokinetic Model-Driven Drug Delivery."
#' #'
#' #' @param Cet Numeric vector of target effect-site concentrations. It is assumed that these occur delta seconds/min
#' #' apart where delta is the frequency with which the algorithm adjusts the target.
#' #' @param pars Named vector of parameters for a 3-compartment model with effect-site.
#' #' @param init Inital concentrations for the 4 compartments.
#' #' @param delta Frequency of TCI updates. Default is 10 seconds, expressed in terms of minutes = 1/6.
#' #' @param max_kR Maximum infusion rate of TCI pump.
#' #' @param plasma_tol Maximum percent difference between predicted plasma concentration and target concentration permitted
#' #' in order to switch to plasma-targeting mode. Plasma-targeting mode is used primarily to increase computational efficiency
#' #' when the effect-site concentration is sufficiently stable and close to the target concentration.
#' #' @param effect_tol Maximum percent difference between predicted effect-site concentration and target concentration
#' #' permitted in order to switch to plasma-targeting mode.
#'
#' TCI_EffectSite <- function(Cet, pars, delta = 1/6, init = c(0,0,0,0), max_kR = 10000, plasma_tol = 0.1, effect_tol = 0.005, x1 = 0.1, x2 = 0.2, tmax_search = 4, grid_space = 1/120) {
#'   infrt <- rep(NA, length(Cet))
#'   begin <- seq(0,(length(Cet)-1)*delta, delta)
#'   end <- seq(delta,length(Cet)*delta, delta)
#'   con <- matrix(NA, nrow = length(Cet)+1, ncol = 4)
#'   con[1,] <- init
#'
#'   for(i in 1:length(Cet)){
#'     con[con[i,]<0] <- 0
#'
#'     # if target is below current concentration, stop infusion
#'     if(Cet[i] < con[i,4]){
#'       infrt[i] <- 0
#'     }
#'
#'     # plasma targeting
#'     plasma_cond <- (con[i,1] > Cet[i] - plasma_tol*Cet[i]) &  (con[i,1] < Cet[i] + plasma_tol*Cet[i]) &
#'       (con[i,4] > Cet[i] - effect_tol*Cet[i]) & (con[i,4] < Cet[i] + effect_tol*Cet[i])
#'
#'     if(plasma_cond & is.na(infrt[i])){
#'       sol_x1 <- pk_basic_solution_3cpt_metab_c(k10 = pars[1],
#'                                                k12 = pars[2],
#'                                                k21 = pars[3],
#'                                                k13 = pars[4],
#'                                                k31 = pars[5],
#'                                                v1  = pars[6],
#'                                                v2  = pars[7],
#'                                                v3  = pars[8],
#'                                                ke0 = pars[9],
#'                                                kR  = x1, c0 = con[i,])
#'       sol_x2 <- pk_basic_solution_3cpt_metab_c(k10 = pars[1],
#'                                                k12 = pars[2],
#'                                                k21 = pars[3],
#'                                                k13 = pars[4],
#'                                                k31 = pars[5],
#'                                                v1  = pars[6],
#'                                                v2  = pars[7],
#'                                                v3  = pars[8],
#'                                                ke0 = pars[9],
#'                                                kR  = x2, c0 = con[i,])
#'       Cp_x1 <- sol_x1[[1]](delta)
#'       Cp_x2 <- sol_x2[[1]](delta)
#'       m <- (Cp_x2 - Cp_x1) / (x2-x1)
#'       b <- Cp_x1 - m*x1
#'       infrt_try <- (Cet[i] - b) / m
#'       infrt[i] <- ifelse(infrt_try < 0, 0, infrt_try)
#'       if(infrt[i] > max_kR) infrt[i] <- max_kR
#'     }
#'     if(is.na(infrt[i])){
#'       infrt[i] <- TCI_basic(Ce = Cet[i], pars = pars, init = con[i,], tmax_search = tmax_search, grid_space = grid_space)
#'     }
#'
#'     sol_x0 <- pk_basic_solution_3cpt_metab_c(k10 = pars[1],
#'                                              k12 = pars[2],
#'                                              k21 = pars[3],
#'                                              k13 = pars[4],
#'                                              k31 = pars[5],
#'                                              v1  = pars[6],
#'                                              v2  = pars[7],
#'                                              v3  = pars[8],
#'                                              ke0 = pars[9],
#'                                              kR  = infrt[i], c0 = con[i,])
#'     con[i+1,] <- c(sol_x0$c_1(delta), sol_x0$c_2(delta), sol_x0$c_3(delta), sol_x0$c_4(delta))
#'   }
#'
#'   if(any(is.nan(infrt))) stop("NaN in infusion rate")
#'   ivt_out <- lapply(1:length(Cet), function(j) list(begin = begin[j], end = end[j], k_R = infrt[j]))
#'   return(ivt_out)
#' }
#'
#'
#'
#'
#'
#' #' @param inf_duration Duration of infusion to be returned (minutes)
#' TCI_basic <- function(Ce, pars, init, inf_duration = 1/6, tmax_search = inf_duration*120, grid_space = inf_duration/10, return_jpeak = F){
#'   # course without infusion - use current concentration
#'   B <- pk_solution_3cpt_metab_singleinf(pars = pars, ivt = list(begin = 0, end = inf_duration, k_R = 0), init = init, ce_only = T)
#'   # course with infusion starting from 0 concentration
#'   E <- pk_solution_3cpt_metab_singleinf(pars = pars, ivt = list(begin = 0, end = inf_duration, k_R = 1), init = c(0,0,0,0), ce_only = T)
#'   # find tpeak
#'   tm_seq <- seq(grid_space,tmax_search,grid_space)
#'   ceproj <- E(tm_seq)
#'   tpeak = tm_seq[which.max(ceproj)]
#'
#'   if(all(init == 0)){
#'     kR <- unname((Ce-B(tpeak)) / E(tpeak))
#'     if(return_jpeak) return(c(kR = kR, jpeak = tpeak, Ce = E(tpeak)))
#'     else return(kR)
#'   } else{
#'     tm_seq2 <- seq(grid_space,tpeak+0.5,grid_space)
#'     jpeak0 = tpeak - 0.1# initialize search
#'     jpeak1 = jpeak0 + 0.1
#'     while(jpeak0 != jpeak1){
#'       jpeak0 = jpeak1
#'       I0 = (Ce - B(jpeak0)) / E(jpeak0)
#'       ceproj = B(tm_seq2) + E(tm_seq2)*I0
#'       jpeak1 = tm_seq2[which.max(ceproj)]
#'     }
#'
#'     kR = unname((Ce-B(jpeak1)) / E(jpeak1))
#'     if(kR < 0) kR = 0
#'     if(return_jpeak) return(c(kR = kR, jpeak = jpeak1, Ce = B(jpeak1) + E(jpeak1)*I0))
#'     else return(kR)
#'   }
#' }
#'
#' #' @examples
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
#' Hinv(pars = pars_pd_current, bis = bis_target, E0 = E0, gamma = gamma_eval)
#' ce_target <- Hinv(pars = pars["c50"], bis = 60, E0 = pars["E0"], gamma = pars["gamma"])
#' TCI_basic(Ce = ce_target, pars = pars, init = c(0,0,0,0))





