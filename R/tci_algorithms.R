# --------------------------------------------------------------------------------------------------------------------------------
# - TCI algorithms ---------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------

# #' TCI algorithm for plasma targeting
# #'
# #' TCI algorithm based on the algorithm described by Jacobs (1990).
# #'
# #' @param Cpt Target plasma concentration
# #' @param pkmod PK model
# #' @param dtm Duration of the infusion
# #' @param maxrt Maximum infusion rate. Defaults to 200 ml/min in reference to the
# #' maximum infusion rate of 1200 ml/h permitted by
# #' existing TCI pumps (e.g. Anestfusor TCI program).
# #' @param ... Arguments passed on to pkmod.
# #' @return Numeric value
# #' @examples
# #' my_mod <- pkmod("1", pars_pk = c(CL = 10, V1 = 10))
# #' tci_plasma(Ct = 2, pkmod = my_mod, dtm = 1)
# #' # update CL parameter
# #' tci_plasma(Ct = 2, pkmod = my_mod, dtm = 1, pars = c(CL = 15))
# #' @export
# tci_plasma <- function(Ct, pkmod, dtm, maxrt = 1200, ...){
#
#   pkmod <- update(pkmod,...)
#
#   Cp1 <- pkmod$pkfn(tm = dtm, kR = 1, ...)
#   Cp2 <- pkmod$pkfn(tm = dtm, kR = 2, ...)
#
#   # for multi-compartment models, use only concentration in compartment 'cmpt'
#   if(!is.null(dim(Cp1))){
#     Cp1 <- Cp1[pkmod$pcmpt,]
#     Cp2 <- Cp2[pkmod$pcmpt,]
#   }
#
#   m <- Cp2 - Cp1
#   b <- Cp1 - m
#   infrt <- (Cpt - b) / m
#   if(infrt < 0)
#     infrt <- 0
#   if(infrt > maxrt)
#     infrt <- maxrt
#   return(unname(infrt))
# }


#' TCI algorithm for plasma targeting
#'
#' TCI algorithm based on the algorithm described by Jacobs (1990).
#'
#' @param Ct Target plasma concentration
#' @param pkmod PK model
#' @param dtm Duration of the infusion
#' @param maxrt Maximum infusion rate. Defaults to 200 ml/min in reference to the
#' maximum infusion rate of 1200 ml/h permitted by
#' existing TCI pumps (e.g. Anestfusor TCI program).
#' @param ... Arguments passed on to update.pkmod.
#' @return Numeric value
#' @examples
#' # plasma targeting
#' my_mod <- pkmod(pars_pk = c(CL = 10, V1 = 10))
#' tci_plasma(Ct = 2, my_mod)
#' # update CL parameter
#' tci_plasma(Ct = 2, my_mod, pars_pk = c(CL = 15))
#' @export
tci_plasma <- function(Ct, pkmod, dtm = 1/6, maxrt = 1200, ...){

  pkmod <- update(pkmod,...)

  Cp1 <- with(pkmod, pkfn(tm = dtm, kR = 1, pars = pars_pk, init = init))
  Cp2 <- with(pkmod, pkfn(tm = dtm, kR = 2, pars = pars_pk, init = init))

  if(!is.null(dim(Cp1))){
    Cp1 <- Cp1[pkmod$pcmpt,]
    Cp2 <- Cp2[pkmod$pcmpt,]
  }

  m <- Cp2 - Cp1
  b <- Cp1 - m
  infrt <- (Ct - b) / m
  if(infrt < 0)
    infrt <- 0
  if(infrt > maxrt)
    infrt <- maxrt
  return(unname(infrt))
}


#' TCI algorithm for effect-site targeting
#'
#' Function for calculating a TCI infusion schedule corresponding to a set of target concentrations.
#' This function makes use of formulas described by Shafer and Gregg (1992) in "Algorithms to rapidly achieve
#' and maintain stable drug concentrations at the site of drug effect with a computer-controlled infusion pump"
#'
#' @param Ct Numeric vector of target effect-site concentrations.
#' @param pkmod PK model
#' @param dtm Frequency of TCI updates. Default is 1/6 minutes = 10 seconds.
#' @param tmax_search Outer bound on times searched to find a maximum concentration
#' following an infusion of duration dtm. Defaults to 20 minutes. May need to be increased
#' if a drug has a slow elimination rate.
#' @param maxrt Maximum infusion rate of TCI pump. Defaults to 1200.
#' @param grid_len Number of time points used to identify time of maximum concentration.
#' Can be increased for more precision.
#' @param ... List or vector of named values to be passed on to update.pkmod.
#' @return Numeric value
#' @examples
#' # three compartment model with effect site
#' my_mod <- pkmod(pars_pk = c(v1 = 8.995, v2 = 17.297, v3 = 120.963,
#' cl = 1.382, q2 = 0.919, q3 = 0.609, ke0 = 1.289))
#' tci_effect_only(Ct = 2, pkmod = my_mod)
#' # update parameters
#' tci_effect_only(Ct = 2, pkmod = my_mod, pars_pk = c(v1 = 12, cl = 2))
#' @export
tci_effect_only <- function(Ct, pkmod, dtm = 1/6, tmax_search = 10,
                       maxrt = 1200, grid_len = 1200, ...){

  # updates to pkmod
  pkmod <- update(pkmod, ...)

  # check that effect site compartment is specified
  if(is.null(pkmod$ecmpt)) stop("ecmpt must be specified in new_pkmod to use tci_effect")

  # infusions corresponding to unit infusion for duration dtm and a null infusion
  unit_inf <- create_inf(times = c(0,dtm,tmax_search), infrt = c(1,0,0))
  null_inf <- create_inf(times = c(0,tmax_search), infrt = c(0,0))

  # predict concentrations with no additional infusions and starting concentrations
  B <- function(tm) predict(pkmod, inf = null_inf, tms = tm)[,pkmod$ecmpt]

  # predict concentrations with unit infusion and no starting concentrations
  E <- function(tm) predict(pkmod, inf = unit_inf, init = rep(0,pkmod$ncmpt), tms = tm)[,pkmod$ecmpt]

  # predict to find the longest time of maximum concentration
  # this will always be shorter when any prior drug has been infused
  grid_tmax <- seq(0,tmax_search,length.out = grid_len)
  con_proj <- E(grid_tmax)
  con_dif <- diff(con_proj)
  while(all(con_dif > 0)){
    tmax_search <- tmax_search*2
    grid_tmax <- seq(0,tmax_search,length.out = grid_len)
    con_proj <- E(grid_tmax)
    con_dif <- diff(con_proj)
  }

  peak_ix <- min(which(diff(con_proj) < 0))

  if(all(pkmod$init == 0)){
    kR <- Ct / con_proj[peak_ix]
  } else{
    tms <- seq(0, grid_tmax[peak_ix]+0.5, length.out = grid_len)
    Bpred <- B(tms)
    Epred <- E(tms)
    peak_ix <- which.max(Epred)
    jpeak1 <- tms[peak_ix]
    jpeak0 <- tms[peak_ix-1]
    iter = 0

    while(jpeak0 != jpeak1){
      if(iter > 100) stop("Effect-site TCI algorithm did not converge.")
      jpeak0 = jpeak1
      I0 = (Ct - Bpred[which(tms == jpeak0)]) / Epred[which(tms == jpeak0)]
      jpeak1 = tms[which.max(Bpred + Epred*I0)]
      iter = iter + 1
    }
    kR = unname((Ct-Bpred[which(tms == jpeak1)]) / Epred[which(tms == jpeak1)])
  }

  if(kR < 0) kR = 0
  if(kR > maxrt) kR = maxrt

  return(unname(kR))
}


# #' Effect-site TCI algorithm with plasma targeting within small range of target
# #'
# #' Modified effect-site TCI algorithm that switches to plasma-targeting when the
# #' plasma concentration is within 20\%
# #' of the target and the effect-site concentration is within 0.5\% of the target.
# #' The modification decreases computation time and prevents oscillatory behavior
# #' in the effect-site concentrations.
# #'
# #' @param Ct Numeric vector of target effect-site concentrations.
# #' @param pkmod PK model
# #' @param cptol Percentage of plasma concentration required to be within to switch
# #' to plasma targeting.
# #' @param cetol Percentage of effect-site concentration required to be within to switch
# #' to plasma targeting.
# #' @param cp_cmpt Position of central compartment. Defaults to first compartment.
# #' @param ce_cmpt Position of effect-site compartment. Defaults to last compartment.
# #' @param ... Arguments passed on to 'tci_plasma' and 'tci_effect' functions.
# #' @return Numeric value
# #' @examples
# #' tci_comb(Ct = 2, pkmod = pkmod2cpt, dtm = 1, pars = c(CL = 15, V1 = 10, Q2 = 10, V2 = 20))
# #' @export
# tci_comb <- function(Ct, pkmod, cptol = 0.2, cetol = 0.05, cp_cmpt = NULL, ce_cmpt = NULL, ...){
#
#   list2env(list(...), envir = environment())
#
#   if(!("init" %in% ls())) init <- eval(formals(pkmod)$init)
#   if(is.null(cp_cmpt)) cp_cmpt <- 1
#   if(is.null(ce_cmpt)) ce_cmpt <- length(init)
#
#   if(Ct <= init[ce_cmpt])
#     return(0)
#
#   if(Ct>0){
#     if(abs((Ct-init[cp_cmpt]) / Ct) <= cptol & abs((Ct-init[ce_cmpt]) / Ct) <= cetol){
#       tci_plasma(Cpt = Ct, pkmod = pkmod, ...)
#     } else{
#       tci_effect(Cet = Ct, pkmod = pkmod, ...)
#     }
#   } else 0
#
# }

#' Effect-site TCI algorithm with plasma targeting within small range of target
#'
#' Modified effect-site TCI algorithm that switches to plasma-targeting when the
#' plasma concentration is within 20\%
#' of the target and the effect-site concentration is within 0.5\% of the target.
#' The modification decreases computation time and prevents oscillatory behavior
#' in the effect-site concentrations.
#'
#' @param Ct Numeric vector of target effect-site concentrations.
#' @param pkmod PK model
#' @param dtm TCI update frequency. Defaults to 1/6, corresponding to 10-second
#' intervals if model parameters are in terms of minutes.
#' @param cptol Percentage of plasma concentration required to be within to switch
#' to plasma targeting.
#' @param cetol Percentage of effect-site concentration required to be within to switch
#' to plasma targeting.
#' @param ... Arguments passed on to 'tci_plasma' and 'tci_effect_only' functions, including
#' to update.pkmod.
#' @return Numeric value
#' @examples
#' my_mod <- pkmod(pars_pk = c(v1 = 8.995, v2 = 17.297, v3 = 120.963, cl = 1.382,
#' q2 = 0.919, q3 = 0.609, ke0 = 1.289))
#' tci_effect(Ct = 2, pkmod = my_mod)
#' # update parameters
#' tci_effect(Ct = 2, pkmod = my_mod, pars_pk = c(v1 = 12, cl = 2))
#' @export
tci_effect <- function(Ct, pkmod, dtm = 1/6, cptol = 0.2, cetol = 0.05, ...){

  pkmod <- update(pkmod, ...)

  # return zero if concentration exceeds target
  if(Ct <= pkmod$init[pkmod$ecmpt])
    return(0)

  if(Ct>0){
    if(abs((Ct-pkmod$init[pkmod$pcmpt]) / Ct) <= cptol & abs((Ct-pkmod$init[pkmod$ecmpt]) / Ct) <= cetol){
      tci_plasma(Ct = Ct, pkmod = pkmod, dtm = dtm)
    } else{
      tci_effect_only(Ct = Ct, pkmod = pkmod, dtm = dtm)
    }
  } else 0
}


#' Apply TCI algorithm
#'
#' Apply a TCI algorithm to a set of targets and a `pkmod` object to calculate infusion
#' rates.
#'
#' Predict method for tciinf objects.
#'
#' Apply a TCI algorithm and pkmod to a set of sequential targets.
#' @param targets A matrix or data frame with columns 'value' and 'time'. Times
#' indicate when the TCI algorithm should begin infusions to reach each target.
#' @param pkmod `pkmod` object created by `pkmod()`.
#' @param type Type of TCI algorithm to be used. Options are plasma- or effect-site
#' targeting.
#' @param dtm TCI update frequency. Defaults to 1/6, corresponding to 10-second
#' intervals if model parameters are in terms of minutes.
#' @param custom_alg Custom TCI algorithm to be used instead of default plasma-
#' or effect-site targeting algorithms. The algorithm should be a function that
#' takes minimum arguments `Ct`, `pkmod`, and `dtm` and returns a single infusion
#' rate. See `tci_plasma` or `tci_effect` for examples and vignette on custom
#' models/algorithms for more details.
#' @param inittm Initial time to start TCI algorithm. Cannot be greater than
#' the minimum value of `targets[,"time"]`.
#' @param ignore_pd Logical. Should the PD component of the pkmod object (if present)
#' be ignored. By default, predict.tciinf will assume that 'value' refers to PD
#' targets if a PD model is specified.
#' @param ... Arguments passed to TCI algorithm
#' @examples
#' # set of targets and times
#' targets = cbind(value = c(2,3,4,4), time = c(0,2,3,10))
#' # 3-compartment model with effect-site
#' my_mod <- pkmod(pars_pk = c(v1 = 8.995, v2 = 17.297, v3 = 120.963, cl = 1.382,
#' q2 = 0.919, q3 = 0.609, ke0 = 1.289))
#' # plasma targeting
#' apply_tci(targets, my_mod, "plasma")
#' # effect-site targeting
#' apply_tci(targets, my_mod, "effect")
#' # incorporate  PD model
#' my_mod_pd <- update(my_mod, pars_pd = c(c50 = 2.8, gamma = 1.47, e0 = 93, emx = 93),
#' pdfn = emax, pdinv = emax_inv)
#' pd_targets = cbind(value = c(70,60,50,50), time = c(0,2,3,10))
#' apply_tci(pd_targets, my_mod_pd, "effect")
#' @rdname apply_tci
#' @export
apply_tci <- function(targets, pkmod, type = c("plasma","effect"), dtm = NULL, custom_alg = NULL, inittm = 0, ignore_pd = FALSE, ...){

  if(!all(c("value","time") %in% colnames(targets)))
    stop("targets must have column names 'value' and 'time'")

  type <- match.arg(type)

  if(is.null(custom_alg)){
    tci_alg <- list(tci_plasma, tci_effect)[[which(type == c("plasma","effect"))]]
  } else{
    if(!all(c("Ct","pkmod","dtm") %in% names(formals(custom_alg))))
      stop('pkmod must contain arguments ("Ct","pkmod","dtm").')
    tci_alg <- custom_alg
  }

  if(is.null(dtm)) dtm <- eval(formals(tci_alg)$dtm)

  tms <- targets[,"time"]
  if(any(inittm > tms)) stop("inittm cannot be greater than any target times")

  if(!is.null(pkmod$pdfn) & !ignore_pd){
    Ct <- pkmod$pdinv(targets[,"value"], pars = pkmod$pars_pd)
  } else{
    Ct <- targets[,"value"]
  }

  # create step function to define targets at any point
  tms <- tms-inittm
  sf <- stepfun(tms, c(0,Ct))
  # updatetms <- seq(0, max(tms), dtm)
  updatetms <- seq(0, max(tms)-dtm, dtm)

  inf <- rep(NA, length(updatetms))
  ini <- matrix(NA, nrow = pkmod$ncmpt, ncol = length(updatetms)+1)
  ini[,1] <- pkmod$init

  # iterate through times
  for(i in 1:length(updatetms)){
    inf[i] <- tci_alg(sf(updatetms[i]), pkmod = pkmod, dtm = dtm, init = ini[,i], ...)
    ini[,i+1] <- pkmod$pkfn(tm = dtm, pars = pkmod$pars_pk, kR = inf[i], init = ini[,i])
  }

  startcon <- matrix(ini[,-ncol(ini)], ncol = nrow(ini), nrow = ncol(ini)-1, byrow = TRUE)
  endcon <- matrix(ini[,-1], ncol = nrow(ini), nrow = ncol(ini)-1, byrow = TRUE)
  dose <- create_inf(times = updatetms+inittm, infrt = inf, duration = dtm)
  out <- cbind(dose, sf(updatetms), startcon, endcon)
  colnames(out) <- c("begin","end","infrt","Ct",paste0("c",1:pkmod$ncmpt, "_start"),
                     paste0("c",1:pkmod$ncmpt, "_end"))

  if(!is.null(pkmod$pdfn) & !ignore_pd){
    pdt <- pkmod$pdfn(out[,"Ct"], pkmod$pars_pd)
    pdresp_start <- pkmod$pdfn(out[,paste0("c",pkmod$ecmpt,"_start")], pkmod$pars_pd)
    pdresp_end <- pkmod$pdfn(out[,paste0("c",pkmod$ecmpt,"_end")], pkmod$pars_pd)

    out <- cbind(out,
                 pdt = pdt,
                 pdresp_start = pdresp_start,
                 pdresp_end = pdresp_end)
  }

  return(out)
}




