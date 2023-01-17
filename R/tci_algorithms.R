# --------------------------------------------------------------------------------------------------------------------------------
# - TCI algorithms ---------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------


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
  inf_rate <- (Ct - b) / m
  if(inf_rate < 0)
    inf_rate <- 0
  if(inf_rate > maxrt)
    inf_rate <- maxrt
  return(unname(inf_rate))
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
  unit_inf <- inf_manual(inf_tms = c(0,dtm,tmax_search), inf_rate = c(1,0,0))
  null_inf <- inf_manual(inf_tms = c(0,tmax_search), inf_rate = c(0,0))

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


#' Apply a TCI algorithm to a `pkmod` object
#'
#' Apply a TCI algorithm to a set of targets and a `pkmod` object to calculate infusion
#' rates.
#' @param pkmod `pkmod` object created by `pkmod()`.
#' @param target_vals A vector of numeric values indicating PK or PD targets for TCI algorithm.
#' @param target_tms A vector of numeric values indicating times at which the TCI algorithm should
#' begin targeting each value.
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
#' the minimum value of `target_tms`.
#' @param ignore_pd Logical. Should the PD component of the pkmod object (if present)
#' be ignored. By default, predict.tciinf will assume that 'value' refers to PD
#' targets if a PD model is specified.
#' @param ... Arguments passed to TCI algorithm
#' @examples
#' # 3-compartment model with effect-site
#' my_mod <- pkmod(pars_pk = c(v1 = 8.995, v2 = 17.297, v3 = 120.963, cl = 1.382,
#' q2 = 0.919, q3 = 0.609, ke0 = 1.289))
#' # plasma targeting
#' apply_tci(my_mod, target_vals = c(2,3,4,4), target_tms = c(0,2,3,10), "plasma")
#' # effect-site targeting
#' apply_tci(my_mod, target_vals = c(2,3,4,4), target_tms = c(0,2,3,10), "effect")
#' # incorporate  PD model
#' my_mod_pd <- update(my_mod, pars_pd = c(c50 = 2.8, gamma = 1.47, e0 = 93, emx = 93),
#' pdfn = emax, pdinv = emax_inv)
#' apply_tci(my_mod_pd, target_vals = c(70,60,50,50), target_tms = c(0,2,3,10), "effect")
#' @rdname apply_tci
#' @export
apply_tci <- function(pkmod, target_vals, target_tms, type = c("plasma","effect"), dtm = NULL, custom_alg = NULL, inittm = 0, ignore_pd = FALSE, ...){

  if(length(target_vals) != length(target_tms))
    stop("'target_vals' and 'target_tms' must have the same length")
  # if(!inherits(pkmod,"pkmod")) stop("pkmod must have class 'pkmod'")

  type <- match.arg(type)

  if(is.null(custom_alg)){
    tci_alg <- list(tci_plasma, tci_effect)[[which(type == c("plasma","effect"))]]
  } else{
    if(!all(c("Ct","pkmod","dtm") %in% names(formals(custom_alg))))
      stop('pkmod must contain arguments ("Ct","pkmod","dtm").')
    tci_alg <- custom_alg
  }

  if(is.null(dtm)) dtm <- eval(formals(tci_alg)$dtm)

  tms <- target_tms
  if(any(inittm > tms)) stop("inittm cannot be greater than any target times")

  if(!is.null(pkmod$pdfn) & !ignore_pd){
    Ct <- pkmod$pdinv(target_vals, pars = pkmod$pars_pd)
  } else{
    Ct <- target_vals
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
  dose <- inf_manual(inf_tms = updatetms+inittm, inf_rate = inf, duration = dtm)
  out <- cbind(dose, sf(updatetms), startcon, endcon)
  colnames(out) <- c("begin","end","inf_rate","Ct",paste0("c",1:pkmod$ncmpt, "_start"),
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


#' Target-controlled infusion
#'
#' Apply a TCI algorithm to a set of targets and a `pkmod` or `poppkmod` object to calculate infusion
#' rates.
#' @param pkmod `pkmod` object created by `pkmod()` or a `poppkmod` object created by `poppkmod()`.
#' @param target_vals A vector of numeric values indicating PK or PD targets for TCI algorithm.
#' @param target_tms A vector of numeric values indicating times at which the TCI algorithm should
#' begin targeting each value.
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
#' the minimum value of `target_tms`.
#' @param ignore_pd Logical. Should the PD component of the pkmod object (if present)
#' be ignored. By default, predict.tciinf will assume that 'value' refers to PD
#' targets if a PD model is specified.
#' @param pop_fn Function applied to the distribution of predicted values. E.g.,
#' `median` will calculate doses such that the median value in the population will
#' obtain the target value. Only applicable to `poppkmod` objects.
#' @param ... Arguments passed to TCI algorithm
#' @examples
#' # 3-compartment model with effect-site
#' my_mod <- pkmod(pars_pk = c(v1 = 8.995, v2 = 17.297, v3 = 120.963, cl = 1.382,
#' q2 = 0.919, q3 = 0.609, ke0 = 1.289))
#' # plasma targeting
#' inf_tci(my_mod, target_vals = c(2,3,4,4), target_tms = c(0,2,3,10), "plasma")
#' # effect-site targeting
#' inf_tci(my_mod, target_vals = c(2,3,4,4), target_tms = c(0,2,3,10), "effect")
#' # poppkmod object
#' data <- data.frame(ID = 1:5, AGE = seq(20,60,by=10), TBW = seq(60,80,by=5),
#' HGT = seq(150,190,by=10), MALE = c(TRUE,TRUE,FALSE,FALSE,FALSE))
#' elvd_mod <- poppkmod(data, drug = "ppf", model = "eleveld")
#' inf_tci(elvd_mod, target_vals = c(2,3,4,4), target_tms = c(0,2,3,10), "effect")
#' @export
inf_tci <- function(pkmod, target_vals, target_tms, type = c("plasma","effect"), dtm = NULL, custom_alg = NULL, inittm = 0, ignore_pd = FALSE, pop_fn = NULL,...){

  if(!(inherits(pkmod,"pkmod")|inherits(pkmod, "poppkmod")))
    stop("pkmod must have class 'pkmod' or 'poppkmod'")
  if(inherits(pkmod,"pkmod")){
    apply_tci(pkmod, target_vals, target_tms, type = match.arg(type), dtm = dtm, custom_alg = custom_alg, inittm = inittm, ignore_pd = ignore_pd, ...)
  } else{
    infs <- lapply(pkmod$pkmods, apply_tci, target_vals = target_vals, target_tms=target_tms, type = match.arg(type), dtm = dtm, custom_alg = custom_alg, inittm = inittm, ignore_pd = ignore_pd, ...)
    cbind(id=rep(pkmod$ids, each = nrow(infs[[1]])),do.call("rbind", infs))
  }
}


#' @name inf_tci_pop
#' @title Calculate infusion rates to reach target values within a population
#' @description This function will calculate infusion rates required to reach targets at a specified
#' position (e.g., median, 25th percentile) of the response distribution in a population.
#' The user supplies a `poppkmod` object, a set of values to be obtained, a set of
#' times at which to reach them, and a function to apply to the distribution of responses.
#' @param mod Object with class `poppkmod`, created by `poppkmod()`.
#' @param target_vals A vector of numeric values indicating PK or PD targets for TCI algorithm.
#' @param target_tms A vector of numeric values indicating times at which the TCI algorithm should
#' begin targeting each value.
#' @param pop_fn Function to apply to distribution of response values. Defaults to setting
#' median value equal to targets.
#' @param dtm TCI update frequency. Defaults to 1/6, corresponding to 10-second
#' intervals if model parameters are in terms of minutes.
#' @param inf_duration Optional parameter to describe duration of infusions. Can be
#' less than or equal to `dtm`. Defaults to `dtm` if unspecified.
#' @param cmpt Compartment used for targeting. Defaults to PD response if applicable and
#' effect-site compartment if not.
#' @param inittm Initial time to start TCI algorithm. Cannot be greater than
#' the minimum value of `target_tms`.
#' @param maxrt Maximum infusion rate.
#' @examples
#' nid = 100
#' data <- data.frame(ID = 1:nid, AGE = sample(10:70, nid, replace = TRUE),
#' TBW = sample(40:90, nid, TRUE), HGT = sample(130:210, nid, TRUE),
#' MALE = sample(c(TRUE,FALSE),nid,TRUE))
#' elvd_mod <- poppkmod(data, drug = "ppf", model = "eleveld")
#' # calculate infusions to keep 90% of patients below 70
#' inf <- inf_tci_pop(elvd_mod, target_vals = c(70,70), target_tms = c(0,50),
#' pop_fn = function(...) quantile(..., 0.9), dtm = 10)
#' tms <- seq(0,50,0.1)
#' resp <- as.data.frame(predict(elvd_mod, inf, tms))
#' ggplot(resp, aes(x = time, y = pdresp, group = id)) +
#' geom_line(alpha = 0.1) +
#' geom_hline(yintercept = 70)
#' @export
inf_tci_pop <- function(mod, target_vals, target_tms, pop_fn = median, dtm = 1/6,
                        inf_duration = NULL, cmpt = NULL, inittm = 0, maxrt = 1e5){

  if(length(target_vals) != length(target_tms))
    stop("'target_vals' and 'target_tms' must have the same length")
  if(!inherits(mod,"poppkmod")) stop("mod must have class 'poppkmod'")

  if(is.null(cmpt)){
    # warning("'cmpt' is not specified. Using PD response or effect site concentration")
    cmpt <- ifelse(is.null(mod$pkmods[[1]]$pdfn),
                   paste0("c",mod$pkmods[[1]]$ecmpt),
                   "pdresp")
  }

  # set infusion duration if not specified
  if(is.null(inf_duration)) inf_duration <- dtm
  if(inf_duration > dtm) stop("'inf_duration' cannot be greater than 'dtm'")

  # check initial time
  if(any(inittm > target_tms)) stop("inittm cannot be greater than any target times")

  # Create step function to define targets at any point
  tms <- target_tms-inittm
  sf <- stepfun(tms, c(0,target_vals))
  updatetms <- seq(0, max(tms)-dtm, dtm)

  # initialize values
  inf_all <- rep(NA, length(updatetms))
  ncmpt <- mod$pkmods[[1]]$ncmpt

  # store original initial values
  init_orig <- sapply(mod$pkmods, `[[`, "init")

  # function to minimize
  fopt <- function(infrt, popmod, tm_pred, val){
    inf <- inf_manual(inf_tms = 0, inf_rate = infrt, duration = inf_duration)
    pred <- predict(popmod, inf, tm_pred)[,cmpt]
    return(pop_fn(pred)-val)
  }

  # iterate through times
  for(i in 1:length(updatetms)){
    # check bounds
    lv <- fopt(0, mod, dtm, sf(updatetms[i]))
    uv <- fopt(maxrt, mod, dtm, sf(updatetms[i]))
    if(sign(lv) == sign(uv)){
      warning("Some targets could not be reached on the interval provided and minimum/maximum
              values were used. Consider increasing 'maxrt' to avoid this.")
      inf_all[i] <- c(0,maxrt)[which.min(abs(c(lv,uv)))]
    } else{
      # calculate infusion rate
      inf_all[i] <- uniroot(fopt, interval = c(0,maxrt), popmod = mod,
                            tm_pred = dtm, val=sf(updatetms[i]))$root
    }

    # predict values at update time
    predi <- predict(mod, inf = inf_manual(0,inf_all[i],inf_duration), tms = dtm)

    # update initial values
    mod$pkmods <- lapply(1:length(mod$pkmods), function(ii){
      update(mod$pkmods[[ii]], init = predi[ii,paste0("c",1:ncmpt)])
    })
  }

  # return infusion rates
  out <- cbind(begin = updatetms, end = updatetms + inf_duration, inf_rate = inf_all)
  return(out)
}


#'  infusion schedule
#'
#' Returns a data frame describing a set of infusions to be administered. Output
#' can be used as argument "inf" in predict_pkmod.
#'
#' @param inf_tms Vector of time to begin infusions. If duration is NULL, times
#' expected to include both infusion start and infusion end times.
#' @param inf_rate Vector of infusion rates. Must lave length equal to either
#' length(inf_tms) or 1. If length(inf_rate)=1, the same infusion rate will be
#' administered at each time. Either inf_rate or target must be specified, but not both.
#' @param duration Optional duration of infusions.
#' @return Matrix of infusion rates, start and end times.
#' @importFrom utils tail
#' @examples
#' # specify start and end times, as well as infusion rates (including 0).
#' inf_manual(inf_tms = c(0,0.5,4,4.5,10), inf_rate = c(100,0,80,0,0))
#' # specify start times, single infusion rate and single duration
#' inf_manual(inf_tms = c(0,4), inf_rate = 100, duration = 0.5)
#' # multiple infusion rates, single duration
#' inf_manual(inf_tms = c(0,4), inf_rate = c(100,80), duration = 0.5)
#' # multiple sequential infusion rates
#' inf_manual(inf_tms = seq(0,1,1/6), inf_rate = 100, duration = 1/6)
#' # single infusion rate, multiple durations
#' inf_manual(inf_tms = c(0,4), inf_rate = 100, duration = c(0.5,1))
#' # multiple infusion rates, multiple durations
#' inf_manual(inf_tms = c(0,4), inf_rate = c(100,80), duration = c(0.5,1))
#' @export
inf_manual <- function(inf_tms, inf_rate, duration = NULL)
{

  if(!is.null(duration) & !(length(duration) %in% c(1, length(inf_tms))))
    stop("duration length must be equal to 1 or length of time vector")
  if(!(length(inf_rate) %in% c(1, length(inf_tms))))
    stop("inf_rate length must be equal to 1 or length of time vector")
  if(!(length(inf_rate) %in% c(1, length(inf_tms))) & tail(inf_rate,1)!=0 & !is.null(duration))
    stop("Final value of inf_rate must be zero or duration must be non-null")
  if(!is.null(duration)){
    endtms <- inf_tms + duration
    alltms <- union(round(inf_tms,5), round(endtms,5))
    ord <- order(alltms)
    if(length(inf_rate) == 1) inf_rate <- c(rep(inf_rate,length(inf_tms)), rep(0,length(endtms)))
    out <- as.matrix(cbind(
      begin = alltms[ord][-length(alltms)],
      end = alltms[ord][-1],
      inf_rate = inf_rate[ord][-length(alltms)]))
    out[is.na(out[,"inf_rate"]),"inf_rate"] <- 0
  } else{
    out <- as.matrix(cbind(begin = inf_tms[-length(inf_tms)],
                           end = inf_tms[-1],
                           inf_rate = inf_rate[-length(inf_rate)]))
  }
  return(out)
}
