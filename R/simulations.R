# --------------------------------------------------------------------------------------------------------------------------------
# - Simulation functions ---------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------

#' Function to simulate data from a specified PK or PK-PD model with a specified infusion schedule.
#' @param inf An infusion rate object outputted from either the 'create_intvl' function or the 'iterate_tci_grid' function
#' @param pkmod PK model
#' @param pars_pk0 "True" parameter estimates used to simulate data observations.
#' @param sigma_add Additive residual error standard deviation.
#' @param sigma_mult Multiplicative residual error standard deviation.
#' @param log_err Logical. Should the error be log-normally distributed?
#' @param init Initial concentrations.
#' @param tms Observation times. Defaults to beginning of each infusion if unspecified.
#' @param pdmod PD model if applicable.
#' @param pars_pd0 PD model parameters if applicable.
#' @param ecmpt Effect-site compartment number. Defaults to last compartment.
#' @param delay Delay between generation and observation of measurements.
#' @param max_pdval Maximum PD value.
#' @param min_pdval Minimum PD value
#' @import stats
#'
#' @export
gen_data <- function(inf, pkmod, pars_pk0,
                     sigma_add = 0, sigma_mult = 0,
                     log_err = FALSE,
                     init = NULL, tms = NULL,
                     pdmod = NULL, pars_pd0 = NULL,
                     ecmpt = NULL, delay = 0,
                     max_pdval = 100, min_pdval = 0){

  if(is.null(tms)){
    tms <- inf[,"begin"]
  }

  if(is.null(init))
    init <- eval(formals(pkmod)$init)

  con0 <- predict(object = pkmod, inf = inf, tms = tms, pars = pars_pk0, init = init)

  # additive and multiplicative errors
  eadd  <- rnorm(nrow(con0),0,sigma_add)
  emult <- rnorm(nrow(con0),0,sigma_mult)
  if(is.null(pdmod)){

    if(log_err)
      con0 <- cbind(con0, cobs = exp(log(con0[,"c1"])*(1+emult) + eadd))
    else con0 <- cbind(con0, cobs = con0[,"c1"]*(1+emult) + eadd)

  } else{
    # possible time delay
    timeobs <- con0[,"time"] + delay

    if(is.null(ecmpt))
      ecmpt <- length(eval(formals(pkmod)$init))

    # pd observations
    pd0 <- pdmod(con0[,paste0("c",ecmpt)], pars_pd0)

    if(log_err)
      pdobs <- exp(log(pd0)*(1+emult) + eadd)
    else pdobs <- pd0*(1+emult) + eadd

    con0 <- cbind(con0,
                  timeobs = timeobs,
                  pd0 = pd0,
                  pdobs = pdobs)

    # replace with max/min values if necessary
    con0[con0[,"pdobs"] > max_pdval,"pdobs"] <- max_pdval
    con0[con0[,"pdobs"] < min_pdval,"pdobs"] <- min_pdval
  }

  out <- list(sim = con0,
              inf = inf,
              init = init,
              pkmod = pkmod,
              pdmod = pdmod,
              pars_pk0 = pars_pk0,
              pars_pd0 = pars_pd0,
              sigma_add = sigma_add,
              sigma_mult = sigma_mult,
              ecmpt = ecmpt,
              delay = delay)
  class(out) <- c(class(out), "datasim")
  return(out)
}





#' Combine simulation outputs
#'
#' Function to merge objects with class datasim from different infusion schedules
#' infusion schedules can be passed directly in or as a list.
#' @param ... Set of datasim objects created from `gen_data` function.
#'
#' @export
combine_sim <- function(...){
  simlist <- list(...)

  # allow for the possibility of passing in null simulations
  simlist[sapply(simlist,is.null)] <- NULL
  if(length(simlist) == 1)
    return(simlist[[1]])

  out <- vector("list", length(simlist[[1]]))
  names(out) <- names(simlist[[1]])

  out$sim <- do.call("rbind", lapply(simlist, `[[`, "sim"))
  out$inf <- do.call("rbind", lapply(simlist, `[[`, "inf"))

  lnames <- names(out[sapply(out, is.null)])
  for(i in 1:length(lnames)){
    if(length(unique(lapply(simlist, `[[`, lnames[i]))) == 1){
      out[[lnames[i]]] <- simlist[[1]][[lnames[i]]]
    } else{
      lapply(simlist, `[[`, lnames[i])
    }
  }

  out$init <- simlist[[1]]$init
  class(out) <- c(class(out), "datasim")
  return(out)
}

#' Apply a population PK model to a data frame
#'
#' Function to apply saved population PK or PK-PD models to a
#' data frame of patient values.
#'
#' @param patient_df Dataframe with patient covariate values. Must have
#' names used by model "mod"
#' @param mod Population PK model to apply to rows of patient_df
#' @param ... Arguments passed on to population PK model.
#' @export
apply_poppk <- function(patient_df, mod = c("marsh","schnider","eleveld"), ...){
  switch(match.arg(mod),
         marsh = marsh_poppk(patient_df, ...),
         schnider = schnider_poppk(patient_df, ...),
         eleveld = eleveld_poppk(patient_df, ...))
}



#' Calculate logged prior value
#'
#' Function to return the prior probability for a set of parameters
#' assuming a log-normal distribution.
#'
#' @param lpr log parameter values to evaluate
#' @param mu mean for model parameters and mean residual error
#' @param sig variance covariance matrix for model parameters
#' @importFrom mvtnorm dmvnorm
#'
#' @export
log_prior <- function(lpr, mu, sig){
  mvtnorm::dmvnorm(lpr, mu, sig, log = TRUE)
}



#' Evaluate log-likelihood
#'
#' Function to evaluate the log likelihood given a set of logged parameter values and a set of observed BIS values.
#' It is assumed that the full set of parameters are given by indices (pk_ix, pd_ix), of which a subset may be
#' fixed (i.e. not updated, but still used to evaluate PK-PD functions).
#'
#' @param lpr Set of logged PK-PD-error parameter values to be updated. The final value of lpr is
#' assumed to be the residual error term.
#' @param dat data frame with columns c("time","bis") corresponding to observed time and bis values
#' @param pk_ix indices of (pars_pk,pars_pd) corresponding to PK function values
#' @param pd_ix indices of (pars_pk,pars_pd) corresponding to PD function values
#' @param fixed_ix indices of (pars_pk,pars_pd) corresponding to PD function values
#' @param fixed_lpr values used by PD function that are not updated.
#' @importFrom truncnorm dtruncnorm
#'
#' @export
log_likelihood <- function(lpr, dat, pk_ix, pd_ix, fixed_ix = NULL, fixed_lpr = NULL){

  # reconstruct vector of PK-PD parameters
  err <- exp(lpr[length(lpr)])
  epars <- rep(NA, length(c(pk_ix,pd_ix)))
  epars[fixed_ix] <- exp(fixed_lpr)
  epars[is.na(epars)] <- exp(lpr[-length(lpr)])
  names(epars)[fixed_ix] <- names(fixed_lpr)
  names(epars)[-fixed_ix] <- names(lpr)[-length(lpr)]

  # restrict times to allow for time lag
  if("timeobs" %in% colnames(dat$sim)){
    tm_all <- dat$sim[,"time"]
    tm_obs <- dat$sim[,"timeobs"]
    tms <- tm_all[tm_obs <= max(tm_all)]
    val_obs <- dat$sim[tm_obs <= max(tm_all),"pdobs"]
  } else{
    tms <- dat$sim[,"time"]
    if("pdobs" %in% colnames(dat$sim))
      val_obs <- dat$sim[,"pdobs"]
    else val_obs <- dat$sim[,"cobs"]
  }

  # predicted initial concentrations from data
  ini <- dat$inf[1,grep("c[0-9]_start",colnames(dat$inf))]

  # predict concentrations at lpr
  cp <- predict(dat$pkmod,
                inf = dat$inf,
                tms = tms,
                pars = epars[pk_ix],
                init = ini)

  # if pd model isn't given, evaluate likelihood at pk observations
  if(is.null(dat$pdmod)){
    return(sum(log(truncnorm::dtruncnorm(x = val_obs, mean = cp[,"c1"], sd = err, a = 0))))
  } else{
    # apply function to lpr if applicable - e.g. have parameters change based on predicted concentration
    econ <- cp[,paste0("c",dat$ecmpt)]
    bisp <- dat$pdmod(ce = econ, pars = epars[pd_ix])
    return(sum(log(truncnorm::dtruncnorm(x = val_obs, mean = bisp, sd = err, a = 0, b = 100))))
  }
}


# log_likelihood <- function(lpr, dat, pk_ix, pd_ix, err_ix){
#
#   # restrict times to allow for time lag
#   if("timeobs" %in% colnames(dat$sim)){
#     tm_all <- dat$sim[,"time"]
#     tm_obs <- dat$sim[,"timeobs"]
#     tms <- tm_all[tm_obs <= max(tm_all)]
#     val_obs <- dat$sim[tm_obs <= max(tm_all),"pdobs"]
#   } else{
#     tms <- dat$sim[,"time"]
#     if("pdobs" %in% colnames(dat$sim))
#       val_obs <- dat$sim[,"pdobs"]
#     else val_obs <- dat$sim[,"cobs"]
#   }
#
#   # predicted initial concentrations from data
#   ini <- dat$inf[1,grep("c[0-9]_start",colnames(dat$inf))]
#
#   # predict concentrations at lpr
#   cp <- predict(dat$pkmod,
#                 inf = dat$inf,
#                 tms = tms,
#                 pars = exp(lpr[pk_ix]),
#                 init = ini)
#
#   # if pd model isn't given, evaluate likelihood at pk observations
#   if(is.null(dat$pdmod)){
#     return(sum(log(truncnorm::dtruncnorm(x = val_obs, mean = cp[,"c1"], sd = exp(lpr[err_ix]), a = 0))))
#   } else{
#     # apply function to lpr if applicable - e.g. have parameters change based on predicted concentration
#     econ <- cp[,paste0("c",dat$ecmpt)]
#     bisp <- dat$pdmod(ce = econ, pars = exp(c(lpr[pd_ix],fixed_lpr)))
#     return(sum(log(truncnorm::dtruncnorm(x = val_obs, mean = bisp, sd = exp(lpr[err_ix]), a = 0, b = 100))))
#   }
# }


#' Function to evaluate the negative log posterior given a set of logged parameter values and observed BIS values.
#' @param lpr logged PK-PD-error parameter values
#' @param dat data frame with columns corresponding to observed time and PD response values.
#' @param mu Mean of prior distribution.
#' @param sig Variance-covariance matrix of prior distribution.
#' @param ... Arguments passed on to log-likelihood.
#'
#' @export
log_posterior_neg <- function(lpr, dat, mu, sig, ...) {
    -1*(log_prior(lpr,mu,sig) + log_likelihood(lpr, dat, ...))
}


#' Closed-loop targets
#'
#' Format data frame of closed-loop targets.
#' @param time Times at which target values are set
#' @param target Response target values
#'
#' @export
cl_targets <- function(time, target){
  data.frame(time = time, target = target)
}


#' Closed-loop updates
#'
#' Set parameters for closed-loop updates.
#'
#' @param time Times at which PK or PK-PD parameters should be updated
#' @param full_data Vector of logical values indicating if all data should be used at update
#' or only data since last update. Once first "FALSE" value is observed, the prior variance-
#' covariance matrix is overwritten. Consequently, any "TRUE" updates will only use data since
#' the last "FALSE" update.
#' @param plot_progress Vector of logical values. Should values be plotted at each update?
#'
#' @export
cl_updates <- function(time, full_data = TRUE, plot_progress = FALSE){
  data.frame(time = time, full_data = full_data, plot_progress = plot_progress)
}


#' Bayesian closed-loop control
#'
#' Function to provide Bayesian closed-loop control.
#'
#' @param targets Data frame with columns ("time","target")
#' @param updates Data frame of times at which closed-loop updates should be conducted and
#' optional variable with logical values named 'full_data' indicating if full updates should
#' be used. Defaults to partial.
#' @param prior List with elements "mu" and "sig" specifying the prior mean and covariance
#' matrices for the logged parameter values.
#' @param true_pars Vector of true patient PK-PD parameters.
#' @param pkmod PK model
#' @param pdmod PD model
#' @param pdinv Inverse PD model
#' @param init0 True initial concentrations
#' @param init_p Predicted initial concentrations
#' @param obs_tms Times at which observations are collected. If null, observations will be
#' made at fixed intervals specified by 'dtm'.
#' @param dt_obs Interval between measurements.
#' @param sim_starttm Start time of simulation
#' @param tci_alg TCI algorithm used. Defaults to effect-site targeting.
#' @param print_progress Logical. Should current update times be printed to the console.
#'
#' @importFrom utils head tail
#' @export
bayes_control <- function(targets, updates, prior, true_pars,
                          pkmod = pkmod3cptm, pdmod = emax_eleveld,
                          pdinv = inv_emax_eleveld,
                          init0 = NULL, init_p = NULL, obs_tms = NULL,
                          dt_obs = 1/6, sim_starttm = 0, tci_alg = "effect",
                          print_progress = FALSE){

  # set observation/measurement times
  if(is.null(obs_tms)) obs_tms <- seq(dt_obs, max(targets$time), dt_obs)

  # set true and predicted initial concentrations if not specified
  ncpt <- length(eval(formals(pkmod)$init))
  if(is.null(init0)) init0 <- rep(0,ncpt)
  if(is.null(init_p)) init_p <- rep(0,ncpt)
  init_start <- init0

  if(is.vector(updates)) updates <- data.frame(time = updates)
  if(!("time" %in% names(updates))) stop("dataframe updates must have column named 'time'")
  if(!("full_data" %in% names(updates))) updates$full_data <- FALSE
  if(!("plot_progress" %in% names(updates))) updates$plot_progress <- FALSE

  # add simulation start time to list of update times
  update_tms <- c(sim_starttm, updates$time)
  update_full <- c(NA, updates$full_data)
  plot_progress <- c(NA, updates$plot_progress)

  # combine update times into set of target times
  na.locf <- function(x) {
    v <- !is.na(x)
    c(NA, x[v])[cumsum(v)+1]
  }
  targets_new <- data.frame(time = sort(union(targets$time, update_tms)))
  targets_new <- merge(targets_new, targets, all.x = TRUE)
  targets_new$target <- na.locf(targets_new$target)

  true_pk <- unlist(true_pars$pars_pkpd[true_pars$pk_ix])
  true_pd <- unlist(true_pars$pars_pkpd[-true_pars$pk_ix])

  prior$pars_pkpd <- unlist(prior$pars_pkpd)
  prior0 <- prior
  dat0 <- NULL
  lpr_all <- NULL

  for(i in 2:length(update_tms)){

    if(print_progress) print(paste("Updating at time t =", update_tms[i]))

    prior_pk <- prior$pars_pkpd[prior$pk_ix]
    prior_pd <- prior$pars_pkpd[prior$pd_ix]

    # subset targets and observation times to period being updated
    targets_sub <- targets_new[targets_new$time <= update_tms[i] & targets_new$time >= update_tms[i-1],]
    obs_tms_sub <- obs_tms[obs_tms <= update_tms[i] & obs_tms > update_tms[i-1]]

    # calculate tci infusions at prior parameter estimates for update period
    inf <- tci_pd(pdresp  = targets_sub$target,
                  tms     = targets_sub$time,
                  pdinv   = pdinv,
                  pdmod   = pdmod,
                  pkmod   = pkmod,
                  pars_pk = prior_pk,
                  pars_pd = prior_pd,
                  init = init_p,
                  tci_alg = tci_alg)

    # generate data under true model
    dat <- gen_data(inf = inf,
                    tms = obs_tms_sub,
                    pkmod = pkmod,
                    pdmod = pdmod,
                    pars_pk0 = true_pk, # true pk parameters
                    pars_pd0 = true_pd, # true pd parameters
                    sigma_add = true_pars$err, # random error
                    delay = true_pars$delay, # bis delay in minutes
                    max_pdval = 100,
                    init = init0)

    # merge sampled dataset with prior observations
    dat0 <- combine_sim(dat0,dat)

    # update parameters based on generated data
    # use prior parameter values as starting point - separate fixed parameters
    if(any(!is.null(prior$fixed_ix))){
      lpr <- log(c(prior$pars_pkpd[-prior$fixed_ix], err = prior$err))
      lpr_fixed <- log(prior$pars_pkpd[prior$fixed_ix])
    } else{
      lpr <- log(c(prior$pars_pkpd, err = prior$err))
      lpk_fixed <- NULL
    }

    lpr_all <- rbind(lpr_all, lpr)

    # indicate if full dataset should be used for updates
    if(update_full[i]){
      dat_eval <- dat0

      # use full dataset and original vcov matrix
      post_est <- optim(par = lpr,
                        fn = log_posterior_neg,
                        dat = dat_eval,
                        mu = lpr,
                        sig = prior0$sig,
                        pk_ix = prior$pk_ix,
                        pd_ix = prior$pd_ix,
                        fixed_ix = prior$fixed_ix,
                        fixed_lpr = lpr_fixed,
                        method = "BFGS",
                        hessian = FALSE)

    } else{
      dat_eval <- dat
      post_est <- optim(par = lpr,
                        fn = log_posterior_neg,
                        dat = dat_eval,
                        mu = lpr,
                        sig = prior$sig,
                        pk_ix = prior$pk_ix,
                        pd_ix = prior$pd_ix,
                        fixed_ix = prior$fixed_ix,
                        fixed_lpr = lpr_fixed,
                        method = "BFGS",
                        hessian = TRUE)

     # update vcov matrix
     prior$sig <- solve(post_est$hessian)
    }

    if(plot_progress[i]){
      print(plot(dat0, lpars_update = post_est$par,
                         lpars_fixed = log(prior$pars_pkpd[prior$fixed_ix])))
    }

    # update prior values
    if(any(!is.null(prior$fixed_ix))){
      prior$pars_pkpd[-prior$fixed_ix] <- exp(head(post_est$par,-1))
      prior$err <- exp(tail(post_est$par,1))
    } else{
      prior$pars_pkpd <- exp(head(post_est$par,-1))
      prior$err <- exp(tail(post_est$par,1))
    }

    # update true and predicted initial values
    init0 <- dat0$sim[nrow(dat0$sim),grep("c[0-9]",colnames(dat0$sim))]
    init_p <- as.numeric(predict(pkmod,
                                 inf = dat0$inf,
                                 tms = update_tms[i],
                                 pars = prior$pars_pkpd[prior$pk_ix],
                                 init = init_start)[-1])

  }

  # save final posterior parameter values
  if(any(!is.null(prior$fixed_ix))){
    lpr <- log(c(prior$pars_pkpd[-prior$fixed_ix], err = prior$err))
    lpr_fixed <- log(prior$pars_pkpd[prior$fixed_ix])
  } else{
    lpr <- log(c(prior$pars_pkpd, err = prior$err))
    lpk_fixed <- NULL
  }

  lpr_all <- rbind(lpr_all, lpr)

  out <- list(dat = dat0,
              lpr = lpr_all,
              posterior = prior,
              prior = prior0,
              true_pars = true_pars)
  class(out) <- "bayessim"
  return(out)
}



#' Sigmoid target function
#'
#' @param lpars Logged parameter values
#' @param tms Times to evaluate sigmoid function
#' @param bis0 BIS value with no drug administered
#' @param ... Arguments passed on to 'restrict_sigmoid' function
#'
#' @export
sigmoid_targetfn <- function(lpars, tms, bis0 = 93, ...)
  emax(tms, restrict_sigmoid(t50 = exp(lpars), BIS0 = bis0, ...))


# #' Cubic target function
# #' Cubic function restricted to pass through BIS = 93 at t=0 and BIS = 50 at t=5
# cubic_targetfn <- function(pars, tms, bis0 = 93, tfinal = 10, bisfinal = 50){
#   b1 = pars[1]
#   b2 = pars[2]
#   b3 = 1/tfinal^3 * (bisfinal - bis0 - tfinal*b1 - tfinal^2*b2)
#   return(bis0 + b1*tms + b2*tms^2 + b3*tms^3)
# }



#' Apply target function to a PK-PD model
#'
#' Function to apply any specified target function to a PK-PD model
#' and TCI algorithm.
#'
#' @param lp Logged parameter values
#' @param tm Time values to evaluate
#' @param targetfn Target function
#' @param prior_pk Prior PK point estimates
#' @param prior_pd Prior PD point estimates
#' @param pkmod PK model to evaluate
#' @param pdmod PD model to evaluate
#' @param pdinv Inverse PD model
#' @param ... Additional arguments passed on to tci_pd
#'
#' @export
apply_targetfn <- function(lp, tm, targetfn, prior_pk, prior_pd,
                           pkmod = pkmod3cptm,
                           pdmod = emax_eleveld,
                           pdinv = inv_emax_eleveld, ...){

  return(tci_pd(pdresp  = targetfn(lp, tm, ...),
                tms     = tm,
                pdinv   = pdinv,
                pdmod   = pdmod,
                pkmod   = pkmod,
                pars_pk = prior_pk,
                pars_pd = prior_pd,
                init = c(0,0,0,0)))
}




# # objective function minimizing weighted combination of target over- and under-shoot
# # over- and under-shoot is relative to the target BIS = 50, NOT an indifference region
# # (e.g. 40 to 60).
# Phi_WgtOvershoot <- function(inf,
#                              pars_pk0,
#                              pars_pd0,
#                              alpha = 0.2,
#                              target = 50,
#                              dtm = 1/20,
#                              init = NULL,
#                              pkmod = pkmod3cptm,
#                              pdmod = emax_eleveld,
#                              pdinv = inv_emax_eleveld,
#                              ecmpt = NULL){
#
#   if(is.null(init))
#     init <- eval(formals(pkmod)$init)
#
#   if(is.null(ecmpt))
#     ecmpt <- length(eval(formals(pkmod)$init))
#
#   # inf <- as.data.frame(inf)
#
#   rng <- range(inf[,c("begin","end")])
#   tms <- seq(dtm,rng[2],dtm)
#
#   # predict concentrations and responses
#   con0 <- predict(pkmod = pkmod, inf = inf, tms = tms, pars = pars_pk0, init = init)
#   pd_pred <- pdmod(con0[,paste0("c",ecmpt)], pars_pd0)
#
#   phi1 <- sum((pd_pred[pd_pred>target] - target))*dtm # integral above target
#   phi2 <- sum((target - pd_pred[pd_pred<target]))*dtm # integral below target
#
#   return(alpha*phi1 + (1-alpha)*phi2)
# }
#
#
# Phi_WgtOvershoot_region <- function(inf,
#                              pars_pk0,
#                              pars_pd0,
#                              alpha = 0.2,
#                              ub = 60,
#                              lb = 40,
#                              dtm = 1/20,
#                              init = NULL,
#                              pkmod = pkmod3cptm,
#                              pdmod = emax_eleveld,
#                              pdinv = inv_emax_eleveld,
#                              ecmpt = NULL){
#
#   if(is.null(init))
#     init <- eval(formals(pkmod)$init)
#
#   if(is.null(ecmpt))
#     ecmpt <- length(eval(formals(pkmod)$init))
#
#   rng <- range(inf[,c("begin","end")])
#   tms <- seq(dtm,rng[2],dtm)
#
#   # predict concentrations and responses
#   con0 <- predict(pkmod = pkmod, inf = inf, tms = tms,
#                   pars = pars_pk0, init = init)[,paste0("c",ecmpt)]
#   pd_pred <- pdmod(con0, pars_pd0)
#
#   phi1 <- sum((pd_pred[pd_pred>ub] - ub))*dtm # integral above target
#   phi2 <- sum((lb - pd_pred[pd_pred<lb]))*dtm # integral below target
#
#   return(alpha*phi1 + (1-alpha)*phi2)
# }
#
#
#
# #' Function to evaluate an objective function phi over a population with infusions defined
# #' by a target function.
# #' This should have an option to sample a set of parameters from a covariance matrix if
# #' specific parameters are not provided.
#
# Phi_ED <- function(lp, targetfn, tms, Phi, alpha,
#                    prior_pars_pk, prior_pars_pd,
#                    eb_pars_pk, eb_pars_pd,
#                    targetfn_args = list(), phi_args = list()){
#
#   library(parallel)
#
#   if(nrow(prior_pars_pk) != nrow(prior_pars_pd)) stop("PK and PD parameters must have the same number of rows.")
#   npars <- nrow(prior_pars_pk)
#
#   # evaluate the infusion schedule at each
#   # phi_i <- mclapply(1:npars, function(i){
#   phi_i <- lapply(1:npars, function(i){
#     targetfn_args$lp <- lp
#     targetfn_args$tm <- tms
#     targetfn_args$targetfn <- targetfn
#     targetfn_args$prior_pk <- unlist(prior_pars_pk[i,])
#     targetfn_args$prior_pd <- unlist(prior_pars_pd[i,])
#
#     # infusion schdule
#     infi <- do.call(apply_targetfn, targetfn_args)
#
#     # calculate objective function
#     Phi(inf = infi,
#         alpha = alpha,
#         pars_pk0 = unlist(eb_pars_pk[i,]),
#         pars_pd0 = unlist(eb_pars_pd[i,]),
#         dtm = 1/200)
#   })
#
#   return(mean(unlist(phi_i)))
# }


# --------------------------------------------------------------------------------------------------------------------------------
# - Old functions ----------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
#
# #' Return an infusion schedule defined by an Emax sigmoid curve.
# #' @param beta Parameters (c50, gamma) for Emax target function
# #' @param theta Vector of PK-PD parameters (theta_PK, theta_PD)
# #' @param ini Initial concentration values
# #' @param t0 Starting time
# #' @param E0 Effect at ce = 0. Used as starting point for sigmoid curve.
# #' @param gamma Slope at c50
# #' @param tmx Final time to evaluate infusion schedule to.
# #' @param bist Target BIS, used to set the maximum effect value of the target sigmoid curve.
# #' @param delta Time interval between TCI updates. Defaults to 1/6 minutes = 10 seconds.
# sig_inf <- function(beta, theta, ini, t0, E0, gamma, tmx = 10, bist = 50, delta = 1/6, ...){
#   if(t0+delta <= tmx) tm_eval <- seq(t0+delta,tmx,delta)
#   else tm_eval <- tmx
#   bis_targets <- E0-(E0-bist)*(tm_eval^beta[2] / (tm_eval^beta[2] + beta[1]^beta[2]))
#   ce_targets <- Hinv(pars = theta[10], bis = bis_targets, E0 = E0, gamma = gamma)
#   kR_future <- TCI_EffectSite(Cet = ce_targets, pars = theta[1:9], init = ini, delta = delta, ...)
#   if(t0 != 0){
#     kR_vec <- unlist(kR_future)
#     kR_vec[names(kR_vec) == "begin"] <- kR_vec[names(kR_vec) == "begin"] + t0
#     kR_vec[names(kR_vec) == "end"] <- kR_vec[names(kR_vec) == "end"] + t0
#     kR_future <- relist(kR_vec, kR_future)
#   }
#   return(kR_future)
# }
#
# #' Objective function for a single infusion - numerically estimates the weighted integral between the time-BIS curve and the target BIS value.
# #' @param kR Infusion schedule to be evaluated.
# #' @param theta Vector of PK-PD parameters (theta_PK, theta_PD)
# #' @param E0 Effect at ce = 0. Used as starting point for sigmoid curve.
# #' @param gamma Hill parameter (slope at c50) for PD model
# #' @param alpha Weight associated with the integral above the target BIS value and below the time-BIS curve.
# #' @param bist Target BIS, used to set the maximum effect value of the target sigmoid curve.
# #' @param tfinal Final time to evaluate infusion schedule to.
# #' @param dtm Resolution of integral
# #' @export
# Phi <- function(kR, theta, E0, gamma, t0, alpha = 0.01, bist = 50, tfinal = 10, dtm = 1e-2, log = F){
#   tms <- seq(0,tfinal,length.out = 1/dtm)
#   ce <- pk_solution_3cpt_metab(pars = theta[1:9], ivt = kR, init = c(0,0,0,0))(tms)[4,]
#   ce[ce < 0] <- 0
#   bisp <- Emax1(theta[10],ce = ce, gamma = gamma, E0 = E0)
#   phi1 <- sum((bisp[bisp>bist] - bist))*dtm*tfinal # integral above target
#   phi2 <- sum((bist - bisp[bisp<bist]))*dtm*tfinal # integral below target
#   if(log) return(log(alpha*phi1 + (1-alpha)*phi2))
#   else return(alpha*phi1 + (1-alpha)*phi2)
# }
#
#
# sim_BIS <- function(pars, ivt, init, gamma, E0, BIS_sampling_freq = 1){
#   pars_pk <- pars[1:9]
#   pars_pd <- pars[10]
#   sampling_window <- c(ivt[[1]]$begin, tail(ivt,1)[[1]]$end)
#   delta_tms_min <- BIS_sampling_freq / 60
#   tm_seq <- seq(sampling_window[1] + delta_tms_min, sampling_window[2], by = delta_tms_min)
#   sol <- pk_solution_3cpt_metab(pars_pk, ivt = ivt, init = init)
#   con_est <- sol(tm_seq)
#   bis_vals <- Emax1(pars = pars_pd, ce = con_est[4,], gamma = gamma, E0 = E0)
#   out <- list(bis = bis_vals, con = con_est, bis_tms = tm_seq)
#   return(out)
# }
#
# quad_pars <- function(lpars, sig, kR, gamma, E0, nquad = 4, ev_dim = 3, prune = 0){
#   b <- function(lpars, ...) sim_BIS(pars = exp(lpars), ...)$bis
#   grad <- fdGrad(pars = lpars[1:10], fun = b, ivt = kR, init = c(0,0,0,0), gamma = gamma, E0 = E0, BIS_sampling_freq = 1)
#   if(any(is.nan(rowSums(grad)))) grad <- grad[-is.nan(rowSums(grad)),]
#   ev <- svd(grad)
#   sig_ev <- t(ev$v[,c(1:ev_dim)]) %*% sig[1:10,1:10] %*% ev$v[,c(1:ev_dim)]
#   quad <- mgauss.hermite(n = nquad, mu = rep(0,ev_dim), sigma = sig_ev, prune = prune)
#   return(list(points = t(t(quad$points %*% t(ev$v[,1:ev_dim])) + lpars[1:10]),
#               weights = quad$weights,
#               var_exp = cumsum(ev$d^2/ sum(ev$d^2))[ev_dim]))
# }
#
# # quasi-monte-carlo
# quad_pars_qmc <- function(lpars, sig, npoints = 200, par_dim = 10){
#   xy = randtoolbox::sobol(npoints,dim=par_dim)
#   z = apply(xy,2,qnorm)
#   ev <- eigen(sig[1:par_dim,1:par_dim])
#   A <- t(t(ev$vectors)*sqrt(ev$values))
#   y <- t(t(A)%*%t(z) + lpars[1:par_dim])
#   return(y)
# }
#
# sample_data <- function(ivt_d, pars_pk_tci, pars_pk0, pars_pd0, sigma_bis0, gamma1, gamma2, E0, start_time, init_pred, init_true, delay, nmin_sample = 1/6, sample_freq = 1/60){
#   # times to be evaluated (starting at 0)
#   tms0 <- seq(sample_freq, nmin_sample, sample_freq)
#
#   # true concentration is used to generate observed BIS
#   sol0 <- pk_solution_3cpt_metab(pars = pars_pk0, ivt = ivt_d, init = init_true)
#   con0 <- sol0(tms0)
#
#   # simulate BIS from true concentrations
#   gamma0 <- ifelse(con0[4,] <= pars_pd0, gamma1, gamma2)
#   bis_t <- Emax1(pars = pars_pd0, ce = con0[4,], gamma = gamma0, E0 = E0)
#   bis_obs <- truncnorm::rtruncnorm(n = length(con0[4,]), mean = bis_t, sd = sigma_bis0, a = 0, b = 100)
#
#   # real times being sampled
#   tms_sample <- seq(start_time + sample_freq, start_time + nmin_sample, sample_freq)
#   tmsobs <- tms_sample + delay/60
#
#   ivt_d[[1]]$begin <- ivt_d[[1]]$begin + start_time
#   ivt_d[[1]]$end <- ivt_d[[1]]$end + start_time
#
#   dat0 <- data.frame(cbind(time = tms_sample, timedelay = tmsobs, bis = bis_obs, bis_t = bis_t,
#                            c1_t = con0[1,], c2_t = con0[2,], c3_t = con0[3,], ce_t = con0[4,]))
#   return(dat0)
# }
#
#
# update_pars <- function(lp, dat, ivt, lpr, gamma, E0){
#   post_est <- nlm(f = log_posterior_neg, p = lp, ivt=ivt, dat=dat,  lhyper = lpr, gamma = gamma, E0 = E0,
#                   hessian = T,
#                   steptol=1e-6, gradtol=1e-6, stepmax = 5,
#                   iterlim = 2000)
#   post_est_pars <- post_est$estimate
#   post_hes <- post_est$hessian
#   lpost <- list(mu = post_est_pars, sig = solve(post_hes))
#   return(lpost)
# }
#
# sigmoid_update_objfn <- function(lgamma, E0, bis_target = 50, nmin = 10, nmin_target = 5, eps = 1, ...){
#   # bis_t50 = nmin*(eps/(bis0 - bis_target + eps))^(1/gamma) # express c50 parameter in terms of gamma
#   lbis_t50 = log(nmin_target) + 1/exp(lgamma) * (log(eps) - log(E0 - bis_target + eps))
#   sig_ED(lbeta = c(lbis_t50, lgamma), E0 = E0, bis_target = bis_target, nmin = nmin, epsilon = 0, ...)
# }
#
# sigmoid_update_objfn2 <- function(lbeta, E0, bis_target = 50, nmin = 10, nmin_target = 5, eps = 1, ...){
#   sig_ED(lbeta = lbeta, E0 = E0, bis_target = bis_target, nmin = nmin, epsilon = 0, ...)
# }
#
#
# sig_ED <- function(lbeta, ltheta_hat, Theta_samples, weights, gamma, E0, sig, ivt0, t0, t1 = NULL, init_p, p = NULL, alpha = 0.01, epsilon = 0.1, bis_target = 50, nmin = 10, dtm = 1e-2){
#   if(is.null(t1)) t1 <- nmin
#   kR_all <- c(ivt0, sig_inf(beta =exp(lbeta), theta = exp(ltheta_hat), ini = init_p, t0 = t0, tmx = t1, E0 = E0, gamma = gamma))
#
#   if(is.null(p)){
#     # evaluate quadrature points at infusion schedule
#     obj_fn_val <- c(t(vapply(1:nrow(Theta_samples), function(l){
#       Phi(kR = kR_all, theta = Theta_samples[l,], alpha = alpha, bist = bis_target, tfinal = nmin,
#           t0 = t0,
#           dtm = dtm, gamma = gamma, E0 = E0)
#     }, FUN.VALUE = numeric(1))) %*% weights)
#   } else{
#     # evaluate quadrature points at infusion schedule
#     obj_fn_val <- quantile(vapply(1:nrow(Theta_samples), function(l){
#       Phi(kR = kR_all, theta = Theta_samples[l,], alpha = alpha, bist = bis_target, tfinal = nmin,
#           t0 = t0,
#           dtm = dtm, gamma = gamma, E0 = E0)
#     }, FUN.VALUE = numeric(1)), probs = p)
#   }
#   return(obj_fn_val)
# }
#
#
#
# sigmoid_induct_eleveld <- function(pars0, lhyper, target_update_tms, alpha = 0.01, p = NULL, epsilon = 0.1, delta = 1/6,
#                                    nmin = 10, nmin_prop = 1/2, beta_init = c(2,3), fixed_bis = F,
#                                    bis_target = 50, plot_path = T, delta_bis = 5, update_pars_logical = T,
#                                    qmc = F, laplace_appx = T, seed = NULL, laplace_update_tm = 4, eps = 1){
#   library(mvtnorm, quietly = T)
#   library(numDeriv, quietly = T)
#   library(nloptr, quietly = T)
#   if(!is.null(seed)) set.seed(seed)
#   pars_pk0 <- pars0$pk
#   pars_pd0 <- pars0$pd
#   sigma_bis0 <- pars0$sigma
#
#   # assumed to not vary by individual
#   lag0 = pars0$lag
#   gamma1 = pars0$gamma
#   gamma2 = pars0$gamma2
#   E0 = pars0$E0
#
#   update_tms <- seq(0, nmin-delta, delta)
#   target_update_ix <- which(update_tms %in% target_update_tms)
#   next_update_tm <- c(target_update_tms[-1],10)
#   dat_obs <- as.data.frame(matrix(NA, ncol = 12, nrow = 0))
#   lpars_prior <- lhyper
#   # lpars_post_list <- matrix(NA, nrow = length(update_tms), ncol = length(lpars_prior$mu))
#
#   for(r in 1:length(update_tms)){
#     print(r)
#     pars_pk_current <- exp(lpars_prior$mu[1:9])
#     pars_pd_current <- exp(lpars_prior$mu[10])
#
#     # Set initial values
#     if(nrow(dat_obs) == 0) {
#       init_0 <- init_p <- c(0,0,0,0) # predicted and true initial concentrations
#       bis_p <- E0
#       t0 = 0
#       beta_est <- beta_init
#       gamma_eval = gamma1
#       ivt_eval <- NULL
#     } else{
#       # true concentrations
#       init_0 <- as.numeric(dat_obs[nrow(dat_obs),c("c1_t","c2_t","c3_t","ce_t")])
#       # current time of evaluation
#       t0 <- ivt_eval[[length(ivt_eval)]]$end
#       # predicted concentration given infusion schedule administered
#       init_p <- c(t(pk_solution_3cpt_metab(pars = pars_pk_current, ivt = ivt_eval, init = c(0,0,0,0))(t0)))
#       # gamma value isn't updated, but switches based on whether the concentration is greater than C50
#       gamma_eval <- unname(ifelse(init_p[4] <= pars_pd_current, gamma1, gamma2))
#       # currently predicted BIS
#       bis_p <- Emax1(pars = pars_pd_current, ce = init_p[4], gamma = gamma_eval, E0 = E0)
#     }
#
#     if(fixed_bis){ # target a fixed value of BIS directly
#       ce_target <- Hinv(pars = pars_pd_current, bis = bis_target, E0 = E0, gamma = gamma_eval)
#       kR_new <- TCI_basic(Ce = ce_target, pars = pars_pk_current, init = init_p)
#       ivt_new <- list(list(begin = t0, end = t0 + delta, k_R = kR_new))
#     } else{
#
#       if(is.null(target_update_tms)) t1 = t0 + delta
#
#       # update beta coefficients
#       if(r %in% target_update_ix){
#         print("Updating target sigmoid function")
#
#         # evaluate infusion schedule up until the next update time
#         t1 <- next_update_tm[which(target_update_ix == r)]
#
#         kR_quad <- c(ivt_eval, sig_inf(beta = beta_init, theta = c(pars_pk_current, pars_pd_current), ini = init_p, t0 = t0, gamma = gamma_eval, E0 = E0, tmx = t1))
#         # }
#
#         if(qmc){
#           samples <- quad_pars_qmc(lpars = lpars_prior$mu, sig = lpars_prior$sig, npoints = 100)
#           quad <- list(points = samples, weights = rep(1,nrow(samples)))
#         } else{
#           quad <- quad_pars(lpars_prior$mu, sig = lpars_prior$sig, kR = kR_quad, gamma = gamma_eval, E0 = E0)
#         }
#
#         if(!is.null(eps)){
#           # optimize version - one parameter constrained to equal bis target + eps at final point
#           # nmin_target_prop = t0 + nmin_prop*(nmin-t0)
#           nmin_target_prop = nmin_prop*nmin
#           opt_res <- optimize(f = sigmoid_update_objfn, interval = c(0,3),
#                               ltheta_hat = lpars_prior$mu,
#                               Theta_samples = exp(quad$points), weights = quad$weights,
#                               sig = lpars_prior$sig,
#                               t0 = t0, ivt0 = ivt_eval, init_p = init_p, t1 = t1,
#                               E0 = E0, gamma = gamma_eval,
#                               alpha = alpha, bis_target = bis_target, p = p,
#                               nmin = nmin,
#                               nmin_target = nmin_target_prop,
#                               eps = eps)
#           lgamma_opt = opt_res$minimum
#           lbis_t50_opt = log(nmin_target_prop) + 1/exp(lgamma_opt) * (log(1) - log(E0 - bis_target + 1))
#           beta_est <- exp(c(lbis_t50_opt, lgamma_opt))
#         } else{
#           optimx_res <- optimx::optimx(par = log(beta_est), fn = sig_ED, ltheta_hat = lpars_prior$mu,
#                                        Theta_samples = exp(quad$points), weights = quad$weights,
#                                        sig = lpars_prior$sig, t0 = t0, ivt0 = ivt_eval, init_p = init_p, E0 = E0, gamma = gamma_eval,
#                                        alpha = alpha, bis_target = bis_target, nmin = nmin,
#                                        method = "nlm", control = list(rel.tol = 1e-4))
#           beta_est <- c(unname(exp(coef(optimx_res))))
#         }
#
#         print(beta_est)
#       }
#
#       ivt_new <- sig_inf(beta = beta_est, theta = exp(lpars_prior$mu), ini = init_p, t0 = t0, E0 = E0, gamma = gamma_eval,
#                          tmx = t1,
#                          bist = bis_target)
#     }
#
#     ivt_eval <- c(ivt_eval, ivt_new[1])
#     dat_i <- sample_data(ivt_d = list(list(begin = 0, end = delta, k_R = ivt_eval[[length(ivt_eval)]]$k_R)),
#                          pars_pk_tci = pars_pk_current, pars_pk0 = pars_pk0, pars_pd0 = pars_pd0, sigma_bis0 = sigma_bis0,
#                          gamma1 = gamma1, gamma2 = gamma2, E0 = E0, delay = lag0, start_time = update_tms[r],
#                          init_pred = init_p, init_true = init_0, nmin_sample = delta)
#     dat_obs <- rbind(dat_obs, dat_i)
#
#     if(update_pars_logical & update_tms[min(r+1, length(update_tms))] >= lag0/60){
#       if(laplace_appx & update_tms[r] > laplace_update_tm){
#         dat_eval = dat_obs[dat_obs$timedelay <=(update_tms[r] + delta) & dat_obs$timedelay > update_tms[r],]
#         lpr_eval = lpars_prior
#       } else{
#         dat_eval = dat_obs[dat_obs$timedelay <= (update_tms[r] + delta),]
#         lpr_eval = lhyper
#       }
#
#       lpars_post <- update_pars(lp = lpars_prior$mu,
#                                 ivt = ivt_eval,
#                                 dat = dat_eval,
#                                 lpr = lpr_eval,
#                                 gamma = gamma_eval,
#                                 E0 = E0)
#     } else{
#       lpars_post <- lpars_prior
#     }
#
#     # save values for comparison of appoximation methods
#     # lpars_post_list[[r]] <- lpars_post
#
#     if(plot_path){
#       tseq <- seq(0, nmin, 1/60)
#       plot(dat_obs$time, dat_obs$bis, xlim = c(0,nmin), ylim = c(0,100),
#            col = ifelse(dat_obs$timedelay <=update_tms[min(r+1, length(update_tms))], rgb(red = 0, green = 0, blue = 1, alpha = 0.2), rgb(red = 0, green = 0, blue = 0, alpha = 0.2)),
#            pch = 19,
#            ylab = "bis", xlab = "min")
#       lines(dat_obs$time, dat_obs$bis_t)
#       lines(tseq, Emax(pars = beta_est, ce = tseq, E0 = E0, Emx = E0-bis_target), col = 2, lty = 2)
#       abline(h = 50)
#
#       # plot predicted bis
#       solp <- pk_solution_3cpt_metab(pars = exp(lpars_prior$mu[1:9]), ivt = c(ivt_eval,ivt_new[-1]), init = c(0,0,0,0))
#       lines(tseq, Emax1(pars = exp(lpars_post$mu[10]), ce = solp(tseq)[4,], gamma = gamma_eval, E0 = E0), col = 4, lty = 2)
#       inf_tms <- sapply(ivt_eval, `[[`, "begin")
#       inf_amt <- sapply(ivt_eval, `[[`, "k_R")
#       for(j in 1:r){
#         rug(x = inf_tms[j], quiet = T, ticksize = -1/8*(inf_amt[j]/max(inf_amt[1:r])),
#             side = 3, lwd = 0.5, col = 4, outer = F, line = 0)
#       }
#     }
#     lpars_prior <- lpars_post
#   }
#   out <- list(lpars_post = lpars_post, dat = dat_obs, ivt = ivt_eval)
#   return(out)
# }
#
#
#
#
# minimax_induct_eleveld <- function(pars0, lhyper, target_update_tms, alpha = 0.05, p = NULL, delta = 1/6,
#                                    nmin = 10, bis_target = 50, bis_lb = 40, plot_path = T, update_pars_logical = T,
#                                    qmc = T, laplace_appx = T, seed = NULL, laplace_update_tm = 4, eps = 1,
#                                    target_range = F){
#   library(mvtnorm, quietly = T)
#   library(numDeriv, quietly = T)
#   library(nloptr, quietly = T)
#   if(!is.null(seed)) set.seed(seed)
#   pars_pk0 <- pars0$pk
#   pars_pd0 <- pars0$pd
#   sigma_bis0 <- pars0$sigma
#
#   # assumed to not vary by individual
#   lag0 = pars0$lag
#   gamma1 = pars0$gamma
#   gamma2 = pars0$gamma2
#   E0 = pars0$E0
#
#   # functions
#   max_ce <- function(pars, I0, init, tmax_search = 20, grid_space = 1/60, B = NULL, E = NULL){
#     if(is.null(B) & is.null(E)){
#       # course without infusion - use current concentration
#       B <- pk_solution_3cpt_metab_singleinf(pars = pars, ivt = list(begin = 0, end = 1/6, k_R = 0), init = init, ce_only = T)
#       # course with infusion starting from 0 concentration
#       E <- pk_solution_3cpt_metab_singleinf(pars = pars, ivt = list(begin = 0, end = 1/6, k_R = 1), init = c(0,0,0,0), ce_only = T)
#     }
#     tms <- seq(grid_space,tmax_search,grid_space)
#     ceproj <- B(tms) + I0*E(tms)
#     maxce <- max(ceproj)
#     jpeak = tms[which(ceproj == maxce)]
#     return(c(ce = maxce, jpeak = jpeak))
#   }
#
#   overshoot_objfn <- function(lkR, esamples, init, pars_fixed, alpha = 0.25, B_samples = NULL, E_samples = NULL){
#
#     kR <- exp(lkR)
#     if(is.null(B_samples) | is.null(E_samples)){
#       B_samples <- apply(esamples, 1, pk_solution_3cpt_metab_singleinf, ivt = list(begin = 0, end = 1/6, k_R = 0), init = init, ce_only = T)
#       E_samples <- apply(esamples, 1, pk_solution_3cpt_metab_singleinf, ivt = list(begin = 0, end = 1/6, k_R = 1), init = c(0,0,0,0), ce_only = T)
#     }
#     ce_samples <- t(sapply(1:nrow(esamples), function(i){
#       max_ce(pars = esamples[i,], I0 = kR, init = init_d, B = B_samples[[i]], E = E_samples[[i]])
#     }))
#
#     min_bis_p <- unlist(sapply(1:nrow(esamples), function(i){
#       unname(Emax1(pars = esamples[i,10], ce = ce_samples[i,1], gamma = pars_fixed["gamma"], E0 = pars_fixed["E0"]))
#     }))
#
#     prop_overshoot <- mean(min_bis_p < bis_lb)
#     return(prop_overshoot - alpha)
#   }
#
#   update_tms <- seq(0, nmin-delta, delta)
#   # target_update_ix <- which(update_tms %in% target_update_tms)
#   # next_update_tm <- c(target_update_tms[-1],10)
#   dat_obs <- as.data.frame(matrix(NA, ncol = 12, nrow = 0))
#   lpars_prior <- lhyper
#   # lpars_post_list <- matrix(NA, nrow = length(update_tms), ncol = length(lpars_prior$mu))
#
#   for(r in 1:length(update_tms)){
#     print(r)
#     pars_pk_current <- exp(lpars_prior$mu[1:9])
#     pars_pd_current <- exp(lpars_prior$mu[10])
#
#     # Set initial values
#     if(nrow(dat_obs) == 0) {
#       init_0 <- init_p <- c(0,0,0,0) # predicted and true initial concentrations
#       bis_p <- E0
#       t0 = 0
#       gamma_eval = gamma1
#       ivt_eval <- NULL
#     } else{
#       # true concentrations
#       init_0 <- as.numeric(dat_obs[nrow(dat_obs),c("c1_t","c2_t","c3_t","ce_t")])
#       # current time of evaluation
#       t0 <- ivt_eval[[length(ivt_eval)]]$end
#       # predicted concentration given infusion schedule administered
#       init_p <- c(t(pk_solution_3cpt_metab(pars = pars_pk_current, ivt = ivt_eval, init = c(0,0,0,0))(t0)))
#       # gamma value isn't updated, but switches based on whether the concentration is greater than C50
#       gamma_eval <- unname(ifelse(init_p[4] <= pars_pd_current, gamma1, gamma2))
#       # currently predicted BIS
#       bis_p <- Emax1(pars = pars_pd_current, ce = init_p[4], gamma = gamma_eval, E0 = E0)
#     }
#
#     samples <- quad_pars_qmc(lpars = lpars_prior$mu, sig = lpars_prior$sig, npoints = 100)
#     quad <- list(points = samples, weights = rep(1,nrow(samples)))
#
#     # calculate infusion necessary to reach 50 at point estimate
#     ce_target <- Hinv(pars = pars_pd_current, bis = bis_target, E0 = E0, gamma = gamma_eval)
#     kR_50 <- TCI_basic(Ce = ce_target, pars = pars_pk_current, init = init_p)
#
#     # if within range, target bis = 50
#     if(bis_p < 60 & bis_p > 40 & target_range) kR_alpha <- kR_50
#     else{
#       if(kR_50 > 0){ # min infusion is 0 --> don't calculate unless necessary.
#         # calculate infusion that overshoots for alpha percent
#         B_samples <- apply(exp(samples), 1, pk_solution_3cpt_metab_singleinf, ivt = list(begin = 0, end = 1/6, k_R = 0), init = init_p, ce_only = T)
#         E_samples <- apply(exp(samples), 1, pk_solution_3cpt_metab_singleinf, ivt = list(begin = 0, end = 1/6, k_R = 1), init = c(0,0,0,0), ce_only = T)
#
#         overshoot_lb <- overshoot_objfn(1, esamples = exp(samples), init = init_p,
#                                         B_samples = B_samples, E_samples = E_samples,
#                                         pars_fixed = c(E0 = E0, gamma = gamma_eval), alpha = alpha)
#         if(overshoot_lb < 0){
#           # f <- function(ld) -0.20489 + 0.039*ld
#           # alpha_test <- f(-log(det(lpars_prior$sig)))
#           # print(paste("alpha:",alpha_test))
#           # print(paste("log|D|:",-log(det(lpars_prior$sig))))
#           alpha_test <- alpha
#           kR_alpha <- exp(uniroot(f = overshoot_objfn, interval = c(1,10), esamples = exp(samples), init = init_p,
#                                   B_samples = B_samples, E_samples = E_samples, pars_fixed = c(E0 = E0, gamma = gamma_eval), alpha = alpha_test)$root)
#         } else kR_alpha <- 0
#       } else kR_alpha <- 0
#     }
#
#
#     ivt_new <- list(list(begin = t0, end = t0 + delta, k_R = min(kR_alpha, kR_50)))
#
#     ivt_eval <- c(ivt_eval, ivt_new)
#     dat_i <- sample_data(ivt_d = list(list(begin = 0, end = delta, k_R = ivt_eval[[length(ivt_eval)]]$k_R)),
#                          pars_pk_tci = pars_pk_current, pars_pk0 = pars_pk0, pars_pd0 = pars_pd0, sigma_bis0 = sigma_bis0,
#                          gamma1 = gamma1, gamma2 = gamma2, E0 = E0, delay = lag0, start_time = update_tms[r],
#                          init_pred = init_p, init_true = init_0, nmin_sample = delta)
#     dat_obs <- rbind(dat_obs, dat_i)
#
#     if(update_pars_logical & update_tms[min(r+1, length(update_tms))] >= lag0/60){
#       if(laplace_appx & update_tms[r] > laplace_update_tm){
#         dat_eval = dat_obs[dat_obs$timedelay <=(update_tms[r] + delta) & dat_obs$timedelay > update_tms[r],]
#         lpr_eval = lpars_prior
#       } else{
#         dat_eval = dat_obs[dat_obs$timedelay <= (update_tms[r] + delta),]
#         lpr_eval = lhyper
#       }
#
#       lpars_post <- update_pars(lp = lpars_prior$mu,
#                                 ivt = ivt_eval,
#                                 dat = dat_eval,
#                                 lpr = lpr_eval,
#                                 gamma = gamma_eval,
#                                 E0 = E0)
#     } else{
#       lpars_post <- lpars_prior
#     }
#
#     if(plot_path){
#       tseq <- seq(0, nmin, 1/60)
#       plot(dat_obs$time, dat_obs$bis, xlim = c(0,nmin), ylim = c(0,100),
#            col = ifelse(dat_obs$timedelay <=update_tms[min(r+1, length(update_tms))], rgb(red = 0, green = 0, blue = 1, alpha = 0.2), rgb(red = 0, green = 0, blue = 0, alpha = 0.2)),
#            pch = 19,
#            ylab = "bis", xlab = "min")
#       lines(dat_obs$time, dat_obs$bis_t)
#       abline(h = 50)
#
#       # plot predicted bis
#       solp <- pk_solution_3cpt_metab(pars = exp(lpars_prior$mu[1:9]), ivt = c(ivt_eval,ivt_new[-1]), init = c(0,0,0,0))
#       lines(tseq, Emax1(pars = exp(lpars_post$mu[10]), ce = solp(tseq)[4,], gamma = gamma_eval, E0 = E0), col = 4, lty = 2)
#       inf_tms <- sapply(ivt_eval, `[[`, "begin")
#       inf_amt <- sapply(ivt_eval, `[[`, "k_R")
#       for(j in 1:r){
#         rug(x = inf_tms[j], quiet = T, ticksize = -1/8*(inf_amt[j]/max(inf_amt[1:r])),
#             side = 3, lwd = 0.5, col = 4, outer = F, line = 0)
#       }
#     }
#     lpars_prior <- lpars_post
#   }
#   out <- list(lpars_post = lpars_post, dat = dat_obs, ivt = ivt_eval)
#   return(out)
# }
#
#
# plot_res <- function(res_list, lpars_fixed){
#   nplots <- length(res_list)
#   par(mar=c(1,1,1,1), mfrow = c(nplots/2,2))
#   pars_fixed <- exp(lpars_fixed)
#
#   plot_person <- function(res, xaxt_val){
#     lpars <- res$lpars_post
#     dat <- res$dat
#     ivt <- res$ivt
#     nmin <- max(dat$time)
#     tms <- dat$time
#
#     sol.pr <- pk_solution_3cpt_metab(pars = exp(lpars$mu[1:9]), ivt = ivt, init = c(0,0,0,0))
#     con.pr <- sol.pr(tms)
#     gamma_eval <- ifelse(con.pr[4,] <= exp(lpars$mu[10]), pars_fixed["gamma"], pars_fixed["gamma2"])
#
#     bis.pr <- Emax1(pars = exp(lpars$mu[10]), ce = con.pr[4,], gamma = gamma_eval, E0 = pars_fixed["E0"])
#
#     plot(dat$time, dat$bis_t, type = "l", xlim = c(0,nmin), ylim = c(0,100), xaxt = xaxt_val); abline(h = 50, lty=2)
#     polygon(x = rep(par()$usr[c(1,2)], each = 2), y = c(40,60,60,40), col=rgb(.75,.75,.75, alpha = 0.1))
#     lines(dat$time, bis.pr, col = 4)
#   }
#
#   for(i in 1:nplots){
#     xaxt_val <- ifelse(i %in% c(nplots, nplots-1), "s", "n")
#     plot_person(res_list[[i]], xaxt_val = xaxt_val)
#     # title(sub = paste("Patient", i), cex = 0.8, adj= 1, line = -1)
#   }
#
#   par(mar=c(5.1,4.1,4.1,2.1))
#   par(mfrow = c(1,1))
# }
#
#
# # extract summary statistics from simulation object
# res_stats <- function(sim, lb = 40, ub = 60){
#   tms <- sim[[1]]$dat$time
#   bis_vals <- sapply(sim, function(x) x$dat$bis_t)
#   first_pass <- apply(bis_vals, 2, function(x) tms[which(x<ub)[1]])
#   overshoot <- apply(bis_vals, 2, function(x) tms[which(x<lb)[1]])
#   max_overshoot <- apply(bis_vals, 2, min)
#   max_overshoot[max_overshoot>lb] <- NA
#   second_pass <- apply(bis_vals, 2, function(x) {
#     tm_overshoot <- tms[which(x<lb)[1]]
#     tms[intersect(which(x>lb), which(tms > tm_overshoot))][1]}
#   )
#   stable_entry <- ifelse(is.na(max_overshoot), first_pass, second_pass)
#   return(cbind(first_pass, max_overshoot, second_pass, stable_entry))
# }
#
# rel_inf <- function(sim_ref, sim_test){
#   inf_ref <- sapply(sim_ref, function(x) sum(sapply(x$ivt, `[[`, "k_R")))
#   inf_test <- sapply(sim_test, function(x) sum(sapply(x$ivt, `[[`, "k_R")))
#   return(inf_test/inf_ref)
# }
