#' ------------------------------------------------------------------------------------
#' Closed-loop ------------------------------------------------------------------------
#' ------------------------------------------------------------------------------------


#' Simulate closed-loop control
#'
#' Simulate closed-loop control using Bayesian updates
#'
#' @param pkmod_prior `pkmod` object describing a PK/PK-PD model that is used to calculate
#' TCI infusion rates and is updated as data are simulated and incorporated. Must have an
#' associated Omega matrix.
#' @param pkmod_true `pkmod` object describing the patient's "true" response. This model
#' will be used to simulate observations.
#' @param target_vals A vector of numeric values indicating PK or PD targets for TCI algorithm.
#' @param target_tms A vector of numeric values indicating times at which the TCI algorithm should
#' begin targeting each value.
#' @param obs_tms Times at which data values should be simulated from `pkmod_true`.
#' @param update_tms Times at which `pkmod_prior` should be updated using all available
#' simulated observations.
#' @param type Type of TCI algorithm to be used. Options are "plasma" and "effect".
#' Defaults to "effect". Will be overwritten if `custom_alg` is non-null.
#' @param custom_alg Custom TCI algorithm to overwrite default plasma- or effect-site targeting.
#' @param resp_bounds Optional vector of two values indicating minimum and maximum values possible for the response.
#' @param delay Optional numeric value indicating a temporal delay between when observations
#' are simulated and when they should be made available for updating `pkmod_prior`. For example,
#' a delay should be set to account for a processing time delay in Bispectral Index measurements
#' or the time required to measure drug concentrations from collected samples.
#' @param seed An integer used to initialize the random number generator.
#' @examples
#' prior_vcov <- matrix(diag(c(0.265,0.346,0.209,0.610,0.565,0.597,0.702,0.463)),
#' 8,8, dimnames = list(NULL,c('cl','q2','q3','v','v2','v3','ke0','sigma_add')))
#' pkmod_prior <- pkmod(pars_pk = c(cl = 10, q2 = 2, q3 =20, v = 15, v2 = 30, v3 = 50, ke0 = 1.2),
#' sigma_add = 0.2, log_response = TRUE, Omega = prior_vcov)
#' pkmod_true  <- pkmod(pars_pk = c(cl = 16, q2 = 4, q3 =10, v = 20, v2 = 20, v3 = 80, ke0 = 0.8),
#' sigma_add = 0.1, log_response = TRUE)
#' target_vals <- c(2,3,4,3,3)
#' target_tms <- c(0,5,10,36,60)
#' obs_tms <- c(1,2,4,8,12,16,24,36,48)
#' update_tms <- c(5,15,25,40,50)
#' sim <- clc(pkmod_prior, pkmod_true, target_vals, target_tms, obs_tms, update_tms, seed = 200)
#' len <- 500
#' tms <- seq(0,60,length.out = len)
#' # true, prior, posterior plasma concentrations
#' df <- data.frame(time = rep(tms,3),
#'                  value = c(predict(pkmod_true, sim$inf,tms)[,1],
#'                  predict(pkmod_prior, sim$inf,tms)[,1],
#'                  predict(sim$pkmod_post, sim$inf, tms)[,1]),
#'                  type = c(rep("true",len),rep("prior",len),rep("posterior",len)))
#' library(ggplot2)
#' ggplot(df, aes(x = time, y = value, color = type)) +
#'   geom_line() +
#'   geom_point(data = sim$obs, aes(x = time, y = obs), inherit.aes = FALSE) +
#'   geom_step(data = data.frame(time = target_tms, value = target_vals),
#'   aes(x = time, y = value), inherit.aes = FALSE)
#'
#'
#' # PK-PD example with observation delay (30 sec)
#' prior_vcov <- matrix(diag(c(0.265,0.346,0.209,0.702,0.242,0.230)),6,6,
#' dimnames = list(NULL,c('cl','q2','q3','ke0','c50','sigma_add')))
#' pkpdmod_prior <- update(pkmod_prior, pars_pd = c(c50 = 2.8, gamma = 1.47, e0 = 93, emx = 93),
#' pdfn = emax, pdinv = emax_inv, sigma_add = 4, log_response = FALSE, Omega = prior_vcov)
#' pkpdmod_true <- update(pkmod_true, pars_pd = c(c50 = 3.4, gamma = 1.47, e0 = 93, emx = 93),
#' pdfn = emax, pdinv = emax_inv, sigma_add = 3, log_response = FALSE)
#' target_vals <- c(75,60,50,50)
#' target_tms <- c(0,3,6,10)
#' obs_tms <- seq(1/6,10,1/6)
#' update_tms <- seq(1,10,0.5)
#' sim_pkpd <- clc(pkpdmod_prior, pkpdmod_true, target_vals, target_tms, obs_tms,
#' update_tms, seed = 201, delay = 0.5)
#' # plot results
#' tms <- seq(0,10,length.out = len)
#' df <- data.frame(time = rep(tms,3),
#'                  value = c(predict(pkpdmod_true, sim_pkpd$inf,tms)[,"pdresp"],
#'                  predict(pkpdmod_prior, sim_pkpd$inf,tms)[,"pdresp"],
#'                  predict(sim_pkpd$pkmod_post, sim_pkpd$inf, tms)[,"pdresp"]),
#'                  type = c(rep("true",len),rep("prior",len),rep("posterior",len)))
#' library(ggplot2)
#' ggplot(df, aes(x = time, y = value, color = type)) +
#'   geom_line() +
#'   geom_point(data = sim_pkpd$obs, aes(x = time, y = obs), inherit.aes = FALSE) +
#'   labs(x = "Hours", y = "Bispectral Index") + theme_bw() +
#'   geom_vline(xintercept = update_tms, linetype = "dotted", alpha = 0.6) +
#'   geom_step(data = data.frame(time = target_tms, value = target_vals),
#'   aes(x = time, y = value), inherit.aes = FALSE)
#' @export
clc <- function(pkmod_prior, pkmod_true, target_vals, target_tms, obs_tms, update_tms,
                         type = c("effect","plasma"), custom_alg = NULL,
                         resp_bounds = NULL, delay = 0, seed = NULL){

  if(is.null(seed)) seed <- sample(1:1e5,1)
  set.seed(seed)
  type <- match.arg(type)

  if(length(target_vals) != length(target_tms))
    stop("'target_vals' and 'target_tms' must have the same length")

  if(is.null(pkmod_prior$Omega))
    stop("'pkmod_prior' must have an associated Omega matrix of random effects")

  update_parnms <- colnames(pkmod_prior$Omega)
  update_tms <- union(c(0,update_tms),max(target_tms))
  lpars <- matrix(NA, nrow = length(update_tms)-1, ncol = length(update_parnms), dimnames = list(NULL, update_parnms))

  init_prior <- matrix(NA, nrow = length(update_tms), ncol = pkmod_prior$ncmpt,
                      dimnames = list(NULL, paste0("c",1:pkmod_prior$ncmpt)))
  init_true <- matrix(NA, nrow = length(update_tms), ncol = pkmod_true$ncmpt,
                      dimnames = list(NULL, paste0("c",1:pkmod_true$ncmpt)))

  init_prior[1,] <- pkmod_prior$init
  init_true[1,] <- pkmod_true$init
  inf_all <- obs_all <- obs_tms_obs <- NULL
  vcovi <- vector("list", length(update_tms))

  # combine update times into set of target times
  na.locf <- function(x) {
    v <- !is.na(x)
    c(NA, x[v])[cumsum(v)+1]
  }
  targets_new <- data.frame(time = sort(union(target_tms, update_tms)))
  targets_new <- merge(targets_new, cbind(time = target_tms, value = target_vals), all.x = TRUE)
  targets_new$value <- na.locf(targets_new$value)

  for(i in 2:length(update_tms)){

    # subset targets and observation times to period being updated
    targets_sub <- targets_new[targets_new$time <= update_tms[i] & targets_new$time >= update_tms[i-1],]
    # infusions up until update time
    infi <- inf_tci(pkmod_prior, targets_sub[,"value"], targets_sub[,"time"], type = type, custom_alg = custom_alg, inittm = update_tms[i-1])[,1:3]
    inf_all <- rbind(inf_all, infi)

    # full infusion schedule in case any(obs_tms_sub < min(infi[,"begin"]))
    obs_tms_sub <- obs_tms[obs_tms+delay <= update_tms[i] & obs_tms+delay > update_tms[i-1]]
    obsi <- simulate(pkmod_true, inf = inf_all, tms = obs_tms_sub, init = init_true[1,], resp_bounds=resp_bounds)
    obs_tms_obs <- c(obs_tms_obs, obs_tms_sub)
    obs_all <- rbind(obs_all, obsi)
    # update model parameters
    lpars[i-1,] <- with(pkmod_prior, log(c(pars_pk,pars_pd,sigma_add=sigma_add,sigma_mult=sigma_mult)))[update_parnms]
    pkmod_prior <- bayes_update(lpars[i-1,], pkmod_prior, inf_all, tms = obs_tms_obs,
                                  obs = obs_all, init = init_prior[1,])

    vcovi[[i]] <- pkmod_prior$Omega

    # update true and predicted initial values
    init_prior[i,] <- predict(pkmod_prior, inf = inf_all, tms = update_tms[i], init = init_prior[1,])[,paste0("c",1:pkmod_prior$ncmpt)]
    init_true[i,]  <- predict(pkmod_true, inf = inf_all, tms = update_tms[i], init = init_true[1,])[,paste0("c",1:pkmod_true$ncmpt)]
    pkmod_prior <- update(pkmod_prior, init = init_prior[i,])
    pkmod_true <- update(pkmod_true, init = init_true[i,])
  }

  out <- list(obs = data.frame(time = obs_tms_obs, obs = obs_all),
              inf = inf_all,
              pars = exp(lpars),
              pkmod_post = update(pkmod_prior, init = init_prior[1,]),
              seed = seed)

  return(out)
}


#' Update PK-PD model parameters from observed data values
#'
#' Function will update parameters of `pkmod` object based on available data.
#' Parameters from `pkmod` are used as initial values in the search.
#'
#' @param lpars Logged parameter values. Can be a subset of the full set of PK or PK-PD parameter values.
#' @param pkmod `pkmod` object. Mean values are a subset of log(pars_pk), log(pars_pd),
#' log(sigma_add), log(sigma_mult). PK-PD parameter values not specified in `lpars` will be inferred from `pkmod`.
#' @param inf Infusion schedule
#' @param tms Times associated with observations
#' @param obs Observed values (concentrations or PD response values)
#' @param update_init Logical. Should initial values be updated in addition to parameter values?
#' @param ... Arguments passed to update.pkmod
#' @return `pkmod` object with PK/PK-PD/sigma parameters updated based on minimizing negative logged posterior.
#' @examples
#' # evaluate negative log posterior at a new set of parameters
#' lpars = log(c(cl=11,q2=3,q3=25,v=20,v2=40,v3=80,ke0=1.15,sigma_add=0.3))
#' prior_vcov <- matrix(diag(c(0.265,0.346,0.209,0.610,0.565,0.597,0.702,0.463)),
#' 8,8, dimnames = list(NULL,names(lpars)))
#' my_mod <- pkmod(pars_pk = c(cl = 10, q2 = 2, q3 =20, v = 15, v2 = 30, v3 = 50, ke0 = 1.2),
#' sigma_add = 0.2, log_response = TRUE, Omega = prior_vcov)
#' inf <- inf_manual(inf_tms = 0, inf_rate = 80, duration = 2)
#' tms <- c(1,2,4,8,12)
#' obs <- simulate(my_mod, inf = inf, tms = tms)
#' bayes_update(lpars, my_mod, inf, tms, obs)
#'
#' # evaluate log-prior for subset of parameters (remove volume parameters)
#' lpars_sub = log(c(cl=11,q2=3,q3=25,ke0=1.15,sigma_add=0.15))
#' my_mod <- update(my_mod, Omega = matrix(diag(c(0.265,0.346,0.209,0.702,0.463)),5,5,
#'  dimnames = list(NULL,names(lpars_sub))))
#' bayes_update(lpars_sub, my_mod, inf, tms, obs)
#'
#' # add a pd response and replace multiplicative error with additive error
#' my_mod_pd <- update(my_mod, pars_pd = c(c50 = 2.8, gamma = 1.47, e0 = 93, emx = 93),
#' pdfn = emax, pdinv = emax_inv, ecmpt = 4, sigma_mult = 0, sigma_add = 4)
#' # simulate observations
#' obs_pd <- simulate(my_mod_pd, inf = inf, tms = seq(0,12,0.5))
#' # evaluate likelihood at new parameters
#' lpars_pd_eval = log(c(cl=11,q3=25,v=15,ke0=1.15,sigma_add=4,c50=5,gamma=1))
#' prior_vcov_pd <- matrix(diag(c(0.265,0.209,0.610,0.702,0.230,0.242,0.1)),7,7,
#' dimnames = list(NULL,names(lpars_pd_eval)))
#' my_mod_pd <- update(my_mod_pd, Omega = prior_vcov_pd)
#' bayes_update(lpars_pd_eval, my_mod_pd, inf, tms, obs_pd)
#' @export
bayes_update <- function(lpars, pkmod, inf, tms, obs, update_init = FALSE, ...){

  pkmod <- update(pkmod, ...)
  opt <- try(optim(par = lpars, fn = log_posterior_neg, method = "BFGS", hessian = TRUE,
                    pkmod = pkmod, inf = inf, tms = tms, obs = obs, control = list(maxit = 5000)))

  if(class(opt) == "try-error") browser()
  if(opt$convergence>0) warning(paste("Convergence:", opt$convergence))

  vcov_est <- solve(opt$hessian)
  if(any(diag(vcov_est)<0)){
    vcov_est <- pkmod$Omega
  }
  pars_est <- exp(opt$par)
  pkmod_new <- update(pkmod,
                      pars_pk = pars_est[names(pars_est) %in% names(pkmod$pars_pk)],
                      Omega = vcov_est)

  if(any(names(pars_est) %in% names(pkmod$pars_pd)))
    pkmod_new <- update(pkmod_new, pars_pd = pars_est[names(pars_est) %in% names(pkmod$pars_pd)])
  if("sigma_add" %in% names(pars_est))
    pkmod_new <- update(pkmod_new, sigma_add = unname(pars_est["sigma_add"]))
  if("sigma_mult" %in% names(pars_est))
    pkmod_new <- update(pkmod_new, sigma_mult = unname(pars_est["sigma_mult"]))

  if(update_init){
    con_nms <- paste0("c",1:pkmod_new$ncmpt)
    ini <- predict(pkmod_new, inf = inf, tms = inf[nrow(inf),"end"])[,con_nms]
    pkmod_new <- update(pkmod_new, init = ini)
  }

  return(pkmod_new)
}

#' Evaluate the log likelihood of a vector of parameter values
#'
#' Evaluate the log liklihood of parameters given observed data. Can be applied to PK or PK-PD models.
#'
#' @param lpars Named vector of logged parameter values to be evaluated. This should include any PK or PD parameters,
#' as well as residual error standard deviations (sigma_add or sigma_mult) that are to be evaluated.
#' @param pkmod `pkmod` object. Mean values are a subset of log(pars_pk), log(pars_pd),
#' log(sigma_add), log(sigma_mult). PK-PD parameter values not specified in `lpars` will be inferred from `pkmod`.
#' @param inf Infusion schedule
#' @param tms Times associated with observations
#' @param obs Observed values (concentrations or PD response values)
#' @return Numeric value of length 1
#' @examples
#' my_mod <- pkmod(pars_pk = c(cl = 10, q2 = 2, q3 =20, v = 15, v2 = 30, v3 = 50,
#'  ke0 = 1.2), sigma_mult = 0.2)
#' inf <- inf_manual(inf_tms = 0, inf_rate = 80, duration = 2)
#' tms <- c(1,2,4,8,12)
#' obs <- simulate(my_mod, inf = inf, tms = tms)
#' # evaluate log-likelihood at a new set of parameters
#' lpars = log(c(cl=11,q2=3,q3=25,v=15,v2=30,v3=50,ke0=1.15,sigma_mult=0.3))
#' log_likelihood(lpars, my_mod, inf, tms, obs)
#'
#' # estimate for a subset of parameters (exclude q2, v2, v3)
#' lpars_sub = log(c(cl=11,q3=25,v=15,ke0=1.15,sigma_mult=0.3))
#' log_likelihood(lpars_sub, my_mod, inf, tms, obs)
#'
#' # add a pd response and replace multiplicative error with additive error
#' my_mod_pd <- update(my_mod, pars_pd = c(c50 = 2.8, gamma = 1.47, e0 = 93,
#' emx = 93), pdfn = emax, pdinv = emax_inv, ecmpt = 4, sigma_mult = 0, sigma_add = 4)
#' # simulate observations
#' obs_pd <- simulate(my_mod_pd, inf = inf, tms = seq(0,12,0.5))
#' # evaluate likelihood at new parameters
#' lpars_pd <- log(c(cl=11,q3=25,v=15,ke0=1.15,sigma_add=4,c50=5,gamma=1))
#' log_likelihood(lpars_pd, my_mod_pd, inf, tms = seq(0,12,0.5), obs_pd)
#' @export
log_likelihood <- function(lpars, pkmod, inf, tms, obs){

  # assign new parameters
  pars <- exp(lpars)
  nms_new <- names(pars)
  pknew <- pars[nms_new %in% names(pkmod$pars_pk)]
  pdnew <- pars[nms_new %in% names(pkmod$pars_pd)]

  sigma_add = ifelse(any(nms_new == "sigma_add"), pars["sigma_add"], pkmod$sigma_add)
  sigma_mult = ifelse(any(nms_new == "sigma_mult"), pars["sigma_mult"], pkmod$sigma_mult)

  # obs_cmpt if PD response
  if(!is.null(pkmod$pdfn)){
    obs_cmpt <- pkmod$ncmpt +1
  } else{
    obs_cmpt <- pkmod$pcmpt
  }

  # predict concentrations at tms
  resp <- predict(pkmod, inf, tms, pars_pk = pknew, pars_pd = pdnew)[,obs_cmpt]
  if(pkmod$log_response){
    obs <- log(obs)
    resp <- log(resp)
  }

  return(sum(dnorm(obs, resp, sd = resp*sigma_mult + sigma_add, log = TRUE)))
}

#' Calculate logged-prior probability for a set of parameters
#'
#' Calculate logged-prior probability for a set of parameters, assuming that parameter values
#' are log-normally distributed. Mean values are set as the logged parameter values in
#' the `pkmod` object. Variances are given by the diagonal elements of `prior_vcov`.
#'
#' @param lpars Logged parameter values. Can be a subset of the full set of PK or PK-PD parameter values.
#' @param pkmod `pkmod` object. Mean values are a subset of log(pars_pk), log(pars_pd),
#' log(sigma_add), log(sigma_mult). PK-PD parameter values not specified in `lpars` will be inferred from `pkmod`.
#' @return Numeric value of length 1
#' @examples
#' # evaluate log-prior for pk parameters + residual
#' lpars = log(c(cl=11,q2=3,q3=25,v=15,v2=30,v3=50,ke0=1.15,sigma_add=0.15))
#' prior_vcov <- matrix(diag(c(0.265,0.346,0.209,0.610,0.565,0.597,0.702,0.463)), 8,8,
#' dimnames = list(NULL,names(lpars)))
#' my_mod <- pkmod(pars_pk = c(cl = 10, q2 = 2, q3 =20, v = 15, v2 = 30, v3 = 50, ke0 = 1.2),
#' sigma_add = 0.2, log_response = TRUE, Omega = prior_vcov)
#' log_prior(lpars, my_mod)
#'
#' # evaluate log-prior for subset of parameters (remove volume parameters)
#' lpars_sub = log(c(cl=11,q2=3,q3=25,ke0=1.15,sigma_add=0.15))
#' prior_vcov_sub <- matrix(diag(c(0.265,0.346,0.209,0.702,0.463)), 5,5,
#' dimnames = list(NULL,names(lpars_sub)))
#' my_mod <- update(my_mod, Omega = prior_vcov_sub)
#' log_prior(lpars_sub, my_mod)
#' @export
log_prior <- function(lpars, pkmod){

  if(is.null(pkmod$Omega)) stop("pkmod$Omega must be specified to use log_prior")
  if(length(lpars) != nrow(pkmod$Omega)) stop("pkmod$Omega must have row and column dimensions equal to length of lpars")
  if(!all(names(lpars) %in% colnames(pkmod$Omega))) stop("Names of 'lpars' and colnames of pkmod$Omega must be the same")

  # construct mean vector in correct order (pkmod$Omega)
  mu_all <- log(c(pkmod$pars_pk, pkmod$pars_pd, sigma_add = pkmod$sigma_add, sigma_mult = pkmod$sigma_mult))
  mu <- mu_all[colnames(pkmod$Omega)]
  lpars <- lpars[colnames(pkmod$Omega)]

  # evaluate log prior
  vc <- pkmod$Omega
  lp <- mvtnorm::dmvnorm(lpars, mu, vc, log = TRUE)
  while(is.infinite(lp)){
    diag(vc) <- diag(vc) + diag(vc)*0.5
    lp <- mvtnorm::dmvnorm(lpars, mu, vc, log = TRUE)
  }

  return(lp)
}


#' Evaluate the negative log posterior value of a parameter vector
#'
#' Evaluate the negative log posterior value of a parameter vector given a set of
#' observations and prior distribution for log-normally distributed PK or PK-PD parameters.
#'
#' @param lpars Logged parameter values. Can be a subset of the full set of PK or PK-PD parameter values.
#' @param pkmod `pkmod` object. Mean values are a subset of log(pars_pk), log(pars_pd),
#' log(sigma_add), log(sigma_mult). PK-PD parameter values not specified in `lpars` will be inferred from `pkmod`.
#' @param inf Infusion schedule
#' @param tms Times associated with observations
#' @param obs Observed values (concentrations or PD response values)
#' @return Numeric value of length 1
#' @examples
#' # evaluate negative log posterior at a new set of parameters
#' lpars = log(c(cl=11,q2=3,q3=25,v=20,v2=40,v3=80,ke0=1.15,sigma_add=0.3))
#' prior_vcov <- matrix(diag(c(0.265,0.346,0.209,0.610,0.565,0.597,0.702,0.463)),
#' 8,8, dimnames = list(NULL,names(lpars)))
#' my_mod <- pkmod(pars_pk = c(cl = 10, q2 = 2, q3 =20, v = 15, v2 = 30, v3 = 50,
#' ke0 = 1.2), sigma_add = 0.2, log_response = TRUE, Omega = prior_vcov)
#' inf <- inf_manual(inf_tms = 0, inf_rate = 80, duration = 2)
#' tms <- c(1,2,4,8,12)
#' obs <- simulate(my_mod, inf = inf, tms = tms)
#' log_posterior_neg(lpars, my_mod, inf, tms, obs)
#'
#' # evaluate log-prior for subset of parameters (remove volume parameters)
#' lpars_sub = log(c(cl=11,q2=3,q3=25,ke0=1.15,sigma_add=0.15))
#' my_mod <- update(my_mod, Omega = matrix(diag(c(0.265,0.346,0.209,0.702,0.463)),
#' 5,5, dimnames = list(NULL,names(lpars_sub))))
#' log_posterior_neg(lpars_sub, my_mod, inf, tms, obs)
#'
#' # add a pd response and replace multiplicative error with additive error
#' # evaluate likelihood at new parameters
#' lpars_pd_eval = log(c(cl=11,q3=25,v=15,ke0=1.15,sigma_add=4,c50=5,gamma=1))
#' prior_vcov_pd <- matrix(diag(c(0.265,0.209,0.610,0.702,0.230,0.242,0.1)),7,7,
#' dimnames = list(NULL,names(lpars_pd_eval)))
#' my_mod_pd <- update(my_mod, pars_pd = c(c50 = 2.8, gamma = 1.47, e0 = 93, emx = 93),
#' pdfn = emax, pdinv = emax_inv, ecmpt = 4, sigma_mult = 0, sigma_add = 4,
#' Omega = prior_vcov_pd)
#' # simulate observations
#' obs_pd <- simulate(my_mod_pd, inf = inf, tms = seq(0,12,0.5))
#' log_posterior_neg(lpars_pd_eval, my_mod_pd, inf, tms, obs_pd)
#' @export
log_posterior_neg <- function(lpars, pkmod, inf, tms, obs){
  -1*(log_prior(lpars, pkmod) + log_likelihood(lpars, pkmod, inf, tms, obs))
}




#' Simulate closed-loop control using Bayesian updates
#'
#' @param pkmod_prior `pkmod` or `poppkmod` object describing a PK/PK-PD model that is used to calculate
#' TCI infusion rates and is updated as data are simulated and incorporated. Must have an
#' associated Omega matrix.
#' @param pkmod_true `pkmod` or `poppkmod` object describing the patient's "true" response. This model
#' will be used to simulate observations.
#' @param target_vals A vector of numeric values indicating PK or PD targets for TCI algorithm.
#' @param target_tms A vector of numeric values indicating times at which the TCI algorithm should
#' begin targeting each value.
#' @param obs_tms Times at which data values should be simulated from `pkmod_true`.
#' @param update_tms Times at which `pkmod_prior` should be updated using all available
#' simulated observations.
#' @param type Type of TCI algorithm to be used. Options are "plasma" and "effect".
#' Defaults to "effect". Will be overwritten if `custom_alg` is non-null.
#' @param custom_alg Custom TCI algorithm to overwrite default plasma- or effect-site targeting.
#' @param resp_bounds Optional vector of two values indicating minimum and maximum values possible for the response.
#' @param delay Optional numeric value indicating a temporal delay between when observations
#' are simulated and when they should be made available for updating `pkmod_prior`. For example,
#' a delay should be set to account for a processing time delay in Bispectral Index measurements
#' or the time required to measure drug concentrations from collected samples.
#' @param seed An integer used to initialize the random number generator.
#' @param verbose Logical. Print progress as simulation is run.
#' @examples
#' data <- data.frame(ID = 1:3, AGE = c(20,30,40), TBW = c(60,70,80),
#' HGT = c(150,160,170), MALE = c(TRUE,FALSE,TRUE))
#' pkmod_prior <- poppkmod(data, drug = "ppf", model = "eleveld")
#' pkmod_true  <- poppkmod(data, drug = "ppf", model = "eleveld", sample = TRUE)
#' obs_tms <- seq(1/6,10,1/6)
#' update_tms <- c(2,4,6,8)
#' target_vals = c(75,60,50,50)
#' target_tms = c(0,3,6,10)
#' sim <- simulate_clc(pkmod_prior, pkmod_true, target_vals, target_tms, obs_tms,
#' update_tms, seed = 200)
#' len <- 500
#' tms <- seq(0,10,length.out = len)
#' resp <- data.frame(rbind(predict(pkmod_true, sim$inf, tms),
#' predict(pkmod_prior, sim$inf, tms),
#' predict(sim$pkmod_post, sim$inf, tms)))
#' resp$type = c(rep("true",len*3),rep("prior",len*3),rep("posterior",len*3))
#' library(ggplot2)
#' ggplot(resp) + geom_line(aes(x = time, y = pdresp, color = factor(id))) +
#' facet_wrap(~type) + labs(x = "Hours", y = "Bispectral Index") +
#'   geom_step(data = data.frame(time = target_tms, value = target_vals),
#'   aes(x = time, y = value), inherit.aes = FALSE)
#' @export
simulate_clc <- function(pkmod_prior, pkmod_true, target_vals, target_tms, obs_tms, update_tms,
                         type = c("effect","plasma"), custom_alg = NULL,
                         resp_bounds = NULL, delay = 0, seed = NULL, verbose = TRUE){

  if(is.null(seed)) seed <- sample(1:1e5,1)
  set.seed(seed)
  type <- match.arg(type)

  if(class(pkmod_prior) != class(pkmod_true) | !class(pkmod_prior) %in% c("pkmod","poppkmod"))
    stop("pkmod_prior and pkmod_true must either both have class 'pkmod' or 'poppkmod'")

  pkmod_post <- pkmod_prior
  if(class(pkmod_prior) == "pkmod"){
    out <- clc(pkmod_prior, pkmod_true, target_vals, target_tms, obs_tms, update_tms, type, custom_alg, resp_bounds, delay, seed)
    ids = as.factor(1)
  } else{
    if(length(pkmod_prior$pkmods) != length(pkmod_true$pkmods))
      stop("'pkmod_prior' and 'pkmod_true' must have the same number of pkmods")

    sim <- lapply(1:length(pkmod_prior$pkmods), function(i){
      if(verbose) print(paste0("Simulating ID=",pkmod_prior$ids[i]))
      clc(pkmod_prior$pkmods[[i]], pkmod_true$pkmods[[i]], target_vals, target_tms, obs_tms, update_tms, type, custom_alg, resp_bounds, delay)
      })

    pkmod_post$pkmods <- lapply(sim, `[[`, "pkmod_post")

    ids <- as.factor(pkmod_prior$ids)
    out <- list(obs = cbind(id = rep(ids, each = nrow(sim[[1]]$obs)),
                            type = "true",
                     do.call("rbind",lapply(sim, `[[`, "obs"))),
         inf = cbind(id = rep(ids, each = nrow(sim[[1]]$inf)),
                     do.call("rbind",lapply(sim, `[[`, "inf"))),
         pars = cbind(id = rep(ids, each = nrow(sim[[1]]$pars)),
                     do.call("rbind",lapply(sim, `[[`, "pars"))),
         pkmod_post = pkmod_post,
         seed = seed)
  }

  tmseq <- seq(min(out$inf[,"begin"]), max(out$inf[,"end"]), length.out = 1e3)
  resp <- list(prior = predict(pkmod_prior, out$inf,tms = tmseq),
               posterior = predict(pkmod_post, out$inf, tms = tmseq),
               true = predict(pkmod_true, out$inf, tms = tmseq))
  if(length(unique(sapply(resp, ncol)))==1){
    resp <- as.data.frame(do.call("rbind", resp))
    resp$type <- c(rep("prior", 1e3*length(ids)),
                   rep("posterior",1e3*length(ids)),
                   rep("true", 1e3*length(ids)))
    if("id" %in% names(out$resp)) out$resp$id <- as.factor(out$resp$id)
  }
  out <- c(out, list(resp = resp))
  out
}




