# --------------------------------------------------------------------------------------------------------------------------------
# - Simulation functions ---------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------

pal <- brewer.pal(5, "BrBG")

#' Function to simulate data from a specified PK or PK-PD model with a specified infusion schedule.
#' @param inf An infusion rate object outputted from either the 'create_intvl' function or the 'iterate_tci_grid' function
gen_data <- function(inf, pkmod, pars_pk0, sigma.add = 0, sigma.mult = 0, init = NULL, tms = NULL, pdmod = NULL, pars_pd0 = NULL, ecmpt = NULL, delay = 0, max_pdval = NULL, min_pdval = NULL){

  if(any(!(c("infrt","begin","end") %in% names(inf)))) stop('Names of argument "inf" must include ("infrt","begin","end").')

  if(is.null(tms)){
    tms <- inf[,"begin"]
  }

  if(is.null(init))
    init <- eval(formals(pkmod)$init)

  con0 <- as.data.frame(predict(pkmod = pkmod, inf = inf, tms = tms, pars = pars_pk0, init = init))

  # additive and multiplicative errors
  eadd <- rnorm(nrow(con0),0,sigma.add)
  emult <- rnorm(nrow(con0),0,sigma.mult)
  if(is.null(pdmod)){
    con0$cobs <- con0[,"c1"]*(1+emult) + eadd
  } else{
    # possible time delay
    con0$timeobs <- con0$time + delay

    if(is.null(ecmpt))
      ecmpt <- length(eval(formals(pkmod)$init))

    # pd observations
    con0$pd0 <- pdmod(con0[,paste0("c",ecmpt)], pars_pd0)
    con0$pdobs <- con0$pd0*(1+emult) + eadd

    # replace with max/min values if necessary
    con0$pdobs[con0$pdobs > max_pdval] <- max_pdval
    con0$pdobs[con0$pdobs < min_pdval] <- min_pdval
  }

  out <- list(sim = con0,
              inf = inf,
              pkmod = pkmod,
              pdmod = pdmod,
              pars_pk = pars_pk0,
              pars_pd = pars_pd0,
              sigma.add = sigma.add,
              sigma.mult = sigma.mult,
              ecmpt = ecmpt,
              delay = delay)
  class(out) <- c(class(out), "datasim")
  return(out)
}


# function to merge datasim objects from different infusion schedules
# infusion schedules can be passed directly in or as a list.
combine_sim <- function(...){
  simlist <- list(...)
  if(length(simlist) == 1) simlist <- simlist[[1]]

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
  class(out) <- c(class(out), "datasim")
  return(out)
}


plot.datasim <- function(datasim){

  r <- range(datasim$sim$time)
  tms <- seq(r[1], r[2], diff(r)/1000)

  cp <- data.frame(predict(pkmod = datasim$pkmod,
                           inf = datasim$inf,
                           tms = tms,
                           pars = datasim$pars_pk,
                           init = unlist(head(datasim$inf[,c("c1_start","c2_start","c3_start","c4_start")],1))))

  if(is.null(datasim$pdmod)){

    out <- ggplot(cp, aes(x = time, y = c1)) +
      geom_line() +
      geom_point(data = datasim$sim, aes(x = time, y = cobs), shape = 1, col = "blue") +
      labs(x = "Time (min)", y = "Concentration")

  } else{
    cp$pdp <- datasim$pdmod(ce = cp[,paste0("c",datasim$ecmpt)], pars = datasim$pars_pd)

    out <- ggplot(cp, aes(x = time, y = pdp)) +
      geom_line() +
      geom_point(data = datasim$sim, aes(x = time, y = pdobs), shape = 1, col = "blue") +
      labs(x = "Time (min)", y = "Effect")

    # if("pdt" %in% names(datasim$inf)) out <- out + geom_line(data = datasim$inf, aes(x = begin, y = pdt), linetype = "dashed", color = "red")
  }

  out
}
.S3method("plot","datasim",plot.datasim)


#' Function to apply saved population PK or PK-PD models to a dataframe of patient values.
apply_poppk <- function(patient_df, mod = c("marsh","schnider","eleveld"), ...){
  switch(match.arg(mod),
         marsh = marsh_poppk(patient_df, ...),
         schnider = schnider_poppk(patient_df, ...),
         eleveld = eleveld_poppk(patient_df, ...))
}




#' Return an infusion schedule defined by an Emax sigmoid curve.
#' @param beta Parameters (c50, gamma) for Emax target function
#' @param theta Vector of PK-PD parameters (theta_PK, theta_PD)
#' @param ini Initial concentration values
#' @param t0 Starting time
#' @param E0 Effect at ce = 0. Used as starting point for sigmoid curve.
#' @param gamma Slope at c50
#' @param tmx Final time to evaluate infusion schedule to.
#' @param bist Target BIS, used to set the maximum effect value of the target sigmoid curve.
#' @param delta Time interval between TCI updates. Defaults to 1/6 minutes = 10 seconds.
sig_inf <- function(beta, theta, ini, t0, E0, gamma, tmx = 10, bist = 50, delta = 1/6, ...){
  if(t0+delta <= tmx) tm_eval <- seq(t0+delta,tmx,delta)
  else tm_eval <- tmx
  bis_targets <- E0-(E0-bist)*(tm_eval^beta[2] / (tm_eval^beta[2] + beta[1]^beta[2]))
  ce_targets <- Hinv(pars = theta[10], bis = bis_targets, E0 = E0, gamma = gamma)
  kR_future <- TCI_EffectSite(Cet = ce_targets, pars = theta[1:9], init = ini, delta = delta, ...)
  if(t0 != 0){
    kR_vec <- unlist(kR_future)
    kR_vec[names(kR_vec) == "begin"] <- kR_vec[names(kR_vec) == "begin"] + t0
    kR_vec[names(kR_vec) == "end"] <- kR_vec[names(kR_vec) == "end"] + t0
    kR_future <- relist(kR_vec, kR_future)
  }
  return(kR_future)
}


# --------------------------------------------------------------------------------------------------------------------------------
# - Old functions ----------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------

#' Return an infusion schedule defined by an Emax sigmoid curve.
#' @param beta Parameters (c50, gamma) for Emax target function
#' @param theta Vector of PK-PD parameters (theta_PK, theta_PD)
#' @param ini Initial concentration values
#' @param t0 Starting time
#' @param E0 Effect at ce = 0. Used as starting point for sigmoid curve.
#' @param gamma Slope at c50
#' @param tmx Final time to evaluate infusion schedule to.
#' @param bist Target BIS, used to set the maximum effect value of the target sigmoid curve.
#' @param delta Time interval between TCI updates. Defaults to 1/6 minutes = 10 seconds.
sig_inf <- function(beta, theta, ini, t0, E0, gamma, tmx = 10, bist = 50, delta = 1/6, ...){
  if(t0+delta <= tmx) tm_eval <- seq(t0+delta,tmx,delta)
  else tm_eval <- tmx
  bis_targets <- E0-(E0-bist)*(tm_eval^beta[2] / (tm_eval^beta[2] + beta[1]^beta[2]))
  ce_targets <- Hinv(pars = theta[10], bis = bis_targets, E0 = E0, gamma = gamma)
  kR_future <- TCI_EffectSite(Cet = ce_targets, pars = theta[1:9], init = ini, delta = delta, ...)
  if(t0 != 0){
    kR_vec <- unlist(kR_future)
    kR_vec[names(kR_vec) == "begin"] <- kR_vec[names(kR_vec) == "begin"] + t0
    kR_vec[names(kR_vec) == "end"] <- kR_vec[names(kR_vec) == "end"] + t0
    kR_future <- relist(kR_vec, kR_future)
  }
  return(kR_future)
}

#' Objective function for a single infusion - numerically estimates the weighted integral between the time-BIS curve and the target BIS value.
#' @param kR Infusion schedule to be evaluated.
#' @param theta Vector of PK-PD parameters (theta_PK, theta_PD)
#' @param E0 Effect at ce = 0. Used as starting point for sigmoid curve.
#' @param gamma Hill parameter (slope at c50) for PD model
#' @param alpha Weight associated with the integral above the target BIS value and below the time-BIS curve.
#' @param bist Target BIS, used to set the maximum effect value of the target sigmoid curve.
#' @param tfinal Final time to evaluate infusion schedule to.
#' @param dt Resolution of integral
Phi <- function(kR, theta, E0, gamma, t0, alpha = 0.01, bist = 50, tfinal = 10, dt = 1e-2, log = F){
  tms <- seq(0,tfinal,length.out = 1/dt)
  ce <- pk_solution_3cpt_metab(pars = theta[1:9], ivt = kR, init = c(0,0,0,0))(tms)[4,]
  ce[ce < 0] <- 0
  bisp <- Emax1(theta[10],ce = ce, gamma = gamma, E0 = E0)
  phi1 <- sum((bisp[bisp>bist] - bist))*dt*tfinal # integral above target
  phi2 <- sum((bist - bisp[bisp<bist]))*dt*tfinal # integral below target
  if(log) return(log(alpha*phi1 + (1-alpha)*phi2))
  else return(alpha*phi1 + (1-alpha)*phi2)
}


sim_BIS <- function(pars, ivt, init, gamma, E0, BIS_sampling_freq = 1){
  pars_pk <- pars[1:9]
  pars_pd <- pars[10]
  sampling_window <- c(ivt[[1]]$begin, tail(ivt,1)[[1]]$end)
  delta_tms_min <- BIS_sampling_freq / 60
  tm_seq <- seq(sampling_window[1] + delta_tms_min, sampling_window[2], by = delta_tms_min)
  sol <- pk_solution_3cpt_metab(pars_pk, ivt = ivt, init = init)
  con_est <- sol(tm_seq)
  bis_vals <- Emax1(pars = pars_pd, ce = con_est[4,], gamma = gamma, E0 = E0)
  out <- list(bis = bis_vals, con = con_est, bis_tms = tm_seq)
  return(out)
}

quad_pars <- function(lpars, sig, kR, gamma, E0, nquad = 4, ev_dim = 3, prune = 0){
  b <- function(lpars, ...) sim_BIS(pars = exp(lpars), ...)$bis
  grad <- fdGrad(pars = lpars[1:10], fun = b, ivt = kR, init = c(0,0,0,0), gamma = gamma, E0 = E0, BIS_sampling_freq = 1)
  if(any(is.nan(rowSums(grad)))) grad <- grad[-is.nan(rowSums(grad)),]
  ev <- svd(grad)
  sig_ev <- t(ev$v[,c(1:ev_dim)]) %*% sig[1:10,1:10] %*% ev$v[,c(1:ev_dim)]
  quad <- mgauss.hermite(n = nquad, mu = rep(0,ev_dim), sigma = sig_ev, prune = prune)
  return(list(points = t(t(quad$points %*% t(ev$v[,1:ev_dim])) + lpars[1:10]),
              weights = quad$weights,
              var_exp = cumsum(ev$d^2/ sum(ev$d^2))[ev_dim]))
}

# quasi-monte-carlo
quad_pars_qmc <- function(lpars, sig, npoints = 200, par_dim = 10){
  xy = randtoolbox::sobol(npoints,dim=par_dim)
  z = apply(xy,2,qnorm)
  ev <- eigen(sig[1:par_dim,1:par_dim])
  A <- t(t(ev$vectors)*sqrt(ev$values))
  y <- t(t(A)%*%t(z) + lpars[1:par_dim])
  return(y)
}

sample_data <- function(ivt_d, pars_pk_tci, pars_pk0, pars_pd0, sigma_bis0, gamma1, gamma2, E0, start_time, init_pred, init_true, delay, nmin_sample = 1/6, sample_freq = 1/60){
  # times to be evaluated (starting at 0)
  tms0 <- seq(sample_freq, nmin_sample, sample_freq)

  # true concentration is used to generate observed BIS
  sol0 <- pk_solution_3cpt_metab(pars = pars_pk0, ivt = ivt_d, init = init_true)
  con0 <- sol0(tms0)

  # simulate BIS from true concentrations
  gamma0 <- ifelse(con0[4,] <= pars_pd0, gamma1, gamma2)
  bis_t <- Emax1(pars = pars_pd0, ce = con0[4,], gamma = gamma0, E0 = E0)
  bis_obs <- truncnorm::rtruncnorm(n = length(con0[4,]), mean = bis_t, sd = sigma_bis0, a = 0, b = 100)

  # real times being sampled
  tms_sample <- seq(start_time + sample_freq, start_time + nmin_sample, sample_freq)
  tmsobs <- tms_sample + delay/60

  ivt_d[[1]]$begin <- ivt_d[[1]]$begin + start_time
  ivt_d[[1]]$end <- ivt_d[[1]]$end + start_time

  dat0 <- data.frame(cbind(time = tms_sample, timedelay = tmsobs, bis = bis_obs, bis_t = bis_t,
                           c1_t = con0[1,], c2_t = con0[2,], c3_t = con0[3,], ce_t = con0[4,]))
  return(dat0)
}


update_pars <- function(lp, dat, ivt, lpr, gamma, E0){
  post_est <- nlm(f = log_posterior_neg, p = lp, ivt=ivt, dat=dat,  lhyper = lpr, gamma = gamma, E0 = E0,
                  hessian = T,
                  steptol=1e-6, gradtol=1e-6, stepmax = 5,
                  iterlim = 2000)
  post_est_pars <- post_est$estimate
  post_hes <- post_est$hessian
  lpost <- list(mu = post_est_pars, sig = solve(post_hes))
  return(lpost)
}

sigmoid_update_objfn <- function(lgamma, E0, bis_target = 50, nmin = 10, nmin_target = 5, eps = 1, ...){
  # bis_t50 = nmin*(eps/(bis0 - bis_target + eps))^(1/gamma) # express c50 parameter in terms of gamma
  lbis_t50 = log(nmin_target) + 1/exp(lgamma) * (log(eps) - log(E0 - bis_target + eps))
  sig_ED(lbeta = c(lbis_t50, lgamma), E0 = E0, bis_target = bis_target, nmin = nmin, epsilon = 0, ...)
}

sigmoid_update_objfn2 <- function(lbeta, E0, bis_target = 50, nmin = 10, nmin_target = 5, eps = 1, ...){
  sig_ED(lbeta = lbeta, E0 = E0, bis_target = bis_target, nmin = nmin, epsilon = 0, ...)
}


sig_ED <- function(lbeta, ltheta_hat, Theta_samples, weights, gamma, E0, sig, ivt0, t0, t1 = NULL, init_p, p = NULL, alpha = 0.01, epsilon = 0.1, bis_target = 50, nmin = 10, dt = 1e-2){
  if(is.null(t1)) t1 <- nmin
  kR_all <- c(ivt0, sig_inf(beta =exp(lbeta), theta = exp(ltheta_hat), ini = init_p, t0 = t0, tmx = t1, E0 = E0, gamma = gamma))

  if(is.null(p)){
    # evaluate quadrature points at infusion schedule
    obj_fn_val <- c(t(vapply(1:nrow(Theta_samples), function(l){
      Phi(kR = kR_all, theta = Theta_samples[l,], alpha = alpha, bist = bis_target, tfinal = nmin,
          t0 = t0,
          dt = dt, gamma = gamma, E0 = E0)
    }, FUN.VALUE = numeric(1))) %*% weights)
  } else{
    # evaluate quadrature points at infusion schedule
    obj_fn_val <- quantile(vapply(1:nrow(Theta_samples), function(l){
      Phi(kR = kR_all, theta = Theta_samples[l,], alpha = alpha, bist = bis_target, tfinal = nmin,
          t0 = t0,
          dt = dt, gamma = gamma, E0 = E0)
    }, FUN.VALUE = numeric(1)), probs = p)
  }
  return(obj_fn_val)
}



sigmoid_induct_eleveld <- function(pars0, lhyper, target_update_tms, alpha = 0.01, p = NULL, epsilon = 0.1, delta = 1/6,
                                   nmin = 10, nmin_prop = 1/2, beta_init = c(2,3), fixed_bis = F,
                                   bis_target = 50, plot_path = T, delta_bis = 5, update_pars_logical = T,
                                   qmc = F, laplace_appx = T, seed = NULL, laplace_update_tm = 4, eps = 1){
  library(mvtnorm, quietly = T)
  library(numDeriv, quietly = T)
  library(nloptr, quietly = T)
  if(!is.null(seed)) set.seed(seed)
  pars_pk0 <- pars0$pk
  pars_pd0 <- pars0$pd
  sigma_bis0 <- pars0$sigma

  # assumed to not vary by individual
  lag0 = pars0$lag
  gamma1 = pars0$gamma
  gamma2 = pars0$gamma2
  E0 = pars0$E0

  update_tms <- seq(0, nmin-delta, delta)
  target_update_ix <- which(update_tms %in% target_update_tms)
  next_update_tm <- c(target_update_tms[-1],10)
  dat_obs <- as.data.frame(matrix(NA, ncol = 12, nrow = 0))
  lpars_prior <- lhyper
  # lpars_post_list <- matrix(NA, nrow = length(update_tms), ncol = length(lpars_prior$mu))

  for(r in 1:length(update_tms)){
    print(r)
    pars_pk_current <- exp(lpars_prior$mu[1:9])
    pars_pd_current <- exp(lpars_prior$mu[10])

    # Set initial values
    if(nrow(dat_obs) == 0) {
      init_0 <- init_p <- c(0,0,0,0) # predicted and true initial concentrations
      bis_p <- E0
      t0 = 0
      beta_est <- beta_init
      gamma_eval = gamma1
      ivt_eval <- NULL
    } else{
      # true concentrations
      init_0 <- as.numeric(dat_obs[nrow(dat_obs),c("c1_t","c2_t","c3_t","ce_t")])
      # current time of evaluation
      t0 <- ivt_eval[[length(ivt_eval)]]$end
      # predicted concentration given infusion schedule administered
      init_p <- c(t(pk_solution_3cpt_metab(pars = pars_pk_current, ivt = ivt_eval, init = c(0,0,0,0))(t0)))
      # gamma value isn't updated, but switches based on whether the concentration is greater than C50
      gamma_eval <- unname(ifelse(init_p[4] <= pars_pd_current, gamma1, gamma2))
      # currently predicted BIS
      bis_p <- Emax1(pars = pars_pd_current, ce = init_p[4], gamma = gamma_eval, E0 = E0)
    }

    if(fixed_bis){ # target a fixed value of BIS directly
      ce_target <- Hinv(pars = pars_pd_current, bis = bis_target, E0 = E0, gamma = gamma_eval)
      kR_new <- TCI_basic(Ce = ce_target, pars = pars_pk_current, init = init_p)
      ivt_new <- list(list(begin = t0, end = t0 + delta, k_R = kR_new))
    } else{

      if(is.null(target_update_tms)) t1 = t0 + delta

      # update beta coefficients
      if(r %in% target_update_ix){
        print("Updating target sigmoid function")

        # evaluate infusion schedule up until the next update time
        t1 <- next_update_tm[which(target_update_ix == r)]

        kR_quad <- c(ivt_eval, sig_inf(beta = beta_init, theta = c(pars_pk_current, pars_pd_current), ini = init_p, t0 = t0, gamma = gamma_eval, E0 = E0, tmx = t1))
        # }

        if(qmc){
          samples <- quad_pars_qmc(lpars = lpars_prior$mu, sig = lpars_prior$sig, npoints = 100)
          quad <- list(points = samples, weights = rep(1,nrow(samples)))
        } else{
          quad <- quad_pars(lpars_prior$mu, sig = lpars_prior$sig, kR = kR_quad, gamma = gamma_eval, E0 = E0)
        }

        if(!is.null(eps)){
          # optimize version - one parameter constrained to equal bis target + eps at final point
          # nmin_target_prop = t0 + nmin_prop*(nmin-t0)
          nmin_target_prop = nmin_prop*nmin
          opt_res <- optimize(f = sigmoid_update_objfn, interval = c(0,3),
                              ltheta_hat = lpars_prior$mu,
                              Theta_samples = exp(quad$points), weights = quad$weights,
                              sig = lpars_prior$sig,
                              t0 = t0, ivt0 = ivt_eval, init_p = init_p, t1 = t1,
                              E0 = E0, gamma = gamma_eval,
                              alpha = alpha, bis_target = bis_target, p = p,
                              nmin = nmin,
                              nmin_target = nmin_target_prop,
                              eps = eps)
          lgamma_opt = opt_res$minimum
          lbis_t50_opt = log(nmin_target_prop) + 1/exp(lgamma_opt) * (log(1) - log(E0 - bis_target + 1))
          beta_est <- exp(c(lbis_t50_opt, lgamma_opt))
        } else{
          optimx_res <- optimx::optimx(par = log(beta_est), fn = sig_ED, ltheta_hat = lpars_prior$mu,
                                       Theta_samples = exp(quad$points), weights = quad$weights,
                                       sig = lpars_prior$sig, t0 = t0, ivt0 = ivt_eval, init_p = init_p, E0 = E0, gamma = gamma_eval,
                                       alpha = alpha, bis_target = bis_target, nmin = nmin,
                                       method = "nlm", control = list(rel.tol = 1e-4))
          beta_est <- c(unname(exp(coef(optimx_res))))
        }

        print(beta_est)
      }

      ivt_new <- sig_inf(beta = beta_est, theta = exp(lpars_prior$mu), ini = init_p, t0 = t0, E0 = E0, gamma = gamma_eval,
                         tmx = t1,
                         bist = bis_target)
    }

    ivt_eval <- c(ivt_eval, ivt_new[1])
    dat_i <- sample_data(ivt_d = list(list(begin = 0, end = delta, k_R = ivt_eval[[length(ivt_eval)]]$k_R)),
                         pars_pk_tci = pars_pk_current, pars_pk0 = pars_pk0, pars_pd0 = pars_pd0, sigma_bis0 = sigma_bis0,
                         gamma1 = gamma1, gamma2 = gamma2, E0 = E0, delay = lag0, start_time = update_tms[r],
                         init_pred = init_p, init_true = init_0, nmin_sample = delta)
    dat_obs <- rbind(dat_obs, dat_i)

    if(update_pars_logical & update_tms[min(r+1, length(update_tms))] >= lag0/60){
      if(laplace_appx & update_tms[r] > laplace_update_tm){
        dat_eval = dat_obs[dat_obs$timedelay <=(update_tms[r] + delta) & dat_obs$timedelay > update_tms[r],]
        lpr_eval = lpars_prior
      } else{
        dat_eval = dat_obs[dat_obs$timedelay <= (update_tms[r] + delta),]
        lpr_eval = lhyper
      }

      lpars_post <- update_pars(lp = lpars_prior$mu,
                                ivt = ivt_eval,
                                dat = dat_eval,
                                lpr = lpr_eval,
                                gamma = gamma_eval,
                                E0 = E0)
    } else{
      lpars_post <- lpars_prior
    }

    # save values for comparison of appoximation methods
    # lpars_post_list[[r]] <- lpars_post

    if(plot_path){
      tseq <- seq(0, nmin, 1/60)
      plot(dat_obs$time, dat_obs$bis, xlim = c(0,nmin), ylim = c(0,100),
           col = ifelse(dat_obs$timedelay <=update_tms[min(r+1, length(update_tms))], rgb(red = 0, green = 0, blue = 1, alpha = 0.2), rgb(red = 0, green = 0, blue = 0, alpha = 0.2)),
           pch = 19,
           ylab = "bis", xlab = "min")
      lines(dat_obs$time, dat_obs$bis_t)
      lines(tseq, Emax(pars = beta_est, ce = tseq, E0 = E0, Emx = E0-bis_target), col = 2, lty = 2)
      abline(h = 50)

      # plot predicted bis
      solp <- pk_solution_3cpt_metab(pars = exp(lpars_prior$mu[1:9]), ivt = c(ivt_eval,ivt_new[-1]), init = c(0,0,0,0))
      lines(tseq, Emax1(pars = exp(lpars_post$mu[10]), ce = solp(tseq)[4,], gamma = gamma_eval, E0 = E0), col = 4, lty = 2)
      inf_tms <- sapply(ivt_eval, `[[`, "begin")
      inf_amt <- sapply(ivt_eval, `[[`, "k_R")
      for(j in 1:r){
        rug(x = inf_tms[j], quiet = T, ticksize = -1/8*(inf_amt[j]/max(inf_amt[1:r])),
            side = 3, lwd = 0.5, col = 4, outer = F, line = 0)
      }
    }
    lpars_prior <- lpars_post
  }
  out <- list(lpars_post = lpars_post, dat = dat_obs, ivt = ivt_eval)
  return(out)
}




minimax_induct_eleveld <- function(pars0, lhyper, target_update_tms, alpha = 0.05, p = NULL, delta = 1/6,
                                   nmin = 10, bis_target = 50, bis_lb = 40, plot_path = T, update_pars_logical = T,
                                   qmc = T, laplace_appx = T, seed = NULL, laplace_update_tm = 4, eps = 1,
                                   target_range = F){
  library(mvtnorm, quietly = T)
  library(numDeriv, quietly = T)
  library(nloptr, quietly = T)
  if(!is.null(seed)) set.seed(seed)
  pars_pk0 <- pars0$pk
  pars_pd0 <- pars0$pd
  sigma_bis0 <- pars0$sigma

  # assumed to not vary by individual
  lag0 = pars0$lag
  gamma1 = pars0$gamma
  gamma2 = pars0$gamma2
  E0 = pars0$E0

  # functions
  max_ce <- function(pars, I0, init, tmax_search = 20, grid_space = 1/60, B = NULL, E = NULL){
    if(is.null(B) & is.null(E)){
      # course without infusion - use current concentration
      B <- pk_solution_3cpt_metab_singleinf(pars = pars, ivt = list(begin = 0, end = 1/6, k_R = 0), init = init, ce_only = T)
      # course with infusion starting from 0 concentration
      E <- pk_solution_3cpt_metab_singleinf(pars = pars, ivt = list(begin = 0, end = 1/6, k_R = 1), init = c(0,0,0,0), ce_only = T)
    }
    tms <- seq(grid_space,tmax_search,grid_space)
    ceproj <- B(tms) + I0*E(tms)
    maxce <- max(ceproj)
    jpeak = tms[which(ceproj == maxce)]
    return(c(ce = maxce, jpeak = jpeak))
  }

  overshoot_objfn <- function(lkR, esamples, init, pars_fixed, alpha = 0.25, B_samples = NULL, E_samples = NULL){

    kR <- exp(lkR)
    if(is.null(B_samples) | is.null(E_samples)){
      B_samples <- apply(esamples, 1, pk_solution_3cpt_metab_singleinf, ivt = list(begin = 0, end = 1/6, k_R = 0), init = init, ce_only = T)
      E_samples <- apply(esamples, 1, pk_solution_3cpt_metab_singleinf, ivt = list(begin = 0, end = 1/6, k_R = 1), init = c(0,0,0,0), ce_only = T)
    }
    ce_samples <- t(sapply(1:nrow(esamples), function(i){
      max_ce(pars = esamples[i,], I0 = kR, init = init_d, B = B_samples[[i]], E = E_samples[[i]])
    }))

    min_bis_p <- unlist(sapply(1:nrow(esamples), function(i){
      unname(Emax1(pars = esamples[i,10], ce = ce_samples[i,1], gamma = pars_fixed["gamma"], E0 = pars_fixed["E0"]))
    }))

    prop_overshoot <- mean(min_bis_p < bis_lb)
    return(prop_overshoot - alpha)
  }

  update_tms <- seq(0, nmin-delta, delta)
  # target_update_ix <- which(update_tms %in% target_update_tms)
  # next_update_tm <- c(target_update_tms[-1],10)
  dat_obs <- as.data.frame(matrix(NA, ncol = 12, nrow = 0))
  lpars_prior <- lhyper
  # lpars_post_list <- matrix(NA, nrow = length(update_tms), ncol = length(lpars_prior$mu))

  for(r in 1:length(update_tms)){
    print(r)
    pars_pk_current <- exp(lpars_prior$mu[1:9])
    pars_pd_current <- exp(lpars_prior$mu[10])

    # Set initial values
    if(nrow(dat_obs) == 0) {
      init_0 <- init_p <- c(0,0,0,0) # predicted and true initial concentrations
      bis_p <- E0
      t0 = 0
      gamma_eval = gamma1
      ivt_eval <- NULL
    } else{
      # true concentrations
      init_0 <- as.numeric(dat_obs[nrow(dat_obs),c("c1_t","c2_t","c3_t","ce_t")])
      # current time of evaluation
      t0 <- ivt_eval[[length(ivt_eval)]]$end
      # predicted concentration given infusion schedule administered
      init_p <- c(t(pk_solution_3cpt_metab(pars = pars_pk_current, ivt = ivt_eval, init = c(0,0,0,0))(t0)))
      # gamma value isn't updated, but switches based on whether the concentration is greater than C50
      gamma_eval <- unname(ifelse(init_p[4] <= pars_pd_current, gamma1, gamma2))
      # currently predicted BIS
      bis_p <- Emax1(pars = pars_pd_current, ce = init_p[4], gamma = gamma_eval, E0 = E0)
    }

    samples <- quad_pars_qmc(lpars = lpars_prior$mu, sig = lpars_prior$sig, npoints = 100)
    quad <- list(points = samples, weights = rep(1,nrow(samples)))

    # calculate infusion necessary to reach 50 at point estimate
    ce_target <- Hinv(pars = pars_pd_current, bis = bis_target, E0 = E0, gamma = gamma_eval)
    kR_50 <- TCI_basic(Ce = ce_target, pars = pars_pk_current, init = init_p)

    # if within range, target bis = 50
    if(bis_p < 60 & bis_p > 40 & target_range) kR_alpha <- kR_50
    else{
      if(kR_50 > 0){ # min infusion is 0 --> don't calculate unless necessary.
        # calculate infusion that overshoots for alpha percent
        B_samples <- apply(exp(samples), 1, pk_solution_3cpt_metab_singleinf, ivt = list(begin = 0, end = 1/6, k_R = 0), init = init_p, ce_only = T)
        E_samples <- apply(exp(samples), 1, pk_solution_3cpt_metab_singleinf, ivt = list(begin = 0, end = 1/6, k_R = 1), init = c(0,0,0,0), ce_only = T)

        overshoot_lb <- overshoot_objfn(1, esamples = exp(samples), init = init_p,
                                        B_samples = B_samples, E_samples = E_samples,
                                        pars_fixed = c(E0 = E0, gamma = gamma_eval), alpha = alpha)
        if(overshoot_lb < 0){
          # f <- function(ld) -0.20489 + 0.039*ld
          # alpha_test <- f(-log(det(lpars_prior$sig)))
          # print(paste("alpha:",alpha_test))
          # print(paste("log|D|:",-log(det(lpars_prior$sig))))
          alpha_test <- alpha
          kR_alpha <- exp(uniroot(f = overshoot_objfn, interval = c(1,10), esamples = exp(samples), init = init_p,
                                  B_samples = B_samples, E_samples = E_samples, pars_fixed = c(E0 = E0, gamma = gamma_eval), alpha = alpha_test)$root)
        } else kR_alpha <- 0
      } else kR_alpha <- 0
    }


    ivt_new <- list(list(begin = t0, end = t0 + delta, k_R = min(kR_alpha, kR_50)))

    ivt_eval <- c(ivt_eval, ivt_new)
    dat_i <- sample_data(ivt_d = list(list(begin = 0, end = delta, k_R = ivt_eval[[length(ivt_eval)]]$k_R)),
                         pars_pk_tci = pars_pk_current, pars_pk0 = pars_pk0, pars_pd0 = pars_pd0, sigma_bis0 = sigma_bis0,
                         gamma1 = gamma1, gamma2 = gamma2, E0 = E0, delay = lag0, start_time = update_tms[r],
                         init_pred = init_p, init_true = init_0, nmin_sample = delta)
    dat_obs <- rbind(dat_obs, dat_i)

    if(update_pars_logical & update_tms[min(r+1, length(update_tms))] >= lag0/60){
      if(laplace_appx & update_tms[r] > laplace_update_tm){
        dat_eval = dat_obs[dat_obs$timedelay <=(update_tms[r] + delta) & dat_obs$timedelay > update_tms[r],]
        lpr_eval = lpars_prior
      } else{
        dat_eval = dat_obs[dat_obs$timedelay <= (update_tms[r] + delta),]
        lpr_eval = lhyper
      }

      lpars_post <- update_pars(lp = lpars_prior$mu,
                                ivt = ivt_eval,
                                dat = dat_eval,
                                lpr = lpr_eval,
                                gamma = gamma_eval,
                                E0 = E0)
    } else{
      lpars_post <- lpars_prior
    }

    if(plot_path){
      tseq <- seq(0, nmin, 1/60)
      plot(dat_obs$time, dat_obs$bis, xlim = c(0,nmin), ylim = c(0,100),
           col = ifelse(dat_obs$timedelay <=update_tms[min(r+1, length(update_tms))], rgb(red = 0, green = 0, blue = 1, alpha = 0.2), rgb(red = 0, green = 0, blue = 0, alpha = 0.2)),
           pch = 19,
           ylab = "bis", xlab = "min")
      lines(dat_obs$time, dat_obs$bis_t)
      abline(h = 50)

      # plot predicted bis
      solp <- pk_solution_3cpt_metab(pars = exp(lpars_prior$mu[1:9]), ivt = c(ivt_eval,ivt_new[-1]), init = c(0,0,0,0))
      lines(tseq, Emax1(pars = exp(lpars_post$mu[10]), ce = solp(tseq)[4,], gamma = gamma_eval, E0 = E0), col = 4, lty = 2)
      inf_tms <- sapply(ivt_eval, `[[`, "begin")
      inf_amt <- sapply(ivt_eval, `[[`, "k_R")
      for(j in 1:r){
        rug(x = inf_tms[j], quiet = T, ticksize = -1/8*(inf_amt[j]/max(inf_amt[1:r])),
            side = 3, lwd = 0.5, col = 4, outer = F, line = 0)
      }
    }
    lpars_prior <- lpars_post
  }
  out <- list(lpars_post = lpars_post, dat = dat_obs, ivt = ivt_eval)
  return(out)
}


plot_res <- function(res_list, lpars_fixed){
  nplots <- length(res_list)
  par(mar=c(1,1,1,1), mfrow = c(nplots/2,2))
  pars_fixed <- exp(lpars_fixed)

  plot_person <- function(res, xaxt_val){
    lpars <- res$lpars_post
    dat <- res$dat
    ivt <- res$ivt
    nmin <- max(dat$time)
    tms <- dat$time

    sol.pr <- pk_solution_3cpt_metab(pars = exp(lpars$mu[1:9]), ivt = ivt, init = c(0,0,0,0))
    con.pr <- sol.pr(tms)
    gamma_eval <- ifelse(con.pr[4,] <= exp(lpars$mu[10]), pars_fixed["gamma"], pars_fixed["gamma2"])

    bis.pr <- Emax1(pars = exp(lpars$mu[10]), ce = con.pr[4,], gamma = gamma_eval, E0 = pars_fixed["E0"])

    plot(dat$time, dat$bis_t, type = "l", xlim = c(0,nmin), ylim = c(0,100), xaxt = xaxt_val); abline(h = 50, lty=2)
    polygon(x = rep(par()$usr[c(1,2)], each = 2), y = c(40,60,60,40), col=rgb(.75,.75,.75, alpha = 0.1))
    lines(dat$time, bis.pr, col = 4)
  }

  for(i in 1:nplots){
    xaxt_val <- ifelse(i %in% c(nplots, nplots-1), "s", "n")
    plot_person(res_list[[i]], xaxt_val = xaxt_val)
    # title(sub = paste("Patient", i), cex = 0.8, adj= 1, line = -1)
  }

  par(mar=c(5.1,4.1,4.1,2.1))
  par(mfrow = c(1,1))
}


# extract summary statistics from simulation object
res_stats <- function(sim, lb = 40, ub = 60){
  tms <- sim[[1]]$dat$time
  bis_vals <- sapply(sim, function(x) x$dat$bis_t)
  first_pass <- apply(bis_vals, 2, function(x) tms[which(x<ub)[1]])
  overshoot <- apply(bis_vals, 2, function(x) tms[which(x<lb)[1]])
  max_overshoot <- apply(bis_vals, 2, min)
  max_overshoot[max_overshoot>lb] <- NA
  second_pass <- apply(bis_vals, 2, function(x) {
    tm_overshoot <- tms[which(x<lb)[1]]
    tms[intersect(which(x>lb), which(tms > tm_overshoot))][1]}
  )
  stable_entry <- ifelse(is.na(max_overshoot), first_pass, second_pass)
  return(cbind(first_pass, max_overshoot, second_pass, stable_entry))
}

rel_inf <- function(sim_ref, sim_test){
  inf_ref <- sapply(sim_ref, function(x) sum(sapply(x$ivt, `[[`, "k_R")))
  inf_test <- sapply(sim_test, function(x) sum(sapply(x$ivt, `[[`, "k_R")))
  return(inf_test/inf_ref)
}
