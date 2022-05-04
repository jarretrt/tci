# --------------------------------------------------------------------------------------------------------------------------------
# - PK-PD model methods ----------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------

## Methods for pkmod objects
#' Create an object with class "pkmod"
#' @param pars_pk Vector or matrix of named PK parameters. If not specified, the pkmod function will be
#' inferred from the parameter names. Print `list_parnms()` for acceptable parameter names.
#' @param init Vector of initial concentrations. Will default to values of zero in all compartments if not specified.
#' @param pkfn PK model function. Functions provided in `tci` include `pkmod1cpt`, `pkmod2cpt`, `pkmod3cpt`, and `pkmod3cptm`.
#' User-defined functions should be specified here.
#' @param pars_pd PD model parameters if a PD model is specified
#' @param pdfn PD model function
#' @param pdinv Inverse PD model function for use in TCI algorithms
#' @param pcmpt Index of plasma compartment. Defaults to first compartment if not specified.
#' @param ecmpt Index of effect-site compartment if a PD model is specified. Will default to last compartment if unspecified.
#' @return A list with class pkmod for which print, plot, predict, and simulate methods exist.
#' @param sigma_add Standard deviation of additive residual error
#' @param sigma_mult Standard deviation of multiplicative residual error
#' @param log_response Logical value indicating if the response should be logged prior to adding error.
#' @param Omega Optional matrix of random effect parameters. Column names should correspond to names of
#' pars_pk, pars_pd, and sigma_add or sigma_mult.
#' @rdname init_pkmod
#' @examples
#' # create a pkmod object for a one compartment model with plasma targeting
#' init_pkmod()
#' @export
init_pkmod <- function(pars_pk = NULL, init = NULL, pkfn = NULL, pars_pd = NULL, pdfn = NULL,
                       pdinv = NULL, pcmpt = NULL, ecmpt = NULL, sigma_add = NULL, sigma_mult = NULL,
                       log_response = NULL, Omega = NULL){

  pkmod_obj <- list(pars_pk = pars_pk,
                    init = init,
                    pkfn = pkfn,
                    pars_pd = pars_pd,
                    pdfn = pdfn,
                    pdinv = pdinv,
                    pcmpt = pcmpt,
                    ecmpt = ecmpt,
                    sigma_add = sigma_add,
                    sigma_mult = sigma_mult,
                    log_response = log_response,
                    Omega = Omega)

  class(pkmod_obj) <- "pkmod"
  return(pkmod_obj)
}


#' pkmod validation checks
#'
#' Function to provide validation checks for a pkmod object
#' @param x Object of class "pkmod"
#' @examples
#' validate_pkmod(init_pkmod(pars_pk = c(CL = 10, V1 = 10,Q2 = 4,V2=20)))
#' @return Returns a list with class "pkmod" if validation checks are passed. Returns an error if not.
#' @export
validate_pkmod <- function(x){

  ## match PK model
  if(is.null(x$pkfn)){

    if(is.null(x$pars_pk))
      stop("pars_pk must be provided")

    if(is.null(names(x$pars_pk)))
      stop("pars_pk must be named. See `pars_pk` argument of ?pkmod for more details.")

    x$pkfn <- infer_pkfn(names(x$pars_pk))
  } else{
    if(!all(c("tm","kR","pars","init") %in% names(formals(x$pkfn))))
      stop('pkfn must contain arguments ("tm","kR","pars","init").')
  }

  x$ncmpt <- length(eval(formals(x$pkfn)$init))

  ## store initial values
  if(is.null(x$init)) x$init <- rep(0,x$ncmpt)
  if(is.null(x$pcmpt)) x$pcmpt = 1
  if(is.null(x$ecmpt)) x$ecmpt = length(x$init)

  ## check Omega matrix if provided
  if(!is.null(x$Omega)){
    if(nrow(x$Omega)!=ncol(x$Omega)) stop("Omega must be a square matrix")
    if(!all(colnames(x$Omega) %in% c(names(x$pars_pk), names(x$pars_pd), "sigma_add","sigma_mult")))
      stop("Column names of Omega must match names of pars_pk, pars_pd, sigma_add, sigma_mult")
  }

  return(x)
}
#' Create a pkmod object
#'
#' User function to create pkmod objects with validation checks. Multiple rows are permitted for
#' arguments `pars_pk`, `pars_pd`, `init`
#'
#' @param pars_pk Vector or matrix of named PK parameters. If not specified, the pkmod function will be
#' inferred from the parameter names. Print `list_parnms()` for acceptable parameter names.
#' @param init Vector of initial concentrations. Will default to values of zero in all compartments if not specified.
#' @param pkfn PK model function. Functions provided in `tci` include `pkmod1cpt`, `pkmod2cpt`, `pkmod3cpt`, and `pkmod3cptm`.
#' User-defined functions should be specified here.
#' @param pars_pd PD model parameters if a PD model is specified
#' @param pdfn PD model function
#' @param pdinv Inverse PD model function for use in TCI algorithms
#' @param pcmpt Index of plasma compartment. Defaults to first compartment if not specified.
#' @param ecmpt Index of effect-site compartment if a PD model is specified. Will default to last compartment if unspecified.
#' @param sigma_add Standard deviation of additive residual error
#' @param sigma_mult Standard deviation of multiplicative residual error
#' @param log_response Logical value indicating if the response should be logged prior to adding error.
#' @param Omega Optional matrix of random effect parameters. Column names should correspond to names of
#' pars_pk, pars_pd, and sigma_add or sigma_mult.
#' @rdname pkmod
#' @return Returns a list with class "pkmod" if validation checks are passed. Returns an error if not.
#' @examples
#' # 1-compartment model
#' pkmod(pars_pk = c(CL = 10, V1 = 10))
#' # 2-compartment model
#' pkmod(pars_pk = c(CL = 10, V1 = 10, Q = 3, v2 = 20))
#' # 2-compartment model with random effects matrix
#' Omega <- matrix(diag(c(0.3,0.2,0,0.4)), 4,4, dimnames = list(NULL,c("CL","V1","Q","v2")))
#' pkmod(pars_pk = c(CL = 10, V1 = 10, Q = 3, v2 = 20), Omega = Omega)
#' @export
pkmod <- function(pars_pk = NULL, init = NULL, pkfn = NULL, pars_pd = NULL, pdfn = NULL,
                  pdinv = NULL, pcmpt = NULL, ecmpt = NULL, sigma_add = 0, sigma_mult = 0,
                  log_response = FALSE, Omega = NULL){

  new_mod <- init_pkmod(pars_pk = pars_pk,
                        init = init,
                        pkfn = pkfn,
                        pars_pd = pars_pd,
                        pdfn = pdfn,
                        pdinv = pdinv,
                        ecmpt = ecmpt,
                        sigma_add = sigma_add,
                        sigma_mult = sigma_mult,
                        log_response = log_response,
                        Omega = Omega)
  validate_pkmod(new_mod)
}


#' Print pkmod
#'
#' Print method for pkmod objects
#' @param x Object with class "pkmod" created by `pkmod()`
#' @param ... Arguments passed to `update.pkmod()`
#' @param digits Number of significant digits to print
#'
#' @examples
#' # 1-parameter model
#' pkmod(pars_pk = c(CL = 10, V1 = 10))
#' @return Prints description of pkmod
#' @export
print.pkmod <- function(x, ..., digits = 3){

  x <- update(x,...)
  cat("--- PK model -------------------------------------------", "\n")
  cat(paste0(x$ncmpt,"-compartment PK model"),"\n")
  if(!is.null(x$pars_pk))
    cat(paste0("PK parameters: ",
               paste(names(x$pars_pk), "=", signif(x$pars_pk,digits),
                     collapse = ", ")),"\n")
  if(!is.null(x$init))
    cat(paste0("Initial concentrations: ",paste0("(",paste(signif(x$init,digits),
                                                           collapse = ","),")")),"\n")

  cat(paste0("Plasma compartment: ",x$pcmpt),"\n")
  if(!is.null(x$ecmpt)) cat(paste0("Effect compartment: ",x$ecmpt),"\n")
  if(!is.null(x$pars_pd) & length(x$pars_pd) > 0){
    cat("--- PD model -------------------------------------------", "\n")
    cat(paste0("PD parameters: ",paste(names(x$pars_pd), "=", signif(x$pars_pd,digits),
                                       collapse = ", ")),"\n")
  }
  if(!is.null(x$sigma_add)|!is.null(x$sigma_mult)|!is.null(x$log_response)){
    cat("--- Simulation -----------------------------------------", "\n")
  }

  if(!is.null(x$sigma_add)) cat(paste0("Additive error SD: ",signif(x$sigma_add,digits)),"\n")
  if(!is.null(x$sigma_mult)) cat(paste0("Multiplicative error SD: ",signif(x$sigma_mult,digits)),"\n")
  if(!is.null(x$log_response)) cat(paste0("Logged response: ",x$log_response),"\n")
  if(!is.null(x$Omega)){
    cat("Random effects","\n")
    print(signif(x$Omega,digits))
  }
}


#' Update method for pkmod
#'
#' Update parameters or initial values of a pkmod object
#' @param object Object with class pkmod
#' @param ... Updated values passed to object.
#' @return Returns the pkmod object with elements replaced
#' @examples
#' # initial pkmod object
#' (my_mod <- pkmod(pars_pk = c(CL = 10, V1 = 10)))
#' # update a subset of parameters and initial values
#' update(my_mod, pars_pk = c(CL = 20), init = 3, sigma_add = 1, log_response = TRUE)
#' @export
update.pkmod <- function(object, ...){

  update_elements <- list(...)
  if(length(update_elements) == 0) return(object)
  if(is.list(update_elements[[1]])) update_elements <- update_elements[[1]]

  # subset to elements with valid names
  update_elements <- update_elements[names(update_elements) %in% c(names(formals(pkmod)))]

  for(nm in names(update_elements)){
    # partial update - if matrix, assumes that the whole matrix will be replaced
    if(nm %in% c("pars_pk","pars_pd") & is.vector(object[[nm]])){
      object[[nm]][names(update_elements[[nm]])] <- update_elements[[nm]]
    } else{
      object[[nm]] <- update_elements[[nm]]
    }
  }

  return(validate_pkmod(object))
}

# update.pkmod <- function(object, ...){
#  update_elements <- list(...)
#  if(length(update_elements) == 0) return(object)
#  if(is.list(update_elements[[1]])) update_elements <- update_elements[[1]]
#  if(!all(names(update_elements) %in% names(object))){
#    nms_mismatch <- names(update_elements)[!names(update_elements) %in% names(object)]
#    warning(paste(paste(nms_mismatch, collapse = ", "), "are not valid names."))
#  }
#  # subset to elements with valid names
#  update_elements <- update_elements[names(update_elements) %in% c(names(object))]
#  for(nm in names(update_elements)){
#    # partial update - if matrix, assumes that the whole matrix will be replaced
#    if(nm %in% c("pars_pk","pars_pd")){
#      # if(is.vector(update_elements[[nm]])){
#      #   update_elements[[nm]] <- matrix(update_elements[[nm]], nrow = 1,
#      #                                   dimnames = list(NULL, names(update_elements[[nm]])))
#      # }
#      # if(dim(update_elements[[nm]])[1] == dim(object[[nm]])[1]){
#        object[[nm]][,names(update_elements[[nm]])] <- update_elements[[nm]]
#      # } else{
#      #   object[[nm]] <- update_elements[[nm]]
#      # }
#    } else{
#      object[[nm]] <- update_elements[[nm]]
#    }
#  }
#  return(validate_pkmod(object))
# }



#' Predict method for pkmod objects
#'
#' Predict concentrations from a pkmod object - can be a user defined function
#'
#' @param object An object with class pkmod.
#' @param inf A matrix with columns "begin","end","infrt" indicating when infusions should
#' be administered. Can be created by `create_inf` or `apply_tci`.
#' @param tms Times at which to calculate predicted concentrations.
#' @param ... List or vector of values to be passed on to update.pkmod.
#' @return Matrix of predicted concentrations associated with a pkmod object and
#' and infusion schedule.
#' @examples
#' # dosing schedule
#' dose <- create_inf(times = c(0,0.5,4,4.5,10), c(100,0,80,0,0))
#' # pkmod object
#' my_mod <- pkmod(pars_pk = c(CL = 15, V1 = 10, Q2 = 10, V2 = 20))
#' # predict at specific times
#' predict(my_mod, inf = dose, tms = c(1.5,2.5,3))
#' # predict with an Emax PD function
#' my_mod_pd <- pkmod(pars_pk = c(v1 = 8.995, v2 = 17.297, v3 = 120.963, cl = 1.382,
#' q2 = 0.919, q3 = 0.609, ke0 = 1.289),
#' pars_pd = c(c50 = 2.8, gamma = 1.47, gamma2 = 1.89, e0 = 93, emx = 93),
#' pdfn = emax, pdinv = emax_inv, ecmpt = 4)
#' predict(my_mod_pd, inf = dose, tms = c(1.5,2.5,3))
#' # predict with a subset of new PK-PD parameters
#' predict(my_mod_pd, inf = dose, tms = c(1.5,2.5,3), pars_pk = c(ke0 = 0.8),
#' pars_pd = c(c50 = 2, e0 = 100))
#' @export
predict.pkmod <- function(object, inf, tms, ...){

  # update pkmod
  object <- update(object, ...)

  # check columns of inf
  if(!all(c("infrt","begin","end") %in% colnames(inf)))
    stop("inf must include 'infrt','begin','end' as column names")

  # apply Rcpp implementation to selected PK models
  if(any(identical(object$pkfn, pkmod1cpt),
         identical(object$pkfn, pkmod2cpt),
         identical(object$pkfn, pkmod3cpt),
         identical(object$pkfn, pkmod3cptm))){

    starttm <- inf[1,"begin"]
    begin <- inf[,"begin"] - starttm
    end <- inf[,"end"] - starttm
    infs <- inf[,"infrt"]
    tms_eval <- tms - starttm
    pars <- format_pars(object$pars_pk, ncmpt = object$ncmpt)

    if(object$ncmpt == 1){
      pred <- t.default(pksol1cpt(tms_eval, pars, begin, end, infs, object$init))
    }
    if(object$ncmpt == 2){
      pred <- pksol2cpt(tms_eval, pars, begin, end, infs, object$init)
    }
    if(object$ncmpt == 3){
      pred <- pksol3cpt(tms_eval, pars, begin, end, infs, object$init)
    }
    if(object$ncmpt == 4){
      pred <- pksol3cptm(tms_eval, pars, begin, end, infs, object$init)
    }

  } else{

    # round times to avoid mismatches in floating-point precision
    tm_digits = 7

    # Times to evaluate concentrations at.
    # round times - this is needed to prevent errors associated with rounding numeric values
    tms <- round(tms, tm_digits)

    if(inf[nrow(inf),"infrt"] !=0 & inf[nrow(inf),"end"]<max(tms)){
      inf <- rbind(inf, c(begin = inf[nrow(inf),"end"], end = max(tms), infrt = 0))
    }

    # if times are provided, predict at those times plus boundaries for initial values
    bnd <- sort(unique(
      round(as.numeric(unlist(inf[,c("begin","end")])),tm_digits)
    ))
    tms_all <- sort(unique(
      round(c(bnd,tms),tm_digits)
    ))
    tms_eval <- split(tms_all, findInterval(tms_all, bnd,
                                            rightmost.closed = TRUE,
                                            left.open = TRUE))
    init <- vector("list", nrow(inf)+1)

    # Pass on initial concentrations to first element of init. Use values if specified, else defaults.
    init[[1]] <- object$init

    # get indexes of times and initialize matrix for predictions
    tm_ix <- lapply(tms_eval, function(x) match(x,tms_all))
    pred <- matrix(NA, nrow = object$ncmpt, ncol = length(tms_all))

    # Predict concentrations and store initial values.
    for(i in 1:nrow(inf)){
      # start predictions from tm=0
      pred[,tm_ix[[i]]] <- object$pkfn(tm = tms_eval[[i]] - inf[i,"begin"],
                                       kR = inf[i,"infrt"],
                                       pars = object$pars_pk,
                                       init = init[[i]])
      init[[i+1]] <- pred[,length(unlist(tm_ix[1:i]))]
    }

    # Replace any negative values
    pred[pred<0] <- 0

    # reduce to predictions only at tms
    pred <- pred[,unique(unlist(tms_eval)) %in% tms]
  }

  pred <- t.default(pred)
  colnames(pred) <- paste0("c",1:object$ncmpt)

  # predict PD component if present
  if(!is.null(object$pdfn)){
    pred <- cbind(pred, pdresp=object$pdfn(ce = pred[,object$ecmpt],
                                           pars = object$pars_pd))
  }

  return(pred)
}


#' Simulate method for pkmod objects
#'
#' Simulate observations from a pkmod object using a set of infusion rates or TCI targets.
#'
#' @param object An object with class "pkmod" generated by `pkmod`
#' @param nsim Number of observations to simulate at each time point. Defaults to 1.
#' @param seed An integer used to initialize the random number generator.
#' @param ... Arguments passed to `update.pkmod`.
#' @param inf A matrix of infusion rates with columns 'begin', 'end', and 'infrt'. This
#' can be created manually, by `create_inf`, or by `apply_tci`.
#' @param tms Times at which to simulate observations.
#' @param obs_cmpt Integer value indicating compartment in which observations are taken.
#' Overridden if a PD model is included.
#' @param resp_bounds Optional vector of two values indicating minimum and maximum values possible for the response.
#' @examples
#' # simulate data from a 2 compartment model with multiplicative error
#' dose <- create_inf(times = c(0,0.5,4,4.5,10), c(100,0,80,0,0))
#' my_mod <- pkmod(pars_pk = c(CL = 10, V1 = 10, Q2 = 4, V2 = 30))
#' tci_targets = cbind(value = c(2,3,4,4), time = c(0,2,3,10))
#' inf <- apply_tci(tci_targets, my_mod, "plasma")
#' simulate(my_mod, nsim = 3, inf = inf, tms = c(1,2,4,6,10), sigma_mult = 0.2)
#' # simulate with PD model
#' my_mod_pd <- pkmod(pars_pk = c(v1 = 8.995, v2 = 17.297, v3 = 120.963, cl = 1.382,
#' q2 = 0.919, q3 = 0.609, ke0 = 1.289),
#' pars_pd = c(c50 = 2.8, gamma = 1.47, gamma2 = 1.89, e0 = 93, emx = 93),
#' pdfn = emax, pdinv = emax_inv, ecmpt = 4)
#' simulate(my_mod_pd, inf = dose, tms = c(1.5,2.5,3))
#' @import stats
#' @export
simulate.pkmod <- function(object, nsim = 1, seed = NULL, ..., inf, tms,
                           obs_cmpt = 1, resp_bounds = NULL){

  if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if(is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }

  object <- update(object, ...)

  # predict at each set of true parameter values
  con <- predict(object = object, inf = inf, tms = tms)

  # PD response if applicable
  if(!is.null(object$pdfn)){
    resp <- con[,object$ncmpt+1]
  } else{
    resp <- con[,paste0("c",obs_cmpt)]
  }

  # add error term
  eadd  <- matrix(rnorm(nsim*length(resp),0,object$sigma_add),ncol = nsim)
  emult <- matrix(rnorm(nsim*length(resp),0,object$sigma_mult),ncol = nsim)

  if(object$log_response)
    obs <- exp(log(resp)*(1+emult) + eadd)
  else
    obs <- resp*(1+emult) + eadd

  if(!is.null(resp_bounds)){
    lb <- min(resp_bounds, na.rm = TRUE)
    ub <- max(resp_bounds, na.rm = TRUE)
    obs[obs<lb] <- lb
    obs[obs>ub] <- ub
  }

  return(obs)
}




## -- Plotting methods ---------------------------------------------------------

#' Plot method for 'pkmod'
#'
#' Will show predicted concentrations in compartments associated with an infusion schedule.
#' @param x An object with class pkmod.
#' @param ... Arguments passed on to predict.pkmod. Could include changing initial concentrations or parameter values through pkmod_update. See ?predict.pkmod
#' @param inf An infusion schedule object with columns "begin","end","infrt".
#' @param endtm Final time to evaluate predictions
#' @param title Plot title
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param plot_pk Logical. Should PK concentrations be displayed in addition to PD, when PD is present.
#' @param ylab_resp y-axis label for PD component, if present.
#' @return ggplot object displaying predicted concentrations for a pkmod object.
#' @rdname plot
#' @examples
#' # dosing schedule
#' dose <- create_inf(times = c(0,0.5,4,4.5,10), c(100,0,80,0,0))
#' # pkmod object
#' my_mod <- pkmod(pars_pk = c(CL = 15, V1 = 10, Q2 = 10, V2 = 20))
#' # plot predicted concentrations
#' plot(my_mod, inf = dose, ylab = "Concentration (mg/L)", xlab = "Minutes")
#' # plot with PD component
#' my_mod_pd <- pkmod(pars_pk = c(v1 = 8.995, v2 = 17.297, v3 = 120.963, cl = 1.382,
#' q2 = 0.919, q3 = 0.609, ke0 = 1.289),
#' pars_pd = c(c50 = 2.8, gamma = 1.47, gamma2 = 1.89, e0 = 93, emx = 93),
#' pdfn = emax, pdinv = emax_inv, ecmpt = 4)
#' plot(my_mod_pd, inf = dose, ylab = "Concentration (mg/L)", ylab_resp = "Bispectral Index",
#' xlab = "Minutes")
#' # plot with alternate PK and PD parameters
#' plot(my_mod_pd, inf = dose, pars = c(cl = 10), pars_pd = c(c50 = 2, emx = 100),
#' ylab = "Concentration (mg/L)", ylab_resp = "Bispectral Index", xlab = "Minutes")
#' @export
plot.pkmod <- function(x, ..., inf, endtm = NULL, title = NULL, xlab = "Time",
                       ylab = "Concentration", plot_pk = TRUE,
                       ylab_resp = "Response"){

  pdresp <- id <- value <- variable <- NULL

  # set of times to predict at
  begintm <- min(inf[,"begin"])
  if(is.null(endtm)) endtm <- max(inf[,"end"])
  tms <- seq(begintm, endtm, length.out = 1000)

  # predict concentrations
  con <- data.frame(predict(x, inf = inf, tms = tms, return_init = FALSE, ...))
  colnames(con) <- gsub("^c", "Compt. ", colnames(con))
  con[,"time"] <- tms

  if(!("id" %in% colnames(con))) con <- cbind(id = rep(1,nrow(con)), con)
  ids <- unique(con$id)

  if(is.null(x$pdfn)){
    tmp <- reshape::melt(con, id = c("id","time"))
    p <- ggplot2::ggplot(tmp,
                         ggplot2::aes(x = time,
                                      y = value,
                                      linetype = variable,
                                      color = variable))
    for(i in ids) p <- p + ggplot2::geom_line(data = tmp[tmp$id ==i,])
    p <- p +
      ggplot2::labs(y = ylab,
                    x = xlab,
                    color = "Compartment",
                    linetype = "Compartment",
                    title = title)
    p
  } else{
    # plot PK model
    if(plot_pk){
      tmp <- reshape::melt(con[,-which(colnames(con)=="pdresp")], id = c("id","time"))
      p1 <- ggplot2::ggplot(tmp,
                            ggplot2::aes(x = time,
                                         y = value,
                                         linetype = variable,
                                         color = variable))

      for(i in ids) p1 <- p1 + ggplot2::geom_line(data = tmp[tmp$id ==i,])
      p1 <- p1 +
        ggplot2::labs(y = ylab,
                      x = xlab,
                      color = "Compartment",
                      linetype = "Compartment") +
        ggplot2::theme(legend.position="bottom")
    } else{
      p1 <- NULL
    }

    p2 <- ggplot2::ggplot(con, ggplot2::aes(x = time, y = pdresp, group = id)) +
      ggplot2::geom_line() +
      ggplot2::labs(y = ylab_resp, x = xlab) +
      ggplot2::lims(y = c(0,100))

    if(!is.null(p1)){
      gridExtra::grid.arrange(p2, p1, nrow = 2, top = title)
    } else{
      gridExtra::grid.arrange(p2, top = title)
    }
  }
}


# #' Plotting method for tciinf objects. Will work for outputs from "iterate_tci_grid" or "tci_pd".
# #' Note: the tci targets function is left-continuous. Consequently, it plots
# #' final target value plotted will be a repetition of the preceding value
# #'
# #' @param x Object with class "tciinf" created by `tciinf`
# #' @param ... Arguments passed to predict.tciinf.
# #' @param targets A matrix or data frame with columns 'value' and 'time'. Times
# #' indicate when the TCI algorithm should begin to reach each target. If the TCI
# #' algorithm should begin immediately, either the first time value must be 0 or
# #' the `inittm` parameter must be set.
# #' @param cmpts Which compartment concentrations should be plotted. Defaults to all compartments.
# #' @param display Logical. Should plots be printed or returned as an arrangeGrob object?
# #' @param title Title of plot.
# #' @param xlab x-axis label
# #' @param ylab_con y-axis label for concentration-time plot
# #' @param ylab_resp y-axis label for response-time plot
# #' @return gtable object using gridExtra::arrangeGrob
# #' @examples
# #' # plasma-targeting
# #' my_mod <- pkmod(pars = c(CL = 10, V1 = 10), pkmod_fun = pkmod1cpt, init = 0)
# #' my_tci <- tciinf(tci_alg = tci_plasma, pkmod = my_mod, dtm = 1/6)
# #' plot(my_tci, targets = cbind(value = c(2,3,4,4), time = c(0,2,3,10)))
# #' # effect-site targeting TCI algorithm with PD model
# #' my_mod_pd <- pkmod(pars = c(v1 = 8.995, v2 = 17.297, v3 = 120.963, cl = 1.382, q2 = 0.919, q3 = 0.609, ke0 = 1.289), pkmod = pkmod3cptm, pars_pd = c(c50 = 2.8, gamma = 1.47, gamma2 = 1.89, e0 = 93, emx = 93, sigma = 8.03, bis_delay = 28.263), pdmod = emax, pdinv = emax_inv, ecmpt = 4, init = c(0,0,0,0))
# #' my_tci_pd <- tciinf(tci_alg = tci_comb, pkmod = my_mod_pd, dtm = 1/6)
# #' plot(my_tci_pd, targets = cbind(value = c(70,60,50,50), time = c(0,2,3,10)))
# #' # effect-site targeting for a matrix of PK-PD parameter values
# #' evd_pkpd_pars <- eleveld_poppk(merge(eleveld_pk, eleveld_pd))[c(1,20,50,70,100),]
# #' evd_pkmod <- pkmod(pars = evd_pkpd_pars[,1:7], pars_pd = evd_pkpd_pars[,8:11], pkmod_fun = pkmod3cptm, ecmpt = 4, pdmod_fun = emax_eleveld, pdinv_fun = emax_inv_eleveld)
# #' my_tci_pd_matrix <- tciinf(tci_alg = tci_comb, pkmod = evd_pkmod, dtm = 1/6)
# #' plot(my_tci_pd_matrix, targets = cbind(value = c(70,60,50,50), time = c(0,2,3,10)))
# #' @rdname plot
# #' @export
# plot.tciinf <- function(x, ..., targets, title = NULL, cmpts = NULL, display = TRUE, xlab = "Time",
#                         ylab_con = "Concentration", ylab_resp = "Response", ylab_inf = "Infusion Rate"){
#
#   begin <- value <- variable <- NULL
#   tciinf <- as.data.frame(predict(x, targets = targets, ...))
#
#   if(!is.null(cmpts)){
#     if(length(setdiff(cmpts,1:x$pkmod$ncmpt))>0) stop("Compartment number specified that doesn't exist")
#   } else{
#     cmpts <- 1:x$pkmod$ncmpt
#   }
#   ncpt <- length(cmpts)
#
#   if("id" %in% names(tciinf)){
#     last_rows <- as.matrix(tciinf[tciinf$begin == max(tciinf$begin),])
#     tmp <- matrix(NA,nrow = length(unique(tciinf$id)), ncol = ncol(tciinf))
#     colnames(tmp) <- colnames(tciinf)
#     tmp[,grep("id|start|begin|Ct|pdt", colnames(tmp))] <- last_rows[,grep("id|end|Ct|pdt", colnames(last_rows))]
#     tciinf <- as.data.frame(rbind(tciinf, tmp))
#   } else{
#     # move end of last row to start column for plotting
#     tciinf <- rbind(tciinf, NA)
#     tciinf[nrow(tciinf),grep("start|begin|Ct|pdt", names(tciinf))] <-
#       tciinf[nrow(tciinf)-1,grep("end|Ct|pdt", names(tciinf))]
#     tciinf$id <- rep(1,nrow(tciinf))
#   }
#
#   ids <- unique(tciinf$id)
#   # scale alpha for readability
#   alph <- emax(length(ids), pars = c(e0 = 1, emx = 0.8, c50 = 10, gamma = 1))
#
#   # data frame for PK plot
#   tciinfm <- reshape::melt(tciinf[,c("id","begin","Ct", paste0("c",cmpts,"_start"))],
#                            id.vars = c("id","begin"))
#   tciinfm$variable <- gsub("_start","",tciinfm$variable)
#   tciinfm$variable <- gsub("c","Compt. ",tciinfm$variable)
#   tciinfm$variable <- gsub("Ct","Target",tciinfm$variable)
#   tciinfm$variable <- factor(tciinfm$variable,
#                              levels = c("Target",grep("Compt. ",unique(tciinfm$variable),
#                                                       value = TRUE)))
#
#   ppk <- ggplot2::ggplot(tciinfm) +
#     ggplot2::geom_step(data = tciinfm[tciinfm$variable == "Target",],
#                        ggplot2::aes(x = begin, y = value, color = variable, linetype = variable, group = id),
#                        size = 1, alpha = alph)
#
#   for(i in ids){
#     ppk <- ppk + ggplot2::geom_line(data = tciinfm[tciinfm$variable != "Target" & tciinfm$id == i,],
#                                     ggplot2::aes(x = begin, y = value, linetype = variable, color = variable),
#                                     size = 1, alpha = alph)
#   }
#   ppk <- ppk + ggplot2::scale_linetype_manual("", values = c(cmpts,1),
#                                               labels = c(paste0("Compt. ",cmpts),"Target")) +
#     ggplot2::scale_color_manual("", values = c(unname(pal[2:(ncpt+1)]),"black"),
#                                 labels = c(paste0("Compt. ",cmpts),"Target")) +
#     ggplot2::labs(x = xlab, y = ylab_con, color = "", linetype = "") +
#     ggplot2::theme(legend.position="bottom")
#
#   if("pdresp_start" %in% names(tciinf)){
#     tciinfm2 <- reshape::melt(tciinf[,c("id","begin","pdt","pdresp_start")],id.vars = c("id","begin"))
#
#     ppd <- ggplot2::ggplot(tciinfm2,
#                            ggplot2::aes(x = begin, y = value)) +
#       ggplot2::geom_step(data = tciinfm2[tciinfm2$variable == "pdt",],
#                          ggplot2::aes(color = "col1", linetype = "solid", group = id), size = 1, alpha = alph)
#     for(i in ids){
#       ppd <- ppd + ggplot2::geom_line(data = tciinfm2[tciinfm2$variable != "pdt" & tciinfm2$id == i,],
#                                       ggplot2::aes(color = "col2", linetype = "solid"), size = 1, alpha = alph)
#     }
#     ppd <- ppd +
#       ggplot2::scale_linetype_manual("", values = c("solid","solid"), labels = c("PD Target","PD Response")) +
#       ggplot2::scale_color_manual("", values = unname(pal[c(1,5)]), labels = c("PD Target","PD Response")) +
#       ggplot2::ylim(c(0,100)) +
#       ggplot2::labs(x = xlab, y = ylab_resp) +
#       ggplot2::theme(legend.position="bottom")
#
#   }
#
#   pinf <- ggplot2::ggplot(tciinf) +
#     ggplot2::geom_step(ggplot2::aes(x = begin, y = infrt, group = id), size = 1,
#                        color = pal[6], alpha = alph) +
#     ggplot2::labs(x = xlab, y = ylab_inf)
#
#
#   if("pdresp_start" %in% names(tciinf)){
#     gb <- gridExtra::arrangeGrob(pinf, ppk, ppd, ncol = 1, top = title)
#   } else{
#     # gb <- ppk + ggplot2::ggtitle(title)
#     gb <- gridExtra::arrangeGrob(pinf, ppk, ncol = 1, top = title)
#   }
#
#   if(display){
#     plot(gb)
#   } else{
#     gb
#   }
# }

# #' Plotting method for tciinf objects. Will work for outputs from "iterate_tci_grid" or "tci_pd".
# #' Note: the tci targets function is left-continuous. Consequently, it plots
# #' final target value plotted will be a repetition of the preceding value
# #'
# #' @param x Object with class "tciinf" created by functions
# #' `iterate_tci_grid` or `tci_pd`
# #' @param ... \dots
# #' @param title Title of plot.
# #' @param display Logical. Should plots be printed or returned as an arrangeGrob object?
# #' @param xlab x-axis label
# #' @param ylab_con y-axis label for concentration-time plot
# #' @param ylab_resp y-axis label for response-time plot
# #' @return gtable object using gridExtra::arrangeGrob
# #' @rdname plot
# #' @export
# plot.tciinf <- function(x, ..., title = NULL, display = TRUE, xlab = "Time",
#                         ylab_con = "Concentration", ylab_resp = "Response"){
#
#   begin <- value <- variable <- NULL
#
#   tciinf <- as.data.frame(x)
#   # move end of last row to start column for plotting
#   tciinf <- rbind(tciinf, NA)
#   ncpt <- length(grep("c[0-9]_start",names(tciinf)))
#   tciinf[nrow(tciinf),grep("start|begin|Ct|pdt", names(tciinf))] <-
#     tciinf[nrow(tciinf)-1,grep("end|Ct|pdt", names(tciinf))]
#
#   # data frame for PK plot
#   tciinfm <- reshape::melt(tciinf[,c("begin","Ct", paste0("c",1:ncpt,"_start"))],
#                            id.vars = "begin")
#   tciinfm$variable <- gsub("_start","",tciinfm$variable)
#   tciinfm$variable <- gsub("c","Compt. ",tciinfm$variable)
#   tciinfm$variable <- gsub("Ct","Target",tciinfm$variable)
#   tciinfm$variable <- factor(tciinfm$variable,
#                              levels = c("Target",grep("Compt. ",unique(tciinfm$variable),
#                                                       value = TRUE)))
#
#   ppk <- ggplot2::ggplot(tciinfm) +
#     ggplot2::geom_step(data = tciinfm[tciinfm$variable == "Target",],
#                        ggplot2::aes(x = begin, y = value, color = variable, linetype = variable),
#               size = 1) +
#     ggplot2::geom_line(data = tciinfm[tciinfm$variable != "Target",],
#                        ggplot2::aes(x = begin, y = value, color = variable, linetype = variable),
#               size = 1) +
#     ggplot2::scale_linetype_manual("", values = c(1:ncpt,1),
#                           labels = c(paste0("Compt. ",1:ncpt),"Target")) +
#     ggplot2::scale_color_manual("", values = c(unname(pal[2:(ncpt+1)]),"black"),
#                        labels = c(paste0("Compt. ",1:ncpt),"Target")) +
#     ggplot2::labs(x = xlab, y = ylab_con, color = "", linetype = "") +
#     ggplot2::theme(legend.position="bottom")
#
#
#   if("pdresp_start" %in% names(tciinf)){
#     tciinfm2 <- reshape::melt(tciinf[,c("begin","pdt","pdresp_start")],id.vars = "begin")
#
#     ppd <- ggplot2::ggplot(tciinfm2,
#                            ggplot2::aes(x = begin, y = value)) +
#       ggplot2::geom_step(data = tciinfm2[tciinfm2$variable == "pdt",],
#                 ggplot2::aes(color = "col1", linetype = "solid"), size = 1) +
#       ggplot2::geom_line(data = tciinfm2[tciinfm2$variable != "pdt",],
#                          ggplot2::aes(color = "col2", linetype = "solid"), size = 1)+
#       ggplot2::scale_linetype_manual("", values = c("solid","solid"), labels = c("PD Target","PD Response")) +
#       ggplot2::scale_color_manual("", values = unname(pal[c(1,5)]), labels = c("PD Target","PD Response")) +
#       ggplot2::ylim(c(0,100)) +
#       ggplot2::labs(x = xlab, y = ylab_resp) +
#       ggplot2::theme(legend.position="bottom")
#   }
#
#   if("pdresp_start" %in% names(tciinf)){
#     gb <- gridExtra::arrangeGrob(ppd, ppk, nrow = 2, top = title)
#   } else{
#     gb <- ppk + ggplot2::ggtitle(title)
#   }
#
#   if(display){
#     plot(gb)
#   } else{
#     gb
#   }
# }




## -- Old functions ------------------------------------------------------------

# #' Plot method for datasim objects
# #'
# #' @param x An object with class datasim, created by function `gen_data`.
# #' @param pars_prior Named vector of prior PK or PK-PD parameters
# #' @param pars_post  Named vector of posterior PK or PK-PD parameters
# #' @param pk_ix Indicies of parameter vector(s) corresponding to PK parameters
# #' @param pd_ix Indicies of parameter vector(s) corresponding to PD parameters
# #' @param xlab x-axis label
# #' @param ylab_con y-axis label for concentration-time plot
# #' @param ylab_resp y-axis label for response-time plot
# #' @param ... \dots
# #' @return ggplot object displaying simulated data
# #' @rdname plot
# #' @importFrom ggplot2 aes
# #' @export
# plot.datasim <- function(x, ..., pars_prior = NULL, pars_post = NULL, pk_ix = NULL, pd_ix = NULL, xlab = "Time", ylab_con = "Concentration", ylab_resp = "Response"){
#
#   if((!is.null(pars_prior) | !is.null(pars_post)) & (is.null(pk_ix) | is.null(pd_ix)))
#     stop("pk_ix and pd_ix must be specified if pars_prior or pars_post are non-null")
#
#   # initialize for CRAN check -- will be created by other functions
#   c1 <- variable <- cobs <- pdp <- pdobs <- NULL
#
#   datasim <- x
#   datasim$sim <- as.data.frame(datasim$sim)
#   datasim$inf <- as.data.frame(datasim$inf)
#   r <- range(datasim$inf[,c("begin","end")])
#   tms <- seq(r[1], r[2], length.out = 1000)
#
#   # predict concentrations at true parameter values
#   cp <- data.frame(predict_pkmod(datasim$pkmod,
#                            inf = datasim$inf,
#                            tms = tms,
#                            pars = datasim$pars_pk0,
#                            init = datasim$init))
#
#   if(!is.null(datasim$pdmod)){
#     cp$pdp <- datasim$pdmod(ce = cp[,paste0("c",datasim$ecmpt)],
#                             pars = datasim$pars_pd0)
#   }
#   cp$variable <- as.factor("Observed")
#
#   # predict concentrations at prior parameters, if specified
#   if(!is.null(pars_prior)){
#     cp_prior <- data.frame(predict_pkmod(datasim$pkmod,
#                                    inf = datasim$inf,
#                                    tms = tms,
#                                    pars = pars_prior[pk_ix],
#                                    init = datasim$init))
#
#     if(!is.null(datasim$pdmod)){
#       cp_prior$pdp <- datasim$pdmod(ce = cp_prior[,paste0("c",datasim$ecmpt)],
#                                     pars = pars_prior[pd_ix])
#     }
#
#     cp_prior$variable <- as.factor("Prior")
#
#   } else{
#     cp_prior <- NULL
#   }
#
#   # predict concentrations at posterior parameters, if specified
#   if(!is.null(pars_post)){
#     cp_post <- data.frame(predict_pkmod(datasim$pkmod,
#                                   inf = datasim$inf,
#                                   tms = tms,
#                                   pars = pars_post[pk_ix],
#                                   init = datasim$init))
#
#     if(!is.null(datasim$pdmod)){
#       cp_post$pdp <- datasim$pdmod(ce = cp_post[,paste0("c",datasim$ecmpt)],
#                                    pars = pars_post[pd_ix])
#     }
#
#     cp_post$variable <- as.factor("Posterior")
#
#   } else{
#     cp_post <- NULL
#   }
#
#   df <- rbind(cp,cp_prior,cp_post)
#
#   # plot for pk simulations
#   if(is.null(datasim$pdmod)){
#     out <- ggplot2::ggplot(df,
#                            ggplot2::aes(x = time,
#                                         y = c1,
#                                         color = variable,
#                                         linetype = variable)) +
#       ggplot2::geom_line() +
#       ggplot2::geom_point(data = datasim$sim,
#                           ggplot2::aes(x = time, y = cobs),
#                           shape = 16, col = pal["navy"],
#                           inherit.aes = FALSE, alpha = 1) +
#       ggplot2::scale_color_manual(values = unname(pal[c(1,4)])) +
#       ggplot2::labs(x = xlab, y = ylab_con, color = "", linetype = "")
#
#   } else{
#     # plot for pd simulations
#
#     # add targets to data frame
#     tmp <- datasim$inf[,c("begin","pdt")]
#     names(tmp) <- c("time","pdp")
#     tmp$variable <- "Target"
#     for(nm in setdiff(names(df), names(tmp))){
#       tmp[,nm] <- NA
#     }
#     df <- rbind(df, tmp)
#     lv <- length(levels(df$variable))
#
#     out <- ggplot2::ggplot(df,
#                            ggplot2::aes(x = time,
#                                         y = pdp,
#                                         color = variable,
#                                         linetype = variable)) +
#       ggplot2::geom_line(ggplot2::aes(x = time, y = pdp, color = variable, linetype = variable),
#                          size = 1) +
#       ggplot2::geom_point(data = as.data.frame(datasim$sim),
#                           ggplot2::aes(x = time, y = pdobs),
#                           color = unname(pal["darkgrey"]),
#                           inherit.aes = FALSE, alpha = 0.5) +
#       ggplot2::scale_linetype_manual("", values = c(1:(lv-1),1)) +
#       ggplot2::scale_color_manual("", values = unname(c(pal[2:lv],pal[1]))) +
#       ggplot2::labs(x = xlab, y = ylab_resp) +
#       ggplot2::theme(legend.position="bottom")
#   }
#
#   return(out)
# }


# #' Plot method for bayessim objects
# #'
# #' Plot output returned by "bayes_control" function.
# #'
# #' @param x Object returned from "bayes_control" function
# #' @param ... \dots
# #' @return ggplot object displaying simulated data.
# #' @rdname plot
# #' @importFrom ggplot2 aes
# #' @export
# plot.bayessim <- function(x, ..., xlab = "Time", ylab_con = "Concentration", ylab_resp = "Response"){
#
#   pars_prior <- x$prior$pars_pkpd
#   pars_post <- rep(NA, length(pars_prior))
#   pars_post[x$prior$fixed_ix] <- pars_prior[x$prior$fixed_ix]
#   pars_post[is.na(pars_post)] <- exp(x$lpr[nrow(x$lpr),-ncol(x$lpr)])
#   names(pars_post) <- names(pars_prior)
#
#   # plot.datasim method
#   plot(x$dat,
#        pars_prior = pars_prior,
#        pars_post = pars_post,
#        pk_ix = x$prior$pk_ix,
#        pd_ix = x$prior$pd_ix,
#        xlab = xlab,
#        ylab_con = ylab_con,
#        ylab_resp = ylab_resp)
# }






## Methods for 'tciinf' -------------------------------------------------------
#
#
# #' Initialize an object with class "tciinf"
# #'
# #' @param tci_alg Function implementing TCI algorithm.
# #' @param pkmod A pkmod object created by `pkmod`.
# #' @return A list with class pkmod for which print, plot, predict, and simulate methods exist.
# #' @rdname init_pkmod
# #' @examples
# #' # create a tciinf object for plasma targeting (tci_plasma function) with a one compartment model
# #' my_mod <- pkmod(pars = c(CL = 10, V1 = 10), pkmod_fun = pkmod1cpt, init = 0)
# #' init_tciinf(tci_alg = tci_plasma, pkmod = my_mod)
# #' @export
# init_tciinf <- function(tci_alg, pkmod, ...){
# tciinf_obj <- c(list(tci_alg = tci_alg, pkmod = pkmod), list(...))
# class(tciinf_obj) <- "tciinf"
# return(tciinf_obj)
# }
#
#
# #' tciinf validation checks
# #'
# #' Function to provide validation checks for a tciinf object
# #' @param x Object of class "tciinf"
# #' @examples
# #' # create a tciinf object for plasma targeting (tci_plasma function) with a one compartment model
# #' my_mod <- pkmod(pars = c(CL = 10, V1 = 10), pkmod_fun = pkmod1cpt, init = 0)
# #' validate_tciinf(init_tciinf(tci_alg = tci_plasma, pkmod = my_mod, dtm = 1/6))
# #' @return Returns a list with class "tciinf" if validation checks are passed. Returns an error if not.
# #' @export
# validate_tciinf <- function(x){

## check that tci_alg is a function
# if(!is.function(x$tci_alg))
#   stop("tci_alg must be a function")
#
# ## check that tci_alg has the required arguments
# if(!all(c("Ct","pkmod","dtm") %in% names(formals(x$tci_alg))))
#   stop("tci_alg must have argument names c('Ct','pkmod','dtm)")
#
# if(is.null(x$dtm)) x$dtm <- eval(formals(x$tci_alg)$dtm)
#
# ## check that pkmod is a 'pkmod' object
# if(!inherits(x$pkmod,"pkmod"))
#   stop("pkmod must have class 'pkmod'")
#
# return(x)
# }
#
# #' Create an object with class "tciinf"
# #'
# #' @param tci_alg Function implementing TCI algorithm
# #' @param pkmod pkmod object created by `pkmod()`
# #' @return A list with class pkmod for which print, plot, predict, and simulate methods exist.
# #' @rdname init_pkmod
# #' @examples
# #' # create a tciinf object for plasma targeting (tci_plasma function) with a one compartment model
# #' my_mod <- pkmod(pars = c(CL = 10, V1 = 10), pkmod_fun = pkmod1cpt, init = 0)
# #' my_tci <- tciinf(tci_alg = tci_plasma, pkmod = my_mod)
# #' # hybrid/effect-site targeting
# #' mod_effect <- pkmod(pars = c(v1 = 8.995, v2 = 17.297, v3 = 120.963, cl = 1.382, q2 = 0.919, q3 = 0.609, ke0 = 1.289), pkmod = pkmod3cptm, pars_pd = c(c50 = 2.8, gamma = 1.47, gamma2 = 1.89, e0 = 93, emx = 93), pdmod = emax, pdinv = emax_inv, ecmpt = 4)
# #' my_tci_effect <- tciinf(tci_alg = tci_comb, pkmod = mod_effect)
# #' @export
# tciinf <- function(tci_alg, pkmod, ...){
# if(missing(tci_alg)){
#   stop("TCI algorithm must be specified. Options include 'tci_plasma', 'tci_effect',
#          ,'tci_comb' or a user-supplied algorithm.")
# }
# if(missing(pkmod)){
#   stop("A pkmod object must be specified through pkmod()")
# }
# tciinf_obj <- init_tciinf(tci_alg = tci_alg, pkmod = pkmod, ...)
# return(validate_tciinf(tciinf_obj))
# }
#
#
# #' Print method for tciinf objects
# #' @param x Object with class "tciinf"
# #' @examples
# #' my_mod <- pkmod(pars = c(CL = 10, V1 = 10), pkmod_fun = pkmod1cpt, init = 0)
# #' tciinf(tci_alg = tci_plasma, pkmod = my_mod, dtm = 1/6)
# #' @return Prints description of pkmod
# #' @export
# print.tciinf <- function(x){
# cat(paste0("TCI algorithm for a ", x$pkmod$ncmpt,"-compartment PK model"),"\n")
# cat(paste0("Update frequency: ",round(x$dtm,3)),"\n")
# }
#
#
#
#
# #' Predict method for tciinf objects.
# #'
# #' Apply a TCI algorithm and pkmod to a set of sequential targets.
# #' @param object An object created by `tciinf`.
# #' @param targets A matrix or data frame with columns 'value' and 'time'. Times
# #' indicate when the TCI algorithm should begin to reach each target. If the TCI
# #' algorithm should begin immediately, either the first time value must be 0 or
# #' the `inittm` parameter must be set.
# #' @param inittm Initial time to start TCI algorithm. Cannot be after target times.
# #' @param ignore_pd Logical. Should the PD component of the pkmod object (if present)
# #' be ignored. By default, predict.tciinf will assume that 'value' refers to PD
# #' targets if a PD model is specified.
# #' @param ... Arguments passed to TCI algorithm (tci_alg)
# #' @examples
# #' # plasma-targeting
# #' my_mod <- pkmod(pars = c(CL = 10, V1 = 10), pkmod_fun = pkmod1cpt, init = 0)
# #' my_tci <- tciinf(tci_alg = tci_plasma, pkmod = my_mod, dtm = 1/6)
# #' predict(my_tci, targets = cbind(value = c(2,3,4,4), time = c(0,2,3,10)))
# #' # effect-site targeting TCI algorithm with PD model
# #' my_mod_pd <- pkmod(pars = c(v1 = 8.995, v2 = 17.297, v3 = 120.963, cl = 1.382, q2 = 0.919, q3 = 0.609, ke0 = 1.289), pkmod = pkmod3cptm, pars_pd = c(c50 = 2.8, gamma = 1.47, gamma2 = 1.89, e0 = 93, emx = 93, sigma = 8.03, bis_delay = 28.263), pdmod = emax, pdinv = emax_inv, ecmpt = 4, init = c(0,0,0,0))
# #' my_tci_pd <- tciinf(tci_alg = tci_comb, pkmod = my_mod_pd, dtm = 1/6)
# #' predict(my_tci_pd, targets = cbind(value = c(70,60,50,50), time = c(0,2,3,10)))
# #' # effect-site targeting for a matrix of PK-PD parameter values
# #' evd_pkpd_pars <- eleveld_poppk(merge(eleveld_pk, eleveld_pd))[c(1,20,50,70,100),]
# #' evd_pkmod <- pkmod(pars = evd_pkpd_pars[,1:7], pars_pd = evd_pkpd_pars[,8:11], pkmod_fun = pkmod3cptm, ecmpt = 4, pdmod_fun = emax_eleveld, pdinv_fun = emax_inv_eleveld)
# #' my_tci_pd_matrix <- tciinf(tci_alg = tci_comb, pkmod = evd_pkmod, dtm = 1/6)
# #' predict(my_tci_pd_matrix, targets = cbind(value = c(70,60,50,50), time = c(0,2,3,10)))
# #' @rdname predict
# #' @export
# predict.tciinf <- function(object, targets, inittm = 0, ignore_pd = FALSE, ...){
#
# if(!all(c("value","time") %in% colnames(targets)))
#   stop("targets must have column names 'value' and 'time'")
#
# tms <- targets[,"time"]
# if(any(inittm > tms)) stop("inittm cannot be greater than any target times")
#
# if(!is.null(object$pkmod$pdmod_fun) & !ignore_pd){
#   if(is.matrix(object$pkmod$pars_pd)){
#     Ct <- apply(object$pkmod$pars_pd, 1, function(x) object$pkmod$pdinv(targets[,"value"],x))
#   } else{
#     Ct <- matrix(object$pkmod$pdinv(targets[,"value"], object$pkmod$pars_pd), ncol = 1)
#   }
# } else{
#   Ct <- matrix(targets[,"value"], ncol = 1)
# }
#
# # create step function to define targets at any point
# tms <- tms-inittm
#
# # define sequence of update times
# updatetms <- seq(0, max(tms), object$dtm)
#
# calc_infs <- function(j){
#
#   sf <- stepfun(tms, c(0,Ct[,j]))
#
#   if(is.matrix(object$pkmod$pars)){
#     pars_pk <- object$pkmod$pars[j,]
#   } else{
#     pars_pk <- object$pkmod$pars
#   }
#   if(is.matrix(object$pkmod$pars_pd)){
#     pars_pd <- object$pkmod$pars_pd[j,]
#   } else{
#     pars_pd <- object$pkmod$pars_pd
#   }
#
#   inf <- rep(NA, length(updatetms))
#   ini <- matrix(NA, nrow = object$pkmod$ncmpt, ncol = length(updatetms)+1)
#   ini[,1] <- object$pkmod$init
#
#   # iterate through times
#   for(i in 1:length(updatetms)){
#     inf[i] <- with(object, tci_alg(sf(updatetms[i]), pkmod = pkmod, dtm = dtm, pars = pars_pk,...))
#     ini[,i+1] <- with(object, pkmod$pkmod_fun(tm = dtm, kR = inf[i], pars = pars_pk, init = ini[,i]))
#     object$pkmod$init <- ini[,i+1]
#   }
#
#   startcon <- matrix(ini[,-ncol(ini)], ncol = nrow(ini), nrow = ncol(ini)-1, byrow = TRUE)
#   endcon <- matrix(ini[,-1], ncol = nrow(ini), nrow = ncol(ini)-1, byrow = TRUE)
#   dose <- create_intvl(times = updatetms+inittm, infrt = inf, duration = object$dtm)
#   out <- cbind(dose, sf(updatetms), startcon, endcon)
#   colnames(out) <- c("infrt","begin","end","Ct",paste0("c",1:object$pkmod$ncmpt, "_start"),
#                      paste0("c",1:object$pkmod$ncmpt, "_end"))
#
#   if(!is.null(object$pkmod$pdmod_fun) & !ignore_pd){
#     pdt <- object$pkmod$pdmod_fun(out[,"Ct"], pars_pd)
#     pdresp_start <- object$pkmod$pdmod_fun(out[,paste0("c",object$pkmod$ecmpt,"_start")], pars_pd)
#     pdresp_end <- object$pkmod$pdmod_fun(out[,paste0("c",object$pkmod$ecmpt,"_end")], pars_pd)
#
#     out <- cbind(out,
#                  pdt = pdt,
#                  pdresp_start = pdresp_start,
#                  pdresp_end = pdresp_end)
#   }
#
#   return(out)
# }
#
# if(is.vector(object$pkmod$pars)){
#   out <- calc_infs(1)
# } else{
#   out <- lapply(1:nrow(object$pkmod$pars), function(j) calc_infs(j))
#   id <- rep(1:nrow(object$pkmod$pars), each = nrow(out[[1]]))
#   out <- cbind(id = id, do.call("rbind", out))
# }
#
# return(out)
# }




# ## Methods for pkmod objects
# #' Create an object with class "pkmod"
# #' @param ncmpt Number of compartments used for mammillary compartment model. Options "one","two","three", and "three_effect"
# #' link to functions `pkmod1cpt`, `pkmod2cpt`, `pkmod3cpt`, and `pkmod3cptm`, respectively. Can be overruled by providing
# #' a user-defined PK function to `pkmod`.
# #' @param pkmod_fun PK model function that returns predicted concentrations for a single infusion rate.
# #' Can be user-defined or a package function: pkmod1cpt, pkmod2cpt, pkmod3cpt,
# #' for 1-, 2-, and 3-compartment models or pkmod3cptm for a 3-compartment metabolite model.
# #' @param pars Vector or matrix of parameters for pkmod object
# #' @param init Initial concentrations. Will default to values of zero in all compartments if not specified.
# #' @param pars_pd PD model parameters if a PD model is specified
# #' @param pdmod_fun PD model function
# #' @param pdinv_fun Inverse PD model for use in TCI algorithms
# #' @param pcmpt Index of plasma compartment. Defaults to first compartment if not specified.
# #' @param ecmpt Index of effect-site compartment if a PD model is specified. Will default to last compartment if unspecified.
# #' @param sigma_add Standard deviation of additive residual error.
# #' @param sigma_mult Standard deviation of multiplicative residual error.
# #' @param log_response Logical. Should the response be natural-log transformed before error is applied? Defaults to FALSE.
# #' @param delay Delay in PD observations. Defaults to 0.
# #' @param omega_pk Optional vector of random effect variances for PK parameters.
# #' Values should be set to zero if parameters are fixed in the population. It is
# #'  assumed that random effects are log-normally distributed.
# #' @param omega_pd Optional vector of random effect variances for PD parameters.
# #' Values should be set to zero if parameters are fixed in the population. It is
# #'  assumed that random effects are log-normally distributed.
# #' @param omega_sigma Variance of random effect on additive residual error.
# #' @param max_pdval Maximum allowable value of PD response.
# #' @param min_pdval Minimum allowable value of PD response.
# #' @return A list with class pkmod for which print, plot, predict, and simulate methods exist.
# #' @rdname init_pkmod
# #' @examples
# #' # create a pkmod object for a one compartment model
# #' init_pkmod(pars = c(CL = 10, V1 = 10), pkmod_fun = pkmod1cpt, init = 0)
# #' # create a joint PK-PD model
# #' init_pkmod(pars = c(v1 = 8.995, v2 = 17.297, v3 = 120.963, cl = 1.382, q2 = 0.919, q3 = 0.609, ke0 = 1.289), pkmod = pkmod3cptm, pars_pd = c(c50 = 2.8, gamma = 1.47, gamma2 = 1.89, e0 = 93, emx = 93), pdmod = emax, pdinv = emax_inv, ecmpt = 4)
# #' @export
# init_pkmod <- function(ncmpt = NULL, pars = NULL, pkmod_fun = NULL, init = NULL, pars_pd = NULL, pdmod_fun = NULL,
#                        pdinv_fun = NULL, pcmpt = NULL, ecmpt = NULL, sigma_add = NULL,
#                        sigma_mult = NULL, log_response = NULL, delay = NULL,
#                        pars_pk0 = NULL, pars_pd0 = NULL, omega_pk = NULL, omega_pd = NULL,
#                        omega_sigma = NULL, max_pdval = NULL, min_pdval = NULL){
#
# if(!is.function(pkmod_fun)) stop("`pkmod_fun` must be a function. Available models in `tci` are 'pkmod1cpt',
#        'pkmod2cpt','pkmod3cpt','pkmod3cmptm'. Custom models can be specified.
#        See vignettes for examples.")
#
# if(!("init" %in% names(formals(pkmod_fun)))) stop("pkmod_fun must include a vector argument
#                                                     'init' with initial concentrations.")
#
# if("ID" %in% toupper(colnames(pars))){
#   ids <- pars[,which(toupper(colnames(pars)) == "ID")]
#   pars <- pars[,-which(toupper(colnames(pars)) == "ID")]
# } else{
#   ids <- 1:ifelse(is.matrix(pars), nrow(pars), 1)
# }
#
# ## check for initial concentrations and set to zero if not specified
# ncmpt <- length(eval(formals(pkmod_fun)$init))
# if(is.null(init)) {
#   init <- rep(0,ncmpt)
# }
# if(is.null(pcmpt)) pcmpt = 1
#
# pkmod_obj <- list(ids = ids,
#                   pars = pars,
#                   pkmod_fun = pkmod_fun,
#                   init = init,
#                   ncmpt = ncmpt,
#                   pars_pd = pars_pd,
#                   pdmod_fun = pdmod_fun,
#                   pdinv_fun = pdinv_fun,
#                   pcmpt = pcmpt,
#                   ecmpt = ecmpt,
#                   sigma_add = sigma_add,
#                   sigma_mult = sigma_mult,
#                   log_response = log_response,
#                   delay = delay,
#                   omega_pk = omega_pk,
#                   omega_pd = omega_pd,
#                   omega_sigma = omega_sigma,
#                   max_pdval = max_pdval,
#                   min_pdval = min_pdval)
#
# class(pkmod_obj) <- "pkmod"
# return(pkmod_obj)
# }
#
#
# #' pkmod validation checks
# #'
# #' Function to provide validation checks for a pkmod object
# #' @param x Object of class "pkmod"
# #' @examples
# #' # parameters must include V1 for one compartment model
# #' my_mod <- init_pkmod(pars = c(CL = 10, V = 10), pkmod_fun = pkmod1cpt, init = 0)
# #' validate_pkmod(my_mod)
# #' my_mod2 <- init_pkmod(pars = c(CL = 10, V1 = 10), pkmod_fun = pkmod1cpt, init = 0)
# #' validate_pkmod(my_mod2)
# #' # Schnider population PK model for propofol
# #' dat <- data.frame(AGE  = c(20,40,65),
# #'                   TBM  = c(50,70,90),
# #'                   HGT  = c(150,170,200),
# #'                   MALE = c(TRUE,FALSE,TRUE))
# #' my_mod3 <- init_pkmod(pars = schnider_poppk(dat), pkmod_fun = pkmod3cptm)
# #' validate_pkmod(my_mod3)
# #' @return Returns a list with class "pkmod" if validation checks are passed. Returns an error if not.
# #' @export
# validate_pkmod <- function(x){

## verify parameters are named correctly (case indifferent)
# if(is.vector(x$pars))
#   parnms <- tolower(names(x$pars))
# else
#   parnms <- tolower(colnames(x$pars))
#
# if(x$ncmpt == 1){
#   if(any(!(c("k10","v1") %in% parnms)) & any(!(c("cl","v1") %in% parnms)))
#     stop("pars must have names ('k10','v1') or ('cl','v1')")
# }
# if(x$ncmpt == 2){
#   if(any(!(c("k10","k12","k21","v1","v2") %in% parnms)) & any(!(c("cl","q2","v1","v2") %in% parnms)))
#     stop("pars must have names ('k10','k12','k21','v1','v2') or ('cl','q2','v1','v2')")
# }
# if(x$ncmpt == 3){
#   if(any(!(c("k10","k12","k21","k13","k31","v1","v2","v3") %in% parnms)) & any(!(c("cl","q2","q3","v1","v2","v3") %in% parnms)))
#     stop("pars must have names ('k10','k12','k21','k13','k31','v1','v2','v3') or ('cl','q2','q3','v1','v2','v3')")
# }
# if(x$ncmpt == 4){
#   if(any(!(c("k10","k12","k21","k13","k31","ke0","v1","v2","v3") %in% parnms)) &
#      any(!(c("cl","q2","q3","ke0","v1","v2","v3") %in% parnms)))
#     stop("pars must have names ('k10','k12','k21','k13','k31','v1','v2','v3','ke0') or ('cl','q2','q3','v1','v2','v3','ke0')")
# }
#
# ## effect-site assignment/warning if PD is specified and ecmpt is not
# if(!is.null(x$pdmod) & is.null(x$ecmpt)){
#   # warning("Effect-site compartment is not specified. Defaulting to last compartment: ecmpt = ncmpt")
#   x$ecmpt <- x$ncmpt
# }
#
# ## check that PD parameters are specified if PD model is
# if(!is.null(x$pdmod) & is.null(x$pars_pd)){
#   stop("pars_pd must be specified if pdmod is non-null")
#   x$ecmpt <- x$ncmpt
#
#   ## check ordering of parameters in PD model
#   if(!all.equal(names(formals(x$pdmod)), c("ce","pars"))){
#     stop("pdmod must have arguments named c('ce','pars')")
#   }
# }
#
# return(x)
# }
#
#
# #' Create pkmod object
# #'
# #' User function to create pkmod objects with validation checks
# #'
# #' @param pars Vector or matrix of parameters for pkmod object
# #' @param pkmod_fun PK model function that returns predicted concentrations for a single infusion rate
# #' @param init initial concentrations
# #' @param pars_pd PD model parameters if a PD model is specified
# #' @param pdmod_fun PD model function
# #' @param pdinv_fun Inverse PD model for use in TCI algorithms
# #' @param pcmpt Index of plasma compartment. Defaults to first compartment if not specified.
# #' @param ecmpt Index of effect-site compartment if a PD model is specified. Will default to last compartment if unspecified.
# #' @param sigma_add Standard deviation of additive residual error.
# #' @param sigma_mult Standard deviation of multiplicative residual error.
# #' @param log_response Logical. Should the response be natural-log transformed before error is applied? Defaults to FALSE.
# #' @param delay Delay in PD observations. Defaults to 0.
# #' @param omega_pk Optional vector of random effect variances for PK parameters. Values should be set to zero if parameters are fixed in the population. It is assumed that random effects are log-normally distributed.
# #' @param omega_pd Optional vector of random effect variances for PD parameters. Values should be set to zero if parameters are fixed in the population. It is  assumed that random effects are log-normally distributed.
# #' @param omega_sigma Variance of random effect on additive residual error.
# #' @param max_pdval Maximum allowable value of PD response.
# #' @param min_pdval Minimum allowable value of PD response.
# #' @return A list with class pkmod for which print, plot, predict, and simulate methods exist.
# #' @rdname pkmod
# #' @return Returns a list with class "pkmod" if validation checks are passed. Returns an error if not.
# #' @examples
# #' # create a one compartment pkmod object
# #' pkmod(pars = c(CL = 10, V1 = 10), pkmod_fun = pkmod1cpt, init = 0)
# #' # create a three compartment pkmod object with effect site and an emax PD component
# #' pkmod(pars = c(v1 = 8.995, v2 = 17.297, v3 = 120.963, cl = 1.382, q2 = 0.919, q3 = 0.609, ke0 = 1.289), pkmod = pkmod3cptm, pars_pd = c(c50 = 2.8, gamma = 1.47, gamma2 = 1.89, e0 = 93, emx = 93), pdmod = emax, pdinv = emax_inv, ecmpt = 4)
# #' @export
# pkmod <- function(pkmod_fun, pars, init = NULL, pars_pd = NULL, pdmod_fun = NULL,
#                   pdinv_fun = NULL, ecmpt = NULL, sigma_add = NULL,
#                   sigma_mult = NULL, log_response = NULL,
#                   delay = NULL, pars_pk0 = NULL, pars_pd0 = NULL,
#                   omega_pk = NULL, omega_pd = NULL, omega_sigma = NULL, max_pdval = NULL,
#                   min_pdval = NULL){
# if(missing(pkmod_fun)){
#   stop("pkmod_fun is missing. Options in package include 'pkmod1cpt', 'pkmod2cpt',
#          'pkmod3cpt', and 'pkmod3cptm'. User-defined functions can also be supplied.")
# }
# new_mod <- init_pkmod(pkmod_fun = pkmod_fun,
#                       pars = pars,
#                       init = init,
#                       pars_pd = pars_pd,
#                       pdmod_fun = pdmod_fun,
#                       pdinv_fun = pdinv_fun,
#                       ecmpt = ecmpt,
#                       sigma_add = sigma_add,
#                       sigma_mult = sigma_mult,
#                       log_response = log_response,
#                       delay = delay,
#                       omega_pk = omega_pk,
#                       omega_pd = omega_pd,
#                       omega_sigma = omega_sigma,
#                       max_pdval = max_pdval,
#                       min_pdval = min_pdval)
# validate_pkmod(new_mod)
# }
#
# #' Print method for pkmod objects
# #' @param x Object with class "pkmod"
# #' @examples
# #' # 1-compartment PK model
# #' my_mod <- pkmod(pars = c(CL = 10, V1 = 10), pkmod = pkmod1cpt, init = 0)
# #' print(my_mod)
# #'
# #' # eleveld population PK-PD model
# #' evd_pkpd_pars <- eleveld_poppk(merge(eleveld_pk, eleveld_pd))
# #' # random effect standard deviations
# #' omega_pk <- c(0.610,0.565,0.597,0.265,0.346,0.209,0.463,0.702)
# #' omega_pd <- c(0.242,0,0,0)
# #' evd_pkmod <- pkmod(pars = evd_pkpd_pars[,1:7], pars_pd = evd_pkpd_pars[,8:11], pkmod_fun = pkmod3cptm, ecmpt = 4, pdmod_fun = emax_eleveld, pdinv_fun = emax_inv_eleveld, sigma_add = 8.03, omega_pk = omega_pk, omega_pd = omega_pd, omega_sigma = 0.230)
# #' print(evd_pkmod)
# #' @return Prints description of pkmod
# #' @export
# print.pkmod <- function(x){
# cat("--- PK model -------------------------------------------", "\n")
# cat(paste0(x$ncmpt,"-compartment PK model"),"\n")
# cat(paste0("Number of observations: ", length(unique(x$ids))),"\n")
# if(is.vector(x$pars)){
#   cat(paste0("PK parameters: ",paste(names(x$pars), "=", x$pars,collapse = ", ")),"\n")
# } else{
#   cat(paste0("PK parameters: ", paste(colnames(x$pars), collapse = ",")),"\n")
# }
# cat(paste0("Initial concentrations: ",paste0("(",paste(x$init,collapse = ","),")")),"\n")
# cat(paste0("Plasma compartment: ",x$pcmpt),"\n")
# if(!is.null(x$ecmpt)) cat(paste0("Effect compartment: ",x$ecmpt),"\n")
# if(!all(sapply(x[c("pars_pd","delay")],is.null)))
#   cat("--- PD model -------------------------------------------", "\n")
#
# if(!is.null(x$pars_pd)){
#   if(is.vector(x$pars_pd)){
#     cat(paste0("PD parameters: ",paste(names(x$pars_pd), "=", x$pars_pd,collapse = ", ")),"\n")
#   } else{
#     cat(paste0("PD parameters: ", paste(colnames(x$pars_pd), collapse = ",")),"\n")
#   }
# }
# if(!is.null(x$delay)) cat(paste0("Observation delay: ",x$delay),"\n")
# if(!is.null(x$min_pdval)) cat(paste0("Minimum PD value: ",x$min_pdval),"\n")
# if(!is.null(x$max_pdval)) cat(paste0("Maximum PD value: ",x$max_pdval),"\n")
#
# if(!all(sapply(x[c("sigma_add","sigma_mult","log_response","omega_pk","omega_pd")],is.null)))
#   cat("--- Simulations ----------------------------------------", "\n")
# if(!is.null(x$sigma_add)) cat(paste0("Additive residual error: ",x$sigma_add),"\n")
# if(!is.null(x$sigma_mult)) cat(paste0("Multiplicative residual error: ",x$sigma_mult),"\n")
# if(!is.null(x$log_response)) cat(paste0("Logged response: ",x$log_response),"\n")
# if(!is.null(x$omega_pk)) cat(paste0("Random effect SD for PK: ",paste0("(",paste(x$omega_pk,collapse = ","),")")),"\n")
# if(!is.null(x$omega_pd)) cat(paste0("Random effect SD for PD: ",paste0("(",paste(x$omega_pd,collapse = ","),")")),"\n")
# if(!is.null(x$omega_sigma)) cat(paste0("Random effect SD for residual error: ",x$omega_sigma),"\n")
# }
#
#
# #' Update method for pkmod
# #'
# #' Update parameters or initial values of a pkmod object
# #' @param object Object with class pkmod
# #' @param ... Updated values passed to object.
# #' @return Returns the pkmod object with elements replaced
# #' @examples
# #' # initial pkmod object
# #' my_mod <- pkmod(pars = c(CL = 10, V1 = 10), pkmod = pkmod1cpt, init = 0)
# #' print(my_mod)
# #' # update a subset of parameters and initial values. Add multiplicative error for simulations.
# #' update(my_mod, pars = c(CL = 20), init = 3, sigma_mult = 0.2)
# #' @export
# update.pkmod <- function(object, ...){#

# update_elements <- list(...)
# if(length(update_elements) == 0) return(object)
# if(is.list(update_elements[[1]])) update_elements <- update_elements[[1]]

# if(any(!names(update_elements) %in% c(names(formals(pkmod))))){
#   nms_mismatch <- names(update_elements)[!names(update_elements) %in% names(formals(pkmod))]
#   warning(paste(paste(nms_mismatch, collapse = ", "), "are not valid names."))
# }

# # subset to elements with valid names
# update_elements <- update_elements[names(update_elements) %in% c(names(formals(pkmod)))]

# for(nm in names(update_elements)){
#   # partial update - if matrix, assumes that the whole matrix will be replaced
#   if(nm %in% c("pars","pars_pd") & is.vector(object[[nm]])){
#     object[[nm]][names(update_elements[[nm]])] <- update_elements[[nm]]
#   } else{
#     object[[nm]] <- update_elements[[nm]]
#   }
# }

# return(validate_pkmod(object))
# }
