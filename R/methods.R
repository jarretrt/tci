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


#' Print poppkmod
#'
#' Print method for poppkmod objects
#' @param x Object with class "pkmod" created by `poppkmod()`
#' @param digits Number of significant digits to print
#' @param ... Additional arguments. Not used
#' @examples
#' data <- data.frame(ID = 1:5,
#' AGE = seq(20,60,by=10),
#' TBW = seq(60,80,by=5),
#' HGT = seq(150,190,by=10),
#' MALE = c(TRUE,TRUE,FALSE,FALSE,FALSE))
#' poppkmod(data, drug = "ppf", model = "eleveld")
#' poppkmod(data, drug = "ppf", model = "eleveld")
#' @return Prints description of pkmod
#' @export
print.poppkmod <- function(x, ..., digits = 3){

  cat("--- Model summary --------------------------------------", "\n")
  drug_type <- ifelse(x$drug == "ppf","propofol","remifentanil")
  mod_name <- x$model
  substr(mod_name,1,1) <- toupper(substr(mod_name,1,1))
  cat(paste(mod_name, "population model for",drug_type), "\n")
  cat(paste("Number of individuals:",length(x$pkmods)), "\n")

  cat("--- PK model -------------------------------------------", "\n")
  pkstr <- ifelse(x$pkmods[[1]]$ncmpt == 4, "3 compartment-effect", paste(x$pkmods[[1]]$ncmpt,"compartment"))
  cat(paste("Structural PK model:", pkstr), "\n")
  cat(paste("PK parameter names:", paste(names(x$pkmods[[1]]$pars_pk), collapse = " ")),"\n")
  tci_type <- ifelse("ke0" %in% tolower(names(x$pkmods[[1]]$pars_pk)), "Effect-site","Plamsa")
  cat(paste("Default TCI targeting:",tci_type),"\n")

  cat("--- PD model -------------------------------------------", "\n")
  cat(paste("Structural PD model:", ifelse(!is.null(x$pkmods[[1]]$pdfn), "Emax", "None")),"\n")
  cat(paste("PD parameter names:", paste(names(x$pkmods[[1]]$pars_pd), collapse = " ")),"\n")
  if(!is.null(x$pkmods[[1]]$pars_pd)){
    vals_pd <- summary(t(sapply(x$pkmods, `[[`, "pars_pd")), digits = digits)[c(1,3,6),]
    knitr::kable(vals_pd, format = "rst")
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
      # object[[nm]][(names(update_elements[[nm]]))] <- update_elements[[nm]]
      object[[nm]][match(tolower(names(update_elements[[nm]])), tolower(names(object[[nm]])))] <- update_elements[[nm]]
    } else{
      object[[nm]] <- update_elements[[nm]]
    }
  }

  return(validate_pkmod(object))
}




#' Predict method for pkmod objects
#'
#' Predict concentrations from a pkmod object - can be a user defined function
#'
#' @param object An object with class pkmod.
#' @param inf A matrix with columns "begin","end","inf_rate" indicating when infusions should
#' be administered. Can be created by `inf_manual` or `inf_tci`.
#' @param tms Times at which to calculate predicted concentrations.
#' @param ... List or vector of values to be passed on to update.pkmod.
#' @return Matrix of predicted concentrations associated with a pkmod object and
#' and infusion schedule.
#' @examples
#' # dosing schedule
#' dose <- inf_manual(inf_tms = c(0,0.5,4,4.5,10), inf_rate = c(100,0,80,0,0))
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
  if(!all(c("inf_rate","begin","end") %in% colnames(inf)))
    stop("inf must include 'inf_rate','begin','end' as column names")

  # apply Rcpp implementation to selected PK models
  if(any(identical(object$pkfn, pkmod1cpt),
         identical(object$pkfn, pkmod2cpt),
         identical(object$pkfn, pkmod3cpt),
         identical(object$pkfn, pkmod3cptm))){

    starttm <- inf[1,"begin"]
    begin <- inf[,"begin"] - starttm
    end <- inf[,"end"] - starttm
    infs <- inf[,"inf_rate"]
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

    if(inf[nrow(inf),"inf_rate"] !=0 & inf[nrow(inf),"end"]<max(tms)){
      inf <- rbind(inf, c(begin = inf[nrow(inf),"end"], end = max(tms), inf_rate = 0))
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
                                       kR = inf[i,"inf_rate"],
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


#' Predict method for pkmod objects
#'
#' Predict concentrations from a pkmod object - can be a user defined function
#'
#' @param object An object with class poppkmod.
#' @param inf A matrix with columns "begin","end","inf_rate" indicating when infusions should
#' be administered. Can be created by `inf_manual` or `inf_tci`.
#' @param tms Times at which to calculate predicted concentrations.
#' @param ... List or vector of values to be passed on to predict.pkmod.
#' @return Matrix of predicted concentrations associated with a poppkmod object and
#' and infusion schedule.
#' @examples
#' # dosing schedule
#' dose <- inf_manual(inf_tms = c(0,0.5,4,4.5,10), inf_rate = c(100,0,80,0,0))
#' # poppkmod object
#' data <- data.frame(ID = 1:5, AGE = seq(20,60,by=10), TBW = seq(60,80,by=5),
#' HGT = seq(150,190,by=10), MALE = c(TRUE,TRUE,FALSE,FALSE,FALSE))
#' elvd_mod <- poppkmod(data, drug = "ppf", model = "eleveld")
#' predict(elvd_mod, inf = dose, tms = c(1.5,2.5,3))
#' @export
predict.poppkmod <- function(object, inf, tms, ...){
  if(is.data.frame(inf)) inf <- as.matrix(inf)
  if("id" %in% colnames(inf)){
    inf <- lapply(split(inf, f = as.factor(inf[,"id"])), matrix, ncol = ncol(inf), dimnames = list(NULL, colnames(inf)))
    pred <- do.call("rbind", lapply(1:length(object$pkmods), function(i) predict(object$pkmods[[i]], inf = inf[[i]], tms = tms, ...)))
  } else{
    pred <- do.call("rbind", lapply(1:length(object$pkmods), function(i) predict(object$pkmods[[i]], inf = inf, tms = tms, ...)))
  }
  cbind(id=rep(object$ids, each = length(tms)),
        time=rep(tms, times = length(object$ids)),
        pred)
}

#' Simulate method for pkmod objects
#'
#' Simulate observations from a pkmod object
#'
#' @param object An object with class "pkmod" generated by `pkmod`
#' @param nsim Number of observations to simulate at each time point. Defaults to 1.
#' @param seed An integer used to initialize the random number generator.
#' @param ... Arguments passed to `update.pkmod`.
#' @param inf A matrix of infusion rates with columns 'begin', 'end', and 'inf_rate'. This
#' can be created manually, by `inf_manual`, or by `inf_tci`.
#' @param tms Times at which to simulate observations.
#' @param obs_cmpt Integer value indicating compartment in which observations are taken.
#' Overridden if a PD model is included.
#' @param resp_bounds Optional vector of two values indicating minimum and maximum values possible for the response.
#' @examples
#' # simulate data from a 2 compartment model with multiplicative error
#' dose <- inf_manual(inf_tms = c(0,0.5,4,4.5,10), inf_rate = c(100,0,80,0,0))
#' my_mod <- pkmod(pars_pk = c(CL = 10, V1 = 10, Q2 = 4, V2 = 30))
#' inf <- inf_tci(my_mod, target_vals = c(2,3,4,4), target_tms = c(0,2,3,10), "plasma")
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



#' Simulate method for poppkmod objects
#'
#' Simulate observations from a poppkmod object
#'
#' @param object An object with class "poppkmod" generated by `poppkmod`
#' @param nsim Number of observations to simulate at each time point. Defaults to 1.
#' @param seed An integer used to initialize the random number generator.
#' @param ... Arguments passed to `update.pkmod`.
#' @param inf A matrix of infusion rates with columns 'begin', 'end', and 'inf_rate'. This
#' can be created manually, by `inf_manual`, or by `inf_tci`.
#' @param tms Times at which to simulate observations.
#' @param obs_cmpt Integer value indicating compartment in which observations are taken.
#' Overridden if a PD model is included.
#' @param resp_bounds Optional vector of two values indicating minimum and maximum values possible for the response.
#' @examples
#' # simulate data from a 2 compartment model with multiplicative error
#' dose <- inf_manual(inf_tms = c(0,0.5,4,4.5,10), inf_rate = c(100,0,80,0,0))
#' # poppkmod object
#' data <- data.frame(ID = 1:5, AGE = seq(20,60,by=10), TBW = seq(60,80,by=5),
#'  HGT = seq(150,190,by=10), MALE = c(TRUE,TRUE,FALSE,FALSE,FALSE))
#' elvd_mod <- poppkmod(data, drug = "ppf", model = "eleveld")
#' inf <- inf_tci(elvd_mod, target_vals = c(2,3,4,4), target_tms = c(0,2,3,10), "plasma")
#' simulate(elvd_mod, nsim = 3, inf = inf, tms = c(1,2,4,6,10))
#' @export
simulate.poppkmod <- function(object, nsim = 1, seed = NULL, ..., inf, tms,
                              obs_cmpt = 1, resp_bounds = NULL){

  if(!is.null(seed)) set.seed(seed)
  inf <- lapply(split(inf, f = as.factor(inf[,"id"])), matrix, ncol = ncol(inf), dimnames = list(NULL, colnames(inf)))
  sims <- do.call("rbind",lapply(1:length(object$pkmods), function(i){
    simulate(object$pkmods[[i]], nsim = nsim, inf = inf[[i]], tms = tms,
             obs_cmpt=obs_cmpt, resp_bounds=resp_bounds, ...)}))
  colnames(sims) <- paste0("sim",1:ncol(sims))
  cbind(id=rep(object$ids, each = length(tms)),
        time=rep(tms, times = length(object$ids)),
        sims)
}



## -- Plotting methods ---------------------------------------------------------

#' Plot method for 'pkmod'
#'
#' Will show predicted concentrations in compartments associated with an infusion schedule.
#' @param x An object with class pkmod.
#' @param ... Arguments passed on to predict.pkmod. Could include changing initial concentrations or parameter values through pkmod_update. See ?predict.pkmod
#' @param inf An infusion schedule object with columns "begin","end","inf_rate".
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
#' dose <- inf_manual(inf_tms = c(0,0.5,4,4.5,10), inf_rate = c(100,0,80,0,0))
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


plot.tcisim <- function(tcisim){

}
