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

  cat("tci pkmod object\n")
  cat("See ?update.pkmod to modify or add elements\n")
  x <- update(x,...)
  cat("\nPK model \n")
  cat(paste0(" ", x$ncmpt,"-compartment PK model"),"\n")
  if(!is.null(x$pars_pk))
    cat(paste0(" PK parameters: ",
               paste(names(x$pars_pk), "=", signif(x$pars_pk,digits),
                     collapse = ", ")),"\n")
  if(!is.null(x$init))
    cat(paste0(" Initial concentrations: ",paste0("(",paste(signif(x$init,digits),
                                                           collapse = ","),")")),"\n")

  cat(paste0(" Plasma compartment: ",x$pcmpt),"\n")
  if(!is.null(x$ecmpt)) cat(paste0(" Effect compartment: ",x$ecmpt),"\n")
  if(!is.null(x$pars_pd) & length(x$pars_pd) > 0){
    cat("\nPD model \n")
    cat(paste0(" PD parameters: ",paste(names(x$pars_pd), "=", signif(x$pars_pd,digits),
                                       collapse = ", ")),"\n")
  }
  if(!is.null(x$sigma_add)|!is.null(x$sigma_mult)|!is.null(x$log_response)){
    cat("\nSimulation\n")
  }

  if(!is.null(x$sigma_add)) cat(paste0(" Additive error SD: ",signif(x$sigma_add,digits)),"\n")
  if(!is.null(x$sigma_mult)) cat(paste0(" Multiplicative error SD: ",signif(x$sigma_mult,digits)),"\n")
  if(!is.null(x$log_response)) cat(paste0(" Logged response: ",x$log_response),"\n")
  if(!is.null(x$Omega)){
    cat(" Random effects","\n")
    print(signif(x$Omega,digits))
  }

  # cat("\nUpdate Method\n")
  # cat(" Modify or add elements via 'update()' \n")
  # cat(" e.g., update(pkmod, pars_pk = c(CL=5, V1=2), init = 1) \n")
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
#' update(my_mod, pars_pk = c(CL = 20, V1 = 2), init = 3, sigma_add = 1, log_response = TRUE)
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
#' @param return_times Logical. Should prediction times be returned along with responses? Defaults to FALSE.
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
predict.pkmod <- function(object, inf, tms, return_times = FALSE, ...){

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

  if(return_times) pred <- cbind(time = tms, pred)

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



## -- poppkmod methods ---------------------------------------------------------

#' @name init_poppkmod
#' @title Initialize a `poppkmod` object.
#'
#' @description Generate a `poppkmod` object from an existing population PK model
#' for propofol or remifentanil using patient covariates. Available models for
#' propofol are the Marsh, Schnider, and Eleveld models. Available models for
#' remifentanil are the Minto, Kim, and Eleveld models. Input is
#' a data frame with rows corresponding to individuals and columns
#' recording patient covariates.
#' @param data Data frame of patient covariates. ID values, if used, should be in a column labeled "id" or "ID"
#' @param drug "ppf" for propofol or "remi" for remifentanil. Defaults to "ppf".
#' @param model Model name. Options are "marsh", "schnider", or "eleveld" if
#' drug = "ppf", or "minto", "kim", or "eleveld" if drug = "remi".
#' @param sample Logical. Should parameter values be sampled from interindividual distribution (TRUE)
#' or evaluated at point estimates for covariates (FALSE)? Defaults to FALSE.
#' @param PD Should the PD component be evaluated for PK-PD models. Defaults to TRUE.
#' @return `poppkmod` object
#' @examples
#' data <- data.frame(ID = 1:5, AGE = seq(20,60,by=10), TBW = seq(60,80,by=5),
#' HGT = seq(150,190,by=10), MALE = c(TRUE,TRUE,FALSE,FALSE,FALSE))
#' init_poppkmod(data, drug = "ppf", model = "eleveld")
#' init_poppkmod(data, drug = "remi", model = "kim")
#' @export
init_poppkmod <- function(data=NULL, drug = NULL, model = NULL, sample = NULL, PD = NULL){
  poppkmod_obj <- list(data = data)
  attr(poppkmod_obj,"drug") <- drug
  attr(poppkmod_obj,"model") <- model
  class(poppkmod_obj) <- "poppkmod"
  return(poppkmod_obj)
}


#' @name validate_poppkmod
#' @title Perform validation checks on a `poppkmod` object
#'
#' @description Perform validation checks on a `poppkmod` object created by `init_poppkmod`.
#' @param x Object with class "poppkmod" created by init_poppkmod
#' @return `poppkmod` object or error
#' @examples
#' data <- data.frame(ID = 1:5, AGE = seq(20,60,by=10), TBW = seq(60,80,by=5),
#' HGT = seq(150,190,by=10), MALE = c(TRUE,TRUE,FALSE,FALSE,FALSE))
#' validate_poppkmod(init_poppkmod(data, drug = "ppf", model = "eleveld"))
#' @export
validate_poppkmod <- function(x){

  if(!inherits(x, "poppkmod"))
    stop("x must have class 'poppkmod'")

  drug <- attr(x,"drug")
  model <- attr(x,"model")

    if(drug == "ppf" & model %in% c("minto","kim"))
      stop("'drug' must be set to 'remi' to use Minto or Kim models")

    if(drug == "remi" & model %in% c("marsh","schnider"))
      stop("'drug' must be set to 'ppf' to use Marsh or Schnider models")

    if("id" %in% tolower(names(x$data))){
      ids <- x$data[,grep("id",names(x$data), ignore.case = TRUE)]
    } else{
      ids <- 1:nrow(x$data)
    }

    if(!is.data.frame(x$data)){
      warning("Converting 'data' to a data frame")
      x$data <- as.data.frame(x$data)
    }

    if(model == "marsh"){
      if(!("TBW" %in% names(x$data))) stop("data must have column 'TBW'")
    }

    if(model == "schnider"){
      if(any(!(c("AGE","HGT") %in% names(x$data))) | !("LBM" %in% names(x$data) | all(c("TBW","MALE") %in% names(x$data))))
        stop("data must have columns 'AGE','HGT' and either 'TBW' and 'MALE' or 'LBW'")
    }

    if(model == "eleveld" & drug == "ppf"){
      if(any(!(c("AGE","TBW","HGT","MALE") %in% names(x$data))))
        stop("data must have columns 'AGE','TBW','HGT' and 'MALE'")
    }

    if(model == "minto"){
      if(any(!(c("AGE") %in% names(x$data))) | !("LBM" %in% names(x$data) | all(c("TBW","MALE","HGT") %in% names(x$data))))
        stop("data must have columns 'AGE','HGT' and either 'TBW' and 'MALE' or 'LBW'")
    }

    if(model == "kim"){
      if(any(!(c("AGE","TBW") %in% names(x$data))) | (all(!c("BMI","HGT") %in% names(x$data)) & !"FFM" %in% names(x$data)))
        stop("data must have columns 'AGE','TBW' and either 'BMI' and 'HGT' or 'FFM'")
    }

    if(model == "eleveld" & drug == "remi"){
      if(any(!(c("AGE","MALE","TBW") %in% names(x$data))) | (all(!c("BMI","HGT") %in% names(x$data)) & !"FFM" %in% names(x$data)))
        stop("data must have columns 'AGE', 'MALE', and 'TBW' and either 'HGT' or 'BMI'")
    }

  return(x)
}


#' @name poppkmod
#' @title Implement a population pharmacokinetic/pharmacodynamic model.
#'
#' @description Create a `poppkmod` object using an existing population PK model
#' for propofol or remifentanil using patient covariates. Available models for
#' propofol are the Marsh, Schnider, and Eleveld models. Available models for
#' remifentanil are the Minto, Kim, and Eleveld models. Input is
#' a data frame with rows corresponding to individuals and columns
#' recording patient covariates. An ID column is optional, but will be generated
#' as 1:nrow(data) if not supplied. Covariates required by each model are
#'
#'  \strong{Propofol}
#'  \itemize{
#'   \item Marsh: TBW
#'   \item Schnider: (AGE, HGT, TBW, MALE) or (AGE, HGT, LBW)
#'   \item Eleveld: AGE, TBW, HGT, MALE
#'  }
#'
#'  \strong{Remifentanil}
#'  \itemize{
#'   \item Minto: (AGE, HGT, TBW, MALE) or (AGE, HGT, LBW)
#'   \item Kim: (AGE, TBW, BMI, HGT) or (AGE, TBW, FFM)
#'   \item Eleveld: (AGE, MALE, TBW, HGT) or (AGE, MALE, TBW, BMI)
#'  }
#'
#' \strong{Abbreviations}
#' \itemize{
#'  \item TBW = Total body weight (kg)
#'  \item LBW = Lean body weight (kg)
#'  \item FFM = Fat-free mass (kg)
#'  \item AGE = Age (years)
#'  \item HGT = Height (cm)
#'  \item MALE = Male (1/0, TRUE/FALSE)
#'  \item BMI = Body mass index (kg/\eqn{m^2})
#'  }
#'
#' @param data Data frame of patient covariates. ID values, if used, should be in a column labeled "id" or "ID"
#' @param drug "ppf" for propofol or "remi" for remifentanil. Defaults to "ppf".
#' @param model Model name. Options are "marsh", "schnider", or "eleveld" if
#' drug = "ppf", or "minto", "kim", or "eleveld" if drug = "remi".
#' @param sample Logical. Should parameter values be sampled from interindividual distribution (TRUE)
#' or evaluated at point estimates for covariates (FALSE)? Defaults to FALSE.
#' @param PD Logical. If applicable, should the PD component be evaluated for PK-PD models. Defaults to TRUE.
#' @param custom_pkmod_fn Function that generates a custom `pkmod` object. If applicable,
#' this should be a function of all covariates used with names corresponding to `data` argument
#' (see ?pkmod_schnider for an example). If no covariates are used, the `data`
#' argument must contain a column labeled "id" indicating the set of IDs.
#' @param ... Arguments passed on to each pkmod object
#' @return `poppkmod` object
#' @examples
#' data <- data.frame(ID = 1:5, AGE = seq(20,60,by=10), TBW = seq(60,80,by=5),
#' HGT = seq(150,190,by=10), MALE = c(TRUE,TRUE,FALSE,FALSE,FALSE))
#' poppkmod(data, drug = "ppf", model = "eleveld")
#' poppkmod(data, drug = "remi", model = "kim")
#' @export
poppkmod <- function(data, drug = c("ppf","remi"), model = c("marsh","schnider","eleveld","minto","kim"),
                     sample = FALSE, PD = TRUE, custom_pkmod_fn = NULL,...){

  if(!is.null(custom_pkmod_fn)){
    # implement a custom poppkmod
    new_mod <- init_poppkmod(data = data,
                             drug = "Custom",
                             model = "Custom")
    new_mod <- validate_poppkmod(new_mod)

    fn <- custom_pkmod_fn
    fn_nms <- names(formals(fn))
    cov_nms <- setdiff(fn_nms,"...")
    dat_nms <- names(new_mod$data)

    new_mod$pkmods <- with(new_mod,
                           lapply(1:nrow(data), function(i)
                             do.call(fn, c(as.list(data[i,dat_nms %in% cov_nms, drop = FALSE]),
                                           list(...)))))
  } else{

    model <- match.arg(model)
    drug <- match.arg(drug)
    # initialize object
    new_mod <- init_poppkmod(data = data,
                             drug = drug,
                             model = model)

    drug <- attr(new_mod,"drug")
    model <- attr(new_mod,"model")

    # perform validation checks
    new_mod <- validate_poppkmod(new_mod)

    # create list of pkmod objects
    pkmod_fn_nm <- paste0("pkmod_", model)
    if(model == "eleveld") pkmod_fn_nm <- paste0(pkmod_fn_nm,"_",drug)

    fn <- get(pkmod_fn_nm)
    fn_nms <- names(formals(fn))
    cov_nms <- setdiff(fn_nms,"...")
    dat_nms <- names(new_mod$data)

    if("PD" %in% fn_nms){
      new_mod$pkmods <- with(new_mod,
                             lapply(1:nrow(data), function(i)
                               do.call(fn, c(as.list(data[i,dat_nms %in% cov_nms, drop = FALSE]),
                                             PD = PD,
                                             list(...))))
      )
    } else{
      new_mod$pkmods <- with(new_mod,
                             lapply(1:nrow(data), function(i)
                               do.call(fn, c(as.list(data[i,dat_nms %in% cov_nms, drop = FALSE]),
                                             list(...)))))
    }
  }

    if(sample){
      new_mod$pkmods <- lapply(new_mod$pkmods, sample_pkmod)
    }

    if("id" %in% tolower(dat_nms)){
      ix <- which(tolower(dat_nms) == "id")
      new_mod$ids <- new_mod$data[,ix]
    } else{
      id <- 1:nrow(new_mod$data)
    }

    return(new_mod)
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
#' poppkmod(data, drug = "remi", model = "kim")
#' @return Prints description of pkmod
#' @export
print.poppkmod <- function(x, ..., digits = 3){

  drug <- attr(x, "drug")
  model <- attr(x, "model")
  cat("tci poppkmod\n")

  cat("\nModel summary\n")
  drug_type <- ifelse(drug == "ppf","propofol","remifentanil")
  mod_name <- model
  substr(mod_name,1,1) <- toupper(substr(mod_name,1,1))
  cat(paste0(" ", mod_name, " population model for ",drug_type), "\n")
  cat(paste(" Number of individuals:",length(x$pkmods)), "\n")

  cat("\nPK model\n")
  pkstr <- ifelse(x$pkmods[[1]]$ncmpt == 4, "3 compartment-effect", paste(x$pkmods[[1]]$ncmpt,"compartment"))
  cat(paste(" Structural PK model:", pkstr), "\n")
  cat(paste(" PK parameter names:", paste(names(x$pkmods[[1]]$pars_pk), collapse = " ")),"\n")
  tci_type <- ifelse("ke0" %in% tolower(names(x$pkmods[[1]]$pars_pk)), "Effect-site","Plasma")
  cat(paste(" Default TCI targeting:",tci_type),"\n")

  if(!is.null(x$pkmods[[1]]$pdfn)){
    cat("\nPD model\n")
    cat(" Structural PD model: Emax \n")
    cat(paste(" PD parameter names:", paste(names(x$pkmods[[1]]$pars_pd), collapse = " ")),"\n")
  }

  if(!is.null(x$pkmods[[1]]$pars_pd)){
    vals_pd <- summary(t(sapply(x$pkmods, `[[`, "pars_pd")), digits = digits)[c(1,3,6),]
    knitr::kable(vals_pd, format = "rst")
  }

  cat("\nSee poppkmod$pkmods for a list of individual models")
}

# #' Update method for poppkmod
# #'
# #' Update parameters or initial values of a poppkmod object
# #' @param object Object with class poppkmod
# #' @param ... Updated values passed to object.
# #' @return Returns the pkmod object with elements replaced
# #' @examples
# #' # initial poppkmod object
# #' data <- data.frame(ID = 1:5,
# #' AGE = seq(20,60,by=10),
# #' TBW = seq(60,80,by=5),
# #' HGT = seq(150,190,by=10),
# #' MALE = c(TRUE,TRUE,FALSE,FALSE,FALSE))
# #' my_mod <- poppkmod(data, drug = "ppf", model = "eleveld")
# #' # update a subset of parameters and initial values
# #' update(my_mod, pars_pk = c(CL = 4))$pkmods
# #' update(my_mod, pars_pk = c(Q2 = c(1:5)))$pkmods
# #'
# #' update(my_mod, pars_pk = list(CL = 1:5, Q2 = 5:1))$pkmods
# #' @export
# update.poppkmod <- function(object, ...){
#
#   update_elements <- list(...)
#   if(length(update_elements) == 0) return(object)
#   if(is.list(update_elements[[1]])) update_elements <- update_elements[[1]]
#
#   # subset to elements with valid names
#   update_elements <- update_elements[names(update_elements) %in% names(formals(pkmod))]
#
#   nid <- length(object$pkmods)
#   for(i in 1:length(update_elements)){
#     if(!length(update_elements[[i]]) %in% c(1,nid))
#       stop("Arguments must have length 1 or length(data$ID)")
#
#     if(length(update_elements[[i]]) == 1){
#       object$pkmods <- lapply(object$pkmods, update, update_elements[i])
#     } else{
#
#       object$pkmods <- lapply(1:length(object$pkmods), function(j){
#         elm <- list(update_elements[[i]][j])
#         names(elm[[1]]) <- gsub(paste0(j,"$"), "", names(elm[[1]]))
#         names(elm) <- names(update_elements)[i]
#         update(object$pkmods[[j]], elm)
#       })
#     }
#
#   }
#
#   return(validate_poppkmod(object))
# }


#' Summary method for `poppkmod` objects
#'
#' Summarize parameter distribution


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



## -- sim_tci methods ---------------------------------------------------------

#' Print method for sim_tci class
#'
#' Print object with class "sim_tci" created by `simulate_tci()`.
#' @param x Object with class "sim_tci" created by `simulate_tci()`
#' @param ... Other arguments. Not currently used.
#' @return Prints a description of the simulation.
#' @examples
#' data <- data.frame(ID = 1:3, AGE = c(20,30,40), TBW = c(60,70,80),
#' HGT = c(150,160,170), MALE = c(TRUE,FALSE,TRUE))
#' pkmod_prior <- poppkmod(data, drug = "ppf", model = "eleveld")
#' pkmod_true  <- poppkmod(data, drug = "ppf", model = "eleveld", sample = TRUE)
#' obs_tms <- seq(1/6,10,1/6)
#' target_vals = c(75,60,50,50)
#' target_tms = c(0,3,6,10)
#'
#' # open-loop simulation (without update_tms)
#' sim_ol <- simulate_tci(pkmod_prior, pkmod_true, target_vals, target_tms, obs_tms,
#' seed = 200)
#'
#' # closed-loop simulation (with update_tms)
#' \dontrun{
#' sim_cl <- simulate_tci(pkmod_prior, pkmod_true, target_vals, target_tms, obs_tms,
#' update_tms = c(2,4,6,8), seed = 200)
#' }
#' @export
print.sim_tci <- function(x, ...){
  cat("tci sim_tci\n")

  cat("\n", "TCI Simulation", "\n", sep = "")
  nid <- ifelse("id" %in% colnames(x$resp), length(unique(x$resp$id)), 1)
  cat("  N individuals: ", nid, "\n", sep = "")
  cat("  Control type: ", x$control, "\n", sep = "")

  resp_type <- ifelse("pdresp" %in% names(x$resp), "PD","PK")
  cat("  Response type: ", resp_type, "\n", sep = "")

  tm_range <- round(range(x$resp[,"time"]),1)
  cat("  Simulation duration: [",tm_range[1],",",tm_range[2],"]" , "\n", sep = "")

  if(x$control == "closed-loop")
    cat("  Update times: ", paste(x$update_tms, collapse = ", "), "\n", sep = "")

  obs_tms <- round(x$obs_tms,1)
  if(length(obs_tms)>10) {
    obs_tms <- c(obs_tms[1:4], obs_tms[(length(obs_tms)-2):length(obs_tms)])
    cat("  Observation times: ", paste(obs_tms[1:4], collapse = ", "),", ... , ",
        paste(obs_tms[5:7], collapse = ", "), sep = "")
  } else{
    cat("  Observation times: ", paste(round(x$obs_tms,1), collapse = ", "), sep = "")
  }
}




#' Plot method for sim_tci class
#'
#' Plot object with class "sim_tci" created by `simulate_tci()`.
#' @param x Object with class "sim_tci" created by `simulate_tci()`
#' @param ... Other arguments. Not currently used.
#' @param yvar Response variable. Options are concentrations ("c1","c2",...) or
#' "pdresp" for a PD response. Only one variable may be plotted at a time.
#' @param id Subset of IDs to plot. Will default to all if unspecified. Can be
#' displayed in separate plots via `wrap_id` argument.
#' @param type Type of response to plot. Options are "prior", "true", or "posterior"
#' if closed-loop control was used.
#' @param show_inf Logical. Display infusion rates alongside response values.
#' @param show_data Logical. Display simulated data values in addition to responses.
#' @param show_updates Logical, for closed-loop only. Show update times along x-axis.
#' @param wrap_id Logical. Separate plots by ID value.
#' @return Plots simulation results
#' @import ggplot2
#' @importFrom reshape melt
#' @examples
#' data <- data.frame(ID = 1:2, AGE = c(30,40), TBW = c(70,80),
#' HGT = c(160,170), MALE = c(FALSE,TRUE))
#' pkmod_prior <- poppkmod(data, drug = "ppf", model = "eleveld")
#' pkmod_true  <- poppkmod(data, drug = "ppf", model = "eleveld", sample = TRUE)
#' obs_tms <- seq(1/6,10,1/6)
#' target_vals = c(75,60,50,50)
#' target_tms = c(0,3,6,10)
#'
#' # open-loop simulation (without update_tms)
#' sim_ol <- simulate_tci(pkmod_prior, pkmod_true, target_vals, target_tms, obs_tms,
#' seed = 200)
#' plot(sim_ol, id = 1, type = "true")
#' plot(sim_ol, yvar = "c4", type = "true")
#' plot(sim_ol, yvar = "c4", type = "true", wrap_id = TRUE, show_inf = TRUE)
#'
#' # closed-loop simulation (with update_tms)
#' \dontrun{
#' sim_cl <- simulate_tci(pkmod_prior, pkmod_true, target_vals, target_tms, obs_tms,
#' update_tms = c(2,4,6), seed = 200)
#' plot(sim_cl, type = "posterior", id = 1, show_inf = TRUE)
#' plot(sim_cl, type = "posterior", wrap_id = TRUE, show_data = TRUE)
#' plot(sim_cl, yvar = "c4", wrap_id = TRUE)
#' }
#' @export
plot.sim_tci <- function(x, ..., yvar = NULL, id = NULL,
                         type = c("true","prior","posterior"),
                         show_inf = FALSE, show_data = FALSE,
                         show_updates = FALSE, wrap_id = FALSE){

  value = NA; variable = NA
  type = match.arg(type)

  if(length(yvar)>1) stop("Only 1 yvar value is allowed")
  if(type == "posterior" & x$control == "open-loop")
    stop("posterior type is not available for open-loop simulations")

  resp_options <- names(x$resp)[names(x$resp) %in% c(paste0("c",1:9),"pdresp")]
  if(!is.null(yvar)){
    if(!yvar %in% resp_options)
      stop(paste0("yvar must be one of ", paste(resp_options, collapse = ", ")))
  }

  if(is.null(yvar) & "pdresp" %in% names(x$resp)) yvar <- "pdresp"
  if(is.null(yvar) & !("pdresp" %in% names(x$resp))) yvar <- "c1"

  if(!"id" %in% colnames(x$resp)){
    x$resp <- cbind(id = 1, x$resp)
    x$inf <- cbind(id = 1, x$inf)
    x$obs <- cbind(id = 1, x$obs)
  }
  if(is.null(id)) id <- unique(x$resp[,"id"])
  target_var <- ifelse(yvar=="pdresp", "pdt","Ct")

  resp <- melt(x$resp, id.vars = c("id","time","type"))
  resp <- resp[resp$variable == yvar,]
  resp <- resp[resp$type == type,]
  resp$type <- "Response"
  resp$id <- as.factor(resp$id)

  # infusion data frame
  inf <- data.frame(id = as.factor(rep(x$inf[,"id"],times = 2)),
                    time = rep(x$inf[,"begin"],times = 2),
                    variable = c(x$inf[,target_var], x$inf[,"inf_rate"]),
                    type = factor(rep(c("Response","Infusion Rate"), each = nrow(x$inf)),
                                  levels = c("Infusion Rate","Response"), ordered = TRUE))

  if(!show_inf){
    inf <- inf[inf$type == "Response",]
  }

  # observation data frame
  obs <- as.data.frame(x$obs)
  obs$type = "Response"

  # subset to selected ids
  resp <- resp[resp$id %in% id,]
  obs <- obs[obs$id %in% id,]
  inf <- inf[inf$id %in% id,]

  p <- ggplot(resp, aes(x = time, y = value, color = id, linetype = id)) +
      geom_line() +
      geom_step(data = inf, aes(x = time, y = variable), color = "black") +
      facet_wrap(~type, nrow = 2, scales = "free")

  if(show_data) p <- p + geom_point(data = obs, aes(x = time, y = obs), alpha = 0.5)

  if(show_updates & x$control == "closed-loop"){
    update_df <- as.data.frame(expand.grid(id = id,
                                           time = x$update_tms))
    update_df$type = "Response"
    p <- p + geom_vline(data = update_df, aes(xintercept = time), linetype = "dotted")
  }

  if(wrap_id) p <- p + facet_grid(type~id, scales = "free")

  p
}

