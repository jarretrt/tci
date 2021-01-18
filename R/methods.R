# --------------------------------------------------------------------------------------------------------------------------------
# - PK-PD model methods ----------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------


#' Predict concentrations from a pkmod object
#'
#' predict method to apply pk model piecewise to infusion schedule
#' @param object An object with class pkmod.
#' @param ... Arguments passed on to pkmod
#' @param inf An infusion schedule object with columns "begin","end","infrt".
#' @param tms Times to evaluate predictions at. Will default to a sequence
#' spanning the infusions at intervals of dtm.
#' @param dtm Interval used for prediction if argument tms is unspecified.
#' @param return_init Logical indicating if concentrations at time 0 should
#' be returned. Defaults to FALSE.
#' @param remove_bounds Logical, indicating if concentrations calculated at
#' changes in infusion rates should be returned if not included in prediction
#' times. Defaults to TRUE, so that only concentrations at specified times
#' are returned.
#'
#' @export
predict.pkmod <- function(object, ..., inf, tms = NULL, dtm = 1/6, return_init = FALSE,
                          remove_bounds = TRUE){

  fptol = 1e-10
  tm_digits = 7

  if(!all(c("infrt","begin","end") %in% colnames(inf)))
    stop("inf must include 'infrt','begin','end' as column names")

  dot.args <- list(...)

  if(!("init" %in% names(formals(object))))
    stop("object must contain argument 'init'")

  if(!("pars" %in% names(dot.args)) &
     is.symbol(formals(object)$pars))
    stop("PK parameters must be passed as 'pars' within predict or set
         as defaults in PK model object")

  # Times to evaluate concentrations at. Defaults to a sequence of values at intervals of dtm.
  if(!is.null(tms)){
    # round times - this is needed to prevent errors associated with rounding numeric values
    tms <- round(tms, tm_digits)
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
  } else{
    # if times are not provided, predict across a grid of points
    tms_eval <- mapply(seq, inf[,"begin"]+dtm, inf[,"end"], by = dtm,
                       SIMPLIFY = FALSE)
    tms_all <- unlist(tms_eval)
  }

  init <- vector("list", nrow(inf)+1)

  # Pass on initial concentrations to first element of init. Use values if specified, else defaults.
  if("init" %in% names(dot.args)){
    init[[1]] <- unlist(dot.args$init)
    dot.args$init <- NULL
  } else {
    init[[1]] <- eval(formals(object)$init)
  }

  if("pars" %in% names(dot.args)){
    pars <- unlist(dot.args$pars)
    dot.args$pars <- NULL
  } else {
    pars <- eval(formals(object)$pars)
  }

  # get indexes of times and initialize matrix for predictions
  tm_ix <- lapply(tms_eval, function(x) match(x,tms_all))
  ncmpt <- length(eval(formals(object)$init))
  pred <- matrix(NA, nrow = ncmpt, ncol = length(tms_all))

  # Predict concentrations and store initial values.
  for(i in 1:nrow(inf)){
    pred[,tm_ix[[i]]] <- do.call("object", c(list(tm = tms_eval[[i]],
                                                 kR = inf[i,"infrt"],
                                                 pars = pars,
                                                 init = init[[i]],
                                                 inittm = inf[i,"begin"]),
                                            dot.args))
    init[[i+1]] <- pred[,length(unlist(tm_ix[1:i]))]
  }

  # Replace any negative values
  pred[pred<0] <- 0

  # Return predicted concentrations
  if(dim(pred)[1] == 1) {
    predtms <- cbind(unique(unlist(tms_eval)), c(pred))
  } else{
    predtms <- cbind(unique(unlist(tms_eval)), t.default(pred))
  }

  # Add on t=0 concentrations
  if(return_init) predtms <- rbind(c(inf[1,"begin"], init[[1]]), predtms)

  # # remove transition concentrations
  if(!is.null(tms) & remove_bounds){
    predtms <- matrix(predtms[which(predtms[,1] %in% tms),],
                      nrow = length(tms), byrow = FALSE)
  }

  colnames(predtms) <- c("time",paste0("c",1:length(init[[1]])))
  return(predtms)
}


#' Plot object with class 'pkmod'
#'
#' Will show predicted concentrations in compartments associated with an infusion schedule.
#' @param x An object with class pkmod.
#' @param ... Arguments passed on to predict.pkmod
#' @param inf An infusion schedule object with columns "begin","end","infrt".
#' @param npts Number of points used to evaluate predicted concentrations.
#'
#' @rdname plot
#' @export
plot.pkmod <- function(x, ..., inf, npts = 1000, title = NULL){

  value <- variable <- NULL

  # set dtm based on range between points
  dtm <- diff(range(inf[,"begin"], inf[,"end"])) / npts
  # predict concentrations
  con <- data.frame(predict(x, inf = inf, dtm = dtm, return_init = TRUE, ...))

  colnames(con) <- gsub("^c", "Cmpt", colnames(con))
  ggplot2::ggplot(reshape::melt(con, id = "time"),
                  ggplot2::aes(x = time,
                               y = value,
                               linetype = variable,
                               color = variable)) +
    ggplot2::geom_line() +
    ggplot2::labs(y = "Concentration",
                  x = "Time",
                  color = "Compartment",
                  linetype = "Compartment",
                  title = title) +
    ggplot2::scale_color_manual(values = unname(pal))
}


#' Plot method for PD models.
#'
#' User can provide a series of effect-site concentrations and a PD model or
#'  an infusion schedule with a PK-PD model.
#'
#' @param x An object with class pdmod.
#' @param ... Arguments passed on to method predict.pkmod.
#' @param pkmod PK model
#' @param inf An infusion schedule object with columns "begin","end","infrt".
#' @param pars_pd Parameters used by pdmod.
#' @param pars_pk Parameters used by pkmod.
#' @param npts Number of points used to evaluate predicted concentrations.
#' @param plot_pk Logical. Should PK concentrations be plotted alongside
#' the PD response. Defaults to TRUE.
#' @param title Title of plot
#' @param ecmpt Effect-site compartment number. Defaults to the last
#' compartment concentration returned by pkmod.
#'
#' @rdname plot
#' @export
plot.pdmod <- function(x, ..., pkmod, inf, pars_pd, pars_pk = NULL, npts = 1000,
                       plot_pk = TRUE, title = NULL, ecmpt = NULL){

  value <- variable <- pdresp <- NULL

  # set dtm based on range between points
  dtm <- diff(range(inf[,"begin"], inf[,"end"])) / npts
  # predict concentrations
  if(is.null(pars_pk)){
    con <- data.frame(predict(pkmod, inf = inf, dtm = dtm, return_init = TRUE, ...))
  } else{
    con <- data.frame(predict(pkmod, inf = inf, dtm = dtm, return_init = TRUE, pars = pars_pk, ...))
  }

  # effect site compartment
  if(is.null(ecmpt))
    ecmpt <- length(eval(formals(pkmod)$init))

  # predict PD values
  pd <- data.frame(time = con$time, pdresp = x(con[,ecmpt], pars_pd))

  if(plot_pk){
    p1 <- ggplot2::ggplot(reshape::melt(con, id = "time"),
                          ggplot2::aes(x = time,
                                       y = value,
                                       linetype = variable,
                                       color = variable)) +
      ggplot2::geom_line() +
      ggplot2::labs(y = "Concentration",
                    x = "Time",
                    color = "Compartment",
                    linetype = "Compartment") +
      ggplot2::theme(legend.position="bottom") +
      ggplot2::scale_color_manual(values = unname(pal))
  } else{
    p1 <- NULL
  }

  p2 <- ggplot2::ggplot(pd, ggplot2::aes(x = time, y = pdresp)) +
    ggplot2::geom_line(color = unname(pal[5])) +
    ggplot2::labs(y = "Response", x = "Time") +
    ggplot2::lims(y = c(0,100))

  if(!is.null(p1)){
    gridExtra::grid.arrange(p2, p1, nrow = 2, top = title)
  } else{
    gridExtra::grid.arrange(p2, top = title)
  }
}


#' Plotting method for tciinf objects. Will work for outputs from "iterate_tci_grid" or "tci_pd".
#' Note: the tci targets function is left-continous. Consequently, it plots
#' final target value plotted will be a repetition of the preceeding value
#'
#' @param x Object with class "tciinf" created by functions
#' `iterate_tci_grid` or `tci_pd`
#' @param ... \dots
#' @param title Title of plot.
#' @param display Logical. Should plots be printed or returned as an arrangeGrob object?
#'
#' @rdname plot
#' @export
plot.tciinf <- function(x, ..., title = NULL, display = TRUE){

  begin <- value <- variable <- NULL

  tciinf <- as.data.frame(x)
  # move end of last row to start column for plotting
  tciinf <- rbind(tciinf, NA)
  ncpt <- length(grep("c[0-9]_start",names(tciinf)))
  tciinf[nrow(tciinf),grep("start|begin|Ct|pdt", names(tciinf))] <-
    tciinf[nrow(tciinf)-1,grep("end|Ct|pdt", names(tciinf))]

  # data frame for PK plot
  tciinfm <- reshape::melt(tciinf[,c("begin","Ct", paste0("c",1:ncpt,"_start"))],
                           id.vars = "begin")
  tciinfm$variable <- gsub("_start","",tciinfm$variable)
  tciinfm$variable <- gsub("c","Cmpt",tciinfm$variable)
  tciinfm$variable <- gsub("Ct","Target",tciinfm$variable)
  tciinfm$variable <- factor(tciinfm$variable,
                             levels = c("Target",grep("Cmpt",unique(tciinfm$variable),
                                                      value = TRUE)))

  ppk <- ggplot2::ggplot(tciinfm) +
    ggplot2::geom_step(data = tciinfm[tciinfm$variable == "Target",],
                       ggplot2::aes(x = begin, y = value, color = variable, linetype = variable),
              size = 1) +
    ggplot2::geom_line(data = tciinfm[tciinfm$variable != "Target",],
                       ggplot2::aes(x = begin, y = value, color = variable, linetype = variable),
              size = 1) +
    ggplot2::scale_linetype_manual("", values = c(1:ncpt,1),
                          labels = c(paste0("Cmpt",1:ncpt),"Target")) +
    ggplot2::scale_color_manual("", values = c(unname(pal[2:(ncpt+1)]),"black"),
                       labels = c(paste0("Cmpt",1:ncpt),"Target")) +
    ggplot2::labs(x = "Time", y = "Concentration", color = "", linetype = "") +
    ggplot2::theme(legend.position="bottom")


  if("pdresp_start" %in% names(tciinf)){
    tciinfm2 <- reshape::melt(tciinf[,c("begin","pdt","pdresp_start")],id.vars = "begin")

    ppd <- ggplot2::ggplot(tciinfm2,
                           ggplot2::aes(x = begin, y = value)) +
      ggplot2::geom_step(data = tciinfm2[tciinfm2$variable == "pdt",],
                ggplot2::aes(color = "col1", linetype = "solid"), size = 1) +
      ggplot2::geom_line(data = tciinfm2[tciinfm2$variable != "pdt",],
                         ggplot2::aes(color = "col2", linetype = "solid"), size = 1)+
      ggplot2::scale_linetype_manual("", values = c("solid","solid"), labels = c("PD Target","PD Response")) +
      ggplot2::scale_color_manual("", values = unname(pal[c(1,5)]), labels = c("PD Target","PD Response")) +
      ggplot2::ylim(c(0,100)) +
      ggplot2::labs(x = "Time", y = "PD response") +
      ggplot2::theme(legend.position="bottom")
  }

  if("pdresp_start" %in% names(tciinf)){
    # gridExtra::grid.arrange(ppd, ppk, nrow = 2, top = title)
    gb <- gridExtra::arrangeGrob(ppd, ppk, nrow = 2, top = title)
  } else{
    # gridExtra::grid.arrange(ppk, top = title)
    gb <- gridExtra::arrangeGrob(ppk, top = title)
  }

  if(display){
    plot(gb)
  } else{
    gb
  }
}


#' Plot method for datasim objects
#'
#' @param x An object with class datasim, created by function `gen_data`.
#' @param pars_prior Named vector of prior PK or PK-PD parameters
#' @param pars_post  Named vector of posterior PK or PK-PD parameters
#' @param pk_ix Indicies of parameter vector(s) corresponding to PK parameters
#' @param pd_ix Indicies of parameter vector(s) corresponding to PD parameters
#' @param ... \dots
#'
#' @rdname plot
#' @importFrom ggplot2 aes
#' @export
plot.datasim <- function(x, ..., pars_prior = NULL, pars_post = NULL, pk_ix = NULL, pd_ix = NULL){

  if((!is.null(pars_prior) | !is.null(pars_post)) & (is.null(pk_ix) | is.null(pd_ix)))
    stop("pk_ix and pd_ix must be specified if pars_prior or pars_post are non-null")

  # initialize for CRAN check -- will be created by other functions
  c1 <- variable <- cobs <- pdp <- pdobs <- NULL

  datasim <- x
  datasim$sim <- as.data.frame(datasim$sim)
  datasim$inf <- as.data.frame(datasim$inf)
  r <- range(datasim$inf[,c("begin","end")])
  tms <- seq(r[1], r[2], length.out = 1000)

  # predict concentrations at true parameter values
  cp <- data.frame(predict(datasim$pkmod,
                           inf = datasim$inf,
                           tms = tms,
                           pars = datasim$pars_pk0,
                           init = datasim$init))

  if(!is.null(datasim$pdmod)){
    cp$pdp <- datasim$pdmod(ce = cp[,paste0("c",datasim$ecmpt)],
                            pars = datasim$pars_pd0)
  }
  cp$variable <- as.factor("Observed")

  # predict concentrations at prior parameters, if specified
  if(!is.null(pars_prior)){
    cp_prior <- data.frame(predict(datasim$pkmod,
                                   inf = datasim$inf,
                                   tms = tms,
                                   pars = pars_prior[pk_ix],
                                   init = datasim$init))

    if(!is.null(datasim$pdmod)){
      cp_prior$pdp <- datasim$pdmod(ce = cp_prior[,paste0("c",datasim$ecmpt)],
                                    pars = pars_prior[pd_ix])
    }

    cp_prior$variable <- as.factor("Prior")

  } else{
    cp_prior <- NULL
  }

  # predict concentrations at posterior parameters, if specified
  if(!is.null(pars_post)){
    cp_post <- data.frame(predict(datasim$pkmod,
                                  inf = datasim$inf,
                                  tms = tms,
                                  pars = pars_post[pk_ix],
                                  init = datasim$init))

    if(!is.null(datasim$pdmod)){
      cp_post$pdp <- datasim$pdmod(ce = cp_post[,paste0("c",datasim$ecmpt)],
                                   pars = pars_post[pd_ix])
    }

    cp_post$variable <- as.factor("Posterior")

  } else{
    cp_post <- NULL
  }

  df <- rbind(cp,cp_prior,cp_post)

  # plot for pk simulations
  if(is.null(datasim$pdmod)){
    out <- ggplot2::ggplot(df,
                           ggplot2::aes(x = time,
                                        y = c1,
                                        color = variable,
                                        linetype = variable)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(data = datasim$sim,
                          ggplot2::aes(x = time, y = cobs),
                          shape = 16, col = pal["navy"],
                          inherit.aes = FALSE, alpha = 1) +
      ggplot2::scale_color_manual(values = unname(pal[c(1,4)])) +
      ggplot2::labs(x = "Time (min)", y = "Concentration", color = "", linetype = "")

  } else{
    # plot for pd simulations

    # add targets to data frame
    tmp <- datasim$inf[,c("begin","pdt")]
    names(tmp) <- c("time","pdp")
    tmp$variable <- "Target"
    for(nm in setdiff(names(df), names(tmp))){
      tmp[,nm] <- NA
    }
    df <- rbind(df, tmp)
    lv <- length(levels(df$variable))

    out <- ggplot2::ggplot(df,
                           ggplot2::aes(x = time,
                                        y = pdp,
                                        color = variable,
                                        linetype = variable)) +
      ggplot2::geom_line(ggplot2::aes(x = time, y = pdp, color = variable, linetype = variable),
                         size = 1) +
      ggplot2::geom_point(data = as.data.frame(datasim$sim),
                          ggplot2::aes(x = time, y = pdobs),
                          color = unname(pal["darkgrey"]),
                          inherit.aes = FALSE, alpha = 0.5) +
      ggplot2::scale_linetype_manual("", values = c(1:(lv-1),1)) +
      ggplot2::scale_color_manual("", values = unname(c(pal[2:lv],pal[1]))) +
      ggplot2::labs(x = "Time", y = "PD response") +
      ggplot2::theme(legend.position="bottom")
  }


  return(out)
}


#' Plot method for bayessim objects
#'
#' Plot output returned by "bayes_control" function.
#'
#' @param x Object returned from "bayes_control" function
#' @param ... \dots
#'
#' @rdname plot
#' @importFrom ggplot2 aes
#' @export
plot.bayessim <- function(x, ...){

  pars_prior <- x$prior$pars_pkpd
  pars_post <- rep(NA, length(pars_prior))
  pars_post[x$prior$fixed_ix] <- pars_prior[x$prior$fixed_ix]
  pars_post[is.na(pars_post)] <- exp(x$lpr[nrow(x$lpr),-ncol(x$lpr)])
  names(pars_post) <- names(pars_prior)

  # plot.datasim method
  plot(x$dat,
       pars_prior = pars_prior,
       pars_post = pars_post,
       pk_ix = x$prior$pk_ix,
       pd_ix = x$prior$pd_ix)
}
