# --------------------------------------------------------------------------------------------------------------------------------
# - PK-PD model helper functions and methods -------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------

# dosing schedule
create_intvl <- function(dose, inittm = 0){
  # b <- cut2(dose$time, breaks = c(inittm,dose$time), include.lowest = TRUE, right = FALSE)
  b <- cut2(dose$time +inittm, breaks = c(inittm,dose$time+inittm), include.lowest = TRUE, right = FALSE)
  ss <- t(sapply(stringr::str_extract_all(levels(b),"-?[0-9.]+"), as.numeric))
  setNames(data.frame(levels(b), dose$infrt, ss), c("intvl","infrt","begin","end"))
}
dose <- data.frame(time = c(0.5,4,4.5,10), infrt = c(100,0,100,0))
create_intvl(dose)

# predict method to apply pk model piecewise to infusion schedule
predict.pkmod <- function(pkmod, inf, tms = NULL, dt = 1/6, len_out = NULL, return_init = FALSE, remove_bounds = TRUE, ...){

  pred <- vector("list", length(inf$infrt))
  init <- vector("list", length(inf$infrt)+1)

  # Times to evaluate concentrations at. Defaults to a sequence of values at intervals of dt.
  if(!is.null(tms)){
    b <- unique(as.numeric(unlist(stringr::str_extract_all(inf$intvl,"-?[0-9.]+"))))
    tms_all <- unique(sort(c(b,tms)))
    # tms_eval <- split(tms_all, cut(tms_all, breaks = b, right = T))
    tms_eval <- split(tms_all, findInterval(tms_all, b, rightmost.closed = T, left.open = T))
  } else{
    if(is.null(len_out)) tms_eval <- mapply(seqby, inf$begin+dt, inf$end, by = dt, SIMPLIFY = F)
    else tms_eval <- mapply(seq, inf$begin+dt, inf$end, length.out = len_out)
  }

  # Pass on initial concentrations to first element of init. Use values if specified, else defaults.
  dot.args <- list(...)
  if("init" %in% names(dot.args)){
    init[[1]] <- unlist(dot.args$init)
    dot.args$init <- NULL
  } else {
    init[[1]] <- eval(formals(pkmod)$init)
  }

  # Predict concentrations and store initial values.
  for(i in 1:nrow(inf)){
    pred[[i]] <- do.call("pkmod", c(list(tm = tms_eval[[i]], kR = inf$infrt[i], init = init[[i]], inittm = inf$begin[i]), dot.args))
    init[[i+1]] <- tail_vec(pred[[i]]) # extract last element OR column
  }
  pred <- lapply(pred, function(x) ifelse(x < 0, 0, x))

  # Return predicted concentrations
  if(is.null(dim(pred[[1]]))) {
    predtms <- cbind(unique(unlist(tms_eval)),  do.call("c", pred))
  } else{
    predtms <- cbind(unique(unlist(tms_eval)),  t(do.call("cbind", pred)))
  }

  # Add on t=0 concentrations
  if(return_init) predtms <- rbind(c(inf$begin[1], init[[1]]), predtms)

  # remove transition concentrations
  if(!is.null(tms) & remove_bounds) predtms <- matrix(predtms[which(predtms[,1] %in% tms),], nrow = length(tms), byrow = F)

  colnames(predtms) <- c("time",paste0("c",1:length(init[[1]])))
  return(predtms)
}

# when calling generic "predict" on an object of class "pkmod", use method "predict.pkmod"
.S3method(generic = "predict", class = "pkmod", method = predict.pkmod)


# Note: ... arguments SHOULD be passed to plotting function instead of predict. Will fix later.
plot.pkmod <- function(pkmod, inf, npts = 1000, title = NULL, ...){
  # set dt based on range between points
  dt <- diff(range(inf$begin, inf$end)) / npts
  # predict concentrations
  con <- data.frame(predict(pkmod, inf, dt = dt, return_init = T, ...))

  ggplot(melt(con, id = "time"), aes(x = time, y = value, linetype = variable, color = variable)) +
    geom_line() +
    labs(y = "Concentration", x = "Time", color = "Compartment", linetype = "Compartment", title = title)
}

.S3method(generic = "plot", class = "pkmod", method = plot.pkmod)



# Plot method for PD models. User can provide a series of effect-site concentrations and a PD model OR an infusion schedule with a PK-PD model
plot.pdmod <- function(pdmod, pkmod, inf, pars_pd, pars_pk, npts = 1000, plot_pk = TRUE, title = NULL, ecmpt = NULL, ...){

  # set dt based on range between points
  dt <- diff(range(inf$begin, inf$end)) / npts
  # predict concentrations
  con <- data.frame(predict(pkmod, inf, dt = dt, return_init = T, pars = pars_pk, ...))

  # effect site comparment
  if(is.null(ecmpt))
    ecmpt <- length(eval(formals(pkmod)$init))

  # predict PD values
  pd <- data.frame(time = con$time, pdresp = pdmod(con[,ecmpt], pars_pd))

  if(plot_pk){
    p1 <- ggplot(melt(con, id = "time"), aes(x = time, y = value, linetype = variable, color = variable)) +
      geom_line() +
      labs(y = "Concentration", x = "Time", color = "Compartment", linetype = "Compartment") +
      theme(legend.position="bottom")
  } else{
    p1 <- NULL
  }

  p2 <- ggplot(pd, aes(x = time, y = pdresp)) +
    geom_line() +
    labs(y = "Response", x = "Time")

  if(!is.null(p1)){
    grid.arrange(p2, p1, nrow = 2, top = title)
  } else{
    grid.arrange(p2, top = title)
  }
}
.S3method("plot","pdmod",plot.pdmod)



#' Function to extend TCI grid to a set of PD targets
tci_pd <- function(pdresp, tms, pkmod, pdmod, pars_pk, pars_pd, pdinv, ecmpt = NULL, ...){
  Ct <- pdinv(pdresp, pars_pd)
  con <- iterate_tci_grid(Ct = Ct, tms = tms, pkmod = pkmod, pars = pars_pk, ...)
  if(is.null(ecmpt))
    ecmpt <- length(eval(formals(pkmod)$init))
  con$pdt <- pdmod(con$Ct, pars_pd)
  con$pdresp <- pdmod(con[,paste0("c",ecmpt,"_start")], pars_pd)
  return(con)
}


#' Plotting method for tciinf objects. Will work for outputs from "iterate_tci_grid" or "tci_pd".
plot.tciinf <- function(tciinf, title = NULL){

  # data frame for PK plot
  tciinfm <- melt(tciinf[,c("infrt","dt","begin","Ct",grep("\\_start",names(tciinf), value = T))], id.vars = c("infrt","dt","begin","Ct"))
  tciinfm$variable <- gsub("\\_start|c", "", tciinfm$variable)
  names(tciinfm)[names(tciinfm) == "variable"] <- "Compartment"

  ppk <- ggplot(tciinf[,c("begin","Ct")], aes(x = begin, y = Ct, fill = "")) +
    geom_line(data = tciinf[,c("begin","Ct")], color = "black", size = 0.8, linetype = "dashed") +
    guides(fill=guide_legend(title="Target"), linetype = "dashed")  +
    geom_line(data = tciinfm, mapping = aes(x = begin, y = value, color = Compartment), size = 1, show.legend = TRUE) +
    ylim(range(tciinfm$value)) +
    labs(x = "Time", y = "Concentration") +
    theme(legend.position="bottom")

  if("pdresp" %in% names(tciinf)){
    ppd <- ggplot(tciinf, aes(x = begin, y = pdt)) +
      geom_line(data = tciinf, color = "black", size = 0.8, linetype = "dashed") +
      # guides(fill=guide_legend(title="Target effect", linetype = "dashed"))  +
      geom_line(data = tciinf,
                mapping = aes(x = begin, y = pdresp), size = 1, show.legend = TRUE) +
      ylim(c(0,100)) +
      labs(x = "Time", y = "PD response") +
      theme(legend.position="bottom")
  }

  if("pdresp" %in% names(tciinf)){
    grid.arrange(ppd, ppk, nrow = 2, top = title)
  } else{
    grid.arrange(ppk, top = title)
  }

}
.S3method("plot","tciinf",plot.tciinf)


#' Function to place restriction on gamma and E50 parameters of target sigmoid
#' such that it passes through point (tfinal, BISfinal+eps)
restrict_sigmoid <- function(t50, tfinal =10, eps = 1, BIS0 = 100, BISfinal = 50){
  gamma <- log((BIS0-BISfinal)/eps - 1, base = tfinal/t50)
  c(c50 = t50, gamma = gamma, E0 = BIS0, Emx = BIS0 - BISfinal)
}


