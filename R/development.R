# Functions in development


#' A function to return the basic solution for a one compartment PK model.
#'
#' @param kR A numeric value specifying an infusion rate.
#' @param pars A named vector of PK parameters with names "k10" "v" or alternately "CL" and "v".
#' @param init Optional initial concentration in central compartment.
#' @return A function that calculates compartment concentrations as a function of time for the specified PK model.
#' @export
pk_basic_solution_1cpt <- function(kR, pars, init = 0){

  # parameters should be named list with values "k10","v" OR "CL","v"
  list2env(as.list(pars), envir = environment())

  if("CL" %in% ls()){
    k10 <- CL / v
  }

  init_amt <- init * v1
  a_1 <- function(t) kR/k10*(1-exp(-t*k10))+init_amt*exp(-t*k10)
  c_1 <- function(t) a_1(t)/v1

  out <- list(c_1=c_1)
  attributes(out) <- c(attributes(out), list(infusion = kR, param = pars, initial = init))
  return(out)
}
attributes(pk_basic_solution_1cpt) <- list(ncompartments = 1)

#' @examples
#' pars_1cpt <- exp(c(k10=-2.496,v1=2.348))
#' pk_basic_solution_1cpt(kR = 10, pars = pars_1cpt)



custom_pk1cpt <- function(kR, pars, init){

}


#' A function to return the basic solution for a two compartment PK model.
#'
#' @param kR A numeric value specifying an infusion rate.
#' @param pars A named list of PK parameters with names ("k10","k12","k21","v1","v2) or alternately ("Cl","Q","v1","v2).
#' Elimination from the second compartment can be permitted by including "k20" within the list of parameters.
#' @param init Optional initial concentration in central compartment. Default is to assume no existing medication in either compartment.
#' @return A function that calculates compartment concentrations as a function of time for the specified PK model.
#' @export
pk_basic_solution_2cpt <- function(kR, pars, init=c(0,0)){

  list2env(as.list(pars), envir = environment())

  if("CL" %in% ls()){
    k10 <- CL / v1
    k12 <- Q / v1
    k21 <- Q / v2
  }
  if(!("k20" %in% ls())){
    k20 <- 0
  }
  E1 <- k10+k12
  E2 <- k21+k20

  #calculate hybrid rate constants
  lambda1 = 0.5*((E1+E2)+sqrt((E1+E2)^2-4*(E1*E2-k12*k21)))
  lambda2 = 0.5*((E1+E2)-sqrt((E1+E2)^2-4*(E1*E2-k12*k21)))

  A1last <- init[1]*v1
  A2last <- init[2]*v2
  Doserate <- kR

  a_1 <- function(t) {
    A1term1 = (((A1last*E2+Doserate+A2last*k21)-A1last*lambda1)*exp(-t*lambda1)-((A1last*E2+Doserate+A2last*k21)-A1last*lambda2)*exp(-t*lambda2))/(lambda2-lambda1)
    A1term2 = Doserate*E2*(1/(lambda1*lambda2)+exp(-t*lambda1)/(lambda1*(lambda1-lambda2))-exp(-t*lambda2)/(lambda2*(lambda1-lambda2)))
    return(A1term1+A1term2)
  }

  a_2 <- function(t) {
    A2term1 = (((A2last*E1+A1last*k12)-A2last*lambda1)*exp(-t*lambda1)-((A2last*E1+A1last*k12)-A2last*lambda2)*exp(-t*lambda2))/(lambda2-lambda1)
    A2term2 = Doserate*k12*(1/(lambda1*lambda2)+exp(-t*lambda1)/(lambda1*(lambda1-lambda2))-exp(-t*lambda2)/(lambda2*(lambda1-lambda2)))
    return(A2term1+A2term2)
  }

  c_1 <- function(t) a_1(t)/v1
  c_2 <- function(t) a_2(t)/v2

  out <- list(c_1=c_1, c_2=c_2)
  attributes(out) <- list(infusion = kR, param = pars, initial = init)

  return(out)
}
attributes(pk_basic_solution_2cpt) <- list(ncompartments = 2)


#' A function to return the basic solution for a three compartment PK model with a metabolite compartment.
#'
#' @param kR A numeric value specifying an infusion rate.
#' @param pars A named list of PK parameters with names ("k10","k12","k21","k13","k31" "v1","v2","v3","ke0").
#' Elimination from the second or third compartments can be permitted by including "k20" or "k30", respectively within the list of parameters.
#' @param init Optional initial concentration in central compartment. Default is to assume no existing medication in any compartment.
#' @return A function that calculates compartment concentrations as a function of time for the specified PK model.
#' @export
pk_basic_solution_3cpt_metab <- function(kR, pars, init=c(0,0,0,0)) {

  # browser()
  list2env(as.list(pars), envir = environment())

  if(!("k20" %in% ls())){
    k20 <- 0
  }
  if(!("k30" %in% ls())){
    k30 <- 0
  }

  kme <- ke0 # k41
  # km  <- kme / 1e5 # k14 Absorption into the effect site is much slower than elimination --> as soon as any drug enters, it is eliminated
  # v4  <- v1 / 1e5
  km = 0
  v4 = 0
  k20 <- 0
  k30 <- 0
  E1 <- k10+k12+k13+km
  E2 <- k21+k20
  E3 <- k31+k30

  a <- E1+E2+E3
  b <- E1*E2+E3*(E1+E2)-k12*k21-k13*k31
  c <- E1*E2*E3-E3*k12*k21-E2*k13*k31

  m <- (3*b - a^2)/3
  n <- (2*a^3 - 9*a*b + 27*c)/27
  Q <- (n^2)/4 + (m^3)/27

  alpha <- sqrt(-1*Q)
  beta <- -1*n/2
  gamma <- sqrt(beta^2+alpha^2)
  theta <- atan2(alpha,beta)

  lambda1 <- a/3 + gamma^(1/3)*(cos(theta/3) + sqrt(3)*sin(theta/3))
  lambda2 <- a/3 + gamma^(1/3)*(cos(theta/3) - sqrt(3)*sin(theta/3))
  lambda3 <- a/3 -(2*gamma^(1/3)*cos(theta/3))

  A1last <- init[1]*v1
  A2last <- init[2]*v2
  A3last <- init[3]*v3
  Amlast <- init[4]*v4
  Doserate <- kR

  B = A2last*k21+A3last*k31
  C = E3*A2last*k21+E2*A3last*k31
  I = A1last*k12*E3-A2last*k13*k31+A3last*k12*k31
  J = A1last*k13*E2+A2last*k13*k21-A3last*k12*k21

  # return concentration functions
  a_1 <- function(t) {
    A1term1 = A1last*(exp(-t*lambda1)*(E2-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E2-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E2-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))
    A1term2 = exp(-t*lambda1)*(C-B*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(B*lambda2-C)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(B*lambda3-C)/((lambda1-lambda3)*(lambda3-lambda2))
    A1term3 = Doserate*((E2*E3)/(lambda1*lambda2*lambda3)-exp(-t*lambda1)*(E2-lambda1)*(E3-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-t*lambda2)*(E2-lambda2)*(E3-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-t*lambda3)*(E2-lambda3)*(E3-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)))
    return(A1term1+A1term2+A1term3)
  }

  a_2 <- function(t) {
    A2term1 = A2last*(exp(-t*lambda1)*(E1-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E1-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E1-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))
    A2term2 = exp(-t*lambda1)*(I-A1last*k12*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(A1last*k12*lambda2-I)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(A1last*k12*lambda3-I)/((lambda1-lambda3)*(lambda3-lambda2))
    A2term3 = Doserate*k12*(E3/(lambda1*lambda2*lambda3)-exp(-t*lambda1)*(E3-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-t*lambda2)*(E3-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-t*lambda3)*(E3-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)))
    return(A2term1+A2term2+A2term3)
  }

  a_3 <- function(t) {
    A3term1 = A3last*(exp(-t*lambda1)*(E1-lambda1)*(E2-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E1-lambda2)*(E2-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E1-lambda3)*(E2-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))
    A3term2 = exp(-t*lambda1)*(J-A1last*k13*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(A1last*k13*lambda2-J)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(A1last*k13*lambda3-J)/((lambda1-lambda3)*(lambda3-lambda2))
    A3term3 = Doserate*k13*(E2/(lambda1*lambda2*lambda3)-exp(-t*lambda1)*(E2-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-t*lambda2)*(E2-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-t*lambda3)*(E2-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)))
    return(A3term1+A3term2+A3term3)
  }

  a_4 <- function(t) {
    Amterm1 = Amlast*exp(-t*kme) +km*A1last*(exp(-t*lambda1)*(E2-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1)*(kme-lambda1))+exp(-t*lambda2)*(E2-lambda2)*(E3-lambda2)/((kme-lambda2)*(lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E2-lambda3)*(E3-lambda3)/((kme-lambda3)*(lambda1-lambda3)*(lambda2-lambda3))+exp(-t*kme)*(E2-kme)*(E3-kme)/((lambda1-kme)*(lambda2-kme)*(lambda3-kme)))
    Amterm2 = km*(exp(-t*lambda1)*(B*lambda1-C)/((lambda1-lambda2)*(lambda1-lambda3)*(lambda1-kme))+exp(-t*lambda2)*(C-B*lambda2)/((lambda1-lambda2)*(lambda2-lambda3)*(lambda2-kme))+exp(-t*lambda3)*(C-B*lambda3)/((lambda1-lambda3)*(lambda3-lambda2)*(lambda3-kme))-exp(-t*kme)*(B*kme-C)/((lambda1-kme)*(kme-lambda2)*(kme-lambda3)))
    Amterm3 = km*Doserate*((E2*E3)/(lambda1*lambda2*lambda3*kme)-exp(-t*lambda1)*(E2-lambda1)*(E3-lambda1)/(lambda1*(kme-lambda1)*(lambda2-lambda1)*(lambda3-lambda1))-exp(-t*lambda2)*(E2-lambda2)*(E3-lambda2)/(lambda2*(kme-lambda2)*(lambda1-lambda2)*(lambda3-lambda2))-exp(-t*lambda3)*(E2-lambda3)*(E3-lambda3)/(lambda3*(kme-lambda3)*(lambda1-lambda3)*(lambda2-lambda3))-exp(-t*kme)*(E2-kme)*(E3-kme)/(kme*(lambda1-kme)*(lambda2-kme)*(lambda3-kme)))
    return(Amterm1+Amterm2+Amterm3)
  }

  c_1 <- function(t) a_1(t)/v1
  c_2 <- function(t) a_2(t)/v2
  c_3 <- function(t) a_3(t)/v3
  c_4 <- function(t) a_4(t)/v4

  out <- list(c_1=c_1, c_2=c_2, c_3=c_3, c_4=c_4)
  attributes(out) <- list(infusion = kR, param = pars, initial = init)

  return(out)
}
attributes(pk_basic_solution_3cpt_metab) <- list(ncompartments = 4)



iterate_sol <- function(ivt, basic_sol, ...)
{
  ncpt <- attr(basic_sol,"ncompartments")
  basic_sol_args <- list(...)

  # browser()
  # if initial values aren't specified, set to 0 for all compartments
  if("init" %in% names(basic_sol_args)) init <- basic_sol_args$init
  else init <- rep(0,ncpt)

  if(!any(sapply(ivt, is.list))) ivt <- list(ivt)

  ibe <- sapply(ivt, `[[`, 'begin')
  ied <- sapply(ivt, `[[`, 'end')
  prd <- sort(unique(c(0, c(ibe,ied), Inf)))
  rits <- list()
  for(i in 1:(length(prd)-1)) {
    civt <- sapply(ivt, function(iv) {
      if(prd[i] >= iv$begin && prd[i] < iv$end) {
        iv$k_R
      } else { 0 }
    })

    rit <- list(begin=prd[i], end=prd[i+1],
                idose=sum(civt))
    if(i == 1) {
      rit$init <- init
    } else {
      rit$init <- unname(sapply(which(sapply(rits[[i-1]], is.function)), function(j)
        rits[[i-1]][[j]](rits[[i-1]]$end-rits[[i-1]]$begin)))
    }

    basic_sol_args$init <- rit$init
    sol <- do.call(basic_sol, c(kR = rit$idose, basic_sol_args))
    rits[[i]] <- c(rit, sol)
  }

  ret <- function(tms) {
    # browser()
    sapply(tms, function(t) {
      val <- rep(NA,ncpt)
      for(rit in rits) {
        if(t >= rit$begin && t <= rit$end) {
          val <- sapply(which(sapply(rit, is.function)), function(j) c(rit[[j]](t-rit$begin)))
          break
        }
      }
      val
    })
  }
  return(ret)
}
#' @examples
#' pars_1cpt <- c(k10=1.5,v1=10)
#' pars_2cpt <- c(k10=1.5,k12=0.15,k21=0.09,v1=10,v2=15)
#' pars_3cpt <- c(k10=1.5,k12=0.15,k21=0.09,k13=0.8,k31=0.8,v1=10,v2=15,v3=100,ke0=4)
#' con1cpt <- pk_basic_solution_1cpt(kR = 1, pars = pars_1cpt)
#' con2cpt <- pk_basic_solution_2cpt(kR = 1, pars = pars_2cpt)
#' con3cpt <- pk_basic_solution_3cpt_metab(kR = 1, pars = pars_3cpt)
#' ivt_test = list(list(begin = 0, end = 1/6, k_R = 1), list(begin = 3, end = 3+1/6, k_R = 1))
#' sol_1cpt <- iterate_sol(ivt = ivt_test, basic_sol = pk_basic_solution_1cpt, pars = pars_1cpt, init = 0)
#' sol_2cpt <- iterate_sol(ivt = ivt_test, basic_sol = pk_basic_solution_2cpt, pars = pars_2cpt, init = c(0,0))
#' sol_3cpt <- iterate_sol(ivt = ivt_test, basic_sol = pk_basic_solution_3cpt_metab, pars = pars_3cpt, init = c(0,0,0,0))
#' tms <- seq(0,10,0.01)
#' plot(tms, sol_1cpt(tms), type = "l")
#' lines(tms, sol_2cpt(tms)[1,], col = 2)
#' lines(tms, sol_3cpt(tms)[1,], col = 3)


# create pkmod s3 class object
create_pkmod <- function(ncpts = c("cmpt1","cmpt2","cmpt3m"), pars, init = NULL, ivt = NULL, vcov = NULL){

  if(!is.numeric(pars)) stop("PK parameters must be passed as a vector of numeric values")

  # identify basic solution corresponding to specified compartment model
  basic_sol_opts <- c("pk_basic_solution_1cpt","pk_basic_solution_2cpt","pk_basic_solution_3cpt_metab")
  basic_sol_name <- basic_sol_opts[c("cmpt1","cmpt2","cmpt3m") == ncpts]
  basic_sol <- get(basic_sol_name)
  ncpt <- c(1,2,3)[c("cmpt1","cmpt2","cmpt3m") == ncpts]

  pkmod <- list(basic_sol = basic_sol, pars = pars, ncpt = ncpt, init = init, ivt = ivt, vcov = vcov)

  # if infusion schedule is provided, add piecewise solution to list
  if(!is.null(ivt)){

    if(is.null(init)){
      print("Calculated solution assumes no prior medication in any compartment")
      init <- formals(basic_sol)$init
    }
    sol <- iterate_sol(ivt, basic_sol, pars = pars, init = init)
    pkmod <- c(pkmod, sol = list(sol))
  }

  class(pkmod) <- "pkmod"
  return(pkmod)
}
#' @examples
#' pk1cpt <- create_pkmod(ncpts = "cmpt1", pars = list(k10=-2.496,v1=2.348), vcov = diag(2), ivt = list(list(begin = 0, end = 1/6, k_R = 1), list(begin = 3, end = 3+1/6, k_R = 1)))
#' class(pk1cpt)


#' create update method that changes elements of the pkmod object and (by default) changes downstream effects
#' e.g. if parameters are updated and an infusion schedule exists,
#' @param pkmod A list with class pkmod
#' @param args_update A named list of arguments to the pkmod object that are to be changed
update.pkmod <- function(pkmod, args_update){

  if(!all(names(args_update) %in% attributes(pkmod)$names)){
    stop(paste("Argument names must be limited to the following:", paste(attributes(pk1cpt)$names, collapse = ", ")))
  }

  # replace values with updated versions
  pkmod[names(args_update)] <- args_update

  # update pk solution if an infusion schedule has been specified - this is downstream of any updates
  if(!is.null(pkmod$ivt)){

    if(is.null(pkmod$init)){
      print("Calculated solution assumes no prior medication in any compartment")
      init <- formals(pkmod$basic_sol)$init
    }

    pkmod$sol <- iterate_sol(pkmod$ivt, pkmod$basic_sol, pars = pkmod$pars, init = pkmod$init)
  }

  return(pkmod)
}
#' @examples
#' pk1cpt <- create_pkmod(ncpts = "cmpt1", pars = c(k10=0.5,v1=5), vcov = diag(2), ivt = list(list(begin = 0, end = 1/6, k_R = 1), list(begin = 3, end = 3+1/6, k_R = 1)))
#' pk1cpt$pars
#' plot(seq(0,8,0.1), pk1cpt$sol(seq(0,8,0.1)), type = "l", xlab = "min", ylab = "concentration", ylim = c(0,0.07))
#' pars2 <- c(k10=0.3,v1=4)
#' pk1cpt <- update(pk1cpt, list(pars = pars2))
#' pk1cpt$pars
#' lines(seq(0,8,0.1), pk1cpt$sol(seq(0,8,0.1)), col = 2)



# create summary method
summary.pkmod <- function(pkmod, alpha = 0.05, ...){
  browser()
  lb <- pkmod$pars + qnorm(alpha/2)*diag(pkmod$vcov)
  ub <- pkmod$pars - qnorm(alpha/2)*diag(pkmod$vcov)
  out <- data.frame(parameters = names(pkmod$pars), value = unlist(pkmod$pars))
  if(!is.null(pkmod$vcov)){
    out <- cbind(out, lb = lb, ub = ub)
  }
  knitr::kable(out, row.names = F, ...)
}
#' @examples
#' pk1cpt <- create_pkmod(ncpts = "cmpt1", pars = k10=-2.496,v1=2.348, vcov = diag(2), ivt = list(list(begin = 0, end = 1/6, k_R = 1), list(begin = 3, end = 3+1/6, k_R = 1)))
#' summary(pk1cpt, digits = 3)


# create predict method
predict.pkmod <- function(pkmod, tms = NULL, cmpt = 1, plot_con = T){

  if(is.null(pkmod$ivt)){
    stop("An infusion schedule has not yet been associated with the pkmod object.
         Use the 'update' function to add an infusion schedule.")
  }

  if(is.null(tms)){
    inf_start <- pkmod$ivt[[1]]$begin
    inf_end <- pkmod$ivt[[length(pkmod$ivt)]]$end
    tms <- seq(inf_start, inf_end + (inf_end - inf_start), length.out = 100)
  }

  con_p <- pkmod$sol(tms)

  # browser()

  f <- function(pars) unname(iterate_sol(pkmod$ivt, pkmod$basic_sol, pars = pars, init = pkmod$init)(tms))

  gd <- fdGrad(pars = pkmod$pars, fun = f)


}
#' @examples
#' pk1cpt <- create_pkmod(ncpts = "cmpt1", pars = list(k10=0.5,v1=5), vcov = diag(2), ivt = list(list(begin = 0, end = 1/6, k_R = 1), list(begin = 3, end = 3+1/6, k_R = 1)))
#' predict(pk1cpt)


# potentially keep population prior as a separate object with specifications of CDF so that





# TCI functions


# TCI function (possibly user-specified) - take in a target concentration and a pk model object and return a single infusion rate
# Note: requires that user define function in terms of basic PK functions available

#' @param Ct Target plasma concentration as scalar or vector. If vector, infusions are independent of eachother.
#' @param pkfn PK function
#' @param cmpt Compartment number for plasma concentration
#' @param dt Duration of infusion
#' @param max_kR Maximum infusion rate administered by TCI device (g/h)
#' @param ... Parameters to be passed on to pkfn. PK parameters are required, intial values can be passed on.
TCI_plasma <- function(Ct, pkfn, cmpt = 1, dt = 1/6, max_kR = 10000, ...){
  # if(is.null(cmpt)) warning("Using first compartment as plasma concentration")
  cmpt <- 1; x1 = 1; x2 = 2
  sol_x1 <- pkfn(kR = x1, ...)
  sol_x2 <- pkfn(kR = x2, ...)
  Cp_x1 <- sol_x1[[cmpt]](dt)
  Cp_x2 <- sol_x2[[cmpt]](dt)
  m <- (Cp_x2 - Cp_x1) / (x2-x1)
  b <- Cp_x1 - m*x1
  kR <- (Ct - b) / m
  kR <- ifelse(kR < 0, 0, kR)
  kR <- ifelse(kR > max_kR, max_kR, kR)
  attr(kR, "dt") <- dt
  return(kR)
}

#' @examples
#' lpars_3cpt <- list(k10=-2.496,k12=-1.992,k21=-2.376,k13=-2.868,k31=-5.244,v1=2.348,v2=2.731,v3=4.724,ke0=1.302)
#' TCI_plasma(3, pk_basic_solution_3cpt_metab, pars = lapply(lpars_3cpt,exp), init = c(1,0,0,0))


#' @param Ct Target effect-site concentration
#' @param pkfn PK function object
#' @param cmpt Compartment number for effect-site concentration (assumed to be last if unspecified)
#' @param dt Duration of infusion
#' @param max_kR Maximum infusion rate administered by TCI device (g/h)
#' @param ... Parameters to be passed on to pkfn. PK parameters are required, intial values can be passed on.
#' @param tmax_search Maximum time (min) until peak is observed following an infusion of duration dt. Defaults to 20 minutes.
#' @param grid_space Time resolution (min) associated with search for peak concentration. Defaults to 1 second.
#' @param return_all Logical. Return time of maximum concentration and vector of concentrations in addition to the infusion rate.
TCI_effect <- function(Ct, pkfn, cmpt = NULL, dt = 1/6, max_kR = 10000, tmax_search = 20, grid_space = 1/60, return_all = F, ...){
  if(is.null(cmpt)) cmpt <- attr(pkfn,"ncompartments")
  pkfn_args_init0 <- list(...)
  if(is.null(pkfn_args_init0$init)) init <- eval(formals(pkfn)$init)
  else init <- pkfn_args_init0$init
  pkfn_args_init0$init <- rep(0, attr(pkfn,"ncompartments"))

  extract_col_vec <- function(x, col = 1){
    if(is.vector(x)) return(x)
    else return(x[,col])
  }

  E <- do.call(iterate_sol, c(ivt = list(list(begin = 0, end = dt, k_R = 1)), basic_sol = pkfn, pkfn_args_init0))
  B <- iterate_sol(ivt = list(list(begin = 0, end = dt, k_R = 0)), basic_sol = pkfn, ...)

  tm_seq <- seq(grid_space,tmax_search,grid_space)
  ceproj <- E(tm_seq)[cmpt,]
  tpeak = tm_seq[which.max(ceproj)]

  if(all(init == 0)){
    kR <- unname((Ct-B(tpeak)[cmpt,]) / E(tpeak)[cmpt,])
    if(return_all) return(c(kR = kR, jpeak = tpeak, Ct = E(tpeak)[cmpt,]))
    else return(kR)
  } else{
    tm_seq2 <- seq(grid_space,tpeak+0.5,grid_space)
    jpeak0 = tpeak - 0.1
    jpeak1 = jpeak0 + 0.1
    while(jpeak0 != jpeak1){
      jpeak0 = jpeak1
      I0 = (Ct - B(jpeak0)[cmpt,]) / E(jpeak0)[cmpt,]
      # ceproj = B(tm_seq2)[cmpt,] + E(tm_seq2)[cmpt,]*I0
      ceproj = B(tm_seq2)[cmpt,] + E(tm_seq2)[cmpt,] %*% t(I0)
      # the peak is always at the same point regardless of infusion rate
      jpeak1 = tm_seq2[which.max(extract_col_vec(ceproj))]
    }
    kR = unname((Ct-B(jpeak1)[cmpt]) / E(jpeak1)[cmpt])
    kR <- ifelse(kR < 0, 0, kR)
    kR <- ifelse(kR > max_kR, max_kR, kR)

    if(return_all) return(c(kR = kR, jpeak = jpeak1, Ct = B(jpeak1)[cmpt] + E(jpeak1)[cmpt]*I0))
    else return(kR)
  }
}

#' @examples
#' lpars_3cpt <- list(k10=-2.496,k12=-1.992,k21=-2.376,k13=-2.868,k31=-5.244,v1=2.348,v2=2.731,v3=4.724,ke0=1.302)
#' TCI_effect(3, pk_basic_solution_3cpt_metab, pars = lapply(lpars_3cpt,exp)) # uses default initial values for pkfn.
#' TCI_effect(c(3,4,5), pk_basic_solution_3cpt_metab, pars = lapply(lpars_3cpt,exp), init = c(1,0,0,0))
#'
#' First compartment is equivalent to plasma-targeting
#' TCI_plasma(c(3,4,5), pk_basic_solution_3cpt_metab, pars = lapply(lpars_3cpt,exp), init = c(1,0,0,0))
#' TCI_effect(c(3,4,5), pk_basic_solution_3cpt_metab, pars = lapply(lpars_3cpt,exp), init = c(1,0,0,0), cmpt = 1)



#' If the infusion rate is provided

#' Function to iterate TCI algorithm
#' We want to minimize the number of PK function calls
#' Process:
#' 1. Initialize values - matrix for concentrations
#' 2. Calculate infusion necessary to reach target
#' 3.

#' Function takes in arguments for concentrations (and possibly times), an algorithm, and a pk function
#' Returns a vector of infusion rates
#' (?) Also return bolus/maintenance infusion rates and times associated with reaching the target.

TCI <- function(Ct, tci_algorithm, pkfn, start = 0, end = NULL, ...){

  pkfn_args_init0 <- list(...)
  cmpt <- formals(tci_algorithm)$cmpt
  if(is.null(cmpt)) cmpt <- attr(pkfn,"ncompartments")
  if(is.null(pkfn_args_init0$init)) init <- eval(formals(pkfn)$init)
  else init <- pkfn_args_init0$init

  dt <- eval(formals(tci_algorithm)$dt)

  if(is.vector(Ct)){
    end <- length(Ct)*dt
    starttms <- seq(start,end,by=dt)
    sf <- stepfun(starttms, c(0))
    endtms <- starttms+dt
  } else{
    sf <- stepfun(Ct$time, c(init[cmpt],Ct$Ct))
  }
  endtm <- ifelse(is.matrix(Ct), )
  tm_seq <- seq(starttm, )


}

#' @examples
#' TCI(c(3,4,5), tci_algorithm = TCI_effect, pkfn = pk_basic_solution_3cpt_metab, pars = lapply(lpars_3cpt,exp), init = c(1,0,0,0))
#' tms <- seq(0,10,1/6)
#' Ct <- c(rep(3,1/(1/6)), rep(4,1/(1/6)), rep(5,4/(1/6)), rep(4,4/(1/6)), 4)
#' df <- data.frame(time = tms, Ct = Ct)
#' TCI(df, tci_algorithm = TCI_effect, pkfn = pk_basic_solution_3cpt_metab, pars = lapply(lpars_3cpt,exp), init = c(1,0,0,0))


TCI_gen <- function(pkfn, tmax_search = 20, grid_space = 1/60, cmpt = NULL, return_jpeak = F, delta = 1/6, max_kR = 10000, ...){

  if(is.null(cmpt)) cmpt <- attr(pkfn,"ncompartments")

  if(cmpt == 1){ # plasma targeting
    x1 = 0.1; x2 = 0.2

    TCI_fn <- function(Ce, init){
      sol_x1 <- pkfn(kR = x1, init = init, ...)
      sol_x2 <- pkfn(kR = x2, init = init, ...)
      Cp_x1 <- sol_x1[[1]](delta)
      Cp_x2 <- sol_x2[[1]](delta)
      m <- (Cp_x2 - Cp_x1) / (x2-x1)
      b <- Cp_x1 - m*x1
      infrt_try <- (Ce - b) / m
      infrt <- ifelse(infrt_try < 0, 0, infrt_try)
      if(infrt > max_kR) infrt <- max_kR
      return(infrt)
    }

  } else{
    # course with infusion starting from 0 concentration
    E <- iterate_sol(ivt = list(list(begin = 0, end = 1/6, k_R = 1)), basic_sol = pkfn, init = rep(0, attr(pkfn,"ncompartments")), ...)

    tm_seq <- seq(grid_space,tmax_search,grid_space)
    ceproj <- E(tm_seq)[cmpt,]
    tpeak = tm_seq[which.max(ceproj)]

    TCI_fn <- function(Ce, init){
      # course without infusion - use current concentration
      B <- iterate_sol(ivt = list(list(begin = 0, end = 1/6, k_R = 0)), basic_sol = pkfn, init = init, ...)

      if(all(init == 0)){
        kR <- unname((Ce-B(tpeak)[cmpt]) / E(tpeak)[cmpt])
        if(return_jpeak) return(c(kR = kR, jpeak = tpeak, Ce = E(tpeak)[cmpt,]))
        else return(kR)
      } else{
        tm_seq2 <- seq(grid_space,tpeak+0.5,grid_space)
        jpeak0 = tpeak - 0.1
        jpeak1 = jpeak0 + 0.1
        while(jpeak0 != jpeak1){
          jpeak0 = jpeak1
          I0 = (Ce - B(jpeak0)[cmpt]) / E(jpeak0)[cmpt]
          ceproj = B(tm_seq2)[cmpt,] + E(tm_seq2)[cmpt,]*I0
          jpeak1 = tm_seq2[which.max(ceproj)]
        }
        kR = unname((Ce-B(jpeak1)[cmpt]) / E(jpeak1)[cmpt])
        if(kR < 0) kR = 0
        if(return_jpeak) return(c(kR = kR, jpeak = jpeak1, Ce = B(jpeak1)[cmpt] + E(jpeak1)[cmpt]*I0))
        else return(kR)
      }
    }
  }
  mostattributes(TCI_fn) <- attributes(pkfn)
  return(TCI_fn)
}

#' @examples
#' lpars_1cpt <- list(k10=-2.496,v1=2.348)
#' lpars_3cpt <- list(k10=-2.496,k12=-1.992,k21=-2.376,k13=-2.868,k31=-5.244,v1=2.348,v2=2.731,v3=4.724,ke0=1.302)
#' TCI_3cpt_effect <- TCI_gen(pkfn = pk_basic_solution_3cpt_metab, pars = lapply(lpars_3cpt,exp), cmpt = 4)
#' TCI_1cpt_central <- TCI_gen(pkfn = pk_basic_solution_1cpt, pars = lapply(lpars_1cpt,exp), cmpt = 1)
#' TCI_3cpt_effect(1,c(0,0,0,0))
#'
#' TCI_1cpt_central(Ce = 5, init = 0)


# take in component TCI function and initial value - return function that will calculate infusions to reach a set of
# targets at (possibly) set times
iterate_tci <- function(tci, init, start_tm = 0, ...){

  tci_out <- function(Ce, tms = NULL, ...){
    if(!is.null(tms)){

    }

    infrt <- rep(NA, length(Cet))
    begin <- seq(start_tm,(length(Cet)-1)*delta, delta)
    end <- seq(delta,length(Cet)*delta, delta)
    con <- matrix(NA, nrow = length(Cet)+1, ncol = attr(tci,"ncompartments"))
    con[1,] <- init


  }

}

#' @examples
#' Ce <- c(1,1.5,2,1.5)
#' tms <- 1:4


