# New test comment
pk_basic_solution_1cpt <- function(kR, pars, init = 0){

  # parameters should be named list with values "k10","v" OR "CL","v"
  list2env(pars, envir = environment())

  if("CL" %in% ls()){
    k10 <- CL / v
  }

  init_amt <- init * v1
  a_1 <- function(t) kR/k10*(1-exp(-t*k10))+init_amt*exp(-t*k10)
  c_1 <- function(t) a_1(t)/v1

  out <- list(c_1=c_1)
  attributes(out) <- list(infusion = kR, param = pars, initial = init)
  return(out)
}
attributes(pk_basic_solution_1cpt) <- list(ncompartments = 1)



# two compartment model
pk_basic_solution_2cpt <- function(kR, pars, init=c(0,0,0,0)){

  list2env(pars, envir = environment())

  if("CL" %in% ls()){
    k10 <- CL / v1
    k12 <- Q / v1
    k21 <- Q / v2
  }
    k20 <- 0
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


pk_basic_solution_3cpt_metab <- function(kR, pars, init=c(0,0,0,0)) {

  list2env(pars, envir = environment())

  kme <- ke0 # k41
  km  <- kme / 1e5 # k14 Absorption into the effect site is much slower than elimination --> as soon as any drug enters, it is eliminated
  v4  <- v1 / 1e5
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

#' examples
# lpars_1cpt <- list(k10=-2.496,v1=2.348)
# lpars_3cpt <- list(k10=-2.496,k12=-1.992,k21=-2.376,k13=-2.868,k31=-5.244,v1=2.348,v2=2.731,v3=4.724,ke0=1.302)
# con1cpt <- pk_basic_solution_1cpt(kR = 1, pars = lapply(lpars_1cpt,exp))
# con3cpt <- pk_basic_solution_3cpt_metab(kR = 1, pars = lapply(lpars_3cpt,exp))
# ivt_test = list(list(begin = 0, end = 1/6, k_R = 1), list(begin = 3, end = 3+1/6, k_R = 1))
# sol_1cpt <- iterate_sol(ivt = ivt_test, basic_sol = pk_basic_solution_1cpt, pars = list(k10 = exp(0.1), v1 = 10), init = 0)
# sol_3cpt <- iterate_sol(ivt = ivt_test, basic_sol = pk_basic_solution_3cpt_metab, pars = lapply(lpars_3cpt,exp), init = c(0,0,0,0))
# tms <- seq(0,10,0.01)
# plot(tms, sol_3cpt(tms)[1,], type = "l", ylim = c(0,0.03))
# lines(tms, sol_1cpt(tms), col = 2)



#' Function to update

