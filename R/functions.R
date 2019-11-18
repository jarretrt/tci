# --------------------------------------------------------------------------------------------------------------------------------
# - Functions --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------

#' 3 compartment IV infusion with first-order absorption between compartments and with an additional effect-site compartment.
#' The analytical solutions implemented in this function are provided in "ADVAN-style analytical solutions for common pharmacokinetic models" by
#' Abuhelwa et al. 2015.
#'
#' This function takes in arguments for each of the absorption and elimination rate constants of a three-compartment model
#' as well as initial concentrations, c0. ke0 gives the rate of elimination from the effect-site compartment into the
#' central compartment (i.e. k41). The rate of absorption into the effect-site compartment is set at 1/10,000 the value of ke0.
#' The function returns a set of functions that calculate the concentration in each of the four compartments as a function of
#' time.
#' @param kR Infusion rate (e.g. ml/min).
#' @param k10 Rate of excretion from central compartment.
#' @param k12 Rate of transfer from compartment 1 to compartment 2.
#' @param k21 Rate of transfer from compartment 2 to compartment 1.
#' @param k13 Rate of transfer from compartment 1 to compartment 3.
#' @param k31 Rate of transfer from compartment 3 to compartment 1.
#' @param ke0 Rate of transfer from effect-site compartment to compartment 1.
#' @param v1 Volume of compartment 1.
#' @param v2 Volume of compartment 2.
#' @param v3 Volume of compartment 3.
pk_basic_solution_3cpt_metab <- function(kR,
                                         k10 = k_10_d,
                                         k12 = k_12_d,
                                         k21 = k_21_d,
                                         k13 = k_13_d,
                                         k31 = k_31_d,
                                         v1 = v_1_d,
                                         v2 = v_2_d,
                                         v3 = v_3_d,
                                         ke0 = k_e0_d,
                                         c0=c(0,0,0,0)) {


  kme <- ke0 # k41
  km  <- ke0 / 1e5 # k14 Absorption into the effect site is much slower than elimination --> as soon as any drug enters, it is eliminated
  v4  = v1 / 1e5
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

  A1last <- c0[1]*v1
  A2last <- c0[2]*v2
  A3last <- c0[3]*v3
  Amlast <- c0[4]*v4
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

  return(list(c_1=c_1, c_2=c_2, c_3=c_3, c_4=c_4))
}

# compiled version
pk_basic_solution_3cpt_metab_c <- compiler::cmpfun(pk_basic_solution_3cpt_metab)

#' This function extends the function pk_basic_solution_3cpt_metab to a specified infusion schedule, rather than a single
#' infusion.
#' @param pars Named vector of parameters for a 3-compartment model with effect-site.
#' @param ivt Infusion schedule given in the form of a named list
#' (e.g. list(list(begin = 0, end = 2, k_R = 1), list(begin = 4, end = 6, k_R = 1)))
#' @param init inital concentrations for the 4 compartments.
pk_solution_3cpt_metab <- function(pars, ivt = ivt_d, init=init_d)
{
  k_10 = pars[1]
  k_12 = pars[2]
  k_21 = pars[3]
  k_13 = pars[4]
  k_31 = pars[5]
  v_1  = pars[6]
  v_2  = pars[7]
  v_3  = pars[8]
  k_e0 = pars[9]

  # k is associted with the rate of change between compartments
  # d indicates the test value
  # ivt = intravenous transfusion (schedule)
  ## create a list of event times
  ibe <- sapply(ivt, `[[`, 'begin') # extract the "begin" item (transfusion start time) from each list
  ied <- sapply(ivt, `[[`, 'end') # extract the transfusion end time from each list
  # browser()
  prd <- sort(unique(c(0, c(ibe,ied), Inf))) # order all time periods
  rits <- list()
  ## compute basic solution in each interval
  for(i in 1:(length(prd)-1)) { # there are length(prd)-1 = 10 intervals, split into half-hours
    civt <- sapply(ivt, function(iv) { # change in concentration due to IV transfusion in first half-hour of each of 5 sessions
      if(prd[i] >= iv$begin && prd[i] < iv$end) { # prd[i] = iv$end in the second half-hour
        iv$k_R
      } else { 0 }
    })

    rit <- list(begin=prd[i], end=prd[i+1],
                idose=sum(civt))
    if(i == 1) {
      rit$init <- init
    } else {
      rit$init <- c(rits[[i-1]]$c_1(rits[[i-1]]$end-rits[[i-1]]$begin),
                    rits[[i-1]]$c_2(rits[[i-1]]$end-rits[[i-1]]$begin),
                    rits[[i-1]]$c_3(rits[[i-1]]$end-rits[[i-1]]$begin),
                    rits[[i-1]]$c_4(rits[[i-1]]$end-rits[[i-1]]$begin))
    }

    sol <- pk_basic_solution_3cpt_metab_c(kR=rit$idose, k10 = k_10, k12 = k_12, k21 = k_21, k13 = k_13,
                                          k31 = k_31, v1 = v_1, v2 = v_2, v3 = v_3, ke0 = k_e0, c0=rit$init)

    rits[[i]] <- c(rit, sol)
  }

  ret <- function(tms) {
    sapply(tms, function(t) {
      val <- rep(NA,4)
      for(rit in rits) {
        if(t >= rit$begin && t <= rit$end) {
          val <- c(rit$c_1(t-rit$begin),
                   rit$c_2(t-rit$begin),
                   rit$c_3(t-rit$begin),
                   rit$c_4(t-rit$begin))
          break
        }
      }
      val
    })
  }
  return(ret)
}

#' Piecewise solution for a single infusion followed by a period with no infusion.
#' This function is similar to pk_solution_3cpt_metab, except that it accepts and
#' implements only the first infusion (ivt[[1]]). This function exists primarily for
#' reducing computational speed when searching for time until maximum concentration.
#'
#' @param pars Named vector of parameters for a 3-compartment model with effect-site.
#' @param ivt Infusion schedule given in the form of a named list
#' (e.g. list(list(begin = 0, end = 2, k_R = 1), list(begin = 4, end = 6, k_R = 1)))
#' @param init Inital concentrations for the 4 compartments.
pk_solution_3cpt_metab_singleinf <- function(pars, ivt, init, ce_only = F){

  k_10 = pars[1]
  k_12 = pars[2]
  k_21 = pars[3]
  k_13 = pars[4]
  k_31 = pars[5]
  v_1  = pars[6]
  v_2  = pars[7]
  v_3  = pars[8]
  k_e0 = pars[9]

  begin <- ivt$begin
  end <- ivt$end
  delta <- end - begin

  sol_inf <- pk_basic_solution_3cpt_metab_c(kR=ivt$k_R, k10 = k_10, k12 = k_12, k21 = k_21, k13 = k_13,
                                            k31 = k_31, v1 = v_1, v2 = v_2, v3 = v_3, ke0 = k_e0, c0=init)
  c_end_inf <- c(sol_inf$c_1(delta), sol_inf$c_2(delta), sol_inf$c_3(delta), sol_inf$c_4(delta))
  sol_coast <- pk_basic_solution_3cpt_metab_c(kR=0, k10 = k_10, k12 = k_12, k21 = k_21, k13 = k_13,
                                              k31 = k_31, v1 = v_1, v2 = v_2, v3 = v_3, ke0 = k_e0, c0=c_end_inf)

  if(ce_only){
    calc_con <- function(tms){
      tms_inf <- tms[tms <= end]
      tms_coast <- tms[tms > end]
      return(c(sol_inf$c_4(tms_inf-begin), sol_coast$c_4(tms_coast-end)))
    }
  } else{
    calc_con <- function(tms){
      tms_inf <- tms[tms <= end]
      tms_coast <- tms[tms > end]
      val_inf <- cbind(sol_inf$c_1(tms_inf-begin), sol_inf$c_2(tms_inf-begin), sol_inf$c_3(tms_inf-begin), sol_inf$c_4(tms_inf-begin))
      val_coast <- cbind(sol_coast$c_1(tms_coast-end), sol_coast$c_2(tms_coast-end), sol_coast$c_3(tms_coast-end), sol_coast$c_4(tms_coast-end))
      return(t(rbind(val_inf, val_coast)))
    }
  }
  return(calc_con)
}


#' Function for calculating a TCI infusion schedule corresponding to a set of target concentrations. Primary inputs are a vector
#' of target effect-site concentrations and a set of parameters for a 3-compartment PK model. Initial concentrations (init)
#' can be set, as can the frequency with which the TCI device updates the target (delta), the maximum infusion rate, and
#' the sensitivity associated with switching to plasma-targeting mode (plasma_tol, effect_tol).
#'
#' This function makes use of formulas described by Shafer and Gregg (1992) in "Algorithms to rapidly achieve
#' and maintain stable drug concentrations at the site of drug effect with a computer-controlled infusion pump", as well as
#' the plasma-targeting algorithm described by Jacobs (1990) in "Algorithm for Optimal Linear Model-Based Control with
#' Application to Pharmacokinetic Model-Driven Drug Delivery."
#'
#' @param Cet Numeric vector of target effect-site concentrations. It is assumed that these occur delta seconds/min
#' apart where delta is the frequency with which the algorithm adjusts the target.
#' @param pars Named vector of parameters for a 3-compartment model with effect-site.
#' @param init Inital concentrations for the 4 compartments.
#' @param delta Frequency of TCI updates. Default is 10 seconds, expressed in terms of minutes = 1/6.
#' @param max_kR Maximum infusion rate of TCI pump.
#' @param plasma_tol Maximum percent difference between predicted plasma concentration and target concentration permitted
#' in order to switch to plasma-targeting mode. Plasma-targeting mode is used primarily to increase computational efficiency
#' when the effect-site concentration is sufficiently stable and close to the target concentration.
#' @param effect_tol Maximum percent difference between predicted effect-site concentration and target concentration
#' permitted in order to switch to plasma-targeting mode.


TCI_basic <- function(Ce, pars, init, tmax_search = 20, grid_space = 1/60, return_jpeak = F){
  # course without infusion - use current concentration
  B <- pk_solution_3cpt_metab_singleinf(pars = pars, ivt = list(begin = 0, end = 1/6, k_R = 0), init = init, ce_only = T)
  # course with infusion starting from 0 concentration
  E <- pk_solution_3cpt_metab_singleinf(pars = pars, ivt = list(begin = 0, end = 1/6, k_R = 1), init = c(0,0,0,0), ce_only = T)
  # find tpeak
  tm_seq <- seq(grid_space,tmax_search,grid_space)
  ceproj <- E(tm_seq)
  tpeak = tm_seq[which.max(ceproj)]

  if(all(init == 0)){
    kR <- unname((Ce-B(tpeak)) / E(tpeak))
    if(return_jpeak) return(c(kR = kR, jpeak = tpeak, Ce = E(tpeak)))
    else return(kR)
  } else{
    tm_seq2 <- seq(grid_space,tpeak+0.5,grid_space)
    jpeak0 = tpeak - 0.1# initialize search
    jpeak1 = jpeak0 + 0.1
    while(jpeak0 != jpeak1){
      jpeak0 = jpeak1
      I0 = (Ce - B(jpeak0)) / E(jpeak0)
      ceproj = B(tm_seq2) + E(tm_seq2)*I0
      jpeak1 = tm_seq2[which.max(ceproj)]
    }

    kR = unname((Ce-B(jpeak1)) / E(jpeak1))
    if(kR < 0) kR = 0
    if(return_jpeak) return(c(kR = kR, jpeak = jpeak1, Ce = B(jpeak1) + E(jpeak1)*I0))
    else return(kR)
  }
}


TCI_EffectSite <- function(Cet, pars, delta = 1/6, init = c(0,0,0,0), max_kR = 10000, plasma_tol = 0.1, effect_tol = 0.005, x1 = 0.1, x2 = 0.2, tmax_search = 4, grid_space = 1/120) {
  infrt <- rep(NA, length(Cet))
  begin <- seq(0,(length(Cet)-1)*delta, delta)
  end <- seq(delta,length(Cet)*delta, delta)
  con <- matrix(NA, nrow = length(Cet)+1, ncol = 4)
  con[1,] <- init

  for(i in 1:length(Cet)){
    con[con[i,]<0] <- 0

    # if target is below current concentration, stop infusion
    if(Cet[i] < con[i,4]){
      infrt[i] <- 0
    }

    # plasma targeting
    plasma_cond <- (con[i,1] > Cet[i] - plasma_tol*Cet[i]) &  (con[i,1] < Cet[i] + plasma_tol*Cet[i]) &
      (con[i,4] > Cet[i] - effect_tol*Cet[i]) & (con[i,4] < Cet[i] + effect_tol*Cet[i])

    if(plasma_cond & is.na(infrt[i])){
      sol_x1 <- pk_basic_solution_3cpt_metab_c(k10 = pars[1],
                                               k12 = pars[2],
                                               k21 = pars[3],
                                               k13 = pars[4],
                                               k31 = pars[5],
                                               v1  = pars[6],
                                               v2  = pars[7],
                                               v3  = pars[8],
                                               ke0 = pars[9],
                                               kR  = x1, c0 = con[i,])
      sol_x2 <- pk_basic_solution_3cpt_metab_c(k10 = pars[1],
                                               k12 = pars[2],
                                               k21 = pars[3],
                                               k13 = pars[4],
                                               k31 = pars[5],
                                               v1  = pars[6],
                                               v2  = pars[7],
                                               v3  = pars[8],
                                               ke0 = pars[9],
                                               kR  = x2, c0 = con[i,])
      Cp_x1 <- sol_x1[[1]](delta)
      Cp_x2 <- sol_x2[[1]](delta)
      m <- (Cp_x2 - Cp_x1) / (x2-x1)
      b <- Cp_x1 - m*x1
      infrt_try <- (Cet[i] - b) / m
      infrt[i] <- ifelse(infrt_try < 0, 0, infrt_try)
      if(infrt[i] > max_kR) infrt[i] <- max_kR
    }
    if(is.na(infrt[i])){
      infrt[i] <- TCI_basic(Ce = Cet[i], pars = pars, init = con[i,], tmax_search = tmax_search, grid_space = grid_space)
    }

    sol_x0 <- pk_basic_solution_3cpt_metab_c(k10 = pars[1],
                                             k12 = pars[2],
                                             k21 = pars[3],
                                             k13 = pars[4],
                                             k31 = pars[5],
                                             v1  = pars[6],
                                             v2  = pars[7],
                                             v3  = pars[8],
                                             ke0 = pars[9],
                                             kR  = infrt[i], c0 = con[i,])
    con[i+1,] <- c(sol_x0$c_1(delta), sol_x0$c_2(delta), sol_x0$c_3(delta), sol_x0$c_4(delta))
  }

  if(any(is.nan(infrt))) stop("NaN in infusion rate")
  ivt_out <- lapply(1:length(Cet), function(j) list(begin = begin[j], end = end[j], k_R = infrt[j]))
  return(ivt_out)
}


#' Function to generate PK parameters for Eleveld model.
#' @param theta Vector of fixed effects
#' @param eta Vector of random effect variances
#' @param patient_vars Named list of observed patient characteristics
#' @param returnQ Logical value indicating if clearance parameters should be returned instead of rate parameters
gen_eleveld_pk_pars <- function(theta, eta, patient_vars, returnQ = F){

  AGE  <- patient_vars$AGE
  PMA  <- patient_vars$PMA
  WGT  <- patient_vars$WGT
  HGT  <- patient_vars$HGT
  M1F2 <- patient_vars$M1F2
  TECH <- patient_vars$TECH
  A1V2 <- patient_vars$A1V2

  MALE <- ifelse(M1F2 == 1, 1, 0)
  OPIATE <- ifelse(TECH == 2, 1, 0)
  ARTERIAL <- ifelse(A1V2 == 1, 1, 0)
  AGEref=35
  WGTref=70
  HGTref=170

  if(is.na(PMA)) PMA = AGE + 40/52
  PMAref = AGEref + 40/52
  BMI = 10000 * WGT / HGT / HGT
  BMIref = 10000 * WGTref / HGTref / HGTref

  V1ref = theta[1]
  V2ref = theta[2]
  V3ref = theta[3]

  faging <- function(x) exp(x*(AGE-AGEref))
  fsigmoid <- function(x,E50,lambda) x^lambda / (x^lambda + E50^lambda)
  fcentral <- function(x) fsigmoid(x,theta[12],1)
  fCLmaturation <- fsigmoid(PMA*52,theta[8],theta[9])
  fCLmaturation_ref = fsigmoid(PMAref*52,theta[8],theta[9])
  fQ3maturation <- fsigmoid(AGE*52+40,theta[14],1)
  fQ3maturation_ref <- fsigmoid(AGEref*52+40,theta[14],1)
  fopiates <- function(x) ifelse(OPIATE,exp(x*AGE),1)

  fAlSallami <-     ifelse(MALE, (0.88+(1-0.88)/(1+(AGE/13.4)^(-12.7)))*(9270*WGT)/(6680+216*BMI),
                           (1.11+(1-1.11)/(1+(AGE/7.1)^(-1.1)))*(9270*WGT)/(8780+244*BMI))
  # fAlSallami_ref <- ifelse(MALE, (0.88+(1-0.88)/(1+(AGEref/13.4)^(-12.7))) *(9270*WGTref)/(6680+216*BMIref),
  # (1.11+(1-1.11)/(1+(AGEref/7.1)^(-1.1)))*(9270*WGTref)/(8780+244*BMIref))

  # apparently the reference case doesn't change with sex
  fAlSallami_ref <- (0.88+(1-0.88)/(1+(AGEref/13.4)^(-12.7))) *(9270*WGTref)/(6680+216*BMIref)
  V1arterial <- theta[1]*fcentral(WGT)/fcentral(WGTref)*exp(eta[1])
  V1venous <- V1arterial*(1+theta[17]*(1-fcentral(WGT)))
  V2 <- theta[2]*WGT/WGTref*faging(theta[10])*exp(eta[2])
  V3 <- theta[3]*(fAlSallami/fAlSallami_ref)*fopiates(theta[13])*exp(eta[3])

  CL <- c(t(c(MALE,1-MALE))%*%theta[c(4,15)]*(WGT/WGTref)^0.75*(fCLmaturation/fCLmaturation_ref)*fopiates(theta[11])*exp(eta[4]))
  Q2arterial <- theta[5]*(V2/V2ref)^(0.75)*(1+theta[16]*(1-fQ3maturation))*exp(eta[5])
  Q2venous <- Q2arterial*theta[18]
  Q3 <- theta[6]*(V3/V3ref)^(0.75)*(fQ3maturation/fQ3maturation_ref)*exp(eta[6])

  V1 = ifelse(ARTERIAL, V1arterial, V1venous)
  K10 = CL/V1
  K12 = ifelse(ARTERIAL, Q2arterial/V1, Q2venous/V1)
  K21 = ifelse(ARTERIAL, Q2arterial/V2, Q2venous/V2)
  K13 = Q3/V1
  K31 = Q3/V3

  if(!returnQ){
    return(c(k10 = K10, k12 = K12, k21 = K21, k13 = K13, k31 = K31, v1 = V1, v2 = V2, v3 = V3))
  } else{
    Q2 <- ifelse(ARTERIAL, Q2arterial, Q2venous)
    return(c(CL = CL, Q2 = Q2, Q3 = Q3, v1 = V1, v2 = V2, v3 = V3))
  }
}

#' R code adapted from nonmem PK file provided in supplementary material of Eleveld et al. Function takes in fixed effect
#' parameter estimates and random effect variance estimates to return parameters for a 3 compartment pk model with an
#' effect site compartment.
#' @param THETA Vector of fixed effects
#' @param ETA Vector of random effect variances
#' @param PATIENT_VARS Named list of patient covariate values
#' @param returnQ Optional logical value to indicate if clearance values should be returned instead of elimination rate constants.
gen_eleveld_pk_pars_nonmem <- function(THETA, ETA, PATIENT_VARS, returnQ = F){

  # convert to logged values - add one when negative
  neg_ix <- which(THETA<0)
  THETA[-neg_ix] <- log(THETA[-neg_ix])
  THETA[neg_ix]  <- log(1+THETA[neg_ix])

  AGE  <- PATIENT_VARS$AGE
  PMA  <- PATIENT_VARS$PMA
  WGT  <- PATIENT_VARS$WGT
  HGT  <- PATIENT_VARS$HGT
  M1F2 <- PATIENT_VARS$M1F2
  TECH <- PATIENT_VARS$TECH
  A1V2 <- PATIENT_VARS$A1V2

  # Al-sallami FFM
  # 1.7 = reference height (m)
  # 70 = reference weight (kg)
  HT2=(HGT/100)*(HGT/100)
  MATM=0.88+((1-0.88)/(1+(AGE/13.4)^(-12.7)))
  MATF=1.11+((1-1.11)/(1+(AGE/7.1)^(-1.1)))
  MATR=0.88+((1-0.88)/(1+(35/13.4)^(-12.7))) # first term in Al-Sallami formula for reference age
  FFMM=MATM*42.92*(HT2)*WGT/(30.93*(HT2)+WGT)
  FFMF=MATF*37.99*(HT2)*WGT/(35.98*(HT2)+WGT)
  FFMR=MATR*42.92*(1.7*1.7)*70/(30.93*(1.7*1.7)+70)
  FFM=FFMM*(2-M1F2) + FFMF*(M1F2-1)
  NFFM=FFM/FFMR

  # maturation
  DV1=1
  DV2=1
  DV3=1

  # sigmoidal maturation of CL
  PMW=PMA*52
  PMR=(35+40/52)*52
  ME50=exp(THETA[8])
  MGAM=exp(THETA[9])
  MCL=(PMW^MGAM)/(PMW^MGAM+ME50^MGAM)
  RCL=(PMR^MGAM)/(PMR^MGAM+ME50^MGAM)
  DCL=MCL/RCL
  DQ2=1

  # sigmoidal maturation of Q3 based on 40 weeks gestation
  PMEW=AGE*52+40
  PMER=35*52+40
  QE50=exp(THETA[14])
  MQ3=PMEW/(PMEW+QE50)
  RQ3=PMER/(PMER+QE50)
  DQ3=MQ3/RQ3

  # aging
  KV1=1
  KV2=exp(THETA[10]*(AGE-35))
  KV3=exp(THETA[13]*(AGE)*(TECH-1))
  KCL=exp(THETA[11]*(AGE)*(TECH-1))
  KQ2=1
  KQ3=1

  # covariate structure
  # V1 scales sigmoid with weight
  VV50=exp(THETA[12])
  CV1=WGT/(WGT+VV50)
  RV1=70/(70+VV50)
  M1 =(CV1/RV1) * KV1 * DV1
  VCV1=(A1V2-1)*(1-CV1)
  V1 =exp(THETA[1]+ETA[1]) * M1 * (1+VCV1*exp(THETA[17]))
  M2 =(WGT/70)^1 * KV2 * DV2
  V2 =exp(THETA[2]+ETA[2]) * M2
  M3 =(NFFM)^1 * KV3 * DV3 # this corresponds to the middle part of the V3 equation
  V3 =exp(THETA[3]+ETA[3]) * M3
  M4 =(WGT/70)^0.75 * KCL * DCL
  CL =exp((2-M1F2)*THETA[4]+(M1F2-1)*THETA[15]+ETA[4]) * M4
  RV2=exp(THETA[2])
  M5 =(V2/RV2)^0.75 * KQ2 * DQ2
  KM5=1+exp(THETA[16])*(1-MQ3)
  Q2 =exp(THETA[5]+ETA[5]+(A1V2-1)*THETA[18]) * M5 * KM5
  RV3=exp(THETA[3])
  M6 =(V3/RV3)^0.75 * KQ3 * DQ3
  Q3 =exp(THETA[6]+ETA[6]) * M6

  K10 = CL/V1
  K12 = Q2/V1
  K21 = Q2/V2
  K13 = Q3/V1
  K31 = Q3/V3

  if(returnQ) return(c(CL = CL, Q2 = Q2, Q3 = Q3, V1 = V1, V2 = V2, V3 = V3))
  else return(c(k10 = K10, k12 = K12, k21 = K21, k13 = K13, k31 = K31, V1 = V1, V2 = V2, V3 = V3))
}


#' Function to generate Emax PD model described by Eleveld et al
#' @param theta Fixed effect parameters on non-logged scale
#' @param eta Random effect
#' @param patient_vars Named list of patient covariate values
gen_eleveld_pd_pars <- function(theta, eta, patient_vars){
  AGE = patient_vars$AGE
  WGT = patient_vars$WGT
  ARTERIAL = ifelse(patient_vars$A1V2 == 1,1,0)
  AGEref=35
  faging <- function(x) exp(x*(AGE-AGEref))
  Ce50 = theta[1]*faging(theta[7])*exp(eta[1])
  ke0 = ifelse(ARTERIAL, theta[2]*(WGT/70)^(-0.25)*exp(eta[2]), theta[8]*(WGT/70)^(-0.25)*exp(eta[2]))
  BISbaseline = theta[3]
  gamma = theta[4] # value of gamma when Ce < Ce50
  gamma2 = theta[9]
  sigma = theta[5]*exp(eta[3])
  bis_delay = 15 + exp(theta[6]*AGE)
  return(c(ke0 = ke0, c50 = Ce50, gamma = gamma, gamma2 = gamma2, E0 = BISbaseline, sigma = sigma, bis_delay = bis_delay))
}


#' Function to return the prior probability for a set of parameters assuming they are log-normally distributed.
#' It is assumed that the last value of of lpr is the prior mean for the variance parameter.
#' @param lpr log parameter values to evaluate
#' @param lhyper hyperparameters for mean and error distributions. List with values "mu", "sig" described below.
#' mu: mean for model parameters and mean of the error distribution
#' sig: variance covariance matrix for model parameters and standard deviation of error distribution
log_prior <- function(lpr, lhyper){
  dmvnorm(lpr, lhyper$mu, lhyper$sig, log = TRUE)
}

#' Function to evaluate the log likelihood given a set of logged parameter values and a set of observed BIS values.
#' @param lpr logged PK-PD-error parameter values
#' @param ivt infusion schedule
#' @param dat data frame with columns c("time","bis") corresponding to observed time and bis values
#' @param ini initial concentrations
log_likelihood <- function(lpr, ivt, dat, gamma, E0, ini = c(0,0,0,0), sum_vals = T){
  epr <- exp(lpr)
  pars_pk <- epr[1:9]; pars_pd <- epr[10]; sig = epr[11]
  sol <- pk_solution_3cpt_metab(pars_pk, ivt = ivt, init = ini)
  con_est <- sol(dat$time)
  con_est[con_est<0] <- 0
  bis_p <- Emax1(pars = pars_pd, ce = con_est[4,], gamma = gamma, E0 = E0)
  # truncated normal distribution
  if(sum_vals) {return(with(dat, sum(log(truncnorm::dtruncnorm(x = bis, mean = bis_p, sd = sig, a = 0, b = 100)))))}
  else{return(with(dat, truncnorm::dtruncnorm(x = bis, mean = bis_p, sd = sig, a = 0, b = 100)))}
}


#' Function to evaluate the negative log posterior given a set of logged parameter values and observed BIS values.
#' @param lpr logged PK-PD-error parameter values
#' @param ivt infusion schedule
#' @param dat data frame with columns corresponding to  observed time and bis values
#' @param lhyper hyperparameter values to be passed to log_prior()
#' @param gamma gamma parameter of PD model (fixed in Eleveld model)
#' @param E0 E0 parameter of PD model (fixed in Eleveld model)
log_posterior_neg <- function(lpr, ivt, dat, lhyper, gamma, E0) {
  dat <- na.omit(dat)
  if(nrow(dat) < 1) {
    -1*log_prior(lpr,lhyper)
  } else {
    -1*(log_prior(lpr,lhyper) + log_likelihood(lpr, ivt, dat, gamma = gamma, E0 = E0))
  }
}





#' Emax function. Assumes maximum value of Emx = 100 and maximum change in value of 100 (i.e. minimum value of 0).
#' @param pars Vector of parameters: (c50,gamma).
#' @param ce Effect-site concentration.
#' @param Emx Maximum value of function.
#' @param E0 Maximum difference between highest and lowest function value.
Emax <- function(pars, ce, Emx = 100, E0 = 100) E0 - Emx*(ce^pars[2] / (ce^pars[2] + pars[1]^pars[2]))


#' Emax function with only c50 to be estimated
#' @param pars c50
#' @param ce Effect site concentrations
#' @param gamma Slope at c50
#' @param E0 Effect at ce = 0
Emax1 <- function(pars, ce, gamma, E0) BISpred <- E0*(1 - ce^gamma/(ce^gamma + pars[1]^gamma))


# Corresponding inverse Emax function
#' @param pars c50
#' @param bis BIS values
#' @param E0 Effect at ce = 0
#' @param gamma Slope at c50
Hinv <- function(pars, bis, E0, gamma) unname(pars[1] * ((E0 - bis) /bis)^(1/gamma))


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


#' Finite difference gradient
#' @param pars parameters of function that gradient is evaluated with respect to.
#' @param fun function to be evaluated.
fdGrad <- function (pars, fun, ...,
                    .relStep = (.Machine$double.eps)^(1/2),
                    minAbsPar = 0) {
  npar <- length(pars)
  incr <- ifelse(abs(pars) <= minAbsPar, .relStep,
                 (abs(pars)-minAbsPar) * .relStep)
  ival <- do.call(fun, list(pars, ...))
  diff <- rep(0,npar)
  sapply(1:npar, function(i) {
    del <- rep(0,npar)
    del[i] <- incr[i]
    (do.call(fun, list(pars+del, ...))-ival)/incr[i]
  })
}

## perform quadrature of multivariate normal

#' Compute Gauss-Hermite quadrature points and weights for a one-dimensional integral.
#' @param points -- number of points
#' @param interlim -- maximum number of Newton-Raphson iterations
hermite <- function (points, z) {
  p1 <- 1/pi^0.4
  p2 <- 0
  for (j in 1:points) {
    p3 <- p2
    p2 <- p1
    p1 <- z * sqrt(2/j) * p2 - sqrt((j - 1)/j) * p3
  }
  pp <- sqrt(2 * points) * p2
  c(p1, pp)
}

gauss.hermite <- function (points, iterlim = 50) {
  x <- w <- rep(0, points)
  m <- (points + 1)/2
  for (i in 1:m) {
    z <- if (i == 1)
      sqrt(2 * points + 1) - 2 * (2 * points + 1)^(-1/6)
    else if (i == 2)
      z - sqrt(points)/z
    else if (i == 3 || i == 4)
      1.9 * z - 0.9 * x[i - 2]
    else 2 * z - x[i - 2]
    for (j in 1:iterlim) {
      z1 <- z
      p <- hermite(points, z)
      z <- z1 - p[1]/p[2]
      if (abs(z - z1) <= 1e-15)
        break
    }
    if (j == iterlim)
      warning("iteration limit exceeded")
    x[points + 1 - i] <- -(x[i] <- z)
    w[i] <- w[points + 1 - i] <- 2/p[2]^2
  }
  r <- cbind(x * sqrt(2), w/sum(w))
  colnames(r) <- c("Points", "Weights")
  r
}


#' Compute multivariate Gaussian quadrature points
#' @param n     - number of points each dimension before pruning
#' @param mu    - mean vector
#' @param sigma - covariance matrix
#' @param prune - NULL - no pruning; (0-1) - fraction to prune
mgauss.hermite <- function(n, mu, sigma, prune=NULL) {
  if(!all(dim(sigma) == length(mu)))
    stop("mu and sigma have nonconformable dimensions")

  dm  <- length(mu)
  gh  <- gauss.hermite(n)
  #idx grows exponentially in n and dm
  idx <- as.matrix(expand.grid(rep(list(1:n),dm)))
  pts <- matrix(gh[idx,1],nrow(idx),dm)
  wts <- apply(matrix(gh[idx,2],nrow(idx),dm), 1, prod)

  ## prune
  if(!is.null(prune)) {
    qwt <- quantile(wts, probs=prune)
    pts <- pts[wts > qwt,]
    wts <- wts[wts > qwt]
  }

  ## rotate, scale, translate points
  eig <- eigen(sigma)
  rot <- eig$vectors %*% diag(sqrt(eig$values))
  pts <- t(rot %*% t(pts) + mu)
  return(list(points=pts, weights=wts))
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




