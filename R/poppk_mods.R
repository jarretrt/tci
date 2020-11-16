#' ------------------------------------------------------------------------------------
#' Population PK and PK-PD functions --------------------------------------------------
#' ------------------------------------------------------------------------------------

#' Marsh population PK model.
#'
#' Takes in a vector of patient weights and returns a data frame of patient PK-PD parameters.
#' KE0 parameter set to 1.2 in accordance with Base Primea system and recommendations
#' from Absalom et al. (2009) "Pharmacokinetic models for propofol- Defining and
#' illuminating the devil in the detail"
#'
#' @param df data frame with column titled "TBM" giving patient total body
#' mass in kg.
#' @param rate Logical. Should elimination rate constants be returned instead
#' of clearance parameters.
#' @export
marsh_poppk <- function(df, rate = TRUE){
  if(!("TBM" %in% names(df))) stop('The data frame must have a column named "TBM"')

  df$V1  = 0.228 * df$TBM
  df$V2  = 0.464 * df$TBM
  df$V3  = 2.89  * df$TBM

  K10 = 0.119
  K12 = 0.112
  K13 = 0.042
  K21 = 0.055
  K31 = 0.0033

  if(rate){
    df$K10 = K10
    df$K12 = K12
    df$K21 = K21
    df$K13 = K13
    df$K31 = K31
  } else{
    df$CL  <- K10 * df$V1
    df$Q12 <- K12 * df$V1
    df$Q21 <- K21 * df$V2
    df$Q13 <- K13 * df$V1
    df$Q31 <- K31 * df$V3
  }

  df$KE0 = 1.2

  class(df) <- c(class(df),"poppk")
  return(df)
}
#' @examples
#' marsh_poppk(data.frame(TBM = c(50,70,90)))


#' Schnider population PK model
#'
#' @name schnider_poppk
#' @title schnider_poppk
#' @param df data frame with variable names "AGE","TBM","HGT","MALE"
#' @param rate Logical. Should rate parameters be returned rather than clearance.
#' Defaults to FALSE
#' @param rand Logical. Should a vector of Monte Carlo samples be returned instead
#' of point estimates at patient covariate values. Defaults to FALSE.
#' @export
schnider_poppk <- function(df, rate = FALSE, rand = FALSE){

  covar_names <- c("AGE","TBM","HGT","MALE")
  if(!all(covar_names %in% names(df))) stop(
    paste('The data frame must have a columns named',
          paste(covar_names, collapse = ", ")))

  # fixed effect estimates
  theta0 <- c(4.27,18.9,238,1.89,1.29,0.836,-0.391,0.0456,-0.0681,0.0264,0.024)
  # random effect variances
  omega <- c(0.278,2.330,34.900,0.059,0.112,0.044,0.070,0.009,0.017,0.009,0.005)

  if(rand){
    # sample random effects
    eta <- mvtnorm::rmvnorm(nrow(df), rep(0,length(omega)), sigma = diag(omega^2))
    theta <- t(theta0 * t(exp(eta)))
  } else{
    theta <- matrix(theta0, nrow = 1)
  }

  # caclulate lean body mass
  df$LBM <- ifelse(df$MALE, 1.1*df$TBM - 128*(df$TBM/df$HGT)^2, 1.07*df$TBM - 148*(df$TBM/df$HGT)^2)

  df$V1  = theta[,1]
  df$V2  = theta[,2] + theta[,7]*(df$AGE-53)
  df$V3  = theta[,3]
  CL <- theta[,4] + (df$TBM-77)*theta[,8] + (df$LBM-59)*(theta[,9]) + (df$HGT-177)*theta[,10]
  Q2 <- theta[,5] + (theta[,11])*(df$AGE-53)
  Q3 <- 0.836

  if(rate){
    df$K10 <- CL / df$V1
    df$K12 <- Q2 / df$V1
    df$K21 <- Q2 / df$V2
    df$K13 <- Q3 / df$V1
    df$K31 <- Q3 / df$V3
  } else{
    df$CL <- CL
    df$Q2 <- Q2
    df$Q3 <- Q3
  }

  df$KE0 <- 0.456

  class(df) <- c(class(df),"poppk")
  return(df)
}
#' @examples
#' dat <- data.frame(AGE  = c(20,40,65),
#'                   TBM  = c(50,70,90),
#'                   HGT  = c(150,170,200),
#'                   MALE = c(TRUE,FALSE,TRUE))
#' schnider_poppk(dat, rand = F, rate = F)
#' schnider_poppk(dat, rand = T, rate = T)


#' Eleveld et al. population PK model.
#'
#' Function takes a data frame of patient covariate values with variable names
#' "AGE","TBM","HGT","MALE" and returns PK parameter values.
#'
#' @name eleveld_poppk
#' @title eleveld_poppk
#'
#' @param df Data frame with variable names "AGE","TBM","HGT","MALE"
#' @param PD Logical. Should PD parameters be returned in addition to PK parameters.
#' Defaults to TRUE.
#' @param rate Logical. Should rate parameters be returned rather than clearance.
#' Defaults to FALSE
#' @param rand Logical. Should a vector of Monte Carlo samples be returned instead
#' of point estimates at patient covariate values. Defaults to FALSE.
#' @export
#'
eleveld_poppk <- function(df, PD = TRUE, rate = FALSE, rand = FALSE){

  # fixed effect estimates
  theta <- c(6.28,25.5,273,1.79,1.75,1.11,0.191,42.3,9.06,-0.0156,-0.00286,33.6,-0.0138,68.3,2.10,1.30,1.42,0.68)
  # random effect variances
  omega <- c(0.610,0.565,0.597,0.265,0.346,0.209,0.463)
  # simulate random effects if rand = TRUE
  if(rand){
    eta <- mvtnorm::rmvnorm(n=nrow(df), mean = rep(0,length(omega)), sigma=diag(omega))
  } else{
    eta <- matrix(0, nrow = nrow(df), ncol = length(omega))
  }

  AGE  <- df$AGE
  PMA  <- df$PMA
  WGT  <- df$WGT
  HGT  <- df$HGT
  M1F2 <- df$M1F2
  TECH <- df$TECH
  A1V2 <- df$A1V2
  MALE <- ifelse(M1F2 == 1, 1, 0)
  OPIATE <- ifelse(TECH == 2, 1, 0)
  ARTERIAL <- ifelse(A1V2 == 1, 1, 0)
  AGEref=35
  WGTref=70
  HGTref=170

  df$PMA[is.na(PMA)] <- AGE[is.na(PMA)] + 40/52
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

  fAlSallami <- ifelse(MALE, (0.88+(1-0.88)/(1+(AGE/13.4)^(-12.7)))*(9270*WGT)/(6680+216*BMI),
                           (1.11+(1-1.11)/(1+(AGE/7.1)^(-1.1)))*(9270*WGT)/(8780+244*BMI))
  fAlSallami_ref <- (0.88+(1-0.88)/(1+(AGEref/13.4)^(-12.7))) *(9270*WGTref)/(6680+216*BMIref)

  V1arterial <- theta[1]*fcentral(WGT)/fcentral(WGTref)*exp(eta[,1])
  V1venous <- V1arterial*(1+theta[17]*(1-fcentral(WGT)))
  V2 <- theta[2]*WGT/WGTref*faging(theta[10])*exp(eta[,2])
  V3 <- theta[3]*(fAlSallami/fAlSallami_ref)*fopiates(theta[13])*exp(eta[,3])

  CL <- (MALE*theta[4] + (1-MALE)*theta[15])*(WGT/WGTref)^0.75*(fCLmaturation/fCLmaturation_ref)*fopiates(theta[11])*exp(eta[,4])
  Q2arterial <- theta[5]*(V2/V2ref)^(0.75)*(1+theta[16]*(1-fQ3maturation))*exp(eta[,5])
  Q2venous <- Q2arterial*theta[18]
  Q3 <- theta[6]*(V3/V3ref)^(0.75)*(fQ3maturation/fQ3maturation_ref)*exp(eta[,6])

  V1 = ifelse(ARTERIAL, V1arterial, V1venous)
  K10 = CL/V1
  K12 = ifelse(ARTERIAL, Q2arterial/V1, Q2venous/V1)
  K21 = ifelse(ARTERIAL, Q2arterial/V2, Q2venous/V2)
  K13 = Q3/V1
  K31 = Q3/V3

  if(rate){
    df$K10 <- K10
    df$K12 <- K12
    df$K21 <- K21
    df$K13 <- K13
    df$K31 <- K31
    df$V1  <- V1
    df$V2  <- V2
    df$V3  <- V3
  } else{
    df$CL <- CL
    df$Q2 <- ifelse(ARTERIAL, Q2arterial, Q2venous)
    df$Q3 <- Q3
    df$V1  <- V1
    df$V2  <- V2
    df$V3  <- V3
  }

  df$LN_SIGMA <- theta[7]*exp(eta[,7])

  if(PD){
    theta_pd <- c(3.08,0.146,93.0,1.47,8.03,0.0517,-0.00635,1.24,1.89)
    omega_pd <- c(0.242,0.702,0.230)

    if(rand){
      eta_pd <- mvtnorm::rmvnorm(n=nrow(df), mean = rep(0,length(omega_pd)), sigma=diag(omega_pd))
    } else{
      eta_pd <- matrix(0, nrow = nrow(df), ncol = length(omega_pd))
    }

    df$CE50 = theta_pd[1]*faging(theta_pd[7])*exp(eta_pd[,1])
    df$KE0 = ifelse(ARTERIAL, theta_pd[2]*(WGT/70)^(-0.25)*exp(eta[2]), theta_pd[8]*(WGT/70)^(-0.25)*exp(eta_pd[,2]))
    df$BIS0 = theta_pd[3]
    df$GAMMA = theta_pd[4] # value of gamma when Ce < Ce50
    df$GAMMA2 = theta_pd[9]
    df$SIGMA = theta_pd[5]*exp(eta_pd[,3])
    df$BIS_DELAY = 15 + exp(theta_pd[6]*AGE)
  }

  class(df) <- c(class(df),"poppk")
  return(df)
}
#' @examples
#' data(eleveld_pk)
#' df <- eleveld_pk[1:10,c("ID","AGE","WGT","HGT","M1F2","PMA","TECH","BMI","FFM","A1V2")]
#' eleveld_poppk(df, rate = TRUE)
#' eleveld_poppk(df, rate = FALSE, PD = TRUE, rand = TRUE)





