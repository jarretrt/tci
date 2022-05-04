#' -----------------------------------------------------------------------------
#' Population PK and PK-PD functions -------------------------------------------
#' -----------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Propofol models --------------------------------------------------------------
# ------------------------------------------------------------------------------

#' Marsh population PK model for propofol
#'
#' Evaluates the Marsh propofol model at patient covariates (total body mass) and
#' returns a `pkmod` object. KE0 parameter set to 1.2 in accordance with recommendations
#' from Absalom et al., 2009 "Pharmacokinetic models for propofol- Defining and
#' illuminating the devil in the detail."
#'
#' @param TBW Weight (kg)
#' @param ... Arguments passed to `pkmod`
#' @return `pkmod` object with Marsh population PK parameters
#' @examples
#' pkmod_marsh(TBW = 50)
#' @export
pkmod_marsh <- function(TBW,...){
  K10 = 0.119
  K12 = 0.112
  K13 = 0.042
  K21 = 0.055
  K31 = 0.0033
  KE0 = 1.2
  V1= 0.228 * TBW
  V2= 0.464 * TBW
  V3= 2.89  * TBW
  pkmod(pars_pk = c(V1=V1,V2=V2,V3=V3,K10=K10,K12=K12,K13=K13,K21=K21,K31=K31,KE0=KE0),...)
}


#' @name pkmod_schnider
#' @title Schnider population PK model for propofol
#'
#' @description Evaluate Schnider population PK model at patient covariate values.
#' Published in Schnider et al. (1998). "The influence of method of administration
#' and covariates on the pharmacokinetics of propofol in adult volunteers."
#' Anesthesiology 88 (5):1170-82.
#' @param AGE Age (years)
#' @param LBM Lean body mass (kg). Can be provided instead of TBW, HGT, and MALE
#' @param TBW Weight (kg). Used to calculate LBM if LBM is not provided.
#' @param HGT Height (cm). Used to calculate LBM if LBM is not provided.
#' @param MALE Logical. Used to calculate LBM if LBM is not provided.
#' @param ... Arguments passed to `pkmod`
#' @return `pkmod` object with Schnider population PK parameters
#' @examples
#' pkmod_schnider(AGE = 40,HGT=170,LBM = 43.9)
#' pkmod_schnider(AGE = 40,HGT=170,TBW=50,MALE=TRUE)
#' @export
pkmod_schnider <- function(AGE, HGT, LBM = NULL,TBW=NULL,MALE=NULL,...){

  if(is.null(LBM) & any(is.null(MALE),is.null(TBW),is.null(HGT)))
    stop("MALE, TBW, and HGT must all be provided if LBM is NULL")

  if(is.null(LBM)){
    # calculate lean body mass
    LBM <- ifelse(MALE, 1.1*TBW - 128*(TBW/HGT)^2,
                         1.07*TBW - 148*(TBW/HGT)^2)
  }

  # fixed effect estimates
  theta <- c(4.27,18.9,238,1.89,1.29,0.836,-0.391,0.0456,-0.0681,0.0264,0.024)

  V1  <- theta[1]
  V2  <- theta[2] + theta[7]*(AGE-53)
  V3  <- theta[3]
  CL  <- theta[4] + (LBM-77)*theta[8] + (LBM-59)*(theta[9]) + (HGT-177)*theta[10]
  Q2  <- theta[5] + (theta[11])*(AGE-53)
  Q3  <- 0.836
  KE0 <- 0.456

  # random effect distribution
  cv <- c(4.04,0,14.35,10.05,0,11.79)
  Omega <- diag(log((cv/100)^2+1))
  colnames(Omega) <- c("V1","V2","V3","CL","Q2","Q3")
  pkmod(pars_pk = c(V1=V1,V2=V2,V3=V3,CL=CL,Q2=Q2,Q3=Q3,KE0=KE0), Omega = Omega, ...)
}




#' @name pkmod_eleveld_ppf
#' @title Eleveld population PK model for propofol
#' @description Function takes patient covariate values required for the Eleveld
#' PK or PK-PD model for propofol and returns a `pkmod` object with the appropriate model
#' parameters.
#' @param AGE Age (years)
#' @param TBW Weight (kg)
#' @param HGT Height (cm)
#' @param MALE Sex, logical
#' @param OPIATE Logical indicating presence of opiates. Defaults to TRUE.
#' @param ARTERIAL PK based on arterial sampling rather than venous. Defaults to TRUE.
#' @param PMA Post-menstrual age. Calculated as AGE + 40 weeks if not provided.
#' @param PD Logical. Should PD parameters be returned in addition to PK parameters.
#' @param ... Arguments passed to `pkmod`
#' @return `pkmod` object with Eleveld population PK or PK-PD parameters
#' @examples
#' pkmod_eleveld_ppf(AGE = 40,TBW = 56,HGT=150,MALE = TRUE)
#' @export
pkmod_eleveld_ppf <- function(AGE,TBW,HGT,MALE,OPIATE=TRUE,ARTERIAL=TRUE,PMA=NULL,PD=TRUE,...){

  # fixed effect estimates
  theta <- c(6.28,25.5,273,1.79,1.75,1.11,0.191,42.3,9.06,-0.0156,-0.00286,33.6,-0.0138,68.3,2.10,1.30,1.42,0.68)
  theta_pd <- c(3.08,0.146,93.0,1.47,8.03,0.0517,-0.00635,1.24,1.89)

  if(is.null(PMA)) PMA <- AGE + 40/52
  # BMI uses HGT in meters
  BMI = 10000 * TBW / HGT / HGT
  AGEref=35
  WGTref=70
  HGTref=170
  PMAref = AGEref + 40/52
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
  fAlSallami <- ifelse(MALE, (0.88+(1-0.88)/(1+(AGE/13.4)^(-12.7)))*(9270*TBW)/(6680+216*BMI),
                       (1.11+(1-1.11)/(1+(AGE/7.1)^(-1.1)))*(9270*TBW)/(8780+244*BMI))
  fAlSallami_ref <- (0.88+(1-0.88)/(1+(AGEref/13.4)^(-12.7))) *(9270*WGTref)/(6680+216*BMIref)

  V1arterial <- theta[1]*fcentral(TBW)/fcentral(WGTref)
  V1venous <- V1arterial*(1+theta[17]*(1-fcentral(TBW)))
  V1 = ifelse(ARTERIAL, V1arterial, V1venous)
  V2 <- theta[2]*TBW/WGTref*faging(theta[10])
  V3 <- theta[3]*(fAlSallami/fAlSallami_ref)*fopiates(theta[13])

  CL <- (MALE*theta[4] + (1-MALE)*theta[15])*(TBW/WGTref)^0.75*(fCLmaturation/fCLmaturation_ref)*fopiates(theta[11])
  Q2arterial <- theta[5]*(V2/V2ref)^(0.75)*(1+theta[16]*(1-fQ3maturation))
  Q2venous <- Q2arterial*theta[18]
  Q2 <- ifelse(ARTERIAL, Q2arterial, Q2venous)
  Q3 <- theta[6]*(V3/V3ref)^(0.75)*(fQ3maturation/fQ3maturation_ref)
  # include KE0 with PK parameters for effect-site targeting
  KE0 <- ifelse(ARTERIAL, theta_pd[2]*(TBW/70)^(-0.25), theta_pd[8]*(TBW/70)^(-0.25))

  # PD parameters
  C50   <- theta_pd[1]*faging(theta_pd[7])
  BIS0   <- theta_pd[3]
  GAMMA <- theta_pd[4] # value of gamma when Ce < Ce50
  GAMMA2 <- theta_pd[9]

  # random effect variances
  omega_pk <- c(0.610,0.565,0.597,0.265,0.346,0.209,0.702,0.463)
  omega_pd <- c(0.242,0.230)

  if(PD){
    # drop random effect on residual PK error
    Omega <- diag(c(omega_pk[-8],omega_pd))
    colnames(Omega) <- c("V1","V2","V3","CL","Q2","Q3","KE0","C50","sigma_add")
  } else{
    Omega <- diag(omega_pk)
    colnames(Omega) <- c("V1","V2","V3","CL","Q2","Q3","KE0","sigma_add")
  }

  # residual error
  sig <- ifelse(PD, 8.03, 0.191)

  # logged response
  lr <- ifelse(PD, FALSE, TRUE)

  if(PD){
    pkmod(pars_pk = c(V1=V1,V2=V2,V3=V3,CL=CL,Q2=Q2,Q3=Q3,KE0=KE0),
          pars_pd = c(C50=C50,E0=BIS0,EMX=BIS0,GAMMA=GAMMA,GAMMA2=GAMMA2),
          pdfn = emax_eleveld,
          pdinv = emax_inv_eleveld,
          sigma_add = sig,
          log_response = lr,
          Omega = Omega, ...)
  } else{
    pkmod(pars_pk = c(V1=V1,V2=V2,V3=V3,CL=CL,Q2=Q2,Q3=Q3,KE0=KE0),
          sigma_add = sig,
          log_response = lr,
          Omega = Omega, ...)
  }
}


# ------------------------------------------------------------------------------
# Remifentanil models -----------------------------------------------------------
# ------------------------------------------------------------------------------


#' @name pkmod_minto
#' @title Minto population PK model for remifentanil
#'
#' @description Evaluate Minto population PK model at patient covariate values.
#' Published in Minto et al. (1997). "Influence of Age and Gender on the
#' Pharmacokinetics and Pharmacodynamics of Remifentanil: I. Model Development."
#' Anesthesiology 86:10-23 doi: https://doi.org/10.1097/00000542-199701000-00004.
#' Residual error standard deviations are taken from Eleveld et al. (2017).
#' "An Allometric Model of Remifentanil Pharmacokinetics and Pharmacodynamics."
#' Anesthesiology Vol. 126, 1005–1018 doi: https://doi.org/10.1097/ALN.0000000000001634.
#' @param AGE Age (years)
#' @param LBM Lean body mass (kg). Can be provided instead of TBW, HGT, and MALE
#' @param TBW Weight (kg). Used to calculate LBM if LBM is not provided.
#' @param HGT Height (cm). Used to calculate LBM if LBM is not provided.
#' @param MALE Logical. Used to calculate LBM if LBM is not provided.
#' @param PD Logical. Should PD parameters be returned in addition to PK parameters.
#' @param ... Arguments passed to `pkmod`
#' @return `pkmod` object with Schnider population PK parameters
#' @examples
#' pkmod_minto(AGE = 40,HGT=170,LBM = 43.9)
#' pkmod_minto(AGE = 40,HGT=170,TBW=50,MALE=TRUE)
#' @export
pkmod_minto <- function(AGE, HGT, LBM = NULL,TBW=NULL,MALE=NULL, PD = TRUE, ...){

  if(is.null(LBM) & any(is.null(MALE),is.null(TBW),is.null(HGT)))
    stop("MALE, TBW, and HGT must all be provided if LBM is NULL")

  if(is.null(LBM)){
    # calculate lean body mass
    LBM <- ifelse(MALE, 1.1*TBW - 128*(TBW/HGT)^2,
                  1.07*TBW - 148*(TBW/HGT)^2)
  }

  # fixed effect estimates for PK model
  theta <- c(5.1,9.82,5.42,2.6,2.05,0.076,
             -0.0201,-0.0811,-0.0162,-0.0301,-0.00113,
             0.072,0.108,0.0191)

  V1  <- theta[1] + theta[7]*(AGE-40) + theta[12]*(LBM-55)
  V2  <- theta[2] + theta[8]*(AGE-40) + theta[13]*(LBM-55)
  V3  <- theta[3]
  CL  <- theta[4] + theta[9]*(AGE-40) + theta[14]*(LBM-55)
  Q2  <- theta[5] + theta[10]*(AGE-40)
  Q3  <- theta[6] + theta[11]*(AGE-40)

  # fixed effect estimates for PD model
  theta_pd <- c(0.595,20,5.5,13.1,2.44,-0.007,-0.148)
  KE0   <- theta_pd[1] + theta_pd[6] * (AGE-40)
  E0    <- theta_pd[2]
  EMX   <- theta_pd[3]
  C50  <- theta_pd[4] + theta_pd[7] * (AGE-40)
  GAMMA <- theta_pd[5]

  # random effect distribution
  omega_pk <- (c(26,29,66,14,36,41)/100)^2
  omega_pd <- (c(60,15,30,45,54)/100)^2

  if(PD){

    Omega <- diag(c(omega_pk,omega_pd))
    colnames(Omega) <- c("V1","V2","V3","CL","Q2","Q3","KE0","E0","EMX","C50","GAMMA")

    pkmod(pars_pk = c(V1=V1,V2=V2,V3=V3,CL=CL,Q2=Q2,Q3=Q3,KE0=KE0),
          pars_pd = c(E0=E0,EMX=EMX,C50=C50,GAMMA=GAMMA),
          pdfn = emax,
          pdinv = emax_inv,
          sigma_add = 1.96,
          log_response = FALSE,
          Omega = Omega, ...)
  } else{

    Omega <- diag(c(omega_pk,omega_pd[1]))
    colnames(Omega) <- c("V1","V2","V3","CL","Q2","Q3","KE0")

    pkmod(pars_pk = c(V1=V1,V2=V2,V3=V3,CL=CL,Q2=Q2,Q3=Q3,KE0=KE0),
          sigma_add = 0.111,
          log_response = TRUE,
          Omega = Omega, ...)
  }
}



#' @name pkmod_kim
#' @title Kim population PK model for remifentanil
#'
#' @description Evaluate Kim population PK model at patient covariate values.
#' Published in Kim et al. (2017). "Disposition of Remifentanil in Obesity: A New
#' Pharmacokinetic Model Incorporating the Influence of Body Mass"
#' Anesthesiology Vol. 126, 1019–1032. doi: https://doi.org/10.1097/ALN.0000000000001635
#' @param AGE Age (years)
#' @param TBW Total body weight (kg).
#' @param BMI Body mass index
#' @param MALE Logical. Used to calculate LBM if LBM is not provided.
#' @param FFM Optional specification of
#' @param ... Arguments passed to `pkmod`
#' @return `pkmod` object with Schnider population PK parameters
#' @examples
#' pkmod_kim(AGE = 40,TBW = 75, BMI = 30, MALE = TRUE)
#' pkmod_kim(AGE = 40,TBW = 75, FFM = 52.83)
#' @export
pkmod_kim <- function(AGE,TBW,BMI=NULL,MALE=NULL,FFM=NULL,...){

  if((is.null(BMI)|is.null(MALE)) & is.null(FFM))
    stop("BMI and MALE must be specified if FFM is NULL")

  if(is.null(FFM)){
    # calculate fat-free mass
    FFM <- ifelse(MALE, 9.27*10^3*TBW / (6.68*10^3 + 216*BMI),
                  9.27*10^3*TBW / (8.78*10^3 + 244*BMI))
  }

  # fixed effect estimates for PK model
  theta <- c(4.76,8.4,4,2.77,1.94,0.197,NA,NA,
             0.658,0.573,0.0936,0.0477,0.336,0.0149,0.0280)

  V1  <- theta[1] * (TBW/74.5)^theta[9]
  V2  <- theta[2] * (FFM/52.3)^theta[10] - theta[11]*(AGE-37)
  V3  <- theta[3] - theta[12]*(AGE-37)
  CL  <- theta[4] * (TBW/74.5)^theta[13] - theta[14]*(AGE-37)
  Q2  <- theta[5] - theta[15]*(AGE-37)
  Q3  <- theta[6]

  # random effect distribution
  cv <- c(53.6,56.3,59,21.2,41.2,59)
  Omega <- diag(log((cv/100)^2+1))
  colnames(Omega) <- c("V1","V2","V3","CL","Q2","Q3")

  pkmod(pars_pk = c(V1=V1,V2=V2,V3=V3,CL=CL,Q2=Q2,Q3=Q3),
        sigma_add = 0.0217,
        sigma_mult = 0.168,
        log_response = FALSE,
        Omega = Omega, ...)
}


#' @name pkmod_eleveld_remi
#' @title Eleveld population PK model for remifentanil
#' @description Function takes patient covariate values required for the Eleveld
#' PK or PK-PD model for propofol and returns a `pkmod` object with the appropriate model
#' parameters.
#' @param AGE Age (years)
#' @param TBW Total body weight (kg)
#' @param HGT Height (cm)
#' @param MALE Sex, logical
#' @param OPIATE Logical indicating presence of opiates. Defaults to TRUE.
#' @param ARTERIAL PK based on arterial sampling rather than venous. Defaults to TRUE.
#' @param PMA Post-menstrual age. Calculated as AGE + 40 weeks if not provided.
#' @param PD Logical. Should PD parameters be returned in addition to PK parameters.
#' @param ... Arguments passed to `pkmod`
#' @return `pkmod` object with Eleveld population PK or PK-PD parameters
#' @examples
#' pkmod_eleveld_remi(AGE = 40,TBW = 56,HGT=150,MALE = TRUE)
#' @export
pkmod_eleveld_remi <- function(AGE, MALE, TBW = NULL,HGT=NULL,BMI = NULL, FFM = NULL,PD=TRUE,...){

  if(is.null(BMI) & is.null(TBW) & is.null(HGT))
    stop("TBW and HGT must be specified if BMI is NULL")

  if(is.null(BMI)) BMI = 10000 * TBW / HGT / HGT

  AGEref=35
  TBWref=70
  HGTref=170
  BMIref = 10000 * TBWref / HGTref / HGTref

  # calculate FFM
  FFM <- ifelse(MALE, (0.88+(1-0.88)/(1+(AGE/13.4)^(-12.7)))*(9270*TBW)/(6680+216*BMI),
                       (1.11+(1-1.11)/(1+(AGE/7.1)^(-1.1)))*(9270*TBW)/(8780+244*BMI))
  FFMref <- (0.88+(1-0.88)/(1+(AGEref/13.4)^(-12.7))) *(9270*TBWref)/(6680+216*BMIref)

  faging <- function(x) exp(x*(AGE-AGEref))
  fsigmoid <- function(x,E50,lambda) x^lambda / (x^lambda + E50^lambda)

  # fixed effect estimates
  theta <- c(2.88,-0.00554,-0.00327,-0.0315,0.470,-0.0260)
  theta_pd <- c(3.08,0.146,93.0,1.47,8.03,0.0517,-0.00635,1.24,1.89)

  V1ref = 5.81
  V2ref = 8.82
  V3ref = 5.03
  CLref = 2.58
  Q2ref = 1.72
  Q3ref = 0.124

  SIZE = FFM/FFMref
  KMAT = fsigmoid(TBW,theta[1],2)
  KMATref = fsigmoid(TBWref,theta[1],2)
  KSEX = ifelse(MALE, 1, 1+theta[5]*fsigmoid(AGE,12,6)*(1-fsigmoid(AGE,45,6)))

  V1 = V1ref * SIZE * faging(theta[2])
  V2 = V2ref * SIZE * faging(theta[3])*KSEX
  V3 = V2ref * SIZE * faging(theta[4])*exp(theta[6]*(TBW-70))
  CL = CLref * SIZE^0.75 * (KMAT/KMATref) * KSEX * faging(theta[3])
  Q2 = Q2ref * (V2/V2ref)^0.75 * faging(theta[2]) * KSEX
  Q3 = Q3ref * (V3/V3ref)^0.75 * faging(theta[2])

  # PD model
  KE0 = 1.09 * faging(-0.0289)
  E0 = 19.9
  EMX = 5.66
  C50 = 12.7
  GAMMA = 2.87

  # random effect variances
  omega_pk <- c(0.104,0.115,0.810,0.0197,0.0547,0.285)
  omega_pd <- c(0.021,0.101,0.391,0.379,0.947)

  if(PD){
    # drop random effect on residual PK error
    Omega <- diag(c(omega_pk,omega_pd))
    colnames(Omega) <- c("V1","V2","V3","CL","Q2","Q3","E0","EMX","C50","GAMMA","KE0")

    pkmod(pars_pk = c(V1=V1,V2=V2,V3=V3,CL=CL,Q2=Q2,Q3=Q3,KE0=KE0),
          pars_pd = c(E0=E0,EMX=EMX,C50=C50,GAMMA=GAMMA),
          pdfn = emax_remi,
          pdinv = emax_inv_remi,
          sigma_add = 1.96,
          log_response = FALSE,
          Omega = Omega, ...)

  } else{
    Omega <- diag(c(omega_pk,omega_pd[5]))
    colnames(Omega) <- c("V1","V2","V3","CL","Q2","Q3","KE0")

    pkmod(pars_pk = c(V1=V1,V2=V2,V3=V3,CL=CL,Q2=Q2,Q3=Q3,KE0=KE0),
          sigma_add = 1.11,
          log_response = TRUE,
          Omega = Omega, ...)
  }
}

# ------------------------------------------------------------------------------
# Extra functions --------------------------------------------------------------
# ------------------------------------------------------------------------------

#' Sample parameters from a population PK model
#'
#' @param pkmod `pkmod` object with associated Omega matrix describing random effect variances
#' @param log_normal Logical. Assumes random effects are log-normally distributed
#' and multiplicative if TRUE, additive and normally distributed if FALSE.
#' @param ... Arguments passed to update.pkmod
#' @examples
#' sample_pkmod(pkmod_schnider(AGE = 40,HGT=170,TBW=50,MALE=TRUE))
#' sample_pkmod(pkmod_eleveld(AGE = 40,TBW = 56,HGT=150,MALE = TRUE, PD = FALSE))
#' @export
sample_pkmod <- function(pkmod, log_normal = TRUE, ...){

  pkmod <- update(pkmod, ...)
  if(is.null(pkmod$Omega)) stop("pkmod must have an Omega matrix")

  eta <- c(mvtnorm::rmvnorm(1, mean = rep(0, length(diag(pkmod$Omega))), sigma = pkmod$Omega))
  names(eta) <- colnames(pkmod$Omega)
  pknms <- names(pkmod$pars_pk)[names(pkmod$pars_pk) %in% names(eta)]
  pdnms <- names(pkmod$pars_pd)[names(pkmod$pars_pd) %in% names(eta)]

  if(log_normal){
    pars_pk_new <- pkmod$pars_pk[pknms] * exp(eta[pknms])
    pars_pd_new <- pkmod$pars_pd[pdnms] * exp(eta[pdnms])
    sigma_add_new <- ifelse("sigma_add" %in% names(eta),
                            pkmod$sigma_add * exp(eta["sigma_add"]),
                            pkmod$sigma_add)

    sigma_mult_new <- ifelse("sigma_mult" %in% names(eta),
                            pkmod$sigma_mult * exp(eta["sigma_mult"]),
                            pkmod$sigma_mult)
  } else{
    pars_pk_new <- pkmod$pars_pk[pknms] + eta[pknms]
    pars_pd_new <- pkmod$pars_pd[pdnms] + eta[pdnms]
    sigma_add_new <- ifelse("sigma_add" %in% names(eta),
                            pkmod$sigma_add + eta["sigma_add"],
                            pkmod$sigma_add)

    sigma_mult_new <- ifelse("sigma_mult" %in% names(eta),
                             pkmod$sigma_mult + eta["sigma_mult"],
                             pkmod$sigma_mult)
  }

  update(pkmod, pars_pk = pars_pk_new, pars_pd = pars_pd_new,
         sigma_add = sigma_add_new, sigma_mult = sigma_mult_new)
}




