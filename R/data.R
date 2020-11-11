#' Eleveld et al. pharmacokinetic data
#'
#' Empirical Bayes (EB) estimates of PK parameters made by the Eleveld et al (2018) PK-PD model.
#'
#' @docType data
#'
#' @usage data(eleveld_pk)
#'
#' @references Eleveld et al. (2018) British Journal of Anesthesia Vol. 120, 5:942-959
#' (\href{https://bjanaesthesia.org/article/S0007-0912(18)30051-5/abstract}{BJA})
#'
#' @format A data frame with 1033 rows and 16 variables:
#' \describe{
#'   \item{ID}{Patient ID}
#'   \item{V1}{EB estimate of first compartment volume}
#'   \item{V2}{EB estimate of second compartment volume}
#'   \item{V3}{EB estimate of third compartment volume}
#'   \item{CL}{EB estimate of clearance for the first compartment}
#'   \item{Q2}{EB estimate of inter-compartmental clearance for second compartment}
#'   \item{Q3}{EB estimate of inter-compartmental clearance for third compartment}
#'   \item{AGE}{Patient's age (years)}
#'   \item{WGT}{Patient's weight (kg)}
#'   \item{HGT}{Patient's height (cm)}
#'   \item{M1F2}{Patient's sex: male = 1, female = 2}
#'   \item{PMA}{Patient's post-menstrual age. Assumed to be age + 40 weeks if not provided}
#'   \item{TECH}{Presence of concomitant anaesthetic techniques (Local anesthetic = 1, Opioids = 2)}
#'   \item{BMI}{Patient's BMI}
#'   \item{FFM}{Patient's fat-free mass (FFM)}
#'   \item{A1V2}{Sampling site: arterial sampling = 1, venous sampling = 2}
#' }
"eleveld_pk"




#' Eleveld et al. pharmacodynamic data
#'
#' Empirical Bayes (EB) estimates of PD parameters made by the Eleveld et al (2018) PK-PD model.
#'
#' @docType data
#'
#' @usage data(eleveld_pd)
#'
#' @references Eleveld et al. (2018) British Journal of Anesthesia Vol. 120, 5:942-959
#' (\href{https://bjanaesthesia.org/article/S0007-0912(18)30051-5/abstract}{BJA})
#'
#' @format A data frame with 122 rows and 15 variables:
#' \describe{
#'   \item{ID}{Patient ID}
#'   \item{E50}{EB estimate of effect-site concentration required to achieve 50 percent response}
#'   \item{KE0}{EB estimate of elimination rate from effect-site compartment}
#'   \item{Emax}{EB estimate of baseline bispectral index (BIS) with no drug administered}
#'   \item{GAM}{EB estimate of Hill parameter when the effect-site concentration is less than E50}
#'   \item{GAM1}{EB estimate of Hill parameter when the effect-site concentration is greater than than E50}
#'   \item{RESD}{EB estimate of residual error term}
#'   \item{ALAG1}{Estimated time lag in BIS measurements due to patient age (fixed-effects only)}
#'   \item{AGE}{Patient's age (years)}
#'   \item{WGT}{Patient's weight (kg)}
#'   \item{HGT}{Patient's height (cm)}
#'   \item{M1F2}{Patient's sex: male = 1, female = 2}
#'   \item{A1V2}{Sampling site: arterial sampling = 1, venous sampling = 2}
#'   \item{PMA}{Patient's post-menstrual age. Assumed to be age + 40 weeks if not provided}
#'   \item{TECH}{Presence of concomitant anaesthetic techniques (Local anesthetic = 1, Opioids = 2)}
#' }
"eleveld_pd"
