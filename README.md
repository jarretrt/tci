README
================

# tci

The `tci` package implements plasma- and effect-site targeting
target-controlled infusion (TCI) algorithms to 1-, 2-, 3-, and
3-compartment/effect-site pharmacokinetic (PK) models with intravenous
drug administration. TCI algorithms can further be applied to
pharmacodynamic (PD) targets when a PD model and its inverse are
specified. Functions are supplied for simulation of patient responses to
TCI infusion schedules under open- or Bayesian closed-loop control. See
the `overview` vignette for further details. Custom user-defined PK
models or TCI algorithms can be specified, as illustrated in the
`custom` vignette.

## Installation

The `tci` package can be installed from
[CRAN](https://cran.rstudio.com/web/packages/tci/index.html) using the
command.

``` r
install.packages("tci")
```

The most recent version can be downloaded from GitHub using the
`devtools` package and loaded as follows.

``` r
devtools::install_github("jarretrt/tci")
```

``` r
library(tci)
library(ggplot2) # for plotting
```

## Examples

### PK/PKPD Models

PK and PKPD models are created using the function `pkmod`. The only
required argument to `pkmod` is a named vector `pars_pk`, indicating the
PK parameters for a specific individual. Acceptable PK parameter names
can be found by calling `list_parnms()`. The number of compartments will
be inferred from the parameter names.

``` r
# 3-compartment PK model with effect site
(mod3ecpt <- pkmod(pars_pk = c(cl = 10, q2 = 2, q3 =20, v = 15, v2 = 30, v3 = 50, ke0 = 1.2)))
```

    ## --- PK model ------------------------------------------- 
    ## 4-compartment PK model 
    ## PK parameters: cl = 10, q2 = 2, q3 = 20, v = 15, v2 = 30, v3 = 50, ke0 = 1.2 
    ## Initial concentrations: (0,0,0,0) 
    ## Plasma compartment: 1 
    ## Effect compartment: 4 
    ## --- Simulation ----------------------------------------- 
    ## Additive error SD: 0 
    ## Multiplicative error SD: 0 
    ## Logged response: FALSE

A PD model can be added by specifying arguments `pdfn`, `pdinv`, and
`pars_pd`, with the first to referring to functions implementing the PD
function and its inverse, and `pars_pd` specifying a named vector of PD
parameters used by `pdfn` and `pdinv`. Arguments can be added to an
existing `pkmod` object using the `update.pkmod` method.

``` r
# 3-compartment PK model with effect site and Emax PD 
(mod3ecpt_pd <- update(mod3ecpt, pdfn = emax, pdinv = emax_inv, 
                 pars_pd = c(e0 = 100, emx = 100, c50 = 3.5, gamma = 2.2)))
```

    ## --- PK model ------------------------------------------- 
    ## 4-compartment PK model 
    ## PK parameters: cl = 10, q2 = 2, q3 = 20, v = 15, v2 = 30, v3 = 50, ke0 = 1.2 
    ## Initial concentrations: (0,0,0,0) 
    ## Plasma compartment: 1 
    ## Effect compartment: 4 
    ## --- PD model ------------------------------------------- 
    ## PD parameters: e0 = 100, emx = 100, c50 = 3.5, gamma = 2.2 
    ## --- Simulation ----------------------------------------- 
    ## Additive error SD: 0 
    ## Multiplicative error SD: 0 
    ## Logged response: FALSE

## Population PK models

Several published population PK models are currently implemented in
`tci`. For propofol, these include the Marsh, Schnider, and Eleveld
models. For remifentanil, they include the Minto, Kim, and Eleveld
models. Each population PK model function takes patient covariates
(e.g., age, total body weight) and returns a `pkmod` object that
implements the population model at the corresponding structural PK/PK-PD
parameter values.

The function `list_pkmods` prints the available models.

``` r
list_pkmods()
```

| Population.model      | Function               | Drug         | Type    | Required.covariates              |
|:----------------------|:-----------------------|:-------------|:--------|:---------------------------------|
| Marsh                 | pkmod\_marsh()         | Propofol     | PK      | TBW                              |
| Schnider              | pkmod\_schnider()      | Propofol     | PK      | AGE, HGT, LBM or (TBW and MALE)  |
| Eleveld (propofol)    | pkmod\_eleveld\_ppf()  | Propofol     | PK/PKPD | AGE, HGT, MALE, TBW              |
| Minto                 | pkmod\_minto()         | Remifentanil | PK/PKPD | AGE, LBM or (MALE, TBW, and HGT) |
| Kim                   | pkmod\_kim()           | Remifentanil | PK      | AGE, TBW, FFM or (MALE and BMI)  |
| Eleveld (remifentanil | pkmod\_eleveld\_remi() | Remifentanil | PK/PKPD | AGE, MALE, BMI or (TBW and HGT)  |

By default, calling population PK model functions will return a model
with parameter values corresponding to the point estimates for the set
of covariates. Parameter values can be sampled, however, by calling the
function `sample_pkmod` on any `pkmod` object with an `Omega` matrix
describing inter-individual variability.

``` r
# pkmod at point estimates for Eleveld propofol model
mod_elvd <- pkmod_eleveld_ppf(AGE = 40,TBW = 56,HGT=150,MALE = TRUE)

# sample parameter values for new simulated individual
set.seed(1)
mod_elvd_sample <- sample_pkmod(mod_elvd)
```

## Infusion schedules

An infusion schedule is required to evaluate a `pkmod` object. This
schedule should be a matrix with column labels “begin,” “end,” and
“infrt,” indicating infusion begin times, end times, and infusion rates.
It can be created directly by the user, or outputted by the `create_inf`
or `apply_tci` functions. In the former function, the user specifies
infusion start times, durations, and infusion rates.

``` r
# Infusions of 100 mg/min for 30 sec at 0, 3, 6 min
(multi_inf <- create_inf(times = c(0,3,6),duration = 0.5,infrt = 100))
```

    ##      begin end infrt
    ## [1,]   0.0 0.5   100
    ## [2,]   0.5 3.0     0
    ## [3,]   3.0 3.5   100
    ## [4,]   3.5 6.0     0
    ## [5,]   6.0 6.5   100

Typically, the `apply_tci` will be used to calculate infusion rates
required to reach specified targets. `apply_tci` requires 1) a set of
target concentrations (or PD response values) and corresponding times at
which the target is set, and 2) a `pkmod` object. It has “plasma” and
“effect” settings, implementing the Jacobs and Shafer algorithms,
respectively. Effect-site targeting will be used by default if the model
has a effect-site equilibrium parameter, *k*<sub>*e*0</sub>, defined.

``` r
# target concentrations for TCI algorithm
plasma_targets <- cbind(value = c(2,3,4,4), time = c(0,2,3,10))

# set ignore_pd = TRUE to indicate that targets are PK, not PD, values
inf_pkvals <- apply_tci(plasma_targets, pkmod = mod_elvd_sample, ignore_pd = TRUE)

# target response using effect-site algorithm
effect_targets <- cbind(value = c(70,60,50,50), time = c(0,2,3,10))
inf_pdvals <- apply_tci(effect_targets, pkmod = mod_elvd)
```

## Predict and simulate methods

`predict.pkmod` and `simulate.pkmod` methods are used to predict and
simulate PK or PD responses at designated times.

``` r
# prediction/observation times
tms_pred <- seq(0,10,0.01)
tms_obs <- seq(0,10,1/3)

# predict values at point estimates
pre <- predict(mod_elvd, inf = inf_pdvals, tms = tms_pred)

# simulate responses under sampled parameter values, representing IIV
obs <- simulate(mod_elvd_sample, seed = 1, inf = inf_pdvals, tms = tms_obs)

# plot results
dat_pre <- data.frame(time = tms_pred, value = pre[,"pdresp"])
dat_obs <- data.frame(time = tms_obs, con = obs)

ggplot(dat_pre, aes(x = time, y = value)) + 
  geom_line() + 
  geom_point(data = dat_obs, aes(x = time, y = con)) +
  xlab("Minutes") + ylab("Bispectral Index")
```

<img src="README_files/figure-gfm/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />
