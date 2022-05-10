README
================

# tci

The `tci` package implements plasma- and effect-site targeting
target-controlled infusion (TCI) algorithms to 1-, 2-, 3-, and
3-compartment/effect-site pharmacokinetic (PK) models with intravenous
drug administration. TCI algorithms can further be applied to
pharmacodynamic (PD) targets when a PD model and its inverse are
specified. Functions are supplied for simulation of individual responses
to TCI infusion schedules under open- or Bayesian closed-loop control.
See the `overview` vignette for further details. Custom user-defined PK
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

PK and PKPD models are created using the function `pkmod` for individual
PK/PKPD models or `poppkmod` to access previously published population
PK models. The only required argument to `pkmod` is a named vector
`pars_pk`, indicating the PK parameters for a specific individual.
Acceptable PK parameter names can be found by calling `list_parnms()`.
The number of compartments will be inferred from the parameter names.

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
models. Each population PK model function takes individual covariates
(e.g., age, total body weight) and returns a `pkmod` object that
implements the population model at the corresponding structural PK/PK-PD
parameter values.

Population PK/PKPD models are created using the function `poppkmod`. The
primary argument to `poppkmod` is a data frame with individual covariate
values to be evaluated by the model. Available models and required
covariates can be viewed by printing `list_pkmods()`.

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

The user additionally selects the appropriate drug (“ppf” for propofol
or “remi” for remifentanil) and model. By setting the argument
`sample=TRUE` PK/PKPD parameters will be randomly sampled from the
distribution described by the model’s interindividual variability.

``` r
data <- data.frame(ID = 1:5, 
                   AGE = seq(20,60,by=10), 
                   TBW = seq(60,80,by=5), 
                   HGT = seq(150,190,by=10), 
                   MALE = c(TRUE,TRUE,FALSE,FALSE,FALSE))
ppf_eleveld <- poppkmod(data, drug = "ppf", model = "eleveld")
```

## Predict and simulate methods

S3 methods for `predict` and `simulate` exist for both `pkmod` and
`poppkmod` objects. `predict` is used to predict concentrations over
time associated with a model object, while `simulate` is used to
simulate observations from the model. Both methods require an infusion
schedule that indicates infusion begin times, end times, and infusion
rates with columns labeled “begin,” “end,” and “inf\_rate.” Infusion
schedules can be created manually via the function `inf_manual` or
calculated by a TCI algorithm to reach designated targets (plasma or
effect-site) with `inf_tci`. By default, `inf_tci` will interpret
targets as PD response values if the model object has a PD component.
This can be overruled by setting the argument `ignore_pd=TRUE`,
indicating that targets are concentrations.

``` r
# Infusions of 100 mg/min for 30 sec at 0, 3, 6 min
(multi_inf <- inf_manual(inf_tms = c(0,3,6), inf_rate = 100, duration = 0.5))
```

    ##      begin end inf_rate
    ## [1,]   0.0 0.5      100
    ## [2,]   0.5 3.0        0
    ## [3,]   3.0 3.5      100
    ## [4,]   3.5 6.0        0
    ## [5,]   6.0 6.5      100

``` r
# Infusions to reach and maintain a target effect-site concentration of 3 mg/L for
# ten minutes, updated each 30 seconds
inf_pk <- inf_tci(pkmod = ppf_eleveld, 
                   target_vals = c(3,3), 
                   target_tms = c(0,10), 
                   type = "effect", 
                   dtm = 0.5, 
                   ignore_pd = TRUE)
# plot results
ggplot(as.data.frame(inf_pk)) + 
  geom_step(aes(x = begin, y = inf_rate, color = factor(id))) +
  xlab("Minutes") + ylab("Infusion rate (mg/min)")
```

<img src="README_files/figure-gfm/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

## Open- and closed-loop simulations

A system is termed “open-loop” if a controlled variable, such as the
infusion rate, is calculated without feedback from the system and
closed-loop if feedback is incorporated. By default, TCI systems are
open-loop: infusion rates are calculated to reach targets according to a
PK or PKPD model that is believed to describe an individual. In
practice, a clinician will adjust infusion rates as needed, thereby
manually “closing the loop.” The `tci` package simulates this process,
absent clinician involvement, using the function `simulate_olc`.

In practice, the dynamics governing an individual’s response will not be
identical to those predicted by a model. To simulate this mismatch,
therefore, two models are required: 1) a “prior” model that is believed
to describe the individual and used to calculate infusion rates, and 2)
a “true” model that is used to generate responses given the infusions
administered. Typically, the prior will be a population model evaluated
at patient covariate values, while the true model parameters may be
sampled from the distributions representing inter- and/or intra-patient
variability in parameters, or may be associated with a completely
different model to simulate model misspecification. In addition to these
two models, `simulate_olc` requires a set of targets, times at which
targets are set, and times at which observations are to be simulated.

``` r
# Sample "true" PKPD parameters from distribution of interindividual variability
set.seed(1)
ppf_eleveld_true <- poppkmod(data, drug = "ppf", model = "eleveld", sample = TRUE)

# Target BIS=50 for 10 minutes with observations collected every 30 seconds
sim_ol <- simulate_olc(pkmod_prior = ppf_eleveld, 
                       pkmod_true = ppf_eleveld_true, 
                       target_vals = c(50,50),
                       target_tms = c(0,10), 
                       obs_tms = seq(1/2,10,1/2))

ggplot(sim_ol$resp) + 
  geom_line(aes(x = time, y = pdresp, color = id)) + 
  geom_point(data = sim_ol$obs, aes(x = time, y = obs, color = id), alpha = 0.3) +
  facet_wrap(~type) +
  labs(x = "Hours", y = "Bispectral Index")
```

<img src="README_files/figure-gfm/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

During closed-loop control, as implemented in `tci`, simulated data are
used to periodically perform Bayesian updates of the prior patient model
to estimate a posterior patient-specific model. This model is then used
to calculate subsequent infusion rates. Closed-loop control is simulated
by the function `simulate_clc`, which takes the same arguments as
`simulate_olc` in addition to a set of times at which updates should be
performed.

``` r
# Target BIS=50 for 10 minutes with observations collected every 30 seconds
sim_cl <- simulate_clc(pkmod_prior = ppf_eleveld, 
                       pkmod_true = ppf_eleveld_true, 
                       target_vals = c(50,50),
                       target_tms = c(0,10), 
                       obs_tms = seq(1/2,10,1/2), 
                       update_tms = 1:10, 
                       verbose = FALSE)

ggplot(sim_cl$resp) + 
  geom_line(aes(x = time, y = pdresp, color = id)) + 
  geom_point(data = sim_cl$obs, aes(x = time, y = obs, color = id), alpha = 0.3) +
  facet_wrap(~type) +
  labs(x = "Hours", y = "Bispectral Index")
```

<img src="README_files/figure-gfm/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />
