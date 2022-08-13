README
================

# tci

The `tci` package implements target-controlled infusion (TCI) algorithms
to calculate infusion schedules for compartmental pharmacokinetic (PK)
and pharmacokinetic-pharmacodynamic (PK-PD) models. Using these infusion
schedules, functions in `tci` can be used to simulate PK or PK-PD
responses with or without model misspecification. Generated responses
can be used to simulate closed-loop control by implementing Bayesian
updates to the PK or PK-PD model based on the “observed” (i.e.,
simulated) data.

Closed-form solutions are provided for one, two, or three compartment
mammillary models (i.e., all peripheral compartments are joined to a
central compartment), as well as a three-compartment model with an
adjoining effect-site. Alternative models, potentially based on ordinary
differential equations (ODE), can be specified using other R packages
and adapted for use with `tci` functions (see the “Custom PK models and
algorithms” vignette).

`tci` implements both plasma-and effect-site targeting TCI algorithms
based on work by Jacobs (1990) and Shafer and Gregg (1992),
respectively. Users can implement alternative user-defined TCI
algorithms, however. Again, see the “Custom PK models and algorithms”
vignette for further details.

Several population PK models commonly used for TCI are implemented in
`tci`. These include the Marsh, Schnider, and Eleveld models for
propofol, and the Minto, Kim, and Eleveld models for remifentanil. Use
of these models is illustrated in the “Population PK models” vignette.

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

Equations implementing 1-, 2-, 3-compartment and 3-compartment-effect
structural PK models are included in the `tci` package. The function
`pkmod` will automatically infer the correct structure based on the
parameter names.

``` r
# 3-compartment model with effect site
(mod3ecpt <- pkmod(pars_pk = c(cl = 10, q2 = 2, q3 =20, v = 15, v2 = 30, v3 = 50, ke0 = 1.2)))
```

    ## tci pkmod object
    ## See ?update.pkmod to modify or add elements
    ## 
    ## PK model 
    ##  4-compartment PK model 
    ##  PK parameters: cl = 10, q2 = 2, q3 = 20, v = 15, v2 = 30, v3 = 50, ke0 = 1.2 
    ##  Initial concentrations: (0,0,0,0) 
    ##  Plasma compartment: 1 
    ##  Effect compartment: 4 
    ## 
    ## Simulation
    ##  Additive error SD: 0 
    ##  Multiplicative error SD: 0 
    ##  Logged response: FALSE

Acceptable parameter names can be viewed by calling `list_parnms()`.
Less-commonly used parameters, such as clearance from a peripheral
compartment, are also permissible.

``` r
# acceptable parameter names
list_parnms()
```

    ## Acceptable names for 'pars_pk' vector (case-insensitive) 
    ## 
    ## First compartment options
    ##  Central volume: 'v','v1' 
    ##  Elimination: 'cl','cl1','k10','ke' 
    ## 
    ## Second compartment options
    ##  Peripheral volume: 'v2' 
    ##  Transfer: 'q','q2','k12','k21' 
    ##  Elimination: 'cl2','k20' 
    ## 
    ## Third compartment options
    ##  Second peripheral volume: 'v3' 
    ##  Transfer: 'q3','k13','k31' 
    ##  Elimination: 'cl3','k30' 
    ## 
    ## Effect-site
    ##  Elimination: 'ke0'

Elements of `pkmod` objects can be modified through an `update.pkmod`
method. Perhaps most usefully, this allows for partial modifications to
PK-PD parameters. For example, the effect-site equilibrium constant can
be easily updated.

``` r
update(mod3ecpt, pars_pk = c(ke0 = 0.9), init = c(1,0.2,0.3,1))
```

    ## tci pkmod object
    ## See ?update.pkmod to modify or add elements
    ## 
    ## PK model 
    ##  4-compartment PK model 
    ##  PK parameters: cl = 10, q2 = 2, q3 = 20, v = 15, v2 = 30, v3 = 50, ke0 = 0.9 
    ##  Initial concentrations: (1,0.2,0.3,1) 
    ##  Plasma compartment: 1 
    ##  Effect compartment: 4 
    ## 
    ## Simulation
    ##  Additive error SD: 0 
    ##  Multiplicative error SD: 0 
    ##  Logged response: FALSE

Most functions in the `tci` package pass additional arguments to
`update.pkmod` allowing for easy modification of `pkmod` objects as
needed.

## TCI Infusion schedules

TCI algorithms are implemented using the function `tci_inf` (manual
infusions are implemented by `inf_manual`). The user supplies a set of
targets, times at which the target is set, and a `pkmod` object. The TCI
algorithm (defaults to `type = "plasma"`) is iteratively applied to
calculate infusion rates required to reach each target in turn. By
default, infusion rates are updated in increments of 1/6, corresponding
to every 10-second intervals if infusions rate units are in amount per
minute. Infusion rates themselves must have the same units as the PK
elimination parameters. If elimination rates are in different units,
such as hours, then the TCI update frequency should be modified by the
argument `dtm`.

``` r
# effect-site targeting for three-compartment effect site model
inf_3ecpt <- inf_tci(target_vals = c(2,3,4,4), target_tms = c(0,2,3,10), 
                     pkmod = mod3ecpt, type = "effect")
head(inf_3ecpt)
```

    ##        begin     end inf_rate Ct c1_start   c2_start  c3_start  c4_start
    ## [1,] 0.00000 0.16667 643.6921  2 0.000000 0.00000000 0.0000000 0.0000000
    ## [2,] 0.16667 0.33333   0.0000  2 6.033672 0.03532348 0.2079682 0.5966263
    ## [3,] 0.33333 0.50000   0.0000  2 4.301795 0.09138590 0.5235366 1.4096671
    ## [4,] 0.50000 0.66667   0.0000  2 3.136188 0.13103476 0.7267367 1.8178340
    ## [5,] 0.66667 0.83333  19.3835  2 2.350071 0.15960595 0.8548449 1.9785480
    ## [6,] 0.83333 1.00000   0.0000  2 2.000000 0.18174014 0.9390928 2.0109479
    ##        c1_end     c2_end    c3_end    c4_end
    ## [1,] 6.033672 0.03532348 0.2079682 0.5966263
    ## [2,] 4.301795 0.09138590 0.5235366 1.4096671
    ## [3,] 3.136188 0.13103476 0.7267367 1.8178340
    ## [4,] 2.350071 0.15960595 0.8548449 1.9785480
    ## [5,] 2.000000 0.18174014 0.9390928 2.0109479
    ## [6,] 1.586605 0.19939667 0.9931813 1.9678507

## Population PK models

Population PK models are implemented by the function `poppkmod`. The
user must supply a data frame with the set of covariates (e.g., weight,
age) required by the model. Several published population PK models are
currently implemented in `tci`. For propofol, these include the Marsh,
Schnider, and Eleveld models. For remifentanil, they include the Minto,
Kim, and Eleveld models. See `?poppkmod` or the population PK model
vignette for details. `list_pkmods()` will list available population PK
models and covariates required by each.

``` r
# data frame of patient covariates
data <- data.frame(ID = 1:5, AGE = seq(20,60,by=10), 
                   TBW = seq(60,80,by=5), HGT = seq(150,190,by=10), 
                   MALE = c(TRUE,TRUE,FALSE,FALSE,FALSE))
# Eleveld population PK model for propofol
pkpd_elvd <- poppkmod(data = data, drug = "ppf", model = "eleveld")
```

As with the `pkmod` class, `poppkmod` objects can be used by `inf_tci`
and have `predict` and `simulate` methods to predict and simulate PK-PD
responses, respectively.

PK-PD parameter values can be drawn at random from the
inter-/intra-individual variability distribution, as described by the
`pkmod`
![Omega](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Omega "Omega")
matrix, by either 1) setting the argument `sample = TRUE` when calling
`poppkmod`, or 2) by using the function `sample_iiv`.

``` r
set.seed(1)
pkpd_elvd_iiv <- sample_iiv(pkpd_elvd)
```

## Simulations

Simulations are best implemented through the function `simulate_tci`,
which allows for model misspecification as well as Bayesian updates to
model parameters based on previously observed data (i.e., “closed-loop”
control). `simulate_tci` can be used for both `pkmod` or `poppkmod`
classes. Required arguments to `simulate_tci` are 1) a prior PK model
(`pkmod_prior`) that is used to calculate infusion rates and may be
updated throughout the simulation if update times are provided, 2) a
true PK model (`pkmod_true`) that is used to simulate observations, 3)
TCI target values, 4) TCI target times, and 5) times to simulate
observations.

``` r
# TCI target values (PD response)
target_vals <- c(75,60,50,50)
# values are in terms of minutes. 1/6 = 10 seconds
# TCI target times
target_tms <- c(0,3,6,10)
# observation times 
obs_tms <- seq(1/6,10,1/6)

# simulate without updates ("open-loop")
sim_ol <- simulate_tci(pkmod_prior = pkpd_elvd, pkmod_true = pkpd_elvd_iiv, 
             target_vals, target_tms, obs_tms, type = "effect", seed = 1)
```

`simulate_tci` returns an object with class `sim_tci` that can be
plotted using the `ggplot2` library.

``` r
plot(sim_ol)
```

<img src="README_files/figure-gfm/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

Modifications can be made to the plot to show a subset of responses,
concentrations instead of PD response values, infusion rates, and
simulated data.

``` r
plot(sim_ol, yvar = "c4", id = c(1,3,5), show_inf = TRUE, wrap_id = TRUE)
```

<img src="README_files/figure-gfm/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

Closed-loop simulations can be implemented by specifying a set of update
times. We illustrate this with updates each minute and a processing
delay of 20 seconds.

``` r
sim_cl <- simulate_tci(pkmod_prior = pkpd_elvd, pkmod_true = pkpd_elvd_iiv, 
             target_vals, target_tms, obs_tms, update_tms = 1:10, delay = 1/3,
               type = "effect", seed = 1)
```

    ## [1] "Simulating ID=1"
    ## [1] "Simulating ID=2"
    ## [1] "Simulating ID=3"
    ## [1] "Simulating ID=4"
    ## [1] "Simulating ID=5"

Since `plot.sim_tci` returns a `ggplot2` object, it is easy to modify
aspects such as titles and axis labels using `ggplot2` functions.

``` r
plot(sim_cl) + 
  xlab("Minutes") + 
  ylab("Bispectral Index") + 
  ggtitle("Closed-loop simulation of Eleveld propofol model", 
          subtitle = "Minute updates, processing delay of 20 seconds")
```

<img src="README_files/figure-gfm/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Jacobs1990" class="csl-entry">

Jacobs, James R. 1990. “<span class="nocase">Algorithm for Optimal
Linear Model-Based Control with Application to Pharmacokinetic
Model-Driven Drug Delivery</span>.” *IEEE Transactions on Biomedical
Engineering* 37 (1): 107–9. <https://doi.org/10.1109/10.43622>.

</div>

<div id="ref-Shafer1992" class="csl-entry">

Shafer, Steven L., and Keith M. Gregg. 1992. “<span
class="nocase">Algorithms to rapidly achieve and maintain stable drug
concentrations at the site of drug effect with a computer-controlled
infusion pump</span>.” *Journal of Pharmacokinetics and
Biopharmaceutics* 20 (2): 147–69. <https://doi.org/10.1007/BF01070999>.

</div>

</div>
