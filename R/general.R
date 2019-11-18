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

