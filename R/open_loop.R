

#' Simulate open-loop control
#'
#' Simulate open-loop control with target-controlled infusion for a `pkmod` object
#'
#' @param pkmod_prior `pkmod` object describing a PK/PK-PD model that is used to calculate
#' TCI infusion rates and is updated as data are simulated and incorporated. Must have an
#' associated Omega matrix.
#' @param pkmod_true `pkmod` object describing the patient's "true" response. This model
#' will be used to simulate observations.
#' @param target_vals A vector of numeric values indicating PK or PD targets for TCI algorithm.
#' @param target_tms A vector of numeric values indicating times at which the TCI algorithm should
#' begin targeting each value.
#' @param obs_tms Times at which data values should be simulated from `pkmod_true`.
#' @param type Type of TCI algorithm to be used. Options are "plasma" and "effect".
#' Defaults to "effect". Will be overwritten if `custom_alg` is non-null.
#' @param custom_alg Custom TCI algorithm to overwrite default plasma- or effect-site targeting.
#' @param resp_bounds Optional vector of two values indicating minimum and maximum values possible for the response.
#' @param seed An integer used to initialize the random number generator.
#' @examples
#' pkmod_prior <- pkmod(pars_pk = c(cl = 10, q2 = 2, q3 =20, v = 15, v2 = 30, v3 = 50, ke0 = 1.2))
#' pkmod_true  <- pkmod(pars_pk = c(cl = 16, q2 = 4, q3 =10, v = 20, v2 = 20, v3 = 80, ke0 = 0.8),
#' sigma_add = 0.1, log_response = TRUE)
#' target_vals <- c(2,3,4,3,3)
#' target_tms <- c(0,5,10,36,60)
#' obs_tms <- c(1,2,4,8,12,16,24,36,48)
#' sim <- olc(pkmod_prior, pkmod_true, target_vals, target_tms, obs_tms)
#' len <- 500
#' tms <- seq(0,60,length.out = len)
#' df <- data.frame(time = rep(tms,2),
#'                  value = c(predict(pkmod_true, sim$inf,tms)[,1],
#'                  predict(pkmod_prior, sim$inf,tms)[,1]),
#'                  type = c(rep("true",len),rep("prior",len)))
#' library(ggplot2)
#' ggplot(df, aes(x = time, y = value, color = type)) +
#'   geom_step(data = data.frame(time = target_tms, value = target_vals),
#'   aes(x = time, y = value), inherit.aes = FALSE) +
#'   geom_line() +
#'   geom_point(data = sim$obs, aes(x = time, y = obs), inherit.aes = FALSE)
#' @export
olc <- function(pkmod_prior, pkmod_true, target_vals, target_tms, obs_tms, type=c("effect","plasma"),
                custom_alg = NULL, resp_bounds = NULL, seed = NULL){

  if(class(pkmod_prior) != "pkmod" | class(pkmod_true) != "pkmod")
    stop("Both 'pkmod_prior' and 'pkmod_true' must be pkmod objects.")

  if(is.null(seed)) seed <- sample(1:1e5,1)
  set.seed(seed)
  type = match.arg(type)

  # calculate infusion based on prior model
  inf <- inf_tci(pkmod_prior, target_vals, target_tms, type = type, custom_alg = custom_alg)

  # simulate response from true model
  obs <- simulate(pkmod_true, tms = obs_tms, inf = inf, resp_bounds=resp_bounds)

  return(list(obs = data.frame(time = obs_tms, obs = obs), inf = inf, seed=seed))
}


#' Simulate open-loop control using TCI
#'
#' Simulate open-loop control using TCI for `pkmod` or `poppkmod` objects.
#'
#' @param pkmod_prior `pkmod` object describing a PK/PK-PD model that is used to calculate
#' TCI infusion rates and is updated as data are simulated and incorporated. Must have an
#' associated Omega matrix.
#' @param pkmod_true `pkmod` object describing the patient's "true" response. This model
#' will be used to simulate observations.
#' @param target_vals A vector of numeric values indicating PK or PD targets for TCI algorithm.
#' @param target_tms A vector of numeric values indicating times at which the TCI algorithm should
#' begin targeting each value.
#' @param obs_tms Times at which data values should be simulated from `pkmod_true`.
#' @param type Type of TCI algorithm to be used. Options are "plasma" and "effect".
#' Defaults to "effect". Will be overwritten if `custom_alg` is non-null.
#' @param custom_alg Custom TCI algorithm to overwrite default plasma- or effect-site targeting.
#' @param resp_bounds Optional vector of two values indicating minimum and maximum values possible for the response.
#' @param seed An integer used to initialize the random number generator.
#' @examples
#' data <- data.frame(ID = 1:5, AGE = seq(20,60,by=10), TBW = seq(60,80,by=5),
#' HGT = seq(150,190,by=10), MALE = c(TRUE,TRUE,FALSE,FALSE,FALSE))
#' pkmod_prior <- poppkmod(data, drug = "ppf", model = "eleveld")
#' pkmod_true  <- poppkmod(data, drug = "ppf", model = "eleveld", sample = TRUE)
#' obs_tms <- seq(1/6,10,1/6)
#' target_vals = c(75,60,50,50)
#' target_tms = c(0,3,6,10)
#' sim <- simulate_olc(pkmod_prior, pkmod_true, target_vals, target_tms, obs_tms)
#' len <- 500
#' tms <- seq(0,10,length.out = len)
#' resp <- data.frame(rbind(predict(pkmod_true, sim$inf, tms),
#' predict(pkmod_prior, sim$inf, tms)))
#' resp$type = c(rep("true",len*5),rep("prior",len*5))
#' library(ggplot2)
#' ggplot(resp) + geom_line(aes(x = time, y = pdresp, color = factor(id))) + facet_wrap(~type) +
#'   labs(x = "Hours", y = "Bispectral Index") + theme_bw() +
#'   geom_step(data = data.frame(time = target_tms, value = target_vals),
#'   aes(x = time, y = value), inherit.aes = FALSE)
#' @export
simulate_olc <- function(pkmod_prior, pkmod_true, target_vals, target_tms, obs_tms, type=c("effect","plasma"),
                         custom_alg = NULL, resp_bounds = NULL, seed = NULL){

  if(is.null(seed)) seed <- sample(1:1e5,1)
  set.seed(seed)
  type <- match.arg(type)

  if(class(pkmod_prior) != class(pkmod_true) | !class(pkmod_prior) %in% c("pkmod","poppkmod"))
    stop("pkmod_prior and pkmod_true must either both have class 'pkmod' or 'poppkmod'")

  if(class(pkmod_prior) == "pkmod"){
    out <- olc(pkmod_prior, pkmod_true, target_vals, target_tms, obs_tms, type, custom_alg, resp_bounds)
    ids <- 1
    out$obs$type = "true"
  } else{
    if(length(pkmod_prior$pkmods) != length(pkmod_true$pkmods))
      stop("'pkmod_prior' and 'pkmod_true' must have the same number of pkmods")

    sim <- lapply(1:length(pkmod_prior$pkmods), function(i){
      olc(pkmod_prior$pkmods[[i]], pkmod_true$pkmods[[i]], target_vals, target_tms, obs_tms, type, custom_alg, resp_bounds)
    })

    ids <- as.factor(pkmod_prior$ids)
    out <- list(obs = cbind(id = rep(ids, each = nrow(sim[[1]]$obs)),
                            type = "true",
                     do.call("rbind",lapply(sim, `[[`, "obs"))),
         inf = cbind(id = rep(ids, each = nrow(sim[[1]]$inf)),
                     do.call("rbind",lapply(sim, `[[`, "inf"))),
         seed = sapply(sim, `[[`, "seed"))
  }

  tmseq <- seq(min(out$inf[,"begin"]), max(out$inf[,"end"]), length.out = 1e3)
  out$resp <- as.data.frame(rbind(
    predict(pkmod_prior, out$inf,tms = tmseq),
    predict(pkmod_true, out$inf, tms = tmseq)))
  out$resp$type <- c(rep("prior", 1e3*length(ids)), rep("true", 1e3*length(ids)))
  if("id" %in% names(out$resp)) out$resp$id <- as.factor(out$resp$id)
  if(!"time" %in% names(out$resp)) out$resp <- cbind(time = tmseq, out$resp)

  out
}
