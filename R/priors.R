#' Flat (improper) prior
#'
#' Equivalent to unpenalised logistic regression via `glm()`.
#'
#' @return A prior specification object of class `c("prior_flat", "bpm_prior")`.
#' @export
flat <- function() {
  structure(list(type = "flat"), class = c("prior_flat", "bpm_prior"))
}

#' Jeffreys prior
#'
#' Fitted via `brglm2::brglm_fit` with `type = "AS_mean"`, which implements the
#' mean bias-reducing adjusted scores (Firth 1993; Kosmidis & Firth 2021).
#' This is equivalent to a Jeffreys prior on the coefficients.
#'
#' @return A prior specification object of class `c("prior_jeffreys", "bpm_prior")`.
#' @export
jeffreys <- function() {
  structure(list(type = "jeffreys"), class = c("prior_jeffreys", "bpm_prior"))
}

#' log-F(m, m) prior
#'
#' Places an independent log-F(m, m) prior on each slope coefficient via data
#' augmentation (Greenland 2001; Hanley & Shapiro 2022). The intercept is left
#' unpenalised. Larger `m` implies stronger shrinkage toward zero.
#'
#' Approximate 95% prior intervals for the odds ratio:
#' * `m = 1`: (1/648, 648)
#' * `m = 2`: (1/39, 39)
#' * `m = 5`: (1/7.1, 7.1)
#'
#' @param m Positive numeric scalar. Pseudo-observation weight per slope.
#'   Default `2`.
#' @return A prior specification object of class `c("prior_logf", "bpm_prior")`.
#' @export
log_f <- function(m = 2) {
  if (!is.numeric(m) || length(m) != 1L || m <= 0)
    stop("`m` must be a single positive number.")
  structure(list(type = "logf", m = m), class = c("prior_logf", "bpm_prior"))
}

#' Bayesian Ridge prior
#'
#' Estimates a global ridge (L2) penalty on all slope coefficients via REML,
#' using `mgcv::gam` with `paraPen`. Predictors are standardised internally;
#' coefficients are returned on the original (unstandardised) scale.
#'
#' @return A prior specification object of class `c("prior_bridge", "bpm_prior")`.
#' @export
bridge <- function() {
  structure(list(type = "bridge"), class = c("prior_bridge", "bpm_prior"))
}

#' @export
print.bpm_prior <- function(x, ...) {
  msg <- switch(x$type,
    flat     = "flat (unpenalised glm)",
    jeffreys = "Jeffreys (brglm2, AS_mean)",
    logf     = sprintf("log-F(m = %g)", x$m),
    bridge   = "Bayesian Ridge (mgcv REML)",
    paste("unknown:", x$type)
  )
  cat("Prior:", msg, "\n")
  invisible(x)
}
