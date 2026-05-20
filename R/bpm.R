#' Fit a Bayesian prediction model
#'
#' Fits a logistic regression model under one of four shrinkage priors and
#' returns a self-contained `bpm` object that can be used with
#' [predict.bpm()] using base R only at deployment time.
#'
#' @param formula A model formula (same syntax as [stats::glm()]).
#' @param data A data frame containing the variables in `formula`.
#' @param family A family object. Only `binomial(link = "logit")` is supported
#'   in v1; any other value raises an informative error.
#' @param prior A prior specification created by [flat()], [jeffreys()],
#'   [log_f()], or [bridge()]. Default: `log_f(m = 2)`.
#' @param projpred Logical. If `TRUE`, compute self-projection coefficients at
#'   fit time and cache them on the returned object. Default `FALSE`. Required
#'   for `predict(..., method = "pm_proj")`.
#' @param model Logical. If `TRUE` (default), store the model frame on the
#'   returned object (used by `predict()` when `newdata` is absent).
#' @param ... Currently unused.
#'
#' @return An object of class `"bpm"` with components:
#' \describe{
#'   \item{`coefficients`}{Named numeric vector of fitted coefficients.}
#'   \item{`vcov`}{Posterior covariance matrix (Laplace approximation).}
#'   \item{`family`}{The family object.}
#'   \item{`prior`}{The prior specification object.}
#'   \item{`formula`}{The model formula.}
#'   \item{`terms`}{Terms object for building prediction design matrices.}
#'   \item{`contrasts`}{Factor contrast encodings from the training data.}
#'   \item{`xlevels`}{Factor level sets from the training data.}
#'   \item{`call`}{The matched call.}
#'   \item{`model`}{The model frame (if `model = TRUE`).}
#'   \item{`projection`}{Self-projection coefficients (if `projpred = TRUE`).}
#'   \item{`fit_method`}{Internal tag identifying the fitting backend.}
#' }
#'
#' @seealso [predict.bpm()], [flat()], [jeffreys()], [log_f()], [bridge()]
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(200), y = rbinom(200, 1, 0.3))
#' fit <- bpm(y ~ x, data = d, prior = log_f(m = 2))
#' predict(fit, data.frame(x = 0.5))
#'
#' @export
bpm <- function(formula, data,
                family   = binomial(link = "logit"),
                prior    = log_f(m = 2),
                projpred = FALSE,
                model    = TRUE,
                ...) {
  cl <- match.call()

  if (!identical(family$family, "binomial") ||
      !identical(family$link,   "logit"))
    stop("BayesCPM v1 supports only `binomial(link = 'logit')`.")

  if (!inherits(prior, "bpm_prior"))
    stop("`prior` must be created by flat(), jeffreys(), log_f(), or bridge().")

  mf        <- model.frame(formula, data = data, na.action = na.omit)
  y         <- model.response(mf)
  terms_obj <- attr(mf, "terms")
  X         <- model.matrix(terms_obj, mf)
  contrasts <- attr(X, "contrasts")
  xlevels   <- .get_xlevels(terms_obj, mf)

  fit_result <- switch(class(prior)[1L],
    prior_flat     = .fit_flat(X, y, family),
    prior_jeffreys = .fit_jeffreys(X, y, family),
    prior_logf     = .fit_logf(X, y, prior$m, family),
    prior_bridge   = .fit_bridge(X, y, family),
    stop("Unknown prior class: ", class(prior)[1L])
  )

  if (!fit_result$converged)
    warning("Model did not converge.", call. = FALSE)

  projection <- if (isTRUE(projpred))
    c(
      .compute_projection(X, X, fit_result$coefficients, fit_result$vcov,
                          family, self = TRUE),
      list(
        terms     = terms_obj,
        contrasts = contrasts,
        xlevels   = xlevels,
        family    = family
      )
    )
  else
    NULL

  structure(
    list(
      coefficients = fit_result$coefficients,
      vcov         = fit_result$vcov,
      family       = family,
      prior        = prior,
      formula      = formula,
      terms        = terms_obj,
      contrasts    = contrasts,
      xlevels      = xlevels,
      call         = cl,
      model      = if (isTRUE(model)) mf else NULL,
      projection = projection,
      fit_method = fit_result$fit_method
    ),
    class = "bpm"
  )
}
