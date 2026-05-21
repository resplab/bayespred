#' Fit a Bayesian prediction model
#'
#' Fits a logistic regression model under one of four shrinkage priors and
#' returns a self-contained `bpmfit` object. The deployable `bpm` posterior
#' (coefficients, covariance, link) is stored as `$posterior` and can be
#' extracted with [posterior()]. Both `bpmfit` and `bpm` objects accept
#' [predict()]; the latter always requires `newdata`.
#'
#' @param formula A model formula (same syntax as [stats::glm()]).
#' @param data A data frame containing the variables in `formula`.
#' @param family A family object. Only `binomial(link = "logit")` is supported
#'   in v1; any other value raises an informative error.
#' @param prior A prior specification created by [flat()], [jeffreys()],
#'   [log_f()], or [bridge()]. Default: `log_f(m = 2)`.
#' @param model Logical. If `TRUE` (default), store the model frame on the
#'   returned object (used by `predict()` when `newdata` is absent).
#' @param ... Currently unused.
#'
#' @return An object of class `"bpmfit"` with components:
#' \describe{
#'   \item{`posterior`}{A `bpm` object (see [posterior()]) containing the
#'     deployable posterior: `coefficients` (MAP estimate), `vcov`
#'     (Laplace-approximated posterior covariance), `family`, `terms`,
#'     `contrasts`, and `xlevels`. This is everything needed for prediction
#'     with base R only.}
#'   \item{`prior`}{The prior specification object.}
#'   \item{`formula`}{The model formula.}
#'   \item{`call`}{The matched call.}
#'   \item{`model`}{The model frame (if `model = TRUE`).}
#'   \item{`fit_method`}{Internal tag identifying the fitting backend.}
#' }
#'
#' @seealso [predict.bpmfit()], [predict.bpm()], [posterior()],
#'   [flat()], [jeffreys()], [log_f()], [bridge()]
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(200), y = rbinom(200, 1, 0.3))
#' fit <- bpmfit(y ~ x, data = d, prior = log_f(m = 2))
#' predict(fit, data.frame(x = 0.5))
#'
#' # Extract the deployable posterior and predict from it
#' post <- posterior(fit)
#' predict(post, data.frame(x = 0.5))
#'
#' @export
bpmfit <- function(formula, data,
                family = binomial(link = "logit"),
                prior  = log_f(m = 2),
                model  = TRUE,
                ...) {
  cl <- match.call()

  if (!identical(family$family, "binomial") ||
      !identical(family$link,   "logit"))
    stop("bayespred v1 supports only `binomial(link = 'logit')`.")

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

  post <- structure(
    list(
      coefficients = fit_result$coefficients,
      vcov         = fit_result$vcov,
      family       = family,
      terms        = terms_obj,
      contrasts    = contrasts,
      xlevels      = xlevels
    ),
    class = "bpm"
  )

  structure(
    list(
      posterior  = post,
      prior      = prior,
      formula    = formula,
      call       = cl,
      model      = if (isTRUE(model)) mf else NULL,
      fit_method = fit_result$fit_method
    ),
    class = "bpmfit"
  )
}
