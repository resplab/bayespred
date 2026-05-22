# Internal: shared prediction core.
# `posterior` is a bpm object; `X_new` is the already-built design matrix.
.predict_bpm_core <- function(posterior, X_new, type, method, se.fit, interval) {
  beta    <- posterior$coefficients
  V       <- posterior$vcov
  eta     <- as.numeric(X_new %*% beta)
  var_eta <- as.numeric(rowSums((X_new %*% V) * X_new))
  sigma   <- sqrt(pmax(var_eta, 0))

  # SE on the same scale as `type`, matching predict.glm() convention.
  # link:     SE of the linear predictor = sigma = sqrt(x'Vx)
  # response: delta-method SE of the probability = mu(1-mu) * sigma
  se_fit <- if (type == "link") {
    sigma
  } else {
    mu_pe <- plogis(eta)
    mu_pe * (1 - mu_pe) * sigma
  }

  # ---- type = "link" ----------------------------------------------------------
  if (type == "link") {
    if (!is.null(interval)) {
      if (!is.numeric(interval) || length(interval) != 1L ||
          interval < 0 || interval >= 1)
        stop("`interval` must be 0 (SE only) or a single number in (0, 1).",
             call. = FALSE)
      if (interval == 0)
        return(data.frame(fit = eta, se.fit = se_fit,
                          row.names = rownames(X_new)))
      z <- qnorm(1 - (1 - interval) / 2)
      return(data.frame(fit    = eta,
                        lwr    = eta - z * sigma,
                        upr    = eta + z * sigma,
                        se.fit = se_fit,
                        row.names = rownames(X_new)))
    }
    if (isTRUE(se.fit))
      return(list(fit = eta, se.fit = se_fit, residual.scale = 1))
    return(eta)
  }

  # ---- type = "response" ------------------------------------------------------
  fit <- switch(method,
    pe        = plogis(eta),
    pm        = .pm_quadrature(eta, sigma),
    pm_mackay = .pm_mackay(eta, sigma)
  )

  if (!is.null(interval)) {
    if (!is.numeric(interval) || length(interval) != 1L ||
        interval < 0 || interval >= 1)
      stop("`interval` must be 0 (SE only) or a single number in (0, 1).",
           call. = FALSE)
    if (interval == 0)
      return(data.frame(fit = fit, se.fit = se_fit,
                        row.names = rownames(X_new)))
    z <- qnorm(1 - (1 - interval) / 2)
    return(data.frame(fit    = fit,
                      lwr    = plogis(eta - z * sigma),
                      upr    = plogis(eta + z * sigma),
                      se.fit = se_fit,
                      row.names = rownames(X_new)))
  }

  if (isTRUE(se.fit))
    return(list(fit = fit, se.fit = se_fit, residual.scale = 1))

  fit
}

# ------------------------------------------------------------------------------

#' Predict from a bpm posterior object
#'
#' Computes predictions from a `bpm` object (the deployable posterior obtained
#' via [posterior()]). Requires `newdata`; no training data is stored on a
#' `bpm` object. For in-sample predictions, call [predict.bpmfit()] on the
#' original `bpmfit` instead.
#'
#' The interface is intentionally compatible with [stats::predict.glm()].
#' See [predict.bpmfit()] for full parameter documentation.
#'
#' @param object A `bpm` object (obtained via [posterior()]).
#' @param newdata A data frame of new observations. Always required.
#' @param type Character: `"response"` (default) or `"link"`.
#' @param method Character: `"pm"` (default), `"pe"`, or `"pm_mackay"`.
#'   Ignored when `type = "link"`.
#' @param se.fit Logical (default `FALSE`). If `TRUE`, return a list with
#'   `fit`, `se.fit`, and `residual.scale`, matching [stats::predict.glm()].
#' @param interval `NULL` (default), `0` (SE only, no bounds), or a coverage
#'   probability in (0, 1). Cannot be combined with `se.fit = TRUE`.
#' @param dispersion Ignored; included for [stats::predict.glm()] compatibility.
#' @param na.action Function to handle `NA`s in `newdata`. Default
#'   [stats::na.pass].
#' @param ... Currently unused.
#'
#' @return See [predict.bpmfit()] for return-value details.
#'
#' @seealso [posterior()], [predict.bpmfit()]
#' @export
predict.bpm <- function(object, newdata,
                        type       = c("response", "link", "terms"),
                        method     = c("pm", "pe", "pm_mackay"),
                        se.fit     = FALSE,
                        interval   = NULL,
                        dispersion = NULL,
                        na.action  = na.pass,
                        ...) {
  if (missing(newdata) || is.null(newdata))
    stop(
      "`newdata` is required for `predict.bpm()`. ",
      "To predict on the training data, call `predict()` on the `bpmfit` ",
      "object instead.",
      call. = FALSE
    )

  type   <- match.arg(type)
  method <- match.arg(method)

  if (type == "terms")
    stop("`type = 'terms'` is not implemented for bpm objects.", call. = FALSE)

  if (!is.null(interval) && isTRUE(se.fit))
    stop("`interval` and `se.fit = TRUE` cannot be used together.", call. = FALSE)

  tt    <- delete.response(object$terms)
  mf_nd <- model.frame(tt, newdata, xlev = object$xlevels, na.action = na.action)
  X_new <- model.matrix(tt, mf_nd, contrasts.arg = object$contrasts)

  .predict_bpm_core(object, X_new, type, method, se.fit, interval)
}

# ------------------------------------------------------------------------------

#' Predict from a bpmfit model
#'
#' Computes predictions from a fitted `bpmfit` object. The interface is
#' intentionally compatible with [stats::predict.glm()]: `type` controls the
#' output scale and `se.fit` requests standard errors. The additional `method`
#' argument selects the Bayesian computation used for `type = "response"`.
#'
#' When `newdata` is omitted, predictions are computed on the training data
#' using the stored model frame (requires `model = TRUE` at fit time). To
#' predict on new data from a portable, self-contained object, extract the
#' posterior with [posterior()] and call [predict.bpm()] directly.
#'
#' @param object A `bpmfit` object.
#' @param newdata A data frame of new observations. If omitted, predictions are
#'   made on the training data (requires `model = TRUE` at fit time).
#' @param type Character string controlling the output scale (same as
#'   [stats::predict.glm()]):
#' \describe{
#'   \item{`"response"`}{(default) Predicted probability, controlled by `method`.}
#'   \item{`"link"`}{Linear predictor \eqn{\eta = x^\top\hat\beta}.
#'     `method` is ignored.}
#'   \item{`"terms"`}{Not implemented in v1.}
#' }
#' @param method Character string selecting the Bayesian computation. Only
#'   relevant when `type = "response"`:
#' \describe{
#'   \item{`"pm"`}{(default) Posterior mean via 30-point Gauss-Hermite quadrature.}
#'   \item{`"pe"`}{Plug-in estimate: `plogis(X %*% beta_hat)`. Equivalent to
#'     `predict.glm(..., type = "response")`.}
#'   \item{`"pm_mackay"`}{MacKay (1992) closed-form approximation to the PM.}
#' }
#' @param se.fit Logical. If `TRUE`, return a list with elements `fit`,
#'   `se.fit`, and `residual.scale`, matching [stats::predict.glm()].
#'   `se.fit` is on the same scale as `type`: for `"link"` it is
#'   \eqn{\sqrt{x^\top V x}}; for `"response"` it is the delta-method
#'   approximation \eqn{\hat\mu(1-\hat\mu)\sqrt{x^\top V x}}. Default `FALSE`.
#' @param interval Either `NULL` (point predictions only, the default), `0`
#'   (return `fit` and `se.fit` without bounds), or a numeric scalar in (0, 1)
#'   giving credible interval coverage (e.g. `0.95`). Intervals are computed on
#'   the linear-predictor scale (normal approximation) and back-transformed via
#'   `plogis` for `type = "response"`. Cannot be combined with `se.fit = TRUE`.
#' @param dispersion Ignored. Included for compatibility with
#'   [stats::predict.glm()]; the binomial dispersion is always 1.
#' @param na.action Function to handle `NA`s in `newdata`. Default
#'   [stats::na.pass], matching [stats::predict.glm()].
#' @param ... Currently unused.
#'
#' @return
#' \itemize{
#'   \item If `interval = NULL` and `se.fit = FALSE`: a numeric vector.
#'   \item If `interval = 0`: a data frame with columns `fit` and `se.fit`
#'     only — no bounds. `se.fit` is on the same scale as `type`.
#'   \item If `interval` is in (0, 1): a data frame with columns `fit`,
#'     `lwr`, `upr`, and `se.fit` (on the same scale as `type`).
#'   \item If `se.fit = TRUE`: a list with elements `fit`, `se.fit`, and
#'     `residual.scale` (1 for binomial).
#' }
#'
#' @details
#' The marginal posterior of the linear predictor is approximated as
#' \deqn{\eta_i \mid \text{data} \approx N(x_i^\top\hat\beta,\; x_i^\top V x_i).}
#' All credible intervals and standard errors are derived from this
#' approximation. For `type = "response"` the PM point estimate lies inside
#' but not at the centre of the interval (Jensen's inequality).
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(200), y = rbinom(200, 1, 0.3))
#' fit <- bpmfit(y ~ x, data = d, prior = log_f(m = 2))
#'
#' # glm-style
#' predict(fit, data.frame(x = 0.5), type = "response")
#' predict(fit, data.frame(x = 0.5), type = "link")
#' predict(fit, data.frame(x = 0.5), type = "link",    se.fit = TRUE)
#' predict(fit, data.frame(x = 0.5), type = "response", se.fit = TRUE)
#'
#' # Bayesian extensions
#' predict(fit, data.frame(x = 0.5), method = "pm",   interval = 0.95)
#' predict(fit, data.frame(x = 0.5), type = "link",   interval = 0.95)
#'
#' @export
predict.bpmfit <- function(object, newdata,
                        type       = c("response", "link", "terms"),
                        method     = c("pm", "pe", "pm_mackay"),
                        se.fit     = FALSE,
                        interval   = NULL,
                        dispersion = NULL,
                        na.action  = na.pass,
                        ...) {
  if (identical(method, "pm_proj"))
    stop(
      '`method = "pm_proj"` is no longer available on bpmfit objects. ',
      "Use `predict(bpmproj_pm(fit), newdata)` instead.",
      call. = FALSE
    )

  type   <- match.arg(type)
  method <- match.arg(method)

  if (type == "terms")
    stop("`type = 'terms'` is not implemented for bpmfit objects.", call. = FALSE)

  if (!is.null(interval) && isTRUE(se.fit))
    stop("`interval` and `se.fit = TRUE` cannot be used together.", call. = FALSE)

  # ---- build prediction design matrix ----------------------------------------
  post <- object$posterior
  if (missing(newdata) || is.null(newdata)) {
    if (is.null(object$model))
      stop("No model frame stored. Refit with `model = TRUE` or supply `newdata`.",
           call. = FALSE)
    X_new <- model.matrix(post$terms, object$model,
                          contrasts.arg = post$contrasts)
  } else {
    tt    <- delete.response(post$terms)
    mf_nd <- model.frame(tt, newdata, xlev = post$xlevels,
                         na.action = na.action)
    X_new <- model.matrix(tt, mf_nd, contrasts.arg = post$contrasts)
  }

  .predict_bpm_core(post, X_new, type, method, se.fit, interval)
}
