#' Predict from a bpm model
#'
#' Computes predictions from a fitted `bpm` object. The interface is
#' intentionally compatible with [stats::predict.glm()]: `type` controls the
#' output scale and `se.fit` requests standard errors. The additional `method`
#' argument selects the Bayesian computation used for `type = "response"`.
#'
#' @param object A `bpm` object.
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
#'   `se.fit`, and `residual.scale`, matching the structure of
#'   [stats::predict.glm()]. For `type = "link"`, `se.fit` is
#'   \eqn{\sqrt{x^\top V x}}. For `type = "response"`, `se.fit` is the
#'   delta-method approximation \eqn{\hat\mu(1-\hat\mu)\sqrt{x^\top V x}}
#'   evaluated at the plug-in mean. Default `FALSE`.
#' @param interval Either `NULL` (point predictions only, the default) or a
#'   numeric scalar in (0, 1) giving credible interval coverage (e.g. `0.95`).
#'   The interval is computed on the linear-predictor scale (normal
#'   approximation) and back-transformed via `plogis` for `type = "response"`.
#'   Cannot be combined with `se.fit = TRUE`.
#' @param dispersion Ignored. Included for compatibility with
#'   [stats::predict.glm()]; the binomial dispersion is always 1.
#' @param na.action Function to handle `NA`s in `newdata`. Default
#'   [stats::na.pass], matching [stats::predict.glm()].
#' @param ... Currently unused.
#'
#' @return
#' \itemize{
#'   \item If `interval = NULL` and `se.fit = FALSE`: a numeric vector.
#'   \item If `interval` is specified: a data frame with columns `fit`,
#'     `lwr`, `upr`, and `se.link` (SE of the linear predictor,
#'     \eqn{\sqrt{x^\top V x}}). `se.link` is included automatically because
#'     it is computed regardless and is the fundamental uncertainty quantity.
#'   \item If `se.fit = TRUE`: a list with elements `fit`, `se.fit`, and
#'     `residual.scale` (1 for binomial), matching [stats::predict.glm()].
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
#' fit <- bpm(y ~ x, data = d, prior = log_f(m = 2))
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
predict.bpm <- function(object, newdata,
                        type       = c("response", "link", "terms"),
                        method     = c("pm", "pe", "pm_mackay"),
                        se.fit     = FALSE,
                        interval   = NULL,
                        dispersion = NULL,
                        na.action  = na.pass,
                        ...) {
  if (identical(method, "pm_proj"))
    stop(
      '`method = "pm_proj"` is no longer available on bpm objects. ',
      "Use `predict(project_pm(fit), newdata)` instead.",
      call. = FALSE
    )

  type   <- match.arg(type)
  method <- match.arg(method)

  if (type == "terms")
    stop("`type = 'terms'` is not implemented for bpm objects.", call. = FALSE)

  if (!is.null(interval) && isTRUE(se.fit))
    stop("`interval` and `se.fit = TRUE` cannot be used together.", call. = FALSE)

  # ---- build prediction design matrix ----------------------------------------
  if (missing(newdata) || is.null(newdata)) {
    if (is.null(object$model))
      stop("No model frame stored. Refit with `model = TRUE` or supply `newdata`.",
           call. = FALSE)
    X_new <- model.matrix(object$terms, object$model,
                          contrasts.arg = object$contrasts)
  } else {
    tt    <- delete.response(object$terms)
    mf_nd <- model.frame(tt, newdata, xlev = object$xlevels,
                         na.action = na.action)
    X_new <- model.matrix(tt, mf_nd, contrasts.arg = object$contrasts)
  }

  beta <- object$coefficients
  V    <- object$vcov

  eta     <- as.numeric(X_new %*% beta)
  var_eta <- as.numeric(rowSums((X_new %*% V) * X_new))
  sigma   <- sqrt(pmax(var_eta, 0))

  # ---- type = "link" ----------------------------------------------------------
  if (type == "link") {
    if (!is.null(interval)) {
      if (!is.numeric(interval) || length(interval) != 1L ||
          interval <= 0 || interval >= 1)
        stop("`interval` must be a single number in (0, 1).", call. = FALSE)
      z <- qnorm(1 - (1 - interval) / 2)
      return(data.frame(fit     = eta,
                        lwr     = eta - z * sigma,
                        upr     = eta + z * sigma,
                        se.link = sigma,
                        row.names = rownames(X_new)))
    }
    if (isTRUE(se.fit))
      return(list(fit = eta, se.fit = sigma, residual.scale = 1))
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
        interval <= 0 || interval >= 1)
      stop("`interval` must be a single number in (0, 1).", call. = FALSE)
    z <- qnorm(1 - (1 - interval) / 2)
    return(data.frame(fit     = fit,
                      lwr     = plogis(eta - z * sigma),
                      upr     = plogis(eta + z * sigma),
                      se.link = sigma,
                      row.names = rownames(X_new)))
  }

  if (isTRUE(se.fit)) {
    # Delta-method SE for the response scale, evaluated at the plug-in mean.
    # Matches predict.glm(type = "response", se.fit = TRUE) for method = "pe".
    mu_pe  <- plogis(eta)
    se_fit <- mu_pe * (1 - mu_pe) * sigma
    return(list(fit = fit, se.fit = se_fit, residual.scale = 1))
  }

  fit
}
