# Compute self-projection coefficients (Sadatsafavi & Riley 2026, Section 4).
# Called by bpm() when projpred = TRUE, and by add_projection().
#
# Strategy: compute the posterior-mean predicted probability at every
# development-sample row, then fit a logistic regression of those soft labels
# on the same design matrix. The resulting coefficients define a simplified
# linear predictor whose inverse-logit approximates the PM predictor.
.compute_projection <- function(X, beta, V, family) {
  mu_dev    <- as.numeric(X %*% beta)
  var_dev   <- rowSums((X %*% V) * X)
  sigma_dev <- sqrt(pmax(var_dev, 0))

  pm_dev <- .pm_quadrature(mu_dev, sigma_dev)

  # Soft labels in (0,1) are intentional; suppress the non-integer warning.
  fit_proj <- suppressWarnings(
    glm.fit(X, pm_dev, family = family)
  )

  list(
    coefficients = fit_proj$coefficients,
    fit_method   = "self_projection"
  )
}

#' Add self-projection to an existing bpm object
#'
#' Computes self-projection coefficients from a fitted `bpm` object without
#' refitting the model. The projection is computed from the already-stored
#' `coefficients` and `vcov`, so this is fast regardless of which prior was
#' used.
#'
#' @param object A `bpm` object.
#' @param data Optional data frame of the development sample. Required when
#'   `object` was fitted with `model = FALSE`; ignored otherwise (the stored
#'   model frame is used).
#'
#' @return The same `bpm` object with the `projection` field populated.
#'   Afterwards, `predict(object, newdata, method = "pm_proj")` and
#'   `coef(object, type = "projection")` will work.
#'
#' @examples
#' set.seed(1)
#' d   <- data.frame(x = rnorm(200), y = rbinom(200, 1, 0.3))
#' fit <- bpm(y ~ x, data = d, prior = log_f(m = 2))   # projpred = FALSE
#' fit <- add_projection(fit)
#' predict(fit, data.frame(x = 0.5), method = "pm_proj")
#'
#' @export
add_projection <- function(object, data = NULL) {
  if (!inherits(object, "bpm"))
    stop("`object` must be a bpm object.", call. = FALSE)

  # Resolve the development-sample design matrix.
  if (!is.null(object$model)) {
    X <- model.matrix(object$terms, object$model,
                      contrasts.arg = object$contrasts)
  } else if (!is.null(data)) {
    tt <- delete.response(object$terms)
    mf <- model.frame(tt, data, xlev = object$xlevels)
    X  <- model.matrix(tt, mf, contrasts.arg = object$contrasts)
  } else {
    stop(
      "No model frame stored on the object. ",
      "Either refit with `model = TRUE` or supply `data`.",
      call. = FALSE
    )
  }

  object$projection <- .compute_projection(
    X, object$coefficients, object$vcov, object$family
  )
  object
}
