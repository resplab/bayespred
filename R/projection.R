#' Project a bpm model onto a simplified linear predictor (PM method)
#'
#' Computes posterior-mean (PM) predictions on the development sample and
#' regresses them onto a design matrix, producing a standalone `bpm_proj_pm`
#' object with no covariance matrix. All three aspects of the projection model
#' default to the main `bpm` fit but can be overridden independently:
#'
#' * **`formula`** — defaults to the main model's formula (self-projection).
#'   Supply a simpler or different formula to project onto a reduced predictor
#'   set, e.g. `~ age + sex`.
#' * **`family`** — defaults to `object$family` (same link as the main model).
#'   Override to project onto a different link function.
#' * **`data`** — the sample used both to compute PM predictions (always via
#'   the main model's coefficients and vcov) and to fit the projection formula.
#'   Can be the original development sample, or a different dataset entirely —
#'   for example, local data from a new site, enabling domain-specific
#'   projection without refitting the main model. When omitted, the model frame
#'   stored on `object` is used (requires `model = TRUE` at fit time).
#'
#' The returned object is self-contained and can be saved and deployed with
#' [predict.bpm_proj_pm()] using base R only — no dependency on the original
#' fit or on `mgcv` / `brglm2`.
#'
#' @param object A `bpm` object.
#' @param formula Formula for the projection model. `NULL` (default) uses the
#'   same formula as `object`. A one- or two-sided formula (LHS is ignored)
#'   overrides this, e.g. `~ age + sex`.
#' @param data Data frame used to compute PM predictions and fit the projection.
#'   Can be the original development sample or a different dataset (e.g. local
#'   data from a new site). When omitted, the model frame stored on `object` is
#'   used (requires `model = TRUE` at fit time).
#' @param family Family object for the projection fit. `NULL` (default) reuses
#'   `object$family`.
#' @param ... Currently unused.
#'
#' @return An object of class `"bpm_proj_pm"` with elements `coefficients`,
#'   `terms`, `contrasts`, `xlevels`, `family`, and `call`.
#'
#' @seealso [predict.bpm_proj_pm()]
#' @export
project_pm <- function(object, ...) UseMethod("project_pm")

#' @export
project_pm.bpm <- function(object, formula = NULL, data = NULL,
                            family = NULL, ...) {
  cl          <- match.call()
  proj_family <- if (is.null(family)) object$family else family

  if (is.null(formula)) {
    # Self-projection -----------------------------------------------------------
    if (!is.null(object$model)) {
      X_dev <- model.matrix(object$terms, object$model,
                            contrasts.arg = object$contrasts)
    } else if (!is.null(data)) {
      tt    <- delete.response(object$terms)
      mf    <- model.frame(tt, data, xlev = object$xlevels)
      X_dev <- model.matrix(tt, mf, contrasts.arg = object$contrasts)
    } else {
      stop(
        "No model frame stored on this object. ",
        "Supply `data` or refit with `model = TRUE`.",
        call. = FALSE
      )
    }
    mu  <- as.numeric(X_dev %*% object$coefficients)
    sig <- sqrt(pmax(as.numeric(rowSums((X_dev %*% object$vcov) * X_dev)), 0))
    pm  <- .pm_quadrature(mu, sig)
    fit <- suppressWarnings(glm.fit(X_dev, pm, family = proj_family))

    structure(
      list(
        coefficients = fit$coefficients,
        terms        = object$terms,
        contrasts    = object$contrasts,
        xlevels      = object$xlevels,
        family       = proj_family,
        call         = cl
      ),
      class = "bpm_proj_pm"
    )

  } else {
    # Custom projection ---------------------------------------------------------
    if (is.null(data))
      stop("Supply `data` for custom projection.", call. = FALSE)

    tt_dev <- delete.response(object$terms)
    mf_dev <- model.frame(tt_dev, data, xlev = object$xlevels)
    X_dev  <- model.matrix(tt_dev, mf_dev, contrasts.arg = object$contrasts)

    mu  <- as.numeric(X_dev %*% object$coefficients)
    sig <- sqrt(pmax(as.numeric(rowSums((X_dev %*% object$vcov) * X_dev)), 0))
    pm  <- .pm_quadrature(mu, sig)

    if (length(formula) == 3L) formula <- formula[-2L]
    proj_terms    <- terms(formula, data = data)
    mf_proj       <- model.frame(proj_terms, data)
    X_proj        <- model.matrix(proj_terms, mf_proj)
    proj_contrasts <- attr(X_proj, "contrasts")
    proj_xlevels  <- .get_xlevels(proj_terms, mf_proj)

    fit <- suppressWarnings(glm.fit(X_proj, pm, family = proj_family))

    structure(
      list(
        coefficients = fit$coefficients,
        terms        = proj_terms,
        contrasts    = proj_contrasts,
        xlevels      = proj_xlevels,
        family       = proj_family,
        call         = cl
      ),
      class = "bpm_proj_pm"
    )
  }
}

# ------------------------------------------------------------------------------

#' Predict from a PM-projected model
#'
#' @param object A `bpm_proj_pm` object returned by [project_pm()].
#' @param newdata A data frame of new observations.
#' @param type `"response"` (default) for predicted probability;
#'   `"link"` for the linear predictor.
#' @param na.action Function to handle `NA`s in `newdata`. Default
#'   [stats::na.pass].
#' @param ... Currently unused.
#' @return A numeric vector.
#' @seealso [project_pm()]
#' @export
predict.bpm_proj_pm <- function(object, newdata,
                                 type      = c("response", "link"),
                                 na.action = na.pass, ...) {
  type <- match.arg(type)
  tt   <- delete.response(object$terms)
  mf   <- model.frame(tt, newdata, xlev = object$xlevels, na.action = na.action)
  X    <- model.matrix(tt, mf, contrasts.arg = object$contrasts)
  eta  <- as.numeric(X %*% object$coefficients)
  if (type == "link") return(eta)
  as.numeric(object$family$linkinv(eta))
}

#' Extract coefficients from a bpm_proj_pm object
#' @param object A `bpm_proj_pm` object.
#' @param ... Ignored.
#' @export
coef.bpm_proj_pm <- function(object, ...) object$coefficients

#' Print a bpm_proj_pm object
#' @param x A `bpm_proj_pm` object.
#' @param ... Ignored.
#' @export
print.bpm_proj_pm <- function(x, ...) {
  cat("PM-Projected Prediction Model\n\n")
  cat("Call:\n"); print(x$call); cat("\n")
  cat("Coefficients:\n")
  print(round(x$coefficients, 4))
  invisible(x)
}
