# Internal: compute projection coefficients.
# X_dev:  design matrix used to compute PM on the development sample.
# X_proj: design matrix the projection coefficients are fitted against.
#         Identical to X_dev for self-projection; different for custom projection.
.compute_projection <- function(X_dev, X_proj, beta, V, family, self = TRUE) {
  mu    <- as.numeric(X_dev %*% beta)
  sigma <- sqrt(pmax(as.numeric(rowSums((X_dev %*% V) * X_dev)), 0))
  pm    <- .pm_quadrature(mu, sigma)
  fit   <- suppressWarnings(glm.fit(X_proj, pm, family = family))
  list(
    coefficients = fit$coefficients,
    fit_method   = if (self) "self_projection" else "custom_projection"
  )
}

#' Compute projection coefficients from a bpm object
#'
#' Returns the projection result without modifying the object. For
#' self-projection (no `formula`), the posterior-mean predictions on the
#' development sample are regressed back on the same design matrix, producing
#' a simplified linear predictor. For custom projection, they are regressed on
#' a user-supplied (typically simpler) formula — useful for building clinical
#' risk scores with fewer or different predictors.
#'
#' To cache the result on the object for use with
#' `predict(..., method = "pm_proj")`, pass the return value to
#' [add_projection()] or use it directly in a pipe:
#' ```r
#' fit <- bpm(...) |> add_projection()
#' fit <- bpm(...) |> add_projection(formula = ~ age + sex, data = dat)
#' ```
#'
#' @param object A `bpm` object.
#' @param formula Optional one- or two-sided formula specifying the projection
#'   model. If `NULL` (default), self-projection onto the original design matrix
#'   is performed.
#' @param data Data frame containing the development sample. Required for
#'   custom projection; also required for self-projection when `object` was
#'   fitted with `model = FALSE`.
#' @param family Family object for the projection fit. Defaults to
#'   `object$family` (same link as the main model).
#' @param ... Currently unused.
#'
#' @return A list with elements `coefficients`, `fit_method`, `terms`,
#'   `contrasts`, `xlevels`, and `family` — everything needed for
#'   `predict(..., method = "pm_proj")` to reconstruct the design matrix on
#'   new data.
#'
#' @seealso [add_projection()], [predict.bpm()]
#' @export
project <- function(object, ...) UseMethod("project")

#' @export
project.bpm <- function(object, formula = NULL, data = NULL, family = NULL, ...) {
  proj_family <- if (is.null(family)) object$family else family

  if (is.null(formula)) {
    # Self-projection ---------------------------------------------------------
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
    raw <- .compute_projection(X_dev, X_dev, object$coefficients, object$vcov,
                               proj_family, self = TRUE)
    c(raw, list(
      terms     = object$terms,
      contrasts = object$contrasts,
      xlevels   = object$xlevels,
      family    = proj_family
    ))

  } else {
    # Custom projection -------------------------------------------------------
    if (is.null(data))
      stop("Supply `data` for custom projection.", call. = FALSE)

    # PM predictions are computed from the main model on `data`.
    tt_dev <- delete.response(object$terms)
    mf_dev <- model.frame(tt_dev, data, xlev = object$xlevels)
    X_dev  <- model.matrix(tt_dev, mf_dev, contrasts.arg = object$contrasts)

    # Projection design matrix from the custom formula.
    if (length(formula) == 3L) formula <- formula[-2L]  # drop LHS if present
    proj_terms    <- terms(formula, data = data)
    mf_proj       <- model.frame(proj_terms, data)
    X_proj        <- model.matrix(proj_terms, mf_proj)
    proj_contrasts <- attr(X_proj, "contrasts")
    proj_xlevels  <- .get_xlevels(proj_terms, mf_proj)

    raw <- .compute_projection(X_dev, X_proj, object$coefficients, object$vcov,
                               proj_family, self = FALSE)
    c(raw, list(
      terms     = proj_terms,
      contrasts = proj_contrasts,
      xlevels   = proj_xlevels,
      family    = proj_family
    ))
  }
}

#' Cache projection coefficients on a bpm object
#'
#' Computes projection coefficients via [project()] and stores them on the
#' object. Pipe-friendly: takes and returns the `bpm` object.
#'
#' ```r
#' fit <- bpm(...) |> add_projection()
#' fit <- bpm(...) |> add_projection(formula = ~ age + sex, data = dat)
#' ```
#'
#' @inheritParams project
#' @return The same `bpm` object with `$projection` populated.
#' @seealso [project()], [predict.bpm()]
#' @export
add_projection <- function(object, ...) UseMethod("add_projection")

#' @export
add_projection.bpm <- function(object, formula = NULL, data = NULL,
                                family = NULL, ...) {
  object$projection <- project(object, formula = formula, data = data,
                                family = family, ...)
  object
}
