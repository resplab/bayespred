#' Project a bpmfit onto a simplified deployable model
#'
#' Unified entry point for all projection methods. Dispatches to the
#' appropriate method based on `type`:
#'
#' * `"pm"` — posterior-mean projection: computes PM soft labels on the
#'   development sample and regresses them onto a design matrix. The result is
#'   a self-contained `bpmproj_pm` object requiring only base R at deployment.
#' * `"full"` — full Bayesian projection (not yet implemented).
#'
#' @param object A `bpmfit` object.
#' @param type Character. Projection method. `"pm"` (default) for
#'   posterior-mean projection; `"full"` is reserved for future use.
#' @param formula Formula for the projection model. `NULL` (default) uses the
#'   same formula as `object` (self-projection). A one- or two-sided formula
#'   (LHS is ignored) overrides this, e.g. `~ age + sex`.
#' @param data Data frame used to compute PM predictions and fit the projection.
#'   Can be the original development sample or a different dataset (e.g. local
#'   data from a new site). When omitted, the model frame stored on `object` is
#'   used (requires `model = TRUE` at fit time).
#' @param family Family object for the projection fit. `NULL` (default) reuses
#'   `object$posterior$family`.
#' @param compute_divergence Logical (default `TRUE`). If `TRUE`, compute the
#'   median KL divergence and median Jensen-Shannon divergence (both in bits)
#'   between the full model's PM predictions and the projected model's fitted
#'   values on the development sample. These are stored as `$kl_median` and
#'   `$js_median` on the returned object.
#' @param ... Currently unused.
#'
#' @return For `type = "pm"`, an object of class `"bpmproj_pm"` with elements
#'   `coefficients`, `terms`, `contrasts`, `xlevels`, `family`, `call`,
#'   `kl_median`, and `js_median`.
#' @seealso [predict.bpmproj_pm()]
#' @export
bpmproject <- function(object, type = c("pm", "full"),
                       formula = NULL, data = NULL, family = NULL,
                       compute_divergence = TRUE, ...) {
  type <- match.arg(type)
  switch(type,
    pm   = .bpmproj_pm(object, formula = formula, data = data,
                       family = family,
                       compute_divergence = compute_divergence, ...),
    full = stop('type = "full" projection is not yet implemented.', call. = FALSE)
  )
}

# ------------------------------------------------------------------------------
# Internal helpers

# Per-observation KL divergence D_KL(p || q) in bits.
.kl_bits <- function(p, q) {
  p * log2(p / q) + (1 - p) * log2((1 - p) / (1 - q))
}

# Per-observation Jensen-Shannon divergence in bits (symmetric, bounded [0,1]).
.js_bits <- function(p, q) {
  m <- (p + q) / 2
  0.5 * .kl_bits(p, m) + 0.5 * .kl_bits(q, m)
}

# Compute median KL and JS divergence; replace non-finite values with NA and
# warn if any are found.
.divergence_summary <- function(pm, q_fitted) {
  kl  <- .kl_bits(pm, q_fitted)
  js  <- .js_bits(pm, q_fitted)
  bad <- !is.finite(kl) | !is.finite(js)
  if (any(bad)) {
    kl[bad] <- NA_real_
    js[bad] <- NA_real_
    warning(
      sum(bad), " observation(s) produced non-finite divergence values ",
      "(projected probabilities outside (0, 1)); set to NA.",
      call. = FALSE
    )
  }
  list(kl_median = median(kl, na.rm = TRUE),
       js_median = median(js, na.rm = TRUE))
}

# ------------------------------------------------------------------------------
# Internal PM projection implementation

.bpmproj_pm <- function(object, formula = NULL, data = NULL, family = NULL,
                         compute_divergence = TRUE, ...) {
  cl          <- match.call()
  post        <- object$posterior
  proj_family <- if (is.null(family)) post$family else family

  if (is.null(formula)) {
    # Self-projection -----------------------------------------------------------
    if (!is.null(object$model)) {
      X_dev <- model.matrix(post$terms, object$model,
                            contrasts.arg = post$contrasts)
    } else if (!is.null(data)) {
      tt    <- delete.response(post$terms)
      mf    <- model.frame(tt, data, xlev = post$xlevels)
      X_dev <- model.matrix(tt, mf, contrasts.arg = post$contrasts)
    } else {
      stop(
        "No model frame stored on this object. ",
        "Supply `data` or refit with `model = TRUE`.",
        call. = FALSE
      )
    }
    mu  <- as.numeric(X_dev %*% post$coefficients)
    sig <- sqrt(pmax(as.numeric(rowSums((X_dev %*% post$vcov) * X_dev)), 0))
    pm  <- .pm_quadrature(mu, sig)
    fit <- suppressWarnings(glm.fit(X_dev, pm, family = proj_family))

    div <- if (isTRUE(compute_divergence))
      .divergence_summary(pm, fit$fitted.values)
    else
      list(kl_median = NA_real_, js_median = NA_real_)

    structure(
      list(
        coefficients = fit$coefficients,
        terms        = post$terms,
        contrasts    = post$contrasts,
        xlevels      = post$xlevels,
        family       = proj_family,
        kl_median    = div$kl_median,
        js_median    = div$js_median,
        call         = cl
      ),
      class = "bpmproj_pm"
    )

  } else {
    # Custom projection ---------------------------------------------------------
    if (is.null(data))
      stop("Supply `data` for custom projection.", call. = FALSE)

    tt_dev <- delete.response(post$terms)
    mf_dev <- model.frame(tt_dev, data, xlev = post$xlevels)
    X_dev  <- model.matrix(tt_dev, mf_dev, contrasts.arg = post$contrasts)

    mu  <- as.numeric(X_dev %*% post$coefficients)
    sig <- sqrt(pmax(as.numeric(rowSums((X_dev %*% post$vcov) * X_dev)), 0))
    pm  <- .pm_quadrature(mu, sig)

    if (length(formula) == 3L) formula <- formula[-2L]
    proj_terms     <- terms(formula, data = data)
    mf_proj        <- model.frame(proj_terms, data)
    X_proj         <- model.matrix(proj_terms, mf_proj)
    proj_contrasts <- attr(X_proj, "contrasts")
    proj_xlevels   <- .get_xlevels(proj_terms, mf_proj)

    fit <- suppressWarnings(glm.fit(X_proj, pm, family = proj_family))

    div <- if (isTRUE(compute_divergence))
      .divergence_summary(pm, fit$fitted.values)
    else
      list(kl_median = NA_real_, js_median = NA_real_)

    structure(
      list(
        coefficients = fit$coefficients,
        terms        = proj_terms,
        contrasts    = proj_contrasts,
        xlevels      = proj_xlevels,
        family       = proj_family,
        kl_median    = div$kl_median,
        js_median    = div$js_median,
        call         = cl
      ),
      class = "bpmproj_pm"
    )
  }
}

# ------------------------------------------------------------------------------

#' Predict from a PM-projected model
#'
#' @param object A `bpmproj_pm` object returned by [bpmproject()].
#' @param newdata A data frame of new observations.
#' @param type `"response"` (default) for predicted probability;
#'   `"link"` for the linear predictor.
#' @param na.action Function to handle `NA`s in `newdata`. Default
#'   [stats::na.pass].
#' @param ... Currently unused.
#' @return A numeric vector.
#' @seealso [bpmproject()]
#' @export
predict.bpmproj_pm <- function(object, newdata,
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

#' Extract coefficients from a bpmproj_pm object
#' @param object A `bpmproj_pm` object.
#' @param ... Ignored.
#' @export
coef.bpmproj_pm <- function(object, ...) object$coefficients

#' Print a bpmproj_pm object
#' @param x A `bpmproj_pm` object.
#' @param ... Ignored.
#' @export
print.bpmproj_pm <- function(x, ...) {
  cat("PM-Projected Prediction Model\n\n")
  cat("Call:\n"); print(x$call); cat("\n")
  cat("Coefficients:\n")
  print(round(x$coefficients, 4))
  if (!is.na(x$kl_median)) {
    cat("\nProjection accuracy (development sample):\n")
    cat(sprintf("  Median KL divergence: %.5f bits\n", x$kl_median))
    cat(sprintf("  Median JS divergence: %.5f bits\n", x$js_median))
  }
  invisible(x)
}
