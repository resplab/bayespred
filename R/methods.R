#' Print a bpm object
#' @param x A `bpm` object.
#' @param ... Ignored.
#' @export
print.bpm <- function(x, ...) {
  cat("Bayesian Prediction Model (bpm)\n\n")
  cat("Call:\n"); print(x$call); cat("\n")
  cat("Prior: "); print(x$prior)
  cat("Fit method:", x$fit_method, "\n")
  if (!is.null(x$model)) {
    y <- model.response(x$model)
    cat(sprintf("n = %d observations, %d events (%.1f%%)\n",
                length(y), sum(y), 100 * mean(y)))
  }
  cat("\nCoefficients:\n")
  print(round(x$coefficients, 4))
  if (!is.null(x$projection))
    cat('\n[Self-projection available: coef(x, type = "projection")]\n')
  invisible(x)
}

#' Summarise a bpm object
#'
#' Returns a coefficient table analogous to `summary.glm`.
#'
#' @param object A `bpm` object.
#' @param ... Ignored.
#' @return An object of class `"summary.bpm"`, printed via `print.summary.bpm`.
#' @export
summary.bpm <- function(object, ...) {
  beta <- object$coefficients
  se   <- sqrt(diag(object$vcov))
  zval <- beta / se
  pval <- 2 * pnorm(-abs(zval))

  coef_table <- cbind(
    Estimate     = beta,
    `Std. Error` = se,
    `z value`    = zval,
    `Pr(>|z|)`   = pval
  )

  y <- if (!is.null(object$model)) model.response(object$model) else NULL

  structure(
    list(
      call         = object$call,
      prior        = object$prior,
      fit_method   = object$fit_method,
      n            = if (!is.null(y)) length(y) else NA_integer_,
      n_events     = if (!is.null(y)) sum(y)    else NA_integer_,
      coefficients = coef_table
    ),
    class = "summary.bpm"
  )
}

#' @export
print.summary.bpm <- function(x, ...) {
  cat("Call:\n"); print(x$call); cat("\n")
  cat("Prior: "); print(x$prior)
  cat("Fit method:", x$fit_method, "\n")
  if (!is.na(x$n))
    cat(sprintf("n = %d observations, %d events\n\n", x$n, x$n_events))
  cat("Coefficients:\n")
  printCoefmat(x$coefficients, digits = 4, signif.stars = TRUE,
               P.values = TRUE, has.Pvalue = TRUE)
  invisible(x)
}

#' Extract coefficients from a bpm object
#'
#' @param object A `bpm` object.
#' @param type Either `"main"` (default) for the fitted coefficients, or
#'   `"projection"` for the self-projection coefficients (requires
#'   `projpred = TRUE` at fit time).
#' @param ... Ignored.
#' @export
coef.bpm <- function(object, type = c("main", "projection"), ...) {
  type <- match.arg(type)
  if (type == "projection") {
    if (is.null(object$projection))
      stop(
        "No projection coefficients found on this object. ",
        "Run `fit <- add_projection(fit)` to compute them without refitting, ",
        "or refit with `projpred = TRUE`.",
        call. = FALSE
      )
    object$projection$coefficients
  } else {
    object$coefficients
  }
}

#' Extract the posterior covariance matrix from a bpm object
#'
#' @param object A `bpm` object.
#' @param ... Ignored.
#' @export
vcov.bpm <- function(object, ...) {
  object$vcov
}
