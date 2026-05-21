#' Print a bpm object
#' @param x A `bpmfit` object.
#' @param ... Ignored.
#' @export
print.bpmfit <- function(x, ...) {
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
  invisible(x)
}

#' Summarise a bpm object
#'
#' Returns a coefficient table analogous to `summary.glm`.
#'
#' @param object A `bpmfit` object.
#' @param ... Ignored.
#' @return An object of class `"summary.bpmfit"`, printed via `print.summary.bpmfit`.
#' @export
summary.bpmfit <- function(object, ...) {
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
    class = "summary.bpmfit"
  )
}

#' @export
print.summary.bpmfit <- function(x, ...) {
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
#' @param object A `bpmfit` object.
#' @param ... Ignored.
#' @export
coef.bpmfit <- function(object, ...) object$coefficients

#' Extract the posterior covariance matrix from a bpm object
#'
#' @param object A `bpmfit` object.
#' @param ... Ignored.
#' @export
vcov.bpmfit <- function(object, ...) {
  object$vcov
}
