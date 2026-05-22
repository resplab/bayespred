# ---- bpm methods -------------------------------------------------------------

#' Print a bpm posterior object
#'
#' Shows the family/link, a coefficient table with posterior standard errors,
#' and a note that the full covariance matrix is available via `summary()`.
#'
#' @param x A `bpm` object.
#' @param ... Ignored.
#' @export
print.bpm <- function(x, ...) {
  cat("Bayesian Prediction Model posterior (bpm)\n")
  cat("Family:", x$family$family, "/", x$family$link, "\n\n")
  cat("Coefficients:\n")
  print(round(x$coefficients, 4))
  cat("\nPosterior covariance matrix:\n")
  print(round(x$vcov, 6))
  invisible(x)
}

#' Summarise a bpm posterior object
#'
#' Returns a structured summary of the `bpm` object for programmatic access.
#' Printing produces the same output as [print.bpm()].
#'
#' @param object A `bpm` object.
#' @param ... Ignored.
#' @return An object of class `"summary.bpm"` with elements `family`,
#'   `coefficients` (table with Estimate and Std. Error), and `vcov`.
#' @export
summary.bpm <- function(object, ...) {
  se  <- sqrt(diag(object$vcov))
  tab <- cbind(Estimate = object$coefficients, `Std. Error` = se)
  structure(
    list(
      family       = object$family,
      coefficients = tab,
      vcov         = object$vcov
    ),
    class = "summary.bpm"
  )
}

#' @export
print.summary.bpm <- function(x, ...) {
  cat("Bayesian Prediction Model posterior (bpm)\n")
  cat("Family:", x$family$family, "/", x$family$link, "\n\n")
  cat("Coefficients:\n")
  print(round(x$coefficients, 4))
  cat("\nPosterior covariance matrix:\n")
  print(round(x$vcov, 6))
  invisible(x)
}

#' Extract coefficients from a bpm object
#' @param object A `bpm` object.
#' @param ... Ignored.
#' @export
coef.bpm <- function(object, ...) object$coefficients

#' Extract the posterior covariance matrix from a bpm object
#' @param object A `bpm` object.
#' @param ... Ignored.
#' @export
vcov.bpm <- function(object, ...) object$vcov

# ---- bpmfit methods ----------------------------------------------------------

#' Print a bpmfit object
#' @param x A `bpmfit` object.
#' @param ... Ignored.
#' @export
print.bpmfit <- function(x, ...) {
  cat("Bayesian Prediction Model (bpmfit)\n\n")
  cat("Call:\n"); print(x$call); cat("\n")
  cat("Prior: "); print(x$prior)
  cat("Fit method:", x$fit_method, "\n")
  if (!is.null(x$model)) {
    y <- model.response(x$model)
    cat(sprintf("n = %d observations, %d events (%.1f%%)\n",
                length(y), sum(y), 100 * mean(y)))
  }
  cat("\nCoefficients:\n")
  print(round(x$posterior$coefficients, 4))
  invisible(x)
}

#' Summarise a bpmfit object
#'
#' Returns a coefficient table analogous to `summary.glm`.
#'
#' @param object A `bpmfit` object.
#' @param ... Ignored.
#' @return An object of class `"summary.bpmfit"`, printed via `print.summary.bpmfit`.
#' @export
summary.bpmfit <- function(object, ...) {
  beta <- object$posterior$coefficients
  se   <- sqrt(diag(object$posterior$vcov))
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

#' Extract coefficients from a bpmfit object
#'
#' @param object A `bpmfit` object.
#' @param ... Ignored.
#' @export
coef.bpmfit <- function(object, ...) object$posterior$coefficients

#' Extract the posterior covariance matrix from a bpmfit object
#'
#' @param object A `bpmfit` object.
#' @param ... Ignored.
#' @export
vcov.bpmfit <- function(object, ...) object$posterior$vcov
