# Generics and methods for likelihood() and posterior().
#
# likelihood() fits the unpenalised (flat-prior) model on-the-fly and returns
# the MLE and observed Fisher information — the data's contribution only,
# independent of the prior.  Useful for multi-centre / federated settings.
#
# posterior() returns the MAP estimate and Laplace-approximated posterior
# covariance (the quantities actually used for prediction).

# ----- generics ---------------------------------------------------------------

#' Extract the likelihood contribution from a model
#'
#' Fits the unpenalised (flat-prior) logistic regression and returns the MLE
#' and observed Fisher information — the data's contribution to inference,
#' independent of the prior. This is the quantity to share in multi-centre or
#' federated settings where sites contribute likelihood information for a later
#' meta-analysis.
#'
#' For the `flat()` prior the likelihood and posterior are identical, so the
#' stored coefficients and vcov are returned directly without refitting.
#'
#' @param object A model object (e.g. `bpm`).
#' @param data Data frame of the development sample. Required when `object` was
#'   fitted with `model = FALSE`; ignored when the model frame is available.
#' @param ... Currently unused.
#' @return A named list with elements `coefficients` (MLE), `vcov`
#'   (inverse observed Fisher information), and `family` (the family object,
#'   carrying the link function and its inverse). The `family` field is
#'   included so that recipients in a meta-analysis know the scale and link on
#'   which the coefficients are expressed, and to support future extension to
#'   other GLMs or survival models.
#' @seealso [posterior()]
#' @export
likelihood <- function(object, ...) UseMethod("likelihood")

#' Extract the posterior from a model
#'
#' Returns the MAP estimate and Laplace-approximated posterior covariance —
#' the quantities used for prediction. These incorporate both the likelihood
#' and the prior.
#'
#' @param object A model object (e.g. `bpm`).
#' @param ... Currently unused.
#' @return A named list with elements `coefficients` (MAP), `vcov`
#'   (posterior covariance), and `family` (the family object, carrying the
#'   link function and its inverse).
#' @seealso [likelihood()]
#' @export
posterior <- function(object, ...) UseMethod("posterior")

# ----- bpm methods ------------------------------------------------------------

#' @rdname likelihood
#' @export
likelihood.bpm <- function(object, data = NULL, ...) {
  # Flat prior: likelihood == posterior, no refitting needed.
  if (inherits(object$prior, "prior_flat"))
    return(list(coefficients = object$coefficients,
                vcov         = object$vcov,
                family       = object$family))

  # Resolve design matrix and response from stored model frame or supplied data.
  if (!is.null(object$model)) {
    X <- model.matrix(object$terms, object$model,
                      contrasts.arg = object$contrasts)
    y <- model.response(object$model)
  } else if (!is.null(data)) {
    tt <- delete.response(object$terms)
    mf <- model.frame(tt, data, xlev = object$xlevels)
    X  <- model.matrix(tt, mf, contrasts.arg = object$contrasts)
    y  <- data[[as.character(object$formula[[2L]])]]
  } else {
    stop(
      "No model frame stored on the object. ",
      "Either refit with `model = TRUE` or supply `data`.",
      call. = FALSE
    )
  }

  flat_fit <- .fit_flat(X, y, object$family)
  list(
    coefficients = flat_fit$coefficients,
    vcov         = flat_fit$vcov,
    family       = object$family
  )
}

#' @export
posterior.bpm <- function(object, ...) {
  list(coefficients = object$coefficients,
       vcov         = object$vcov,
       family       = object$family)
}
