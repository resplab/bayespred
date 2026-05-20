# Generic + methods for likelihood(), posterior(), and add_likelihood().
#
# likelihood() returns the unpenalised MLE and observed Fisher information
# (the data's contribution only — shareable for meta-analysis / federated
# learning without exposing the prior choice).
#
# posterior() returns the MAP estimate and Laplace-approximated posterior
# covariance (prior x likelihood — the quantities used for prediction).
#
# add_likelihood() is the pipe-friendly mutator that fits the unpenalised model,
# caches the result on the bpm object, and returns the updated object.

# ----- generics ---------------------------------------------------------------

#' Extract the likelihood contribution from a model
#'
#' Returns the unpenalised MLE and observed Fisher information — the data's
#' contribution to inference, independent of the prior. This is the quantity
#' to share in multi-centre or federated settings where sites contribute
#' likelihood information for a later meta-analysis.
#'
#' For the `flat()` prior the likelihood and posterior are identical, so no
#' extra fitting is needed.  For all other priors the unpenalised fit must be
#' cached first via [add_likelihood()].
#'
#' @param object A model object (e.g. `bpm`).
#' @param ... Currently unused.
#' @return A named list with elements `coefficients` (MLE), `vcov`
#'   (inverse observed Fisher information), and `family` (the family object,
#'   carrying the link function and its inverse). The `family` field is
#'   included so that recipients in a meta-analysis know the scale and link on
#'   which the coefficients are expressed, and to support future extension to
#'   other GLMs or survival models.
#' @seealso [add_likelihood()], [posterior()]
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
#' @seealso [likelihood()], [add_likelihood()]
#' @export
posterior <- function(object, ...) UseMethod("posterior")

#' Cache the likelihood fit on a bpm object
#'
#' Fits the unpenalised (flat-prior) logistic regression, caches the MLE and
#' observed Fisher information on the object, and returns the updated object.
#' After calling this, [likelihood()] works instantly with no recomputation.
#'
#' This function is pipe-friendly: the object is the first argument and the
#' updated object is returned, so it slots naturally into a `|>` chain:
#'
#' ```r
#' fit <- bpm(...) |> add_likelihood() |> add_projection()
#' ```
#'
#' For the `flat()` prior this is a no-op (likelihood == posterior already);
#' the object is returned unchanged.
#'
#' @param object A `bpm` object.
#' @param data Optional data frame of the development sample. Required when
#'   `object` was fitted with `model = FALSE`; ignored otherwise.
#' @param ... Currently unused.
#' @return The same `bpm` object with the `likelihood_fit` field populated.
#' @seealso [likelihood()], [posterior()], [add_projection()]
#' @export
add_likelihood <- function(object, ...) UseMethod("add_likelihood")

# ----- bpm methods ------------------------------------------------------------

#' @export
likelihood.bpm <- function(object, ...) {
  # Flat prior: likelihood == posterior, no extra fit needed.
  if (inherits(object$prior, "prior_flat"))
    return(list(coefficients = object$coefficients,
                vcov         = object$vcov,
                family       = object$family))

  if (is.null(object$likelihood_fit))
    stop(
      "Likelihood fit not cached on this object. ",
      "Call `fit <- add_likelihood(fit)` first.",
      call. = FALSE
    )

  object$likelihood_fit
}

#' @export
posterior.bpm <- function(object, ...) {
  list(coefficients = object$coefficients,
       vcov         = object$vcov,
       family       = object$family)
}

#' @export
add_likelihood.bpm <- function(object, data = NULL, ...) {
  # No-op for flat prior: likelihood and posterior are already the same.
  if (inherits(object$prior, "prior_flat"))
    return(object)

  # Resolve design matrix from stored model frame or supplied data.
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
  object$likelihood_fit <- list(
    coefficients = flat_fit$coefficients,
    vcov         = flat_fit$vcov,
    family       = object$family
  )
  object
}
