#' Invariant coordinate selection for Compositional Data
#'
#' Applies Invariant coordinate selection (ICS) to compositional data,
#' with either isometric log-ratio (ilr) or additive log-ratio (alr) transformations.
#'
#' @param x A compositional dataset of class \code{acomp} (from the \pkg{compositions} package).
#' @param V An ilr basis matrix. Defaults to \code{compositions::ilrBase(x)}.
#' @param transf Character string specifying the log-ratio transformation to use.
#'   Must be one of \code{"ilr"} (default) or \code{"alr"}.
#' @param ... Additional arguments passed to the \code{ICS} function.
#'
#' @return An object of class \code{ICS}, with the unmixing matrix \code{W} back-transformed to the
#' original compositional space and appropriately labeled.
#'
#' @details
#' This function transforms the compositional data using either ilr or alr transformation before
#' applying the \code{ICS} method. After ICA is performed, the unmixing matrix is transformed back to
#' the compositional space using the inverse transformation (\code{ilrInv} or \code{alrInv}).
#'
#' @importFrom compositions ilr ilrBase ilrInv clr alrInv
#' @importFrom ICS ICS
#'
#' @examples
#' \dontrun{
#' library(compositions)
#' data(SimulatedAmounts)
#' result <- ICS.acomp(acomp(SimulatedAmounts), transf = "ilr")
#' }
#'
#' @export
ICS.acomp <- function(x, V = compositions::ilrBase(x),
                      transf = c("alr", "ilr"), ...) {
  # Match transformation argument
  transf <- match.arg(transf)

  # ILR transformation
  if (transf == "ilr") {
    x_ilr <- compositions::ilr(x, V)
    x_ics <- ICS(x_ilr, ...)
    x_ics$W <- compositions::ilrInv(x_ics$W, V)
  }

  # ALR transformation
  if (transf == "alr") {
    x_alr <- compositions::clr(x)[, -ncol(x), drop = FALSE]
    x_ics <- ICS(x_alr, ...)
    x_ics$W <- compositions::alrInv(x_ics$W)
  }

  # Naming consistency
  colnames(x_ics$W) <- colnames(x)
  rownames(x_ics$W) <- colnames(x_ics$scores)

  return(x_ics)
}
