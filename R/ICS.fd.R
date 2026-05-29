#' Invariant Coordinate Selection (generic)
#'
#' Generic dispatching on the class of `x`. The default method is
#' [ICS::ICS()]; methods are provided for `fd` and `acomp` objects.
#'
#' @param x The input data. The default method forwards to [ICS::ICS()].
#' @param ... Additional arguments passed to methods.
#' @return An object of class `ICS`.
#' @seealso [ICS::ICS()], [ICS.fd()], [ICS.acomp()]
#' @export
ICS <- function(x, ...) UseMethod("ICS")

#' @export
#' @rdname ICS
ICS.default <- function(x, ...) ICS::ICS(x, ...)

#' Invariant coordinate selection for Functional Data
#'
#' Applies Invariant coordinate selection (ICA) to functional data objects
#' (\code{fd} objects from the \pkg{fda} package) using the \code{ICS} method.
#'
#' @param x A functional data object of class \code{fd}.
#' @param slow Logical; if \code{TRUE}, computes the Gram matrix using \code{gram()}.
#'   If \code{FALSE} (default), uses \code{fda::inprod()} for efficiency.
#' @param ... Additional arguments passed to the \code{ICS} function.
#'
#' @return An object of class \code{ICS} and \code{fd}, with fields:
#'   \describe{
#'     \item{\code{W}}{Unmixing matrix in coefficient space}
#'     \item{\code{H}}{Basis functions (eigenfunctions)}
#'     \item{\code{H_dual}}{Dual eigenbasis of functions}
#'     \item{\code{scores}}{Scores of the components (inherited from \code{ICS})}
#'   }
#'
#' @details
#' The function projects the functional data into an approximately orthonormal Z-basis
#' using a transformation matrix. ICA is performed on the transformed coefficients,
#' and the results are projected back into function space as eigenfunctions.
#'
#' The dual eigenbasis \code{H_dual} corresponds to the basis of functions used to reconstruct the original data from the invariant coordinates.
#' in the metric induced by the Gram matrix.
#'
#' @importFrom fda fd inprod
#' @importFrom ICS ICS
#' @seealso \code{\link[fda]{fd}}, \code{\link[ICS]{ICS}}
#'
#' @export
ICS.fd <- function(x, slow = FALSE, ...) {
  fdobj <- x
  # Change of basis between B-splines and ZB-splines
  change_mat <- to_zbsplines(basis = fdobj$basis, inv = TRUE)

  # Compute Gram matrix of ZB-spline basis
  gram <- if (slow) {
    t(change_mat) %*% gram(fdobj$basis) %*% change_mat
  } else {
    t(change_mat) %*% fda::inprod(fdobj$basis, fdobj$basis) %*% change_mat
  }

  # Get transformed coefficients in the ZB-spline basis
  coef_z <- to_zbsplines(fdobj)

  # Apply ICS on crossproduct of coefficients with Gram matrix
  icsobj <- ICS::ICS(crossprod(coef_z, gram), ...)

  # Extract unmixing matrix and compute eigenfunctions
  W <- icsobj$W
  basis <- fdobj$basis
  icsobj$H <- fda::fd(
    to_zbsplines(coefs = t(W), basis = basis, inv = TRUE),
    basis
  )

  # Compute dual eigenbasis for reconstruction
  icsobj$H_dual <- fda::fd(
    to_zbsplines(coefs = solve(W %*% gram), basis = basis, inv = TRUE),
    basis
  )

  class(icsobj) <- c("ICS", "fd")
  return(icsobj)
}
