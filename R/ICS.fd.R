#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @rdname ICS
#' @export
ICS <- function(...) UseMethod("ICS")

#' @export
ICS.default <- ICS::ICS

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param fdobj PARAM_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @seealso
#'  \code{\link[fda]{inprod}}, \code{\link[fda]{fd}}
#'  \code{\link[ICS]{ICS-S3}}, \code{\link[ICS]{ics}}, \code{\link[ICS]{ICS}}
#' @rdname ICS
#' @export
#' @importFrom fda inprod fd
#' @importFrom ICS ICS
#' @import ICS
ICS.fd <- function(fdobj, slow = FALSE, ...) {
  # Fix lack of orthonormality of the basis
  changemat <- to_zbsplines(basis = fdobj$basis, inv = TRUE)
  if (slow) {
    gram <- t(changemat) %*% gram(fdobj$basis) %*% changemat
  } else {
    gram <- t(changemat) %*% fda::inprod(fdobj$basis, fdobj$basis) %*% changemat
  }
  # Multivariate ICS
  icsobj <- ICS::ICS(crossprod(to_zbsplines(fdobj), gram), ...)
  # Compute eigenobjects and their dual
  W <- icsobj$W
  icsobj$H <- fda::fd(
    to_zbsplines(coefs = t(W), basis = fdobj$basis, inv = TRUE),
    fdobj$basis
  )
  icsobj$H_dual <- fda::fd(
    to_zbsplines(coefs = solve(W %*% gram), basis = fdobj$basis, inv = TRUE),
    fdobj$basis
  )
  class(icsobj) <- c("ICS", "fd")
  icsobj
}
