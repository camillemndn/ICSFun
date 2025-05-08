gram <- memoise::memoise(function(bobj) {
  rval <- bobj$rangeval
  p <- bobj$nbasis
  f <- function(t) t(fda::eval.basis(t, bobj)) %*% fda::eval.basis(t, bobj)

  result_matrix <- matrix(0, nrow = p, ncol = p)

  for (i in 1:p) {
    for (j in 1:p) {
      element_function <- Vectorize(function(x) {
        f(x)[i, j]
      })
      result_matrix[i, j] <- stats::integrate(
        element_function,
        rval[1], rval[2]
      )$value
    }
  }
  result_matrix
})

to_zbsplines <- function(
    fdobj = NULL,
    coefs = fdobj$coefs, basis = fdobj$basis, inv = FALSE) {
  rangeval <- basis$rangeval
  knots <- basis$params
  g <- length(knots)
  p <- basis$nbasis
  d <- p - g
  a <- min(rangeval)
  b <- max(rangeval)
  extknots <- c(rep(a, d), knots, rep(b, d))

  dinvmat <- diag(diff(extknots, lag = d)) / d
  dmat <- solve(dinvmat)

  kinvmat <- diag(p)[-p, ]
  kinvmat[lower.tri(kinvmat)] <- 1
  kmat <- diag(p)[, -p]
  kmat[cbind(2:p, 1:(p - 1))] <- -1

  changemat <- if (inv) dmat %*% kmat else kinvmat %*% dinvmat
  if (is.null(coefs)) changemat else changemat %*% coefs
}
