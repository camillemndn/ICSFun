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
