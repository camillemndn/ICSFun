.mdist_simu_test_memo <- memoise::memoise(ICSOutlier::dist_simu_test)

#' Memoised version of [ICSOutlier::dist_simu_test()]
#'
#' A thin wrapper around [ICSOutlier::dist_simu_test()] that memoises results
#' so repeated calls with identical arguments are returned from cache.
#'
#' @param ... Arguments passed to [ICSOutlier::dist_simu_test()].
#' @return The value returned by [ICSOutlier::dist_simu_test()].
#' @seealso [ICSOutlier::dist_simu_test()]
#' @importFrom memoise memoise
#' @export
mdist_simu_test <- function(...) .mdist_simu_test_memo(...)

#' ICS outlier detection (generic)
#'
#' Generic dispatching on the class of `x`. The default method is
#' [ICSOutlier::ICS_outlier()]; a method is provided for `fd` objects.
#'
#' @param x The input data.
#' @param ... Additional arguments passed to methods.
#' @return An object describing detected outliers (class depends on method).
#' @seealso [ICSOutlier::ICS_outlier()], [ICS_outlier.fd()]
#' @export
ICS_outlier <- function(x, ...) UseMethod("ICS_outlier")

#' @export
#' @rdname ICS_outlier
ICS_outlier.default <- function(x, ...) ICSOutlier::ICS_outlier(x, ...)

#' ICS-based outlier detection for functional data
#'
#' Applies [ICS::ICS()] to functional data (objects of class `fd`) and
#' performs outlier detection in the resulting invariant coordinates.
#'
#' @param x A functional data object of class `fd`.
#' @param S1,S2 Scatter functionals.
#' @param S1_args,S2_args Lists of extra arguments passed to `S1` and `S2`.
#' @param ICS_algorithm Algorithm passed to [ICS::ICS()].
#' @param index Optional integer vector of component indices to use.
#' @param method Component-selection method: `"norm_test"` or `"simulation"`.
#' @param test Normality test used when `method = "norm_test"`.
#' @param n_eig Number of replications for normality-test based selection.
#' @param level_test Significance level for component selection.
#' @param adjust Logical; adjust for multiple testing.
#' @param level_dist Significance level for the outlier cutoff.
#' @param n_dist Number of replications used to compute the cutoff.
#' @param type Type of selection; currently `"smallprop"`.
#' @param n_cores,iseed,pkg,q_type Passed to [ICSOutlier::comp_simu_test()] /
#'   [ICSOutlier::dist_simu_test()].
#' @param ... Additional arguments passed to inner routines.
#' @return An object of class `ICS_Out_fd`.
#' @seealso [plot.ICS_Out_fd()], [ICSOutlier::ICS_outlier()]
#' @method ICS_outlier fd
#' @export
#' @importFrom ggplot2 ggplot geom_point aes geom_line geom_hline ggplotGrob guides
#' @importFrom GGally ggpairs ggmatrix_gtable
ICS_outlier.fd <- function(
  x,
  S1 = ICS::ICS_cov,
  S2 = ICS::ICS_cov4,
  S1_args = list(),
  S2_args = list(),
  ICS_algorithm = c("whiten", "standard", "QR"),
  index = NULL,
  method = "norm_test",
  test = "agostino.test",
  n_eig = 10000,
  level_test = 0.05,
  adjust = TRUE,
  level_dist = 0.025,
  n_dist = 10000,
  type = "smallprop",
  n_cores = NULL,
  iseed = NULL,
  pkg = "ICSOutlier",
  q_type = 7,
  ...
) {
  X <- x
  # Step 1: Checks and apply ICS if necessary
  algorithm <- match.arg(ICS_algorithm)
  method <- match.arg(method, c("norm_test", "simulation"))
  if (!(inherits(X, "fd"))) {
    stop("'x' must be of class 'fd'")
  }
  if (!is.function(S1)) {
    stop(paste("S1 must be specified as a function"))
  }
  if (!is.function(S2)) {
    stop(paste("S2 must be specified as a function"))
  }
  if (inherits(X, "ICS")) {
    warning(paste("'X' already has class 'ICS', not applying ICS"))
    object <- X
  } else {
    object <- tryCatch(
      {
        ICS(
          X,
          S1 = S1,
          S2 = S2,
          S1_args = S1_args,
          S2_args = S2_args,
          algorithm = algorithm,
          center = TRUE,
          fix_signs = "scores"
        )
      },
      warning = function(w) stop(w),
      error = function(e) stop(e)
    )
  }
  rownames(object$scores) <- rownames(X)
  row_names <- rownames(object$scores)

  # Step 2: Select the components
  n <- nrow(object$scores)
  p <- ncol(object$scores)
  type <- match.arg(type, c("smallprop"))
  if (is.null(index)) {
    res_method <- switch(
      method,
      norm_test = {
        ICSOutlier::comp_norm_test(
          object,
          test = test,
          level = level_test,
          adjust = adjust,
          type = type
        )
      },
      simulation = {
        ICSOutlier::comp_simu_test(
          object,
          S1 = S1,
          S2 = S2,
          S1_args = S1_args,
          S2_args = S2_args,
          m = n_eig,
          level = level_test,
          adjust = adjust,
          type = type,
          n_cores = n_cores,
          iseed = iseed,
          pkg = pkg,
          q_type = q_type,
          ...
        )
      }
    )
  } else {
    res_method <- list(index = index)
  }

  # Step 3: Detecting the outliers
  if (sum(res_method$index < 0.5)) {
    outliers <- rep(0L, n)
    names(outliers) <- row_names
    IC_distances <- rep(0L, n)
    names(IC_distances) <- row_names
    IC_distances_quantile <- rep(0, n)
  } else {
    empty_object <- list(scores = matrix(nrow = n, ncol = p))
    class(empty_object) <- "ICS"
    IC_distances_quantile <- mdist_simu_test(
      empty_object,
      S1 = S1,
      S2 = S2,
      S1_args = S1_args,
      S2_args = S2_args,
      m = n_dist,
      index = res_method$index,
      level = level_dist,
      n_cores = n_cores,
      iseed = iseed,
      pkg = pkg,
      q_type = q_type,
      ...
    )
    IC_distances <- ICSOutlier::ics_distances(object, index = res_method$index)
    outliers <- as.integer(IC_distances > IC_distances_quantile)
    names(outliers) <- row_names
  }

  combined <- data.frame(
    index = names(object$gen_kurtosis),
    gen_kurtosis = object$gen_kurtosis,
    selected = seq_len(p) %in% res_method$index
  )
  combined$H <- as.list(object$H)
  combined$H_dual <- as.list(object$H_dual)

  res <- list(
    result = combined,
    gen_kurtosis = object$gen_kurtosis,
    X = X,
    scores = object$scores,
    outliers = outliers,
    ics_distances = IC_distances,
    ics_dist_cutoff = IC_distances_quantile,
    level_dist = level_dist,
    level_test = level_test,
    method = method,
    index = res_method$index,
    test = res_method$test,
    criterion = res_method$criterion,
    adjust = res_method$adjust,
    type = res_method$type,
    n_dist = as.integer(n_dist),
    n_eig = as.integer(n_eig),
    S1_label = object$S1_label,
    S2_label = object$S2_label
  )
  class(res) <- "ICS_Out_fd"
  res
}
