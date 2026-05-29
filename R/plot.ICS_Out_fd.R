#' Plot an ICS outlier-detection result for functional data
#'
#' Diagnostic plot for an `ICS_Out_fd` object produced by
#' [ICS_outlier.fd()]: scree plot of generalized kurtoses, selected
#' eigenfunctions, IC distances with cutoff, pair plot of scores, and the
#' raw functional data colored by outlier status.
#'
#' @param x An `ICS_Out_fd` object.
#' @param ... Additional aesthetic arguments forwarded to ggplot layers.
#' @return Invisibly the combined `grid.arrange` grob; called for its side
#'   effect of drawing plots.
#' @export
#' @method plot ICS_Out_fd
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_hline ggplotGrob
#' @importFrom GGally ggpairs ggmatrix_gtable
#' @importFrom gridExtra grid.arrange
#' @importFrom rlang .data
plot.ICS_Out_fd <- function(x, ...) {
  object <- x
  X <- object$X
  IC_distances <- object$ics_distances
  n <- nrow(object$scores)
  p <- ncol(object$scores)

  # ---- Scree plot of generalized kurtoses ----
  screeplot_dat <- data.frame(
    IC = seq_len(p),
    gen_kurtosis = object$gen_kurtosis,
    selected = factor(
      ifelse(seq_len(p) %in% object$index, "selected", "not selected"),
      levels = c("not selected", "selected")
    )
  )
  g1 <- ggplot2::ggplot(
    screeplot_dat, ggplot2::aes(.data$IC, .data$gen_kurtosis)
  ) +
    ggplot2::geom_line(alpha = 0.5) +
    ggplot2::geom_point(ggplot2::aes(color = .data$selected, size = .data$selected))
  plot(g1)

  # ---- Selected eigenfunctions ----
  eig_idx <- object$index
  eigenfun_list <- as.list(object$H_dual)[eig_idx]
  rangeval <- range(vapply(eigenfun_list, function(f) f$basis$rangeval, numeric(2)))
  grid <- seq(rangeval[1], rangeval[2], length.out = 401)
  eig_long <- do.call(rbind, lapply(seq_along(eigenfun_list), function(k) {
    data.frame(
      x = grid,
      y = as.numeric(fda::eval.fd(grid, eigenfun_list[[k]])),
      IC = factor(eig_idx[k])
    )
  }))
  g2 <- ggplot2::ggplot(
    eig_long, ggplot2::aes(.data$x, .data$y, color = .data$IC, group = .data$IC)
  ) +
    ggplot2::geom_line()
  if (inherits(X, "dd")) {
    g2 <- g2 + ggplot2::geom_hline(yintercept = 1 / diff(X$basis$rangeval))
  }
  plot(g2)

  # ---- IC distances with cutoff ----
  outlier_fact <- factor(
    ifelse(object$outliers == 1, "outlier", "not outlier"),
    levels = c("not outlier", "outlier")
  )
  distances_dat <- data.frame(
    Index = seq_len(n),
    IC_distances = IC_distances,
    outlier = outlier_fact
  )
  g3 <- ggplot2::ggplot(
    distances_dat,
    ggplot2::aes(.data$Index, .data$IC_distances, color = .data$outlier)
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = object$ics_dist_cutoff)

  # ---- Pair plot of selected scores ----
  pairs_dat <- data.frame(object$scores, outlier = outlier_fact)
  g3b <- GGally::ggpairs(
    pairs_dat,
    diag = "blank",
    columns = seq_len(max(max(object$index), 2)),
    ggplot2::aes(shape = NA, color = .data$outlier)
  )
  print(g3b)

  # ---- Raw functional data colored by outlier status (skip if X is ICS) ----
  g4 <- NULL
  if (!inherits(X, "ICS")) {
    X_list <- as.list(X)
    rangeval_X <- X[[1]]$basis$rangeval
    grid_X <- seq(rangeval_X[1], rangeval_X[2], length.out = 401)
    X_long <- do.call(rbind, lapply(seq_along(X_list), function(k) {
      data.frame(
        x = grid_X,
        y = as.numeric(fda::eval.fd(grid_X, X_list[[k]])),
        id = factor(k),
        outlier = outlier_fact[k]
      )
    }))
    g4 <- ggplot2::ggplot(
      X_long,
      ggplot2::aes(
        .data$x, .data$y,
        group = .data$id, color = .data$outlier, alpha = .data$outlier
      )
    ) +
      ggplot2::geom_line()
    plot(g4)
  }

  g <- gridExtra::grid.arrange(
    ggplot2::ggplotGrob(g1),
    ggplot2::ggplotGrob(g2),
    GGally::ggmatrix_gtable(g3b),
    if (!is.null(g4)) ggplot2::ggplotGrob(g4) else NULL,
    ncol = 2
  )
  invisible(g)
}
