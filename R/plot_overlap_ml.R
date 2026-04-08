#' Plot propensity score overlap between treatment groups
#'
#' Creates a density overlap plot of estimated propensity scores by treatment
#' group. Substantial overlap supports the positivity assumption. The plot can
#' be shown overall or faceted by cluster.
#'
#' @param x An \code{ml_ps} object from \code{\link{ml_ps}}, or a plain
#'   \code{data.frame} (in which case \code{treatment}, \code{cluster}, and
#'   \code{ps} must be supplied).
#' @param treatment Character string. Treatment variable name. Required when
#'   \code{x} is a \code{data.frame}.
#' @param cluster Character string. Cluster variable name. Required when
#'   \code{x} is a \code{data.frame}.
#' @param ps Character string. Propensity score variable name. Required when
#'   \code{x} is a \code{data.frame}.
#' @param facet_clusters Logical. Facet by cluster? Default \code{FALSE}.
#' @param top_n_clusters Positive integer. Maximum clusters shown when
#'   \code{facet_clusters = TRUE}. Default \code{12}.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @seealso \code{\link{ml_ps}}
#'
#' @examples
#' dat <- simulate_ml_data(n_clusters = 5, cluster_size = 10, seed = 1)
#' ps  <- ml_ps(dat, "z", c("x1", "x2", "x3"), "school_id")
#' plot_overlap_ml(ps)
#'
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_density labs theme_minimal facet_wrap
#' @importFrom stats as.formula
#' @export
plot_overlap_ml <- function(x,
                            treatment      = NULL,
                            cluster        = NULL,
                            ps             = NULL,
                            facet_clusters = FALSE,
                            top_n_clusters = 12L) {
  if (inherits(x, "ml_ps")) {
    d         <- x$data
    treatment <- x$treatment
    cluster   <- x$cluster
    ps        <- ".ps"
  } else {
    if (!is.data.frame(x))
      stop("'x' must be an 'ml_ps' object or a data.frame.", call. = FALSE)
    if (is.null(treatment) || is.null(cluster) || is.null(ps))
      stop("When 'x' is a data.frame, 'treatment', 'cluster', and 'ps' ",
           "must all be supplied.", call. = FALSE)
    d <- x
    .check_required_columns(d, c(treatment, cluster, ps))
  }

  d[[treatment]] <- factor(d[[treatment]], levels = c(0, 1),
                           labels = c("Control", "Treated"))

  if (facet_clusters) {
    tab  <- sort(table(d[[cluster]]), decreasing = TRUE)
    keep <- names(tab)[seq_len(min(as.integer(top_n_clusters), length(tab)))]
    d    <- d[d[[cluster]] %in% keep, , drop = FALSE]
  }

  p <- ggplot2::ggplot(d, ggplot2::aes(x    = .data[[ps]],
                                        fill = .data[[treatment]])) +
    ggplot2::geom_density(alpha = 0.35) +
    ggplot2::labs(x = "Propensity score", y = "Density", fill = "Group") +
    ggplot2::theme_minimal()

  if (facet_clusters)
    p <- p + ggplot2::facet_wrap(
      stats::as.formula(paste("~", cluster)), scales = "free_y"
    )
  p
}
