#' Simulate clustered observational data
#'
#' Generates a multilevel dataset for prototyping and validating
#' \pkg{MLCausal} workflows. The data-generating process (DGP) includes a
#' cluster-level random effect that simultaneously drives treatment selection
#' and the outcome (i.e. unmeasured cluster-level confounding when
#' \code{method = "single"} is used in \code{\link{ml_ps}}), and
#' heterogeneous treatment effects across clusters.
#'
#' @section True estimands:
#' Because treatment effect heterogeneity correlates with the cluster-level
#' random effect, the true ATT \eqn{\neq} ATE:
#' \itemize{
#'   \item \strong{ATE} \eqn{\approx 0.5} (marginal mean of
#'     \eqn{\tau_j = 0.5 + 0.3 u_j}).
#'   \item \strong{ATT} \eqn{> 0.5} because treated units are
#'     over-represented in high-\eqn{u_j} clusters where \eqn{\tau_j > 0.5}.
#' }
#'
#' @param n_clusters Integer. Number of clusters (schools). Default \code{30}.
#' @param cluster_size Integer. Expected cluster size drawn from
#'   \eqn{\mathrm{Poisson}(\lambda = \code{cluster\_size} - 1) + 1}, then
#'   lower-bounded at \code{n_min}.
#' @param n_min Integer. Minimum guaranteed cluster size. Prevents tiny
#'   clusters that lose entire treatment arms after matching. Default \code{10}.
#' @param seed Optional integer. Random seed for reproducibility.
#'
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{school_id}{Integer cluster identifier (1 to \code{n_clusters}).}
#'   \item{x1}{Continuous individual-level covariate, \eqn{N(0,1)}.}
#'   \item{x2}{Binary individual-level covariate, \eqn{\mathrm{Bernoulli}(0.5)}.}
#'   \item{x3}{Continuous covariate correlated with the cluster random effect.}
#'   \item{z}{Binary treatment indicator (1 = treated, 0 = control).}
#'   \item{y}{Continuous outcome.}
#' }
#'
#' @examples
#' dat <- simulate_ml_data(n_clusters = 5, cluster_size = 10, seed = 1)
#' head(dat)
#' table(dat$z)
#'
#' @export
simulate_ml_data <- function(n_clusters   = 30L,
                             cluster_size = 20L,
                             n_min        = 10L,
                             seed         = NULL) {
  if (!is.numeric(n_clusters) || n_clusters < 2L)
    stop("'n_clusters' must be an integer >= 2.", call. = FALSE)
  if (!is.numeric(cluster_size) || cluster_size < 1L)
    stop("'cluster_size' must be a positive integer.", call. = FALSE)
  if (!is.numeric(n_min) || n_min < 1L)
    stop("'n_min' must be a positive integer.", call. = FALSE)

  if (!is.null(seed)) set.seed(seed)

  sizes     <- pmax(stats::rpois(n_clusters, lambda = cluster_size - 1) + 1L,
                    as.integer(n_min))
  school_id <- rep(seq_len(n_clusters), times = sizes)
  n         <- length(school_id)

  u       <- stats::rnorm(n_clusters, 0, 0.7)
  clust_u <- u[school_id]

  x1 <- stats::rnorm(n)
  x2 <- stats::rbinom(n, 1L, 0.5)
  x3 <- stats::rnorm(n, mean = clust_u)

  lin_ps <- -0.2 + 0.7 * x1 - 0.4 * x2 + 0.5 * x3 + 0.8 * clust_u
  p      <- stats::plogis(lin_ps)
  z      <- stats::rbinom(n, 1L, p)

  tau_cluster <- 0.5 + 0.3 * clust_u
  y <- 0.5 + tau_cluster * z +
       0.8 * x1 - 0.3 * x2 + 0.6 * x3 +
       clust_u + stats::rnorm(n)

  data.frame(school_id = school_id,
             x1 = x1, x2 = x2, x3 = x3,
             z  = z,  y  = y,
             stringsAsFactors = FALSE)
}
