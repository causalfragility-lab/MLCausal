#' Cluster-aware propensity score estimation
#'
#' Estimates propensity scores for clustered observational data using one of
#' three specifications to account for between-cluster confounding.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{"mundlak"}}{Augments the propensity score model with cluster
#'     means of each covariate (Mundlak terms). Softly absorbs between-cluster
#'     confounding without requiring within-cluster treatment variation in
#'     every cluster. \strong{Recommended default.}}
#'   \item{\code{"fixed"}}{Adds cluster indicators as a factor. Fully absorbs
#'     between-cluster confounding but requires within-cluster treatment
#'     variation in every cluster.}
#'   \item{\code{"single"}}{Standard logistic regression with no cluster
#'     adjustment. Useful only as a naive baseline.}
#' }
#'
#' @param data A \code{data.frame} containing all required variables.
#' @param treatment Character string. Name of the binary treatment variable
#'   (must be coded \code{0}/\code{1}).
#' @param covariates Character vector of covariate names.
#' @param cluster Character string. Name of the cluster identifier variable.
#' @param method Character string. One of \code{"mundlak"} (default),
#'   \code{"fixed"}, or \code{"single"}.
#' @param estimand Character string. Target estimand: \code{"ATT"} (default)
#'   or \code{"ATE"}. Stored for downstream use by \code{\link{ml_weight}}
#'   and \code{\link{ml_match}}.
#' @param family A GLM family object. Default \code{stats::binomial()}.
#'
#' @return An object of class \code{ml_ps}, a named list with elements:
#' \describe{
#'   \item{call}{The matched call.}
#'   \item{model}{The fitted \code{\link[stats]{glm}} object.}
#'   \item{data}{Working data frame with a \code{.ps} column of clipped
#'     propensity scores and, when \code{method = "mundlak"}, \code{_cm}
#'     cluster-mean columns.}
#'   \item{ps}{Numeric vector of clipped propensity scores.}
#'   \item{treatment}{Name of the treatment variable.}
#'   \item{covariates}{Names of the covariates used.}
#'   \item{cluster}{Name of the cluster variable.}
#'   \item{method}{PS estimation method used.}
#'   \item{estimand}{Target estimand.}
#' }
#'
#' @seealso \code{\link{ml_weight}}, \code{\link{ml_match}},
#'   \code{\link{plot_overlap_ml}}, \code{\link{balance_ml}}
#'
#' @examples
#' dat <- simulate_ml_data(n_clusters = 5, cluster_size = 10, seed = 1)
#' ps  <- ml_ps(dat, treatment = "z", covariates = c("x1", "x2", "x3"),
#'              cluster = "school_id", method = "mundlak", estimand = "ATT")
#' print(ps)
#'
#' @importFrom stats ave binomial fitted glm as.formula
#' @export
ml_ps <- function(data,
                  treatment,
                  covariates,
                  cluster,
                  method    = c("mundlak", "fixed", "single"),
                  estimand  = c("ATT", "ATE"),
                  family    = stats::binomial()) {
  if (!is.data.frame(data))
    stop("'data' must be a data.frame.", call. = FALSE)
  if (!is.character(treatment) || length(treatment) != 1L)
    stop("'treatment' must be a single character string.", call. = FALSE)
  if (!is.character(covariates) || length(covariates) == 0L)
    stop("'covariates' must be a non-empty character vector.", call. = FALSE)
  if (!is.character(cluster) || length(cluster) != 1L)
    stop("'cluster' must be a single character string.", call. = FALSE)

  method   <- match.arg(method)
  estimand <- match.arg(estimand)

  .check_required_columns(data, c(treatment, covariates, cluster))
  .check_treatment(data[[treatment]], treatment)
  .check_covariates(data, covariates)
  .check_cluster_sizes(data, cluster, min_size = 2L)

  d   <- data
  rhs <- covariates

  if (method == "fixed") {
    d[[".cluster_fe"]] <- factor(d[[cluster]])
    rhs <- c(rhs, ".cluster_fe")
  }

  if (method == "mundlak") {
    for (v in covariates) {
      mname      <- paste0(v, "_cm")
      d[[mname]] <- ave(d[[v]], d[[cluster]],
                        FUN = function(x) mean(x, na.rm = TRUE))
      rhs <- c(rhs, mname)
    }
  }

  fml <- stats::as.formula(
    paste(treatment, "~", paste(rhs, collapse = " + "))
  )
  fit <- stats::glm(fml, data = d, family = family)
  ps  <- .clip01(stats::fitted(fit))
  d[[".ps"]] <- ps

  out <- list(
    call       = match.call(),
    model      = fit,
    data       = d,
    treatment  = treatment,
    covariates = covariates,
    cluster    = cluster,
    method     = method,
    estimand   = estimand,
    ps         = ps
  )
  class(out) <- "ml_ps"
  out
}

#' @export
print.ml_ps <- function(x, ...) {
  cat("MLCausal propensity score model\n")
  cat("  Method:   ", x$method,   "\n", sep = "")
  cat("  Estimand: ", x$estimand, "\n", sep = "")
  cat("  N:        ", nrow(x$data), "\n", sep = "")
  cat("  Clusters: ", length(unique(x$data[[x$cluster]])), "\n", sep = "")
  invisible(x)
}
