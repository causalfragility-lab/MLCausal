#' Treatment effect estimation with cluster-robust standard errors
#'
#' Fits a (weighted) linear outcome model and reports the treatment effect
#' estimate with cluster-robust (HC1) standard errors via the sandwich
#' estimator. When \code{covariates} and \code{weights} are both supplied the
#' estimator is doubly robust.
#'
#' @param data A \code{data.frame}. Typically the output of
#'   \code{\link{ml_weight}} or the \code{data_matched} element from
#'   \code{\link{ml_match}}.
#' @param outcome Character string. Name of the continuous outcome variable.
#' @param treatment Character string. Name of the binary treatment variable
#'   (\code{0}/\code{1}).
#' @param cluster Character string. Name of the cluster identifier variable.
#' @param covariates Optional character vector of regression covariates for
#'   doubly robust adjustment. \code{NULL} (default) fits treatment only.
#' @param weights Optional character string. Name of a weight variable in
#'   \code{data} (e.g. \code{"weights"} or \code{"match_weight"}).
#'   \code{NULL} (default) fits an unweighted model.
#'
#' @return An object of class \code{estimate_att_ml}:
#' \describe{
#'   \item{call}{The matched call.}
#'   \item{model}{Fitted \code{\link[stats]{lm}} object.}
#'   \item{vcov}{Cluster-robust variance-covariance matrix (HC1).}
#'   \item{coef_table}{\code{\link[lmtest]{coeftest}} coefficient table.}
#'   \item{estimate}{Point estimate for the treatment coefficient.}
#'   \item{se}{Cluster-robust standard error.}
#'   \item{p_value}{Two-sided p-value.}
#' }
#'
#' @seealso \code{\link{ml_weight}}, \code{\link{ml_match}},
#'   \code{\link{sens_ml}}
#'
#' @examples
#' dat   <- simulate_ml_data(n_clusters = 5, cluster_size = 10, seed = 1)
#' ps    <- ml_ps(dat, "z", c("x1", "x2", "x3"), "school_id")
#' dat_w <- ml_weight(ps)
#' est   <- estimate_att_ml(dat_w, outcome = "y", treatment = "z",
#'                          cluster = "school_id", weights = "weights")
#' print(est)
#'
#' @importFrom stats lm as.formula
#' @importFrom sandwich vcovCL
#' @importFrom lmtest coeftest
#' @export
estimate_att_ml <- function(data,
                            outcome,
                            treatment,
                            cluster,
                            covariates = NULL,
                            weights    = NULL) {
  if (!is.data.frame(data))
    stop("'data' must be a data.frame.", call. = FALSE)
  if (!is.character(outcome)   || length(outcome)   != 1L)
    stop("'outcome' must be a single character string.", call. = FALSE)
  if (!is.character(treatment) || length(treatment) != 1L)
    stop("'treatment' must be a single character string.", call. = FALSE)
  if (!is.character(cluster)   || length(cluster)   != 1L)
    stop("'cluster' must be a single character string.", call. = FALSE)
  if (!is.null(covariates) && !is.character(covariates))
    stop("'covariates' must be a character vector or NULL.", call. = FALSE)
  if (!is.null(weights) &&
      (!is.character(weights) || length(weights) != 1L))
    stop("'weights' must be a single character string or NULL.", call. = FALSE)

  .check_required_columns(data, c(outcome, treatment, cluster))
  if (!is.null(covariates)) .check_required_columns(data, covariates)
  if (!is.null(weights))    .check_required_columns(data, weights)
  .check_treatment(data[[treatment]], treatment)

  rhs <- c(treatment, covariates)
  fml <- stats::as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))
  w   <- if (is.null(weights)) NULL else data[[weights]]

  fit <- stats::lm(fml, data = data, weights = w)
  vc  <- sandwich::vcovCL(fit, cluster = data[[cluster]], type = "HC1")
  ct  <- lmtest::coeftest(fit, vcov. = vc)

  p_col <- grep("^Pr", colnames(ct), value = TRUE)[1L]

  out <- list(
    call       = match.call(),
    model      = fit,
    vcov       = vc,
    coef_table = ct,
    estimate   = unname(ct[treatment, "Estimate"]),
    se         = unname(ct[treatment, "Std. Error"]),
    p_value    = unname(ct[treatment, p_col])
  )
  class(out) <- "estimate_att_ml"
  out
}

#' @export
print.estimate_att_ml <- function(x, ...) {
  cat("MLCausal treatment effect estimate\n")
  cat("  Estimate: ", round(x$estimate, 4), "\n", sep = "")
  cat("  SE:       ", round(x$se,       4), "\n", sep = "")
  cat("  p-value:  ", signif(x$p_value, 4), "\n", sep = "")
  invisible(x)
}
