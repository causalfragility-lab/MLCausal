#' Multilevel inverse probability weights
#'
#' Computes ATT or ATE inverse probability weights from propensity scores
#' estimated by \code{\link{ml_ps}}. Stabilised weights (the default) are
#' recommended to reduce variance without sacrificing consistency.
#'
#' @section Weight formulae:
#' Let \eqn{\hat{e}_i} be the estimated propensity score and
#' \eqn{\bar{p} = \Pr(Z=1)} the marginal treatment prevalence.
#'
#' \strong{Unstabilised:}
#' \itemize{
#'   \item ATT: \eqn{w=1} (treated), \eqn{w=\hat{e}/(1-\hat{e})} (control)
#'   \item ATE: \eqn{w=1/\hat{e}} (treated), \eqn{w=1/(1-\hat{e})} (control)
#' }
#' \strong{Stabilised:}
#' \itemize{
#'   \item ATT: \eqn{w=1} (treated),
#'     \eqn{w=\bar{p}\hat{e}/[(1-\bar{p})(1-\hat{e})]} (control)
#'   \item ATE: \eqn{w=\bar{p}/\hat{e}} (treated),
#'     \eqn{w=(1-\bar{p})/(1-\hat{e})} (control)
#' }
#'
#' @param ps_fit An object of class \code{ml_ps} from \code{\link{ml_ps}}.
#' @param estimand Character string. \code{"ATT"} (default) or \code{"ATE"}.
#'   Defaults to the estimand stored in \code{ps_fit}.
#' @param stabilize Logical. \code{TRUE} (default) uses stabilised weights.
#' @param trim Optional positive numeric. Weights above this value are
#'   Winsorised. \code{NULL} (default) applies no trimming.
#'
#' @return The working data frame from \code{ps_fit} with an appended
#'   \code{weights} column.
#'
#' @seealso \code{\link{ml_ps}}, \code{\link{balance_ml}},
#'   \code{\link{estimate_att_ml}}
#'
#' @examples
#' dat   <- simulate_ml_data(n_clusters = 5, cluster_size = 10, seed = 1)
#' ps    <- ml_ps(dat, "z", c("x1", "x2", "x3"), "school_id")
#' dat_w <- ml_weight(ps, estimand = "ATT", stabilize = TRUE, trim = 10)
#' summary(dat_w$weights)
#'
#' @export
ml_weight <- function(ps_fit,
                      estimand  = c("ATT", "ATE"),
                      stabilize = TRUE,
                      trim      = NULL) {
  if (!inherits(ps_fit, "ml_ps"))
    stop("'ps_fit' must be an object of class 'ml_ps' from ml_ps().",
         call. = FALSE)
  if (!is.logical(stabilize) || length(stabilize) != 1L)
    stop("'stabilize' must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(trim) && (!is.numeric(trim) || trim <= 0))
    stop("'trim' must be a positive number or NULL.", call. = FALSE)

  estimand <- match.arg(estimand)
  d        <- ps_fit$data
  z        <- d[[ps_fit$treatment]]
  ps       <- .clip01(ps_fit$ps)
  p_treat  <- mean(z, na.rm = TRUE)

  if (stabilize) {
    w <- switch(estimand,
      ATT = ifelse(z == 1, 1,
                   p_treat * ps / ((1 - p_treat) * (1 - ps))),
      ATE = ifelse(z == 1, p_treat / ps,
                   (1 - p_treat) / (1 - ps))
    )
  } else {
    w <- switch(estimand,
      ATT = ifelse(z == 1, 1, ps / (1 - ps)),
      ATE = ifelse(z == 1, 1 / ps, 1 / (1 - ps))
    )
  }

  if (!is.null(trim)) {
    n_trimmed <- sum(w > trim, na.rm = TRUE)
    if (n_trimmed > 0)
      message(n_trimmed, " weight(s) Winsorised to upper bound of ", trim, ".")
    w <- pmin(w, trim)
  }

  d$weights <- w
  d
}
