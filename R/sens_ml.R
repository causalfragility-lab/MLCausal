#' Tipping-point sensitivity analysis for omitted cluster confounding
#'
#' Reports how strong a hypothetical omitted confounder (in units of the
#' cluster-robust SE) would need to be to nullify the estimated treatment
#' effect (push \eqn{|z_{\text{adj}}|} below 1.96).
#'
#' @section Statistical model:
#' The analysis assumes a linear bias model: an omitted variable \eqn{U} with
#' confounder strength \eqn{\Gamma} shifts the estimate by
#' \eqn{\delta = \Gamma \times \mathrm{SE}} and the z-statistic by
#' \eqn{\Gamma}. This is transparent and easy to communicate, but is not a
#' full formalisation. For rigorous sensitivity analysis see Cinelli & Hazlett
#' (2020) or Rosenbaum bounds.
#'
#' A message is printed if no tipping point is found within the examined range.
#'
#' @param estimate Numeric scalar. Treatment effect estimate from
#'   \code{\link{estimate_att_ml}}.
#' @param se Numeric scalar. Cluster-robust standard error of the estimate.
#' @param q Numeric vector. Confounder strengths \eqn{\Gamma} to examine.
#'   Default \code{seq(0, 5, by = 0.1)}.
#'
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{confounder_strength}{The \eqn{\Gamma} value examined.}
#'   \item{adjusted_estimate}{Estimate after subtracting assumed bias.}
#'   \item{original_z}{Observed \eqn{|\text{estimate}/\text{SE}|}.}
#'   \item{adjusted_z}{\eqn{\max(0,\, z_{\text{obs}} - \Gamma)}.}
#'   \item{crosses_null}{\code{TRUE} when \eqn{|\text{adjusted\_z}| < 1.96}.}
#' }
#'
#' @references
#' Cinelli, C. & Hazlett, C. (2020). Making sense of sensitivity: Extending
#' omitted variable bias. \emph{JRSS-B}, 82(1), 39--67.
#' \doi{10.1111/rssb.12348}
#'
#' @seealso \code{\link{estimate_att_ml}}
#'
#' @examples
#' dat   <- simulate_ml_data(n_clusters = 5, cluster_size = 10, seed = 1)
#' ps    <- ml_ps(dat, "z", c("x1", "x2", "x3"), "school_id")
#' dat_w <- ml_weight(ps)
#' est   <- estimate_att_ml(dat_w, "y", "z", "school_id",
#'                          weights = "weights")
#' sens  <- sens_ml(est$estimate, est$se)
#' head(sens)
#'
#' @export
sens_ml <- function(estimate,
                    se,
                    q = seq(0, 5, by = 0.1)) {
  if (!is.numeric(estimate) || length(estimate) != 1L || !is.finite(estimate))
    stop("'estimate' must be a single finite numeric value.", call. = FALSE)
  if (!is.numeric(se) || length(se) != 1L || se <= 0 || !is.finite(se))
    stop("'se' must be a single positive finite numeric value.", call. = FALSE)
  if (!is.numeric(q) || any(q < 0))
    stop("'q' must be a numeric vector of non-negative values.", call. = FALSE)

  z_obs <- abs(estimate / se)

  out <- data.frame(
    confounder_strength = q,
    adjusted_estimate   = (abs(estimate) - q * se) * sign(estimate),
    original_z          = z_obs,
    adjusted_z          = pmax(0, z_obs - q),
    crosses_null        = pmax(0, z_obs - q) < 1.96,
    stringsAsFactors    = FALSE
  )

  if (!any(out$crosses_null))
    message("The effect does not cross the null within the examined range ",
            "(max Gamma = ", max(q), "). The result appears robust to ",
            "omitted confounding of this magnitude. Consider extending 'q' ",
            "if a tipping point is still required.")

  out
}
