#' Multilevel covariate balance diagnostics
#'
#' Computes standardised mean differences (SMDs) at two levels simultaneously:
#' the individual level and the cluster-mean level. Cluster-mean SMDs capture
#' between-cluster imbalance that individual-level diagnostics miss in
#' hierarchical data, and are the key diagnostic for assessing whether
#' \code{\link{ml_match}}'s dual-balance optimisation has succeeded.
#'
#' @details
#' A common rule of thumb is \eqn{|\text{SMD}| < 0.10} for adequate balance.
#' Cluster-mean SMDs above this threshold indicate residual between-cluster
#' confounding; in that case increase \code{lambda} in \code{\link{ml_match}}
#' or switch to \code{method = "fixed"} in \code{\link{ml_ps}}.
#'
#' When a cluster-mean SMD cannot be computed (e.g. because the matched sample
#' contains only one treatment arm across clusters), the \code{smd} column
#' returns the string
#' \code{"Cluster-level SMD not estimable due to insufficient within-cluster variation after matching."}
#' rather than \code{NA}. This makes the diagnostic failure explicit and
#' actionable rather than silently missing.
#'
#' @param data A \code{data.frame}.
#' @param treatment Character string. Name of the binary treatment variable
#'   (\code{0}/\code{1}).
#' @param covariates Character vector of covariate names.
#' @param cluster Character string. Name of the cluster identifier variable.
#' @param weights Optional character string. Name of a weight variable in
#'   \code{data} (e.g. \code{"weights"} from \code{\link{ml_weight}} or
#'   \code{"match_weight"} from \code{\link{ml_match}}). \code{NULL} (default)
#'   computes unweighted SMDs.
#'
#' @return An object of class \code{balance_ml}, a named list:
#' \describe{
#'   \item{overall}{Data frame with individual-level SMDs, one row per
#'     covariate. The \code{smd} column is numeric.}
#'   \item{cluster_means}{Data frame with cluster-mean SMDs, one row per
#'     covariate. The \code{smd} column is character: a numeric value
#'     formatted to 4 decimal places, or a descriptive message when not
#'     estimable.}
#'   \item{summary}{Row-bound combination of \code{overall} and
#'     \code{cluster_means}.}
#' }
#'
#' @seealso \code{\link{ml_match}}, \code{\link{ml_weight}},
#'   \code{\link{ml_ps}}
#'
#' @examples
#' dat   <- simulate_ml_data(n_clusters = 5, cluster_size = 10, seed = 1)
#' ps    <- ml_ps(dat, "z", c("x1", "x2", "x3"), "school_id")
#' dat_w <- ml_weight(ps)
#' bal   <- balance_ml(dat_w, treatment = "z",
#'                     covariates = c("x1", "x2", "x3"),
#'                     cluster = "school_id", weights = "weights")
#' print(bal)
#'
#' @importFrom stats aggregate setNames
#' @export
balance_ml <- function(data,
                       treatment,
                       covariates,
                       cluster,
                       weights = NULL) {
  # --- input validation ------------------------------------------------------
  if (!is.data.frame(data))
    stop("'data' must be a data.frame.", call. = FALSE)
  if (!is.character(treatment) || length(treatment) != 1L)
    stop("'treatment' must be a single character string.", call. = FALSE)
  if (!is.character(covariates) || length(covariates) == 0L)
    stop("'covariates' must be a non-empty character vector.", call. = FALSE)
  if (!is.character(cluster) || length(cluster) != 1L)
    stop("'cluster' must be a single character string.", call. = FALSE)
  if (!is.null(weights) &&
      (!is.character(weights) || length(weights) != 1L))
    stop("'weights' must be a single character string or NULL.", call. = FALSE)

  .check_required_columns(data, c(treatment, covariates, cluster))
  if (!is.null(weights)) .check_required_columns(data, weights)

  z <- data[[treatment]]
  w <- if (is.null(weights)) NULL else data[[weights]]

  # --- Individual-level SMDs (numeric) ---------------------------------------
  overall <- do.call(rbind, lapply(covariates, function(v) {
    data.frame(
      variable = v,
      level    = "individual",
      smd      = .compute_smd(data[[v]], z, w, label = v),
      stringsAsFactors = FALSE
    )
  }))

  # --- Cluster-mean SMDs (character, never silent NA) ------------------------
  # aggregate() gives exactly one row per cluster, avoiding the
  # ave()-length-mismatch bug with unequal cluster sizes.
  agg_vars   <- c(treatment, covariates)
  cluster_df <- stats::aggregate(
    stats::setNames(data[agg_vars], agg_vars),
    by  = stats::setNames(list(data[[cluster]]), cluster),
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  zc <- as.integer(cluster_df[[treatment]] >= 0.5)

  not_estimable_msg <- paste0(
    "Cluster-level SMD not estimable due to insufficient ",
    "within-cluster variation after matching."
  )

  cluster_balance <- do.call(rbind, lapply(covariates, function(v) {
    smd_val <- tryCatch({
      val <- withCallingHandlers(
        .compute_smd(cluster_df[[v]], zc, NULL,
                     label = paste0(v, " (cluster mean)")),
        warning = function(w) {
          invokeRestart("muffleWarning")
        }
      )
      # Check for NA after muffling internal warning
      if (is.na(val)) not_estimable_msg else formatC(val, digits = 4,
                                                      format = "f")
    }, error = function(e) {
      not_estimable_msg
    })

    data.frame(
      variable = v,
      level    = "cluster_mean",
      smd      = smd_val,
      stringsAsFactors = FALSE
    )
  }))

  # Combine: overall$smd is numeric, cluster_balance$smd is character.
  # Convert overall$smd to character for a uniform summary table.
  overall_chr        <- overall
  overall_chr$smd    <- formatC(overall$smd, digits = 4, format = "f")
  summary_df         <- rbind(overall_chr, cluster_balance)

  out <- list(
    overall       = overall,          # smd column is numeric
    cluster_means = cluster_balance,  # smd column is character
    summary       = summary_df        # smd column is character (uniform)
  )
  class(out) <- "balance_ml"
  out
}

#' @export
print.balance_ml <- function(x, ...) {
  cat("MLCausal balance diagnostics\n")
  cat("\nIndividual-level SMDs:\n")
  tmp <- x$overall
  tmp$smd <- round(tmp$smd, 4)
  print(tmp, row.names = FALSE)
  cat("\nCluster-mean SMDs:\n")
  print(x$cluster_means, row.names = FALSE)
  invisible(x)
}
