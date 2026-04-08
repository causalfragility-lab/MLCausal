#' Within-cluster nearest-neighbour matching with dual-balance optimisation
#'
#' Performs greedy nearest-neighbour matching within clusters using a composite
#' distance that simultaneously targets balance on individual-level propensity
#' scores \emph{and} cluster-mean covariates. This dual-balance approach is the
#' core methodological innovation of \pkg{MLCausal}: standard within-cluster
#' PS matching achieves individual-level balance but often leaves substantial
#' between-cluster (cluster-mean) imbalance. The composite distance penalises
#' matches that worsen cluster-mean balance.
#'
#' @section Dual-balance distance:
#' For each treated unit \eqn{i} and candidate control \eqn{j} in the same
#' cluster \eqn{g}, the composite matching distance is
#'
#' \deqn{d_{ij} = |\mathrm{logit}(\hat{e}_i) - \mathrm{logit}(\hat{e}_j)| +
#'   \lambda \sum_k \frac{|\bar{x}_{k,g}^{(1)} - \bar{x}_{k,g}^{(0)}|}{s_k}}
#'
#' where \eqn{\hat{e}} is the propensity score, \eqn{\bar{x}_{k,g}^{(t)}} is
#' the running cluster-\eqn{g} mean of covariate \eqn{k} for treatment arm
#' \eqn{t} after each match is tentatively accepted, and \eqn{s_k} is the
#' pooled SD of covariate \eqn{k} in the full sample. The tuning parameter
#' \eqn{\lambda} (\code{lambda}) controls the weight given to cluster-mean
#' balance relative to individual PS proximity. Setting \code{lambda = 0}
#' recovers standard within-cluster PS matching.
#'
#' @param ps_fit An object of class \code{ml_ps} from \code{\link{ml_ps}}.
#' @param ratio Positive integer. Controls matched per treated unit. Default 1.
#' @param replace Logical. Allow control reuse? Default \code{FALSE}.
#' @param caliper Optional positive numeric. Maximum allowable logit-PS
#'   distance. \code{NULL} (default) applies no caliper.
#' @param lambda Non-negative numeric. Weight on the cluster-mean balance
#'   penalty in the composite distance. \code{lambda = 0} gives standard
#'   PS matching; \code{lambda = 1} (default) gives equal weight to PS
#'   proximity and cluster-mean balance.
#'
#' @return An object of class \code{ml_match}, a named list:
#' \describe{
#'   \item{call}{The matched call.}
#'   \item{data_matched}{\code{data.frame} of matched units with columns
#'     \code{match_weight} (1 for treated, \eqn{1/k} for controls) and
#'     \code{pair_id}.}
#'   \item{pairs}{\code{data.frame} of matched pairs with \code{pair_id},
#'     \code{treated_row}, \code{control_row}, \code{ps_distance}, and
#'     \code{composite_distance}.}
#'   \item{treatment, cluster, ratio, replace, caliper, lambda}{Stored
#'     arguments.}
#'   \item{n_unmatched}{Number of treated units that could not be matched.}
#' }
#'
#' @seealso \code{\link{ml_ps}}, \code{\link{balance_ml}},
#'   \code{\link{estimate_att_ml}}
#'
#' @examples
#' dat     <- simulate_ml_data(n_clusters = 5, cluster_size = 10, seed = 1)
#' ps      <- ml_ps(dat, "z", c("x1", "x2", "x3"), "school_id")
#' matched <- ml_match(ps, ratio = 1, caliper = 0.5, lambda = 1)
#' print(matched)
#'
#' @export
ml_match <- function(ps_fit,
                     ratio   = 1L,
                     replace = FALSE,
                     caliper = NULL,
                     lambda  = 1) {
  # --- input checks ----------------------------------------------------------
  if (!inherits(ps_fit, "ml_ps"))
    stop("'ps_fit' must be an object of class 'ml_ps' from ml_ps().",
         call. = FALSE)
  if (!is.numeric(ratio) || ratio < 1L)
    stop("'ratio' must be a positive integer.", call. = FALSE)
  if (!is.logical(replace) || length(replace) != 1L)
    stop("'replace' must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(caliper) && (!is.numeric(caliper) || caliper <= 0))
    stop("'caliper' must be a positive number or NULL.", call. = FALSE)
  if (!is.numeric(lambda) || length(lambda) != 1L || lambda < 0)
    stop("'lambda' must be a non-negative number.", call. = FALSE)

  d         <- ps_fit$data
  treatment <- ps_fit$treatment
  cl        <- ps_fit$cluster
  covs      <- ps_fit$covariates
  d$.id     <- seq_len(nrow(d))
  d$.lps    <- .logit(ps_fit$ps)   # match on logit-PS scale

  # Pre-compute pooled SDs for covariate normalisation in penalty term
  pooled_sd <- vapply(covs, function(v) {
    sd_v <- sqrt(.weighted_var(d[[v]]))
    if (is.na(sd_v) || sd_v == 0) 1 else sd_v   # avoid division by zero
  }, numeric(1))

  matched_rows <- list()
  pair_rows    <- list()
  pair_id      <- 1L
  n_unmatched  <- 0L

  for (g in unique(d[[cl]])) {
    dg    <- d[d[[cl]] == g, , drop = FALSE]
    t_idx <- which(dg[[treatment]] == 1)
    c_idx <- which(dg[[treatment]] == 0)
    if (length(t_idx) == 0 || length(c_idx) == 0) next

    available_controls <- c_idx

    # Running cluster-arm means for the dual-balance penalty.
    # These are updated after each match is accepted.
    cm_treated <- vapply(covs,
                         function(v) mean(dg[[v]][t_idx], na.rm = TRUE),
                         numeric(1))
    cm_control <- vapply(covs,
                         function(v) mean(dg[[v]][c_idx], na.rm = TRUE),
                         numeric(1))

    for (ti in t_idx) {
      if (length(available_controls) == 0) {
        n_unmatched <- n_unmatched + 1L
        next
      }

      # -- PS component of composite distance (logit scale) ------------------
      ps_dist <- abs(dg$.lps[available_controls] - dg$.lps[ti])

      # -- Cluster-mean balance penalty --------------------------------------
      # For each candidate control j, compute how accepting j would change
      # the running cluster-mean difference across covariates.
      cm_penalty <- if (lambda == 0 || length(available_controls) == 0) {
        rep(0, length(available_controls))
      } else {
        n_ctrl_so_far <- length(c_idx) - length(available_controls)
        vapply(seq_along(available_controls), function(ji) {
          j   <- available_controls[ji]
          # Tentative new control cluster mean if j is added
          new_cm_ctrl <- (cm_control * (n_ctrl_so_far + length(c_idx)) +
                          dg[j, covs, drop = FALSE]) /
                         (n_ctrl_so_far + length(c_idx) + 1)
          sum(abs(cm_treated - as.numeric(new_cm_ctrl)) / pooled_sd)
        }, numeric(1))
      }

      composite_dist <- ps_dist + lambda * cm_penalty

      # -- Sort by composite distance ----------------------------------------
      sort_ord      <- order(composite_dist)
      cand_sorted   <- available_controls[sort_ord]
      ps_sorted     <- ps_dist[sort_ord]
      comp_sorted   <- composite_dist[sort_ord]

      # -- Apply caliper on logit-PS distance --------------------------------
      if (!is.null(caliper)) {
        keep        <- ps_sorted <= caliper
        cand_sorted <- cand_sorted[keep]
        ps_sorted   <- ps_sorted[keep]
        comp_sorted <- comp_sorted[keep]
      }

      if (length(cand_sorted) == 0) {
        n_unmatched <- n_unmatched + 1L
        next
      }

      n_chosen    <- min(as.integer(ratio), length(cand_sorted))
      chosen      <- cand_sorted[seq_len(n_chosen)]
      chosen_ps   <- ps_sorted[seq_len(n_chosen)]
      chosen_comp <- comp_sorted[seq_len(n_chosen)]

      # -- Update running control cluster mean --------------------------------
      n_ctrl_so_far <- length(c_idx) - length(available_controls)
      for (jj in chosen) {
        cm_control <- (cm_control * (n_ctrl_so_far + length(c_idx)) +
                       as.numeric(dg[jj, covs, drop = FALSE])) /
                      (n_ctrl_so_far + length(c_idx) + 1)
        n_ctrl_so_far <- n_ctrl_so_far + 1
      }

      pair_rows[[length(pair_rows) + 1L]] <- data.frame(
        pair_id            = pair_id,
        treated_row        = dg$.id[ti],
        control_row        = dg$.id[chosen],
        ps_distance        = chosen_ps,
        composite_distance = chosen_comp
      )

      ids               <- c(dg$.id[ti], dg$.id[chosen])
      block             <- d[d$.id %in% ids, , drop = FALSE]
      block$match_weight <- ifelse(block[[treatment]] == 1, 1, 1 / n_chosen)
      block$pair_id     <- pair_id
      matched_rows[[length(matched_rows) + 1L]] <- block

      if (!replace) available_controls <- setdiff(available_controls, chosen)
      pair_id <- pair_id + 1L
    }
  }

  if (length(matched_rows) == 0)
    stop("No matches found. Consider widening the caliper, reducing lambda, ",
         "or using ml_weight() instead.", call. = FALSE)

  if (n_unmatched > 0)
    message(n_unmatched, " treated unit(s) could not be matched and were ",
            "excluded from the matched sample.")

  data_matched <- do.call(rbind, matched_rows)
  # Drop internal helper columns before returning
  data_matched$.id  <- NULL
  data_matched$.lps <- NULL

  pairs <- do.call(rbind, pair_rows)

  out <- list(
    call         = match.call(),
    data_matched = data_matched,
    pairs        = pairs,
    treatment    = treatment,
    cluster      = cl,
    ratio        = as.integer(ratio),
    replace      = replace,
    caliper      = caliper,
    lambda       = lambda,
    n_unmatched  = n_unmatched
  )
  class(out) <- "ml_match"
  out
}

#' @export
print.ml_match <- function(x, ...) {
  cat("MLCausal matched sample (dual-balance)\n")
  cat("  Matched rows:    ", nrow(x$data_matched), "\n", sep = "")
  cat("  Matched sets:    ", length(unique(x$data_matched$pair_id)), "\n",
      sep = "")
  cat("  Clusters used:   ", length(unique(x$data_matched[[x$cluster]])),
      "\n", sep = "")
  cat("  Unmatched (trt): ", x$n_unmatched, "\n", sep = "")
  cat("  Lambda (balance weight): ", x$lambda, "\n", sep = "")
  if (!is.null(x$caliper))
    cat("  Caliper (logit-PS):      ", x$caliper, "\n", sep = "")
  invisible(x)
}
