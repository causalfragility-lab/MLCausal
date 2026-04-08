# Internal utility functions for MLCausal.
# None of these are exported. Do NOT add roxygen2 #' tags here — doing so
# causes devtools::document() to generate .Rd stubs for internal functions,
# which triggers "undocumented arguments" warnings in R CMD check.

# ---------------------------------------------------------------------------
# Input validation helpers
# ---------------------------------------------------------------------------

.check_required_columns <- function(data, cols) {
  miss <- setdiff(cols, names(data))
  if (length(miss) > 0)
    stop("Missing required columns: ", paste(miss, collapse = ", "),
         call. = FALSE)
}

.check_treatment <- function(z, treatment_name) {
  if (anyNA(z))
    stop("Treatment variable '", treatment_name, "' contains missing values. ",
         "Please remove or impute before proceeding.", call. = FALSE)
  vals <- unique(z)
  if (!all(vals %in% c(0, 1)))
    stop("Treatment variable '", treatment_name, "' must be coded 0/1. ",
         "Found values: ", paste(sort(vals), collapse = ", "), call. = FALSE)
  if (length(vals) < 2L)
    stop("Treatment variable '", treatment_name,
         "' has only one unique value. Both treated (1) and control (0) ",
         "units are required.", call. = FALSE)
}

.check_covariates <- function(data, covariates) {
  for (v in covariates) {
    n_miss <- sum(is.na(data[[v]]))
    if (n_miss > 0)
      warning("Covariate '", v, "' contains ", n_miss,
              " missing value(s). Affected rows will be dropped by glm/lm.",
              call. = FALSE)
  }
}

.check_cluster_sizes <- function(data, cluster, min_size = 2L) {
  sizes <- table(data[[cluster]])
  small <- names(sizes)[sizes < min_size]
  if (length(small) > 0)
    warning(length(small), " cluster(s) have fewer than ", min_size,
            " observation(s): ", paste(small, collapse = ", "), ". ",
            "This may cause instability in cluster-level diagnostics.",
            call. = FALSE)
}

# ---------------------------------------------------------------------------
# Statistical helpers
# ---------------------------------------------------------------------------

.weighted_mean <- function(x, w = NULL) {
  if (is.null(w)) return(mean(x, na.rm = TRUE))
  ok <- is.finite(x) & is.finite(w)
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

.weighted_var <- function(x, w = NULL) {
  if (is.null(w)) return(stats::var(x, na.rm = TRUE))
  ok <- is.finite(x) & is.finite(w)
  x <- x[ok]; w <- w[ok]
  if (length(x) < 2L || sum(w) <= 0) return(NA_real_)
  mu <- sum(w * x) / sum(w)
  sum(w * (x - mu)^2) / sum(w)
}

# Standardised mean difference. Returns NA_real_ with an informative warning
# (not silent) when not estimable.
.compute_smd <- function(x, z, w = NULL, label = NULL) {
  z  <- as.integer(z)
  wt <- if (is.null(w)) rep(1, length(z)) else w
  m1 <- .weighted_mean(x[z == 1], wt[z == 1])
  m0 <- .weighted_mean(x[z == 0], wt[z == 0])
  v1 <- .weighted_var(x[z == 1],  wt[z == 1])
  v0 <- .weighted_var(x[z == 0],  wt[z == 0])
  sd_pool <- sqrt((v1 + v0) / 2)
  if (!is.finite(sd_pool) || sd_pool == 0) {
    lbl <- if (!is.null(label)) paste0(" for '", label, "'") else ""
    warning("SMD", lbl, " not estimable: pooled SD is zero or missing. ",
            "This usually means both treatment arms are not represented ",
            "in this subset. Check cluster composition after matching.",
            call. = FALSE)
    return(NA_real_)
  }
  (m1 - m0) / sd_pool
}

# Clip propensity scores strictly away from 0 and 1.
.clip01 <- function(x, eps = 1e-6) pmin(pmax(x, eps), 1 - eps)

# Logit transformation (safe: clips first).
.logit <- function(p) {
  p <- .clip01(p)
  log(p / (1 - p))
}
