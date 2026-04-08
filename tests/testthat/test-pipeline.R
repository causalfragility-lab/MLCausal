dat <- simulate_ml_data(n_clusters = 8, cluster_size = 12, n_min = 8, seed = 42)
ps  <- ml_ps(dat, treatment = "z", covariates = c("x1", "x2", "x3"),
             cluster = "school_id")

# ── ml_weight ──────────────────────────────────────────────────────────────
test_that("ml_weight returns data frame with weights column", {
  dat_w <- ml_weight(ps)
  expect_true("weights" %in% names(dat_w))   # Fix 3: column is 'weights'
  expect_true(all(dat_w$weights > 0))
  expect_error(ml_weight(list()))
  expect_error(ml_weight(ps, trim = -1))
})

# ── balance_ml ─────────────────────────────────────────────────────────────
test_that("balance_ml returns balance_ml object with correct structure", {
  dat_w <- ml_weight(ps)
  bal   <- balance_ml(dat_w, treatment = "z",    # Fix 3: treatment
                      covariates = c("x1", "x2", "x3"),
                      cluster = "school_id", weights = "weights")
  expect_s3_class(bal, "balance_ml")
  expect_named(bal, c("overall", "cluster_means", "summary"))
  expect_equal(nrow(bal$overall), 3L)
  expect_true(is.numeric(bal$overall$smd))       # Fix 1: overall$smd numeric
  expect_true(is.character(bal$cluster_means$smd)) # Fix 1: cluster smd character
})

test_that("balance_ml cluster-mean SMD is descriptive string when not estimable", {
  # Create a degenerate case: all treated in cluster 1, all control in cluster 2
  d2 <- data.frame(school_id = c(rep(1, 5), rep(2, 5)),
                   x1 = rnorm(10), z = c(rep(1, 5), rep(0, 5)),
                   stringsAsFactors = FALSE)
  # Should not return NA silently; should return descriptive message
  bal2 <- balance_ml(d2, treatment = "z", covariates = "x1",
                     cluster = "school_id")
  expect_true(is.character(bal2$cluster_means$smd))
  expect_true(grepl("not estimable", bal2$cluster_means$smd))  # Fix 1
})

# ── estimate_att_ml ────────────────────────────────────────────────────────
test_that("estimate_att_ml returns correct structure", {
  dat_w <- ml_weight(ps)
  est   <- estimate_att_ml(dat_w, outcome = "y", treatment = "z",  # Fix 3
                           cluster = "school_id", weights = "weights")
  expect_s3_class(est, "estimate_att_ml")
  expect_true(is.numeric(est$estimate))
  expect_true(est$se > 0)
  expect_true(est$p_value >= 0 & est$p_value <= 1)
})

test_that("estimate_att_ml validates inputs", {
  expect_error(estimate_att_ml(list(), "y", "z", "school_id"))
  expect_error(estimate_att_ml(dat, "nonexistent", "z", "school_id"))
})

# ── ml_match (dual-balance) ────────────────────────────────────────────────
test_that("ml_match returns ml_match object with lambda stored", {
  matched <- ml_match(ps, ratio = 1, caliper = 0.5, lambda = 1)
  expect_s3_class(matched, "ml_match")
  expect_true("match_weight" %in% names(matched$data_matched))
  expect_true("pair_id"      %in% names(matched$data_matched))
  expect_equal(matched$lambda, 1)
  expect_true("ps_distance"        %in% names(matched$pairs))
  expect_true("composite_distance" %in% names(matched$pairs))
})

test_that("ml_match lambda=0 recovers standard PS matching", {
  m0 <- ml_match(ps, ratio = 1, caliper = 0.5, lambda = 0)
  m1 <- ml_match(ps, ratio = 1, caliper = 0.5, lambda = 1)
  # With lambda=0, composite distance == ps_distance
  expect_equal(m0$pairs$ps_distance, m0$pairs$composite_distance)
  # Different lambda can produce different matched sets
  expect_s3_class(m0, "ml_match")
  expect_s3_class(m1, "ml_match")
})

test_that("ml_match validates inputs", {
  expect_error(ml_match(list()))
  expect_error(ml_match(ps, ratio = 0))
  expect_error(ml_match(ps, caliper = -1))
  expect_error(ml_match(ps, lambda = -0.1))
})

# ── sens_ml ────────────────────────────────────────────────────────────────
test_that("sens_ml returns correct structure", {
  dat_w <- ml_weight(ps)
  est   <- estimate_att_ml(dat_w, "y", "z", "school_id", weights = "weights")
  sens  <- sens_ml(est$estimate, est$se)
  expect_named(sens, c("confounder_strength", "adjusted_estimate",
                        "original_z", "adjusted_z", "crosses_null"))
  expect_true(is.logical(sens$crosses_null))
  expect_error(sens_ml("a", 0.1))
  expect_error(sens_ml(0.5, -1))
})
