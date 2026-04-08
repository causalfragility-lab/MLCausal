dat <- simulate_ml_data(n_clusters = 5, cluster_size = 10, seed = 1)

test_that("ml_ps returns ml_ps object with all methods", {
  for (m in c("mundlak", "fixed", "single")) {
    ps <- ml_ps(dat, treatment = "z", covariates = c("x1", "x2"),
                cluster = "school_id", method = m)
    expect_s3_class(ps, "ml_ps")
    expect_true(all(ps$ps > 0 & ps$ps < 1))
    expect_equal(length(ps$ps), nrow(dat))
    expect_equal(ps$treatment, "z")   # Fix 3: field is 'treatment'
  }
})

test_that("ml_ps validates inputs", {
  expect_error(ml_ps(list(), "z", "x1", "school_id"))
  expect_error(ml_ps(dat, "nonexistent", "x1", "school_id"))
  dat2      <- dat
  dat2$z    <- 2L
  expect_error(ml_ps(dat2, "z", "x1", "school_id"))
  dat3      <- dat
  dat3$z[1] <- NA
  expect_error(ml_ps(dat3, "z", "x1", "school_id"))
})
