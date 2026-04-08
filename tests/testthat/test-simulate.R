test_that("simulate_ml_data returns correct structure", {
  dat <- simulate_ml_data(n_clusters = 5, cluster_size = 10, seed = 1)
  expect_s3_class(dat, "data.frame")
  expect_named(dat, c("school_id", "x1", "x2", "x3", "z", "y"))
  expect_true(all(dat$z %in% c(0, 1)))
  expect_equal(length(unique(dat$school_id)), 5)
})

test_that("n_min is respected", {
  dat   <- simulate_ml_data(n_clusters = 5, cluster_size = 5, n_min = 8, seed = 1)
  sizes <- table(dat$school_id)
  expect_true(all(sizes >= 8))
})

test_that("simulate_ml_data validates inputs", {
  expect_error(simulate_ml_data(n_clusters = 1))
  expect_error(simulate_ml_data(cluster_size = 0))
  expect_error(simulate_ml_data(n_min = 0))
})
