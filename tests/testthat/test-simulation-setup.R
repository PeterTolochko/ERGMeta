test_that("simulate_ergmeta_data returns expected structure", {
  out <- simulate_ergmeta_data(
    n_networks = 12,
    coef_terms = c("edges", "mutual"),
    seed = 123
  )

  expect_true(is.list(out))
  expect_true(all(c("estimates", "covariance", "truth", "moderators", "diagnostics", "parameters") %in% names(out)))
  expect_true(all(c("network", "coefs", "term", "estimate", "ses") %in% names(out$estimates)))
  expect_true(all(c("network", "coefs_row", "coefs_col", "covariance") %in% names(out$covariance)))
  expect_true(all(c("network", "coefs", "term", "true_theta") %in% names(out$truth)))
  expect_true(nrow(out$estimates) > 0)
})

test_that("simulate_ergmeta_data is reproducible with fixed seed", {
  out1 <- simulate_ergmeta_data(n_networks = 10, seed = 77)
  out2 <- simulate_ergmeta_data(n_networks = 10, seed = 77)

  expect_identical(out1$estimates, out2$estimates)
  expect_identical(out1$covariance, out2$covariance)
  expect_identical(out1$truth, out2$truth)
})

test_that("run_ergmeta_simulation returns summary for selected methods", {
  bench <- run_ergmeta_simulation(
    n_sims = 4,
    n_networks = 16,
    methods = c("sens_random", "cov_random", "cov_random_mod"),
    progress = FALSE,
    seed = 202
  )

  expect_true(is.list(bench))
  expect_true(all(c("results", "summary", "settings") %in% names(bench)))
  expect_true(all(c("method", "coefs", "bias", "rmse", "coverage_95", "mean_tau2") %in% names(bench$summary)))
  expect_true(all(c("sens_random", "cov_random", "cov_random_mod") %in% bench$summary$method))
})
