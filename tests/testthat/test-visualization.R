test_that("plot_sensitivity returns a ggplot object", {
  skip_if_not_installed("ggplot2")

  sim <- simulate_ergmeta_data(n_networks = 12, seed = 123)
  sens <- run_exclusion_sensitivity(
    df = sim$estimates,
    method = "random",
    include_leave_one_out = TRUE
  )

  p <- plot_sensitivity(sens, value = "delta_from_all")
  expect_s3_class(p, "ggplot")

  roadmap_like <- list(sensitivity = sens)
  p2 <- plot_sensitivity(roadmap_like, value = "pooled_estimate")
  expect_s3_class(p2, "ggplot")
})

test_that("plot_tau supports sensitivity and covariance-fit inputs", {
  skip_if_not_installed("ggplot2")

  sim <- simulate_ergmeta_data(n_networks = 14, seed = 456)
  sens <- run_exclusion_sensitivity(
    df = sim$estimates,
    method = "random",
    include_leave_one_out = FALSE
  )
  p_sens <- plot_tau(sens, scale = "tau2")
  expect_s3_class(p_sens, "ggplot")

  fit <- meta_fit_covariance(
    estimates_df = sim$estimates,
    covariance_df = sim$covariance,
    method = "random"
  )
  p_cov <- plot_tau(fit, scale = "tau")
  expect_s3_class(p_cov, "ggplot")

  roadmap_like <- list(covariance_meta = fit)
  p_cov2 <- plot_tau(roadmap_like, scale = "tau2")
  expect_s3_class(p_cov2, "ggplot")
})

test_that("plot_simulation_summary returns a ggplot object", {
  skip_if_not_installed("ggplot2")

  bench <- run_ergmeta_simulation(
    n_sims = 3,
    n_networks = 16,
    methods = c("sens_random", "cov_random"),
    progress = FALSE,
    seed = 789
  )

  p <- plot_simulation_summary(bench, metric = "rmse")
  expect_s3_class(p, "ggplot")
})

test_that("plot helpers validate required columns", {
  skip_if_not_installed("ggplot2")

  expect_error(
    plot_sensitivity(data.frame(scenario = "all"), value = "tau2"),
    "missing required columns"
  )

  expect_error(
    plot_tau(data.frame(coefs = "edges", term = "edges"), scale = "tau2"),
    "missing required columns"
  )

  expect_error(
    plot_simulation_summary(data.frame(method = "m1"), metric = "bias"),
    "missing required columns"
  )
})
