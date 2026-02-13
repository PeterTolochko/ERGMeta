test_that("prepare_meta_inputs keeps diagnostics and original network ids", {
  skip_if_not_installed("ergm")
  skip_if_not_installed("network")

  set.seed(101)

  make_fit <- function(p, study_id) {
    net <- network::network(
      matrix(rbinom(100, 1, p), nrow = 10),
      directed = TRUE
    )
    network::set.network.attribute(net, "study", study_id)
    ergm::ergm(
      net ~ edges,
      control = ergm::control.ergm(
        MCMLE.maxit = 1,
        MCMC.interval = 100,
        MCMC.samplesize = 1000
      )
    )
  }

  fit_1 <- make_fit(0.20, "A")
  fit_3 <- make_fit(0.30, "B")
  out <- prepare_meta_inputs(
    fit_list = list(fit_1, NULL, fit_3),
    grouping_variable = "study"
  )

  expect_identical(out$included_networks, c(1L, 3L))
  expect_identical(out$excluded_networks, 2L)
  expect_identical(sort(unique(out$estimates$network)), c(1L, 3L))
  expect_identical(sort(unique(out$covariance$network)), c(1L, 3L))
  expect_identical(
    out$diagnostics$included_in_meta,
    c(TRUE, FALSE, TRUE)
  )
})

test_that("extract_within_covariance returns tidy covariance output", {
  skip_if_not_installed("ergm")
  skip_if_not_installed("network")

  set.seed(202)

  make_fit <- function(p, study_id) {
    net <- network::network(
      matrix(rbinom(100, 1, p), nrow = 10),
      directed = TRUE
    )
    network::set.network.attribute(net, "study", study_id)
    ergm::ergm(
      net ~ edges,
      control = ergm::control.ergm(
        MCMLE.maxit = 1,
        MCMC.interval = 100,
        MCMC.samplesize = 1000
      )
    )
  }

  cov_df <- extract_within_covariance(
    fit_list = list(make_fit(0.20, "A"), make_fit(0.30, "B")),
    grouping_variable = "study"
  )

  expect_true(all(
    c(
      "network", "coefs_row", "coefs_col", "term_row", "term_col",
      "covariance", "is_observed", "group_var"
    ) %in% names(cov_df)
  ))
  expect_true(all(cov_df$is_observed))
  expect_identical(sort(unique(cov_df$group_var)), c("A", "B"))
})

test_that("run_exclusion_sensitivity returns baseline and leave-one-out scenarios", {
  skip_if_not_installed("ergm")
  skip_if_not_installed("network")

  set.seed(303)

  make_fit <- function(p) {
    net <- network::network(
      matrix(rbinom(100, 1, p), nrow = 10),
      directed = TRUE
    )
    ergm::ergm(
      net ~ edges,
      control = ergm::control.ergm(
        MCMLE.maxit = 1,
        MCMC.interval = 100,
        MCMC.samplesize = 1000
      )
    )
  }

  est_df <- extract_estimates(list(make_fit(0.20), make_fit(0.30)))
  sens <- run_exclusion_sensitivity(
    df = est_df,
    exclude_sets = list(drop_first = 1L),
    include_leave_one_out = TRUE
  )

  expect_true("all" %in% sens$scenario)
  expect_true("leave_out_1" %in% sens$scenario)
  expect_true("leave_out_2" %in% sens$scenario)
  expect_true("drop_first" %in% sens$scenario)

  baseline_rows <- sens[sens$scenario == "all", , drop = FALSE]
  expect_true(all(abs(baseline_rows$delta_from_all) < 1e-12))
  expect_true(all(!is.na(baseline_rows$tau2)))
  expect_true(all(abs(baseline_rows$delta_tau2_from_all) < 1e-12))
})

test_that("meta_fit_covariance supports moderators and returns tau estimates", {
  skip_if_not_installed("ergm")
  skip_if_not_installed("network")

  set.seed(404)

  make_fit <- function(p, density_mod) {
    net <- network::network(
      matrix(rbinom(100, 1, p), nrow = 10),
      directed = TRUE
    )
    network::set.network.attribute(net, "density_mod", density_mod)
    ergm::ergm(
      net ~ edges,
      control = ergm::control.ergm(
        MCMLE.maxit = 1,
        MCMC.interval = 100,
        MCMC.samplesize = 1000
      )
    )
  }

  fit_a <- make_fit(0.20, 0.1)
  fit_b <- make_fit(0.25, 0.2)
  fit_c <- make_fit(0.30, 0.3)

  estimates <- extract_estimates(list(fit_a, fit_b, fit_c), grouping_variable = "density_mod")
  cov_df <- extract_within_covariance(list(fit_a, fit_b, fit_c), grouping_variable = "density_mod")

  mod_df <- unique(estimates[, c("network", "group_var"), drop = FALSE])
  names(mod_df)[2] <- "density_mod"
  mod_df$density_mod <- as.numeric(mod_df$density_mod)

  fit <- meta_fit_covariance(
    estimates_df = estimates,
    covariance_df = cov_df,
    method = "random",
    moderators = "density_mod",
    moderator_df = mod_df
  )

  expect_true(is.list(fit))
  expect_true(all(c("coefficients", "tau", "reference_summary", "fitted") %in% names(fit)))
  expect_true("density_mod" %in% fit$coefficients$predictor)
  expect_true(all(fit$tau$tau2 >= 0))
})

test_that("run_exclusion_sensitivity supports moderator meta-regression", {
  skip_if_not_installed("ergm")
  skip_if_not_installed("network")

  set.seed(505)

  make_fit <- function(p, density_mod) {
    net <- network::network(
      matrix(rbinom(100, 1, p), nrow = 10),
      directed = TRUE
    )
    network::set.network.attribute(net, "density_mod", density_mod)
    ergm::ergm(
      net ~ edges,
      control = ergm::control.ergm(
        MCMLE.maxit = 1,
        MCMC.interval = 100,
        MCMC.samplesize = 1000
      )
    )
  }

  fit_a <- make_fit(0.20, 0.1)
  fit_b <- make_fit(0.25, 0.2)
  fit_c <- make_fit(0.30, 0.3)

  estimates <- extract_estimates(list(fit_a, fit_b, fit_c), grouping_variable = "density_mod")
  mod_df <- unique(estimates[, c("network", "group_var"), drop = FALSE])
  names(mod_df)[2] <- "density_mod"
  mod_df$density_mod <- as.numeric(mod_df$density_mod)

  sens <- run_exclusion_sensitivity(
    df = estimates,
    method = "random",
    moderators = "density_mod",
    moderator_df = mod_df,
    include_leave_one_out = TRUE
  )

  expect_true("density_mod=0.2" %in% sens$reference_profile)
  expect_true(all(is.finite(sens$tau2[!is.na(sens$tau2)])))
})
