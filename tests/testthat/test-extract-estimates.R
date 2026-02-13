test_that("extract_estimates keeps term alignment and network-level groups", {
  skip_if_not_installed("ergm")
  skip_if_not_installed("network")

  set.seed(123)

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

  fit_a <- make_fit(0.20, "A")
  fit_b <- make_fit(0.30, "B")

  out <- extract_estimates(
    fit_list = list(fit_a, fit_b),
    grouping_variable = "study"
  )

  expect_true(all(c("network", "coefs", "term", "estimate", "ses", "group_var") %in% names(out)))
  expect_true(all(out$term == "edges"))
  expect_identical(sort(unique(out$coefs)), "edges")
  expect_identical(out$group_var, c("A", "B"))
})
