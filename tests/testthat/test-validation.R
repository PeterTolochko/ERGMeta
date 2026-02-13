test_that("create_mvbf_formula validates estimate/ses pairs", {
  skip_if_not_installed("brms")

  df <- data.frame(
    estimate_edges = c(0.1, 0.2),
    ses_edges = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )

  f <- create_mvbf_formula(df)
  expect_s3_class(f, "mvbrmsformula")

  bad_df <- data.frame(estimate_edges = c(0.1, 0.2), stringsAsFactors = FALSE)
  expect_error(
    create_mvbf_formula(bad_df),
    "Missing corresponding standard-error columns"
  )
})

test_that("default meta priors match available brms classes", {
  skip_if_not_installed("brms")

  df <- data.frame(
    estimate_edges = c(0.1, 0.2, 0.3),
    ses_edges = c(0.01, 0.02, 0.03),
    estimate_mutual = c(0.2, 0.1, 0.4),
    ses_mutual = c(0.02, 0.03, 0.04),
    stringsAsFactors = FALSE
  )

  current_formula <- create_mvbf_formula(df)
  pri <- .default_meta_priors(current_formula, df)
  expect_s3_class(pri, "brmsprior")
})

test_that("meta_fit fails fast on duplicate network/coef rows", {
  df <- tibble::tibble(
    network = c(1, 1),
    coefs = c("edges", "edges"),
    estimate = c(-1.0, -1.1),
    ses = c(0.1, 0.1)
  )

  expect_error(
    meta_fit(df),
    "duplicate rows for the same network/coefficient combination"
  )
})

test_that("meta_fit fails fast on incomplete network/coef grid", {
  df <- tibble::tibble(
    network = c(1, 1, 2),
    coefs = c("edges", "mutual", "edges"),
    estimate = c(-1.0, -0.4, -0.8),
    ses = c(0.1, 0.2, 0.1)
  )

  expect_error(
    meta_fit(df),
    "Missing estimate/ses values after reshaping"
  )
})

test_that("sequential_estimation rejects empty network list", {
  expect_error(
    sequential_estimation(list(), ~edges),
    "non-empty list of network objects"
  )
})
