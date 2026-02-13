#' Simulate ERGMeta-Style Meta-Analysis Inputs
#'
#' @param n_networks Number of networks to simulate.
#' @param coef_terms Character vector of coefficient term names.
#' @param beta_intercept True coefficient means at moderator value 0.
#' @param beta_moderator True moderator slopes for each coefficient.
#' @param tau Between-network standard deviations for each coefficient.
#' @param tau_cor Common correlation used for between-network covariance.
#' @param within_sd Baseline within-network standard deviations for each coefficient.
#' @param within_cor Common correlation used for within-network covariance.
#' @param within_sd_jitter Log-normal jitter scale for network-specific within SDs.
#' @param moderator_name Name of the simulated network-level moderator.
#' @param moderator_mean Mean of the moderator distribution.
#' @param moderator_sd Standard deviation of the moderator distribution.
#' @param exclusion_prob Probability that a network is excluded (to emulate failed estimation).
#' @param seed Optional integer seed.
#'
#' @returns A list with tidy \code{estimates}, tidy \code{covariance},
#'   \code{truth}, \code{moderators}, \code{diagnostics}, and \code{parameters}.
#' @export
simulate_ergmeta_data <- function(n_networks = 40,
                                  coef_terms = c("edges", "mutual"),
                                  beta_intercept = c(-2.2, 0.4),
                                  beta_moderator = c(0.3, -0.2),
                                  tau = c(0.25, 0.15),
                                  tau_cor = 0,
                                  within_sd = c(0.25, 0.20),
                                  within_cor = 0.25,
                                  within_sd_jitter = 0.15,
                                  moderator_name = "moderator_x",
                                  moderator_mean = 0,
                                  moderator_sd = 1,
                                  exclusion_prob = 0,
                                  seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (!is.numeric(n_networks) || length(n_networks) != 1 || n_networks < 2) {
    stop("n_networks must be a single number >= 2.")
  }
  n_networks <- as.integer(n_networks)

  if (!is.character(coef_terms) || length(coef_terms) == 0) {
    stop("coef_terms must be a non-empty character vector.")
  }
  p <- length(coef_terms)
  coef_ids <- make.names(coef_terms)

  .check_len <- function(x, nm) {
    if (!is.numeric(x) || length(x) != p) {
      stop(sprintf("%s must be numeric with length equal to coef_terms.", nm))
    }
  }
  .check_len(beta_intercept, "beta_intercept")
  .check_len(beta_moderator, "beta_moderator")
  .check_len(tau, "tau")
  .check_len(within_sd, "within_sd")

  if (!is.numeric(tau_cor) || length(tau_cor) != 1 || tau_cor <= -0.99 || tau_cor >= 0.99) {
    stop("tau_cor must be a single number in (-0.99, 0.99).")
  }
  if (!is.numeric(within_cor) || length(within_cor) != 1 || within_cor <= -0.99 || within_cor >= 0.99) {
    stop("within_cor must be a single number in (-0.99, 0.99).")
  }
  if (!is.numeric(within_sd_jitter) || length(within_sd_jitter) != 1 || within_sd_jitter < 0) {
    stop("within_sd_jitter must be a single non-negative number.")
  }
  if (!is.character(moderator_name) || length(moderator_name) != 1 || moderator_name == "") {
    stop("moderator_name must be a single non-empty character string.")
  }
  if (!is.numeric(moderator_mean) || !is.numeric(moderator_sd) || moderator_sd <= 0) {
    stop("moderator_mean must be numeric and moderator_sd must be > 0.")
  }
  if (!is.numeric(exclusion_prob) || length(exclusion_prob) != 1 || exclusion_prob < 0 || exclusion_prob >= 1) {
    stop("exclusion_prob must be in [0, 1).")
  }

  .make_cov <- function(sd_vec, cor_value) {
    if (length(sd_vec) == 1) {
      return(matrix(sd_vec^2, nrow = 1, ncol = 1))
    }
    cor_mat <- matrix(cor_value, nrow = length(sd_vec), ncol = length(sd_vec))
    diag(cor_mat) <- 1
    cov_mat <- diag(sd_vec, nrow = length(sd_vec)) %*% cor_mat %*% diag(sd_vec, nrow = length(sd_vec))
    cov_mat
  }

  .draw_mvnorm_zero <- function(cov_mat) {
    l <- tryCatch(chol(cov_mat), error = function(e) NULL)
    if (is.null(l)) {
      l <- tryCatch(chol(cov_mat + diag(1e-10, nrow(cov_mat))), error = function(e) NULL)
    }
    if (is.null(l)) {
      stop("Could not draw from covariance matrix; matrix is not positive definite.")
    }
    as.vector(t(l) %*% rnorm(nrow(cov_mat)))
  }

  between_cov <- .make_cov(tau, tau_cor)
  moderator_values <- rnorm(n_networks, mean = moderator_mean, sd = moderator_sd)
  network_ids <- seq_len(n_networks)

  excluded_networks <- which(runif(n_networks) < exclusion_prob)
  if (length(excluded_networks) >= n_networks) {
    excluded_networks <- excluded_networks[-length(excluded_networks)]
  }
  included_networks <- setdiff(network_ids, excluded_networks)

  estimates_list <- vector("list", n_networks)
  covariance_list <- vector("list", n_networks)
  truth_list <- vector("list", n_networks)

  for (i in network_ids) {
    x_i <- moderator_values[[i]]
    between_draw <- .draw_mvnorm_zero(between_cov)
    theta_true <- beta_intercept + beta_moderator * x_i + between_draw

    se_i <- within_sd * exp(rnorm(p, mean = 0, sd = within_sd_jitter))
    within_cov_i <- .make_cov(se_i, within_cor)
    error_draw <- .draw_mvnorm_zero(within_cov_i)
    theta_obs <- theta_true + error_draw

    truth_list[[i]] <- tibble::tibble(
      network = i,
      coefs = coef_ids,
      term = coef_terms,
      true_theta = theta_true,
      moderator_value = x_i
    )

    if (!(i %in% included_networks)) {
      next
    }

    estimates_list[[i]] <- tibble::tibble(
      network = i,
      coefs = coef_ids,
      term = coef_terms,
      estimate = theta_obs,
      ses = se_i,
      moderator_value = x_i
    )

    idx <- expand.grid(row = seq_len(p), col = seq_len(p))
    covariance_list[[i]] <- tibble::tibble(
      network = i,
      coefs_row = coef_ids[idx$row],
      coefs_col = coef_ids[idx$col],
      term_row = coef_terms[idx$row],
      term_col = coef_terms[idx$col],
      covariance = within_cov_i[cbind(idx$row, idx$col)],
      is_observed = TRUE
    )
  }

  estimates <- dplyr::bind_rows(estimates_list)
  covariance <- dplyr::bind_rows(covariance_list)
  truth <- dplyr::bind_rows(truth_list)

  if (nrow(estimates) == 0) {
    stop("Simulation produced no included networks. Reduce exclusion_prob.")
  }

  moderators <- tibble::tibble(
    network = network_ids,
    moderator_value = moderator_values
  )
  names(moderators)[2] <- moderator_name

  diagnostics <- tibble::tibble(
    network = network_ids,
    included_in_meta = network_ids %in% included_networks,
    status = ifelse(network_ids %in% included_networks, "estimated", "excluded_simulated")
  )

  out <- list(
    estimates = estimates,
    covariance = covariance,
    truth = truth,
    moderators = moderators,
    diagnostics = diagnostics,
    parameters = list(
      n_networks = n_networks,
      coef_terms = coef_terms,
      beta_intercept = beta_intercept,
      beta_moderator = beta_moderator,
      tau = tau,
      tau_cor = tau_cor,
      within_sd = within_sd,
      within_cor = within_cor,
      within_sd_jitter = within_sd_jitter,
      moderator_name = moderator_name,
      moderator_mean = moderator_mean,
      moderator_sd = moderator_sd,
      exclusion_prob = exclusion_prob
    )
  )
  class(out) <- c("ergmeta_simulated_data", class(out))
  out
}


.collect_sensitivity_baseline <- function(sens_df) {
  base <- sens_df[sens_df$scenario == "all", , drop = FALSE]
  if (nrow(base) == 0) {
    return(tibble::tibble(
      coefs = NA_character_,
      term = NA_character_,
      estimate = NA_real_,
      se = NA_real_,
      tau2 = NA_real_,
      converged = FALSE
    ))
  }
  tibble::tibble(
    coefs = base$coefs,
    term = base$term,
    estimate = base$pooled_estimate,
    se = base$pooled_se,
    tau2 = base$tau2,
    converged = base$converged
  )
}


.collect_covariance_fit <- function(fit_obj) {
  base <- fit_obj$reference_summary[, c("coefs", "term", "pooled_estimate", "pooled_se"), drop = FALSE]
  tau <- fit_obj$tau[, c("coefs", "tau2"), drop = FALSE]
  out <- dplyr::left_join(base, tau, by = "coefs")
  tibble::tibble(
    coefs = out$coefs,
    term = out$term,
    estimate = out$pooled_estimate,
    se = out$pooled_se,
    tau2 = out$tau2,
    converged = fit_obj$converged
  )
}


.summarise_simulation_results <- function(results_df) {
  split_key <- paste(results_df$method, results_df$coefs, sep = "||")
  split_idx <- split(seq_len(nrow(results_df)), split_key)

  summary_list <- lapply(split_idx, function(idx) {
    d <- results_df[idx, , drop = FALSE]
    finite_est <- is.finite(d$estimate)
    finite_err <- is.finite(d$estimate) & is.finite(d$target)
    finite_se <- is.finite(d$se)
    finite_cov <- !is.na(d$coverage_95)
    finite_tau2 <- is.finite(d$tau2)

    tibble::tibble(
      method = d$method[[1]],
      coefs = d$coefs[[1]],
      term = d$term[[1]],
      n_runs = nrow(d),
      n_success = sum(finite_est),
      mean_estimate = if (any(finite_est)) mean(d$estimate[finite_est]) else NA_real_,
      mean_target = if (any(is.finite(d$target))) mean(d$target[is.finite(d$target)]) else NA_real_,
      bias = if (any(finite_err)) mean(d$estimate[finite_err] - d$target[finite_err]) else NA_real_,
      rmse = if (any(finite_err)) sqrt(mean((d$estimate[finite_err] - d$target[finite_err])^2)) else NA_real_,
      mean_se = if (any(finite_se)) mean(d$se[finite_se]) else NA_real_,
      coverage_95 = if (any(finite_cov)) mean(d$coverage_95[finite_cov]) else NA_real_,
      mean_tau2 = if (any(finite_tau2)) mean(d$tau2[finite_tau2]) else NA_real_,
      convergence_rate = mean(d$converged, na.rm = TRUE),
      error_rate = mean(!is.na(d$error))
    )
  })

  dplyr::bind_rows(summary_list)
}


#' Run a Simulation Benchmark for ERGMeta Methods
#'
#' @param n_sims Number of simulation replications.
#' @param n_networks Number of networks per replication.
#' @param coef_terms Character vector of coefficient term names.
#' @param beta_intercept True coefficient means at moderator value 0.
#' @param beta_moderator True moderator slopes for each coefficient.
#' @param tau Between-network standard deviations for each coefficient.
#' @param tau_cor Common correlation used for between-network covariance.
#' @param within_sd Baseline within-network standard deviations for each coefficient.
#' @param within_cor Common correlation used for within-network covariance.
#' @param within_sd_jitter Log-normal jitter scale for network-specific within SDs.
#' @param moderator_name Name of the simulated network-level moderator.
#' @param moderator_mean Mean of the moderator distribution.
#' @param moderator_sd Standard deviation of the moderator distribution.
#' @param exclusion_prob Probability that a network is excluded (to emulate failed estimation).
#' @param methods Character vector selecting methods to run. Options:
#'   \code{"sens_fixed"}, \code{"sens_random"}, \code{"sens_random_mod"},
#'   \code{"cov_random"}, \code{"cov_random_mod"}.
#' @param progress Logical. If \code{TRUE}, prints per-iteration progress.
#' @param seed Optional integer seed.
#'
#' @returns A list with per-run \code{results}, aggregated \code{summary}, and \code{settings}.
#' @export
run_ergmeta_simulation <- function(n_sims = 50,
                                   n_networks = 40,
                                   coef_terms = c("edges", "mutual"),
                                   beta_intercept = c(-2.2, 0.4),
                                   beta_moderator = c(0.3, -0.2),
                                   tau = c(0.25, 0.15),
                                   tau_cor = 0,
                                   within_sd = c(0.25, 0.20),
                                   within_cor = 0.25,
                                   within_sd_jitter = 0.15,
                                   moderator_name = "moderator_x",
                                   moderator_mean = 0,
                                   moderator_sd = 1,
                                   exclusion_prob = 0,
                                   methods = c("sens_fixed", "sens_random", "sens_random_mod", "cov_random", "cov_random_mod"),
                                   progress = TRUE,
                                   seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (!is.numeric(n_sims) || length(n_sims) != 1 || n_sims < 1) {
    stop("n_sims must be a single number >= 1.")
  }
  n_sims <- as.integer(n_sims)

  allowed_methods <- c("sens_fixed", "sens_random", "sens_random_mod", "cov_random", "cov_random_mod")
  methods <- unique(as.character(methods))
  invalid <- setdiff(methods, allowed_methods)
  if (length(invalid) > 0) {
    stop(sprintf("Unknown methods: %s", paste(invalid, collapse = ", ")))
  }
  if (length(methods) == 0) {
    stop("At least one method must be selected.")
  }
  if (!is.logical(progress) || length(progress) != 1 || is.na(progress)) {
    stop("progress must be a single TRUE/FALSE value.")
  }

  coef_ids <- make.names(coef_terms)
  if (length(beta_intercept) != length(coef_ids) || length(beta_moderator) != length(coef_ids)) {
    stop("beta_intercept and beta_moderator must match coef_terms length.")
  }

  results_list <- vector("list", n_sims)

  for (sim in seq_len(n_sims)) {
    if (isTRUE(progress)) {
      message(sprintf("Simulation %s/%s", sim, n_sims))
    }

    sim_data <- simulate_ergmeta_data(
      n_networks = n_networks,
      coef_terms = coef_terms,
      beta_intercept = beta_intercept,
      beta_moderator = beta_moderator,
      tau = tau,
      tau_cor = tau_cor,
      within_sd = within_sd,
      within_cor = within_cor,
      within_sd_jitter = within_sd_jitter,
      moderator_name = moderator_name,
      moderator_mean = moderator_mean,
      moderator_sd = moderator_sd,
      exclusion_prob = exclusion_prob
    )

    estimates <- sim_data$estimates
    covariance <- sim_data$covariance
    moderators <- sim_data$moderators
    included_networks <- sort(unique(estimates$network))
    moderators_included <- moderators[moderators$network %in% included_networks, , drop = FALSE]

    mean_moderator <- mean(moderators_included[[moderator_name]])
    target_values <- beta_intercept + beta_moderator * mean_moderator
    target_df <- tibble::tibble(
      coefs = coef_ids,
      term = coef_terms,
      target = target_values
    )

    template <- tibble::tibble(
      coefs = coef_ids,
      term = coef_terms
    )

    run_method <- function(method_name) {
      method_result <- tryCatch(
        {
          if (method_name == "sens_fixed") {
            sens <- run_exclusion_sensitivity(
              df = estimates,
              include_leave_one_out = FALSE,
              method = "fixed"
            )
            .collect_sensitivity_baseline(sens)
          } else if (method_name == "sens_random") {
            sens <- run_exclusion_sensitivity(
              df = estimates,
              include_leave_one_out = FALSE,
              method = "random"
            )
            .collect_sensitivity_baseline(sens)
          } else if (method_name == "sens_random_mod") {
            sens <- run_exclusion_sensitivity(
              df = estimates,
              include_leave_one_out = FALSE,
              method = "random",
              moderators = moderator_name,
              moderator_df = moderators_included
            )
            .collect_sensitivity_baseline(sens)
          } else if (method_name == "cov_random") {
            fit <- meta_fit_covariance(
              estimates_df = estimates,
              covariance_df = covariance,
              method = "random"
            )
            .collect_covariance_fit(fit)
          } else if (method_name == "cov_random_mod") {
            fit <- meta_fit_covariance(
              estimates_df = estimates,
              covariance_df = covariance,
              method = "random",
              moderators = moderator_name,
              moderator_df = moderators_included
            )
            .collect_covariance_fit(fit)
          } else {
            stop(sprintf("Unknown method: %s", method_name))
          }
        },
        error = function(e) {
          tibble::tibble(
            coefs = coef_ids,
            term = coef_terms,
            estimate = NA_real_,
            se = NA_real_,
            tau2 = NA_real_,
            converged = FALSE,
            error = as.character(e$message)
          )
        }
      )

      if (!("error" %in% names(method_result))) {
        method_result$error <- NA_character_
      }
      method_result <- dplyr::left_join(template, method_result, by = c("coefs", "term"))
      method_result$method <- method_name
      method_result
    }

    method_rows <- dplyr::bind_rows(lapply(methods, run_method))
    method_rows <- dplyr::left_join(method_rows, target_df, by = c("coefs", "term"))
    method_rows$sim <- sim
    method_rows$n_included_networks <- length(included_networks)
    method_rows$ci_low <- method_rows$estimate - 1.96 * method_rows$se
    method_rows$ci_high <- method_rows$estimate + 1.96 * method_rows$se
    method_rows$coverage_95 <- ifelse(
      is.finite(method_rows$ci_low) & is.finite(method_rows$ci_high) & is.finite(method_rows$target),
      method_rows$ci_low <= method_rows$target & method_rows$ci_high >= method_rows$target,
      NA
    )

    results_list[[sim]] <- method_rows
  }

  results_df <- dplyr::bind_rows(results_list)
  summary_df <- .summarise_simulation_results(results_df)

  out <- list(
    results = results_df,
    summary = summary_df,
    settings = list(
      n_sims = n_sims,
      n_networks = n_networks,
      coef_terms = coef_terms,
      beta_intercept = beta_intercept,
      beta_moderator = beta_moderator,
      tau = tau,
      tau_cor = tau_cor,
      within_sd = within_sd,
      within_cor = within_cor,
      within_sd_jitter = within_sd_jitter,
      moderator_name = moderator_name,
      moderator_mean = moderator_mean,
      moderator_sd = moderator_sd,
      exclusion_prob = exclusion_prob,
      methods = methods
    )
  )
  class(out) <- c("ergmeta_simulation_benchmark", class(out))
  out
}
