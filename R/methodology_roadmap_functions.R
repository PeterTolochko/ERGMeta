#' Extract Within-Network Covariance for ERGM Estimates
#'
#' @param fit_list A non-empty list of fitted \code{ergm} model objects.
#' @param grouping_variable (Optional) A network-level attribute used for grouping.
#' @param strict_terms Logical. If \code{TRUE}, all models must have identical coefficient terms.
#'
#' @returns A tibble with one row per covariance cell per network, including
#'   coefficient ids, original terms, covariance values, and an observation flag.
#' @export
extract_within_covariance <- function(fit_list,
                                      grouping_variable = NULL,
                                      strict_terms = TRUE) {
  if (!is.list(fit_list) || length(fit_list) == 0) {
    stop("fit_list must be a non-empty list of fitted models.")
  }
  if (!all(vapply(fit_list, inherits, logical(1), "ergm"))) {
    stop("All elements in fit_list must be ergm objects.")
  }
  if (!is.logical(strict_terms) || length(strict_terms) != 1 || is.na(strict_terms)) {
    stop("strict_terms must be a single TRUE/FALSE value.")
  }

  model_info <- lapply(seq_along(fit_list), function(i) {
    current_fit <- fit_list[[i]]
    current_vcov <- tryCatch(stats::vcov(current_fit), error = function(e) NULL)

    if (is.null(current_vcov) || !is.matrix(current_vcov)) {
      stop(sprintf("Model %s does not provide a valid covariance matrix.", i))
    }
    if (nrow(current_vcov) != ncol(current_vcov)) {
      stop(sprintf("Model %s covariance matrix is not square.", i))
    }

    current_terms <- colnames(current_vcov)
    if (is.null(current_terms)) {
      current_terms <- rownames(current_vcov)
    }
    if (is.null(current_terms)) {
      current_terms <- names(stats::coef(current_fit))
    }
    if (is.null(current_terms) || length(current_terms) != nrow(current_vcov)) {
      current_terms <- paste0("theta_", seq_len(nrow(current_vcov)) - 1L)
    }

    dimnames(current_vcov) <- list(current_terms, current_terms)

    group_value <- NA_character_
    if (!is.null(grouping_variable)) {
      group_value <- network::`%n%`(current_fit$network, grouping_variable)
      if (is.null(group_value) || length(group_value) == 0) {
        group_value <- NA_character_
      } else {
        group_value <- as.character(group_value[[1]])
      }
    }

    list(
      vcov = current_vcov,
      terms = current_terms,
      group_value = group_value
    )
  })

  term_sets <- lapply(model_info, function(x) x$terms)
  if (isTRUE(strict_terms)) {
    reference_terms <- term_sets[[1]]
    mismatched <- which(!vapply(term_sets, function(x) identical(x, reference_terms), logical(1)))
    if (length(mismatched) > 0) {
      stop(sprintf(
        "Coefficient terms differ across fitted models (first mismatch at model %s). Set strict_terms = FALSE to align by term union.",
        mismatched[[1]]
      ))
    }
    target_terms <- reference_terms
  } else {
    target_terms <- Reduce(union, term_sets)
  }

  result <- purrr::map_dfr(seq_along(model_info), function(i) {
    current_vcov <- model_info[[i]]$vcov
    current_terms <- model_info[[i]]$terms

    aligned_vcov <- matrix(
      NA_real_,
      nrow = length(target_terms),
      ncol = length(target_terms),
      dimnames = list(target_terms, target_terms)
    )
    aligned_vcov[current_terms, current_terms] <- current_vcov

    current_df <- as.data.frame(as.table(aligned_vcov), stringsAsFactors = FALSE)
    colnames(current_df) <- c("term_row", "term_col", "covariance")
    current_df$network <- i
    current_df$coefs_row <- make.names(current_df$term_row)
    current_df$coefs_col <- make.names(current_df$term_col)
    current_df$covariance <- as.numeric(current_df$covariance)
    current_df$is_observed <- !is.na(current_df$covariance)
    current_df$group_var <- model_info[[i]]$group_value

    tibble::as_tibble(
      current_df[
        ,
        c(
          "network",
          "coefs_row",
          "coefs_col",
          "term_row",
          "term_col",
          "covariance",
          "is_observed",
          "group_var"
        )
      ]
    )
  })

  if (!is.null(grouping_variable) && any(is.na(result$group_var))) {
    stop("The specified grouping_variable is not found as a network attribute in all fitted models.")
  }

  if (is.null(grouping_variable)) {
    result$group_var <- NULL
  }

  return(result)
}


.extract_moderator_table <- function(df, moderators, moderator_df = NULL) {
  network_ids <- sort(unique(df$network))

  if (is.null(moderators) || length(moderators) == 0) {
    return(tibble::tibble(network = network_ids))
  }
  if (!is.character(moderators)) {
    stop("moderators must be NULL or a character vector.")
  }

  if (is.null(moderator_df)) {
    missing_cols <- setdiff(c("network", moderators), colnames(df))
    if (length(missing_cols) > 0) {
      stop(sprintf(
        "Moderator columns not found in df: %s",
        paste(missing_cols, collapse = ", ")
      ))
    }
    source_df <- df[, c("network", moderators), drop = FALSE]
  } else {
    if (!is.data.frame(moderator_df)) {
      stop("moderator_df must be NULL or a data frame.")
    }
    missing_cols <- setdiff(c("network", moderators), colnames(moderator_df))
    if (length(missing_cols) > 0) {
      stop(sprintf(
        "Moderator columns not found in moderator_df: %s",
        paste(missing_cols, collapse = ", ")
      ))
    }
    source_df <- moderator_df[, c("network", moderators), drop = FALSE]
  }

  moderator_table <- unique(source_df)
  if (any(duplicated(moderator_table$network))) {
    stop("Moderator values must be unique per network.")
  }
  if (!all(network_ids %in% moderator_table$network)) {
    missing_networks <- setdiff(network_ids, moderator_table$network)
    stop(sprintf(
      "Missing moderator values for networks: %s",
      paste(missing_networks, collapse = ", ")
    ))
  }

  moderator_table <- moderator_table[match(network_ids, moderator_table$network), , drop = FALSE]
  if (any(!stats::complete.cases(moderator_table))) {
    stop("Moderator table contains missing values.")
  }

  tibble::as_tibble(moderator_table)
}


.build_design_matrix <- function(moderator_table,
                                 moderators = NULL,
                                 network_ids = NULL,
                                 design_cols = NULL) {
  if (is.null(network_ids)) {
    network_ids <- moderator_table$network
  }

  current_df <- moderator_table[match(network_ids, moderator_table$network), , drop = FALSE]
  if (any(is.na(current_df$network))) {
    stop("Some requested networks are missing from the moderator table.")
  }

  if (is.null(moderators) || length(moderators) == 0) {
    x <- matrix(1, nrow = nrow(current_df), ncol = 1)
    colnames(x) <- "(Intercept)"
    return(list(
      X = x,
      design_cols = colnames(x),
      network_ids = network_ids,
      formula = NULL,
      design_df = current_df[, "network", drop = FALSE]
    ))
  }

  design_formula <- stats::as.formula(
    paste("~", paste(moderators, collapse = " + "))
  )
  x_raw <- stats::model.matrix(
    design_formula,
    data = current_df[, moderators, drop = FALSE]
  )

  if (is.null(design_cols)) {
    x <- x_raw
    design_cols <- colnames(x)
  } else {
    x <- matrix(0, nrow = nrow(x_raw), ncol = length(design_cols))
    colnames(x) <- design_cols
    common <- intersect(colnames(x_raw), design_cols)
    x[, common] <- x_raw[, common, drop = FALSE]
  }

  list(
    X = x,
    design_cols = design_cols,
    network_ids = network_ids,
    formula = design_formula,
    design_df = current_df[, c("network", moderators), drop = FALSE]
  )
}


.build_reference_vector <- function(moderator_table,
                                    moderators = NULL,
                                    design_cols = "(Intercept)",
                                    reference_moderators = NULL) {
  if (is.null(moderators) || length(moderators) == 0) {
    x <- matrix(1, nrow = 1, ncol = 1)
    colnames(x) <- "(Intercept)"
    return(list(
      x = x,
      reference_profile = "intercept_only"
    ))
  }

  if (is.null(reference_moderators)) {
    reference_values <- lapply(moderators, function(m) {
      v <- moderator_table[[m]]
      if (is.numeric(v)) {
        mean(v)
      } else if (is.factor(v)) {
        levels(v)[1]
      } else {
        unique(v)[1]
      }
    })
    names(reference_values) <- moderators
  } else if (is.data.frame(reference_moderators)) {
    reference_values <- as.list(reference_moderators[1, moderators, drop = FALSE])
  } else if (is.list(reference_moderators)) {
    reference_values <- reference_moderators
  } else {
    stop("reference_moderators must be NULL, a named list, or a one-row data frame.")
  }

  missing_names <- setdiff(moderators, names(reference_values))
  if (length(missing_names) > 0) {
    stop(sprintf(
      "reference_moderators is missing: %s",
      paste(missing_names, collapse = ", ")
    ))
  }

  ref_df <- moderator_table[1, moderators, drop = FALSE]
  for (m in moderators) {
    if (is.factor(ref_df[[m]])) {
      ref_df[[m]] <- factor(
        as.character(reference_values[[m]]),
        levels = levels(moderator_table[[m]])
      )
    } else {
      ref_df[[m]] <- reference_values[[m]]
    }
  }

  design_formula <- stats::as.formula(
    paste("~", paste(moderators, collapse = " + "))
  )
  x_raw <- stats::model.matrix(design_formula, data = ref_df)
  x <- matrix(0, nrow = 1, ncol = length(design_cols))
  colnames(x) <- design_cols
  common <- intersect(colnames(x_raw), design_cols)
  x[1, common] <- x_raw[1, common, drop = TRUE]

  profile_text <- paste(
    paste0(moderators, "=", unlist(ref_df[1, moderators, drop = FALSE])),
    collapse = ", "
  )

  list(
    x = x,
    reference_profile = profile_text
  )
}


.safe_inverse <- function(v) {
  chol_obj <- tryCatch(chol(v), error = function(e) NULL)
  if (is.null(chol_obj)) {
    return(NULL)
  }
  list(
    inv = chol2inv(chol_obj),
    logdet = 2 * sum(log(diag(chol_obj)))
  )
}


.stabilize_and_invert <- function(v, jitter = 1e-10, max_attempts = 8) {
  for (k in 0:max_attempts) {
    offset <- jitter * (10^k)
    v_try <- v + diag(offset, nrow(v))
    out <- .safe_inverse(v_try)
    if (!is.null(out)) {
      return(out)
    }
  }
  return(NULL)
}


.fit_univariate_meta <- function(y, s2, x, method = c("random", "fixed")) {
  method <- match.arg(method)
  n <- length(y)

  evaluate <- function(tau2) {
    v <- s2 + tau2
    if (any(!is.finite(v) | v <= 0)) {
      return(NULL)
    }
    w <- 1 / v
    xtwx <- crossprod(x, x * w)
    xtwy <- crossprod(x, y * w)

    beta <- tryCatch(as.vector(solve(xtwx, xtwy)), error = function(e) NULL)
    if (is.null(beta)) {
      return(NULL)
    }
    vcov_beta <- tryCatch(solve(xtwx), error = function(e) NULL)
    if (is.null(vcov_beta)) {
      return(NULL)
    }

    residual <- y - as.vector(x %*% beta)
    nll <- 0.5 * sum(log(v) + (residual^2) / v + log(2 * pi))

    list(
      beta = beta,
      vcov_beta = vcov_beta,
      tau2 = tau2,
      nll = nll
    )
  }

  if (method == "fixed" || n <= 1) {
    out <- evaluate(0)
    if (is.null(out)) {
      stop("Failed to fit fixed-effects model.")
    }
    out$converged <- TRUE
    return(out)
  }

  start_tau2 <- max(1e-8, stats::var(y) - mean(s2))
  opt <- nlminb(
    start = log(start_tau2),
    objective = function(log_tau2) {
      fit <- evaluate(exp(log_tau2))
      if (is.null(fit)) {
        return(Inf)
      }
      fit$nll
    }
  )

  out <- evaluate(exp(opt$par))
  if (is.null(out)) {
    stop("Failed to fit random-effects model.")
  }
  out$converged <- opt$convergence == 0
  out
}


.fit_multivariate_meta <- function(y_list, s_list, x, method = c("random", "fixed")) {
  method <- match.arg(method)

  n <- length(y_list)
  p <- length(y_list[[1]])
  q <- ncol(x)

  evaluate <- function(tau2_vec) {
    if (any(!is.finite(tau2_vec) | tau2_vec < 0)) {
      return(NULL)
    }

    xtvix <- matrix(0, nrow = p * q, ncol = p * q)
    xtviy <- numeric(p * q)
    logdet_sum <- 0
    v_inv_list <- vector("list", n)

    for (i in seq_len(n)) {
      v_i <- s_list[[i]] + diag(tau2_vec, p)
      inv_info <- .stabilize_and_invert(v_i)
      if (is.null(inv_info)) {
        return(NULL)
      }
      v_inv_list[[i]] <- inv_info$inv
      logdet_sum <- logdet_sum + inv_info$logdet

      x_i <- kronecker(diag(p), x[i, , drop = FALSE])
      xtvix <- xtvix + crossprod(x_i, inv_info$inv %*% x_i)
      xtviy <- xtviy + crossprod(x_i, inv_info$inv %*% y_list[[i]])
    }

    beta_vec <- tryCatch(as.vector(solve(xtvix, xtviy)), error = function(e) NULL)
    if (is.null(beta_vec)) {
      return(NULL)
    }
    vcov_beta <- tryCatch(solve(xtvix), error = function(e) NULL)
    if (is.null(vcov_beta)) {
      return(NULL)
    }

    quad <- 0
    for (i in seq_len(n)) {
      x_i <- kronecker(diag(p), x[i, , drop = FALSE])
      residual <- y_list[[i]] - as.vector(x_i %*% beta_vec)
      quad <- quad + as.numeric(crossprod(residual, v_inv_list[[i]] %*% residual))
    }

    nll <- 0.5 * (n * p * log(2 * pi) + logdet_sum + quad)

    list(
      beta_vec = beta_vec,
      vcov_beta = vcov_beta,
      tau2 = tau2_vec,
      nll = nll
    )
  }

  if (method == "fixed" || n <= 1) {
    out <- evaluate(rep(0, p))
    if (is.null(out)) {
      stop("Failed to fit fixed-effects covariance model.")
    }
    out$converged <- TRUE
    return(out)
  }

  diag_means <- Reduce(`+`, lapply(s_list, diag)) / length(s_list)
  start_tau2 <- pmax(1e-8, diag_means * 0.05)

  opt <- nlminb(
    start = log(start_tau2),
    objective = function(log_tau2) {
      fit <- evaluate(exp(log_tau2))
      if (is.null(fit)) {
        return(Inf)
      }
      fit$nll
    }
  )

  out <- evaluate(exp(opt$par))
  if (is.null(out)) {
    stop("Failed to fit random-effects covariance model.")
  }
  out$converged <- opt$convergence == 0
  out
}


.build_meta_arrays <- function(estimates_df, covariance_df, network_ids, coef_ids) {
  term_map <- tapply(estimates_df$term, estimates_df$coefs, function(x) unique(x)[1])
  term_map <- term_map[coef_ids]

  y_list <- vector("list", length(network_ids))
  s_list <- vector("list", length(network_ids))

  for (i in seq_along(network_ids)) {
    current_network <- network_ids[[i]]

    current_est <- estimates_df[
      estimates_df$network == current_network & estimates_df$coefs %in% coef_ids,
      ,
      drop = FALSE
    ]
    if (any(duplicated(current_est$coefs))) {
      stop(sprintf("Duplicate coefficient rows found for network %s.", current_network))
    }
    if (!all(coef_ids %in% current_est$coefs)) {
      stop(sprintf(
        "Missing coefficient estimates for network %s: %s",
        current_network,
        paste(setdiff(coef_ids, current_est$coefs), collapse = ", ")
      ))
    }
    current_est <- current_est[match(coef_ids, current_est$coefs), , drop = FALSE]
    y_list[[i]] <- as.numeric(current_est$estimate)

    current_cov <- covariance_df[
      covariance_df$network == current_network &
        covariance_df$coefs_row %in% coef_ids &
        covariance_df$coefs_col %in% coef_ids,
      ,
      drop = FALSE
    ]
    cov_matrix <- matrix(
      NA_real_,
      nrow = length(coef_ids),
      ncol = length(coef_ids),
      dimnames = list(coef_ids, coef_ids)
    )
    if (nrow(current_cov) > 0) {
      row_idx <- match(current_cov$coefs_row, coef_ids)
      col_idx <- match(current_cov$coefs_col, coef_ids)
      cov_matrix[cbind(row_idx, col_idx)] <- current_cov$covariance
    }

    if (any(is.na(cov_matrix))) {
      stop(sprintf(
        "Covariance matrix is incomplete for network %s and selected coefficients.",
        current_network
      ))
    }

    cov_matrix <- (cov_matrix + t(cov_matrix)) / 2
    s_list[[i]] <- cov_matrix
  }

  list(
    y_list = y_list,
    s_list = s_list,
    term_map = term_map
  )
}


#' Prepare Meta-Analysis Inputs with Diagnostics
#'
#' @param fit_list A non-empty list that may include fitted \code{ergm} models and failed fits (\code{NULL}).
#' @param grouping_variable (Optional) A network-level attribute used for grouping.
#' @param strict_terms Logical. Passed to \code{extract_within_covariance()}.
#'
#' @returns A list containing extracted estimates, covariance data, inclusion diagnostics,
#'   and included/excluded network indices.
#' @export
prepare_meta_inputs <- function(fit_list,
                                grouping_variable = NULL,
                                strict_terms = TRUE) {
  if (!is.list(fit_list) || length(fit_list) == 0) {
    stop("fit_list must be a non-empty list.")
  }

  network_ids <- seq_along(fit_list)
  included_networks <- which(vapply(fit_list, inherits, logical(1), "ergm"))
  excluded_networks <- setdiff(network_ids, included_networks)

  if (length(included_networks) == 0) {
    stop("fit_list does not contain any ergm objects.")
  }

  included_fits <- fit_list[included_networks]
  estimates <- extract_estimates(included_fits, grouping_variable = grouping_variable)
  covariance <- extract_within_covariance(
    included_fits,
    grouping_variable = grouping_variable,
    strict_terms = strict_terms
  )

  estimates$network <- included_networks[estimates$network]
  covariance$network <- included_networks[covariance$network]

  diagnostics <- tibble::tibble(
    network = network_ids,
    status = ifelse(network_ids %in% included_networks, "estimated", "excluded"),
    included_in_meta = network_ids %in% included_networks
  )

  out <- list(
    estimates = estimates,
    covariance = covariance,
    diagnostics = diagnostics,
    included_networks = included_networks,
    excluded_networks = excluded_networks
  )
  class(out) <- c("ergmeta_meta_inputs", class(out))

  return(out)
}


#' Fit a Covariance-Aware Multivariate Meta-Regression
#'
#' @param estimates_df A data frame produced by \code{extract_estimates()}.
#' @param covariance_df A data frame produced by \code{extract_within_covariance()}.
#' @param method Estimation mode: \code{"random"} (default) or \code{"fixed"}.
#' @param moderators Optional character vector of network-level moderator names.
#' @param moderator_df Optional data frame with one row per network and moderator columns.
#' @param reference_moderators Optional named list or one-row data frame for reference predictions.
#' @param coef_subset Optional character vector of coefficient ids to include.
#' @param network_subset Optional integer vector of network ids to include.
#'
#' @returns A list containing coefficient-level meta-regression estimates,
#'   heterogeneity estimates (\code{tau}), and fitted values by network/term.
#' @export
meta_fit_covariance <- function(estimates_df,
                                covariance_df,
                                method = c("random", "fixed"),
                                moderators = NULL,
                                moderator_df = NULL,
                                reference_moderators = NULL,
                                coef_subset = NULL,
                                network_subset = NULL) {
  method <- match.arg(method)

  est_required <- c("network", "coefs", "term", "estimate")
  cov_required <- c("network", "coefs_row", "coefs_col", "covariance")
  missing_est <- setdiff(est_required, colnames(estimates_df))
  missing_cov <- setdiff(cov_required, colnames(covariance_df))

  if (length(missing_est) > 0) {
    stop(sprintf(
      "estimates_df is missing required columns: %s",
      paste(missing_est, collapse = ", ")
    ))
  }
  if (length(missing_cov) > 0) {
    stop(sprintf(
      "covariance_df is missing required columns: %s",
      paste(missing_cov, collapse = ", ")
    ))
  }

  common_networks <- intersect(
    sort(unique(estimates_df$network)),
    sort(unique(covariance_df$network))
  )
  if (!is.null(network_subset)) {
    common_networks <- intersect(common_networks, as.integer(network_subset))
  }
  if (length(common_networks) == 0) {
    stop("No common networks available between estimates_df and covariance_df.")
  }

  available_coefs <- sort(unique(estimates_df$coefs))
  if (!is.null(coef_subset)) {
    coef_ids <- as.character(coef_subset)
    missing_coefs <- setdiff(coef_ids, available_coefs)
    if (length(missing_coefs) > 0) {
      stop(sprintf(
        "Requested coef_subset not found: %s",
        paste(missing_coefs, collapse = ", ")
      ))
    }
  } else {
    coef_ids <- available_coefs
  }

  est_sub <- estimates_df[
    estimates_df$network %in% common_networks &
      estimates_df$coefs %in% coef_ids,
    ,
    drop = FALSE
  ]
  cov_sub <- covariance_df[
    covariance_df$network %in% common_networks &
      covariance_df$coefs_row %in% coef_ids &
      covariance_df$coefs_col %in% coef_ids,
    ,
    drop = FALSE
  ]

  arrays <- .build_meta_arrays(
    estimates_df = est_sub,
    covariance_df = cov_sub,
    network_ids = common_networks,
    coef_ids = coef_ids
  )

  moderator_table <- .extract_moderator_table(
    df = est_sub,
    moderators = moderators,
    moderator_df = moderator_df
  )
  design <- .build_design_matrix(
    moderator_table = moderator_table,
    moderators = moderators,
    network_ids = common_networks
  )

  fit <- .fit_multivariate_meta(
    y_list = arrays$y_list,
    s_list = arrays$s_list,
    x = design$X,
    method = method
  )

  p <- length(coef_ids)
  q <- ncol(design$X)
  beta_mat <- matrix(fit$beta_vec, nrow = p, byrow = TRUE)
  rownames(beta_mat) <- coef_ids
  colnames(beta_mat) <- colnames(design$X)

  beta_se <- sqrt(diag(fit$vcov_beta))
  beta_se_mat <- matrix(beta_se, nrow = p, byrow = TRUE)
  rownames(beta_se_mat) <- coef_ids
  colnames(beta_se_mat) <- colnames(design$X)

  coefficients <- purrr::map_dfr(seq_len(p), function(j) {
    tibble::tibble(
      coefs = coef_ids[[j]],
      term = arrays$term_map[[coef_ids[[j]]]],
      predictor = colnames(design$X),
      estimate = as.numeric(beta_mat[j, ]),
      se = as.numeric(beta_se_mat[j, ])
    )
  })

  reference <- .build_reference_vector(
    moderator_table = moderator_table,
    moderators = moderators,
    design_cols = colnames(design$X),
    reference_moderators = reference_moderators
  )

  reference_summary <- purrr::map_dfr(seq_len(p), function(j) {
    idx <- ((j - 1) * q + 1):(j * q)
    vcov_block <- fit$vcov_beta[idx, idx, drop = FALSE]
    mean_ref <- as.numeric(reference$x %*% beta_mat[j, ])
    se_ref <- sqrt(as.numeric(reference$x %*% vcov_block %*% t(reference$x)))

    tibble::tibble(
      coefs = coef_ids[[j]],
      term = arrays$term_map[[coef_ids[[j]]]],
      reference_profile = reference$reference_profile,
      pooled_estimate = mean_ref,
      pooled_se = se_ref
    )
  })

  fitted_values <- purrr::map_dfr(seq_along(common_networks), function(i) {
    x_i <- design$X[i, , drop = FALSE]
    x_big <- kronecker(diag(p), x_i)
    fitted_i <- as.vector(x_big %*% fit$beta_vec)
    observed_i <- arrays$y_list[[i]]

    tibble::tibble(
      network = common_networks[[i]],
      coefs = coef_ids,
      term = as.character(arrays$term_map[coef_ids]),
      observed = observed_i,
      fitted = fitted_i,
      residual = observed_i - fitted_i
    )
  })

  tau <- tibble::tibble(
    coefs = coef_ids,
    term = as.character(arrays$term_map[coef_ids]),
    tau2 = as.numeric(fit$tau2),
    tau = sqrt(as.numeric(fit$tau2))
  )

  out <- list(
    method = method,
    converged = isTRUE(fit$converged),
    n_networks = length(common_networks),
    predictors = colnames(design$X),
    coefficients = coefficients,
    reference_summary = reference_summary,
    tau = tau,
    fitted = fitted_values
  )
  class(out) <- c("ergmeta_covariance_fit", class(out))

  return(out)
}


#' Run Exclusion Sensitivity for Pooled ERGM Coefficients
#'
#' @param df A data frame produced by \code{extract_estimates()}.
#' @param exclude_sets Optional named list of integer vectors. Each vector is a set of network ids to exclude.
#' @param include_leave_one_out Logical. If \code{TRUE}, automatically add leave-one-network-out scenarios.
#' @param method Estimation mode: \code{"random"} (default) or \code{"fixed"}.
#' @param moderators Optional character vector of network-level moderator names.
#' @param moderator_df Optional data frame with one row per network and moderator columns.
#' @param reference_moderators Optional named list or one-row data frame for reference predictions.
#'
#' @returns A tibble with pooled estimates and heterogeneity summaries by scenario and coefficient.
#' @export
run_exclusion_sensitivity <- function(df,
                                      exclude_sets = NULL,
                                      include_leave_one_out = TRUE,
                                      method = c("random", "fixed"),
                                      moderators = NULL,
                                      moderator_df = NULL,
                                      reference_moderators = NULL) {
  method <- match.arg(method)

  required_columns <- c("network", "coefs", "estimate", "ses")
  missing_columns <- setdiff(required_columns, colnames(df))
  if (length(missing_columns) > 0) {
    stop(sprintf(
      "df is missing required columns: %s",
      paste(missing_columns, collapse = ", ")
    ))
  }
  if (!is.logical(include_leave_one_out) ||
      length(include_leave_one_out) != 1 ||
      is.na(include_leave_one_out)) {
    stop("include_leave_one_out must be a single TRUE/FALSE value.")
  }
  if (any(!is.finite(df$ses) | df$ses <= 0)) {
    stop("All ses values must be finite and strictly positive for pooling.")
  }

  scenario_sets <- list(all = integer(0))
  network_ids <- sort(unique(df$network))

  if (isTRUE(include_leave_one_out) && length(network_ids) > 1) {
    loo_sets <- lapply(network_ids, function(x) as.integer(x))
    names(loo_sets) <- paste0("leave_out_", network_ids)
    scenario_sets <- c(scenario_sets, loo_sets)
  }

  if (!is.null(exclude_sets)) {
    if (!is.list(exclude_sets)) {
      stop("exclude_sets must be NULL or a named list of integer vectors.")
    }
    if (is.null(names(exclude_sets)) || any(names(exclude_sets) == "")) {
      names(exclude_sets) <- paste0("custom_", seq_along(exclude_sets))
    }
    if (any(names(exclude_sets) %in% names(scenario_sets))) {
      stop("exclude_sets scenario names overlap with built-in scenario names.")
    }

    exclude_sets <- lapply(exclude_sets, function(x) {
      if (length(x) == 0) {
        return(integer(0))
      }
      as.integer(unique(x))
    })
    scenario_sets <- c(scenario_sets, exclude_sets)
  }

  moderator_table <- .extract_moderator_table(
    df = df,
    moderators = moderators,
    moderator_df = moderator_df
  )
  global_design <- .build_design_matrix(
    moderator_table = moderator_table,
    moderators = moderators,
    network_ids = network_ids
  )
  reference <- .build_reference_vector(
    moderator_table = moderator_table,
    moderators = moderators,
    design_cols = global_design$design_cols,
    reference_moderators = reference_moderators
  )

  scenario_summary <- purrr::map_dfr(names(scenario_sets), function(scenario_name) {
    excluded <- scenario_sets[[scenario_name]]
    current_df <- df[!(df$network %in% excluded), , drop = FALSE]

    if (nrow(current_df) == 0) {
      return(tibble::tibble(
        scenario = scenario_name,
        method = method,
        coefs = NA_character_,
        term = NA_character_,
        n_networks = 0L,
        pooled_estimate = NA_real_,
        pooled_se = NA_real_,
        tau2 = NA_real_,
        tau = NA_real_,
        converged = FALSE,
        reference_profile = reference$reference_profile,
        n_excluded = length(excluded),
        excluded_networks = paste(sort(excluded), collapse = ",")
      ))
    }

    current_split <- split(current_df, current_df$coefs)
    purrr::map_dfr(names(current_split), function(coef_id) {
      coef_df <- current_split[[coef_id]]
      coef_df <- coef_df[
        is.finite(coef_df$estimate) & is.finite(coef_df$ses) & coef_df$ses > 0,
        ,
        drop = FALSE
      ]
      if (nrow(coef_df) == 0) {
        return(tibble::tibble(
          scenario = scenario_name,
          method = method,
          coefs = coef_id,
          term = NA_character_,
          n_networks = 0L,
          pooled_estimate = NA_real_,
          pooled_se = NA_real_,
          tau2 = NA_real_,
          tau = NA_real_,
          converged = FALSE,
          reference_profile = reference$reference_profile,
          n_excluded = length(excluded),
          excluded_networks = paste(sort(excluded), collapse = ",")
        ))
      }

      if (any(duplicated(coef_df$network))) {
        stop(sprintf("Duplicate network rows found for coefficient %s.", coef_id))
      }

      coef_networks <- sort(unique(coef_df$network))
      coef_df <- coef_df[match(coef_networks, coef_df$network), , drop = FALSE]

      design <- .build_design_matrix(
        moderator_table = moderator_table,
        moderators = moderators,
        network_ids = coef_networks,
        design_cols = global_design$design_cols
      )

      fit <- .fit_univariate_meta(
        y = coef_df$estimate,
        s2 = coef_df$ses^2,
        x = design$X,
        method = method
      )

      pooled <- as.numeric(reference$x %*% fit$beta)
      pooled_se <- sqrt(as.numeric(reference$x %*% fit$vcov_beta %*% t(reference$x)))
      current_term <- unique(coef_df$term)
      current_term <- current_term[!is.na(current_term)]

      tibble::tibble(
        scenario = scenario_name,
        method = method,
        coefs = coef_id,
        term = if (length(current_term) > 0) current_term[[1]] else NA_character_,
        n_networks = length(coef_networks),
        pooled_estimate = pooled,
        pooled_se = pooled_se,
        tau2 = fit$tau2,
        tau = sqrt(fit$tau2),
        converged = isTRUE(fit$converged),
        reference_profile = reference$reference_profile,
        n_excluded = length(excluded),
        excluded_networks = paste(sort(excluded), collapse = ",")
      )
    })
  })

  baseline <- scenario_summary[
    scenario_summary$scenario == "all",
    c("coefs", "pooled_estimate", "tau2")
  ]
  if (nrow(baseline) > 0) {
    colnames(baseline) <- c("coefs", "baseline_estimate", "baseline_tau2")
    scenario_summary <- dplyr::left_join(scenario_summary, baseline, by = "coefs")
    scenario_summary$delta_from_all <- scenario_summary$pooled_estimate - scenario_summary$baseline_estimate
    scenario_summary$delta_tau2_from_all <- scenario_summary$tau2 - scenario_summary$baseline_tau2
    scenario_summary$baseline_estimate <- NULL
    scenario_summary$baseline_tau2 <- NULL
  } else {
    scenario_summary$delta_from_all <- NA_real_
    scenario_summary$delta_tau2_from_all <- NA_real_
  }

  return(scenario_summary)
}


#' Run a Concrete Methodological Roadmap
#'
#' @param fit_list A non-empty list that may include fitted \code{ergm} models and failed fits (\code{NULL}).
#' @param grouping_variable (Optional) A network-level attribute used for grouping.
#' @param strict_terms Logical. Passed to \code{extract_within_covariance()}.
#' @param include_leave_one_out Logical. Passed to \code{run_exclusion_sensitivity()}.
#' @param exclude_sets Optional named list of integer vectors passed to \code{run_exclusion_sensitivity()}.
#' @param sensitivity_method Estimation mode for sensitivity: \code{"random"} (default) or \code{"fixed"}.
#' @param moderators Optional character vector of network-level moderator names.
#' @param moderator_df Optional data frame with one row per network and moderator columns.
#' @param reference_moderators Optional named list or one-row data frame for reference predictions.
#' @param fit_covariance_model Logical. If \code{TRUE}, fit \code{meta_fit_covariance()}.
#' @param covariance_method Estimation mode for covariance-aware fitting: \code{"random"} (default) or \code{"fixed"}.
#'
#' @returns A list with prepared meta-analysis inputs, exclusion-sensitivity summaries,
#'   and an optional covariance-aware multivariate fit.
#' @export
methodology_roadmap <- function(fit_list,
                                grouping_variable = NULL,
                                strict_terms = TRUE,
                                include_leave_one_out = TRUE,
                                exclude_sets = NULL,
                                sensitivity_method = c("random", "fixed"),
                                moderators = NULL,
                                moderator_df = NULL,
                                reference_moderators = NULL,
                                fit_covariance_model = TRUE,
                                covariance_method = c("random", "fixed")) {
  sensitivity_method <- match.arg(sensitivity_method)
  covariance_method <- match.arg(covariance_method)

  prepared <- prepare_meta_inputs(
    fit_list = fit_list,
    grouping_variable = grouping_variable,
    strict_terms = strict_terms
  )

  sensitivity <- run_exclusion_sensitivity(
    df = prepared$estimates,
    exclude_sets = exclude_sets,
    include_leave_one_out = include_leave_one_out,
    method = sensitivity_method,
    moderators = moderators,
    moderator_df = moderator_df,
    reference_moderators = reference_moderators
  )

  covariance_meta <- NULL
  if (isTRUE(fit_covariance_model)) {
    covariance_meta <- meta_fit_covariance(
      estimates_df = prepared$estimates,
      covariance_df = prepared$covariance,
      method = covariance_method,
      moderators = moderators,
      moderator_df = moderator_df,
      reference_moderators = reference_moderators
    )
  }

  out <- list(
    estimates = prepared$estimates,
    covariance = prepared$covariance,
    diagnostics = prepared$diagnostics,
    included_networks = prepared$included_networks,
    excluded_networks = prepared$excluded_networks,
    sensitivity = sensitivity,
    covariance_meta = covariance_meta
  )
  class(out) <- c("ergmeta_methodology_roadmap", class(out))

  return(out)
}
