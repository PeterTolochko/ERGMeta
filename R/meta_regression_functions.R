#' Function to extract estimates and standard errors from a list of fitted models
#'
#' @param fit_list A list of fitted models (e.g., ergm objects).
#' @param grouping_variable (Optional) The variable used for grouping the network objects.
#' 
#' @returns A tibble containing the estimates and the standard errors for each coefficient.
#' @export
extract_estimates <- function(fit_list, grouping_variable = NULL) {
  # Input validation
  if (!is.list(fit_list) || length(fit_list) == 0) {
    stop("fit_list must be a non-empty list of fitted models.")
  }
  if (!all(vapply(fit_list, inherits, logical(1), "ergm"))) {
    stop("All elements in fit_list must be ergm objects.")
  }

  result <- purrr::map_dfr(seq_along(fit_list), function(i) {
    current_fit <- fit_list[[i]]
    current_summary <- summary(current_fit)$coefficients
    if (is.null(current_summary) || ncol(current_summary) < 2) {
      stop(sprintf(
        "Model %s does not contain a valid coefficient matrix with estimates and standard errors.",
        i
      ))
    }

    terms <- rownames(current_summary)
    if (is.null(terms)) {
      terms <- paste0("theta_", seq_len(nrow(current_summary)) - 1L)
    }

    group_value <- NA_character_
    if (!is.null(grouping_variable)) {
      group_value <- network::`%n%`(current_fit$network, grouping_variable)
      if (is.null(group_value) || length(group_value) == 0) {
        group_value <- NA_character_
      } else {
        group_value <- as.character(group_value[[1]])
      }
    }

    tibble::tibble(
      network = i,
      coefs = make.names(terms),
      term = terms,
      estimate = unname(current_summary[, 1]),
      ses = unname(current_summary[, 2]),
      group_var = group_value
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




#' Function to create a multivariate Bayesian formula for a given data frame.
#' This is an internal function.
#' @param df The data frame containing the estimates and SEs
#' @param group_vars Optional vector of grouping variables
#' @param sigma (Optional) Logical indicating whether to include sigma in the formula
#' @returns The multivariate formula for the brm function from brms package
#' @export
create_mvbf_formula <- function(df, group_vars = NULL, sigma = FALSE) {
  if (!is.logical(sigma) || length(sigma) != 1 || is.na(sigma)) {
    stop("sigma must be a single TRUE/FALSE value.")
  }
  if (!is.null(group_vars) && !is.character(group_vars)) {
    stop("group_vars must be NULL or a character vector.")
  }

  estimate_cols <- grep("^estimate_", colnames(df), value = TRUE)
  if (length(estimate_cols) == 0) {
    stop("df must contain one or more columns prefixed with 'estimate_'.")
  }

  ses_cols <- sub("^estimate_", "ses_", estimate_cols)
  missing_ses <- setdiff(ses_cols, colnames(df))
  if (length(missing_ses) > 0) {
    stop(sprintf(
      "Missing corresponding standard-error columns: %s",
      paste(missing_ses, collapse = ", ")
    ))
  }

  if (!is.null(group_vars) && length(group_vars) > 0) {
    missing_group_vars <- setdiff(group_vars, colnames(df))
    if (length(missing_group_vars) > 0) {
      stop(sprintf(
        "Grouping variables not found in df: %s",
        paste(missing_group_vars, collapse = ", ")
      ))
    }
  }

  sigma_flag <- if (isTRUE(sigma)) "TRUE" else "FALSE"
  formula_list <- lapply(seq_along(estimate_cols), function(i) {
    estimate_col <- estimate_cols[[i]]
    ses_col <- ses_cols[[i]]
    bf_str <- sprintf(
      "`%s` | se(`%s`, sigma = %s) ~ 0 + Intercept",
      estimate_col,
      ses_col,
      sigma_flag
    )

    if (!is.null(group_vars) && length(group_vars) > 0) {
      random_effects <- paste(sprintf("(1 | `%s`)", group_vars), collapse = " + ")
      bf_str <- paste(bf_str, "+", random_effects)
    }

    brms::bf(as.formula(bf_str))
  })

  combined_formula <- Reduce(`+`, formula_list)
  out_formula <- brms::mvbf(combined_formula, rescor = TRUE)
  return(out_formula)
}


#' Function to fit a meta-analysis model using Bayesian regression
#'
#' @param df A data frame containing the data for the meta-analysis
#' @param group_vars (optional) Variables indicating the grouping structure of the data
#' @param chains The number of Markov chains to run (default: 4)
#' @param cores The number of CPU cores to use for parallel computing (default: 4)
#' @param iter The number of iterations for each chain (default: 4000)
#' @param priors The prior distribution for the regression coefficients (default: normal(0, 2))
#' @param backend The backend for running the Bayesian regression (default: "cmdstanr")
#' @param sigma (Optional) Logical indicating whether to include sigma in the formula
#' @returns The fitted model object
#' @export
meta_fit <- function(df,
                     group_vars = NULL,
                     chains = 4,
                     cores = 4,
                     iter = 4e3,
                     priors = prior("normal(0, 2)", class = "b"),
                     backend = "cmdstanr",
                     sigma = FALSE) {
  required_columns <- c("network", "coefs", "estimate", "ses")
  missing_columns <- setdiff(required_columns, colnames(df))
  if (length(missing_columns) > 0) {
    stop(sprintf(
      "df is missing required columns: %s",
      paste(missing_columns, collapse = ", ")
    ))
  }
  if (!is.null(group_vars) && !all(group_vars %in% colnames(df))) {
    stop(sprintf(
      "group_vars not found in df: %s",
      paste(setdiff(group_vars, colnames(df)), collapse = ", ")
    ))
  }
  if (any(duplicated(df[, c("network", "coefs")]))) {
    stop("df contains duplicate rows for the same network/coefficient combination.")
  }

  df_wide <- df %>%
    pivot_wider(values_from = c("estimate", "ses"),
                names_from = "coefs")

  current_formula <- create_mvbf_formula(df_wide,
                                         group_vars = group_vars,
                                         sigma = sigma)

  fit <- brms::brm(current_formula,
             prior = priors,
             iter = iter,
             data = df_wide,
             chains = chains,
             cores = cores,
             backend = backend)
  
  return(fit)
}
