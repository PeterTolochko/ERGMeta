#' Function to extract estimates and standard errors from a list of fitted models
#'
#' @param fit_list A list of fitted models (e.g., ergm objects).
#' @param grouping_variable (Optional) The variable used for grouping the network objects.
#' 
#' @returns A tibble containing the estimates and the standard errors for each coefficient.
#' @export
extract_estimates <- function(fit_list, grouping_variable = NULL) {
  # Input validation
  if (length(fit_list) == 0 || !is.list(fit_list)) {
    stop("fit_list must be a non-empty list of fitted models.")
  }
  if (!all(sapply(fit_list, inherits, "ergm"))) {
    stop("All elements in fit_list must be ergm objects.")
  }
  # if (!is.null(grouping_variable) && !all(sapply(fit_list, function(x) grouping_variable %in% names(x$network)))) {
  #   stop("The specified grouping_variable is not found in all network objects.")
  # }
  
  # Extract the group from network objects
  if (!is.null(grouping_variable)) {
    group_vector <- sapply(fit_list, function(x) unique(x$network %v% grouping_variable))
  } else {
    group_vector <- rep(NA, length(fit_list))
  }
  
  # Extract estimates and standard errors from each fitted model
  estimates <- t(sapply(fit_list, function(x) summary(x)$coefficients[ ,1]))
  ses <- t(sapply(fit_list, function(x) summary(x)$coefficients[ ,2]))
  
  # Create column names for estimates and standard errors
  names <- paste0("theta_", 0:(ncol(estimates)-1))
  colnames(estimates) <- names
  colnames(ses) <- names
  
  # Convert estimates to tibble and add grouping variables
  estimates <- estimates %>% as_tibble() %>%
    mutate(
      network = 1:nrow(estimates), # grouping variable
      group_var = group_vector
    ) %>%
    pivot_longer(
      cols = starts_with("theta_"),
      names_to = "coefs"
    ) %>%
    mutate(estimate = value) %>%
    select(-value)
  
  # Convert standard errors to tibble and add grouping variables
  ses <- ses %>% as_tibble() %>%
    mutate(
      network = 1:nrow(ses),
      group_var = group_vector
    ) %>%
    pivot_longer(
      cols = starts_with("theta_"),
      names_to = "coefs"
    ) %>%
    mutate(ses = value) %>%
    select(-value)
  
  # remove group column if no group supplied
  if (is.null(grouping_variable)) {
    estimates <- estimates %>%
      select(-group_var)
    ses <- ses %>%
      select(-group_var)
  }
  
  # Join estimates and standard errors
  result <- left_join(estimates, ses)
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
  num_pairs <- length(grep("^estimate", colnames(df)))
  
  if (is.null(group_vars)) {
    formula_list <- lapply(0:(num_pairs - 1), function(i) {
      estimate_col <- paste("estimate_theta_", i, sep = "")
      ses_col <- paste("ses_theta_", i, sep = "")
      bf_str <- sprintf("%s | se(%s, sigma = %s) ~ 0 + Intercept", estimate_col, ses_col, sigma)
      return(bf_str)
    })
    
  } else {
    
    formula_list <- lapply(0:(num_pairs - 1), function(i) {
      estimate_col <- paste("estimate_theta_", i, sep = "")
      ses_col <- paste("ses_theta_", i, sep = "")
      bf_str <- sprintf("%s | se(%s, sigma = %s) ~ 0 + Intercept ", estimate_col, ses_col, sigma)
      
      for (g in 1:length(group_vars)) {
        bf_str <- paste0(bf_str, sprintf("+ (1 | %s)", group_vars[g]))
      }
      
      return(bf_str)
    })
  }
  
  formula_list <- lapply(1:length(formula_list), function(x) bf(as.formula(formula_list[[x]])))
  combined_formula <- Reduce(`+`, formula_list)
  out_formula <- brms::mvbf(combined_formula,
                      rescor = TRUE)
  return(out_formula)
}


#' Function to fit a meta-analysis model using Bayesian regression
#'
#' @param df A data frame containing the data for the meta-analysis
#' @param group_var (optional) The variable indicating the grouping structure of the data
#' @param chains The number of Markov chains to run (default: 4)
#' @param cores The number of CPU cores to use for parallel computing (default: 4)
#' @param iter The number of iterations for each chain (default: 4000)
#' @param priors The prior distribution for the regression coefficients (default: normal(0, 2))
#' @param backend The backend for running the Bayesian regression (default: "cmdstanr")
#' @returns The fitted model object
#' @export
meta_fit <- function(df,
                     group_vars = NULL,
                     chains = 4,
                     cores = 4,
                     iter = 4e3,
                     priors = prior(normal(0, 2), class = "b"),
                     backend = "cmdstanr") {
  
  df_wide <- df %>%
    pivot_wider(values_from = c(estimate, ses),
                names_from = coefs)
  
  if (is.null(group_vars)) {
    current_formula <- create_mvbf_formula(df_wide)
  } else {
    current_formula <- create_mvbf_formula(df_wide, group_vars)
  }
  
  print("Current Formula:")
  print(current_formula)
  
  fit <- brms::brm(current_formula,
             prior = priors,
             iter = iter,
             data = df_wide,
             chains = chains,
             cores = cores,
             backend = backend)
  
  return(fit)
}
