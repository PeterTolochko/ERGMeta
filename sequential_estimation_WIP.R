require(tidyverse)
require(statnet)
require(brms)


rewrite_formula <- function(formula, network) {
  network_name <- deparse(substitute(network))
  terms <- toString(formula[-c(1:2)]) # everything but "network_name ~"
  # Construct the new formula as a character string
  new_formula <- as.formula(paste(network_name, "~", terms))
  return(new_formula)
}

ergm_estimation <- function(network,
                            formula,
                            control_settings = NULL) {
  current_formula <- rewrite_formula(formula, network)
  out <- tryCatch(
    {
      if (is.null(control_settings)) {
        ergm(current_formula)
      } else {
        ergm(current_formula, control = control_settings)
      }
    },
    error=function(cond) {
      message("\n\nModel does not converge!!!\n\n")
      message(cond)
      return(NA)
    } 
  )
  return(out)
}

sequential_estimation <- function(network_list,
                                  formula,
                                  control_settings = NULL) {
  out_models <- lapply(network_list,
                       function(net) ergm_estimation(net,
                                                     formula,
                                                     control_settings))
  non_null_models <- !is.na(out_models)
  
  if (sum(non_null_models) != length(out_models)) {
    out_models <- out_models[non_null_models]
    cat(paste0("Network ", which(non_null_models == FALSE), " could not be estimated."))
    cat("\nThese networks were removed.\n")
  } else {
    cat("All networks were estimated successfully!")
  }
  
  non_infinite_est <- sapply(out_models, function(x) !is.infinite(summary(x)$coefficients[ ,1]))
  non_infinite_est <- colSums(non_infinite_est)
  infinite_est <- non_infinite_est < max(non_infinite_est)
  
  
  if (sum(infinite_est) == 0) {
    out_models <- out_models
  } else {
    out_models <- out_models[!infinite_est]
    cat(paste0("\nIn addition, network ", which(infinite_est == TRUE), " contains infinite estimates and was removed."))
    cat("\nConsider re-parametrising the networks if you want to keep more networks in.")
  }
  return(out_models)
}




extract_estimates <- function(fit_list, grouping_variable) {
  
  # extract the group from network objects
  group_vector <- sapply(fit_list, function(x) unique(x$network %v% grouping_variable))
  
  
  estimates <- t(sapply(fit_list, function(x) summary(x)$coefficients[ ,1]))
  ses <- t(sapply(fit_list, function(x) summary(x)$coefficients[ ,2]))
  names <- paste0("theta_", 0:(ncol(estimates)-1))
  colnames(estimates) <- names
  colnames(ses) <- names
  
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
  
  left_join(estimates, ses)
}



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
  out_formula <- mvbf(combined_formula,
                      rescor = TRUE)
  return(out_formula)
}





meta_fit <- function(df,
                     group_var = NULL,
                     chains = 4,
                     cores = 4,
                     iter = 4e3,
                     priors = prior(normal(0, 2), class = "b"),
                     backend = "cmdstanr") {
  
  df_wide <- df %>%
    pivot_wider(values_from = c(estimate, ses),
                names_from = coefs)
  
  if (is.null(group_var)) {
    current_formula <- create_mvbf_formula(df_wide)
  } else {
    current_formula <- create_mvbf_formula(df_wide, "group_var")
  }
  
  print("Current Formula:")
  print(current_formula)

  
  fit <- brm(current_formula,
             prior = priors, #TODO!
             iter = iter,
             data = df_wide,
             chains = chains,
             cores = cores,
             backend = backend)
  
  return(fit)
}
