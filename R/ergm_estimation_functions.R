#' Function to rewrite a formula by replacing the response variable with a specified network
#'
#' This function takes an `ergm` formula and a network object and rewrites the formula
#' by replacing the response variable with the specified network.
#'
#' @param formula The original formula to be rewritten. Must be a valid `ergm` formula.
#' @param network The network object to be used as the new response variable.
#' @returns The rewritten formula with the network as the response variable.
#' @export
rewrite_formula <- function(formula, network) {
  # Input validation
  if (!inherits(formula, "formula")) {
    stop("formula must be a valid ergm formula.")
  }
  if (!inherits(network, "network")) {
    stop("network must be a valid network object.")
  }
  
  # Extract the network name
  network_name <- deparse(substitute(network))
  
  # Extract the terms from the original formula
  terms <- attr(terms(formula), "term.labels")
  
  # Construct the new formula using reformulate()
  new_formula <- reformulate(termlabels = terms, response = network_name)
  
  return(new_formula)
}


#' Function to estimate an Exponential Random Graph Model (ERGM) on a given network
#'
#' This function takes a network object, an `ergm` formula, and optional control settings,
#' and estimates an ERGM model using the `ergm` function from the `ergm` package.
#' It handles potential convergence issues and returns the estimated model object or an informative error message.
#'
#' @param network The input network object.
#' @param formula Formula specifying the ERGM model.
#' @param control_settings Optional control settings for the ERGM estimation. Default is control.ergm().
#' @returns The estimated ERGM model object or NULL if estimation fails or produces degenerate results.
#' @export
ergm_estimation <- function(network, formula, control_settings = control.ergm()) {
  # Input validation
  if (!inherits(network, "network")) {
    stop("network must be a valid network object.")
  }
  if (!inherits(formula, "formula")) {
    stop("formula must be a valid ergm formula.")
  }
  
  # Rewrite the formula with the network
  current_formula <- rewrite_formula(formula, network)
  
  # Estimate the ERGM model using tryCatch to capture errors and warnings
  model <- tryCatch(
    {
      ergm::ergm(current_formula, control = control_settings)
    },
    error = function(e) {
      message("Model estimation failed with the following error:")
      message(e$message)
      return(NULL)
    },
    warning = function(w) {
      message("Model estimation generated the following warning:")
      message(w$message)
      return(NULL)
    }
  )
  
  # If the model was estimated, check for convergence and infinite coefficients
  if (!is.null(model)) {
    # Check for convergence issues
    if (isTRUE(model$failure)) {
      message("Model did not converge; skipping network.")
      return(NULL)
    }
    
    # Check if any estimated coefficients are infinite
    # Using ergm::coef() to extract the coefficients
    coefs <- tryCatch(ergm::coef(model), error = function(e) NULL)
    if (!is.null(coefs) && any(is.infinite(coefs))) {
      message("Model estimated but contains infinite coefficient(s); skipping network.")
      return(NULL)
    }
    
    # If everything is fine, return the model
    return(model)
  }
  
  return(NULL)
}


#' Function for sequential estimation of network models
#'
#' This function takes a list of network objects, a formula specifying the model, and optional control settings.
#' It performs estimation of the specified model on each network in the list using the ergm_estimation function.
#' The function checks for any networks that could not be estimated and removes them from the output.
#' It also checks for any networks with infinite estimates and removes them from the output.
#'
#' @param network_list A list of network objects.
#' @param formula A formula specifying the model to be estimated.
#' @param control_settings Optional control settings for the estimation process. Default is control.ergm().
#' @param verbose Logical. If TRUE, prints the index of each network being estimated.
#' @returns A list containing the estimated network models, removed networks, and estimation information.
#'
#' @examples
#' # networks <- list(net1, net2, net3)
#' # result <- sequential_estimation(networks, ~edges + nodematch("gender"), verbose = TRUE)
#' @export
sequential_estimation <- function(network_list, formula, control_settings = control.ergm(), verbose = TRUE) {
  # Input validation
  if (!is.list(network_list) || any(sapply(network_list, function(x) !inherits(x, "network")))) {
    stop("network_list must be a list of network objects.")
  }
  if (!inherits(formula, "formula")) {
    stop("formula must be a valid ergm formula.")
  }
  
  n_networks <- length(network_list)
  estimated_models <- vector("list", n_networks)
  removed_networks <- c()
  
  for (i in seq_along(network_list)) {
    if (verbose) cat("Estimating network", i, "of", n_networks, "\n")
    model <- ergm_estimation(network_list[[i]], formula, control_settings)
    if (is.null(model)) {
      removed_networks <- c(removed_networks, i)
    } else {
      estimated_models[[i]] <- model
    }
  }
  
  # Filter out NULL elements (networks that failed estimation)
  valid_idx <- which(!sapply(estimated_models, is.null))
  estimated_models <- estimated_models[valid_idx]
  
  # Check for networks with infinite estimates or no variance
  infinite_estimates <- sapply(estimated_models, function(x) any(is.infinite(summary(x)$coefficients[, 1])))
  no_variance <- sapply(estimated_models, function(x) any(is.na(summary(x)$coefficients[, 2])))
  remove_indices <- infinite_estimates | no_variance
  
  if (any(remove_indices)) {
    removed_idx <- valid_idx[which(remove_indices)]
    removed_networks <- unique(c(removed_networks, removed_idx))
    estimated_models <- estimated_models[!remove_indices]
    cat("Networks", paste0(removed_idx, collapse = " "), 
        "contain infinite estimates or no variance and were removed.\n")
    cat("Consider re-parametrizing the networks if you want to keep more networks in.\n")
  }
  
  # Prepare the output object
  result <- list(
    estimated_models = estimated_models,
    removed_networks = removed_networks,
    estimation_info = list(
      total_networks = n_networks,
      successfully_estimated = length(estimated_models),
      removed_infinite_estimates = sum(infinite_estimates),
      removed_no_variance = sum(no_variance)
    )
  )
  return(result)
}