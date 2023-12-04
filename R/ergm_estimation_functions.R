#' Function to rewrite a formula by replacing the response variable with a specified network
#' This is mostly for internal use.
#' 
#' @param formula The original formula to be rewritten. Must be `ergm`-formula.
#' @param network The network to be used as the new response variable.
#' @returns The rewritten formula with the network as the response variable
#' @export
rewrite_formula <- function(formula, network) {
  network_name <- deparse(substitute(network))
  terms <- toString(formula[-c(1:2)]) # everything but "network_name ~"
  # Construct the new formula as a character string
  new_formula <- as.formula(paste(network_name, "~", terms))
  return(new_formula)
}


#' Function to estimate an Exponential Random Graph Model (ERGM) on a given network
#' @param network The input network object.
#' @param formula Formula for ERGM estimation.
#' @param contro_settings: Optional control settings for the ERGM estimation.
#' @returns The estimated ERGM model object.
#' @export
ergm_estimation <- function(network,
                            formula,
                            control_settings = NULL) {
  current_formula <- rewrite_formula(formula, network)
  out <- tryCatch(
    {
      if (is.null(control_settings)) {
        ergm::ergm(current_formula)
      } else {
        ergm::ergm(current_formula, control = control_settings)
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


#' Function for sequential estimation of network models
#'
#' This function takes a list of network objects, a formula specifying the model, and optional control settings.
#' It performs estimation of the specified model on each network in the list using the ergm_estimation function.
#' The function checks for any networks that could not be estimated and removes them from the output.
#' It also checks for any networks with infinite estimates and removes them from the output.
#'
#' @param network_list A list of network objects.
#' @param formula A formula specifying the model to be estimated
#' @param control_settings Optional control settings for the estimation process.
#' @returns A list of estimated network models.
#'
#' @examples
#' # networks <- list(net1, net2, net3)
#' # models <- sequential_estimation(networks, ~edges + nodematch("gender"))
#' @export
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
    cat("\nConsider re-parametrizing the networks if you want to keep more networks in.")
  }
  return(out_models)
}