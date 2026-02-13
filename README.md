# ERGMeta

ERGMeta is an R package for meta-analysis of Exponential Random Graph Models (ERGMs) across multiple networks. It provides functions for estimating ERGMs on a list of networks, extracting the estimates and standard errors, fitting meta-analytic models, and visualizing key diagnostics from sensitivity and simulation workflows.

This package is based on the paper: [Tolochko, P., & Boomgaarden, H. G. (2024). Same but different: A comparison of estimation approaches for exponential random graph models for multiple networks. Social Networks, 76, 1-11](https://www.sciencedirect.com/science/article/pii/S0378873323000357).

## Installation

You can install the development version of ERGMeta from GitHub using the `devtools` package:

```R
# Install devtools if not already installed
# install.packages("devtools")

# Install ERGMeta from GitHub
devtools::install_github("PeterTolochko/ERGMeta")
```

To use ERGMeta, you'll need to have [`brms`](https://github.com/paul-buerkner/brms) installed. Additionally, for optimal performance, it is highly recommend installing [`cmdstanr`](https://mc-stan.org/cmdstanr/) as it serves as the backend for brms. Please refer to the respective links for detailed installation instructions.

```r
# Install brms
install.packages("brms")

# Install cmdstanr (optional but recommended)
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

# Install CmdStan toolchain (one-time)
cmdstanr::install_cmdstan()

# Check installation
cmdstanr::cmdstan_version()
```


## Functions

### `sequential_estimation`

Function for sequential estimation of network models.

- `network_list`: A list of network objects.
- `formula`: A formula specifying the model to be estimated.
- `control_settings`: Optional control settings for the estimation process. Default is `ergm::control.ergm()`.

Returns a list containing the estimated network models, removed networks, and estimation information.

### `extract_estimates`

Function to extract estimates and standard errors from a list of fitted models.

- `fit_list`: A list of fitted models (e.g., `ergm` objects).
- `grouping_variable`: (Optional) A network-level attribute used for grouping.

Returns a tibble containing: network id, coefficient id (`coefs`), original coefficient term (`term`), estimate, and standard error.

### `meta_fit`

Function to fit a meta-analysis model using Bayesian regression.

- `df`: A data frame containing the data for the meta-analysis.
- `group_vars`: (optional) A character vector indicating the grouping variables.
- `chains`: The number of Markov chains to run (default: 4).
- `cores`: The number of CPU cores to use for parallel computing (default: 4).
- `iter`: The number of iterations for each chain (default: 4000).
- `priors`: The prior distribution for the regression coefficients (default: `normal(0, 2)`).
- `backend`: The backend for running the Bayesian regression (default: "cmdstanr").
- `sigma`: (Optional) Logical indicating whether to include sigma in the formula (default: FALSE).

Returns the fitted model object.

### Methodological roadmap helpers

- `extract_within_covariance`: extracts per-network ERGM covariance matrices in tidy form.
- `prepare_meta_inputs`: creates estimates, covariance, and inclusion diagnostics in one object.
- `meta_fit_covariance`: fits covariance-aware multivariate meta-regression using full within-network covariance.
- `run_exclusion_sensitivity`: computes fixed/random pooled estimates, including heterogeneity (`tau`), across baseline and exclusion scenarios.
- `methodology_roadmap`: runs the full workflow end-to-end with optional covariance-aware fitting.

### Simulation helpers

- `simulate_ergmeta_data`: generates synthetic estimates/covariances/truth in ERGMeta format.
- `run_ergmeta_simulation`: runs a Monte Carlo benchmark and returns bias, RMSE, coverage, and heterogeneity summaries.

### Visualization helpers

- `plot_sensitivity`: plots scenario-level changes in pooled effects and heterogeneity.
- `plot_tau`: visualizes heterogeneity (`tau` or `tau^2`) from sensitivity runs or covariance-aware fits.
- `plot_simulation_summary`: compares simulation benchmark metrics across methods.

## Vignettes

- End-to-end workflow vignette:
  `vignette("end-to-end-workflow", package = "ERGMeta")`

## Examples

Here are some examples of how to use the functions in the ERGMeta package:

```R
# Load the required packages
library(ERGMeta)
library(network)

# Create a list of network objects
network_list <- list(
  network(matrix(rbinom(100, 1, 0.3), nrow = 10), directed = TRUE),
  network(matrix(rbinom(100, 1, 0.4), nrow = 10), directed = TRUE),
  network(matrix(rbinom(100, 1, 0.5), nrow = 10), directed = TRUE)
)

# Specify ERGM terms.
# The response on the left-hand side is a placeholder and is replaced internally
# with each network in network_list.
formula <- net_placeholder ~ edges + mutual

# Perform sequential estimation on the network list
estimated_models <- sequential_estimation(network_list, formula)

# Extract estimates and standard errors from the fitted models
estimates_df <- extract_estimates(estimated_models$estimated_models)

# Fit a meta-analysis model using Bayesian regression
meta_analysis_model <- meta_fit(estimates_df)

# Print the summary of the meta-analysis model
summary(meta_analysis_model)

```

In this example, we first create a list of network objects using the `network::network` function. We then specify the ERGM formula using the `edges` and `mutual` terms.

Next, we use the `sequential_estimation` function to perform sequential estimation of the ERGM on each network in the list. The `extract_estimates` function is then used to extract the estimates and standard errors from the fitted models.

Finally, we use the `meta_fit` function to fit a meta-analysis model using Bayesian regression on the extracted estimates and standard errors. The summary of the meta-analysis model can be obtained using the `summary` function.

### Minimal working pipeline (fast smoke run)

```R
library(ERGMeta)
library(network)
library(ergm)

set.seed(123)

network_list <- list(
  network(matrix(rbinom(64, 1, 0.25), nrow = 8), directed = TRUE),
  network(matrix(rbinom(64, 1, 0.30), nrow = 8), directed = TRUE),
  network(matrix(rbinom(64, 1, 0.35), nrow = 8), directed = TRUE)
)

formula <- net_placeholder ~ edges + mutual

# Keep this lightweight for a quick test
ctrl <- control.ergm(
  MCMLE.maxit = 1,
  MCMC.interval = 100,
  MCMC.samplesize = 1000
)

fit_out <- sequential_estimation(network_list, formula, control_settings = ctrl, verbose = FALSE)
estimates_df <- extract_estimates(fit_out$estimated_models)

# Small iteration count for smoke testing
meta_model <- meta_fit(
  estimates_df,
  iter = 1000,
  chains = 2,
  cores = 2,
  backend = "cmdstanr"
)

summary(meta_model)
```

### Concrete roadmap example

```R
# fit_list may contain both successful ergm objects and failed fits (NULL)
roadmap <- methodology_roadmap(
  fit_list = estimated_models$estimated_models,
  grouping_variable = NULL,
  include_leave_one_out = TRUE,
  sensitivity_method = "random",
  fit_covariance_model = TRUE
)

# Which networks were included or excluded?
roadmap$diagnostics

# Term-aligned within-network covariance for downstream multivariate pooling
head(roadmap$covariance)

# Covariance-aware multivariate fit output
roadmap$covariance_meta$reference_summary

# Sensitivity of pooled estimates and tau to network exclusions
roadmap$sensitivity
```

### Simulation benchmark example

```R
sim_benchmark <- run_ergmeta_simulation(
  n_sims = 20,
  n_networks = 30,
  methods = c("sens_fixed", "sens_random", "sens_random_mod", "cov_random", "cov_random_mod"),
  seed = 2026
)

sim_benchmark$summary
```

### Visualization example

```R
if (requireNamespace("ggplot2", quietly = TRUE)) {
  sim_data <- simulate_ergmeta_data(n_networks = 20, seed = 99)
  sens <- run_exclusion_sensitivity(
    df = sim_data$estimates,
    method = "random",
    include_leave_one_out = TRUE
  )
  cov_fit <- meta_fit_covariance(
    estimates_df = sim_data$estimates,
    covariance_df = sim_data$covariance,
    method = "random"
  )
  plot_sensitivity(sens, value = "delta_from_all")
  plot_tau(cov_fit, scale = "tau2")
  plot_simulation_summary(sim_benchmark, metric = "rmse")
}
```

## Runtime Notes

- ERGM estimation and Bayesian meta-fitting can be computationally expensive.
- For quick checks, use small networks, fewer ERGM MCMC iterations, and lower `iter/chains` in `meta_fit`.
- For substantive analysis, increase ERGM control settings and `meta_fit` sampling settings, and run sensitivity checks (`run_exclusion_sensitivity`) plus covariance-aware fitting (`meta_fit_covariance`).

## Contributing

Contributions to the ERGMeta package are welcome! If you find any bugs, have suggestions for improvements, or would like to add new features, please open an issue or submit a pull request.

## License

The ERGMeta package is licensed under the MIT License.

## References

Tolochko, P., & Boomgaarden, H. G. (2024). Same but different: A comparison of estimation approaches for exponential random graph models for multiple networks. Social Networks, 76, 1-11.

Bürkner, P. C. (2017). brms: An R package for Bayesian multilevel models using Stan. Journal of Statistical Software, 80(1), 1-28.
