# ERGMeta

ERGMeta is an R package for meta-analysis of Exponential Random Graph Models (ERGMs) across multiple networks. It provides functions for estimating ERGMs on a list of networks, extracting the estimates and standard errors, and performing meta-analysis using Bayesian regression using the [`brms` (Bürkner, 2017)](https://github.com/paul-buerkner/brms) package. In future releases, the plan is to incorporate visualization features for enhanced model interpretation.

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
install.packages("cmdstanr")
```


## Functions

### `sequential_estimation`

Function for sequential estimation of network models.

- `network_list`: A list of network objects.
- `formula`: A formula specifying the model to be estimated.
- `control_settings`: Optional control settings for the estimation process. Default is `NULL`.

Returns a list containing the estimated network models, removed networks, and estimation information.

### `extract_estimates`

Function to extract estimates and standard errors from a list of fitted models.

- `fit_list`: A list of fitted models (e.g., `ergm` objects).
- `grouping_variable`: (Optional) The variable used for grouping the network objects.

Returns a tibble containing the estimates and the standard errors for each coefficient.

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

## Examples

Here are some examples of how to use the functions in the ERGMeta package:

```R
# Load the required packages
library(ERGMeta)
library(ergm)

# Create a list of network objects
network_list <- list(
  network1 = ergm::network(matrix(rbinom(100, 1, 0.3), nrow = 10)),
  network2 = ergm::network(matrix(rbinom(100, 1, 0.4), nrow = 10)),
  network3 = ergm::network(matrix(rbinom(100, 1, 0.5), nrow = 10))
)

# Specify the ERGM formula
formula <- ~ edges + triangle

# Perform sequential estimation on the network list
estimated_models <- sequential_estimation(network_list, formula)

# Extract estimates and standard errors from the fitted models
estimates_df <- extract_estimates(estimated_models)

# Fit a meta-analysis model using Bayesian regression
meta_analysis_model <- meta_fit(estimates_df)

# Print the summary of the meta-analysis model
summary(meta_analysis_model)
```

In this example, we first create a list of network objects using the `ergm::network` function. We then specify the ERGM formula using the `edges` and `triangle` terms.

Next, we use the `sequential_estimation` function to perform sequential estimation of the ERGM on each network in the list. The `extract_estimates` function is then used to extract the estimates and standard errors from the fitted models.

Finally, we use the `meta_fit` function to fit a meta-analysis model using Bayesian regression on the extracted estimates and standard errors. The summary of the meta-analysis model can be obtained using the `summary` function.

## Contributing

Contributions to the ERGMeta package are welcome! If you find any bugs, have suggestions for improvements, or would like to add new features, please open an issue or submit a pull request on the GitHub repository.

## License

The ERGMeta package is licensed under the MIT License.

## References

Tolochko, P., & Boomgaarden, H. G. (2024). Same but different: A comparison of estimation approaches for exponential random graph models for multiple networks. Social Networks, 76, 1-11.