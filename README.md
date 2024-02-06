# ERGMeta

ERGMeta is an R package designed to streamline the estimation of multiple Exponential Random Graph Models (ERGMs). The package provides essential functions for sequentially estimating multiple networks, extracting network coefficients, and conducting hierarchical Bayesian meta-analysis models using the [`brms` (Bürkner, 2017)](https://github.com/paul-buerkner/brms) package. In future releases, the plan is to incorporate visualization features for enhanced model interpretation.

Mostly based on the [Same but different: A comparison of estimation approaches for exponential random graph models for multiple networks](https://www.sciencedirect.com/science/article/pii/S0378873323000357#aep-article-footnote-id1). The goal is to extend and simplify workflows that deal with estimation of multiple exponential random graph models.


## Intallation
To use ERGMeta, you'll need to have [`brms`](https://github.com/paul-buerkner/brms) installed. Additionally, for optimal performance, it is highly recommend installing [`cmdstanr`](https://mc-stan.org/cmdstanr/) as it serves as the backend for brms. Please refer to the respective links for detailed installation instructions.

```r
# Install brms
install.packages("brms")

# Install cmdstanr (optional but recommended)
install.packages("cmdstanr")
```
After installing the dependencies, you can install `ERGMeta` from GitHub using the following:

```r
install.packages("devtools")
devtools::install_github("PeterTolochko/ERGMeta")
```

# Usage

To use `ERGMeta` for estimating multiple ERGMs, extracting coefficients, and conducting hierarchical Bayesian meta-analysis, check out these examples (TODO!)

