#!/usr/bin/env Rscript

# Example simulation benchmark for ERGMeta package methods.

library(ERGMeta)

set.seed(2026)

benchmark <- run_ergmeta_simulation(
  n_sims = 20,
  n_networks = 30,
  coef_terms = c("edges", "mutual"),
  beta_intercept = c(-2.2, 0.4),
  beta_moderator = c(0.35, -0.15),
  tau = c(0.25, 0.15),
  tau_cor = 0.1,
  within_sd = c(0.25, 0.20),
  within_cor = 0.25,
  methods = c("sens_fixed", "sens_random", "sens_random_mod", "cov_random", "cov_random_mod"),
  progress = TRUE
)

print(benchmark$summary)
