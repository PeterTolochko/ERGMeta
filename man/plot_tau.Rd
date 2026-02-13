% Manual page for plot_tau
\name{plot_tau}
\alias{plot_tau}
\title{Plot Heterogeneity (tau / tau\eqn{^2})}
\usage{
plot_tau(x, scale = c("tau2", "tau"), include_baseline = TRUE)
}
\arguments{
\item{x}{Either a data frame from \code{run_exclusion_sensitivity()}
containing scenario-level \code{tau}/\code{tau2}, or an object returned by
\code{meta_fit_covariance()} (or \code{methodology_roadmap()} \code{$covariance_meta}).}

\item{scale}{Quantity to plot: \code{"tau2"} (default) or \code{"tau"}.}

\item{include_baseline}{Logical used for scenario-level sensitivity input.
If \code{FALSE}, drops \code{"all"}.}
}
\value{
A \code{ggplot2} object.
}
\description{
Visualize heterogeneity estimates either by scenario (sensitivity output) or by
coefficient (covariance-aware fit).
}
\examples{
if (requireNamespace("ggplot2", quietly = TRUE)) {
  sim <- simulate_ergmeta_data(n_networks = 14, seed = 202)
  sens <- run_exclusion_sensitivity(
    df = sim$estimates,
    method = "random",
    include_leave_one_out = FALSE
  )
  p1 <- plot_tau(sens, scale = "tau2")
  p1

  fit <- meta_fit_covariance(
    estimates_df = sim$estimates,
    covariance_df = sim$covariance,
    method = "random"
  )
  p2 <- plot_tau(fit, scale = "tau")
  p2
}
}
