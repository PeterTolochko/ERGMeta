% Manual page for plot_simulation_summary
\name{plot_simulation_summary}
\alias{plot_simulation_summary}
\title{Plot Simulation Benchmark Summary}
\usage{
plot_simulation_summary(
  x,
  metric = c("bias", "rmse", "coverage_95", "mean_tau2", "convergence_rate",
    "error_rate")
)
}
\arguments{
\item{x}{Either an object returned by \code{run_ergmeta_simulation()} or a
compatible summary data frame.}

\item{metric}{Metric to plot from the simulation summary.}
}
\value{
A \code{ggplot2} object.
}
\description{
Create a compact method-comparison plot from simulation benchmark summaries.
}
\examples{
if (requireNamespace("ggplot2", quietly = TRUE)) {
  bench <- run_ergmeta_simulation(
    n_sims = 3,
    n_networks = 14,
    methods = c("sens_random", "cov_random"),
    progress = FALSE,
    seed = 303
  )
  p <- plot_simulation_summary(bench, metric = "rmse")
  p
}
}
