% Manual page for plot_sensitivity
\name{plot_sensitivity}
\alias{plot_sensitivity}
\title{Plot Exclusion-Sensitivity Results}
\usage{
plot_sensitivity(
  sensitivity_df,
  value = c("delta_from_all", "pooled_estimate", "tau2"),
  include_baseline = TRUE
)
}
\arguments{
\item{sensitivity_df}{Either a data frame returned by
\code{run_exclusion_sensitivity()} or an object returned by
\code{methodology_roadmap()}.}

\item{value}{Quantity to visualize: \code{"delta_from_all"},
\code{"pooled_estimate"}, or \code{"tau2"}.}

\item{include_baseline}{Logical. If \code{FALSE}, drops the \code{"all"}
scenario from the plot.}
}
\value{
A \code{ggplot2} object.
}
\description{
Convenience plotting helper for exclusion-sensitivity outputs.
}
\examples{
if (requireNamespace("ggplot2", quietly = TRUE)) {
  sim <- simulate_ergmeta_data(n_networks = 12, seed = 101)
  sens <- run_exclusion_sensitivity(
    df = sim$estimates,
    method = "random",
    include_leave_one_out = TRUE
  )
  p <- plot_sensitivity(sens, value = "delta_from_all")
  p
}
}
