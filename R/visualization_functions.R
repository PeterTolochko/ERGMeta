#' Lightweight Visualization Helpers for ERGMeta Workflows
#'
#' The plotting helpers in this file are intentionally lightweight wrappers
#' around key ERGMeta outputs. They use \pkg{ggplot2} only when called.
#'
#' @keywords internal
.require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for plotting. Install it with install.packages('ggplot2').",
      call. = FALSE
    )
  }
}


#' Plot Exclusion-Sensitivity Results
#'
#' @param sensitivity_df Either a data frame returned by \code{run_exclusion_sensitivity()}
#'   or an object returned by \code{methodology_roadmap()}.
#' @param value Quantity to visualize: \code{"delta_from_all"}, \code{"pooled_estimate"}, or \code{"tau2"}.
#' @param include_baseline Logical. If \code{FALSE}, drops the \code{"all"} scenario from the plot.
#'
#' @returns A \code{ggplot2} object.
#' @export
plot_sensitivity <- function(sensitivity_df,
                             value = c("delta_from_all", "pooled_estimate", "tau2"),
                             include_baseline = TRUE) {
  .require_ggplot2()
  value <- match.arg(value)

  if (is.list(sensitivity_df) &&
      !is.data.frame(sensitivity_df) &&
      !is.null(sensitivity_df$sensitivity) &&
      is.data.frame(sensitivity_df$sensitivity)) {
    sensitivity_df <- sensitivity_df$sensitivity
  }

  if (!is.data.frame(sensitivity_df)) {
    stop("sensitivity_df must be a sensitivity data frame or a roadmap object with $sensitivity.")
  }
  if (!is.logical(include_baseline) || length(include_baseline) != 1 || is.na(include_baseline)) {
    stop("include_baseline must be a single TRUE/FALSE value.")
  }

  required_cols <- c("scenario", "coefs", "term", value)
  missing_cols <- setdiff(required_cols, colnames(sensitivity_df))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "sensitivity_df is missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  plot_df <- sensitivity_df[, required_cols, drop = FALSE]
  if (!isTRUE(include_baseline)) {
    plot_df <- plot_df[plot_df$scenario != "all", , drop = FALSE]
  }
  if (nrow(plot_df) == 0) {
    stop("No rows available to plot after applying filters.")
  }

  plot_df$metric_value <- as.numeric(plot_df[[value]])
  if (!any(is.finite(plot_df$metric_value))) {
    stop(sprintf("Column '%s' does not contain finite values to plot.", value))
  }

  scenario_levels <- unique(as.character(plot_df$scenario))
  if ("all" %in% scenario_levels) {
    scenario_levels <- c("all", setdiff(scenario_levels, "all"))
  }
  plot_df$scenario <- factor(plot_df$scenario, levels = scenario_levels)

  value_label <- switch(
    value,
    delta_from_all = "Change From Baseline",
    pooled_estimate = "Pooled Estimate",
    tau2 = "Heterogeneity (tau^2)"
  )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$scenario,
      y = .data$metric_value,
      color = .data$coefs,
      group = .data$coefs
    )
  ) +
    ggplot2::geom_point(na.rm = TRUE, size = 1.8, alpha = 0.9) +
    ggplot2::geom_line(na.rm = TRUE, alpha = 0.75) +
    ggplot2::facet_wrap(~term, scales = "free_y") +
    ggplot2::labs(
      x = "Scenario",
      y = value_label,
      color = "Coefficient",
      title = "Exclusion Sensitivity"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (value == "delta_from_all") {
    p <- p + ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35)
  }

  return(p)
}


#' Plot Heterogeneity (\code{tau} / \code{tau^2})
#'
#' @param x Either:
#'   \itemize{
#'   \item a data frame from \code{run_exclusion_sensitivity()} containing scenario-level \code{tau}/\code{tau2}, or
#'   \item an object returned by \code{meta_fit_covariance()} (or \code{methodology_roadmap()} \code{$covariance_meta}).
#'   }
#' @param scale Quantity to plot: \code{"tau2"} (default) or \code{"tau"}.
#' @param include_baseline Logical used for scenario-level sensitivity input. If \code{FALSE}, drops \code{"all"}.
#'
#' @returns A \code{ggplot2} object.
#' @export
plot_tau <- function(x,
                     scale = c("tau2", "tau"),
                     include_baseline = TRUE) {
  .require_ggplot2()
  scale <- match.arg(scale)

  if (!is.logical(include_baseline) || length(include_baseline) != 1 || is.na(include_baseline)) {
    stop("include_baseline must be a single TRUE/FALSE value.")
  }

  if (is.list(x) &&
      !is.data.frame(x) &&
      !is.null(x$covariance_meta) &&
      is.list(x$covariance_meta) &&
      !is.null(x$covariance_meta$tau) &&
      is.data.frame(x$covariance_meta$tau)) {
    tau_df <- x$covariance_meta$tau
  } else if (is.list(x) && !is.data.frame(x) && !is.null(x$tau) && is.data.frame(x$tau)) {
    tau_df <- x$tau
  } else if (is.data.frame(x)) {
    tau_df <- x
  } else {
    stop("x must be a sensitivity data frame or a covariance-fit object with a $tau data frame.")
  }

  required_cols <- c("coefs", "term", scale)
  missing_cols <- setdiff(required_cols, colnames(tau_df))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Input is missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  if ("scenario" %in% colnames(tau_df)) {
    plot_df <- tau_df[, c("scenario", "coefs", "term", scale), drop = FALSE]
    if (!isTRUE(include_baseline)) {
      plot_df <- plot_df[plot_df$scenario != "all", , drop = FALSE]
    }
    if (nrow(plot_df) == 0) {
      stop("No rows available to plot after applying filters.")
    }
    plot_df$metric_value <- as.numeric(plot_df[[scale]])
    if (!any(is.finite(plot_df$metric_value))) {
      stop(sprintf("Column '%s' does not contain finite values to plot.", scale))
    }

    scenario_levels <- unique(as.character(plot_df$scenario))
    if ("all" %in% scenario_levels) {
      scenario_levels <- c("all", setdiff(scenario_levels, "all"))
    }
    plot_df$scenario <- factor(plot_df$scenario, levels = scenario_levels)

    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x = .data$scenario,
        y = .data$metric_value,
        color = .data$coefs,
        group = .data$coefs
      )
    ) +
      ggplot2::geom_point(na.rm = TRUE, size = 1.8, alpha = 0.9) +
      ggplot2::geom_line(na.rm = TRUE, alpha = 0.75) +
      ggplot2::facet_wrap(~term, scales = "free_y") +
      ggplot2::labs(
        x = "Scenario",
        y = if (scale == "tau2") "tau^2" else "tau",
        color = "Coefficient",
        title = "Heterogeneity Across Scenarios"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid.minor = ggplot2::element_blank()
      )
    return(p)
  }

  plot_df <- tau_df[, c("coefs", "term", scale), drop = FALSE]
  plot_df$metric_value <- as.numeric(plot_df[[scale]])
  if (!any(is.finite(plot_df$metric_value))) {
    stop(sprintf("Column '%s' does not contain finite values to plot.", scale))
  }

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$coefs, y = .data$metric_value, fill = .data$coefs)
  ) +
    ggplot2::geom_col(width = 0.7, alpha = 0.85, show.legend = FALSE, na.rm = TRUE) +
    ggplot2::facet_wrap(~term, scales = "free_y") +
    ggplot2::labs(
      x = "Coefficient",
      y = if (scale == "tau2") "tau^2" else "tau",
      title = "Estimated Heterogeneity"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  return(p)
}


#' Plot Simulation Benchmark Summary
#'
#' @param x Either an object returned by \code{run_ergmeta_simulation()} or a
#'   compatible summary data frame.
#' @param metric Metric to plot from the simulation summary.
#'
#' @returns A \code{ggplot2} object.
#' @export
plot_simulation_summary <- function(x,
                                    metric = c(
                                      "bias",
                                      "rmse",
                                      "coverage_95",
                                      "mean_tau2",
                                      "convergence_rate",
                                      "error_rate"
                                    )) {
  .require_ggplot2()
  metric <- match.arg(metric)

  if (is.list(x) && !is.data.frame(x) && !is.null(x$summary) && is.data.frame(x$summary)) {
    summary_df <- x$summary
  } else if (is.data.frame(x)) {
    summary_df <- x
  } else {
    stop("x must be a simulation benchmark object or a compatible summary data frame.")
  }

  required_cols <- c("method", "coefs", "term", metric)
  missing_cols <- setdiff(required_cols, colnames(summary_df))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Simulation summary is missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  plot_df <- summary_df[, required_cols, drop = FALSE]
  plot_df$metric_value <- as.numeric(plot_df[[metric]])
  if (!any(is.finite(plot_df$metric_value))) {
    stop(sprintf("Column '%s' does not contain finite values to plot.", metric))
  }

  y_label <- switch(
    metric,
    bias = "Bias",
    rmse = "RMSE",
    coverage_95 = "Coverage (95%)",
    mean_tau2 = "Mean tau^2",
    convergence_rate = "Convergence Rate",
    error_rate = "Error Rate"
  )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$method, y = .data$metric_value, fill = .data$method)
  ) +
    ggplot2::geom_col(width = 0.72, alpha = 0.85, show.legend = FALSE, na.rm = TRUE) +
    ggplot2::facet_wrap(~term, scales = "free_y") +
    ggplot2::labs(
      x = "Method",
      y = y_label,
      title = "Simulation Benchmark Summary"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (metric == "bias") {
    p <- p + ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35)
  } else if (metric == "coverage_95") {
    p <- p + ggplot2::geom_hline(yintercept = 0.95, linetype = 2, linewidth = 0.35)
  }

  return(p)
}
