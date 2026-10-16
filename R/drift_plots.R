library(data.table)
library(ggplot2)
source("R/palette.R")

#' Panel plot of obs/clearsky ratio over time, one panel per lake
#'
#' @param drift_dt data.table from run_all_lakes() with date, ratio, lake, clarity_score
#' @return ggplot object
plot_drift_per_lake <- function(drift_dt, min_smooth = 10) {
  # Only add smoother for lakes with enough points
  lakes_ok <- drift_dt[, .N, by = lake][N >= min_smooth, lake]
  smooth_dt <- drift_dt[lake %in% lakes_ok]

  ggplot(drift_dt, aes(x = date, y = ratio)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_point(aes(colour = clarity_score), size = 1.5, alpha = 0.8) +
    geom_smooth(data = smooth_dt, method = "loess", se = FALSE,
                colour = "firebrick", linewidth = 0.7, span = 0.75) +
    scale_colour_viridis_c(name = "Clarity\nscore", direction = -1) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    facet_wrap(~lake, scales = "free_x", ncol = 4) +
    labs(
      x = NULL, y = "Observed / Clear-sky ratio",
      title = "Solar radiation sensor drift",
      subtitle = "Ratio of observed to FAO56 clear-sky radiation on clearest days"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
    )
}

#' Multi-lake summary line plot of drift ratio over time
#'
#' @param drift_dt data.table with date, ratio, lake
#' @return ggplot object
plot_drift_summary <- function(drift_dt, min_smooth = 10) {
  dt <- drift_dt[lake %in% drift_dt[, .N, by = lake][N >= min_smooth, lake]]
  ggplot(dt, aes(x = date, y = ratio, colour = lake)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.9, span = 0.75) +
    scale_colour_limno() +
    labs(
      x = NULL, y = "Observed / Clear-sky ratio",
      title = "Solar radiation sensor drift — all lakes",
      colour = "Lake"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

#' Detailed drift timeseries for a single lake
#'
#' Shows observed mean vs clearsky mean on clear days, plus the ratio.
#'
#' @param drift_dt data.table with date, obs_mean, clearsky_mean, ratio, lake
#' @param lake_name character, lake to plot
#' @return ggplot object
plot_drift_timeseries <- function(drift_dt, lake_name = NULL) {
  if (!is.null(lake_name)) {
    dt <- drift_dt[lake == lake_name]
  } else {
    dt <- copy(drift_dt)
  }

  dt_long <- melt(dt, id.vars = c("date", "lake"),
                  measure.vars = c("obs_mean", "clearsky_mean"),
                  variable.name = "source", value.name = "wm2")
  dt_long[, source := fifelse(source == "obs_mean", "Observed", "Clear-sky")]

  ttl <- if (!is.null(lake_name)) {
    paste(lake_name, "— observed vs clear-sky radiation")
  } else {
    "Observed vs clear-sky radiation — all lakes"
  }

  ggplot(dt_long, aes(x = date, y = wm2, colour = source)) +
    geom_point(size = 1, alpha = 0.6) +
    geom_line(alpha = 0.3, linewidth = 0.4) +
    scale_colour_manual(values = c("Observed" = "steelblue", "Clear-sky" = "orange")) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    facet_wrap(~lake, ncol = 2, scales = "free_x") +
    labs(
      x = NULL, y = expression(paste("Mean radiation (W/", m^2, ")")),
      title = ttl,
      subtitle = "Clear days only, \u00b12 hrs around solar noon",
      colour = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 9),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      panel.spacing = unit(0.3, "lines"),
      legend.position = "bottom"
    )
}

#' Diagnostic plot: single clear day observed vs clearsky curve
#'
#' @param met_dt data.table with DateTime and RadSWD
#' @param clearsky_dt data.table with DateTime and clearsky_wm2
#' @param target_date Date, the day to plot
#' @param solar_noon_dt data.table with date and solar_noon
#' @param hours_window numeric, hours either side of solar noon
#' @return ggplot object
plot_clearday_example <- function(met_dt, clearsky_dt, target_date,
                                   solar_noon_dt, hours_window = 3) {
  target_date <- as.Date(target_date)
  noon <- solar_noon_dt[date == target_date, solar_noon]

  obs <- met_dt[as.Date(DateTime) == target_date]
  cs <- clearsky_dt[as.Date(DateTime) == target_date]

  merged <- merge(obs[, .(DateTime, RadSWD)], cs, by = "DateTime", all = TRUE)
  merged[, hrs_from_noon := as.numeric(difftime(DateTime, noon, units = "hours"))]
  merged <- merged[abs(hrs_from_noon) <= hours_window]

  dt_long <- melt(merged, id.vars = "DateTime",
                  measure.vars = c("RadSWD", "clearsky_wm2"),
                  variable.name = "source", value.name = "wm2")
  dt_long[, source := fifelse(source == "RadSWD", "Observed", "Clear-sky")]

  ggplot(dt_long, aes(x = DateTime, y = wm2, colour = source)) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = noon, linetype = "dotted", colour = "grey40") +
    scale_colour_manual(values = c("Observed" = "steelblue", "Clear-sky" = "orange")) +
    labs(
      x = NULL, y = expression(paste("Radiation (W/", m^2, ")")),
      title = paste("Clear day example:", target_date),
      colour = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}
