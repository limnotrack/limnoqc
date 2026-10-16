library(data.table)
library(ggplot2)
source("R/palette.R")

#' Load humidity data for a single lake with monthly summary stats
#'
#' Computes monthly max (proxy for 100% saturation point), median, and 0.01
#' quantile of HumRel. These track sensor drift:
#'   - Max drifting below 100% suggests sensor degradation
#'   - Low quantile drifting up/down indicates bias shift
#'   - Median shift indicates overall calibration drift
#'
#' @param lake_name character
#' @param rds_dir character
#' @param instruments data.table with lake, date_from, date_to, instrument
#' @return data.table with monthly summaries
load_humidity <- function(lake_name, rds_dir = "data/rds",
                           instruments = NULL) {
  path <- file.path(rds_dir, paste0("full_", lake_name, ".rds"))
  x <- readRDS(path)
  met <- as.data.table(x[["raw"]][["met"]])
  met <- met[!is.na(HumRel), .(DateTime, HumRel)]

  if (nrow(met) == 0) return(data.table())

  met[, `:=`(
    date      = as.Date(DateTime),
    year_month = as.Date(format(DateTime, "%Y-%m-01"))
  )]

  # Monthly stats
  monthly <- met[, .(
    rh_max    = max(HumRel),
    rh_median = median(HumRel),
    rh_q01    = quantile(HumRel, 0.01),
    rh_q99    = quantile(HumRel, 0.99),
    pct_100   = 100 * mean(HumRel >= 100),
    n         = .N
  ), by = year_month]

  monthly[, lake := lake_name]

  # Add instrument type if available
  if (!is.null(instruments)) {
    inst <- instruments[lake == lake_name]
    if (nrow(inst) > 0) {
      monthly[, instrument := "unknown"]
      for (i in seq_len(nrow(inst))) {
        monthly[year_month >= inst$date_from[i] & year_month <= inst$date_to[i],
                instrument := inst$instrument[i]]
      }
    } else {
      monthly[, instrument := "unknown"]
    }
  } else {
    monthly[, instrument := "unknown"]
  }

  setorder(monthly, year_month)
  monthly
}

#' Load humidity summaries for all lakes
#'
#' @param rds_dir character
#' @param instruments data.table from detect_instruments()
#' @return data.table with monthly humidity stats for all lakes
run_all_humidity <- function(rds_dir = "data/rds", instruments = NULL) {
  files <- list.files(rds_dir, pattern = "^full_.*[.]rds$", full.names = FALSE)
  lake_names <- sub("^full_", "", sub("[.]rds$", "", files))

  rbindlist(lapply(lake_names, function(lk) {
    cat("Loading", lk, "...\n")
    tryCatch(load_humidity(lk, rds_dir, instruments),
             error = function(e) {
               cat("  ERROR:", conditionMessage(e), "\n")
               data.table()
             })
  }), fill = TRUE)
}

#' Panel plot of RH drift metrics over time, one panel per lake
#'
#' Shows monthly max (100% point), median, and q01 as lines.
#'
#' @param hum_dt data.table from run_all_humidity()
#' @return ggplot object
plot_humidity_drift <- function(hum_dt) {
  dt <- melt(hum_dt, id.vars = c("lake", "year_month", "instrument"),
             measure.vars = c("rh_max", "rh_median", "rh_q01"),
             variable.name = "metric", value.name = "value")
  dt[, metric := fcase(
    metric == "rh_max",    "Max (100% point)",
    metric == "rh_median", "Median",
    metric == "rh_q01",    "Q0.01"
  )]
  dt[, metric := factor(metric, levels = c("Max (100% point)", "Median", "Q0.01"))]

  ggplot(dt, aes(x = year_month, y = value, colour = metric)) +
    geom_line(linewidth = 0.5, alpha = 0.8) +
    geom_hline(yintercept = 100, linetype = "dashed", colour = "grey50",
               linewidth = 0.3) +
    scale_colour_manual(values = c("Max (100% point)" = "firebrick",
                                    "Median" = "steelblue",
                                    "Q0.01" = "forestgreen")) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    facet_wrap(~lake, ncol = 4, scales = "free_x") +
    labs(
      x = NULL, y = "Relative humidity (%)",
      title = "Relative humidity sensor drift",
      subtitle = "Monthly max, median, and 1st percentile over time",
      colour = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom"
    )
}

#' Summary plot comparing Lufft vs Vaisala RH distributions
#'
#' Boxplots of max, median, and q01 grouped by instrument type.
#'
#' @param hum_dt data.table from run_all_humidity()
#' @return ggplot object
plot_humidity_by_instrument <- function(hum_dt) {
  dt <- hum_dt[instrument %in% c("Lufft", "Vaisala")]
  dt_long <- melt(dt, id.vars = c("lake", "year_month", "instrument"),
                  measure.vars = c("rh_max", "rh_median", "rh_q01"),
                  variable.name = "metric", value.name = "value")
  dt_long[, metric := fcase(
    metric == "rh_max",    "Max (100% point)",
    metric == "rh_median", "Median",
    metric == "rh_q01",    "Q0.01"
  )]

  ggplot(dt_long, aes(x = instrument, y = value, fill = instrument)) +
    geom_boxplot(alpha = 0.6, outlier.size = 0.5) +
    scale_fill_manual(values = c("Lufft" = "steelblue", "Vaisala" = "coral")) +
    facet_wrap(~metric, scales = "free_y", ncol = 3) +
    labs(
      x = NULL, y = "Relative humidity (%)",
      title = "RH sensor characteristics: Lufft vs Vaisala",
      subtitle = "Monthly statistics across all lakes",
      fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}

#' Multi-lake line plot of monthly max RH (100% drift indicator)
#'
#' @param hum_dt data.table from run_all_humidity()
#' @return ggplot object
plot_humidity_max_summary <- function(hum_dt) {
  ggplot(hum_dt, aes(x = year_month, y = rh_max, colour = lake)) +
    geom_hline(yintercept = 100, linetype = "dashed", colour = "grey50") +
    geom_line(linewidth = 0.5, alpha = 0.7) +
    scale_colour_limno() +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(
      x = NULL, y = "Monthly max RH (%)",
      title = "RH sensor 100% saturation point drift",
      subtitle = "Monthly maximum — should reach 100% regularly",
      colour = "Lake"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}
