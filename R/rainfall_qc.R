library(data.table)
library(ggplot2)

#' Detect whether rain data is cumulative (Vaisala) or rate (Lufft)
#'
#' Uses run-length encoding of non-decreasing sequences. Cumulative gauges
#' produce long non-decreasing runs (counter only increases between resets).
#' Rate gauges fluctuate constantly — runs are short.
#'
#' A reading is classified as cumulative if it belongs to a non-decreasing run
#' of at least `min_run` readings AND the run's median value > `min_value`.
#'
#' @param rain numeric vector of PpRain values (ordered by time)
#' @param min_run integer, minimum consecutive non-decreasing readings to
#'   classify as cumulative (default 100, ~8 hours at 5-min intervals)
#' @param min_value numeric, minimum median of non-zero PpRain in a run to
#'   qualify (default 5). Filters out flat-zero stretches.
#' @return logical vector, TRUE = cumulative, same length as rain
detect_cumulative <- function(rain, min_run = 100, min_value = 5) {
  n <- length(rain)
  if (n < 2) return(rep(FALSE, n))

  # Is each reading >= the previous? (non-decreasing)
  # Exclude zero-to-zero transitions — flat zeros are not evidence of cumulative
  d <- diff(rain)
  nondec <- c(FALSE, d >= 0 & !(rain[-1] == 0 & rain[-n] == 0))

  # Run-length encode non-decreasing stretches
  rle_nd <- rle(nondec)
  is_cum <- rep(FALSE, n)

  # Map runs back to indices
  end_idx <- cumsum(rle_nd$lengths)
  start_idx <- end_idx - rle_nd$lengths + 1

  for (i in seq_along(rle_nd$values)) {
    if (rle_nd$values[i] && rle_nd$lengths[i] >= min_run) {
      idx <- start_idx[i]:end_idx[i]
      # Check median of non-zero values in this run
      nonzero <- rain[idx][rain[idx] > 0]
      if (length(nonzero) > 0 && median(nonzero) > min_value) {
        is_cum[idx] <- TRUE
      }
    }
  }

  is_cum
}

#' Convert cumulative rain to interval rain
#'
#' Takes diff of cumulative readings. Negative diffs (counter resets) are
#' set to NA. Large jumps (> max_interval_mm) are also set to NA.
#'
#' @param rain numeric vector of PpRain values
#' @param is_cum logical vector from detect_cumulative()
#' @param max_interval_mm numeric, cap on single-interval rain (default 50 mm)
#' @return numeric vector of per-interval rain (mm), same length as rain
convert_cumulative <- function(rain, is_cum, max_interval_mm = 50) {
  result <- rain  # default: keep raw (rate data)

  if (!any(is_cum)) return(result)

  # For cumulative stretches, compute diff within each contiguous block
  cum_rle <- rle(is_cum)
  end_idx <- cumsum(cum_rle$lengths)
  start_idx <- end_idx - cum_rle$lengths + 1

  for (i in seq_along(cum_rle$values)) {
    if (cum_rle$values[i]) {
      idx <- start_idx[i]:end_idx[i]
      d <- c(0, diff(rain[idx]))
      # Counter resets (negative) -> NA
      d[d < 0] <- NA_real_
      # Unreasonable jumps -> NA
      d[!is.na(d) & d > max_interval_mm] <- NA_real_
      result[idx] <- d
    }
  }

  result
}

#' Apply QC thresholds to rain rate data
#'
#' Based on template QC: cap negatives, flag values > max_rate, flag
#' single-step changes > max_change.
#'
#' @param rain_mmhr numeric vector of rain rate (mm/hr)
#' @param max_rate numeric, maximum plausible rain rate (default 150 mm/hr)
#' @param max_change numeric, maximum plausible single-step change (default 100 mm/hr)
#' @return numeric vector with bad values set to NA
qc_rain_rate <- function(rain_mmhr, max_rate = 150, max_change = 100) {
  rain_mmhr[rain_mmhr < 0] <- NA_real_
  rain_mmhr[!is.na(rain_mmhr) & rain_mmhr > max_rate] <- NA_real_
  d <- c(0, abs(diff(rain_mmhr)))
  rain_mmhr[!is.na(d) & d > max_change] <- NA_real_
  rain_mmhr
}

#' Load rainfall data for a single lake, converting to common units (mm/hr)
#'
#' Automatically detects cumulative (Vaisala) vs rate (Lufft) periods using
#' run-length analysis of non-decreasing sequences. Converts cumulative
#' readings to interval rain via diff(), then to mm/hr.
#'
#' @param lake_name character
#' @param rds_dir character
#' @return data.table with DateTime, PpRain_mmhr, is_cumulative, lake, etc.
load_rainfall <- function(lake_name, rds_dir = "data/rds") {
  path <- file.path(rds_dir, paste0("full_", lake_name, ".rds"))
  x <- readRDS(path)
  met <- as.data.table(x[["raw"]][["met"]])
  met <- met[, .(DateTime, PpRain)]
  met <- met[!is.na(PpRain)]

  if (nrow(met) == 0) return(data.table())

  setorder(met, DateTime)

  # Compute timestep in hours
  met[, dt_hrs := c(NA_real_, as.numeric(diff(DateTime), units = "hours"))]
  med_dt <- median(met$dt_hrs, na.rm = TRUE)
  met[is.na(dt_hrs), dt_hrs := med_dt]
  met[dt_hrs <= 0 | dt_hrs > 24, dt_hrs := med_dt]

  # Detect and convert cumulative periods
  met[, is_cumulative := detect_cumulative(PpRain)]
  n_cum <- sum(met$is_cumulative)
  if (n_cum > 0) {
    cum_pct <- round(100 * n_cum / nrow(met), 1)
    cat("  Cumulative rain detected:", n_cum, "readings (", cum_pct, "%)\n")
  }

  met[, PpRain_interval := convert_cumulative(PpRain, is_cumulative)]

  # Convert interval rain (mm) to rate (mm/hr)
  # Cumulative: interval_mm / dt_hrs
  # Rate: raw value is already mm/hr
  met[, PpRain_mmhr := fifelse(is_cumulative,
                                PpRain_interval / dt_hrs,
                                PpRain)]

  # Flag remaining high-value periods as unknown/bad data.
  # If monthly median PpRain (when positive) > 50 and not already handled
  # as cumulative, these are raw sensor values in unknown units — set to NA.
  met[, ym := format(DateTime, "%Y-%m")]
  bad_months <- met[PpRain_mmhr > 0 & !is_cumulative,
                    .(med = median(PpRain_mmhr)), by = ym][med > 50, ym]
  if (length(bad_months) > 0) {
    cat("  Bad rain months (unknown units):", length(bad_months), "\n")
    met[ym %in% bad_months, PpRain_mmhr := NA_real_]
  }

  # QC thresholds
  met[, PpRain_mmhr := qc_rain_rate(PpRain_mmhr)]
  met[is.na(PpRain_mmhr), PpRain_mmhr := 0]

  # Add time fields
  met[, `:=`(
    lake       = lake_name,
    hydro_year = ifelse(as.integer(format(DateTime, "%m")) >= 7,
                        as.integer(format(DateTime, "%Y")),
                        as.integer(format(DateTime, "%Y")) - 1L),
    month      = as.integer(format(DateTime, "%m")),
    season     = fcase(
      as.integer(format(DateTime, "%m")) %in% c(12, 1, 2), "Summer",
      as.integer(format(DateTime, "%m")) %in% c(3, 4, 5),  "Autumn",
      as.integer(format(DateTime, "%m")) %in% c(6, 7, 8),  "Winter",
      as.integer(format(DateTime, "%m")) %in% c(9, 10, 11), "Spring"
    )
  )]

  met[, .(DateTime, PpRain_mmhr, dt_hrs, is_cumulative, lake,
          hydro_year, month, season)]
}

#' Compute annual rainfall with seasonal adjustment for missing data
#'
#' @param rain_dt data.table from load_rainfall()
#' @return data.table with hydro_year, raw_total_mm, adjusted_total_mm,
#'   coverage_frac, lake
calc_annual_rainfall <- function(rain_dt) {
  if (nrow(rain_dt) == 0) return(data.table())

  rain_dt[, rain_mm := PpRain_mmhr * dt_hrs]

  season_days <- c(Summer = 90, Autumn = 92, Winter = 92, Spring = 91)
  season_hrs <- season_days * 24

  seasonal <- rain_dt[, .(
    total_mm  = sum(rain_mm, na.rm = TRUE),
    obs_hours = sum(dt_hrs, na.rm = TRUE)
  ), by = .(lake, hydro_year, season)]

  seasonal[, expected_hrs := season_hrs[season]]
  seasonal[, coverage := pmin(obs_hours / expected_hrs, 1)]
  seasonal[, adjusted_mm := fifelse(coverage > 0.1,
                                     total_mm / coverage,
                                     NA_real_)]

  annual <- seasonal[, .(
    raw_total_mm      = sum(total_mm, na.rm = TRUE),
    adjusted_total_mm = fifelse(all(!is.na(adjusted_mm)),
                                 sum(adjusted_mm, na.rm = TRUE),
                                 NA_real_),
    coverage_frac     = mean(coverage)
  ), by = .(lake, hydro_year)]

  setorder(annual, lake, hydro_year)
  annual
}

#' Load and compute rainfall for all lakes
#'
#' @param rds_dir character
#' @return list with rain_annual and rain_raw data.tables
run_all_rainfall <- function(rds_dir = "data/rds") {
  files <- list.files(rds_dir, pattern = "^full_.*[.]rds$", full.names = FALSE)
  lake_names <- sub("^full_", "", sub("[.]rds$", "", files))

  all_rain <- rbindlist(lapply(lake_names, function(lk) {
    cat("Loading", lk, "...\n")
    tryCatch(load_rainfall(lk, rds_dir),
             error = function(e) {
               cat("  ERROR:", conditionMessage(e), "\n")
               data.table()
             })
  }), fill = TRUE)

  annual <- calc_annual_rainfall(all_rain)
  list(raw = all_rain, annual = annual)
}

#' Bar plot of annual rainfall by hydrological year, per lake
#'
#' @param annual_dt data.table from calc_annual_rainfall()
#' @return ggplot object
plot_annual_rainfall <- function(annual_dt) {
  dt <- melt(annual_dt, id.vars = c("lake", "hydro_year", "coverage_frac"),
             measure.vars = c("raw_total_mm", "adjusted_total_mm"),
             variable.name = "type", value.name = "total_mm")
  dt[, type := fifelse(type == "raw_total_mm", "Raw", "Adjusted")]
  dt[, hy_label := factor(hydro_year)]

  ggplot(dt[!is.na(total_mm)], aes(x = hy_label, y = total_mm, fill = type)) +
    geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
    scale_fill_manual(values = c("Raw" = "steelblue", "Adjusted" = "coral")) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    facet_wrap(~lake, scales = "free_y", ncol = 4) +
    labs(
      x = "Hydrological year (Jul-Jun)",
      y = "Annual rainfall (mm)",
      title = "Annual rainfall by hydrological year",
      subtitle = "Raw totals vs seasonally adjusted for missing data",
      fill = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 8)
    )
}
