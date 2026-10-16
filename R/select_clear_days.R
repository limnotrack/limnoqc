library(data.table)

#' Select the clearest days from met data based on radiation variability
#'
#' Ranks days by a "clarity score" — the sum of absolute rate-of-change of
#' observed RadSWD, normalised by the daily clear-sky total in the solar noon
#' window. Lower score = smoother = clearer sky.
#'
#' @param met_dt data.table with DateTime and RadSWD columns
#' @param clearsky_dt data.table with DateTime and clearsky_wm2 columns
#'   (already filtered to the solar noon window by calc_clearsky_window)
#' @param solar_noon_dt data.table with date and solar_noon columns
#' @param hours_window numeric, hours either side of solar noon (default 3)
#' @param n_per_month integer, number of clear days to select per month (default 1)
#' @param clarity_threshold numeric or NULL, optional upper bound on clarity
#'   score. Days with score > threshold are excluded. NULL = no threshold.
#' @param min_readings integer, minimum number of readings in the window for a
#'   day to be considered (default 12, i.e. ~1 hour of 5-min data)
#' @param min_obs_fraction numeric, minimum ratio of mean_obs / mean_clearsky
#'   for a day to qualify (default 0.5). Filters out overcast days that happen
#'   to have smooth radiation curves.
#' @return data.table with columns: date, year_month, clarity_score, n_readings,
#'   mean_obs, mean_clearsky
select_clear_days <- function(met_dt, clearsky_dt, solar_noon_dt,
                               hours_window = 3, n_per_month = 1,
                               clarity_threshold = NULL, min_readings = 12,
                               min_obs_fraction = 0.5) {

  # Merge observed with clearsky
  merged <- merge(
    met_dt[, .(DateTime, RadSWD)],
    clearsky_dt[, .(DateTime, clearsky_wm2)],
    by = "DateTime"
  )
  merged[, date := as.Date(DateTime)]

  # Add solar noon and filter to window
  merged <- merge(merged, solar_noon_dt, by = "date")
  merged[, hrs_from_noon := as.numeric(difftime(DateTime, solar_noon, units = "hours"))]
  merged <- merged[abs(hrs_from_noon) <= hours_window]

  if (nrow(merged) == 0) return(data.table())

  # Remove NaN/NA and values < 5 W/m2 (noise floor)
  merged <- merged[!is.na(RadSWD) & RadSWD >= 5]

  if (nrow(merged) == 0) return(data.table())

  # Order by time within each day

  setorder(merged, date, DateTime)

  # Per-day clarity score
  daily <- merged[, {
    n <- .N
    if (n < 2) {
      list(clarity_score = NA_real_, n_readings = n,
           mean_obs = mean(RadSWD),
           mean_clearsky = mean(clearsky_wm2))
    } else {
      d_rad <- abs(diff(RadSWD))
      sum_drad <- sum(d_rad)
      obs_total <- sum(RadSWD)
      score <- if (obs_total > 0) sum_drad / obs_total else NA_real_
      list(clarity_score = score, n_readings = n,
           mean_obs = mean(RadSWD),
           mean_clearsky = mean(clearsky_wm2))
    }
  }, by = date]

  # Filter
  daily <- daily[n_readings >= min_readings & !is.na(clarity_score)]
  daily <- daily[mean_obs / mean_clearsky >= min_obs_fraction]
  if (!is.null(clarity_threshold)) {
    daily <- daily[clarity_score <= clarity_threshold]
  }

  if (nrow(daily) == 0) return(data.table())

  # Select best n_per_month per calendar month
  daily[, year_month := format(date, "%Y-%m")]
  setorder(daily, year_month, clarity_score)
  selected <- daily[, head(.SD, n_per_month), by = year_month]

  # Remove clarity score outliers (IQR upper fence)
  q1 <- quantile(selected$clarity_score, 0.25)
  q3 <- quantile(selected$clarity_score, 0.75)
  upper_fence <- q3 + 1.5 * (q3 - q1)
  n_before <- nrow(selected)
  selected <- selected[clarity_score <= upper_fence]
  if (nrow(selected) < n_before) {
    cat("  Removed", n_before - nrow(selected),
        "clarity score outliers (fence:", round(upper_fence, 4), ")\n")
  }

  setorder(selected, date)
  selected
}
