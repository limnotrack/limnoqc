#' Detect DO sensor gain drift
#'
#' Analyses DO data for gain drift by comparing surface readings during
#' well-mixed periods (where DO should be near 100% saturation) against
#' expected values. Returns a time series of estimated gain factors.
#'
#' @param datetime POSIXct vector
#' @param do_psat DO in percent saturation (or mg/L if temp provided)
#' @param temp Water temperature in °C (required if do_psat is in mg/L)
#' @param depth Depth in metres (for profiling data — uses shallow readings only)
#' @param input_unit "psat" or "mgl" (default "psat")
#' @param pressure Barometric pressure in hPa (for mg/L → %sat conversion)
#' @param max_depth Maximum depth to use for surface readings (default 5m)
#' @param window_days Rolling window for gain estimation (default 30)
#' @param expected_psat Expected surface %sat during well-mixed conditions
#'   (default 98). Can be a single value or a named list with monthly values.
#' @return Data frame with columns:
#'   - `date`: date
#'   - `gain`: estimated sensor gain (1.0 = correct, <1 = reading low)
#'   - `n`: number of surface readings used
#'   - `median_psat`: median observed surface %sat
#' @export
do_detect_drift <- function(datetime, do_psat, temp = NULL, depth = NULL,
                             input_unit = "psat",
                             pressure = 1013.25,
                             max_depth = 5,
                             window_days = 30,
                             expected_psat = 98) {

  # Convert to %sat if needed
  if (input_unit == "mgl") {
    if (is.null(temp)) stop("temp required when input_unit = 'mgl'")
    do_psat <- do_mgl_to_psat(do_psat, temp, pressure)
  }

  # Filter to surface readings
  if (!is.null(depth)) {
    mask <- !is.na(depth) & depth <= max_depth & !is.na(do_psat)
  } else {
    mask <- !is.na(do_psat)
  }

  dt <- datetime[mask]
  psat <- do_psat[mask]

  if (length(dt) < 10) {
    warning("Too few surface DO readings for drift detection")
    return(NULL)
  }

  # Rolling window analysis
  dates <- seq(as.Date(min(dt)), as.Date(max(dt)), by = 1)
  dt_num <- as.numeric(dt)
  window_secs <- window_days * 86400

  results <- lapply(dates, function(d) {
    d_num <- as.numeric(as.POSIXct(d))
    in_window <- abs(dt_num - d_num) <= window_secs / 2
    if (sum(in_window) < 5) return(NULL)

    med_psat <- stats::median(psat[in_window], na.rm = TRUE)

    # Expected %sat (can vary by month)
    if (is.list(expected_psat)) {
      m <- as.numeric(format(d, "%m"))
      exp_val <- expected_psat[[as.character(m)]]
      if (is.null(exp_val)) exp_val <- 98
    } else {
      exp_val <- expected_psat
    }

    data.frame(
      date = d,
      gain = med_psat / exp_val,
      n = sum(in_window),
      median_psat = med_psat
    )
  })

  do.call(rbind, results)
}


#' Apply gain correction to DO data
#'
#' Divides DO values by a gain factor to correct for sensor drift.
#' Supports constant gain, step-change gain, or time-varying gain
#' from [do_detect_drift()].
#'
#' @param do_value DO values (mg/L or %sat)
#' @param datetime POSIXct vector (required for time-varying correction)
#' @param gain Either:
#'   - A single numeric value (constant gain applied to all data)
#'   - A data frame with `date` and `gain` columns (from [do_detect_drift()])
#'   - A list with `date` (POSIXct) and `gain` (numeric) for step changes
#' @param method How to apply time-varying gain:
#'   - "step": apply the nearest gain value (step function)
#'   - "interpolate": linearly interpolate between gain values (default)
#' @return Corrected DO values (same units as input)
#' @export
#'
#' @examples
#' # Constant gain correction
#' do_corrected <- do_correct_gain(do_raw, gain = 0.75)
#'
#' # Step change at a known date
#' do_corrected <- do_correct_gain(do_raw, datetime = dt,
#'   gain = list(date = as.POSIXct("2024-10-01"), gain = 0.75))
do_correct_gain <- function(do_value, datetime = NULL, gain,
                             method = "interpolate") {

  if (is.numeric(gain) && length(gain) == 1) {
    # Constant gain
    return(do_value / gain)
  }

  if (is.null(datetime)) stop("datetime required for time-varying gain correction")

  if (is.data.frame(gain)) {
    gain_dates <- as.numeric(as.POSIXct(gain$date))
    gain_values <- gain$gain
  } else if (is.list(gain)) {
    gain_dates <- as.numeric(gain$date)
    gain_values <- gain$gain
  } else {
    stop("gain must be numeric, data.frame, or list")
  }

  dt_num <- as.numeric(datetime)

  if (method == "interpolate") {
    gain_at_t <- stats::approx(gain_dates, gain_values, xout = dt_num,
                                rule = 2)$y
  } else {
    # Step function: use nearest preceding gain value
    idx <- findInterval(dt_num, gain_dates)
    idx[idx == 0] <- 1
    gain_at_t <- gain_values[idx]
  }

  do_value / gain_at_t
}
