#' Detect step changes (abrupt persistent level shifts)
#'
#' Compares the mean of a trailing window against a leading window at
#' each point. A step change is where the difference exceeds a threshold
#' and persists (unlike a spike, which returns immediately).
#'
#' Useful for detecting sensor swaps, calibration drift onset, or
#' biofouling transitions.
#'
#' @param x Numeric vector.
#' @param window Integer. Half-window size (number of points on each side).
#'   Default 50.
#' @param threshold Numeric. Minimum absolute difference between trailing
#'   and leading window means to flag as a step. If \code{NULL}, auto-set
#'   to 3x the rolling MAD of \code{x}.
#' @return Data frame with columns:
#'   \describe{
#'     \item{index}{Row index of the step}
#'     \item{magnitude}{Signed difference (leading mean - trailing mean)}
#'   }
#'   Empty data frame if no steps detected.
#'
#' @export
#'
#' @examples
#' x <- c(rnorm(200, 10), rnorm(200, 15), rnorm(200, 10))
#' detect_step(x, window = 30, threshold = 2)
detect_step <- function(x, window = 50L, threshold = NULL) {
  n <- length(x)
  w <- as.integer(window)
  if (n < 2L * w + 1L)
    return(data.frame(index = integer(), magnitude = numeric()))

  # Compute trailing and leading means using cumsum (O(n), no loops)
  cs <- cumsum(ifelse(is.na(x), 0, x))
  cn <- cumsum(!is.na(x))

  trail_sum <- cs[(w + 1L):(n - w)] - c(0, cs[1:(n - 2L * w)])
  trail_n   <- cn[(w + 1L):(n - w)] - c(0L, cn[1:(n - 2L * w)])
  lead_sum  <- cs[(2L * w + 1L):n] - cs[(w + 1L):(n - w)]
  lead_n    <- cn[(2L * w + 1L):n] - cn[(w + 1L):(n - w)]

  trail_mean <- ifelse(trail_n > 0, trail_sum / trail_n, NA_real_)
  lead_mean  <- ifelse(lead_n > 0, lead_sum / lead_n, NA_real_)

  diff_val <- lead_mean - trail_mean

  # Auto-threshold from MAD if not supplied
  if (is.null(threshold)) {
    overall_mad <- stats::mad(x, na.rm = TRUE)
    threshold <- 3 * overall_mad
    if (threshold == 0) threshold <- stats::sd(x, na.rm = TRUE) * 3
  }

  # Find where absolute difference exceeds threshold
  idx_offset <- w  # results are offset by w from the start
  candidates <- which(!is.na(diff_val) & abs(diff_val) > threshold)

  if (length(candidates) == 0L)
    return(data.frame(index = integer(), magnitude = numeric()))

  # Keep only the peak of each cluster (local maxima of |diff|)
  abs_d <- abs(diff_val)
  peaks <- candidates[c(TRUE, diff(candidates) > w)]  # at least w apart

  data.frame(
    index     = peaks + idx_offset,
    magnitude = diff_val[peaks]
  )
}


#' Detect loss of diurnal cycle
#'
#' For variables with an expected diurnal signal (solar radiation, air
#' temperature, surface water temperature), flags periods where the daily
#' range drops below a threshold. A sensor that stops showing day/night
#' variation has likely failed.
#'
#' @param x Numeric vector.
#' @param datetime POSIXct vector.
#' @param min_daily_range Numeric. Minimum expected daily range. If the
#'   daily max - min falls below this, the day is flagged.
#' @param n_days Integer. Flag only after \code{n_days} consecutive days
#'   below the threshold. Default 3 (avoids false positives from single
#'   overcast days).
#' @return Logical vector, same length as \code{x}. \code{TRUE} = flagged
#'   (within a period of lost diurnal cycle).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Solar radiation should have > 50 W/m² daily range
#' lost <- detect_diurnal_loss(met$r_swd_wm2, met$datetime,
#'                              min_daily_range = 50, n_days = 3)
#' }
detect_diurnal_loss <- function(x, datetime, min_daily_range,
                                n_days = 3L) {
  n <- length(x)
  flagged <- rep(FALSE, n)
  if (n < 10L) return(flagged)

  # Compute daily range
  day <- as.Date(datetime)
  udays <- sort(unique(day))

  daily_range <- vapply(udays, function(d) {
    vals <- x[day == d]
    vals <- vals[!is.na(vals)]
    if (length(vals) < 4L) return(NA_real_)
    max(vals) - min(vals)
  }, numeric(1))

  # Find runs of low-range days
  low <- !is.na(daily_range) & daily_range < min_daily_range
  r <- rle(low)
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L

  bad_days <- logical(length(udays))
  for (i in seq_along(r$lengths)) {
    if (isTRUE(r$values[i]) && r$lengths[i] >= n_days)
      bad_days[starts[i]:ends[i]] <- TRUE
  }

  # Map back to rows
  bad_day_set <- udays[bad_days]
  flagged <- day %in% bad_day_set
  flagged[is.na(x)] <- FALSE
  flagged
}


#' Detect density inversions in profiles
#'
#' Flags scans where temperature implies an unstable density profile
#' (lighter water beneath heavier water). In freshwater, density is
#' almost entirely determined by temperature, with maximum density at
#' ~3.98 deg C.
#'
#' The test computes freshwater density from temperature at each scan
#' and flags where density decreases with depth within a profile.
#'
#' @param df Data frame of profile data.
#' @param col_profile Profile ID column (default \code{"profile_id"}).
#' @param col_depth Depth column (default \code{"z_meas_m"}).
#' @param col_temp Temperature column (default \code{"t_wtr_degc"}).
#' @param threshold Minimum density decrease (kg/m³) to flag. Default
#'   0.001 (filters out noise-level inversions).
#' @return Logical vector, same length as \code{nrow(df)}.
#'   \code{TRUE} = density inversion at this scan.
#'
#' @details
#' Freshwater density from temperature uses the Kell (1975) equation:
#' \deqn{\rho = 999.842594 + 6.793952e{-2}T - 9.095290e{-3}T^2 +
#'   1.001685e{-4}T^3 - 1.120083e{-6}T^4 + 6.536336e{-9}T^5}
#'
#' Scans are compared to the next-deeper scan within each profile.
#' Both temperature and depth must be non-NA.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' inv <- detect_density_inversion(pro)
#' cat("Inversions:", sum(inv), "of", nrow(pro), "scans\n")
#' }
detect_density_inversion <- function(df,
                                     col_profile = "profile_id",
                                     col_depth   = "z_meas_m",
                                     col_temp    = "t_wtr_degc",
                                     threshold   = 0.001) {
  n <- nrow(df)
  inversion <- rep(FALSE, n)

  temp  <- df[[col_temp]]
  depth <- df[[col_depth]]
  pid   <- df[[col_profile]]

  # Freshwater density (Kell 1975)
  rho <- .fw_density(temp)

  pids <- unique(pid[!is.na(pid)])

  for (p in pids) {
    idx <- which(pid == p)
    if (length(idx) < 2L) next

    # Sort by depth within profile
    ord <- idx[order(depth[idx])]
    d <- rho[ord]

    # Flag where density decreases with depth
    dd <- diff(d)
    bad <- which(!is.na(dd) & dd < -threshold)
    if (length(bad) > 0L) {
      # Flag the deeper scan of each inversion pair
      inversion[ord[bad + 1L]] <- TRUE
    }
  }
  inversion
}


#' Check depth monotonicity within profiles
#'
#' Flags scans where measured depth decreases from the previous scan
#' within the same profile. A profiler cast is always downward —
#' decreasing depth indicates a winch fault, depth sensor problem,
#' or slipping drum.
#'
#' @param df Data frame of profile data.
#' @param col_profile Profile ID column (default \code{"profile_id"}).
#' @param col_depth Depth column (default \code{"z_meas_m"}).
#' @param tolerance Numeric. Allow depth to decrease by up to this many
#'   metres without flagging (accounts for sensor noise). Default 0.1.
#' @return Logical vector, same length as \code{nrow(df)}.
#'   \code{TRUE} = depth decreased from previous scan in same profile.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' nonmono <- detect_depth_nonmonotonic(pro)
#' cat("Non-monotonic:", sum(nonmono), "scans\n")
#' }
detect_depth_nonmonotonic <- function(df,
                                      col_profile = "profile_id",
                                      col_depth   = "z_meas_m",
                                      tolerance   = 0.1) {
  n <- nrow(df)
  flagged <- rep(FALSE, n)

  depth <- df[[col_depth]]
  pid   <- df[[col_profile]]

  pids <- unique(pid[!is.na(pid)])

  for (p in pids) {
    idx <- which(pid == p)
    if (length(idx) < 2L) next

    d <- depth[idx]
    dd <- diff(d)
    bad <- which(!is.na(dd) & dd < -tolerance)
    if (length(bad) > 0L)
      flagged[idx[bad + 1L]] <- TRUE
  }
  flagged
}


#' Flag non-photochemical quenching (NPQ) in chlorophyll fluorescence
#'
#' Surface chlorophyll fluorescence is suppressed at midday under high
#' light due to non-photochemical quenching. This is real biology, not
#' a sensor fault — but the depressed readings don't represent true
#' chlorophyll concentration and should be flagged so users don't
#' interpret the midday dip as a real decrease.
#'
#' Flags surface scans during high-radiation periods. If radiation data
#' is available, uses a radiation threshold. Otherwise, uses solar noon
#' ± a time window.
#'
#' @param df Data frame with datetime, depth, and fluorescence columns.
#' @param col_datetime Datetime column (default \code{"datetime"}).
#' @param col_depth Depth column (default \code{"z_meas_m"}).
#' @param col_fluor Fluorescence column (default \code{"f_chl_rfu"}).
#' @param col_rad Radiation column (default \code{"r_swd_wm2"}).
#'   If not present in \code{df}, time-of-day method is used instead.
#' @param max_depth Numeric. Maximum depth affected by NPQ (metres).
#'   Default 5.
#' @param rad_threshold Numeric. Minimum radiation (W/m²) for quenching
#'   to occur. Default 200.
#' @param noon_hours Numeric. Half-width of the midday window (hours
#'   from solar noon) used when radiation data is unavailable. Default 3.
#' @param solar_noon_hour Numeric. Approximate hour of solar noon in
#'   local time. Default 12.5 (NZST, ~solar noon for NZ).
#' @return Logical vector, same length as \code{nrow(df)}.
#'   \code{TRUE} = likely quenched.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' npq <- detect_quenching(pro)
#' cat("Quenched scans:", sum(npq), "\n")
#' }
detect_quenching <- function(df,
                             col_datetime = "datetime",
                             col_depth    = "z_meas_m",
                             col_fluor    = "f_chl_rfu",
                             col_rad      = "r_swd_wm2",
                             max_depth    = 5,
                             rad_threshold = 200,
                             noon_hours   = 3,
                             solar_noon_hour = 12.5) {
  n <- nrow(df)
  quenched <- rep(FALSE, n)

  if (!(col_fluor %in% names(df))) return(quenched)

  depth <- df[[col_depth]]
  fluor <- df[[col_fluor]]
  dt    <- df[[col_datetime]]

  # Only surface
  shallow <- !is.na(depth) & depth <= max_depth & !is.na(fluor)

  if (col_rad %in% names(df)) {
    # Radiation-based: flag surface scans during high radiation
    rad <- df[[col_rad]]
    # Radiation might be on a different time grid — match nearest
    # For profile data, radiation is typically from the met table.
    # If it's in df (merged), use it directly.
    high_rad <- !is.na(rad) & rad > rad_threshold
    quenched <- shallow & high_rad
  } else {
    # Time-of-day fallback: flag around solar noon
    hour_of_day <- as.numeric(format(dt, "%H")) +
                   as.numeric(format(dt, "%M")) / 60
    midday <- abs(hour_of_day - solar_noon_hour) <= noon_hours
    quenched <- shallow & midday
  }

  quenched
}


# -- internal helpers --------------------------------------------------------

#' Freshwater density from temperature (Kell 1975)
#'
#' @param temp Temperature in degrees Celsius
#' @return Density in kg/m³
#' @noRd
.fw_density <- function(temp) {
  999.842594 +
    6.793952e-2  * temp -
    9.095290e-3  * temp^2 +
    1.001685e-4  * temp^3 -
    1.120083e-6  * temp^4 +
    6.536336e-9  * temp^5
}
