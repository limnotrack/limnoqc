#' Detect spikes using rate-of-change test
#'
#' A spike is a value that departs sharply from its neighbours and
#' immediately returns. The test checks whether the absolute difference
#' between a value and \emph{both} its neighbours exceeds a threshold.
#' This distinguishes genuine spikes (which return) from step changes
#' (which persist).
#'
#' @param x Numeric vector.
#' @param threshold Numeric scalar. Maximum allowable single-step change.
#'   Values where \code{abs(x[i] - x[i-1]) > threshold AND
#'   abs(x[i] - x[i+1]) > threshold} are flagged as spikes.
#' @return Logical vector, same length as \code{x}. \code{TRUE} = spike.
#'   First and last elements are always \code{FALSE}.
#'
#' @details
#' This is the QARTOD "spike test": a value is suspect if it deviates
#' from both neighbours by more than the threshold. A step change only
#' deviates from one neighbour, so it passes.
#'
#' For a stricter test, compare against the mean of neighbours:
#' \code{abs(x[i] - (x[i-1] + x[i+1])/2) > threshold}. This catches
#' smaller spikes but is more sensitive to noise.
#'
#' @export
#'
#' @examples
#' x <- c(10, 10.1, 10, 50, 10.2, 10, 10.1)
#' detect_spikes(x, threshold = 5)
#' # FALSE FALSE FALSE  TRUE FALSE FALSE FALSE
detect_spikes <- function(x, threshold) {
  n <- length(x)
  if (n < 3L) return(rep(FALSE, n))

  spike <- rep(FALSE, n)
  d_prev <- c(NA, abs(diff(x)))
  d_next <- c(abs(diff(x)), NA)

  spike <- !is.na(d_prev) & !is.na(d_next) &
           d_prev > threshold & d_next > threshold

  # NA values are not spikes
  spike[is.na(x)] <- FALSE
  spike
}


#' Detect stuck sensor (flat-line) periods
#'
#' Finds runs of consecutive identical values exceeding a minimum length.
#' A sensor reporting the same value repeatedly has likely failed or
#' frozen. The test uses exact equality by default, or a tolerance band
#' for noisy sensors.
#'
#' @param x Numeric vector.
#' @param n_consec Integer. Minimum number of consecutive identical values
#'   to flag as stuck. Default 10.
#' @param tolerance Numeric. Values within \code{tolerance} of each other
#'   are considered identical. Default 0 (exact match). Set to the sensor's
#'   resolution for noisy channels.
#' @param keep_first Logical. If \code{TRUE} (default), the first value
#'   in each stuck run is not flagged — it's the last real reading before
#'   the sensor froze. If \code{FALSE}, the entire run is flagged.
#' @return Logical vector, same length as \code{x}. \code{TRUE} = stuck.
#'
#' @details
#' Uses run-length encoding for speed — no loops. For a 2M-row vector
#' this completes in milliseconds.
#'
#' NA values break runs (a sensor reporting NA is not "stuck", it's
#' offline).
#'
#' @export
#'
#' @examples
#' x <- c(1:5, rep(7.0, 20), 8:10)
#' detect_stuck(x, n_consec = 10)
detect_stuck <- function(x, n_consec = 10L, tolerance = 0, keep_first = TRUE) {
  n <- length(x)
  stuck <- rep(FALSE, n)
  if (n < n_consec) return(stuck)

  if (tolerance == 0) {
    # Exact match — use rle directly
    r <- rle(x)
  } else {
    # Tolerance band — quantise then rle
    # Round to tolerance to group near-identical values
    xq <- round(x / tolerance) * tolerance
    r <- rle(xq)
  }

  # Flag runs >= n_consec (skip NA runs)
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L

  for (i in seq_along(r$lengths)) {
    if (is.na(r$values[i])) next
    if (r$lengths[i] >= n_consec) {
      flag_start <- if (keep_first) starts[i] + 1L else starts[i]
      stuck[flag_start:ends[i]] <- TRUE
    }
  }
  stuck
}


#' Detect attenuated signal (biofouling / sensor degradation)
#'
#' Computes rolling variance (or range) over a time window and flags
#' periods where the signal variability drops below a fraction of the
#' baseline. Declining variance over weeks is a hallmark of biofouling
#' on optical sensors (chlorophyll, phycocyanin, turbidity).
#'
#' @param x Numeric vector. The variable to check.
#' @param datetime POSIXct vector. Timestamps.
#' @param window_days Numeric. Rolling window width in days. Default 14.
#' @param min_var_frac Numeric. Flag when rolling variance drops below
#'   this fraction of the overall variance. Default 0.1 (10%).
#' @param min_points Integer. Minimum observations in a window to compute
#'   variance. Default 50.
#' @return Logical vector, same length as \code{x}. \code{TRUE} = attenuated.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' atten <- detect_attenuation(pro$f_chl_rfu, pro$datetime,
#'                              window_days = 14, min_var_frac = 0.1)
#' }
detect_attenuation <- function(x, datetime, window_days = 14,
                               min_var_frac = 0.1, min_points = 50L) {
  n <- length(x)
  attenuated <- rep(FALSE, n)
  if (n < min_points) return(attenuated)

  # Overall variance (baseline)
  var_all <- stats::var(x, na.rm = TRUE)
  if (is.na(var_all) || var_all == 0) return(attenuated)

  threshold <- var_all * min_var_frac
  window_secs <- window_days * 86400
  dt_num <- as.numeric(datetime)

  # Sliding window variance using a two-pointer approach
  left <- 1L
  for (right in seq_len(n)) {
    # Advance left pointer to maintain window
    while (left < right && (dt_num[right] - dt_num[left]) > window_secs)
      left <- left + 1L

    if ((right - left + 1L) < min_points) next

    chunk <- x[left:right]
    chunk <- chunk[!is.na(chunk)]
    if (length(chunk) < min_points) next

    v <- stats::var(chunk)
    if (!is.na(v) && v < threshold)
      attenuated[right] <- TRUE
  }
  attenuated
}


#' Check profile completeness
#'
#' For each profile, checks whether it reached all expected target depths
#' and had sufficient scans per hold. Returns a summary table rather than
#' a flag vector — profile-level QC, not row-level.
#'
#' @param df Data frame of profile data.
#' @param col_profile Profile ID column (default \code{"profile_id"}).
#' @param col_depth Measured depth column (default \code{"z_meas_m"}).
#' @param col_target Target depth column (default \code{"z_target_m"}).
#'   If absent, only scan count and depth range are checked.
#' @param expected_targets Numeric vector of expected target depths.
#'   If \code{NULL}, inferred from the most common target set.
#' @param min_scans_per_hold Minimum scans at each target depth.
#'   Default 5.
#' @param depth_tolerance How close \code{z_meas_m} must be to
#'   \code{z_target_m} to count as "reached" (metres). Default 2.
#' @return Data frame with one row per profile:
#'   \describe{
#'     \item{profile_id}{Profile identifier}
#'     \item{n_scans}{Total scans in profile}
#'     \item{n_targets}{Number of unique targets reached}
#'     \item{n_expected}{Number of expected targets}
#'     \item{max_depth}{Deepest measurement}
#'     \item{complete}{Logical: all targets reached with enough scans}
#'     \item{missing_targets}{Comma-separated list of missed targets}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' completeness <- check_profile(pro, expected_targets = seq(1, 31, by = 2))
#' incomplete <- completeness[!completeness$complete, ]
#' }
check_profile <- function(df,
                          col_profile = "profile_id",
                          col_depth   = "z_meas_m",
                          col_target  = "z_target_m",
                          expected_targets = NULL,
                          min_scans_per_hold = 5L,
                          depth_tolerance = 2) {

  pids <- unique(df[[col_profile]])
  pids <- pids[!is.na(pids)]

  has_target <- col_target %in% names(df) && any(!is.na(df[[col_target]]))

  # Infer expected targets from the most common set
  if (is.null(expected_targets) && has_target) {
    targets <- df[[col_target]][!is.na(df[[col_target]])]
    expected_targets <- sort(unique(round(targets, 1)))
  }

  results <- vector("list", length(pids))

  for (i in seq_along(pids)) {
    pid <- pids[i]
    rows <- df[[col_profile]] == pid & !is.na(df[[col_profile]])
    sub <- df[rows, , drop = FALSE]

    n_scans <- nrow(sub)
    max_depth <- max(sub[[col_depth]], na.rm = TRUE)

    if (has_target && !is.null(expected_targets)) {
      # Check which targets were reached with enough scans
      targets_hit <- character()
      for (tgt in expected_targets) {
        near_target <- abs(sub[[col_depth]] - tgt) <= depth_tolerance
        if (sum(near_target, na.rm = TRUE) >= min_scans_per_hold)
          targets_hit <- c(targets_hit, as.character(tgt))
      }
      n_targets <- length(targets_hit)
      n_expected <- length(expected_targets)
      missed <- setdiff(as.character(expected_targets), targets_hit)
      complete <- length(missed) == 0
      missing_str <- if (length(missed)) paste(missed, collapse = ",") else ""
    } else {
      n_targets <- length(unique(round(sub[[col_depth]], 1)))
      n_expected <- NA_integer_
      complete <- NA
      missing_str <- ""
    }

    results[[i]] <- data.frame(
      profile_id      = pid,
      n_scans         = n_scans,
      n_targets       = n_targets,
      n_expected      = n_expected,
      max_depth       = max_depth,
      complete        = complete,
      missing_targets = missing_str,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}
