#' Detect periods where the profiler is stuck (not profiling)
#'
#' Scans a profile dataset for periods where the sensor depth stays
#' within a narrow band for an extended duration — indicating the
#' profiler winch has stalled, the probe is parked, or the system
#' is in a fault state. Returns date ranges of stuck periods.
#'
#' The detection works by computing the rolling depth range over a
#' time window. Normal profiling produces large depth variation
#' (surface to bottom every few hours). When the depth range drops
#' below a threshold for longer than \code{min_duration_hours}, the
#' period is flagged as stuck.
#'
#' @param df Data frame with datetime and depth columns.
#' @param col_datetime Datetime column (default \code{"datetime"}).
#' @param col_depth Depth column (default \code{"z_meas_m"}).
#' @param depth_range_threshold Numeric. If the depth range within a
#'   window stays below this (metres), the profiler is considered stuck.
#'   Default 2 — normal casts span 10-30m, stuck probe varies < 2m.
#' @param window_hours Numeric. Rolling window width for computing
#'   depth range (hours). Default 6 (two typical cast cycles).
#' @param min_duration_hours Numeric. Minimum duration (hours) of low
#'   depth range to flag as stuck. Short pauses (maintenance, reboot)
#'   are normal. Default 12.
#' @param gap_hours Numeric. Merge stuck periods separated by less than
#'   this many hours. Default 6.
#' @return Data frame with columns:
#'   \describe{
#'     \item{date_from}{POSIXct start of stuck period}
#'     \item{date_to}{POSIXct end of stuck period}
#'     \item{duration_hours}{Duration in hours}
#'     \item{mean_depth}{Mean depth during stuck period (where it's parked)}
#'     \item{depth_range}{Depth range during stuck period}
#'     \item{n_readings}{Number of readings in the period}
#'   }
#'   Empty data frame if no stuck periods detected.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' stuck <- detect_stuck_profiler(pro)
#' cat(nrow(stuck), "stuck periods found\n")
#' print(stuck)
#' }
detect_stuck_profiler <- function(df,
                                  col_datetime = "datetime",
                                  col_depth    = "z_meas_m",
                                  depth_range_threshold = 2,
                                  window_hours = 6,
                                  min_duration_hours = 12,
                                  gap_hours = 6) {

  empty <- data.frame(
    date_from      = as.POSIXct(character()),
    date_to        = as.POSIXct(character()),
    duration_hours = numeric(),
    mean_depth     = numeric(),
    depth_range    = numeric(),
    n_readings     = integer(),
    stringsAsFactors = FALSE
  )

  dt    <- df[[col_datetime]]
  depth <- df[[col_depth]]

  valid <- !is.na(dt) & !is.na(depth)
  if (sum(valid) < 10L) return(empty)

  dt    <- dt[valid]
  depth <- depth[valid]

  # Sort by time
  ord   <- order(dt)
  dt    <- dt[ord]
  depth <- depth[ord]

  n <- length(dt)
  dt_num <- as.numeric(dt)
  window_secs <- window_hours * 3600

  # Rolling depth range using two-pointer window
  is_flat <- rep(FALSE, n)
  left <- 1L

  for (right in seq_len(n)) {
    while (left < right && (dt_num[right] - dt_num[left]) > window_secs)
      left <- left + 1L

    if ((right - left) < 2L) next

    d_chunk <- depth[left:right]
    d_range <- max(d_chunk, na.rm = TRUE) - min(d_chunk, na.rm = TRUE)

    if (d_range < depth_range_threshold)
      is_flat[right] <- TRUE
  }

  if (!any(is_flat)) return(empty)

  # Convert flat flags to contiguous runs
  r <- rle(is_flat)
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L

  periods <- list()
  for (i in seq_along(r$lengths)) {
    if (!r$values[i]) next

    s <- starts[i]
    e <- ends[i]
    dur <- as.numeric(difftime(dt[e], dt[s], units = "hours"))

    if (dur < min_duration_hours) next

    d_slice <- depth[s:e]
    periods[[length(periods) + 1L]] <- data.frame(
      date_from      = dt[s],
      date_to        = dt[e],
      duration_hours = round(dur, 1),
      mean_depth     = round(mean(d_slice, na.rm = TRUE), 2),
      depth_range    = round(max(d_slice, na.rm = TRUE) - min(d_slice, na.rm = TRUE), 2),
      n_readings     = e - s + 1L,
      stringsAsFactors = FALSE
    )
  }

  if (length(periods) == 0L) return(empty)

  result <- do.call(rbind, periods)

  # Merge nearby periods (separated by less than gap_hours)
  if (nrow(result) > 1L) {
    merged <- list(result[1L, ])
    for (i in 2:nrow(result)) {
      prev <- merged[[length(merged)]]
      gap <- as.numeric(difftime(result$date_from[i], prev$date_to, units = "hours"))

      if (gap <= gap_hours) {
        # Merge: extend the previous period
        s_idx <- which(dt >= prev$date_from & dt <= result$date_to[i])
        d_slice <- depth[s_idx]
        prev$date_to        <- result$date_to[i]
        prev$duration_hours <- round(as.numeric(difftime(prev$date_to, prev$date_from, units = "hours")), 1)
        prev$mean_depth     <- round(mean(d_slice, na.rm = TRUE), 2)
        prev$depth_range    <- round(max(d_slice, na.rm = TRUE) - min(d_slice, na.rm = TRUE), 2)
        prev$n_readings     <- length(s_idx)
        merged[[length(merged)]] <- prev
      } else {
        merged[[length(merged) + 1L]] <- result[i, ]
      }
    }
    result <- do.call(rbind, merged)
  }

  rownames(result) <- NULL
  result
}
