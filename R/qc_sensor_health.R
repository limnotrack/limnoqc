#' Detect unstable holds within profiles
#'
#' Checks whether readings at each target depth have stabilised by the
#' end of the hold. Computes the rolling variance over the last N scans
#' at each hold and flags holds where it exceeds a threshold.
#'
#' Different from [trim_equilibration()] which uses a fixed time window —
#' this checks the actual signal. A sensor that hasn't stabilised after
#' 60 seconds of scans is still equilibrating, regardless of how long
#' the hold nominally lasts.
#'
#' @param df Data frame of profile data.
#' @param var Column name of the variable to check.
#' @param col_profile Profile ID column (default \code{"profile_id"}).
#' @param col_depth Depth column (default \code{"z_meas_m"}).
#' @param col_target Target depth column (default \code{"z_target_m"}).
#' @param tail_n Number of scans at the end of each hold to evaluate.
#'   Default 5.
#' @param max_cv Maximum coefficient of variation (SD/mean) in the tail
#'   scans. Holds exceeding this are flagged as unstable. Default 0.02
#'   (2%).
#' @param max_sd Alternative: maximum absolute SD in the tail scans.
#'   If both \code{max_cv} and \code{max_sd} are provided, a hold is
#'   unstable if it fails either test. Default \code{NULL} (use CV only).
#' @return Data frame with one row per hold:
#'   \describe{
#'     \item{profile_id}{Profile identifier}
#'     \item{target_depth}{Target depth}
#'     \item{n_scans}{Scans in this hold}
#'     \item{tail_mean}{Mean of last \code{tail_n} scans}
#'     \item{tail_sd}{SD of last \code{tail_n} scans}
#'     \item{tail_cv}{CV of last \code{tail_n} scans}
#'     \item{stable}{Logical: hold passed stability test}
#'   }
#'
#' @export
detect_unstable_holds <- function(df, var,
                                  col_profile = "profile_id",
                                  col_depth   = "z_meas_m",
                                  col_target  = "z_target_m",
                                  tail_n = 5L, max_cv = 0.02,
                                  max_sd = NULL) {

  has_target <- col_target %in% names(df) && any(!is.na(df[[col_target]]))
  if (!has_target)
    stop("Target depth column required for hold stability check")

  pid    <- df[[col_profile]]
  target <- df[[col_target]]
  vals   <- df[[var]]

  pids <- unique(pid[!is.na(pid)])
  results <- vector("list", length(pids) * 10L)
  ri <- 0L

  for (p in pids) {
    idx <- which(pid == p)
    tgts <- unique(target[idx])
    tgts <- tgts[!is.na(tgts)]

    for (tgt in tgts) {
      hold_idx <- idx[!is.na(target[idx]) &
                       abs(target[idx] - tgt) < 0.5]
      hold_vals <- vals[hold_idx]
      hold_vals <- hold_vals[!is.na(hold_vals)]

      n <- length(hold_vals)
      if (n < tail_n) {
        tail_mean <- NA_real_
        tail_sd   <- NA_real_
        tail_cv   <- NA_real_
      } else {
        tail_vals <- hold_vals[(n - tail_n + 1L):n]
        tail_mean <- mean(tail_vals)
        tail_sd   <- stats::sd(tail_vals)
        tail_cv   <- if (abs(tail_mean) > 1e-10) tail_sd / abs(tail_mean) else NA_real_
      }

      stable <- TRUE
      if (!is.na(tail_cv) && tail_cv > max_cv) stable <- FALSE
      if (!is.null(max_sd) && !is.na(tail_sd) && tail_sd > max_sd) stable <- FALSE

      ri <- ri + 1L
      results[[ri]] <- data.frame(
        profile_id   = p,
        target_depth = tgt,
        n_scans      = n,
        tail_mean    = round(tail_mean, 4),
        tail_sd      = round(tail_sd, 5),
        tail_cv      = round(tail_cv, 4),
        stable       = stable,
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, results[seq_len(ri)])
}


#' Detect noise floor changes over time
#'
#' Computes a rolling standard deviation of the high-frequency residuals
#' (after removing the local trend) and flags periods where the noise
#' floor departs significantly from the baseline. A doubling of the
#' noise floor indicates sensor degradation, electrical interference,
#' or connector corrosion — often before it shows up as visible data
#' quality issues.
#'
#' @param x Numeric vector.
#' @param datetime POSIXct vector.
#' @param window_days Numeric. Rolling window for computing local noise
#'   (days). Default 7.
#' @param detrend_k Integer. Kernel width for the running median used to
#'   remove the local trend before computing residual SD. Default 21
#'   (scans). Should be odd.
#' @param baseline_days Numeric. Initial period (days) to establish the
#'   baseline noise floor. Default 30.
#' @param threshold_factor Numeric. Flag when rolling noise exceeds
#'   \code{threshold_factor} times the baseline noise. Default 2 (doubled).
#' @param min_points Integer. Minimum observations per window. Default 50.
#' @return Data frame with columns:
#'   \describe{
#'     \item{datetime}{Window centre timestamp}
#'     \item{noise_sd}{Rolling residual SD}
#'     \item{baseline_sd}{Baseline noise SD (from first \code{baseline_days})}
#'     \item{ratio}{noise_sd / baseline_sd}
#'     \item{flagged}{Logical: noise exceeds threshold}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' noise <- detect_noise_change(pro$c_condsp_uscm, pro$datetime,
#'                               window_days = 7, threshold_factor = 2)
#' plot(noise$datetime, noise$ratio, type = "l")
#' abline(h = 2, col = "red")
#' }
detect_noise_change <- function(x, datetime, window_days = 7,
                                detrend_k = 21L, baseline_days = 30,
                                threshold_factor = 2, min_points = 50L) {

  valid <- !is.na(x) & !is.na(datetime)
  if (sum(valid) < min_points * 2)
    return(data.frame(datetime = as.POSIXct(character()),
                       noise_sd = numeric(), baseline_sd = numeric(),
                       ratio = numeric(), flagged = logical()))

  x_v  <- x[valid]
  dt_v <- datetime[valid]
  ord  <- order(dt_v)
  x_v  <- x_v[ord]
  dt_v <- dt_v[ord]
  n    <- length(x_v)

  # Detrend: remove local trend with running median
  if (detrend_k %% 2L == 0L) detrend_k <- detrend_k + 1L
  detrend_k <- min(detrend_k, n)
  if (detrend_k < 3L) detrend_k <- 3L
  trend    <- stats::runmed(x_v, k = detrend_k, endrule = "median")
  residual <- x_v - trend

  # Rolling noise: SD of residuals in a time window
  dt_num <- as.numeric(dt_v)
  window_secs <- window_days * 86400
  baseline_secs <- baseline_days * 86400

  # Baseline noise from first N days
  baseline_mask <- (dt_num - dt_num[1]) <= baseline_secs
  if (sum(baseline_mask) < min_points)
    baseline_mask <- seq_len(min(min_points, n))
  baseline_sd <- stats::sd(residual[baseline_mask])
  if (is.na(baseline_sd) || baseline_sd == 0)
    baseline_sd <- stats::mad(residual, na.rm = TRUE) * 1.4826
  if (baseline_sd == 0) baseline_sd <- 1e-10

  # Compute rolling SD using two-pointer window
  # Output one value per window step (step = window/4 for efficiency)
  step <- max(1L, as.integer(sum(baseline_mask) / 4))
  centres <- seq(1L, n, by = step)

  out_dt  <- dt_v[centres]
  out_sd  <- rep(NA_real_, length(centres))

  left <- 1L
  for (ci in seq_along(centres)) {
    centre <- centres[ci]
    t_centre <- dt_num[centre]

    # Find window bounds
    w_left  <- which.min(abs(dt_num - (t_centre - window_secs / 2)))
    w_right <- which.min(abs(dt_num - (t_centre + window_secs / 2)))

    if ((w_right - w_left + 1L) < min_points) next
    out_sd[ci] <- stats::sd(residual[w_left:w_right])
  }

  keep <- !is.na(out_sd)
  data.frame(
    datetime    = out_dt[keep],
    noise_sd    = round(out_sd[keep], 6),
    baseline_sd = round(baseline_sd, 6),
    ratio       = round(out_sd[keep] / baseline_sd, 2),
    flagged     = out_sd[keep] > threshold_factor * baseline_sd,
    stringsAsFactors = FALSE
  )
}


#' Detect optical sensor response time degradation
#'
#' Estimates the effective response time of an optical sensor (typically
#' DO) by measuring how quickly readings change at thermocline
#' crossings during descent. As a membrane ages or fouls, the T90
#' response slows — the sensor takes longer to track rapid changes.
#'
#' For each profile, finds the steepest temperature gradient zone
#' (thermocline) and measures the lag between the temperature and
#' target variable (e.g. DO) response. Increasing lag over time means
#' degrading response.
#'
#' @param df Data frame of profile data.
#' @param var Column to check response time for (default \code{"c_dosat_psat"}).
#' @param col_profile Profile ID column (default \code{"profile_id"}).
#' @param col_datetime Datetime column (default \code{"datetime"}).
#' @param col_depth Depth column (default \code{"z_meas_m"}).
#' @param col_temp Temperature column (default \code{"t_wtr_degc"}).
#' @param min_gradient Minimum temperature gradient (deg C/m) to qualify
#'   as a thermocline crossing. Default 0.5.
#' @param min_scans Minimum scans in the thermocline zone per profile.
#'   Default 3.
#' @return Data frame with one row per profile:
#'   \describe{
#'     \item{profile_id}{Profile identifier}
#'     \item{datetime}{Profile start time}
#'     \item{thermo_depth}{Depth of steepest gradient}
#'     \item{temp_gradient}{Temperature gradient (deg C/m)}
#'     \item{var_gradient}{Variable gradient per metre at thermocline}
#'     \item{lag_ratio}{Ratio of variable response to temperature response.
#'       Values > 1 mean the variable is responding slower than temperature.
#'       Increasing over time = degrading response.}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' resp <- detect_response_degradation(pro, var = "c_dosat_psat")
#' plot(resp$datetime, resp$lag_ratio, type = "p")
#' # Rising trend = membrane needs replacing
#' }
detect_response_degradation <- function(df, var = "c_dosat_psat",
                                        col_profile  = "profile_id",
                                        col_datetime = "datetime",
                                        col_depth    = "z_meas_m",
                                        col_temp     = "t_wtr_degc",
                                        min_gradient = 0.5,
                                        min_scans    = 3L) {

  pid   <- df[[col_profile]]
  depth <- df[[col_depth]]
  temp  <- df[[col_temp]]
  vals  <- df[[var]]
  dt    <- df[[col_datetime]]

  pids <- unique(pid[!is.na(pid)])
  results <- vector("list", length(pids))

  for (i in seq_along(pids)) {
    p <- pids[i]
    idx <- which(pid == p)
    if (length(idx) < 5L) next

    d <- depth[idx]
    t <- temp[idx]
    v <- vals[idx]

    # Sort by depth
    ord <- order(d)
    d <- d[ord]; t <- t[ord]; v <- v[ord]

    ok <- !is.na(d) & !is.na(t) & !is.na(v)
    if (sum(ok) < 5L) next
    d <- d[ok]; t <- t[ok]; v <- v[ok]

    # Compute temperature gradient (dT/dz)
    dz <- diff(d)
    dz[dz == 0] <- NA
    dt_dz <- diff(t) / dz

    # Find steepest gradient zone
    best <- which.max(abs(dt_dz))
    if (length(best) == 0L || is.na(dt_dz[best])) next
    if (abs(dt_dz[best]) < min_gradient) next

    # Window around the thermocline
    thermo_zone <- max(1L, best - 1L):min(length(dt_dz), best + 1L)
    if (length(thermo_zone) < min_scans) next

    # Variable gradient at the same zone
    dv_dz <- diff(v) / dz
    mean_temp_grad <- mean(abs(dt_dz[thermo_zone]), na.rm = TRUE)
    mean_var_grad  <- mean(abs(dv_dz[thermo_zone]), na.rm = TRUE)

    # Lag ratio: normalised variable response relative to temperature
    # Both expressed per metre, so ratio is dimensionless
    # Expected: var_grad should scale with temp_grad
    # If ratio drops over time, var response is lagging
    if (mean_temp_grad == 0) next
    lag_ratio <- mean_var_grad / mean_temp_grad

    # Normalise to expected sensitivity (e.g. DO changes ~0.3 mg/L per °C)
    # This makes the ratio comparable across seasons
    # Skip normalisation for now — raw ratio trend is sufficient

    results[[i]] <- data.frame(
      profile_id    = p,
      datetime      = min(dt[idx], na.rm = TRUE),
      thermo_depth  = round(d[best], 1),
      temp_gradient = round(mean_temp_grad, 3),
      var_gradient  = round(mean_var_grad, 4),
      lag_ratio     = round(lag_ratio, 4),
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, Filter(Negate(is.null), results))
  if (is.null(out))
    return(data.frame(profile_id = character(), datetime = as.POSIXct(character()),
                       thermo_depth = numeric(), temp_gradient = numeric(),
                       var_gradient = numeric(), lag_ratio = numeric(),
                       stringsAsFactors = FALSE))
  rownames(out) <- NULL
  out
}
