#' Solar radiation QC: clear-sky comparison and drift analysis
#'
#' Runs the full radiation QC pipeline: compute theoretical clear-sky
#' radiation, identify cloudless days, analyse sensor drift, and flag
#' exceedances. Returns structured results for the app and scripts.
#'
#' @param met Data frame with datetime and radiation columns.
#' @param lat Latitude in decimal degrees (negative for southern hemisphere).
#' @param lng Longitude in decimal degrees.
#' @param elev Elevation in metres above sea level. Default 0.
#' @param rad_col Radiation column name. If \code{NULL}, auto-detects from
#'   \code{r_swd_wm2}, \code{r_swd}, \code{RadSWD}.
#' @param dt_col Datetime column name. Default \code{"datetime"}.
#' @param Lz Longitude of the centre of the local time zone (degrees west
#'   of Greenwich). Default 180 for NZST (UTC+12).
#' @param max_roc_sd Maximum normalised rate-of-change SD for cloudless
#'   day detection. Default 0.02.
#' @param exceedance Ratio above which \code{exceeds_clearsky} is
#'   \code{TRUE}. Default 1.05 (5% above theoretical).
#' @param plot Logical; produce diagnostic plots? Default \code{FALSE}.
#' @param correct Logical; apply correction to the radiation column?
#'   Default \code{FALSE}.
#' @param correct_method Correction method passed to [qc_rad_correct()].
#'   Default \code{"auto"}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{clear_days}{Data frame, one row per clear day: \code{date},
#'       \code{obs_peak}, \code{theo_peak}, \code{ratio},
#'       \code{daily_total_obs}, \code{daily_total_theo}, \code{roc_sd}.}
#'     \item{daily_curves}{Data frame, sub-daily resolution for clear days:
#'       \code{datetime}, \code{observed}, \code{theoretical}, \code{date}.}
#'     \item{all_days}{Data frame, one row per day: \code{date},
#'       \code{is_clear}, \code{obs_peak}, \code{theo_peak}, \code{ratio},
#'       \code{roc_sd}, \code{n_obs}.}
#'     \item{exceeds_clearsky}{Logical vector, same length as input.
#'       \code{TRUE} where observed > theoretical × \code{exceedance}.}
#'     \item{sensor_health}{Data frame: \code{drift_pct_per_year},
#'       \code{n_clear_days}, \code{mean_ratio}, \code{recommendation}.}
#'     \item{met}{The (possibly corrected) met data frame with \code{RadClr},
#'       \code{Date}, \code{ClearDay} columns added.}
#'     \item{analysis}{Output from [qc_rad_analyse()] (trend model, ratios).}
#'   }
#'
#' @export
qc_rad <- function(met, lat, lng, elev = 0,
                   rad_col = NULL, dt_col = "datetime",
                   Lz = 180, max_roc_sd = 0.02,
                   exceedance = 1.05,
                   plot = FALSE, correct = FALSE,
                   correct_method = "auto") {

  # -- Auto-detect radiation column ----------------------------------------
  if (is.null(rad_col)) {
    candidates <- c("r_swd_wm2", "r_swd", "RadSWD", "r_swd_wm2_avg")
    rad_col <- intersect(candidates, names(met))[1]
    if (is.na(rad_col))
      stop("No radiation column found. Supply rad_col explicitly.")
  }

  # -- Normalise column names for internal pipeline ------------------------
  # The sub-functions expect DateTime and RadSWD
  met_internal <- met
  if (dt_col != "DateTime" && dt_col %in% names(met_internal))
    names(met_internal)[names(met_internal) == dt_col] <- "DateTime"
  if (rad_col != "RadSWD" && rad_col %in% names(met_internal))
    names(met_internal)[names(met_internal) == rad_col] <- "RadSWD"

  # -- 1. Clear-sky radiation ----------------------------------------------
  met_internal <- qc_rad_clearsky(met_internal, lat = lat, lng = lng,
                                   elev = elev, Lz = Lz)

  # -- 2. Identify cloudless days ------------------------------------------
  met_internal <- qc_rad_clear_days(met_internal, max_roc_rmse = max_roc_sd)

  # -- 3. Build all_days summary -------------------------------------------
  dates <- unique(met_internal$Date)
  all_days <- do.call(rbind, lapply(dates, function(d) {
    idx <- met_internal$Date == d
    day <- met_internal[idx, ]
    daytime <- !is.na(day$RadClr) & day$RadClr > 20
    obs_dt <- day$RadSWD[daytime]
    clr_dt <- day$RadClr[daytime]
    obs_dt <- obs_dt[!is.na(obs_dt)]
    n_obs <- length(obs_dt)

    obs_peak  <- if (n_obs > 0) max(obs_dt) else NA_real_
    theo_peak <- if (any(daytime)) max(clr_dt, na.rm = TRUE) else NA_real_
    ratio     <- if (!is.na(obs_peak) && !is.na(theo_peak) && theo_peak > 0)
                   obs_peak / theo_peak else NA_real_

    # ROC SD (normalised)
    roc_sd <- NA_real_
    if (n_obs > 5 && obs_peak > 20) {
      roc_obs <- diff(obs_dt / obs_peak)
      roc_sd  <- stats::sd(roc_obs)
    }

    data.frame(date = d, is_clear = any(day$ClearDay),
               obs_peak = obs_peak, theo_peak = theo_peak,
               ratio = ratio, roc_sd = roc_sd, n_obs = as.integer(n_obs),
               stringsAsFactors = FALSE)
  }))

  # -- 4. Build clear_days detail ------------------------------------------
  clear_dates <- all_days$date[all_days$is_clear]
  clear_days <- all_days[all_days$is_clear, , drop = FALSE]

  # Add daily totals (MJ/m²/day)
  dt_sec <- as.numeric(names(sort(table(diff(
    as.numeric(met_internal$DateTime))), decreasing = TRUE))[1])

  clear_days$daily_total_obs  <- NA_real_
  clear_days$daily_total_theo <- NA_real_
  for (i in seq_len(nrow(clear_days))) {
    d <- clear_days$date[i]
    idx <- met_internal$Date == d
    obs <- met_internal$RadSWD[idx]
    clr <- met_internal$RadClr[idx]
    obs[is.na(obs)] <- 0
    clr[is.na(clr)] <- 0
    clear_days$daily_total_obs[i]  <- round(sum(obs) * dt_sec / 1e6, 3)
    clear_days$daily_total_theo[i] <- round(sum(clr) * dt_sec / 1e6, 3)
  }
  clear_days$is_clear <- NULL  # redundant for this table

  # -- 5. Build daily_curves for clear days --------------------------------
  daily_curves <- do.call(rbind, lapply(clear_dates, function(d) {
    idx <- met_internal$Date == d
    data.frame(
      datetime    = met_internal$DateTime[idx],
      observed    = met_internal$RadSWD[idx],
      theoretical = met_internal$RadClr[idx],
      date        = d,
      stringsAsFactors = FALSE
    )
  }))

  # -- 6. Exceeds clearsky (per-row flag) ----------------------------------
  exceeds <- rep(FALSE, nrow(met_internal))
  has_both <- !is.na(met_internal$RadSWD) & !is.na(met_internal$RadClr) &
              met_internal$RadClr > 20
  exceeds[has_both] <- met_internal$RadSWD[has_both] >
                       met_internal$RadClr[has_both] * exceedance

  # -- 7. Analysis (trend) -------------------------------------------------
  analysis <- tryCatch(
    qc_rad_analyse(met_internal, plot = plot),
    error = function(e) {
      message("Radiation trend analysis: ", e$message)
      list(ratios = data.frame(), mean_ratio = NA_real_,
           trend = NULL, stable = NA, drift_per_year = NA_real_)
    }
  )

  # Compute drift_per_year from trend coefficients
  drift_per_year <- NA_real_
  if (!is.null(analysis$trend)) {
    drift_per_year <- as.numeric(stats::coef(analysis$trend)[2]) * 365.25
  }
  analysis$drift_per_year <- drift_per_year

  # -- 8. Sensor health summary --------------------------------------------
  drift_pct <- if (!is.na(drift_per_year) && !is.na(analysis$mean_ratio) &&
                   analysis$mean_ratio > 0)
    round(drift_per_year / analysis$mean_ratio * 100, 1) else NA_real_

  recommendation <- if (is.na(drift_pct)) {
    "insufficient data"
  } else if (abs(drift_pct) < 2 && abs(1 - analysis$mean_ratio) < 0.05) {
    "OK"
  } else if (abs(drift_pct) < 5 && abs(1 - analysis$mean_ratio) < 0.10) {
    "check calibration"
  } else {
    "replace sensor"
  }

  sensor_health <- data.frame(
    drift_pct_per_year = drift_pct,
    n_clear_days       = nrow(clear_days),
    mean_ratio         = round(analysis$mean_ratio, 4),
    recommendation     = recommendation,
    stringsAsFactors   = FALSE
  )

  # -- 9. Correction (optional) --------------------------------------------
  if (correct && !is.null(analysis$trend)) {
    met_internal <- qc_rad_correct(met_internal, analysis,
                                    method = correct_method)
  }

  # -- 10. Plots (optional) ------------------------------------------------
  if (plot) {
    tryCatch({
      qc_rad_plot(met_internal, type = "overlay")
      qc_rad_plot(met_internal, type = "daily")
      qc_rad_plot(met_internal, type = "scatter")
    }, error = function(e) message("Plot error: ", e$message))
  }

  # -- Return ---------------------------------------------------------------
  list(
    clear_days        = clear_days,
    daily_curves      = daily_curves,
    all_days          = all_days,
    exceeds_clearsky  = exceeds,
    sensor_health     = sensor_health,
    met               = met_internal,
    analysis          = analysis
  )
}
