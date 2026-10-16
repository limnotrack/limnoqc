#' Identify cloudless days from radiation time series
#'
#' A day is classified as cloudless when the shape of the observed daytime
#' radiation closely matches the clear-sky model. Specifically, the observed
#' rate-of-change (successive differences) must closely track the clear-sky
#' rate-of-change for the same day.
#'
#' Only days with complete data resolution (no missing timestamps) are
#' considered. Additionally, at least one day per season per year is guaranteed
#' by selecting the best-scoring day for any season-year that has no days
#' passing the RMSE threshold.
#'
#' @param met Data frame with `DateTime`, `RadSWD`, and `RadClr` columns.
#'   Run [qc_rad_clearsky()] first to add `RadClr`.
#' @param min_rad Minimum radiation (W/m2) to consider a timestep as
#'   "daytime". Default 20.
#' @param max_roc_rmse Maximum RMSE between observed and clear-sky
#'   rate-of-change (normalised by clear-sky peak). Default 0.012.
#' @param min_ratio Minimum ratio of observed peak to clear-sky peak for a
#'   day to be considered. Filters out stuck/failed sensors. Default 0.4.
#'
#' @return `met` with additional columns:
#'   \describe{
#'     \item{Date}{Date (Date class) extracted from `DateTime`.}
#'     \item{ClearDay}{Logical; `TRUE` for timesteps belonging to a cloudless
#'       day.}
#'   }
#' @export
qc_rad_clear_days <- function(met, min_rad = 20, max_roc_rmse = 0.012,
                              min_ratio = 0.4) {
  met$Date <- as.Date(format(met$DateTime, "%Y-%m-%d"))
  dates <- unique(met$Date)

  # Detect expected records per day from modal timestep
  diffs <- as.numeric(diff(met$DateTime), units = "mins")
  dt_min <- as.numeric(names(sort(table(diffs), decreasing = TRUE)[1]))
  expected_per_day <- as.integer(round(1440 / dt_min))

  # Score every candidate day
  day_scores <- data.frame(date = character(0), roc_rmse = numeric(0),
                           stringsAsFactors = FALSE)

  for (d in as.character(dates)) {
    idx <- which(met$Date == as.Date(d))
    day <- met[idx, ]

    # Full resolution check: require exactly the expected number of records
    if (nrow(day) != expected_per_day) next

    # Daytime: where clear-sky > min_rad
    daytime <- day$RadClr > min_rad
    if (sum(daytime, na.rm = TRUE) < 10) next

    obs <- day$RadSWD[daytime]
    clr <- day$RadClr[daytime]
    if (any(is.na(obs))) next

    # Peak ratio check
    obs_peak <- max(obs, na.rm = TRUE)
    clr_peak <- max(clr, na.rm = TRUE)
    if (obs_peak < min_rad) next
    if (obs_peak / clr_peak < min_ratio) next

    # Compare rate-of-change: normalise each series by its own peak so
    # shape comparison is independent of sensor amplitude/degradation
    roc_obs <- diff(obs / obs_peak)
    roc_clr <- diff(clr / clr_peak)
    roc_rmse <- sqrt(mean((roc_obs - roc_clr)^2))

    day_scores <- rbind(day_scores, data.frame(date = d, roc_rmse = roc_rmse,
                                               stringsAsFactors = FALSE))
  }

  if (nrow(day_scores) == 0) {
    met$ClearDay <- FALSE
    return(met)
  }

  # Days passing the threshold
  clear_dates <- day_scores$date[day_scores$roc_rmse < max_roc_rmse]

  # Assign season-year to every scored day
  scored_dates <- as.Date(day_scores$date)
  months <- as.integer(format(scored_dates, "%m"))
  years <- as.integer(format(scored_dates, "%Y"))
  # Southern hemisphere seasons: DJF=summer, MAM=autumn, JJA=winter, SON=spring
  season <- ifelse(months %in% c(12, 1, 2), "DJF",
            ifelse(months %in% c(3, 4, 5), "MAM",
            ifelse(months %in% c(6, 7, 8), "JJA", "SON")))
  # For Dec, assign to next year's summer
  season_year <- ifelse(months == 12, years + 1L, years)
  day_scores$season_year <- paste(season_year, season, sep = "-")

  # Ensure at least one day per season-year
  all_sy <- unique(day_scores$season_year)
  for (sy in all_sy) {
    sy_rows <- day_scores[day_scores$season_year == sy, ]
    # Check if any day from this season-year already passed
    if (any(sy_rows$date %in% clear_dates)) next
    # Add the best (lowest RMSE) day from this season-year
    best <- sy_rows$date[which.min(sy_rows$roc_rmse)]
    clear_dates <- c(clear_dates, best)
  }

  clear_dates <- unique(clear_dates)

  met$ClearDay <- met$Date %in% as.Date(clear_dates)
  met
}
