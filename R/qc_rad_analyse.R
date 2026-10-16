#' Analyse offset and drift in solar radiation sensor
#'
#' For each summer (solstice window), selects up to two cloudless days closest
#' to the solstice and computes the ratio of observed to clear-sky peak
#' radiation. Drift is assessed by fitting a linear trend to the ratio over
#' time.
#'
#' @param met Data frame with `DateTime`, `RadSWD`, `RadClr`, `Date`, and
#'   `ClearDay` columns. Run [qc_rad_clear_days()] first.
#' @param solstice_window Number of days either side of the summer solstice
#'   to include. Default 60 (~2 months). Set to `Inf` to use all cloudless
#'   days.
#' @param max_per_summer Maximum number of days to keep per summer year.
#'   The closest-to-solstice days are preferred. Default 2.
#' @param plot Logical; if `TRUE` (default), plot the ratio time series with
#'   trend.
#'
#' @return A list with components:
#'   \describe{
#'     \item{ratios}{Data frame with `Date`, `obs_max`, `clr_max`, `ratio`,
#'       and `fitted` for each cloudless day used.}
#'     \item{mean_ratio}{Mean ratio across used cloudless days.}
#'     \item{trend}{Linear model (`lm`) of ratio ~ date_num.}
#'     \item{stable}{Logical; `TRUE` if the range of the fitted trend is
#'       less than 0.10 (i.e., drift < 10\% of clear-sky).}
#'   }
#' @export
qc_rad_analyse <- function(met, solstice_window = 60, max_per_summer = 2,
                           plot = TRUE) {
  clear_dates <- unique(met$Date[met$ClearDay])

  # Determine hemisphere from peak clear-sky DOY
  clr_by_date <- tapply(met$RadClr[met$ClearDay], met$Date[met$ClearDay],
                        max, na.rm = TRUE)
  peak_doy <- as.numeric(format(
    as.Date(names(which.max(clr_by_date))), "%j"))
  southern <- peak_doy > 300 || peak_doy < 80

  if (southern) {
    solstice_doy <- 355
  } else {
    solstice_doy <- 172
  }

  # Compute distance to solstice for each clear date
  doy_all <- as.numeric(format(clear_dates, "%j"))
  if (southern) {
    dist <- pmin(abs(doy_all - solstice_doy),
                 365 - abs(doy_all - solstice_doy))
  } else {
    dist <- abs(doy_all - solstice_doy)
  }

  # Assign summer year
  all_months <- as.integer(format(clear_dates, "%m"))
  all_years <- as.integer(format(clear_dates, "%Y"))
  if (southern) {
    summer_year <- ifelse(all_months >= 7, all_years + 1L, all_years)
  } else {
    summer_year <- all_years
  }

  # Filter to solstice window
  if (is.finite(solstice_window)) {
    in_window <- dist <= solstice_window
  } else {
    in_window <- rep(TRUE, length(clear_dates))
  }

  # Check data coverage per summer year
  met_dates <- unique(met$Date)
  met_doy <- as.numeric(format(met_dates, "%j"))
  met_months <- as.integer(format(met_dates, "%m"))
  met_years <- as.integer(format(met_dates, "%Y"))
  if (southern) {
    met_summer_year <- ifelse(met_months >= 7, met_years + 1L, met_years)
    met_dist <- pmin(abs(met_doy - solstice_doy),
                     365 - abs(met_doy - solstice_doy))
  } else {
    met_summer_year <- met_years
    met_dist <- abs(met_doy - solstice_doy)
  }
  met_in_window <- met_dist <= solstice_window
  window_days_per_sy <- tapply(met_in_window, met_summer_year, sum)
  expected_window_days <- 2 * solstice_window

  # Build candidate table
  candidates <- data.frame(
    date = clear_dates,
    dist = dist,
    summer_year = summer_year,
    in_window = in_window,
    stringsAsFactors = FALSE
  )

  # For each summer year: pick up to max_per_summer days closest to solstice
  selected <- character(0)
  for (sy in sort(unique(summer_year))) {
    sy_cand <- candidates[candidates$summer_year == sy, ]

    # First try days within the solstice window
    sy_in <- sy_cand[sy_cand$in_window, ]

    if (nrow(sy_in) > 0) {
      # Pick the closest to solstice
      sy_in <- sy_in[order(sy_in$dist), ]
      n_pick <- min(max_per_summer, nrow(sy_in))
      selected <- c(selected, as.character(sy_in$date[1:n_pick]))
    } else {
      # No days in window — force-include best day only if >80% coverage
      sy_coverage <- window_days_per_sy[as.character(sy)]
      if (!is.na(sy_coverage) && sy_coverage / expected_window_days >= 0.8) {
        # Pick the single closest day from outside the window
        sy_cand <- sy_cand[order(sy_cand$dist), ]
        selected <- c(selected, as.character(sy_cand$date[1]))
      }
    }
  }

  clear_dates <- sort(as.Date(unique(selected)))

  if (length(clear_dates) < 2) {
    stop("Need at least 2 cloudless days near the solstice for analysis. ",
         "Try increasing solstice_window or relaxing max_roc_rmse in ",
         "qc_rad_clear_days().")
  }

  # For each cloudless day, get max observed and max clearsky
  ratios <- data.frame(
    Date = as.Date(character(0)),
    obs_max = numeric(0),
    clr_max = numeric(0),
    ratio = numeric(0)
  )

  for (d in as.character(clear_dates)) {
    idx <- met$Date == as.Date(d) & met$ClearDay
    obs_max <- max(met$RadSWD[idx], na.rm = TRUE)
    clr_max <- max(met$RadClr[idx], na.rm = TRUE)
    if (clr_max > 0) {
      ratios <- rbind(ratios, data.frame(
        Date = as.Date(d),
        obs_max = obs_max,
        clr_max = clr_max,
        ratio = obs_max / clr_max
      ))
    }
  }

  # Remove outlier ratios using MAD (robust to single extreme values)
  med_ratio <- median(ratios$ratio)
  mad_ratio <- mad(ratios$ratio)
  if (mad_ratio > 0) {
    keep <- abs(ratios$ratio - med_ratio) < 3 * mad_ratio
    n_removed <- sum(!keep)
    if (n_removed > 0) {
      message("Removed ", n_removed, " outlier day(s): ",
              paste(ratios$Date[!keep], collapse = ", "))
    }
    ratios <- ratios[keep, ]
  }

  if (nrow(ratios) < 2) {
    stop("Fewer than 2 days remain after outlier removal.")
  }

  # Linear trend (simple, avoids overfitting)
  ratios$date_num <- as.numeric(ratios$Date - ratios$Date[1])
  trend <- lm(ratio ~ date_num, data = ratios)
  trend_vals <- predict(trend)

  ratios$fitted <- trend_vals
  drift_range <- diff(range(trend_vals))

  result <- list(
    ratios = ratios[, c("Date", "obs_max", "clr_max", "ratio", "fitted")],
    mean_ratio = mean(ratios$ratio),
    trend = trend,
    stable = drift_range < 0.10
  )

  if (plot) {
    plot(ratios$Date, ratios$ratio, pch = 16, col = "steelblue",
         xlab = "Date", ylab = "Obs / Clear-sky ratio",
         main = "Radiation Sensor Offset & Drift (solstice window)")
    abline(h = 1, col = "grey50", lty = 2)
    pred_dates <- seq(min(ratios$Date), max(ratios$Date), by = "day")
    pred_num <- as.numeric(pred_dates - ratios$Date[1])
    pred_trend <- predict(trend, newdata = data.frame(date_num = pred_num))
    lines(pred_dates, pred_trend, col = "red", lwd = 2)
    legend("topright",
           legend = c(
             sprintf("Mean ratio: %.3f", result$mean_ratio),
             sprintf("Linear trend (%.4f/yr)", coef(trend)[2] * 365.25)),
           pch = c(16, NA), lty = c(NA, 1),
           col = c("steelblue", "red"), bty = "n", cex = 0.8)
  }

  result
}
