library(data.table)
library(MeTo)

#' Calculate solar noon for each date given site longitude
#'
#' @param dates Date vector or POSIXct (unique dates extracted)
#' @param lon numeric, site longitude in degrees East
#' @param tz character, timezone string (default "Etc/GMT-12" for NZST)
#' @return data.table with columns: date, solar_noon (POSIXct)
calc_solar_noon <- function(dates, lon, tz = "Etc/GMT-12") {
  Lz <- 180  # centre longitude of UTC+12
  dates <- as.Date(dates)
  dates <- sort(unique(dates))
  # Solar noon offset from 12:00 local time (hours)
  # Positive lon west of Lz => noon is later; east => earlier
  offset_hrs <- (Lz - lon) / 15
  noon_posix <- as.POSIXct(paste(dates, "12:00:00"), tz = tz) +
    offset_hrs * 3600
  data.table(date = dates, solar_noon = noon_posix)
}

#' Calculate clear-sky radiation for a site at given timestamps
#'
#' Uses MeTo::Rso (FAO56) and converts from MJ/m2/period to W/m2.
#' Builds a lookup table for unique (day-of-year, hour, minute) combinations
#' so Rso() is only called once per unique time-of-year slot, then maps
#' results back to all timestamps. This is much faster for multi-year data.
#'
#' @param datetime POSIXct vector of timestamps
#' @param lat numeric, latitude in degrees (negative for southern hemisphere)
#' @param lon numeric, longitude in degrees East
#' @param elev numeric, elevation in metres above sea level
#' @param tl_hours numeric, calculation period length in hours (default 1)
#' @return data.table with columns: DateTime, clearsky_wm2
calc_clearsky <- function(datetime, lat, lon, elev, tl_hours = 1) {
  dt <- data.table(DateTime = datetime)
  dt[, `:=`(
    doy  = as.integer(format(DateTime, "%j")),
    hour = as.integer(format(DateTime, "%H")),
    min  = as.integer(format(DateTime, "%M"))
  )]
  dt[, slot := paste(doy, hour, min, sep = "-")]

  # Build a 30-min lookup grid for one reference year, then interpolate
  tz <- attr(datetime[1], "tzone")
  if (is.null(tz) || tz == "") tz <- "Etc/GMT-12"
  ref_year <- as.integer(format(datetime[1], "%Y"))

  unique_doys <- sort(unique(dt$doy))
  min_hour <- min(dt$hour)
  max_hour <- max(dt$hour) + 1

  # 30-min grid: each DOY × half-hour slots across the window
  half_hours <- seq(min_hour, max_hour, by = 0.5)
  grid <- CJ(doy = unique_doys, fhour = half_hours)

  ref_dates <- as.Date(grid$doy - 1, origin = paste0(ref_year, "-01-01"))
  ref_times <- as.POSIXct(
    paste(ref_dates, sprintf("%02d:%02d:00",
                             as.integer(grid$fhour),
                             as.integer((grid$fhour %% 1) * 60))),
    tz = tz
  )

  cat("  Rso lookup grid:", nrow(grid), "points\n")

  rso_mj <- vapply(ref_times, function(t) {
    MeTo::Rso(
      x        = t,
      tl       = tl_hours,
      elev     = elev,
      lat.deg  = lat,
      long.deg = lon,
      control  = list(Lz = 180)
    )
  }, numeric(1))

  grid[, clearsky_wm2 := pmax(rso_mj * 1e6 / (tl_hours * 3600), 0)]

  # Interpolate to actual timestamps
  dt[, fhour := hour + min / 60]
  dt[, clearsky_wm2 := {
    g <- grid[doy == .BY$doy]
    approx(g$fhour, g$clearsky_wm2, xout = fhour, rule = 2)$y
  }, by = doy]

  dt[clearsky_wm2 < 0, clearsky_wm2 := 0]
  setorder(dt, DateTime)
  dt[, .(DateTime, clearsky_wm2)]
}

#' Calculate clear-sky radiation for the solar noon window across all dates
#'
#' Computes clearsky at the native met data timestamps, but only within
#' the ±hours_window around solar noon to avoid unnecessary computation.
#'
#' @param met_dt data.table with DateTime column
#' @param lat,lon,elev site coordinates and elevation
#' @param solar_noon_dt data.table from calc_solar_noon() with date, solar_noon
#' @param hours_window numeric, hours either side of solar noon (default 3)
#' @param tl_hours numeric, Rso calculation period (default 1)
#' @return data.table with DateTime, clearsky_wm2 (only rows within window)
calc_clearsky_window <- function(met_dt, lat, lon, elev, solar_noon_dt,
                                  hours_window = 3, tl_hours = 1) {
  dt <- copy(met_dt[, .(DateTime)])
  dt[, date := as.Date(DateTime)]
  dt <- merge(dt, solar_noon_dt, by = "date", all.x = FALSE)
  dt[, hrs_from_noon := as.numeric(difftime(DateTime, solar_noon, units = "hours"))]
  dt <- dt[abs(hrs_from_noon) <= hours_window]

  if (nrow(dt) == 0) return(data.table(DateTime = as.POSIXct(character(0)), clearsky_wm2 = numeric(0)))

  cat("  Computing clear-sky for", nrow(dt), "timestamps...\n")
  cs <- calc_clearsky(dt$DateTime, lat, lon, elev, tl_hours)
  cs
}
