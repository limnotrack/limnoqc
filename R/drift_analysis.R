library(data.table)

#' Load lake data from an RDS file
#'
#' @param lake_name character, e.g. "Rotoiti"
#' @param rds_dir character, path to directory containing full_*.rds files
#' @return list with elements: lake (name), lat, lon, elev, met (data.table)
load_lake_data <- function(lake_name, rds_dir = "data/rds") {
  path <- file.path(rds_dir, paste0("full_", lake_name, ".rds"))
  x <- readRDS(path)

  meta <- x[["meta"]]
  get_meta <- function(attr) as.numeric(meta$value[meta$attribute == attr])

  met <- as.data.table(x[["raw"]][["met"]])

  list(
    lake = lake_name,
    lat  = get_meta("lat"),
    lon  = get_meta("lon"),
    elev = get_meta("masl"),
    met  = met
  )
}

#' Compute drift metrics for selected clear days
#'
#' For each clear day, computes mean observed and clear-sky radiation in the
#' solar noon window, and their ratio.
#'
#' @param met_dt data.table with DateTime and RadSWD
#' @param clearsky_dt data.table with DateTime and clearsky_wm2
#' @param clear_days data.table from select_clear_days() with date column
#' @param solar_noon_dt data.table with date and solar_noon
#' @param hours_window numeric, hours either side of solar noon
#' @return data.table with: date, obs_mean, clearsky_mean, ratio, bias, clarity_score
calc_drift <- function(met_dt, clearsky_dt, clear_days, solar_noon_dt,
                        hours_window = 3) {

  merged <- merge(
    met_dt[, .(DateTime, RadSWD)],
    clearsky_dt[, .(DateTime, clearsky_wm2)],
    by = "DateTime"
  )
  merged[, date := as.Date(DateTime)]
  merged <- merge(merged, solar_noon_dt, by = "date")
  merged[, hrs_from_noon := as.numeric(difftime(DateTime, solar_noon, units = "hours"))]
  merged <- merged[abs(hrs_from_noon) <= hours_window]
  merged <- merged[!is.na(RadSWD) & RadSWD >= 5]

  # Keep only selected clear days
  merged <- merged[date %in% clear_days$date]

  if (nrow(merged) == 0) return(data.table())

  result <- merged[, .(
    obs_mean      = mean(RadSWD, na.rm = TRUE),
    clearsky_mean = mean(clearsky_wm2, na.rm = TRUE)
  ), by = date]

  result[, ratio := obs_mean / clearsky_mean]
  result[, bias := obs_mean - clearsky_mean]

  # Attach clarity scores
  result <- merge(result, clear_days[, .(date, clarity_score)], by = "date")
  setorder(result, date)
  result
}

#' Run drift analysis for a single lake
#'
#' @param lake_name character
#' @param rds_dir character
#' @param hours_window numeric
#' @param n_per_month integer
#' @param clarity_threshold numeric or NULL
#' @return data.table with drift results including lake column
run_lake <- function(lake_name, rds_dir = "data/rds", hours_window = 2,
                      n_per_month = 1, clarity_threshold = 0.01) {
  cat("Processing", lake_name, "...\n")

  ld <- load_lake_data(lake_name, rds_dir)
  cat("  Met data:", nrow(ld$met), "rows,",
      format(range(ld$met$DateTime, na.rm = TRUE), "%Y-%m-%d"), "\n")

  solar_noon_dt <- calc_solar_noon(ld$met$DateTime, ld$lon)
  clearsky_dt <- calc_clearsky_window(ld$met, ld$lat, ld$lon, ld$elev,
                                       solar_noon_dt, hours_window)

  clear_days <- select_clear_days(ld$met, clearsky_dt, solar_noon_dt,
                                   hours_window, n_per_month, clarity_threshold)

  if (nrow(clear_days) == 0) {
    cat("  No clear days found.\n")
    return(data.table())
  }

  cat("  Selected", nrow(clear_days), "clear days\n")

  drift <- calc_drift(ld$met, clearsky_dt, clear_days, solar_noon_dt, hours_window)
  drift[, lake := lake_name]
  drift
}

#' Run drift analysis across all lakes
#'
#' @param rds_dir character
#' @param hours_window numeric
#' @param n_per_month integer
#' @param clarity_threshold numeric or NULL
#' @return data.table with drift results for all lakes
run_all_lakes <- function(rds_dir = "data/rds", hours_window = 2,
                           n_per_month = 1, clarity_threshold = 0.01) {
  files <- list.files(rds_dir, pattern = "^full_.*\\.rds$", full.names = FALSE)
  lake_names <- sub("^full_", "", sub("\\.rds$", "", files))

  results <- rbindlist(lapply(lake_names, function(lk) {
    tryCatch(
      run_lake(lk, rds_dir, hours_window, n_per_month, clarity_threshold),
      error = function(e) {
        cat("  ERROR:", conditionMessage(e), "\n")
        data.table()
      }
    )
  }), fill = TRUE)

  results
}
