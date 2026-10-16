library(data.table)

#' Detect instrument type periods for each lake
#'
#' Classifies instrument as Vaisala (wind in km/h, rain in mm/timestamp) or
#' Lufft (wind in m/s, rain in mm/hr) based on monthly median wind speed.
#' Detects step changes where the median shifts by >factor of 2.
#'
#' @param rds_dir character, path to RDS directory
#' @param ws_threshold numeric, wind speed threshold (m/s vs km/h boundary).
#'   Monthly medians above this are classified as Vaisala (default 4.5)
#' @param change_factor numeric, ratio threshold to detect step changes (default 1.8)
#' @return data.table with columns: lake, date_from, date_to, instrument,
#'   median_ws, note
detect_instruments <- function(rds_dir = "data/rds", ws_threshold = 4.5,
                                change_factor = 1.8) {
  files <- list.files(rds_dir, pattern = "^full_.*[.]rds$", full.names = TRUE)

  results <- rbindlist(lapply(files, function(f) {
    lk <- sub("full_", "", sub("[.]rds$", "", basename(f)))
    x <- readRDS(f)
    met <- as.data.table(x[["raw"]][["met"]])
    cat("Processing", lk, "...\n")

    # Monthly median wind speed
    met[, ym := format(DateTime, "%Y-%m")]
    monthly <- met[!is.na(WndSpd) & WndSpd > 0,
                   .(med_ws = median(WndSpd),
                     date_from = min(DateTime),
                     date_to = max(DateTime)),
                   by = ym]
    setorder(monthly, ym)

    if (nrow(monthly) == 0) {
      return(data.table(lake = lk, date_from = as.Date(NA), date_to = as.Date(NA),
                        instrument = "unknown", median_ws = NA_real_, note = "no data"))
    }

    # Detect breakpoints: where ratio of consecutive medians exceeds threshold
    monthly[, prev_ws := shift(med_ws)]
    monthly[, ratio := med_ws / prev_ws]
    break_idx <- which(!is.na(monthly$ratio) &
                         (monthly$ratio > change_factor | monthly$ratio < 1/change_factor))

    # Build periods from breakpoints
    starts <- c(1, break_idx)
    ends <- c(break_idx - 1, nrow(monthly))

    periods <- rbindlist(lapply(seq_along(starts), function(i) {
      rows <- monthly[starts[i]:ends[i]]
      med <- median(rows$med_ws)
      inst <- if (med > ws_threshold) "Vaisala" else "Lufft"
      note_parts <- character(0)
      if (i > 1) {
        note_parts <- c(note_parts, sprintf("step change (ratio=%.1f)",
                                             monthly$ratio[break_idx[i-1]]))
      }
      data.table(
        lake       = lk,
        date_from  = as.Date(min(rows$date_from)),
        date_to    = as.Date(max(rows$date_to)),
        instrument = inst,
        median_ws  = round(med, 1),
        note       = paste(note_parts, collapse = "; ")
      )
    }))
    periods
  }), fill = TRUE)

  # Merge consecutive periods with same instrument type
  results[, grp := rleid(lake, instrument)]
  merged <- results[, .(
    date_from  = min(date_from),
    date_to    = max(date_to),
    instrument = instrument[1],
    median_ws  = round(median(median_ws), 1),
    note       = paste(note[note != ""], collapse = "; ")
  ), by = .(lake, grp)]
  merged[, grp := NULL]

  setorder(merged, lake, date_from)
  merged
}
