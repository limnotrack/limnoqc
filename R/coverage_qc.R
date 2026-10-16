library(data.table)
library(ggplot2)
source("R/palette.R")

#' Compute monthly data coverage for met and profile data
#'
#' For each lake and month, counts the number of days with data in
#' raw$met and raw$pro, plus per-variable coverage for key met and
#' profile variables.
#'
#' @param rds_dir character
#' @return data.table with lake, year_month, data_type, variable, days_with_data
calc_coverage <- function(rds_dir = "data/rds") {
  files <- list.files(rds_dir, pattern = "^full_.*[.]rds$", full.names = TRUE)

  rbindlist(lapply(files, function(f) {
    lk <- sub("full_", "", sub("[.]rds$", "", basename(f)))
    cat("Loading", lk, "...\n")
    x <- readRDS(f)

    results <- list()

    # Met coverage (use raw to see actual data presence)
    # A day counts if >= 10 non-NA TmpAir measurements
    met <- x[["raw"]][["met"]]
    if (!is.null(met)) {
      met <- as.data.table(met)
      met[, `:=`(date = as.Date(format(DateTime, "%Y-%m-%d")),
                 ym   = as.Date(format(DateTime, "%Y-%m-01")))]
      met_vars <- setdiff(names(met), c("DateTime", "date", "ym"))
      for (v in met_vars) {
        day_n <- met[!is.na(get(v)), .N, by = .(ym, date)]
        cov <- day_n[N >= 10, .(days = .N), by = ym]
        cov[, `:=`(lake = lk, data_type = "Met", variable = v)]
        results <- c(results, list(cov))
      }
      # Overall met: days with >= 10 TmpAir readings
      if ("TmpAir" %in% names(met)) {
        day_n <- met[!is.na(TmpAir), .N, by = .(ym, date)]
        cov_all <- day_n[N >= 10, .(days = .N), by = ym]
      } else {
        cov_all <- met[, .(days = uniqueN(date)), by = ym]
      }
      cov_all[, `:=`(lake = lk, data_type = "Met", variable = "ANY")]
      results <- c(results, list(cov_all))
    }

    # Profile coverage (use raw to see actual data presence)
    # A day counts if >= 1 non-NA TmpWtr measurement
    pro <- x[["raw"]][["pro"]]
    if (!is.null(pro)) {
      pro <- as.data.table(pro)
      pro[, `:=`(date = as.Date(format(DateTime, "%Y-%m-%d")),
                 ym   = as.Date(format(DateTime, "%Y-%m-01")))]
      pro_vars <- intersect(names(pro),
        c("TmpWtr", "DOconc", "DOpsat", "pH", "ORP", "Cond", "SpCond",
          "TurbRT", "FlChlr", "FlPhyc"))
      for (v in pro_vars) {
        cov <- pro[!is.na(get(v)), .(days = uniqueN(date)), by = ym]
        cov[, `:=`(lake = lk, data_type = "Profile", variable = v)]
        results <- c(results, list(cov))
      }
      # Overall profile: days with >= 1 TmpWtr reading
      if ("TmpWtr" %in% names(pro)) {
        cov_all <- pro[!is.na(TmpWtr), .(days = uniqueN(date)), by = ym]
      } else {
        cov_all <- pro[, .(days = uniqueN(date)), by = ym]
      }
      cov_all[, `:=`(lake = lk, data_type = "Profile", variable = "ANY")]
      results <- c(results, list(cov_all))
    }

    rbindlist(results, fill = TRUE)
  }), fill = TRUE)
}

#' Convert days to % coverage and fill gaps with 0
#'
#' For each lake × data_type, fills all months from first record to today
#' with 0% where no data exist, so gaps show as red rather than white.
#' @param cov_dt data.table from calc_coverage()
#' @return data.table with pct column (0–100)
fill_coverage <- function(cov_dt) {
  dt <- copy(cov_dt)

  # Days in each month
  dt[, days_in_month := as.integer(as.Date(format(ym + 32, "%Y-%m-01")) - ym)]
  dt[, pct := round(100 * pmin(days, days_in_month) / days_in_month)]

  # Fill from first record to today with 0% for missing months
  today_ym <- as.Date(format(Sys.Date(), "%Y-%m-01"))
  filled <- dt[, {
    first_ym <- min(ym)
    all_ym <- seq(first_ym, today_ym, by = "month")
    missing <- setdiff(as.character(all_ym), as.character(ym))
    if (length(missing) > 0) {
      gaps <- data.table(ym = as.Date(missing), days = 0L,
                         days_in_month = as.integer(
                           as.Date(format(as.Date(missing) + 32, "%Y-%m-01")) -
                           as.Date(missing)),
                         pct = 0L)
      rbind(.SD, gaps, fill = TRUE)
    } else {
      .SD
    }
  }, by = .(lake, data_type, variable)]

  filled
}

#' Tile plot of temporal coverage — met vs profile, per lake
#'
#' X-axis = month, Y-axis = lake, fill = % days with data.
#' One panel for Met, one for Profile.
#'
#' @param cov_dt data.table from calc_coverage()
#' @return ggplot object
plot_coverage_overview <- function(cov_dt) {

  # Regional council lookup for lake ordering
  council <- c(
    Karapiro = "WRC", Tutira = "HBRC", Wairarapa = "GWRC", Pounui = "GWRC",
    Rotoaira = "WDC", Rotorua = "BOPRC", Rotoiti = "BOPRC", Rotoehu = "BOPRC",
    Okareka = "BOPRC", Okaro = "BOPRC", Okataina = "BOPRC",
    Tarawera = "BOPRC", Rerewhakaaitu = "BOPRC",
    Hayes = "ORC", Wakatipu = "ORC", Wanaka = "ORC",
    Benmore = "ECan", Manapouri = "ES"
  )

  # Fill gaps and convert to %
  dt <- fill_coverage(cov_dt[variable == "ANY"])
  dt[, council := council[lake]]
  dt[is.na(council), council := "Other"]

  # Order lakes by council then lake name
  lake_order <- unique(dt[order(council, lake), lake])
  dt[, lake := factor(lake, levels = lake_order)]

  # Use geom_rect for gap-free tiles
  dt[, `:=`(xmin = ym,
            xmax = as.Date(format(ym + 32, "%Y-%m-01")))]

  ggplot(dt, aes(xmin = xmin, xmax = xmax,
                 ymin = as.numeric(lake) - 0.45,
                 ymax = as.numeric(lake) + 0.45,
                 fill = pct)) +
    geom_rect() +
    scale_fill_gradientn(
      name = "% days\nwith data",
      colours = c("#d73027", "#fc8d59", "#fee08b", "#d9ef8b", "#91cf60", "#1a9850"),
      limits = c(0, 100), oob = scales::squish, na.value = "grey90"
    ) +
    scale_y_continuous(breaks = seq_along(lake_order), labels = lake_order,
                       expand = c(0.02, 0)) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    facet_wrap(~data_type, ncol = 1) +
    labs(
      x = NULL, y = NULL,
      title = "Temporal data coverage",
      subtitle = "% of days with data per month — Met station and Water quality profiles"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold", size = 12),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

#' Tile plot of per-variable coverage for a single data type
#'
#' X-axis = month, Y-axis = variable, faceted by lake.
#'
#' @param cov_dt data.table from calc_coverage()
#' @param type character, "Met" or "Profile"
#' @return ggplot object
#' Compute daily presence/absence for met and profile dataok pivoting back to the 
#'
#' For each lake and day, returns TRUE/FALSE for met and profile presence.
#' Met present = >= 10 non-NA TmpAir readings. Profile present = >= 1 TmpWtr.
#' Fills all days from first record to today; missing days = FALSE.
#'
#' @param rds_dir character
#' @return data.table with lake, date, data_type, present (logical)
calc_daily_coverage <- function(rds_dir = "data/rds") {
  files <- list.files(rds_dir, pattern = "^full_.*[.]rds$", full.names = TRUE)

  rbindlist(lapply(files, function(f) {
    lk <- sub("full_", "", sub("[.]rds$", "", basename(f)))
    cat("Loading", lk, "...\n")
    x <- readRDS(f)

    results <- list()

    # Met: day present if >= 10 non-NA TmpAir readings
    met <- x[["raw"]][["met"]]
    if (!is.null(met) && "TmpAir" %in% names(met)) {
      met <- as.data.table(met)
      met[, date := as.Date(format(DateTime, "%Y-%m-%d"))]
      day_n <- met[!is.na(TmpAir), .N, by = date]
      met_days <- day_n[N >= 10, date]

      # Fill from first to today
      all_days <- seq(min(met_days), Sys.Date(), by = "day")
      dt <- data.table(lake = lk, date = all_days, data_type = "Met",
                        present = all_days %in% met_days)
      results <- c(results, list(dt))
    }

    # Profile: day present if >= 1 non-NA TmpWtr reading
    pro <- x[["raw"]][["pro"]]
    if (!is.null(pro) && "TmpWtr" %in% names(pro)) {
      pro <- as.data.table(pro)
      pro[, date := as.Date(format(DateTime, "%Y-%m-%d"))]
      pro_days <- pro[!is.na(TmpWtr), unique(date)]

      all_days <- seq(min(pro_days), Sys.Date(), by = "day")
      dt <- data.table(lake = lk, date = all_days, data_type = "Profile",
                        present = all_days %in% pro_days)
      results <- c(results, list(dt))
    }

    rbindlist(results)
  }))
}

#' Daily presence/absence plot — met vs profile, per lake
#'
#' X-axis = date, Y-axis = lake, coloured green (data) / red (no data).
#' White before start of record. One panel for Met, one for Profile.
#'
#' @param daily_dt data.table from calc_daily_coverage()
#' @return ggplot object
plot_daily_coverage <- function(daily_dt) {

  council <- c(
    Karapiro = "WRC", Tutira = "HBRC", Wairarapa = "GWRC", Pounui = "GWRC",
    Rotoaira = "WDC", Rotorua = "BOPRC", Rotoiti = "BOPRC", Rotoehu = "BOPRC",
    Okareka = "BOPRC", Okaro = "BOPRC", Okataina = "BOPRC",
    Tarawera = "BOPRC", Rerewhakaaitu = "BOPRC",
    Hayes = "ORC", Wakatipu = "ORC", Wanaka = "ORC",
    Benmore = "ECan", Manapouri = "ES"
  )

  dt <- copy(daily_dt)
  dt[, council := council[lake]]
  dt[is.na(council), council := "Other"]

  lake_order <- unique(dt[order(council, lake), lake])
  dt[, lake := factor(lake, levels = lake_order)]

  ggplot(dt, aes(x = date, y = lake, fill = present)) +
    geom_raster() +
    scale_fill_manual(
      name = "Data",
      values = c("TRUE" = "#1a9850", "FALSE" = "#d73027"),
      labels = c("TRUE" = "Present", "FALSE" = "Missing")
    ) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    facet_wrap(~data_type, ncol = 1) +
    labs(
      x = NULL, y = NULL,
      title = "Daily data presence / absence",
      subtitle = "Met: >= 10 air temp readings/day | Profile: >= 1 water temp reading/day"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold", size = 12),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      panel.background = element_rect(fill = "white", colour = NA),
      legend.position = "bottom"
    )
}

#' Summarise daily coverage by hydro year (July–June)
#'
#' For the first hydro year of each lake × data_type, the denominator is
#' days from first record to 30 June (not the full 365/366). Subsequent
#' years use the full year. Current (incomplete) year uses days up to today.
#'
#' @param daily_dt data.table from calc_daily_coverage()
#' @return data.table with lake, data_type, hydro_year (e.g. "2022-23"),
#'         days_present, days_possible, pct
calc_hydro_year_coverage <- function(daily_dt) {
  dt <- copy(daily_dt)

  # Hydro year: Jul 2022 – Jun 2023 = "2022-23"
  dt[, hy := ifelse(as.integer(format(date, "%m")) >= 7,
                    as.integer(format(date, "%Y")),
                    as.integer(format(date, "%Y")) - 1L)]
  dt[, hydro_year := paste0(hy, "-", substr(hy + 1L, 3, 4))]

  # Summarise
  out <- dt[, .(days_present  = sum(present),
                days_possible = .N),
            by = .(lake, data_type, hy, hydro_year)]

  out[, pct := round(100 * days_present / days_possible)]

  # Mark current (incomplete) hydro year
  current_hy <- ifelse(as.integer(format(Sys.Date(), "%m")) >= 7,
                       as.integer(format(Sys.Date(), "%Y")),
                       as.integer(format(Sys.Date(), "%Y")) - 1L)
  out[hy == current_hy, hydro_year := paste0(hydro_year, "*")]

  setorder(out, lake, data_type, hy)
  out
}

#' Tile plot of hydro-year coverage — lake × year, filled by % coverage
#'
#' @param hy_dt data.table from calc_hydro_year_coverage()
#' @return ggplot object
plot_hydro_year_coverage <- function(hy_dt) {

  council <- c(
    Karapiro = "WRC", Tutira = "HBRC", Wairarapa = "GWRC", Pounui = "GWRC",
    Rotoaira = "WDC", Rotorua = "BOPRC", Rotoiti = "BOPRC", Rotoehu = "BOPRC",
    Okareka = "BOPRC", Okaro = "BOPRC", Okataina = "BOPRC",
    Tarawera = "BOPRC", Rerewhakaaitu = "BOPRC",
    Hayes = "ORC", Wakatipu = "ORC", Wanaka = "ORC",
    Benmore = "ECan", Manapouri = "ES"
  )

  dt <- copy(hy_dt)
  dt[, council := council[lake]]
  dt[is.na(council), council := "Other"]

  lake_order <- unique(dt[order(council, lake), lake])
  dt[, lake := factor(lake, levels = lake_order)]

  # Order hydro years chronologically
  hy_levels <- sort(unique(dt$hydro_year))
  dt[, hydro_year := factor(hydro_year, levels = hy_levels)]

  ggplot(dt, aes(x = hydro_year, y = lake, fill = pct)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = pct), size = 2.8) +
    scale_fill_gradientn(
      name = "% days\nwith data",
      colours = c("#d73027", "#fc8d59", "#fee08b", "#d9ef8b", "#91cf60", "#1a9850"),
      limits = c(0, 100), oob = scales::squish, na.value = "white"
    ) +
    facet_wrap(~data_type, ncol = 1) +
    labs(
      x = NULL, y = NULL,
      title = "Annual data coverage (hydro years, Jul\u2013Jun)",
      subtitle = "% of days with data per hydro year"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold", size = 12),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

plot_coverage_by_variable <- function(cov_dt, type = "Met") {
  dt <- fill_coverage(cov_dt[data_type == type & variable != "ANY"])

  dt[, `:=`(xmin = ym,
            xmax = as.Date(format(ym + 32, "%Y-%m-01")))]

  ggplot(dt, aes(xmin = xmin, xmax = xmax,
                 ymin = as.numeric(factor(variable)) - 0.45,
                 ymax = as.numeric(factor(variable)) + 0.45,
                 fill = pct)) +
    geom_rect() +
    scale_fill_gradientn(
      name = "% days",
      colours = c("#d73027", "#fc8d59", "#fee08b", "#d9ef8b", "#91cf60", "#1a9850"),
      limits = c(0, 100), oob = scales::squish, na.value = "grey90"
    ) +
    scale_y_continuous(
      breaks = seq_along(levels(factor(dt$variable))),
      labels = levels(factor(dt$variable)),
      expand = c(0.02, 0)
    ) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    facet_wrap(~lake, ncol = 2, scales = "free_x") +
    labs(
      x = NULL, y = NULL,
      title = paste(type, "— per-variable temporal coverage"),
      subtitle = "% of days with valid data per month"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      strip.text = element_text(face = "bold", size = 9),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 8),
      panel.spacing = unit(0.5, "lines"),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}
