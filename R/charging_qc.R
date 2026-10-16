library(data.table)
library(ggplot2)
source("R/palette.R")

#' Extract all RadSWD and VBatLg data from lake RDS files
#'
#' Merges met and ctrl data by DateTime for all lakes, saves to
#' data/charging_raw.rds.
#'
#' @param rds_dir character path to RDS directory
#' @param out_path character path to save
#' @return data.table with lake, DateTime, RadSWD, VBatLg
build_charging_raw <- function(rds_dir = "data/rds",
                                out_path = "data/charging_raw.rds") {
  files <- list.files(rds_dir, pattern = "^full_.*[.]rds$", full.names = TRUE)

  raw <- rbindlist(lapply(files, function(f) {
    lk <- sub("full_", "", sub("[.]rds$", "", basename(f)))
    cat("Loading", lk, "...\n")
    x <- readRDS(f)

    met  <- x[["raw"]][["met"]]
    ctrl <- x[["raw"]][["ctrl"]]
    if (is.null(met) || is.null(ctrl)) return(NULL)
    if (!("RadSWD" %in% names(met)) || !("VBatLg" %in% names(ctrl))) return(NULL)

    met  <- as.data.table(met)[, .(DateTime, RadSWD)]
    ctrl <- as.data.table(ctrl)[, .(DateTime, VBatLg)]

    dt <- merge(met, ctrl, by = "DateTime", all = TRUE)
    dt[, lake := lk]
    dt[, .(lake, DateTime, RadSWD, VBatLg)]
  }), fill = TRUE)

  saveRDS(raw, out_path)
  cat("Saved", nrow(raw), "rows to", out_path, "\n")
  raw
}

#' Build daily charging summaries from raw data
#'
#' Filters to daylight hours (08:00-16:00), removes RadSWD outliers
#' (< 0 or > 1500), computes daily totals. Saves to data/daily_charging.rds.
#'
#' @param raw data.table from build_charging_raw() or loaded from RDS
#' @param out_path character path to save
#' @return data.table with lake, date, rad_sum, v_first, v_last, delta_v, v_max, n_vbat
build_daily_charging <- function(raw, out_path = "data/daily_charging.rds") {
  dt <- copy(raw)
  dt[, `:=`(date = as.Date(format(DateTime, "%Y-%m-%d")),
            hour = as.integer(format(DateTime, "%H")))]

  # Remove RadSWD outliers
  dt[RadSWD < 0 | RadSWD > 1500, RadSWD := NA]

  # For each lake-day, find first and last timestamp where RadSWD > 5
  # to define the actual daylight window
  daylight <- dt[!is.na(RadSWD) & RadSWD > 5,
                  .(dawn = min(DateTime), dusk = max(DateTime)),
                  by = .(lake, date)]

  # Trim to daylight window only
  dt <- dt[daylight, on = .(lake, date)][DateTime >= dawn & DateTime <= dusk]

  # Daily mean radiation (W/m2 averaged over 24h)
  # Sum radiation during daylight, then divide by number of intervals in 24h
  # Detect interval from data
  rad <- dt[, .(rad_sum = sum(RadSWD, na.rm = TRUE)),
            by = .(lake, date)]

  # Voltage at 5am (pre-dawn baseline)
  dt_all <- copy(raw)
  dt_all[, `:=`(date = as.Date(format(DateTime, "%Y-%m-%d")),
                hour = as.integer(format(DateTime, "%H")))]
  vbat_5am <- dt_all[!is.na(VBatLg) & hour == 5,
                      .(v_batt = median(VBatLg)),
                      by = .(lake, date)]

  # Voltage at end of day (last reading before midnight)
  vbat_night <- dt_all[!is.na(VBatLg) & hour >= 22,
                        .(v_night = VBatLg[which.max(DateTime)]),
                        by = .(lake, date)]

  # Max voltage during daylight
  vbat_max <- dt[!is.na(VBatLg),
                  .(v_max = max(VBatLg), n_vbat = .N),
                  by = .(lake, date)]

  vbat <- merge(vbat_5am, vbat_night, by = c("lake", "date"))
  vbat <- merge(vbat, vbat_max, by = c("lake", "date"))
  vbat[, delta_v := v_night - v_batt]

  daily <- merge(rad, vbat, by = c("lake", "date"))
  daily <- daily[rad_sum < 150000]

  saveRDS(daily, out_path)
  cat("Saved", nrow(daily), "rows to", out_path, "\n")
  daily
}

#' Compute monthly charging index from daily data
#'
#' @param daily data.table from build_daily_charging() or loaded from RDS
#' @param min_vbat_obs integer minimum VBatLg readings during daylight per day
#' @param vbat_float numeric exclude days where max voltage exceeds this
#' @param min_days_per_month integer minimum sunny days to compute monthly index
#' @return data.table with lake, ym, med_dv, n_days, rad_med, ci
calc_charging_index <- function(daily, min_vbat_obs = 3,
                                 vbat_float = 13.0, min_days_per_month = 3) {
  dt <- copy(daily)
  dt <- dt[n_vbat >= min_vbat_obs & v_max < vbat_float]

  rbindlist(lapply(unique(dt$lake), function(lk) {
    d <- dt[lake == lk]
    if (nrow(d) == 0) return(NULL)

    # Filter to sunny days: top 50% radiation for this site
    rad_median <- median(d$rad_sum, na.rm = TRUE)
    sunny <- d[rad_sum >= rad_median]
    if (nrow(sunny) == 0) return(NULL)

    # Monthly aggregation
    sunny[, ym := as.Date(format(date, "%Y-%m-01"))]
    monthly <- sunny[, .(med_dv  = median(delta_v, na.rm = TRUE),
                          n_days  = .N,
                          rad_med = median(rad_sum, na.rm = TRUE)),
                      by = ym]
    monthly <- monthly[n_days >= min_days_per_month]
    if (nrow(monthly) == 0) return(NULL)

    # Normalise: 1.0 = this site's typical performance
    monthly[, ci := med_dv / median(med_dv, na.rm = TRUE)]
    monthly[, lake := lk]
    monthly[, .(lake, ym, med_dv, n_days, rad_med, ci)]
  }), fill = TRUE)
}

#' Faceted time-series of charging index per lake
#'
#' @param ci_dt data.table from calc_charging_index()
#' @return ggplot object
plot_charging_index <- function(ci_dt) {
  dt <- copy(ci_dt)
  dt[, status := fifelse(ci < 0.5, "poor", "ok")]

  ggplot(dt, aes(x = ym, y = ci)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = 0.5, linetype = "dotted", colour = "#d73027") +
    geom_line(colour = "grey30", linewidth = 0.4) +
    geom_point(aes(colour = status), size = 1.5) +
    scale_colour_manual(values = c(ok = "#1a9850", poor = "#d73027"),
                         guide = "none") +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    coord_cartesian(ylim = c(0, 2.5)) +
    facet_wrap(~lake, ncol = 3, scales = "free_x") +
    labs(
      x = NULL, y = "Charging index (1.0 = typical)",
      title = "Solar panel charging performance",
      subtitle = "Monthly median voltage gain on sunny days, normalised per site. Red = below 50% of typical."
    ) +
    theme_minimal(base_size = 10) +
    theme(
      strip.text = element_text(face = "bold", size = 9),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7)
    )
}

#' Heatmap overview of charging index — all lakes x months
#'
#' @param ci_dt data.table from calc_charging_index()
#' @return ggplot object
plot_charging_overview <- function(ci_dt) {

  council <- c(
    Karapiro = "WRC", Tutira = "HBRC", Wairarapa = "GWRC", Pounui = "GWRC",
    Rotoaira = "WDC", Rotorua = "BOPRC", Rotoiti = "BOPRC", Rotoehu = "BOPRC",
    Okareka = "BOPRC", Okaro = "BOPRC", Okataina = "BOPRC",
    Tarawera = "BOPRC", Rerewhakaaitu = "BOPRC",
    Hayes = "ORC", Wakatipu = "ORC", Wanaka = "ORC",
    Benmore = "ECan", Manapouri = "ES"
  )

  dt <- copy(ci_dt)
  dt[, council := council[lake]]
  dt[is.na(council), council := "Other"]
  lake_order <- unique(dt[order(council, lake), lake])
  dt[, lake := factor(lake, levels = lake_order)]

  dt[, `:=`(xmin = ym,
            xmax = as.Date(format(ym + 32, "%Y-%m-01")))]

  ggplot(dt, aes(xmin = xmin, xmax = xmax,
                 ymin = as.numeric(lake) - 0.45,
                 ymax = as.numeric(lake) + 0.45,
                 fill = pmin(ci, 2))) +
    geom_rect() +
    scale_fill_gradientn(
      name = "Charging\nindex",
      colours = c("#d73027", "#fc8d59", "#fee08b", "#d9ef8b", "#91cf60", "#1a9850"),
      limits = c(0, 2), oob = scales::squish, na.value = "grey90"
    ) +
    scale_y_continuous(breaks = seq_along(lake_order), labels = lake_order,
                       expand = c(0.02, 0)) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(
      x = NULL, y = NULL,
      title = "Charging performance overview",
      subtitle = "Monthly charging index — red = poor charging on sunny days"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}
